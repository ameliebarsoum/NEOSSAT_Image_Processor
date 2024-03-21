import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
import cv2
from scipy.ndimage import shift, median_filter
from skimage.transform import resize
import random
from astropy.stats import mad_std
from ccdproc import CCDData, trim_image, combine
from astride import Streak
import warnings 
warnings.filterwarnings("ignore")
import logging
logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)

INPUT_PATH = '../n1fits/mission/image/outgoing/ASTRO/'
ALIGNED_PATH = 'aligned/'
DIFFERENCED_PATH = 'differenced/'
FLAGGED_PATH = 'flagged/'
STACKED_PATH = 'stacked/stacked_img.fits'
FLAGGED_ORIGINAL_PATH = 'flagged_original/'
CRH_PATH = 'cosmic_ray_hits/'
"""
Input: Path to fits file, and optionally: threshold, fwhm, sharplo, and roundlo for detection via DAOStarFinder
Description: Detects objects in the image using the DAOStarFinder algorithm
"""
def detect_objects(filepath, threshold=15, fwhm=10.0, sharplo=0.2, roundlo=-1.0, border=False, filter=True):

    try:
        with fits.open(filepath) as hdul:
            data = hdul[0].data  # Accessing the data from the primary HDU
    except FileNotFoundError:
        return None
    
    # Apply median filter to fits
    if filter:
        data = median_filter(data, 7)

    # Calculate background statistics
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)

    # DAOStarFinder algorithm for star detection
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold*std, sharplo=sharplo, roundlo=roundlo, exclude_border=border)
    sources = daofind(data - median)

    # Exclude sources along border for alignment
    if border and sources is not None:
        # Get image dimensions
        ny, nx = data.shape
        # Filter out sources within 5 pixels of the edge
        mask = (sources['xcentroid'] > 5) & (sources['xcentroid'] < nx - 5) & \
               (sources['ycentroid'] > 5) & (sources['ycentroid'] < ny - 5)
        filtered_sources = sources[mask]
        return filtered_sources

    return sources

"""
Parameters: Two lists of tuples representing x, y coordinates of sources
Returns: List 1 without any sources that are within a given number of pixels of any source in List 2
"""
def remove_duplicates(source1, source2, threshold=0.2):
    for i in range(len(source1)):
        for j in range(len(source2)):
            if np.sqrt((source1[i][0] - source2[j][0]) ** 2 + (source1[i][1] - source2[j][1]) ** 2) < threshold:
                source1.pop(i)
                break
    return source1

"""
List all distances between each base source and each img source
Returns list of tuples of the distance, id of base source, id of img source
"""
def find_all_distances_between(base_sources, img_sources, max_distance=15):
    distances = []
    for x_base, y_base, id_base in [(base_source['xcentroid'], base_source['ycentroid'], base_source['id']) for base_source in base_sources]:
        for x_img, y_img, id_img in [(img_source['xcentroid'], img_source['ycentroid'], img_source['id']) for img_source in img_sources]:
            distance = np.sqrt((x_base - x_img) ** 2 + (y_base - y_img) ** 2)
            if distance <= max_distance: 
                distances.append((distance, id_base, id_img))
    distances.sort(key=lambda x: x[0])
    return distances

"""
Find the closest pair of stars between two images such that
the distance between these two stars is the same as the distance between at least one unique other pair of stars.
This is done to ensure that the shift found can be validated against another pair, too.
"""
def find_closest_pair(tuples):
    closest_pair = None
    min_distance_pair = 0.1
    
    for i in range(len(tuples)):
        distance1, id1, id2 = tuples[i]
        for j in range(i + 1, len(tuples)):
            distance2, id3, id4 = tuples[j]
            if id3 not in (id1, id2) and id4 not in (id1, id2):
                distance_pair = abs(distance2 - distance1)
                # Shift should "match" at least one other pair to a closeness of min_distance_pair
                if distance_pair < min_distance_pair:
                    closest_pair = (id1, id2, id3, id4)
                    min_distance_pair = distance_pair

    return closest_pair

"""
@base_sources: sources from an image with which the other image will be aligned to
@img_sources: sources from the image to be aligned
@img_filename: filename of the image to be aligned
description: Aligns one image to another based on the brightest star, such that their stars overlap.
"""
def align(base_sources, img_sources, img_filename):
    
    distances = find_all_distances_between(base_sources, img_sources)
    
    if len(base_sources) == 1 or len(img_sources) == 1: # just take the smallest distance
        print("Only one star detected in one of the images. Using the closest star for alignment.")
        base_idx, img_idx, _, _ = distances[0]
    else:
        tuple = find_closest_pair(distances)
        if tuple is None:
            print("No matching stars found for alignment of " + img_filename + ". Not moving on with this file.")
            return -1

        base_idx, img_idx, second_base_idx, second_img_idx = tuple
    
    base_row = base_sources[base_sources['id'] == base_idx]
    img_row = img_sources[img_sources['id'] == img_idx]

    base_row2 = base_sources[base_sources['id'] == second_base_idx]
    img_row2 = img_sources[img_sources['id'] == second_img_idx]
    ref_star_coords_image1_d = base_row2['xcentroid'][0], base_row2['ycentroid'][0]
    ref_star_coords_image2_d = img_row2['xcentroid'][0], img_row2['ycentroid'][0]

    ref_star_coords_image1 = base_row['xcentroid'][0], base_row['ycentroid'][0]
    ref_star_coords_image2 = img_row['xcentroid'][0], img_row['ycentroid'][0]

    # Calculate the shift needed to align the images
    shift_x = ref_star_coords_image1[0] - ref_star_coords_image2[0]
    shift_y = ref_star_coords_image1[1] - ref_star_coords_image2[1]

    # Save png with stars used for alignment circled
    save_flagged_image(INPUT_PATH + img_filename, ALIGNED_PATH + img_filename.split(".")[0] + ".png", [ref_star_coords_image1_d+ref_star_coords_image2_d +ref_star_coords_image1+ ref_star_coords_image2])

    # Shift one image relative to the other using interpolation
    with fits.open(INPUT_PATH + img_filename) as hdul: 
        shifted_data = shift(hdul[0].data, (shift_y, shift_x), mode='nearest')
        hdul[0].data = shifted_data
        
        hdu = fits.PrimaryHDU(shifted_data, header=hdul[0].header)
        hdul_out = fits.HDUList([hdu])
        hdul_out.writeto(ALIGNED_PATH + img_filename, overwrite=True)

    print(f"Aligning images with a shift of ({shift_x}, {shift_y}) pixels as {ALIGNED_PATH + img_filename}.")
    return (shift_x, shift_y) 

"""
Crops all images in the directory to the minimum sized image of the set
"""
def crop_all(ALIGNED_PATH):
    # Get all the images in the directory and load them as CCDData objects
    ccds = [CCDData.read(ALIGNED_PATH + filename, unit='adu') for filename in os.listdir(ALIGNED_PATH) if filename.lower().endswith('.fits')]
    
    # Identify the smallest image size within the specified range
    sizes = [ccd.shape for ccd in ccds]
    
    # Find minimum width, height for cropping
    min_size = (min(sizes, key=lambda x: x[0])[0], min(sizes, key=lambda x: x[1])[1])
    print(f"Cropping images to the minimum size of the set: {min_size}")
    
    # Replace original images with new cropped images
    for filename in os.listdir(ALIGNED_PATH):
        if filename.lower().endswith('.fits'):
            # Open as CCDData to maintain metadata
            ccd = CCDData.read(ALIGNED_PATH + filename, unit='adu')
            trimmed_ccd = trim_image(ccd, fits_section=f'[1:{min_size[1]}, 1:{min_size[0]}]')
            
            # Write the trimmed image back to disk
            trimmed_ccd.write(ALIGNED_PATH + filename, overwrite=True)
    print(f"Cropped files to the minimum size of the set: {trimmed_ccd.shape}")

"""
Perform image stacking on a set of FITS files to reduce noise for pixel differencing
"""
def image_stacker(ALIGNED_PATH):
    paths = [ALIGNED_PATH + filename for filename in os.listdir(ALIGNED_PATH) if filename.lower().endswith('.fits')]
    ccds = [CCDData.read(path, unit="adu") for path in paths]
    
    combined_image = combine(ccds,
                             output_file=STACKED_PATH, # Comment out if you don't want to save the unblurred stacked image
                             method='average',
                             sigma_clip=True,
                             sigma_clip_low_thresh=3,
                             sigma_clip_high_thresh=3,
                             sigma_clip_func=np.median,
                             sigma_clip_dev_func=mad_std,
                             mem_limit=350e6,
                             overwrite_output=True)
##### Comment/Uncomment to apply Gaussian blur to the combined image #####
    # combined_image_np = combined_image.data.astype(np.float32)
    # blurred_image = cv2.GaussianBlur(combined_image_np, (5, 5), 2)
    # hdu = fits.PrimaryHDU(blurred_image, header=ccds[0].header)
    # hdul_out = fits.HDUList([hdu])
    # hdul_out.writeto(STACKED_PATH, overwrite=True)
    
    # return blurred_image
# #####
    return combined_image

"""
Perform pixel differencing on two FITS files
"""
def pixel_difference(image1_file, image2_file, idx):
    try:
        with fits.open(image1_file) as hdul1, fits.open(image2_file) as hdul2:
            data1 = hdul1[0].data
            data2 = hdul2[0].data

            # Resize images to a common size if they are not the same size
            if (data1.shape[0] != data2.shape[0] or data1.shape[1] != data2.shape[1]):
                print("Adjusting size for pixel differencing")
                min_shape = (min(data1.shape[0], data2.shape[0]), min(data1.shape[1], data2.shape[1]))
                data1 = resize(data1, min_shape, mode='constant')
                data2 = resize(data2, min_shape, mode='constant')

            # Perform pixel-wise difference
            diff_data = np.abs(data1 - data2)  # Calculates absolute difference

            # Save the difference image to a new FITS file
            output_filename = DIFFERENCED_PATH + image1_file.split("/")[-1].split(".")[0] + "_" + idx + ".fits"
            hdu = fits.PrimaryHDU(diff_data, header=hdul1[0].header)
            hdul_out = fits.HDUList([hdu])
            hdul_out.writeto(output_filename, overwrite=True)

    except FileNotFoundError:
        print(f"Could not find {image1_file} or {image2_file}. Previous steps may not have been completed.")
        return
    
"""
Performs pixel differencing with the stacked image found at STACKED_PATH
"""
def pixel_difference_with_stack(image1_filename):
    image1_file = os.path.join(ALIGNED_PATH, image1_filename)
    try:
        with fits.open(image1_file) as hdul1, fits.open(STACKED_PATH) as hdul2:
            data1 = hdul1[0].data
            data2 = hdul2[0].data

            # Images are not always the same size. Not sure if this is a result of the cleaning algorithm or not.
            if (data1.shape[0] != data2.shape[0] or data1.shape[1] != data2.shape[1]):
                print("Adjusting size for pixel differencing")
                min_shape = (min(data1.shape[0], data2.shape[0]), min(data1.shape[1], data2.shape[1]))
                data1 = resize(data1, min_shape, mode='edge') 
                data2 = resize(data2, min_shape, mode='edge')

            # Perform pixel-wise difference
            diff_data = data1 - data2 

            # Save the difference image to a new FITS file
            output_filename = DIFFERENCED_PATH + image1_filename
            hdu = fits.PrimaryHDU(diff_data, header=hdul1[0].header)
            hdul_out = fits.HDUList([hdu])
            hdul_out.writeto(output_filename, overwrite=True)

    except FileNotFoundError:
        print(f"Could not find {image1_file} or stacked image. Previous steps may not have been completed.")
        return

"""
Applies Gaussian blur
"""
def blur(path):
    # Open the FITS file
    with fits.open(path) as hdul:
        image = hdul[0].data  # Accessing the data from the primary HDU
        blurred_image = cv2.GaussianBlur(image.astype(np.float32), (5, 5), 2)

        # Save the background-subtracted and darkened image to a new FITS file
        hdu = fits.PrimaryHDU(blurred_image, header=hdul[0].header)
        hdul_out = fits.HDUList([hdu])
        hdul_out.writeto(path, overwrite=True)

"""
Use astride to detect cosmic ray hits in the image
"""
def streak_detection(path):

    max_length = 500
    min_length = 5

    try:
        output_path = CRH_PATH + path.split("/")[1].split(".")[0]
        streak = Streak(path, min_points=40, area_cut=40, connectivity_angle=0.01, output_path=output_path)
    except FileNotFoundError:
        return None
    streak.detect()

    # Remove very short, very long, very thick streaks
    for streak_instance in streak.streaks: 
        x_min = streak_instance.get('x_min')
        x_max = streak_instance.get('x_max')
        y_min = streak_instance.get('y_min')
        y_max = streak_instance.get('y_max')
        length = ((x_max - x_min) ** 2 + (y_max - y_min) ** 2) ** 0.5

        if length > max_length or length < min_length or streak_instance.get('area') > streak_instance.get('perimeter'):
            streak.streaks.remove(streak_instance)

    if len(streak.streaks) == 0:
        return None
    
    print(f"Detected {len(streak.streaks)} cosmic ray hit(s) in {path}.")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    streak.write_outputs()
    streak.plot_figures()

    # Return list of tuples of x_center, y_center of cosmic ray hits
    return [(streak_instance.get('x_center'), streak_instance.get('y_center')) for streak_instance in streak.streaks]

"""
Parameters: input path to fits image, output path, sources
Save image with circled sources to flagged directory
"""
def save_flagged_image(input_path, output_path, sources):
    with fits.open(input_path) as hdul:
        image_data = hdul[0].data

    # Plot the image data
    plt.figure(figsize=(10, 10))
    plt.imshow(image_data, cmap='gray', origin='lower', norm=plt.Normalize(vmin=np.percentile(image_data, 5), vmax=np.percentile(image_data, 95)))

    for i in range(0,len(sources)):
        plt.scatter(sources[i][0], sources[i][1], s=100, facecolors='none', edgecolors='r')

    # Save the image to the output path as PNG
    plt.axis('off')  # Optional: Turn off the axis
    plt.savefig(output_path, bbox_inches='tight', pad_inches=0)
    plt.close()

"""
Helper for choosing random base image for pixel diffencing 
"""
def get_random_files(directory, num_files):
    files = [filename for filename in os.listdir(directory) if filename.lower().endswith('.fits')]
    return random.sample(files, min(num_files, len(files)))


if __name__ == "__main__":

    fits_files = [filename for filename in os.listdir(INPUT_PATH) if filename.lower().endswith('.fits')]
    # Delete prev test files with os.remove here if needed
    if len(fits_files) < 2:
        print("Not enough images to perform pixel differencing.")
        exit()
  
    if not os.path.exists(ALIGNED_PATH):
        os.makedirs(ALIGNED_PATH)
    if not os.path.exists(DIFFERENCED_PATH):
        os.makedirs(DIFFERENCED_PATH) 
    if not os.path.exists(FLAGGED_PATH):
        os.makedirs(FLAGGED_PATH) 
    if not os.path.exists("stacked"):
        os.makedirs("stacked")
    if not os.path.exists(FLAGGED_ORIGINAL_PATH):
        os.makedirs(FLAGGED_ORIGINAL_PATH)
    if not os.path.exists(CRH_PATH):
        os.makedirs(CRH_PATH)

    # Cosmic ray hit detection (works best on raw images)
    for filename in fits_files:
        streak_coords = streak_detection(INPUT_PATH + filename)       

    # Alignment
    source_tuples = []
    base_source = None
    max_stars = -1
    for filename in fits_files:
        num_stars = 0
        cur_source = detect_objects(INPUT_PATH + filename, threshold=25, border=True) # Too many points = inaccurate alignment
        num_stars = len(cur_source)
        source_tuples.append((cur_source, filename))
        print(f"Detected {num_stars} objects in {filename}.")
    
    # Select the source with the median detected objects as the base for alignment
    source_tuples.sort(key=lambda x: len(x[0]))
    base_source = source_tuples[len(source_tuples) // 2][0]
    print(f"Using {source_tuples[len(source_tuples) // 2][1]} having {len(base_source)} detected sources as the base for alignment .")

    for source in source_tuples:
        align(base_source, source[0], source[1])

    # Crop to same size
    crop_all(ALIGNED_PATH)
    
    # Pixel Differencing
    if len(fits_files) == 2 or len(fits_files) == 3:
        for i, base_filename in enumerate(fits_files):
            for j, comparison_filename in enumerate(fits_files):
                if i != j: 
                    base_path = os.path.join(ALIGNED_PATH, base_filename)
                    comparison_path = os.path.join(ALIGNED_PATH, comparison_filename)
                    pixel_difference(base_path, comparison_path, str(i) + str(j))
    else: 
        # Create a stacked image
        stacked_image = image_stacker(ALIGNED_PATH=ALIGNED_PATH)
        stacked_image_sources = detect_objects(STACKED_PATH)
        stacked_image_coordinates = [(source['xcentroid'], source['ycentroid']) for source in stacked_image_sources] # this is probs wrong

        # Perform pixel differencing with all aligned images
        for filename in fits_files:
            # Apply Gaussian blur here if needed
            pixel_difference_with_stack(filename)
            sources = detect_objects(DIFFERENCED_PATH + filename, threshold=50, fwhm=25, sharplo=0.5)
            if sources is None or len(sources) < 1:
                print(f"No objects detected in {filename}.")
                continue
            print(f"Detected {len(sources)} objects in {filename}.")
            print(f"Coordinates: {[(source['xcentroid'], source['ycentroid']) for source in sources]}")

            sources_coordinates = [(source['xcentroid'], source['ycentroid']) for source in sources]
            filtered_sources = remove_duplicates(sources_coordinates, stacked_image_coordinates)
        
            if filtered_sources is not None and len(sources) > 0:
                print(f"Flagged {filename} for classification.")            
                save_flagged_image(DIFFERENCED_PATH + filename, FLAGGED_PATH + filename.split(".")[0] + ".png", sources_coordinates)
                save_flagged_image(INPUT_PATH + filename, FLAGGED_ORIGINAL_PATH + filename.split(".")[0] + ".png", sources_coordinates)     
