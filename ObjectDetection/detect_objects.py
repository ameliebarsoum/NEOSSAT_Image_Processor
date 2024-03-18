import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
import cv2
from scipy.ndimage import shift
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
def detect_objects(filepath, threshold=15, fwhm=10.0, sharplo=0.2, roundlo=-1.0):

    with fits.open(filepath) as hdul:
        data = hdul[0].data  # Accessing the data from the primary HDU

    # Calculate background statistics
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)

    # DAOStarFinder algorithm for star detection
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold*std, sharplo=sharplo, roundlo=roundlo)
    sources = daofind(data - median)
    return sources

"""
Helper to deal with DAOStarFinder QTable
"""
def qtable_to_tuples(qtable):
    tuples_list = []
    for row in qtable:
        x_centroid = row['xcentroid']
        y_centroid = row['ycentroid']
        obj_id = row['id']
        tuples_list.append((x_centroid, y_centroid, obj_id))
    return tuples_list

"""
List all distances between each base source and each img source
Returns list of tuples of the distance, id of base source, id of img source
"""
def find_all_distances_between(base_sources, img_sources):
    base_tuples = qtable_to_tuples(base_sources)
    img_tuples = qtable_to_tuples(img_sources)
    distances = []
    for base in base_tuples:
        for img in img_tuples:
            distance = np.sqrt((base[0] - img[0]) ** 2 + (base[1] - img[1]) ** 2)
            distances.append((distance, base[2], img[2]))
    distances.sort(key=lambda x: x[0])
    return distances

"""
Find the closest pair of stars between two images such that
the distance between these two stars is the same as the distance between at least one unique other pair of stars.
This is done to ensure that the shift found can be validated against another pair, too.
"""
def find_closest_pair(tuples):
    closest_pair = None
    min_distance = float('inf')
    
    for i in range(len(tuples)):
        distance1, id1, id2 = tuples[i]
        
        for j in range(i + 1, len(tuples)):
            distance2, id3, id4 = tuples[j]
            
            if abs(distance2 - distance1) < 1 and not (id1 == id3 or id1 == id4 or id2 == id3 or id2 == id4):
                if distance1 < min_distance:
                    min_distance = distance1
                    closest_pair = (distance1, id1, id2)
                if distance2 < min_distance:
                    min_distance = distance2
                    closest_pair = (distance2, id3, id4)
                    
    return closest_pair

"""
@base_sources: sources from an image with which the other image will be aligned t
@img_sources: sources from the image to be aligned
@img_filename: filename of the image to be aligned
description: Aligns one image to another based on the brightest star, such that their stars overlap.
"""
def align(base_sources, img_sources, img_filename):
    
    distances = find_all_distances_between(base_sources, img_sources)
    
    if len(base_sources) == 1 or len(img_sources) == 1: # just take the smallest distance
        print("Only one star detected in one of the images. Using the closest star for alignment.")
        min_distance, base_idx, img_idx = distances[0]
    else:
        tuple = find_closest_pair(distances)
        if tuple is None:
            print("No matching stars found for alignment of " + img_filename + ". Not moving on with this file.")
            return -1

        min_distance, base_idx, img_idx = tuple
    
    if min_distance > 100:
        print("No matching stars found for alignment of " + img_filename + ". Not moving on with this file.")
        return -1
    
    base_row = base_sources[base_sources['id'] == base_idx]
    img_row = img_sources[img_sources['id'] == img_idx]

    ref_star_coords_image1 = base_row['xcentroid'][0], base_row['ycentroid'][0]
    ref_star_coords_image2 = img_row['xcentroid'][0], img_row['ycentroid'][0]

    # Calculate the shift needed to align the images
    shift_x = ref_star_coords_image1[0] - ref_star_coords_image2[0]
    shift_y = ref_star_coords_image1[1] - ref_star_coords_image2[1]

    # Shift one image relative to the other using interpolation
    with fits.open(INPUT_PATH + img_filename) as hdul:  # Change to GAUSSIANBLUR_PATH if using blur_background
        shifted_data = shift(hdul[0].data, (shift_y, shift_x), mode='nearest')
        hdul[0].data = shifted_data
        
        hdu = fits.PrimaryHDU(shifted_data, header=hdul[0].header)
        hdul_out = fits.HDUList([hdu])
        hdul_out.writeto(ALIGNED_PATH + img_filename, overwrite=True)

    print(f"Aligning images with a shift of ({shift_x}, {shift_y}) pixels as {ALIGNED_PATH + img_filename}.")
    return 0

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
Applies Gaussian blur to the background (non-sources) of image located at @filename based on peaks given by @sources
"""
def blur_background(filename, sources):
    # Open the FITS file
    with fits.open(INPUT_PATH + filename) as hdul:
        image = hdul[0].data  # Accessing the data from the primary HDU

        # Create a mask to exclude areas around bright points
        mask = np.ones_like(image, dtype=np.float32)  # Initialize mask with ones (full image)
        blur_radius = 4  # Define the radius around bright points to exclude from blurring

        x_positions = sources['xcentroid']
        y_positions = sources['ycentroid']
        for i in range(len(sources)):
            x = x_positions[i]
            y = y_positions[i]
            # Exclude a circular region around each bright point from the mask
            yy, xx = np.ogrid[:image.shape[0], :image.shape[1]]
            mask_value = np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * blur_radius ** 2))
            mask *= (1 - mask_value)

        # Apply Gaussian blur to the image using the inverted mask
        inverted_mask = 1 - mask
        blurred_image = cv2.GaussianBlur(image.astype(np.float32), (5, 5), 2)
        
        # Combine the original image with the darkened blurred image using the masks
        final_image = image * inverted_mask + blurred_image * mask
        final_image=  image
        # Save the background-subtracted and darkened image to a new FITS file
        output_filename = ALIGNED_PATH + filename
        hdu = fits.PrimaryHDU(final_image, header=hdul[0].header)
        hdul_out = fits.HDUList([hdu])
        hdul_out.writeto(output_filename, overwrite=True)

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
#####
    return combined_image

"""
Perform pixel differencing on two FITS files
"""
def pixel_difference(image1_filename, image2_filename):

    image1_file = os.path.join(ALIGNED_PATH, image1_filename)
    image2_file = os.path.join(ALIGNED_PATH, image2_filename)

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
            output_filename = DIFFERENCED_PATH + image2_filename
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
            diff_data = np.abs(data1 - data2)  # Calculate absolute difference

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

    straightness_threshold = 0.8
    max_length = 500
    min_length = 15

    output_path = CRH_PATH + path.split("/")[1].split(".")[0]
    streak = Streak(path, min_points=40, area_cut=40, connectivity_angle=0.01, output_path=output_path)
    streak.detect()

    # Remove very short, very long, very thick, and not-straight streaks
    for streak_instance in streak.streaks: 
        x_min = streak_instance.get('x_min')
        x_max = streak_instance.get('x_max')
        y_min = streak_instance.get('y_min')
        y_max = streak_instance.get('y_max')
        length = ((x_max - x_min) ** 2 + (y_max - y_min) ** 2) ** 0.5

        # Calculate straightness ratio
        direct_distance = max(x_max - x_min, y_max - y_min)
        straightness = direct_distance / length if length > 0 else 0

        if length > max_length or length < min_length or straightness < straightness_threshold \
          or streak_instance.get('area') > streak_instance.get('perimeter'):
            streak.streaks.remove(streak_instance)

    if len(streak.streaks) == 0:
        return None
    
    print(f"Detected {len(streak.streaks)} cosmic ray hit(s) in {path}.")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    streak.write_outputs()
    streak.plot_figures()
    return streak

"""
Visualizer for sources on a FITS image
"""
def sources_with_visualizer(path, sources):
    # Open the FITS file
    with fits.open(path) as hdul:
        data = hdul[0].data  # Accessing the data from the primary HDU

    # Display the image
    # First plot: FITS Image
    plt.figure(figsize=(12, 6))

    # Create subplot 1 (1 row, 2 columns, first plot)
    plt.subplot(1, 2, 1)
    plt.imshow(data, cmap='viridis', origin='lower')
    plt.colorbar(label='Pixel Value')
    plt.title(f'FITS file: {path}')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')

    # Second plot: Detected Stars
    # Create subplot 2 (1 row, 2 columns, second plot)
    plt.subplot(1, 2, 2)
    plt.imshow(data, cmap='viridis', origin='lower', vmax=500)
    plt.colorbar(label='Pixel Value')
    plt.scatter(sources['xcentroid'], sources['ycentroid'], s=50, edgecolor='red', facecolor='none', marker='o')
    plt.title('Detected Stars')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')

    plt.tight_layout()

    # Inspect pixel values interactively
    def onclick(event):
        x, y = int(event.xdata), int(event.ydata)
        print(f"Pixel value at (x={x}, y={y}): {data[y, x]}")  # Accessing pixel value

    fig, ax = plt.subplots()
    ax.imshow(data, cmap='viridis', origin='lower')
    plt.title('Click on the image to get pixel value')

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()

    return sources

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
    for filename in os.listdir(INPUT_PATH):
        if filename.lower().endswith('.fits'):
            streak_detection(INPUT_PATH + filename)

    # Alignment
    source_tuples = []
    base_source = None
    for filename in fits_files:
        num_stars = 0
        threshold = 15
        # Find the file with the most detected objects to use as the base for optimal alignment
        max_stars = -1
        while num_stars < 1 and threshold >= 1:
            cur_source = detect_objects(INPUT_PATH + filename)
            num_stars = len(cur_source)
            threshold -= 3
        source_tuples.append((cur_source, filename))
        if num_stars > max_stars:
            max_stars = num_stars
            base_source = cur_source

    for source in source_tuples:
        align(base_source, source[0], source[1])

    if base_source is None:
        print("No FITS files found in the specified directory.")
        exit()

    # Crop to same size
    crop_all(ALIGNED_PATH)
    
    # Pixel Differencing
    if len(fits_files) < 4:
        fits_files = [filename for filename in os.listdir(INPUT_PATH) if filename.lower().endswith('.fits')]
        # Iterate through each file and compare it to every other file
        for i, base_filename in enumerate(fits_files):
            for j, comparison_filename in enumerate(fits_files):
                if i != j:  # Avoid comparing a file to itself
                    base_path = os.path.join(ALIGNED_PATH, base_filename)
                    comparison_path = os.path.join(ALIGNED_PATH, comparison_filename)
                    pixel_difference(base_path, comparison_path)
    else: 
        # Create a stacked image
        stacked_image = image_stacker(ALIGNED_PATH=ALIGNED_PATH)
        # Perform pixel differencing with all blurred & aligned images
        for filename in os.listdir(ALIGNED_PATH):
            if filename.lower().endswith('.fits'):
                pixel_difference_with_stack(filename)
                # Gaussian Blur of differenced images before object detection
                blur(DIFFERENCED_PATH + filename)

    # Object detection on differenced images
    for filename in os.listdir(DIFFERENCED_PATH):
        if filename.lower().endswith('.fits'):
            sources = detect_objects(DIFFERENCED_PATH + filename, threshold=60, fwhm=17, sharplo=0.5, roundlo=-0.8)
            # Move files with stars detected into flagged folder
            if sources is not None and len(sources) > 0:
                print(f"Flagged {filename} for classification.")
                os.rename(DIFFERENCED_PATH + filename, FLAGGED_PATH + filename)
                os.rename(INPUT_PATH + filename, FLAGGED_ORIGINAL_PATH + filename)      
