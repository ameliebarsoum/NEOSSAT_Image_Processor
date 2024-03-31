import os
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils.detection import DAOStarFinder
from collections import defaultdict
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import shift, median_filter
from skimage.transform import resize
from astropy.stats import mad_std
from ccdproc import CCDData, trim_image, combine
from matplotlib.path import Path
from astride import Streak
import shutil
import warnings 
warnings.filterwarnings("ignore")
import logging
logger = logging.getLogger()

INPUT_PATH = "/Users/rebeccamizrahi/Documents/McGill/Capstone/NEOSSAT_Image_Processor/ObjectDetection/13P_copy/"#'n1fits/mission/image/outgoing/ASTRO/' # Relative path from root directory
ALIGNED_PATH = 'aligned/'
DIFFERENCED_PATH = 'differenced/'
FLAGGED_PATH = 'flagged/'
STACKED_PATH = 'stacked/'
FLAGGED_ORIGINAL_PATH = 'flagged_original/'
FLAGGED_ORIGINAL_FITS_PATH = 'flagged_original_fits/'
CRH_PATH = 'cosmic_ray_hits/'
FILTERED_FLAGGED_ORIGINAL_PATH = 'filtered_flagged_original/'

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
Parameters: 
- x, y coordinates of a point
- A list of sources from DAOStarFinder representing x, y coordinates of sources
- threshold: The maximum distance between the point and any source in the list for the function to return True
Returns: True of x, y coordinates are within threshold distance of any source in the list, False otherwise
"""
def is_duplicate(x_coord, y_coord, source_list, threshold=0.05):
    for source in source_list:
        if np.sqrt((x_coord - source['xcentroid']) ** 2 + (y_coord - source['ycentroid']) ** 2) < threshold:
            return True
    return False

"""
List all distances between each base source and each img source
Returns list of tuples of the distance, id of base source, id of img source
"""
def find_all_distances_between(base_sources, img_sources, max_distance=200):
    distances = []
    for x_base, y_base, id_base in [(base_source['xcentroid'], base_source['ycentroid'], base_source['id']) for base_source in base_sources]:
        for x_img, y_img, id_img in [(img_source['xcentroid'], img_source['ycentroid'], img_source['id']) for img_source in img_sources]:
            distance = np.sqrt((x_base - x_img) ** 2 + (y_base - y_img) ** 2)
            if distance <= max_distance: 
                distances.append((distance, id_base, id_img))
    return distances

"""
Find the closest pair of stars between two images such that
the distance between these two stars is the same as the 
distance between at least one unique other pair of stars.
This is done to ensure that the shift found can be validated against another pair, too.
"""
def find_closest_pair(tuples):
    closest_pair = None
    min_distance_pair = 0.1
    max_distance = 15
    
    for i in range(len(tuples)):
        distance1, id1, id2 = tuples[i]
        for j in range(i + 1, len(tuples)):
            distance2, id3, id4 = tuples[j]
            if id3 not in (id1, id2) and id4 not in (id1, id2) and distance2 < max_distance and distance1 < max_distance:
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
TODO: improvement. Currently only works on images taken in the same hour. If 'background' stars move too much, this will not work.
"""
def align(base_sources, img_sources, img_filename):

    if len(base_sources) == 1 or len(img_sources) == 1: # Taking smallest distance often results in bad alignment
            print("Only one star found for alignment of " + img_filename + ". Skipping this image.")
            return -1
    
    with fits.open(INPUT_PATH + img_filename) as hdul:         
        distances = find_all_distances_between(base_sources, img_sources)
        tuple = find_closest_pair(distances)

        if tuple is None:
            print("No pair of stars found for alignment of " + img_filename + ". Skipping this image.")
            return -1
        else:
            base_id, img_id, second_base_id, second_img_id = tuple
            print(f"Aligning {img_filename} to {base_id} and {img_id}.")
    
        base_row = base_sources[base_sources['id'] == base_id]
        img_row = img_sources[img_sources['id'] == img_id]

        ref_star_coords_image1 = base_row['xcentroid'][0], base_row['ycentroid'][0]
        ref_star_coords_image2 = img_row['xcentroid'][0], img_row['ycentroid'][0]

        # Calculate the shift needed to align the images
        shift_x = ref_star_coords_image1[0] - ref_star_coords_image2[0]
        shift_y = ref_star_coords_image1[1] - ref_star_coords_image2[1]

        # Shift one image relative to the other using interpolation
        shifted_data = shift(hdul[0].data, (shift_y, shift_x), mode='nearest')
        hdul[0].data = shifted_data
        
        
        hdu = fits.PrimaryHDU(shifted_data, header=hdul[0].header)
        hdul_out = fits.HDUList([hdu])
        hdul_out.writeto(ALIGNED_PATH + img_filename, overwrite=True)

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
    
    # Replace original images with new cropped images
    for filename in os.listdir(ALIGNED_PATH):
        if filename.lower().endswith('.fits'):
            # Open as CCDData to maintain metadata
            ccd = CCDData.read(ALIGNED_PATH + filename, unit='adu')
            trimmed_ccd = trim_image(ccd, fits_section=f'[1:{min_size[1]}, 1:{min_size[0]}]')
            
            # Write the trimmed image back to disk
            trimmed_ccd.write(ALIGNED_PATH + filename, overwrite=True)

"""
Perform image stacking on a set of FITS files to reduce noise for pixel differencing
"""
def image_stacker(input_path, output_path):
    paths = [input_path + filename for filename in os.listdir(input_path) if filename.lower().endswith('.fits')]
    ccds = [CCDData.read(path, unit="adu") for path in paths]
    
    combined_image = combine(ccds,
                             output_file=output_path, 
                             method='average',
                             sigma_clip=True,
                             sigma_clip_low_thresh=3,
                             sigma_clip_high_thresh=3,
                             sigma_clip_func=np.median,
                             sigma_clip_dev_func=mad_std,
                             mem_limit=350e6,
                             overwrite_output=True)

    return combined_image
    
"""
Performs pixel differencing with the stacked image found at STACKED_PATH
"""
def pixel_difference_with_stack(image1_path, stacked_image_path):

    try:
        with fits.open(image1_path) as hdul1, fits.open(stacked_image_path) as hdul2:
            data1 = hdul1[0].data
            data2 = hdul2[0].data

            # Images are not always the same size. Not sure if this is a result of the cleaning algorithm or not.
            if (data1.shape[0] != data2.shape[0] or data1.shape[1] != data2.shape[1]):
                min_shape = (min(data1.shape[0], data2.shape[0]), min(data1.shape[1], data2.shape[1]))
                data1 = resize(data1, min_shape, mode='edge') 
                data2 = resize(data2, min_shape, mode='edge')

            # Perform pixel-wise difference
            diff_data = data1 - data2 

            # Save the difference image to a new FITS file
            output_filename = DIFFERENCED_PATH + image1_path.split("/")[-1].split(".")[0] + ".fits"
            hdu = fits.PrimaryHDU(diff_data, header=hdul1[0].header)
            hdul_out = fits.HDUList([hdu])
            hdul_out.writeto(output_filename, overwrite=True)

    except FileNotFoundError:
        print(f"Could not find {image1_path} or stacked image. Previous steps may not have been completed.")
        return

"""
Use astride to detect cosmic ray hits in the image
"""
def streak_detection(path):
    
    max_length = 500
    min_length = 5

    try:
        with fits.open(path) as hdul:
            data = hdul[0].data  # Accessing the data from the primary HDU

        # Flip image to match orientation of the streaks
        data = np.flip(data, axis=0)

        hdu = fits.PrimaryHDU(data, header=hdul[0].header)
        hdul_out = fits.HDUList([hdu])
        hdul_out.writeto(path, overwrite=True)

        output_path = CRH_PATH + path.split("/")[-1].split(".")[0]
        streak = Streak(path, min_points=30, area_cut=40, connectivity_angle=0.01, output_path=output_path)
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

    return streak.streaks

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

    print(f"Flagged image saved to {output_path}.")

"""
Description: Check if a point is inside any of the detected streaks.

Parameters:
- x_centroid, y_centroid: Coordinates of the source detected by DAOStarFinder.
- raw_borders: List of borders for each detected streak from ASTRiDE, 
    where each border is a list of (x, y) tuples that defines the polygon.
    
Returns: True if the source is within any streak, False otherwise.
"""
def is_source_in_streak(x_centroid, y_centroid, raw_borders):

    for border in raw_borders:
        # Create a polygon path from streak borders
        x_coords = border['x']
        y_coords = border['y']
        border_points = list(zip(x_coords, y_coords)) # List of (x, y) tuples
        path = Path(border_points)
        # Check if the source's centroid is inside the path
        if path.contains_point((x_centroid, y_centroid)):
            return True  # Source is inside a streak
    return False 

def find_moving_objects(SOURCES, closeness_threshold=5):
    """
    Groups detected sources across images by suspected objects based on spatial consistency.
    
    Parameters:
    - source_detections: Dict, where each key is an image identifier and each value is a list of (x, y) tuples for detected sources.
    - closeness_threshold: Maximum distance between sources to consider them as part of the same object.
    """
    objects = defaultdict(list) # Map object ID to list of (filename, image_id, x, y) tuples
    object_id = 0
    euclidean_distance = lambda x1, y1, x2, y2: np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

    for fname, sources in SOURCES.items():
        if sources is None or len(sources) < 1:
            continue
        for source in sources:
            id, x, y = source['id'], source['xcentroid'], source['ycentroid']
            found = False
            for obj, coords in objects.items():
                if any(euclidean_distance(x, y, x2, y2) < closeness_threshold for _, _, x2, y2 in coords):
                    objects[obj].append((fname, id, x, y))
                    found = True
                    break
            if not found:
                objects[object_id].append((fname, id, x, y))
                object_id += 1

    # Filter out objects only detected in one image
    moving_objects = {obj: coords for obj, coords in objects.items() if len(set(fname for fname, _, _, _ in coords)) > 1}

    # Filter out objects that don't move, i.e. largest distance is less then 0.05
    for obj, coords in moving_objects.items():
        if max(euclidean_distance(x1, y1, x2, y2) for _, _, x1, y1 in coords for _, _, x2, y2 in coords) < 0.05:
            del moving_objects[obj]

    # Filter SOURCES to only include moving objects
    for fname in SOURCES:
        if SOURCES[fname] is None or len(SOURCES[fname]) < 1:
            continue
        SOURCES[fname] = [source for source in SOURCES[fname] if any(source['id'] == id for obj, coords in moving_objects.items() for _, id, _, _ in coords)]

    return SOURCES


def remove_lines(SOURCES, line_threshold=3, noise_threshold=5):
    """
    Filter out sources detected within noisy horizontal lines in the image.
    
    Args:
    - sources: Table of detected sources, as returned by detect_objects.
    - line_threshold: The thickness threshold above which a horizontal region is considered a line.
    - noise_threshold: The intensity threshold to detect noise lines; adjust based on your data.
    
    Returns:
    A filtered table of sources that are not within the detected noisy lines.
    """
    if SOURCES is None:
        return None
    
    for filename in SOURCES:
        if SOURCES[filename] is None or len(SOURCES[filename]) < 1:
            continue
    
        sources = SOURCES[filename]

        with fits.open(ALIGNED_PATH + filename) as hdul:
            data = hdul[0].data
    
        # Calculate the mean intensity of each row
        row_means = np.mean(data, axis=1)
        row_std = np.std(data, axis=1)
        
        # Identify rows that significantly deviate from the mean
        noise_rows = []
        for i, (mean, std) in enumerate(zip(row_means, row_std)):
            if abs(mean - np.mean(row_means)) > noise_threshold * np.mean(row_std):
                noise_rows.append(i)
        
        # Expand identified rows into regions based on line_threshold
        noise_regions = []
        start = None
        for i in noise_rows:
            if start is None:
                start = i
            elif i - start > line_threshold:
                noise_regions.append((start, i))
                start = i
            else:
                continue
        if start is not None:
            noise_regions.append((start, noise_rows[-1] + 1))  # Include the last line
        
        # Filter out sources within these regions
        filtered_sources = sources.copy()
        for start, end in noise_regions:
            filtered_sources = filtered_sources[~((filtered_sources['ycentroid'] >= start) & (filtered_sources['ycentroid'] <= end))]

        SOURCES[filename] = filtered_sources
    
    return SOURCES

"""
Parameters: list of filenames of fits files, the date of observation of the images.
Description: Runs the entire detection pipeline on the input images.
Outputs: Dictionary of filenames and their corresponding detected sources, including of the stacked image.
"""
def DETECTION_PIPELINE(fits_files, date_time_obs):  

    # Init dictionary to store sources
    SOURCES = {}
        
    # Alignment
    source_tuples = []
    base_source = None
    for filename in fits_files:
        cur_source = detect_objects(INPUT_PATH + filename, threshold=15, border=True) # Too many points = inaccurate alignment
        # save detected objects image
        if not os.path.exists("detected/"):
            os.makedirs("detected/")
        save_flagged_image(INPUT_PATH + filename, "detected/" + filename.split(".")[0] + "_d.png", [(source['xcentroid'], source['ycentroid']) for source in cur_source])
        source_tuples.append((cur_source, filename))
    
    # Select the source with the median detected objects as the base for alignment
    source_tuples.sort(key=lambda x: len(x[0]))
    base_source = source_tuples[len(source_tuples) // 2][0]
    # print(f"Using {source_tuples[len(source_tuples) // 2][1]} having {len(base_source)} detected sources as the base for alignment .")

    for source in source_tuples:
        shift = align(base_source, source[0], source[1])

    # Crop to same size
    crop_all(ALIGNED_PATH)
    
    # Create a stacked image
    stacked_image_path = STACKED_PATH + date_time_obs + ".fits"
    image_stacker(input_path=ALIGNED_PATH, output_path=stacked_image_path)
    stacked_image_sources = detect_objects(stacked_image_path)

    # Add stacked image sources to SOURCES
    SOURCES[stacked_image_path] = stacked_image_sources

    # Perform pixel differencing with all aligned images
    for filename in fits_files:
        error = pixel_difference_with_stack(os.path.join(ALIGNED_PATH, filename), stacked_image_path)
        if error == -1:
            continue
        sources = detect_objects(DIFFERENCED_PATH + filename, threshold=40, fwhm=24, sharplo=0.4)
        if sources is None or len(sources) < 1:
            SOURCES[filename] = None
            continue

        # Add sources to SOURCES
        SOURCES[filename] = sources

    return SOURCES


def FILTER_SOURCES(SOURCES, COSMIC_RAY_HITS):

    # Get sources from stacked images 
    stacked_sources = [SOURCES[filename] for filename in SOURCES if "stacked" in filename]

    # Iterate across dictionary
    for filename in SOURCES:
        if SOURCES[filename] is None or "stacked" in filename:
            continue
        
        sources = SOURCES[filename]
        for source in sources:
            x_coord = source['xcentroid']
            y_coord = source['ycentroid']
            rows_to_keep = []
            # Check if the source is in the stacked image, i.e. false positive because it has not moved.
            for stacked_source in stacked_sources:
                if stacked_source is not None and is_duplicate(x_coord, y_coord, stacked_source):
                    print(f"An object flagged for classification is a false positive in {filename}. Removing from unclassified objects.")
                else:
                    rows_to_keep.append(source)

            # Check if source is inside any detected streak
            if filename in COSMIC_RAY_HITS:
                if is_source_in_streak(x_coord, y_coord, COSMIC_RAY_HITS[filename]):
                    print(f"An object flagged for classification is actually a cosmic ray hit in {filename}. Removing from unclassified objects.")
                else:
                    rows_to_keep.append(source)

            SOURCES[filename] = rows_to_keep

    SOURCES = find_moving_objects(SOURCES)
    SOURCES = remove_lines(SOURCES)

    # Save flagged images
    for filename in SOURCES:
        if "stacked" in filename:
            continue
        sources = SOURCES[filename]
        if sources is not None:
            sources_coordinates = [(source['xcentroid'], source['ycentroid']) for source in sources]
            save_flagged_image(INPUT_PATH + filename, FILTERED_FLAGGED_ORIGINAL_PATH + filename.split(".")[0] + ".png", sources_coordinates)

    return SOURCES


def UPDATE_HEADER(SOURCES):

    if SOURCES is None or len(SOURCES) == 0:
        print("No sources detected. Exiting.")
        return -1
    
    for filename in SOURCES:
        if "stacked" in filename:
            continue

        sources = SOURCES[filename]
        if sources is None or len(sources) == 0:
            continue

        with fits.open(INPUT_PATH + filename, mode='update') as hdul:
            header = hdul[0].header
        
            # Add detected source information to the header for each source
            for i, source in enumerate(sources):

                idx = source['id']
                header[f'ID_{idx}'] = (source['id'], f'Identifier of detected source {i}')
                header[f'XCENT_{idx}'] = (source['xcentroid'], f'X centroid of detected source {i}')
                header[f'YCENT_{idx}'] = (source['ycentroid'], f'Y centroid of detected source {i}')
                header[f'MAG_{idx}'] = (source['mag'], f'Magnitude of detected source {i}')

                # header[f'SHARP_{i}'] = (source['sharpness'], f'Sharpness of detected source {i}')
                # header[f'RND1_{i}'] = (source['roundness1'], f'Roundness1 of detected source {i}')
                # header[f'RND2_{i}'] = (source['roundness2'], f'Roundness2 of detected source {i}')
                # header[f'NPIX_{i}'] = (source['npix'], f'Number of pixels of detected source {i}')
                # header[f'SKY_{i}'] = (source['sky'], f'Sky level near detected source {i}')
                # header[f'PEAK_{i}'] = (source['peak'], f'Peak value of detected source {i}')
                # header[f'FLUX_{i}'] = (source['flux'], f'Total flux of detected source {i}')
        
            # Save changes to header
            hdul.flush()  # Writes the updated header to the file
        
        # Now make a copy of the updated file to add the flagged folder
        shutil.copy(INPUT_PATH + filename, FLAGGED_ORIGINAL_FITS_PATH)
        print(f"Flagged {filename} for classification.")
        return 0


if __name__ == "__main__":
    
    if not os.path.exists(INPUT_PATH) or len(os.listdir(INPUT_PATH)) == 0:
        print("Input path does not exist or is empty. Exiting.")
        exit()
    
    # Unzip any files in the input path that end in .gz (cleaner compresses)
    for filename in os.listdir(INPUT_PATH):
        if filename.endswith('.gz'):
            os.system(f"gunzip {INPUT_PATH + filename}")

    # Get all fits files in the input path
    INPUT_FILES = [filename for filename in os.listdir(INPUT_PATH) if filename.lower().endswith('.fits')]

    if not os.path.exists(ALIGNED_PATH):
        os.makedirs(ALIGNED_PATH)
    if not os.path.exists(DIFFERENCED_PATH):
        os.makedirs(DIFFERENCED_PATH) 
    if not os.path.exists(FLAGGED_PATH):
        os.makedirs(FLAGGED_PATH) 
    if not os.path.exists(STACKED_PATH):
        os.makedirs(STACKED_PATH)
    if not os.path.exists(FLAGGED_ORIGINAL_PATH):
        os.makedirs(FLAGGED_ORIGINAL_PATH)
    if not os.path.exists(CRH_PATH):
        os.makedirs(CRH_PATH)
    if not os.path.exists(FILTERED_FLAGGED_ORIGINAL_PATH):
        os.makedirs(FILTERED_FLAGGED_ORIGINAL_PATH)
    if not os.path.exists(FLAGGED_ORIGINAL_FITS_PATH):
        os.makedirs(FLAGGED_ORIGINAL_FITS_PATH)

    # Separate files in input path by date, hour
    fits_files_by_time = {}
    for fits_file in INPUT_FILES:
        with fits.open(INPUT_PATH + fits_file) as hdul:
            header = hdul[0].header
            date_time_obs = header['DATE-OBS']
            date_obs = date_time_obs.split("T")[0]
            hour_obs = date_time_obs.split("T")[1][:2]
            key = date_obs + "_" + hour_obs
            if key not in fits_files_by_time:
                fits_files_by_time[key] = []
            fits_files_by_time[key].append(fits_file)
        
    # Run streak detection on each file
    COSMIC_RAY_HITS = {}
    for filename in INPUT_FILES:
        streaks = streak_detection(INPUT_PATH + filename) # Works best on raw images
        if streaks is not None:
            COSMIC_RAY_HITS[filename] = streaks    

    SOURCES = {}
    # Run object detection pipeline on each date
    for date_time_obs in fits_files_by_time:
        if len(fits_files_by_time[date_time_obs]) < 3:
            continue
        print(f"Running pipeline for {date_time_obs} on {len(fits_files_by_time[date_time_obs])} images.")
        SOURCES.update(DETECTION_PIPELINE(fits_files_by_time[date_time_obs], date_time_obs))
        print(f"Finished pipeline for {date_time_obs}\n")

    # Filter sources
    FILTERED_SOURCES = FILTER_SOURCES(SOURCES, COSMIC_RAY_HITS)
    
    # Update headers
    UPDATE_HEADER(SOURCES)

    # Save flagged images
    if FILTERED_SOURCES is not None and len(FILTERED_SOURCES) > 0:
        #print(f"Flagged {filename} for classification.")    
        sources_coordinates = [(source['xcentroid'], source['ycentroid']) for source in FILTERED_SOURCES]
        save_flagged_image(DIFFERENCED_PATH + filename, FLAGGED_PATH + filename.split(".")[0] + ".png", sources_coordinates)
        save_flagged_image(INPUT_PATH + filename, FLAGGED_ORIGINAL_PATH + filename.split(".")[0] + ".png", sources_coordinates)


    print(f"Finished object detection pipeline on {INPUT_PATH}.")