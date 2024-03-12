import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
import numpy as np
import cv2
from scipy.ndimage import shift
from skimage.transform import resize
import random
from astropy.nddata import CCDData
from astropy.stats import mad_std
from ccdproc import CCDData, trim_image, combine
import warnings 
warnings.filterwarnings("ignore")

INPUT_PATH = 'comet_leonard/'
ALIGNED_PATH = 'aligned/'
GAUSSIANBLUR_PATH = 'gaussian_blur/'
DIFFERENCED_PATH = 'differenced/'
FLAGGED_PATH = 'flagged/'
STACKED_PATH = 'stacked/stacked_img.fits'

"""
Input: Path to fits file
Description: Detects objects in the image using the DAOStarFinder algorithm
"""
def detect_objects(filepath):
    # Open the FITS file
    with fits.open(filepath) as hdul:
        data = hdul[0].data  # Accessing the data from the primary HDU

    # Calculate background statistics
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)

    # DAOStarFinder algorithm for star detection
    daofind = DAOStarFinder(fwhm=10.0, threshold=30.*std)
    sources = daofind(data - median)

    return sources

"""
Helper function to find the closest star to a reference star
"""
def find_closest_star(ref_coords, candidates):
    min_distance = np.inf
    closest_index = -1
    for i, candidate in enumerate(candidates):
        distance = np.sqrt((candidate['xcentroid'] - ref_coords[0]) ** 2 + (candidate['ycentroid'] - ref_coords[1]) ** 2)
        if distance < min_distance:
            min_distance = distance
            closest_index = i
    if min_distance > 10:
        return -1
    return closest_index, min_distance

"""
@base_sources: sources from an image with which the other image will be aligned to
@img_sources: sources from the image to be aligned
@img_filename: filename of the image to be aligned
description: Aligns one image to another based on the brightest star, such that their stars overlap.
"""
def align(base_sources, img_sources, img_filename):

    # If no matching stars are found, try to find the closest star to the brightest star in the base image
    sorted_indices = np.argsort(base_sources['flux'])[::-1]
    # Ignore brightest, it is the most likely of being a cosmic ray/comet/etc which moves.
    sorted_indices = sorted_indices[1:]
    for idx in sorted_indices:
        img_star_idx, distance = find_closest_star((base_sources[idx]['xcentroid'], base_sources[idx]['ycentroid']), img_sources)
        if img_star_idx != -1 and distance < 0.1:
            base_star_idx = idx
            break
    
    if img_star_idx == -1:
        print(f"Failed to find a matching star for alignment. Not moving forward with {img_filename}.")
        return -1

    # Get the coordinates of the selected star in both images
    ref_star_coords_image1 = (base_sources['xcentroid'][base_star_idx], base_sources['ycentroid'][base_star_idx])
    ref_star_coords_image2 = (img_sources['xcentroid'][img_star_idx], img_sources['ycentroid'][img_star_idx])

    # Calculate the shift needed to align the images
    shift_x = ref_star_coords_image1[0] - ref_star_coords_image2[0]
    shift_y = ref_star_coords_image1[1] - ref_star_coords_image2[1]

    # Shift one image relative to the other using interpolation
    with fits.open(GAUSSIANBLUR_PATH + img_filename, mode='update') as hdul:
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
        blurred_image = cv2.GaussianBlur(image.astype(np.float32), (15, 15), 4)
        
        # Combine the original image with the darkened blurred image using the masks
        final_image = image * inverted_mask + blurred_image * mask

        # Save the background-subtracted and darkened image to a new FITS file
        output_filename = GAUSSIANBLUR_PATH + filename
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
                             output_file=STACKED_PATH,
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
Perform pixel differencing on two FITS files
"""
def pixel_difference(image1_filename, image2_filename, idx):

    image1_file = os.path.join(ALIGNED_PATH, image1_filename)
    image2_file = os.path.join(ALIGNED_PATH, image2_filename)

    # Load the aligned FITS files
    try:
        with fits.open(image1_file) as hdul1, fits.open(image2_file) as hdul2:
            data1 = hdul1[0].data
            data2 = hdul2[0].data

            # Resize images to a common size using interpolation
            # Images are not all the same size. Not sure if this is a result of the cleaning algorithm or not.
            min_shape = (min(data1.shape[0], data2.shape[0]), min(data1.shape[1], data2.shape[1]))
            data1_resized = resize(data1, min_shape, mode='constant')
            data2_resized = resize(data2, min_shape, mode='constant')

            # Perform pixel-wise difference
            diff_data = np.abs(data1_resized - data2_resized)  # Calculate absolute difference

            # Save the difference image to a new FITS file
            output_filename = DIFFERENCED_PATH + str(idx) + image2_filename
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
Call this to see a plot of the detected objects
"""
def detect_objects_with_visualizer(filename):
    # Open the FITS file
    with fits.open(INPUT_PATH + filename) as hdul:
        data = hdul[0].data  # Accessing the data from the primary HDU

    # Calculate background statistics
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)

    # DAOStarFinder algorithm for star detection
    daofind = DAOStarFinder(fwhm=10.0, threshold=30.*std)
    sources = daofind(data - median)

    # Display the image
    # First plot: FITS Image
    plt.figure(figsize=(12, 6))

    # Create subplot 1 (1 row, 2 columns, first plot)
    plt.subplot(1, 2, 1)
    plt.imshow(data, cmap='viridis', origin='lower')
    plt.colorbar(label='Pixel Value')
    plt.title('FITS Image')
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

    # Adjust layout for better appearance
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
    if not os.path.exists(GAUSSIANBLUR_PATH):
        os.makedirs(GAUSSIANBLUR_PATH)
    if not os.path.exists(DIFFERENCED_PATH):
        os.makedirs(DIFFERENCED_PATH) 
    if not os.path.exists(FLAGGED_PATH):
        os.makedirs(FLAGGED_PATH) 
    if not os.path.exists("stacked"):
        os.makedirs("stacked")

    # Alignment
    base_source = None
    for filename in fits_files:
        if base_source is None:
            base_source = detect_objects(INPUT_PATH + filename)
            blur_background(filename, base_source)
        else:
            cur_source = detect_objects(INPUT_PATH + filename)
            blur_background(filename, cur_source)
            align(base_source, cur_source, filename)

    if base_source is None:
        print("No FITS files found in the specified directory.")
        exit()

    crop_all(ALIGNED_PATH)
    
    # Pixel Differencing
    if len(fits_files) < 4:
        fits_files = [filename for filename in os.listdir(INPUT_PATH) if filename.lower().endswith('.fits')]
        # Iterate through each file and compare it to every other file
        for i, base_filename in enumerate(fits_files):
            for j, comparison_filename in enumerate(fits_files):
                if i != j:  # Avoid comparing a file to itself
                    base_path = os.path.join(INPUT_PATH, base_filename)
                    comparison_path = os.path.join(INPUT_PATH, comparison_filename)
                    pixel_difference(base_path, comparison_path, j)
    else: 
        stacked_image = image_stacker(ALIGNED_PATH=ALIGNED_PATH)

        # Perform pixel differencing with all aligned images
        for filename in os.listdir(ALIGNED_PATH):
            if filename.lower().endswith('.fits'):
                pixel_difference_with_stack(filename)


    # Object detection on difference images
    for filename in os.listdir(DIFFERENCED_PATH):
        if filename.lower().endswith('.fits'):
            sources = detect_objects(DIFFERENCED_PATH + filename)
            # Move files with stars detected into flagged folder
            if sources is None or len(sources) != 0:
                os.rename(DIFFERENCED_PATH + filename, FLAGGED_PATH + filename)
                print(f"Flagged {filename} for classification.")

