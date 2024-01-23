import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
import numpy as np
import cv2
from skimage.transform import resize
import random

"""
Pre-processing algorithm that applies peak finding, background subtraction
with a masked gaussian blur filter + darkening, and alignment based on peaks.
Preparation for pixel differencing.
Takes in images from the same section of the sky that have gone through 
photometry cleaning algorithm.
"""

INPUT_PATH = 'mission_images/'
ALIGNED_PATH = 'aligned/'
GAUSSIANBLUR_PATH = 'gaussian_blur/'
DIFFERENCED_PATH = 'differenced/'
FLAGGED_PATH = 'flagged/'

def align(sources_image1, sources_image2, image2_filename):

    # Choose a star present in both images for alignment - assume brightest star is present in both
    reference_star_index = np.argmax(sources_image1['flux']) 

    # Get the coordinates of the selected star in both images
    ref_star_coords_image1 = (sources_image1['xcentroid'][reference_star_index], sources_image1['ycentroid'][reference_star_index])
    ref_star_coords_image2 = (sources_image2['xcentroid'][reference_star_index], sources_image2['ycentroid'][reference_star_index])

    # Calculate the shift needed to align the images
    shift_x = ref_star_coords_image1[0] - ref_star_coords_image2[0]
    shift_y = ref_star_coords_image1[1] - ref_star_coords_image2[1]

    # Shift one image relative to the other using interpolation
    with fits.open(GAUSSIANBLUR_PATH + image2_filename, mode='update') as hdul:
        shifted_data = shift(hdul[0].data, (shift_y, shift_x), mode='nearest')
        hdul[0].data = shifted_data
        
        hdu = fits.PrimaryHDU(shifted_data, header=hdul[0].header)
        hdul_out = fits.HDUList([hdu])
        hdul_out.writeto(ALIGNED_PATH + image2_filename, overwrite=True)
    #print(f"Aligning images with a shift of ({shift_x}, {shift_y}) pixels as {ALIGNED_PATH + image2_filename}.")


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
Subtract backgrounds of image located at @filename
based on peaks given by @sources
"""
def subtract_background(filename, sources):
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
            output_filename = DIFFERENCED_PATH + idx + image2_filename
            hdu = fits.PrimaryHDU(diff_data, header=hdul1[0].header)
            hdul_out = fits.HDUList([hdu])
            hdul_out.writeto(output_filename, overwrite=True)
    except FileNotFoundError:
        print(f"Could not find {image1_file} or {image2_file}. Previous steps may not have been completed.")
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
helper for main 
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

    # Alignment
    base_source = None
    for filename in fits_files:
        if base_source is None:
            base_source = detect_objects(INPUT_PATH + filename)
            subtract_background(filename, base_source)
        else:
            cur_source = detect_objects(INPUT_PATH + filename)
            subtract_background(filename, cur_source)
            align(base_source, cur_source, filename)

    if base_source is None:
        print("No FITS files found in the specified directory.")
        exit()

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
        for filename in os.listdir(INPUT_PATH):
            if filename.lower().endswith('.fits'):                
                # Get 3 random files for comparison
                random_files = get_random_files(INPUT_PATH, 3)
                
                for i, random_filename in enumerate(random_files):
                    pixel_difference(filename, random_filename, i)

    # Object detection on difference images
    for filename in os.listdir(DIFFERENCED_PATH):
        if filename.lower().endswith('.fits'):
            sources = detect_objects(DIFFERENCED_PATH + filename)
            # Move files with stars detected into flagged folder
            if sources is None or len(sources) != 0:
                os.rename(DIFFERENCED_PATH + filename, FLAGGED_PATH + filename)
                print(f"Flagged {filename} for classification.")

    