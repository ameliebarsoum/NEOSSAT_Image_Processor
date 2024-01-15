import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
import numpy as np
import cv2
from scipy.ndimage import shift, median_filter
from scipy.ndimage import gaussian_filter
from photutils.utils import make_random_cmap
from photutils import detect_threshold, detect_sources
"""
Pre-processing algorithm that applies peak finding, background subtraction
with a masked gaussian blur filter + darkening, and alignment based on peaks.
Preparation for pixel differencing.
Takes in images from the same section of the sky that have gone through 
photometry cleaning algorithm.
"""

INPUT_PATH = 'mission_images/'
OUTPUT_PATH = 'outgoing/'


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
    with fits.open(INPUT_PATH + 'bg_processed_' + image2_filename, mode='update') as hdul:
        print(f"dimensions of image: {hdul[0].data.shape}")
        shifted_data = shift(hdul[0].data, (shift_y, shift_x), mode='nearest')
        hdul[0].data = shifted_data
        
        hdu = fits.PrimaryHDU(shifted_data, header=hdul[0].header)
        hdul_out = fits.HDUList([hdu])
        hdul_out.writeto(OUTPUT_PATH + 'aligned_' + image2_filename, overwrite=True)
        print(f"dimensions of shifted image: {hdul_out[0].data.shape}")
    print(f"Aligning images with a shift of ({shift_x}, {shift_y}) pixels.")


def detect_objects(filename):
    # Open the FITS file
    with fits.open(INPUT_PATH + filename) as hdul:
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
        
        # Darkening factor for the blurred regions
        darkening_factor = 0.5  # Adjust this value as needed

        # Darken the blurred areas before blending with the original image
        darkened_blurred = blurred_image * darkening_factor
        
        # Combine the original image with the darkened blurred image using the masks
        final_image = image * inverted_mask + darkened_blurred * mask

        # Save the background-subtracted and darkened image to a new FITS file
        output_filename = INPUT_PATH + 'bg_processed_' + filename
        hdu = fits.PrimaryHDU(final_image, header=hdul[0].header)
        hdul_out = fits.HDUList([hdu])
        hdul_out.writeto(output_filename, overwrite=True)

        print(f"Background-subtracted and darkened image saved as {output_filename}.")

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



if __name__ == "__main__":
    if not os.path.exists(OUTPUT_PATH):
        os.makedirs(OUTPUT_PATH)

    base_source = None
    for filename in os.listdir(INPUT_PATH):
        # Check if the filename ends with '.fits' (case insensitive)
        if filename.lower().endswith('.fits'):
            if 'processed' in filename:  # Check if 'processed' is in the filename
                os.remove(os.path.join(INPUT_PATH, filename))  # While testing - delete the file 
            else:
                if base_source is None:
                    base_source = detect_objects(filename)
                    subtract_background(filename, base_source)
                else:
                    cur_source = detect_objects(filename)
                    subtract_background(filename, cur_source)
                    align(base_source, cur_source, filename)

    if base_source is None:
        print("No FITS files found in the specified directory.")
        exit()
