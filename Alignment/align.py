from astropy.convolution import interpolate_replace_nans
import os
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
import numpy as np
import matplotlib.pyplot as plt
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from scipy.spatial import distance
from scipy.ndimage import shift
from skimage.transform import AffineTransform, warp
import itertools


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
    with fits.open(INPUT_PATH + image2_filename, mode='update') as hdul:
        shifted_data = shift(hdul[0].data, (shift_y, shift_x), mode='nearest')
        hdul[0].data = shifted_data
        
        hdu = fits.PrimaryHDU(shifted_data, header=hdul[0].header)
        hdul_out = fits.HDUList([hdu])
        hdul_out.writeto(OUTPUT_PATH + 'aligned_' + filename, overwrite=True)

    print(f"Aligning images with a shift of ({shift_x}, {shift_y}) pixels.")

    # Update header information ?


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


## Call this to see a plot of the detected objects
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
            if base_source is None:
                base_source = detect_objects(filename)
            else:
                align(base_source, detect_objects(filename), filename)

    if base_source is None:
        print("No FITS files found in the specified directory.")
        exit()
