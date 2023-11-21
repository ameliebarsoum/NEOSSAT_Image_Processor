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


INPUT_PATH = 'mission_images/'
# Just working with one image for now
OUTPUT_PATH = 'outgoing/'


def align(base_source, source2, filename):

    matched_indices = []
    tolerance = 1.0  # Tolerance for matching detected objects

    for i, star1 in enumerate(base_source):
        min_dist = np.inf
        min_idx = None
        for j, star2 in enumerate(source2):
            dist = distance.euclidean((star1['xcentroid'], star1['ycentroid']), (star2['xcentroid'], star2['ycentroid']))
            if dist < min_dist and j not in matched_indices:
                min_dist = dist
                min_idx = j
        if min_dist < tolerance and min_idx is not None:
            matched_indices.append(min_idx)

    # Get coordinates of matched stars in both images
    matched_stars_img1 = base_source
    matched_stars_img2 = [source2[i] for i in matched_indices]

    # Calculate the translation needed to align the stars in both images
    x_shift = np.mean([star2['xcentroid'] - star1['xcentroid'] for star1, star2 in zip(matched_stars_img1, matched_stars_img2)])
    y_shift = np.mean([star2['ycentroid'] - star1['ycentroid'] for star1, star2 in zip(matched_stars_img1, matched_stars_img2)])

    with fits.open(INPUT_PATH + filename) as hdul:
        img2 = hdul[0].data

    # Apply the calculated shift to align the second image with the first
    aligned_img2 = shift(img2, (y_shift, x_shift))
    fits.writeto(OUTPUT_PATH + 'shifted_' + filename, aligned_img2, overwrite=True)
    
    return True


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

if __name__ == "__main__":
    if not os.path.exists(OUTPUT_PATH):
        os.makedirs(OUTPUT_PATH)

    base_file = None
    input_files = []

    # Collect all FITS files in the INPUT_PATH directory
    for filename in os.listdir(INPUT_PATH):
        # Check if the filename ends with '.fits' (case insensitive)
        if filename.lower().endswith('.fits'):
            if base_file is None:
                base_file = filename
            else:
                input_files.append(filename)

    if base_file is not None:
        base_source = detect_objects(os.path.join(INPUT_PATH, base_file))

        for filename in input_files:
            source2 = detect_objects(os.path.join(INPUT_PATH, filename))
            align(base_source, source2, filename)
        if len(input_files) == 0:
             print("Only one FITS file found; not enough to find anomalies.")
            
    else:
        print("No FITS files found in the specified directory.")



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