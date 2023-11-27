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


## Rectification with triangle alignment
## Probably unnecessary for our purposes
def get_triangle(objects):
    """
    Inefficient way of obtaining the lengths of each triangle's side.
    Normalized so that the minimum length is 1.
    """

    # 'stars' should be a list of detected stars' coordinates
    # For example, stars = [(x1, y1), (x2, y2), (x3, y3), ...]

    # Find all possible combinations of three stars
    possible_triangles = list(itertools.combinations(range(len(objects)), 3))

    for triangle in possible_triangles:
        # Extract coordinates of three points for each triangle combination
        p1, p2, p3 = objects[triangle[0]], objects[triangle[1]], objects[triangle[2]]

        # Calculate distances between points (you can use different criteria here)
        dist1 = ((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)**0.5
        dist2 = ((p3[0] - p2[0])**2 + (p3[1] - p2[1])**2)**0.5
        dist3 = ((p1[0] - p3[0])**2 + (p1[1] - p3[1])**2)**0.5

        # You can set conditions here to determine a valid triangle
        # For example, you might check for specific angle conditions or ratios between distances
        # If the conditions are met, return the indices of the points forming the triangle
        # Here's a placeholder condition for demonstration purposes
        if dist1 < dist2 + dist3 and dist2 < dist1 + dist3 and dist3 < dist1 + dist2:
            return triangle  # Return indices of points forming a valid triangle

    # Return None if no valid triangle is found
    return None


def main_for_rectification():
# if __name__ == "__main__":
    if not os.path.exists(OUTPUT_PATH):
        os.makedirs(OUTPUT_PATH)

    # detect_objects_with_visualizer('NEOS_SCI_2023138092727.fits')

    file_sources = {}  # Dictionary to store detected objects for each file
    base_source = None  # The file with the most detected objects will be used as the base for optimal alignment
    max_stars = -1

    # Collect all FITS files in the INPUT_PATH directory and detect objects
    for filename in os.listdir(INPUT_PATH):
        # Check if the filename ends with '.fits' (case insensitive)
        if filename.lower().endswith('.fits'):
            file_sources[filename] = detect_objects(filename)
            num_stars = len(file_sources[filename])

            if num_stars > max_stars:
                max_stars = num_stars
                base_source = file_sources[filename] 

    if base_source is None:
        print("No FITS files found in the specified directory.")
        exit()
    
    base_triangle = get_triangle(base_source)
    if base_triangle is None:
        print("No image with enough stars detected to serve as a base for triangle rectification.")
        ### WRITE code for simpler alignment ###
        exit()

    for filename, sources in file_sources.items():
        if sources is base_source:
            continue

        if len(sources) < 3:
            print("Not enough stars found in " + filename + " to align with triangle rectification.")
            ### Write code for simpler alignment ###
            continue

        triangle = get_triangle(sources)
        if triangle is None:
            align(base_source, sources, filename)

    if max_stars == -1:
        print("No FITS files found in the specified directory.")
    elif max_stars < 3:
        print("No image with enough stars detected to serve as a base for triangle rectification.")
