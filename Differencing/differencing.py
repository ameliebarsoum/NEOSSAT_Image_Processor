from astropy.io import fits
import numpy as np
import os
from skimage.transform import resize

INPUT_PATH = '../Alignment/outgoing/'
OUTPUT_PATH = 'outgoing/'


def pixel_difference(image1_filename, image2_filename):
    # Load the aligned FITS files
    with fits.open(image1_filename) as hdul1, fits.open(image2_filename) as hdul2:
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
        output_filename = f"{OUTPUT_PATH}diff_{os.path.basename(image2_filename)[:-5]}.fits"
        hdu = fits.PrimaryHDU(diff_data, header=hdul1[0].header)
        hdul_out = fits.HDUList([hdu])
        hdul_out.writeto(output_filename, overwrite=True)

        print(f"Difference image saved as {output_filename}.")


if __name__ == "__main__":
    if not os.path.exists(OUTPUT_PATH):
        os.makedirs(OUTPUT_PATH)

    base_filename = None
    for filename in os.listdir(INPUT_PATH):
        # Check if the filename ends with '.fits' (case insensitive)
        if filename.lower().endswith('.fits'):
            if base_filename is None:
                base_filename = os.path.join(INPUT_PATH, filename)
            else:
                comparison_filename = os.path.join(INPUT_PATH, filename)
                pixel_difference(base_filename, comparison_filename)

    if base_filename is None:
        print("No FITS files found in the specified directory.")
