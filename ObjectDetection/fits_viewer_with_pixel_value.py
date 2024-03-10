from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# Replace 'your_file.fits' with the path to your FITS file
fits_file = 'mission_images/NEOS_SCI_2023138093527.fits'

# Open the FITS file
with fits.open(fits_file) as hdul:
    data = hdul[0].data  # Accessing the data from the primary HDU

# Display the image
plt.figure(figsize=(8, 6))
plt.imshow(data, cmap='viridis', origin='lower')
plt.colorbar(label='Pixel Value')
plt.title('FITS Image')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.show()

# Inspect pixel values interactively
def onclick(event):
    x, y = int(event.xdata), int(event.ydata)
    print(f"Pixel value at (x={x}, y={y}): {data[y, x]}")  # Accessing pixel value

fig, ax = plt.subplots()
ax.imshow(data, cmap='viridis', origin='lower')
plt.title('Click on the image to get pixel value')

cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()
