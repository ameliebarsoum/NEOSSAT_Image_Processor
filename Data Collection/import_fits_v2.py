#   -- NOTES --
#
#   meant to run on google colab, the script filters the fits images then downloads them
#   to google drive folder called "/content/drive/MyDrive/NEOSsat Images" and places them into 
#   folders, either "object name" or "dark" => so that we can then feed it to jason's cleaning
#   script for dark images.
#
#   left to do...
#   - understand target coordinates
#   - not sure if it's filtering the date properly, check that
#   - lots of duplicate images aren't getting downloaded, check if 'observationID' is what we should be 
#   labeling the fits files
#
#   how to run...
#
#   1. new cell for downloads, run cell
#
#   !pip install astroquery
#   !pip install astropy
#
#   2. new cell with code bellow, run cell (will ask you to connect drive to colab, accept)

import os
from datetime import datetime
from astroquery.cadc import Cadc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy as np
from google.colab import drive
import requests

drive.mount('/content/drive')

# Initialize the CADC client
cadc = Cadc()

# Define the target coordinates and search radius (i dont know need to change)
coords = SkyCoord(240, -30, unit='deg')  # RA, DEC, Unit
radius = 2 * u.deg

# Get today's date in the required format
today_date = datetime.now().strftime('%Y-%m-%d')

# Query CADC for data within the specified region, collection, and date
results = cadc.query_region(coords, radius, collection='NEOSSAT')

# Filter results based on 'time_exposure' and 'instrument_keywords'
accepted_modes = ['16-FINE_POINT', '14-FINE_SLEW']

# Filter results to only include rows where the 'instrument_keywords' column contains either '16-FINE_POINT' or '14-FINE_SLEW'
mask_modes = np.array([keyword in accepted_modes for keyword in results['instrument_keywords']])

# And combine it with your exposure time condition
mask_exposure = results['time_exposure'] > 50.0

# Apply both masks to filter the results
filtered_results = results[mask_modes & mask_exposure]

# Get a list of image URLs based on the filtered results
image_list = cadc.get_image_list(filtered_results, coords, radius)

# Download directory
download_directory = '/content/drive/MyDrive/NEOSsat Images'
if not os.path.exists(download_directory):
    os.makedirs(download_directory)

# Download FITS files and categorize them
for idx, url in enumerate(image_list):

    # Categorize the file based on 'target_name' and move to 'DarkImages' if 'type' is 'dark'
    target_name = filtered_results['target_name'][idx]
    image_type = filtered_results['type'][idx]

    # Check if the image is a dark image then add to dark directory
    if image_type == 'dark':
        target_directory = os.path.join(download_directory, 'DarkImages')

    # otherwise add to corresponding object folder
    else:
        target_directory = os.path.join(download_directory, target_name)

    # make directory if directory doesn't exist
    if not os.path.exists(target_directory):
        os.makedirs(target_directory)

    # name the files based off their observationID
    observation_id = filtered_results['observationID'][idx]  
    filename = str(observation_id) + ".fits"
    file_path = os.path.join(target_directory, filename)

    # download files to their corresponding directories
    if not os.path.exists(file_path):

        # Download the file using HTTP GET request
        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with open(file_path, 'wb') as f:
                f.write(response.content)
            print("Downloaded: " + filename)
        else:
            print("Failed to download " + filename + ": HTTP Status Code " + response.status_code)
    else:
        print("File already exists: " + filename)