## 
## Credit to https://github.com/asc-csa/NEOSSAT_Tutorial/blob/main/Code/Notebook%201_%20Extracting%20Data%20and%20Visualization.ipynb
##

import os
from datetime import datetime
from astroquery.cadc import Cadc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy as np
import requests

def fetch_images(RA=240, DEC=0, radius=None):
    print("Fetching images at ", RA, " ", DEC)
     
    # Initialize the CADC client
    cadc = Cadc()

    # Define the target coordinates and search radius (i dont know need to change)
    coords = SkyCoord(RA, DEC, unit='deg')#SkyCoord(RA, DEC, unit='deg')  # RA, DEC, Unit
    if not radius:
        radius = 2 * u.deg
    else:
        radius = radius * u.deg

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
    download_directory = './FITSImages_' + str(RA) + "_" + str(DEC)
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
        filename = "NEOS_SCI_" + str(observation_id) + ".fits"  # Image cleaner expects this format
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

if __name__ == "__main__":
    RA = input("Enter the right ascension (defaults to 240): ")
    DEC = input("Enter the declination (defaults to 0): ")
    radius = input("Enter the radius range (defaults to 2 degrees): ")
    
    fetch_images(RA, DEC, radius)