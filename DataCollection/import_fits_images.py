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
import sys

def main(args):

    if len(args) < 3 or len(args) > 4: 
        print(usage())
        sys.exit(1) 

    RA = args[1]
    DEC = args[2] if len(args) > 2 else None
    radius = args[3] if len(args) > 3 else None

    print("Fetching images at", RA, DEC, "along a radius of", radius)
     
    # Initialize the CADC client
    cadc = Cadc()

    # Define the target coordinates and search radius (i dont know need to change)
    coords = SkyCoord(RA, DEC, unit='deg')#SkyCoord(RA, DEC, unit='deg')  # RA, DEC, Unit
    if not radius:
        radius = 1 * u.deg
    else:
        radius = float(radius) * u.deg

    # Get today's date in the required format
    today_date = datetime.now().strftime('%Y-%m-%d')

    # Query CADC for data within the specified region, collection, and date
    results = cadc.query_region(coordinates=coords, radius=radius, collection='NEOSSAT')
    print("Found", len(results), "results at " + str(coords))

    # Filter results based on 'time_exposure' and 'instrument_keywords'
    accepted_modes = ['16-FINE_POINT', '14-FINE_SLEW']

    # Filter results to only include rows where the 'instrument_keywords' column contains either '16-FINE_POINT' or '14-FINE_SLEW'
    mask_modes = np.array([keyword in accepted_modes for keyword in results['instrument_keywords']])

    # And combine it with your exposure time condition
    mask_exposure = results['time_exposure'] > 50.0

    # Apply both masks to filter the results
    filtered_results = results[mask_modes & mask_exposure]
    print("Filtered to", len(filtered_results), "images")
    
    if len(filtered_results) == 0:
        print("No images found")
        return

    # Get a list of image URLs based on the filtered results
    image_list = cadc.get_image_list(query_result=filtered_results, coordinates=coords, radius=radius)

    # Download directory - assuming that this is being executed from the root directory (i.e., by run_all.sh)
    download_directory = 'DataCollection/FITSImages_' + str(RA) + "_" + str(DEC)
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
                print("Downloaded " + str(idx) + '/' + str(len(filtered_results)) + ": " + filename)
            else:
                print("Failed to download " + filename + ": HTTP Status Code " + response.status_code)
        else:
            print("File already exists: " + filename)

def usage():
    return """
Usage: ./import_fits_images.py <right ascension> <declination> <radius>
Description: This script that will fetch images from the CADC to perform object detection on.
            Right ascension and declination are celestial coordinates, like longitude and latitude. Both are expressed in degrees.
            Radius represents the range of RA/DEC that is being considered. It is also expresssed in degrees. A float between 0.2 and 1 is recommended.
"""

if __name__ == "__main__":
    main(sys.argv)
