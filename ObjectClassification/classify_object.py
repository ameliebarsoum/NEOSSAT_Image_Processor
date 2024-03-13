from astropy.io import fits
from astropy.wcs import WCS
import requests
from astropy.time import Time
from datetime import datetime

# default date value is the current date 
time_image_taken = datetime.now().strftime('%Y-%m-%d')

def get_centroid_values_from_hdu(filepath):
    # Open the FITS file and extract header information
    with fits.open(filepath) as hdul:
        header = hdul[0].header
        wcs = WCS(header)
        time_image_taken = header['DATE-OBS']

        # Initialize a dictionary to hold RA and Dec with source IDs
        source_coordinates = {}

        # Extract all XCENT and YCENT values and their IDs
        for key in header.keys():
            if key.startswith('XCENT'):
                id = key.split('_')[1]  # Assumes key format is "XCENT_{id}"
                xcent = header[key]
                ycent = header.get(f'YCENT_{id}', None)  # Attempt to get the corresponding YCENT value
                if ycent is not None:  # Ensure there is a corresponding YCENT value

                    # Convert pixel coordinates to world coordinates (RA, Dec)
                    ra, dec = wcs.all_pix2world(xcent, ycent, 0)
                    # Save the RA and Dec with the source ID
                    source_coordinates[id] = {'RA': ra, 'Dec': dec}

        # Print the stored coordinates
        for id, coords in source_coordinates.items():
            print(f"Source ID: {id}, RA: {coords['RA']}, Dec: {coords['Dec']}")

def query_horizons_for_object(ra, dec, time_image_taken):
    """
    Query the JPL Horizons system for objects near a given RA and Dec on a specific date.
    
    Parameters:
    - ra (float): Right Ascension in degrees.
    - dec (float): Declination in degrees.
    - date (str): Date the image was taken, in 'YYYY-MM-DD' format.
    
    Returns:
    - JSON response from the Horizons system with information about any matched objects.
    """
    
    # Convert the date to Julian Date, which is often used in astronomical databases
    jd = Time(time_image_taken).jd
    
    # Construct the query URL
    # Note: This URL structure is conceptual. You will need to refer to the actual Horizons API documentation
    # to construct a proper query. The following URL and parameters are illustrative only.
    url = f"https://ssd-api.jpl.nasa.gov/horizons.api?COMMAND='DES=;';CENTER='500@10';MAKE_EPHEM='YES';TABLE_TYPE='OBSERVER';START_TIME='{jd}';STOP_TIME='{jd + 1}';STEP_SIZE='1 d';OUT_UNITS='AU-D';OBJ_DATA='YES';"
    
    # Make the HTTP GET request
    response = requests.get(url)
    
    # Check if the request was successful
    if response.status_code == 200:
        # Parse and return the JSON data
        return response.json()
    else:
        print("Failed to query the Horizons system:", response.status_code)
        return None