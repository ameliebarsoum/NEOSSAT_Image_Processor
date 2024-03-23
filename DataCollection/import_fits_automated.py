## 
## Credit to https://github.com/asc-csa/NEOSSAT_Tutorial/blob/main/Code/Notebook%201_%20Extracting%20Data%20and%20Visualization.ipynb
##

import os
from datetime import datetime, timedelta
from astropy.time import Time
from astroquery.cadc import Cadc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy as np
import requests
import pyvo
from urllib.parse import urlencode

def fetch_images():
     
    # Initialize the CADC client
    cadc = Cadc()
        
    # Calculate the time frame
    end_time = datetime.now()                      # Current time
    start_time = end_time - timedelta(hours=1)     # One hour ago

    ### FOR TESTING: FIXED TIME INTERVAL #################################

    # # Hard coded start and end times
    # start_time_str = "2021-11-12 23:00:00"
    # end_time_str = "2021-11-13 00:00:00"

    # # Convert string times to datetime objects
    # start_time = datetime.strptime(start_time_str, '%Y-%m-%d %H:%M:%S')
    # end_time = datetime.strptime(end_time_str, '%Y-%m-%d %H:%M:%S')

    #####################################################################

    # Convert to Astropy Time object for easy conversion to MJD
    start_time_astropy = Time(start_time)
    end_time_astropy = Time(end_time)

    # Convert the datetime objects to MJD for the ADQL query
    start_mjd = start_time_astropy.mjd
    end_mjd = end_time_astropy.mjd

    # Formulate the ADQL query
    adql_query = f"""
        SELECT *
        FROM caom2.Plane AS Plane
        JOIN caom2.Observation AS Observation ON Plane.obsID = Observation.obsID 
        WHERE (Observation.collection = 'NEOSSAT'
        AND INTERSECTS(INTERVAL({start_mjd}, {end_mjd}), Plane.time_bounds_samples) = 1
        AND (Plane.quality_flag IS NULL OR Plane.quality_flag != 'junk'))
        AND Plane.time_exposure > 50.0
        AND (Observation.instrument_keywords LIKE '%16-FINE_POINT%' OR Observation.instrument_keywords LIKE '%14-FINE_SLEW%')
        """

    # Execute the query
    results = cadc.exec_sync(adql_query)
    if results is None or len(results) < 1:
        print("No images found")
        exit()

    # Editted function to ignore Coordinates and radius (we want all angles...)
    # Source: https://astroquery.readthedocs.io/en/latest/_modules/astroquery/cadc/core.html
    def get_image_list_modified(cadc, query_result):
        """
        Function to map the results of a CADC query into URLs to
        corresponding data that can be later downloaded, without applying
        cutout operations based on coordinates and radius.

        Parameters
        ----------
        query_result : A `~astropy.table.Table` object
            Result returned by `query_region` or `query_name`, or any
            CADC TAP query that contains the 'publisherID' column.

        Returns
        -------
        list : A list of URLs to data.
        """
        
        if not query_result:
            raise AttributeError('Missing query_result argument')

        try:
            publisher_ids = query_result['publisherID']
        except KeyError:
            raise AttributeError('publisherID column missing from query_result argument')

        result = []

        # Send datalink requests in batches of 20 publisher ids
        batch_size = 20

        # Iterate through list of sublists to send datalink requests in batches
        for pid_sublist in (publisher_ids[pos:pos + batch_size] for pos in
                            range(0, len(publisher_ids), batch_size)):
            datalink = pyvo.dal.adhoc.DatalinkResults.from_result_url(
                '{}?{}'.format(cadc.data_link_url,
                            urlencode({'ID': pid_sublist}, True)),
                session=cadc.cadcdatalink._session)
            
            # Fetch all data URLs without filtering for cutouts
            for service_def in datalink:
                if service_def.semantics == '#this':
                    result.append(service_def.access_url)

        return result
    
    # Get a list of image URLs based on the filtered results
    image_list = get_image_list_modified(cadc, results)

    # Formulate the ADQL query
    CrossRef_query = f"""
        SELECT 
        Observation.observationID AS "Obs. ID",
	    COORD1(CENTROID(Plane.position_bounds)) AS "RA (J2000.0)",
	    COORD2(CENTROID(Plane.position_bounds)) AS "Dec. (J2000.0)"
        FROM caom2.Plane AS Plane
        JOIN caom2.Observation AS Observation ON Plane.obsID = Observation.obsID 
        WHERE (Observation.collection = 'NEOSSAT'
        AND INTERSECTS(INTERVAL({start_mjd}, {end_mjd}), Plane.time_bounds_samples) = 1
        AND (Plane.quality_flag IS NULL OR Plane.quality_flag != 'junk'))
        AND Plane.time_exposure > 50.0
        AND (Observation.instrument_keywords LIKE '%16-FINE_POINT%' OR Observation.instrument_keywords LIKE '%14-FINE_SLEW%')
        """

    # List to cross referance the RA and DEC for folder name
    RA_DEC_List = cadc.exec_sync(CrossRef_query)

    # Create a mapping from Observation ID to RA and DEC
    obs_id_to_radec = {
        row['"Obs. ID"']: (row['"RA (J2000.0)"'], row['"Dec. (J2000.0)"'])
        for row in RA_DEC_List
    }

    # Iterate through each image URL in image_list
    for idx, url in enumerate(image_list):
        observation_id = results['observationID'][idx]  # Assuming this is how you get each observation ID

        # Find corresponding RA and DEC for the observation ID
        if observation_id in obs_id_to_radec:
            RA, DEC = obs_id_to_radec[observation_id]
        else:
            print(f"RA and DEC for observation ID {observation_id} not found.")
            continue

        # Format RA and DEC for file naming (rounded to no decimal places)
        RA_formatted = f"{int(round(RA))}".replace('.', 'p')
        DEC_formatted = f"{int(round(DEC))}".replace('.', 'p').replace('-', 'm')

        # Download directory based on RA and DEC
        download_directory = f'./FITSImages_{RA_formatted}_{DEC_formatted}'
        if not os.path.exists(download_directory):
            os.makedirs(download_directory)

        target_name = results['target_name'][idx]
        image_type = results['type'][idx]

        # Categorize the file based on 'target_name' and move to 'DarkImages' if 'type' is 'dark'
        if image_type == 'dark':
            target_directory = os.path.join(download_directory, 'DarkImages')
        else:
            target_directory = os.path.join(download_directory, target_name)

        # Make directory if it doesn't exist
        if not os.path.exists(target_directory):
            os.makedirs(target_directory)

        # Name the files based off their observationID
        filename = f"NEOS_SCI_{observation_id}.fits"
        file_path = os.path.join(target_directory, filename)

        # Download files to their corresponding directories
        if not os.path.exists(file_path):

            # Download the file using HTTP GET request
            response = requests.get(url, stream=True)
            if response.status_code == 200:
                with open(file_path, 'wb') as f:
                    f.write(response.content)
                print(f"Downloaded {idx + 1}/{len(image_list)}: {filename}")
            else:
                print(f"Failed to download {filename}: HTTP Status Code {response.status_code}")
        else:
            print(f"File already exists: {filename}")

if __name__ == "__main__":
    fetch_images()