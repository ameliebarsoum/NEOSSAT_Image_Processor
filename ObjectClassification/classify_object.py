import os
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import subprocess
import warnings
from astropy.utils.exceptions import AstropyWarning

# Ignore specific FITSFixedWarnings related to 'datfix' and 'unitfix'
warnings.filterwarnings('ignore', category=AstropyWarning, message=".*'datfix' made the change.*")
warnings.filterwarnings('ignore', category=AstropyWarning, message=".*'unitfix' made the change.*")

INPUT_PATH = '../ObjectDetection/flagged_original_fits/'

# Function to call astcheck with the specified input file
def call_astcheck(input_file_path):
    # Path to the astcheck executable
    astcheck_executable_path = './astcheck'  # Adjust this as necessary
    
    # Build the command as a list
    command = [astcheck_executable_path, input_file_path]
    
    # Call astcheck and capture output
    result = subprocess.run(command, capture_output=True, text=True)
    print(result.stdout)
    
    # Count the lines in the output
    num_lines = len(result.stdout.split('\n'))
    
    # Check if there are more than 14 lines, more than 14 means motion exists
    if num_lines > 14:
        return True
    
    # If only 14 lines are emitied, no object found
    else:
        return False

def format_number_to_exact_4_chars(num):
    if num < 0:
        # For negative numbers, allow 1 digit before and 2 digits after the decimal
        formatted = f"{num:.2f}"
        if len(formatted) > 4:  # Truncate if more than 4 chars (rare due to our format)
            formatted = formatted[:4]
    else:
        # For positive numbers, we allow 1 digit before and 2 digits after the decimal
        formatted = f"{num:4.2f}"  # This ensures the total length is 4
        if len(formatted) > 4:  # Adjust if somehow exceeds 4 characters
            formatted = formatted[:4]
    return formatted

def format_mpc_observation(ra_deg, dec_deg, mag, obs_date, obs_code="C53"):
    """
    Formats observation data into the 80-column MPC format.
    
    Parameters:
    - ra_deg, dec_deg: Right Ascension and Declination in decimal degrees.
    - mag: Magnitude of the object.
    - obs_date: Observation date in the format YYYY MM DD.dddddd.
    - obs_code: Observatory code, changed to C53 for NEOSSat (https://www.minorplanetcenter.net/iau/lists/ObsCodes.html)
    
    Returns:
    - Formatted string in the 80-column MPC format.
    """
    # Convert RA and Dec from decimal degrees to sexagesimal format
    coord = SkyCoord(ra=ra_deg * u.degree, dec=dec_deg * u.degree, frame='icrs')
    ra_sexagesimal = coord.ra.to_string(unit=u.hour, sep=' ', pad=True, precision=2, alwayssign=False)
    dec_sexagesimal = coord.dec.to_string(unit=u.degree, sep=' ', pad=True, alwayssign=True, precision=1)

    # Convert obs_date to the required MPC format
    t = Time(obs_date)

    # Truncate to obtain YYYY MM DD.ddddd format
    obs_date_formatted = t.datetime.strftime('%Y %m %d.%f')[:16]  

    # Format the magnitude correctly, handling negative magnitudes
    if mag is not None:

        # Filters in the beam. Fixed to ‘CLEAR’ for NEOSSat
        mag_formatted = format_number_to_exact_4_chars(mag)
        mag_formatted = mag_formatted + " C"

    # Provide a placeholder of the correct width for when magnitude is not available
    else:
        # Six spaces to align correctly in the MPC format
        mag_formatted = "      "

    # Construct the 80-column observation line according to MPC guidelines
    formatted_observation = f"             "                                    # Placeholder for object designation, columns 1-13 are usually for the object identifier
    formatted_observation += f"  "                                              # Placeholder for additional notes, column 14-15
    formatted_observation += f"{obs_date_formatted} "                           # Observation date, occupies columns 16-32
    formatted_observation += f"{ra_sexagesimal} {dec_sexagesimal}"              # RA and Dec, right-justified within their field, column 33-56
    formatted_observation += f"          "                                      # Blank, for column 57-65
    formatted_observation += f"{mag_formatted}"                                 # Magnitude and band, for columns 66-71
    formatted_observation += f"      "                                          # Blank, for columns 72-77
    formatted_observation += f"{obs_code}"                                      # Observatory code, columns 78-80

    return formatted_observation

def findOrbit(observations_file_path):
    return None

if __name__ == "__main__":

    fits_files = [filename for filename in os.listdir(INPUT_PATH) if filename.lower().endswith('.fits')]

    for file in fits_files:

        # Hold all observations to find orbit
        observations_file_path = 'observations_temp.txt'

        # Open the FITS file
        with fits.open(INPUT_PATH + file, mode='update') as hdul:

            # Extract header data
            header = hdul[0].header
            wcs = WCS(header)

            # Date/time parameter
            time_image_taken = header['DATE-OBS']

            # Field to hold the ids of all the unclassified objects inside the fits header
            unclassified_ids = []

            # Extract all XCENT and YCENT values and their IDs
            for key in header.keys():

                # Assumes key format is "XCENT_{id}"
                if key.startswith('XCENT'):
                    id = key.split('_')[1] 
                    
                    # Example of extracting one source's details - this would be looped for multiple sources
                    xcentroid = header.get(f'XCENT_{id}', None)
                    ycentroid = header.get(f'YCENT_{id}', None)
                    mag = header.get(f'MAG_{id}', None)
            
                    # Convert pixel coordinates (xcentroid, ycentroid) to celestial coordinates (RA, Dec)
                    ra, dec = wcs.all_pix2world(xcentroid, ycentroid, 0)

                    # Format and print the example observation
                    formatted_observation = format_mpc_observation(ra, dec, mag, time_image_taken)
                    print(formatted_observation)

                    # Define the output file path
                    output_file_path = 'formatted_observations.txt'

                    # Open the output file in write mode, wiping any existing content and write the new observation
                    
                    with open(output_file_path, 'w') as output_file:
                        output_file.write(formatted_observation)

                    # Append each observation to the file
                    with open(observations_file_path, 'a') as observations_file:
                        observations_file.write(formatted_observation + "\n")

                    # run astcheck to find objects
                    isClassified = call_astcheck(output_file_path)

                    if isClassified:

                        # Remove fields from header if no objects found
                        del header[f'ID_{id}']
                        del header[f'XCENT_{id}']
                        del header[f'YCENT_{id}']
                        del header[f'MAG_{id}']

                    else:

                        # Add 'id' to unclassified ids list if objects found
                        unclassified_ids.append(id)

            # Update the header with unclassified IDs if any
            if unclassified_ids:

                # add ids to header
                header['UNCLASID'] = ','.join(unclassified_ids)

                # If the object isn't found, find orbit and add to header
                orbit = findOrbit(observations_file_path)
                header[f'ORBIT_{id}'] = (orbit, f'Orbit of unidentified object {id}')

                # Check if the 'unclassified' directory exists, create it if not
                unclassified_path = 'unclassified/'
                if not os.path.exists(unclassified_path):
                    os.makedirs(unclassified_path, exist_ok=True)

                # Move the FITS file to the 'unclassified' directory
                os.rename(INPUT_PATH + file, unclassified_path + file)

            else:
                # Delete file if no objects found
                if os.path.exists(observations_file_path):
                    os.remove(observations_file_path)

            # Save changes to the FITS file
            hdul.flush()