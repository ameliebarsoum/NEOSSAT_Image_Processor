import os
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import subprocess
import warnings
from astropy.utils.exceptions import AstropyWarning
from playwright.sync_api import sync_playwright

# Using Playwright with Python: https://playwright.dev/python/docs/intro
# pip install playwright
# playwright install

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

    with sync_playwright() as playwright:

        # Setup playwright browser
        # Set headless=False if you want to see the browser
        browser = playwright.chromium.launch(headless=True)
        page = browser.new_page()

        # Navigate to the website
        page.goto("https://www.projectpluto.com/fo.htm")

        # Read from observations file
        with open(observations_file_path, 'r') as file:
            text_to_submit = file.read()

        # Fill the text box with content from the text file
        page.fill("#TextArea", text_to_submit)

        # Click the "submit" and wait for the navigation to complete
        with page.expect_navigation():
            page.wait_for_selector("input[type='submit'][value=' Compute orbit and ephemerides ']", state="attached")
            page.click("input[type='submit'][value=' Compute orbit and ephemerides ']")

        # Extract all text from the next page
        result_text = page.inner_text("body")

        # Cleanup
        browser.close()

        # Define the start and end markers
        start_marker = "Orbital elements:"
        end_marker = "Residuals in arcseconds:"

        # Find the position of the markers
        start_pos = result_text.find(start_marker) + len(start_marker)
        end_pos = result_text.find(end_marker)

        # Extract the text between the markers to get the orbital elements
        extracted_text = result_text[start_pos:end_pos].strip()

        # if orbital elements exist, return them
        if extracted_text is not None:
            return extracted_text
        
        # otherwise return N/A
        else:
            return "N/A"

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
                    detra = ra.item() # detected ra as scalar
                    detdec = dec.item() # detected dec as scalar

                    header[f'DETRA_{id}'] = (detra, f'Right Ascension of detected object {id}')
                    header[f'DETDEC_{id}'] = (detdec, f'Declination of detected object {id}')                    

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
                        # COMMENT THIS OUT if you want to receive emails for all detected objects
                        del header[f'ID_{id}'] 
                        del header[f'XCENT_{id}']
                        del header[f'YCENT_{id}']
                        del header[f'MAG_{id}']
                        del header[f'DETRA_{id}']
                        del header[f'DETDEC_{id}']

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
