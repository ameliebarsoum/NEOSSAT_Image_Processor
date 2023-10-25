import requests
from bs4 import BeautifulSoup
import os

BASE_URL = "https://data.asc-csa.gc.ca/users/OpenData_DonneesOuvertes/pub/NEOSSAT/ASTRO/2023/001/"
SAVE_DIR = "fits"  # Update with your path

if not os.path.exists(SAVE_DIR):
    os.makedirs(SAVE_DIR)

response = requests.get(BASE_URL)

# Check if the request was successful
if response.status_code == 200:
    soup = BeautifulSoup(response.content, 'html.parser')
    
    # This assumes the FITS files are directly linked on the page. If there's a directory structure, more navigation is needed.
    for link in soup.find_all('a'):
        file_link = link.get('href')
        
        # Check if the link ends with .fits indicating a FITS file. Adjust as needed.
        if file_link.endswith('.fits'):
            file_url = os.path.join(BASE_URL, file_link)
            
            # Download the FITS file
            file_response = requests.get(file_url, stream=True)
            with open(os.path.join(SAVE_DIR, file_link), 'wb') as file:
                for chunk in file_response.iter_content(chunk_size=8192):
                    file.write(chunk)
