# n1sort.py
#
# REVISION HISTORY:
#  1.00 -  			   - 
#  1.10 - M. Hamza      - 2019 Jun 02 - Changed the working and save directory
#  3.00 - V. Abbasi     - 2019 Sep 20 - Updates for integration
#  4.00 - A. Fagarasanu - 2020 Jan 26 - File processing logic, library directories auto creation, compression, logging  
#  4.00 - V. Abbasi     - 2021 Jun 15 - Avoid creating extra image copy in fits_processor/outgoing/*/previous; reorganized code to ensure headers update and compression/move to FTP_DIR for INCOMPLETE and other skipped files
#
from shutil import copy, move, copyfile
from os import listdir,makedirs,path,walk, mkdir, remove
import gzip
from astropy.io import fits
import datetime
import sys
import os
import re
import logging 
import fnmatch 
import configparser
config = configparser.ConfigParser()
config.read(path.join(path.abspath(path.dirname(__file__)), '../fits_cleaner/config/n1clean.config'))

#GENERAL CONFIGURATION
COMPRESS_RAW = config.getboolean('fitsconfig', 'COMPRESS_RAW')

# DIRECTORIES CONFIGURATION 
WORK_DIR = "fits_processor/outgoing/"
FTP_DIR = "outgoing/"
SAVE_DIR = "fits_cleaner/data/images/"
DARK_LIBRARY = "fits_cleaner/data/dark_library/"
CLEAN_INCOMING = "fits_cleaner/incoming/"

# CLEANER SOFTWARE VERSION
CLEANER_SW_VERSION = "1.0"

# VARIABLE USED TO KEEP TRACK OF WHAT OBSERVER IS BEING PROCESSED
OBSERVER = None

# VARIABLES TO COUNT THE NUMBER OF FITS FILES AND DARK FILES 
FITS_NUM = 0 
DARK_NUM = 0 

#DISPLAY MESSAGE COUNTER
msg = 0

#Creating logger to append sub-info in the main logger file from pipeline launcher file 
	
logger = logging.getLogger("n1 Image Processor")
logger.setLevel(logging.INFO)

#Handler creation and addition to the logger 

log_dir = str(sys.argv[3])

fh = logging.FileHandler(log_dir)

def create_libraries(lib_list):
	for lib in lib_list:
		if not os.path.exists(lib):
			makedirs(lib)

def sort_images(imgdir, outdir):

	formatter = logging.Formatter('%(message)s')
	fh.setFormatter(formatter)

	logger.addHandler( fh )

	"""Creates subfolders in outdir for each size of file in imgdir and copies the files over. Requires outdir to be a valid folder."""
	
	if not imgdir[-1] == '/': imgdir += '/'
	if not outdir[-1] == '/': outdir += '/'

	logger.info("\n========================================")
	logger.info("image dir is " + imgdir)
	logger.info("save dir is " + outdir)
	logger.info("dealing with: ")

	for file in listdir(imgdir):
		logger.info(file)
		
		if file.endswith('.fits'):

			global OBSERVER, FITS_NUM, DARK_NUM, msg

			try:
				if ( (file.endswith('_cord.fits')) | (file.endswith('_cor.fits'))) :
					# Assume that headers must have already been applied in "cord" and "cor" products; do not update headers
					raise TypeError("Already cleaned (cor|cord)")

				dark_image = 0
					
				# OPEN FITS FILE AND READ HEADER	
				hdulist = fits.open(imgdir + file, mode="update")
				xsc,ysc = [-i[0]+i[1]+1 for i in [list(map(int,i.split(':'))) for i in hdulist[0].header["CCDSEC"].strip()[1:-1].split(',')]]
				header = hdulist[0].header
				exp = header["EXPOSURE"]
				size = str(xsc) + 'x' + str(ysc) + '_{:0.2f}s'.format(round(exp, 2))
					
				# ADD META DATA
				header.set("CAL_INFO", "", "##### CALIBRATION")
				header.set("ARCHIVE", "NEOSSat", "Data archive")

				# OBSERVATION TYPE
				if int(header['SHUTTER'][0]) != 0: 
					dark_image = 1
					header.set("OBSTYPE", "DARK", "Observation type (DARK | OBJECT)")
				else:
					dark_image = 0
					header.set("OBSTYPE", "OBJECT", "Observation type (DARK | OBJECT)")
						
				if OBSERVER == "ASTRO":
					header.set("RELEASE", header["DATE-OBS"], "Release date of image")


				if (file.endswith('_clean.fits')):
					# Although 'clean' is obsolete, they might still appear from archives; apply some basic headers to them but do not re-clean
					header.set("CAL_LVL", "CALIBRATED", "Calibration level")
					header.set("PRODUCT", "clean", "raw | cor | cord | clean")
					file_id = file.split("_clean.fits")[0]
					hdulist.flush()
					hdulist.close()
					raise TypeError("Already cleaned (clean)")
							
				else:
					header.set("CAL_LVL", "RAW_STANDARD", "Calibration level")
					header.set("PRODUCT", "raw", "raw | cor | cord | clean")
					file_id = file.split(".fits")[0]
					header.set("OBS_ID", file_id, "Observation ID")

					FITS_NUM += 1
							
							
				hdulist.flush()
				hdulist.close()

				if ( header['IMGSTATE'] == 'INCOMPLETE' ):
					# We want to avoid cleaning any INCOMPLETE image states as they can cause problems in the cleaning stages due to banding
					raise TypeError("Image is marked INCOMPLETE")

					
				#GET YEAR AND DAY FROM FILE TO CREATE RELEVANT DIRECTORIES 
				file_name = str(file)
				date = file_name.split("_")[2]
				year = date[0:4]
				day = date[4:7]

				# CREATE A CLEAN LIST OF ALL IMAGES SENT TO CLEANER
				cleanlist = open(CLEAN_INCOMING + "current.cleanlist", "a+")

				
				if dark_image == 1:
						
					DARK_NUM += 1
							
					# CREATE THE DIRECTORIES WITH THE RASTER SIZE IF NOT PRESENT 
					if not size in listdir(DARK_LIBRARY): 
						mkdir(DARK_LIBRARY + size)
						
					# CREATE DIRECTORIES IF THE YEAR IS NOT PRESENT 
					if not year in listdir(DARK_LIBRARY + size): 
						mkdir(DARK_LIBRARY + size + "/" + year)
						
					# CREATE DIRECTORIES IF THE DAY OF YEAR IS NOT PRESENT	 					
					if not day in listdir(DARK_LIBRARY + size + "/" + year):
						mkdir(DARK_LIBRARY + size + "/" + year + "/" + day)

					# ONCE ALL THE DIRECTORIES ARE CREATED THEN COPY THE FILE TO ITS APPROPRIATE SORTED LOCATION
					sortloc_dark = DARK_LIBRARY + size + "/" + year + "/" + day + "/" + "." .join(file.split('.')[:-1]) + ".fits"	
					copy(imgdir + file, sortloc_dark)

					# APPEND TO CLEAN LIST
					cleanlist.write(f'{sortloc_dark}\t\tDARK\n')	
				else:

					# CREATE DAY OF YEAR DIRECTORIES IF ITS NOT PRESENT UNDER THE YEAR DIRECTORY 
					if year in listdir(outdir): 
						if not day in listdir(outdir + year):
							mkdir(outdir + year + "/" + day)
									
					else:
						# CREATE YEAR AND DAY OF YEAR DIRECTORY IF BOTH ARE NOT PRESENT 
						mkdir(outdir + year)
						mkdir(outdir + year + "/" + day)
						
					# ONCE ALL THE DIRECTORIES ARE CREATED THEN COPY THE FILE TO ITS APPROPRIATE SORTED LOCATION	
					sortloc_light = outdir + year + "/" + day + "/" + "." .join(file.split('.')[:-1]) + ".fits"					
					copy(imgdir + file, sortloc_light)

					# ADD LIGHT TO CORRESPONDING LIGHTLIST. CREATE IF NOT PRESENT
					lightlist = open(DARK_LIBRARY + "lightlists/" + size + ".lightlist", "a+")
					lightlist.write(f'{outdir}{year}/{day}/{file_id}_cor.fits\n')
					lightlist.close()

					# APPEND TO CLEAN LIST
					cleanlist.write(f'{sortloc_light}\t\t\t\tLIGHT\n')
					
			except TypeError as e:
                            logger.info("n1sort.py - Skipping cleaner copy for " + file + " due to: " + "TypeError: " + str(e))

			finally:

				# CREATE PREVIOUS DIRECTORY # VA 15-June-2021 removed use of this directory
				if not os.path.exists(imgdir + "previous/"):
					mkdir(imgdir + "previous/")
					
				# IF THE OBSERVER DIRECTORY IS NOT CREATED IN THE FTP DIRECTORY THEN CREATE ONE
				if not os.path.exists(FTP_DIR + OBSERVER):
					mkdir (FTP_DIR + OBSERVER)
					mkdir (FTP_DIR + OBSERVER + "/previous")

                                # ONCE WE HAVE SORTED THE FILE, WE WILL COPY THE FILE TO ITS "previous" DIRECTORY  # VA 15-June-2021: this copy seems unnecessary
				#copy(imgdir + file, imgdir + "previous/" + file)
					
				# MOVE THE RAW FILES TO THE FTP DIRECTORY AND IF COMPRESSION IS ON, COMPRESS BEFORE MOVING  
				if COMPRESS_RAW:
					
					if msg == 0:
						print("\nCompressing RAW files...")
					
					msg += 1

					input = open(imgdir + file, 'rb')
					s = input.read()
					input.close()

					output = gzip.GzipFile(FTP_DIR + OBSERVER + "/" + file + ".gz", 'wb')
					output.write(s)
					output.close()

					remove(imgdir + file)
				else:
					move(imgdir + file, FTP_DIR + OBSERVER + "/" + file)
				pass
		else:
			if file == "ASTRO" or file == "HEOSS":
				if not file in "previous":
					if not file in listdir(outdir): 
						logging.info("creating " + outdir + file)
						logging.info(" and " + outdir + file + "/previous")
						mkdir(outdir + file)
					OBSERVER = file
					sort_images(imgdir + file, outdir + file)				

		
def main():
	global WORK_DIR, FTP_DIR
	global SAVE_DIR, DARK_LIBRARY
	global CLEAN_INCOMING
	global FITS_NUM, DARK_NUM
	
	base_dir = str(sys.argv[2])
	
	WORK_DIR = base_dir + WORK_DIR
	FTP_DIR = base_dir + FTP_DIR
	SAVE_DIR = base_dir + SAVE_DIR
	DARK_LIBRARY = base_dir + DARK_LIBRARY
	CLEAN_INCOMING = base_dir + CLEAN_INCOMING

	lib_list = [WORK_DIR, (WORK_DIR + "ASTRO/"), (WORK_DIR + "HEOSS/"), FTP_DIR, SAVE_DIR, DARK_LIBRARY, \
		(base_dir + "fits_cleaner/incoming/"), (base_dir + "fits_cleaner/logs/"), (base_dir + "fits_cleaner/outgoing/ASTRO/"), \
		(base_dir + "fits_cleaner/outgoing/HEOSS"), (base_dir + "incoming/"), (base_dir + "fits_stats/"), DARK_LIBRARY + "/lightlists/"]
	create_libraries(lib_list)

	sort_images(WORK_DIR, SAVE_DIR)

	logger.info("\n")

	formatter = logging.Formatter('%(asctime)s - %(message)s')
	fh.setFormatter(formatter)

	logger.addHandler( fh )
	
	logger.info( "Finished. " + str(FITS_NUM) + " files sorted. " + str(DARK_NUM) + " dark(s) added.\n" )

	if FITS_NUM == 0:
		print("\nn1sort: No incoming files; ending...\n")
	
	
if __name__ == "__main__":
	main()



