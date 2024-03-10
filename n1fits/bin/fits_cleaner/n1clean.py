# n1clean.py
#
# REVISION HISTORY:
#  1.00 - M. Hamza      - 2019 Aug 19 - Initial release of the cleaner program. -
#  3.00 - V. Abbasi     - 2019 Sep 20 - Miscellaneous fixes in support of integration. -
#  3.10	- G. Ko         - 2019        - Metadata trimming adjustments, dark library. -		  
#  4.00 - A. Fagarasanu - 2020 Jan 26 - Major processing logic overhaul, Mprocessing, HThreading,
#                                       WCS reference adjustments, error catching, logging, various fixes and changes. -
#  4.00 - V. Abbasi     - 2021 June 15 - Ensure getimage_dim() is called when processing new darks to fix issue when mixing processed and unprocessed darks in a list
#                                        + touch-ups in logging and reduction of extra saved image copies
#  4.00 - L. Kuhn       - 2021 Dec 16 - Various fixes and changes including implementing try catches to several key library functions, adapting code to work with J. Rowe image cleaning library, addition of the darkmode parameter, and other minor touch-ups. 
#
from os import listdir, mkdir, path, makedirs
from astropy.wcs import WCS, FITSFixedWarning, Wcsprm
from astropy.nddata import Cutout2D
import fnmatch
import multiprocessing as mp
from multiprocessing import Pool, pool
from tqdm import tqdm , notebook #For monitoring process of big jobs
from tqdm.contrib.concurrent import process_map, thread_map 
from astropy.io import fits
import re #to extract trim section for FITS header
import matplotlib #ploting
import math
#matplotlib.use("Agg") #workaround to stop python icon bouncing in os/x
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm #for better display of FITS images
import time as t 
import datetime
#import neossatlib as neo
import neossat as neo # used for the new library format
import warnings
import numpy as np
from statistics import median
from copy import deepcopy, copy
import datetime
import shutil
import sys
#import ast #For converting string representations of lists into literal lists
import json
import os
import logging
import subprocess
import configparser
import itertools
import traceback
import platform
from platform import uname
config = configparser.ConfigParser()
config.read(path.join(path.abspath(path.dirname(__file__)), 'config/n1clean.config'))

# --- CONFIGURATION ---
VERSION = config.get('fitsconfig', 'VERSION')
SW_REF = config.get('fitsconfig', 'SW_REF')
SAVE_DIR = 			"outgoing/"	 					# DIRECTORY TO SAVE CLEANED IMAGES
IMAGES_DIR = 		"data/images/"					# DIRECTORY WHERE ALL SORTED IMAGES ARE
DARKS_DIR = 		"data/dark_library/"
LOG_DIR = 			"logs/"	
MASTER_DARK = 		None

SAVE_MDARK = config.getboolean('darkprocess', 'SAVE_MDARK')
MIN_DARKS = config.get('darkprocess', 'MIN_DARKS')
MAX_DARKS = config.get('darkprocess', 'MAX_DARKS')
MAX_DARK_SPAN = config.get('darkprocess', 'MAX_DARK_SPAN')
DAY_RANGE_LIMIT = config.get('darkprocess', 'DAY_RANGE_LIMIT')
NODARKS_SPAN = config.get('darkprocess', 'NODARKS_SPAN')
SAA_NW_BOUND = config.get('darkprocess', 'SAA_NW_BOUND')
SAA_SE_BOUND = config.get('darkprocess', 'SAA_SE_BOUND')
darkmode = config.get('darkprocess', 'DARK_MODE')
pixfrac = float(config.get('darkprocess', 'PIX_FRAC'))

SCALE_FACTOR_CORD = config.get('scaling', 'SCALE_FACTOR_CORD')
SCALE_FACTOR_COR = config.get('scaling', 'SCALE_FACTOR_COR')
SCALE_FACTOR_DARKS = config.get('scaling', 'SCALE_FACTOR_DARKS')
SCALE_COR = config.getboolean('scaling', 'SCALE_COR')
SCALE_DARKS = config.getboolean('scaling', 'SCALE_DARKS')
PROGRESS_BAR = config.getboolean('fitsconfig', 'PROGRESS_BAR')
MPROCESSING = config.getboolean('fitsconfig', 'MPROCESSING') # Whether Multi Processing is enabled, pulls choice from config file
nprocesses = mp.cpu_count() 								# Gives number of processing cores available
HTHREADING = config.getboolean('fitsconfig', 'HTHREADING') # Whether Hyper Threading is enabled, pulls choice from config file
nthreads = 64

MEMMAP_ERROR_OS = (config.get('memmap', 'MEMMAP_ERROR_OS')).split(',')

#--- TRACKING FILES ---
INCOMING = "incoming/current.cleanlist"		#File for listing incooming files to clean from n1sort.py
INCOMING_LIGHTS = []						#Stores all of the incoming light frames from the incoming cleanlist
INCOMING_DARKS = []							#Stores all of the incoming dark frames from the incoming cleanlist			
SAA_DARKS = 0								#Number of darks that were found to be inside the define SAA region

# --- FITS IMAGE PARAMETERS ---
#only for plotting the images
bpix = -1.0e10                   			#value to mark bad pixels. Any pixel *below* bpix is considered invalid.
sigscalel = 1.0                 			#low bounds for display clipping.  Keep small to keep background 'dark'
sigscaleh = 1.0                  			#high bounds for display clipping.

# --- EARTH PARAMETERS ---
a_Earth = 6378137							#Apprx. Semi-Major Axis in meters of Earth (radius at equator where it is larger) in WGS84
b_Earth = 6356752							#Apprx. Semi-Minor Axis in meters of Earth (radius at poles where it is smaller) in WGS84

e1_Earth = math.sqrt( (math.pow(a_Earth, 2) - math.pow(b_Earth, 2)) ) / math.pow(a_Earth, 2) #Eccentricity of Earth at equator

e2_Earth = math.sqrt( (math.pow(a_Earth, 2) - math.pow(b_Earth, 2)) ) / math.pow(b_Earth, 2) #Eccentricity of Earth at equator

p_Sat = 0.0  								#semi-latus rectum of Satellite's Location, calculated later
theta_Sat = 0.0   						    #angle of satellite                             

#--- PARAMETERS FOPR FOURIER DECOMPOSITION ---
snrcut = 10.0					 			#When to stop decomposition
fmax = 2 						 			#Maximum number of frequencies in model
xoff = 0 						 			#pixel offsets
yoff = 0
T = 8 							 			#Oversampling
info = 0						 			# 0 - no plots, no output, 1 - no plots, output, 2 - plots, output

#--- PARAMETER FOR EITHER cor OR cord ---
CLEAN_TYPE = None

#--- DARK BINNING PARAMETERS FOR CHOOSING IDEAL MASTERDARK
temp_low = float(config.get('darkprocess', 'TEMP_DARK_LOW'))	#low end of temperature binning (MIN POTENTIAL, Kelvin)
temp_high = float(config.get('darkprocess', 'TEMP_DARK_HIGH')) #high end of temperature binning (MAX POTENTIAL, Kelvin)
standard_bin_width = float(config.get('darkprocess', 'STANDARD_DARK_BIN_WIDTH')) #temperature width of an ideal, 'standard' bin in degrees kelvin

#--- DICTIONARY OF BINNED DARK_LIST MASTERDARKS + MED TEMP
MASTERDARK_DICTIONARIES = {}

#--- DICTIONARY TO KEEP TRACK OF DARK LIST BINNING USED PER LIGHT
DARKS = {}

#--- GLOBAL VARIABLE USED TO KEEP TRACK OF CURRENT OBSERVER
CURRENT_OBSERVER = "ASTRO"

#--- STORING TERMINAL STDOUT ---
terminal = sys.stdout

#--- CLEANER LOG FILE ---
cleaner_log = None

#--- IDENTIFIES USER'S OPERATING SYSTEM AND SETS MEMMAP PARAMETER ACCORDINGLY ---
if any(OS in uname().release for OS in MEMMAP_ERROR_OS):
	MEMMAP_PARAMETER = False
else:
	MEMMAP_PARAMETER = True
        
#--- IDENTIFIES USER'S CURRENT WORKING DIRECTORY ---
workdir = os.getcwd()


def add_cord_header(file, dark_list_used, DARKTMIN, DARKTMAX, DARKTMED):

	global VERSION
        
	# OPEN THE FITS FILE TO RETIREVE THE HEADER
	hdul = fits.open(file, mode="update", memmap= MEMMAP_PARAMETER)
	header = hdul[0].header


	# ADD HEADER INFORMATION RELATED TO "cord" FITS
	header.set("CAL_LVL", "CALIBRATED", "Calibrated level")
	header.set("PRODUCT", "cord", "Clipping, overscan & dark correction applied")
	header.set("RUNID",str(datetime.datetime.utcnow()).split(".")[0].replace(" ", "T"), "Timestamp of data processing run")
	header.set("SW_REF", SW_REF, "Data processing software")


	# APPEND CLEANER SW NAME AND VERSION TO THE HEADER, IF NOT ALREADY DONE
	if "n1clean.py" not in header["CONV_SW"]:
		header.set("CONV_SW", header["CONV_SW"] + ", n1clean.py" )  # Library name changed
		header.set("CONV_VER", VERSION)


	# ADD HEADER FOR TEMPERATURE RANGE USED FOR DARK SUBTRACTION
	header.set("DARKTMIN", DARKTMIN, "[K] Minimum dark temperature used")
	header.set("DARKTMAX", DARKTMAX, "[K] Maximum dark temperature used")
	header.set("DARKTMED", DARKTMED, "[K] Median dark temperature used")

	# ADD HEADER FOR DARK MODE USED FOR DARK PROCESS 
	header.set("DARKMODE", darkmode, "Method to scale masterdark to science image")
    
	dark_list_used
	incrementor = 1

	# ADD ALL THE DARK FITS TO THE IMAGE HEADER
	for dark in dark_list_used:
		file_id = os.path.splitext(os.path.basename(dark))[0].split(".fits")[0]
		if incrementor > 999:
			header.set("DARK" + str(incrementor).zfill(4), file_id)
		else:
			header.set("DARK_" + str(incrementor).zfill(3), file_id)
		incrementor += 1

	# MARK THE END OF THE PRIMARY HEADER
	#header.append('END','')

	# SAVE THE HEADER CHANGES TO THE FITS FILE.
	hdul.flush()

	# CLOSE THE FILE ONCE THE INFORMATION IS PERSISTED
	hdul.close()


def add_cor_header(file):

	# OPEN THE FITS FILE TO RETIREVE THE HEADER
	hdul = fits.open(file, mode="update")
	header = hdul[0].header

	# ADD HEADER INFORMATION RELATED TO "cor" FITS
	header.set("CAL_LVL", "CALIBRATED")
	header.set("PRODUCT", "cor", "Clipping & overscan correction applied")
	header.set("RUNID",str(datetime.datetime.utcnow()).split(".")[0].replace(" ", "T"), "Timestamp of data processing run")
	header.set("SW_REF", SW_REF, "Data processing software")

	# APPEND CLEANER SW NAME AND VERSION TO THE HEADER, IF NOT ALREADY DONE
	if "n1clean.py" not in header["CONV_SW"]:
		header.set("CONV_SW", header["CONV_SW"] + ", n1clean.py" )  # Library name changed
		header.set("CONV_VER", VERSION)

	# WE ARE GOING TO CLIP OVERSCAN, SO WE WILL CHANGE OVERSCAN TO 0
	header.set("OVERSCAN", 0)

	# MARK THE END OF THE PRIMARY HEADER
	#header.append("END",'')

	# SAVE THE HEADER CHANGES TO THE FITS FILE
	hdul.flush()

	# CLOSE THE FILE ONCE THE INFORMATION IS PERSISTED
	hdul.close()


def get_size(file):
	# OPEN THE FITS FILE TO RETRIEVE HEADER
	hdulist = fits.open(file)
	hdr = hdulist[0].header
	exp = hdr["EXPOSURE"]
	hdulist.close()
	try:
		trim,btrim,ysc,xsc,xov,yov = neo.getimage_dim(file)
	except Exception as e:
		print ("Could not get size for", str(file), "due to error:", e)
		print (traceback.format_exc())
	else:  
		size = str(xsc) + 'x' + str(ysc) + '_{:0.2f}s'.format(round(exp, 2))
		return size

def get_temp(file):
	#OPEN THE FITS FILE TO RETRIEVE HEADER
	hdulist = fits.open(file)
	hdr = hdulist[0].header
	temp = hdr["TEMP_CCD"]
	hdulist.close()

	return temp



def cropped_metadata(file, trim):
	'''
		CHANGE DISPLAY METADATA TO FIT CROPPED IMAGE WITHOUT OVERSCAN
		(Doesn't change the actual image, just takes care of a display 
		error some fits viewers would produce)
	'''

	data = neo.read_fitsdata(file) # For changing display metadata
	sh = data.shape

	hdul = fits.open(file, mode = "update")
	hdr = hdul[0].header

	hdr.set("BIASSEC", '', "No Overscan Region Remaining") #NO OVERSCAN REGION ANYMORE
	hdr["TRIMSEC"] = f"[1:{sh[1]},1:{sh[0]}]"
	hdr["DATASEC"] = f"[1:{sh[1]},1:{sh[0]}]"

	VERSION = float(hdr["CONV_VER"])

	if 'CRPIX1'in hdr:
		#Updating WCS pixel reference 
		hdr['CRPIX1'] = hdr['CRPIX1'] - trim[0] + 1
		hdr['CRPIX2'] = hdr['CRPIX2'] - trim[2] + 1

	hdul.flush()
	hdul.close()


def is_in_SAA(file):
	''' 
		This function returns a boolean evaluation of whether the provided image file 
		was taken within the configurable definition of the South Atlantic Anomaly (SAA).
		True if coordinates are within defined region
		False if coordinates are not within defined region
	'''
			
	global p_Sat, theta_Sat, SAA_NW_BOUND, SAA_SE_BOUND

	SAA_NW_BOUND = json.loads(str(SAA_NW_BOUND))
	SAA_SE_BOUND = json.loads(str(SAA_SE_BOUND))

	#OPEN THE FITS FILE TO RETRIEVE HEADER
	hdulist = fits.open(file)
	hdr = hdulist[0].header
	
	if 'GEO_LAT' in hdr and 'GEO_LONG' in hdr:

		longitude = float(hdr['GEO_LONG'])
		latitude = float(hdr['GEO_LAT'])

	else:
		x_ECEF = hdr['EPOS2_1'] * 1000
		y_ECEF = hdr['EPOS2_2'] * 1000
		z_ECEF = hdr['EPOS2_3'] * 1000

		p_Sat = math.sqrt( (math.pow(x_ECEF, 2)) + (math.pow(y_ECEF, 2)) )
		theta_Sat = (z_ECEF * a_Earth) / (p_Sat * b_Earth)

		longitude = (math.atan2(y_ECEF, x_ECEF)) * (180/math.pi)
		latitude = math.atan2((z_ECEF + math.pow(e2_Earth, 2) * b_Earth * math.pow(math.sin(theta_Sat), 2)), \
			(p_Sat - math.pow(e1_Earth, 2) * a_Earth * math.pow(math.cos(theta_Sat), 3))) * (180/math.pi)

	#--- Return the boolean comparison of LONG/LAT with SAA definitions
	return longitude >= SAA_NW_BOUND[1] and longitude <= SAA_SE_BOUND[1] and \
		latitude >= SAA_NW_BOUND[0] and latitude <= SAA_SE_BOUND[0]


def clean_fits(file):
	timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")
	#print (timestamp + "Starting to clean " + os.path.basename(file) )

	size = get_size(file)

	global CURRENT_OBSERVER

	# IF CLEAN_TYPE = 0, THEN WE ARE PROCESSING cor
	if CLEAN_TYPE == 0:

		# SET UP THE FILE NAME FOR THE CLEANED IMAGE
		base = os.path.basename(file)
		x = os.path.splitext(base)
		file_name = x[0] + "_cor.fits"
		new_file = SAVE_DIR + CURRENT_OBSERVER + "/" + file_name
		directory = file.split(base)[0]

                # COPY THE ORIGINAL FILE TO /previous (VA: 15-June-2021 - Removed extra copy, as n1sort has already copied original file)
		#shutil.copy(file, directory + "previous/" + base)

		# ADD NECESSARY HEADERS TO THE FILE
		add_cor_header(file)

		try:
			trim,btrim,xsc,ysc,xov,yov = neo.getimage_dim(file)
			# PERFROM THE CLEANING USING THE neossat LIBRARY
			timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")
			print (timestamp + "Cleaning (COR) file: " + base)
			scidata_cor = neo.clean_sciimage(file,[],xsc,ysc,xov,yov,snrcut,fmax,xoff,yoff,T,info,bpix,darkmode,pixfrac=pixfrac)
		except Exception as e:
			timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")
			print (timestamp + "Skipped", base, "(COR) due to error:", e)
			print (traceback.format_exc())
			os.remove(file)
		else:      
			header = fits.getheader(file)
			if info > 0:
				print("Corrected images, sorting them to their proper directories and statuses")
	
			header['BSCALE']=1, ' ' #make sure BZERO and BSCALE are set
			header['BZERO']=0, 'zero def of pixel (for unsigned short)'  

			# SCALING
			if (SCALE_COR):
				if (SCALE_FACTOR_COR == 'int16'):
					header['BZERO']=32768, 'zero def of pixel (for unsigned short)'
				hdu = fits.PrimaryHDU(scidata_cor)
				hdu.scale(SCALE_FACTOR_COR) #Scaling
				i=0
				for h in header:
					if h!='SIMPLE' and h!='BITPIX' and h!='NAXIS' and h!='NAXIS1' and h!='NAXIS2' and h!='EXTEND':
						hdu.header.append((h,header[i]))
					hdu.header.comments[h] = header.comments[h]
					i+=1
				hdu.writeto(new_file, overwrite=True)
			else:
				# WRITE THE cor FILE TO cleaner/outgoing AND ALSO WRITE IT TO cord_pending TO PROCESS IT LATER TO PRODUCE CORDS
				fits.writeto(new_file, scidata_cor, header, overwrite=True)

			cropped_metadata(new_file, trim)


			fits.writeto(directory + file_name, scidata_cor, header, overwrite=True)

			# ALSO WRITE CHANGES TO /previous
			fits.writeto(directory + "previous/" + file_name, scidata_cor, header, overwrite=True)

			# REMOVE THE ORIGINAL FILE SINCE IT HAS ALREADY BEEN COPIED TO /previous
			os.remove(file)


	elif CLEAN_TYPE == 1:
		base = os.path.basename(file)
		timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")
		print (timestamp + "Cleaning (CORD) file: " + base )

		# RETRIEVE THE MASTER DARK FOR THIS FILE USING THE RASTER SIZE
		#-------------------------------------------------------------
		prev_master_dark = ""
		best_master_dark = ""

		prev_master_dark_temp = 0.0
		best_master_dark_temp = 0.0

		light_temp = get_temp(file)
		temp_diff = 100.0
		bin_num = 0

		prev_bin_used = 0
		best_bin_used = 0

		for masterdark_list in MASTERDARK_DICTIONARIES[size]:
			if len(masterdark_list) != 0:
				if abs(light_temp - masterdark_list[0]) < temp_diff:

					temp_diff = abs(light_temp - masterdark_list[0])

					prev_master_dark = best_master_dark
					best_master_dark = masterdark_list[1]

					prev_master_dark_temp = best_master_dark_temp
					best_master_dark_temp = masterdark_list[0]

					prev_bin_used = best_bin_used
					best_bin_used = bin_num

			bin_num += 1
		
		master_dark = ""
		master_dark_temp = ""
		bin_used = 0

		if (best_master_dark_temp > light_temp) and (prev_master_dark_temp < light_temp) \
			and ((light_temp - prev_master_dark_temp) < standard_bin_width):

			master_dark = prev_master_dark
			master_dark_temp = prev_master_dark_temp
			bin_used = prev_bin_used

		else:

			master_dark = best_master_dark
			master_dark_temp = best_master_dark_temp
			bin_used = best_bin_used

		#-------------------------------------------------------------

		len_bin_used = len(DARKS[size][bin_used])
		first_dark_used = DARKS[size][bin_used][0]
		last_dark_used = DARKS[size][bin_used][len_bin_used-1]

		dark_bin = DARKS[size][bin_used]

		DARKTMIN = get_temp(first_dark_used)
		DARKTMAX = get_temp(last_dark_used) 
		DARKTMED = master_dark_temp

		# ADD "cord" related header
		add_cord_header(file, dark_bin, DARKTMIN, DARKTMAX, DARKTMED)

		# CLEAN THE IMAGE BY APPLYING DARK CORRECTION USING THE neossat library

		try:
			trim,btrim,xsc,ysc,xov,yov = neo.getimage_dim(file)
			# PERFROM THE CLEANING USING THE neossat LIBRARY
			scidata_cord = neo.clean_sciimage(file,master_dark,xsc,ysc,xov,yov,snrcut,fmax,xoff,yoff,T,info,bpix,darkmode,pixfrac=pixfrac)
		except Exception as e:
			timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")
			print (timestamp + "Skipped", base, "(CORD) due to error:", e)
			print (traceback.format_exc())
		else:
			# SET UP THE FILE NAME FOR THE CLEANED IMAGE
			x = os.path.splitext(base)[0]
			x = x.split("_cor")[0]
			file_name = x + "_cord.fits"
			newfile = SAVE_DIR + CURRENT_OBSERVER + "/" + file_name
			directory = file.split("cord_pending/" + base)[0]
            
			#CHANGE DISPLAY METADATA FOR THE CROPPED IMAGE
			cropped_metadata(file, trim)

			header = fits.getheader(file)

			header.insert('IMAGE', ('BZERO',0,'zero def of pixel (for unsigned short)')) #make sure BZERO and BSCALE are set
			header.insert('BZERO', ('BSCALE',1,' '))
            
			# SCALING
			if (SCALE_FACTOR_CORD == 'int16'):
				header['BZERO']=32768, 'zero def of pixel (for unsigned short)'''
			hdu = fits.PrimaryHDU(scidata_cord)
			hdu.scale(SCALE_FACTOR_CORD)
			i=0
			for h in header:
				if h!='SIMPLE' and h!='BITPIX' and h!='NAXIS' and h!='NAXIS1' and h!='NAXIS2' and h!='EXTEND':
					hdu.header.append((h,header[i]))
				hdu.header.comments[h] = header.comments[h]
				i+=1
            
			# WRITE THE CORD FILE TO fits_cleaner/outgoing
			hdu.writeto(newfile,overwrite=True)
            
			timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")
			print (timestamp + "Cleaned", base)
            
		# WRITE THE CORD FILE TO PREVIOUS
		#fits.writeto(directory + "previous/" + file_name, scidata_cord, header, overwrite=True)

		# REMOVE THE "cor" FILE ONCE WE ARE DONE PROCESSING IT
		os.remove(file)


def get_master_dark(dark_list):
	'''
		Returns a master dark created from the provided dark list 
		(or provided dark list bin)
	'''

	global CURRENT_OBSERVER

	if (MPROCESSING):
		pool = mp.Pool(processes=nprocesses)
	elif (HTHREADING):
		pool = mp.pool.ThreadPool(processes=nthreads)

	alldarkdata = []
	if (PROGRESS_BAR):
		with tqdm(total=len(dark_list), miniters=1, mininterval=0, leave=True) as pbar:
			for darkfile in dark_list:

				if darkfile.endswith("_cor.fits"): # Addition for using already processed darks
					print("Using existing processed dark: ", os.path.basename(darkfile) ) 
					darkdata = neo.read_fitsdata(darkfile)
					print(f": darkmin:{np.min(darkdata)},darkmax:{np.max(darkdata)},darkmean:{np.mean(darkdata)},darkmed:{np.median(darkdata)} \n")
					alldarkdata.append(darkdata)

				else:
					print("Processing new dark:", darkfile,'\n')
                   
					try:
						trim,btrim,xsc,ysc,xov,yov=neo.getimage_dim(darkfile)
						if (MPROCESSING or HTHREADING):
							dark_cor = pool.apply_async(neo.darkprocess, args=(workdir,darkfile,xsc,ysc,xov,yov,snrcut,fmax,xoff,yoff,T,bpix))
							alldarkdata.append(dark_cor.get())
							print(f"darkmin:{np.min(dark_cor.get())},darkmax:{np.max(dark_cor.get())},darkmean:{np.mean(alldarkdata)},darkmed:{np.median(alldarkdata)} \n")
						else:
							dark_cor = neo.darkprocess(workdir,darkfile,xsc,ysc,xov,yov,snrcut,fmax,xoff,yoff,T,bpix)
							alldarkdata.append(dark_cor)
							print(f"darkmin:{np.min(dark_cor)},darkmax:{np.max(dark_cor)},darkmean:{np.mean(alldarkdata)},darkmed:{np.median(alldarkdata)} \n") 
					except Exception as e:
						print ("Skipped", darkfile, "due to error:", e)
						print (traceback.format_exc())
					else:
						# Saving processed darks for future processing                    
						base = os.path.basename(darkfile)
						x = os.path.splitext(base)[0]
						filename = x + "_cor.fits"
						newfile = SAVE_DIR + CURRENT_OBSERVER + "/" + filename
						header = fits.getheader(darkfile)

						header['BSCALE']=1, ' ' # Make sure BZERO and BSCALE are set
						header['BZERO']=0, 'zero def of pixel (for unsigned short)'
                        
						if (SCALE_DARKS):
							# SCALING
							if (SCALE_FACTOR_DARKS == 'int16'):
								header['BZERO']=32768, 'zero def of pixel (for unsigned short)'
							if (MPROCESSING or HTHREADING):
								hdu = fits.PrimaryHDU(dark_cor.get())
							else:
								hdu = fits.PrimaryHDU(dark_cor)
							hdu.scale(SCALE_FACTOR_DARKS) 
							i=0
							for h in header:
								if h!='SIMPLE' and h!='BITPIX' and h!='NAXIS' and h!='NAXIS1' and h!='NAXIS2' and h!='EXTEND':
									hdu.header.append((h,header[i]))
								hdu.header.comments[h] = header.comments[h]
								i+=1
							hdu.writeto(newfile,overwrite=True)
						else:
							if (MPROCESSING or HTHREADING):
								fits.writeto(newfile,dark_cor.get(),header,overwrite=True)
							else:
								fits.writeto(newfile,dark_cor,header,overwrite=True)

						add_cor_header(newfile)
						cropped_metadata(newfile, trim)

						date = filename.split("_")[2]
						year = date[0:4]
						day = date[4:7]

						copy_dest = DARKS_DIR + get_size(darkfile) + "/" + year + "/" + day + "/PROCESSED_DARKS/"
						if not os.path.isdir(copy_dest):
							makedirs(copy_dest)
						shutil.copy(newfile, copy_dest)

						pbar.update()
						pbar.clear()                      

	else:

		# TQDM PROGRESS BAR IS NOT ITERATED THROUGH IN THIS INSTANCE

		for darkfile in dark_list:


			if darkfile.endswith("_cor.fits"): # Addition for using already processed darks
				print(f"using previously processed dark: {os.path.basename(darkfile)} \n") 
				darkdata = neo.read_fitsdata(darkfile)
				print(f"darkmin:{np.min(darkdata)},darkmax:{np.max(darkdata)},darkmean:{np.mean(darkdata)},darkmed:{np.median(darkdata)} \n") #gko
				alldarkdata.append(darkdata)

			else:
				print("Processing new dark:", darkfile,'\n')

				try:
					trim,btrim,xsc,ysc,xov,yov=neo.getimage_dim(darkfile)
					if (MPROCESSING or HTHREADING):
						dark_cor = pool.apply_async(neo.darkprocess, args=(workdir,darkfile,xsc,ysc,xov,yov,snrcut,fmax,xoff,yoff,T,bpix))
						alldarkdata.append(dark_cor.get())
						print(f"darkmin:{np.min(dark_cor.get())},darkmax:{np.max(dark_cor.get())},darkmean:{np.mean(alldarkdata)},darkmed:{np.median(alldarkdata)} \n")
					else:
						dark_cor = neo.darkprocess(workdir,darkfile,xsc,ysc,xov,yov,snrcut,fmax,xoff,yoff,T,bpix)
						alldarkdata.append(dark_cor)
						print(f"darkmin:{np.min(dark_cor)},darkmax:{np.max(dark_cor)},darkmean:{np.mean(alldarkdata)},darkmed:{np.median(alldarkdata)} \n") 
				except Exception as e:
					print ("Skipped", darkfile, "due to error:", e)
					print (traceback.format_exc())
				else:   
					# Saving processed darks for future processing #gko

					base = os.path.basename(darkfile)
					x = os.path.splitext(base)[0]
					filename = x + "_cor.fits"
					newfile = SAVE_DIR + CURRENT_OBSERVER + "/" + filename
					header = fits.getheader(darkfile)

					header['BSCALE']=1, ' ' #make sure BZERO and BSCALE are set
					header['BZERO']=0, 'zero def of pixel (for unsigned short)' 

					if (SCALE_DARKS):
						#SCALING
						if (SCALE_FACTOR_DARKS == 'int16'):
							header['BZERO']=32768, 'zero def of pixel (for unsigned short)'
						if (MPROCESSING or HTHREADING):
							hdu = fits.PrimaryHDU(dark_cor.get())
						else:
							hdu = fits.PrimaryHDU(dark_cor)
						hdu.scale(SCALE_FACTOR_DARKS) 
						i=0
						for h in header:
							if h!='SIMPLE' and h!='BITPIX' and h!='NAXIS' and h!='NAXIS1' and h!='NAXIS2' and h!='EXTEND':
								hdu.header.append((h,header[i]))
							hdu.header.comments[h] = header.comments[h]
							i+=1
						hdu.writeto(newfile,overwrite=True)
					else:
						if (MPROCESSING or HTHREADING):
							fits.writeto(newfile,dark_cor.get(),header,overwrite=True)
						else:
							fits.writeto(newfile,dark_cor,header,overwrite=True)

					add_cor_header(newfile)
					cropped_metadata(newfile, trim)

					date = filename.split("_")[2]
					year = date[0:4]
					day = date[4:7]

					copy_dest = DARKS_DIR + get_size(darkfile) + "/" + year + "/" + day + "/PROCESSED_DARKS/"
					shutil.copy(newfile, copy_dest)

	# Close instance of mp and combine processing pool 
	if (MPROCESSING or HTHREADING):
		pool.close()
		pool.join()
	try:
		darkavg = neo.combinedarks(alldarkdata, darkmode=darkmode, pixfrac=pixfrac)
	except Exception as e:
		print ("Unable to combine darks due to error:", e)
		print (traceback.format_exc())
	else: 
		if (SAVE_MDARK):
			mastername = "masterdark" + get_size(darkavg) + datetime.datetime.utcnow().strftime("%Y%j%H%M%S") + ".fits"
			mastername = SAVE_DIR + CURRENT_OBSERVER + '/' + mastername
			hdu = fits.PrimaryHDU(darkavg)
			hdu.writeto(mastername, overwrite=True)

		# Clearing the memory
		del alldarkdata

		return darkavg

    
def produce_cor(YEAR, DAY):
	print(f"starting produce_cor({YEAR}, {DAY}): \n")
	global CLEAN_TYPE
	CLEAN_TYPE = 0

	# LOOP THROUGH EACH OF THE OBSERVERS
	for observer in listdir(IMAGES_DIR):

		light_list = []
		print ("looking in " + IMAGES_DIR + observer + "/" + YEAR + "/" + DAY)

		if os.path.isfile(IMAGES_DIR + observer + "/" + YEAR + "/" + DAY + ".tar.lrz"):

			sys.stdout = terminal
			print(f"\nDecompressing -- images/" + observer + "/" + DAY + ".tar.lrz")
			sys.stdout = cleaner_log
			print("\nDecompressing -- " + IMAGES_DIR + observer + "/" + YEAR + "/" + DAY + ".tar.lrz")

			os.chdir(IMAGES_DIR + observer + "/" + YEAR )
			subprocess.call(["lrzuntar", IMAGES_DIR + observer + "/" + YEAR + "/" + DAY + ".tar.lrz"])
			os.remove(IMAGES_DIR + observer + "/" + YEAR + "/" + DAY + ".tar.lrz")
		

		# CHECK IF ANY THE "DAY" DIRECTORY EXISTS
		if os.path.exists(IMAGES_DIR + observer + "/" + YEAR + "/" + DAY):

			# LOOP THROUGH ALL IMAGES FROM "DAY"
			for light_fits in listdir(IMAGES_DIR + observer + "/" + YEAR + "/" + DAY):

				# SAFE CHECK TO ENSURE WE ARE DEALING WITH .fits FILES AND NOT ANY CLEANED FILES
				if light_fits.endswith(".fits") and not light_fits.endswith("cor.fits"):

					# APPEND THE LIGHT FITS FILE TO A LIST
					light_list.append(IMAGES_DIR + observer + "/" + YEAR + "/" + DAY + "/" + light_fits)

					# CREATE cord_pending DIRECTORY IF IT DOES NOT ALREADY EXIST
					if not os.path.exists(IMAGES_DIR + observer + "/" + YEAR + "/" + DAY + "/cord_pending"):
						mkdir(IMAGES_DIR + observer + "/" + YEAR + "/" + DAY + "/cord_pending")

					# CREATE previous DIRECTORY IF IT DOES NOT ALREADY EXIST
					if not os.path.exists(IMAGES_DIR + observer + "/" + YEAR + "/" + DAY + "/previous"):
						mkdir(IMAGES_DIR + observer + "/" + YEAR + "/" + DAY + "/previous")

			# CHECK IF THERE ARE ANY IMAGES TO PROCESS
			if len(light_list) > 0:

				global CURRENT_OBSERVER, nprocesses

				# ONCE ALL IMAGES ARE ADDED TO THE LIST, USE MULTIPROCESSING AND CLEAN EACH IMAGE IN THE LIST
				timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")
				print (timestamp + "All files collected, starting multiprocessing pool to process images")

				CURRENT_OBSERVER = observer

				sys.stdout = terminal
				print(f"\nStarting COR for {observer}/{YEAR}/{DAY}")
				sys.stdout = cleaner_log
				print(f"\n{timestamp} Producing _COR for {observer}/{YEAR}/{DAY}")
				
				if (PROGRESS_BAR):

					if (MPROCESSING):
						#USE MULTIPROCESSING TO CLEAN IMAGES
						process_map(clean_fits, light_list, max_workers=nprocesses)
					elif (HTHREADING):
						#USE HYPERTHREADING TO CLEAN IMAGES
						thread_map(clean_fits, light_list, max_workers=nthreads)
					else:
						for image in light_list:
							clean_fits(image) #No Multiprocessing or Hyperthreading
				else:

					if (MPROCESSING):
						#USE MULTIPROCESSING TO CLEAN IMAGES
						pool = mp.Pool(processes=nprocesses)
					if (HTHREADING):
						#USE HYPERTHREADING TO CLEAN IMAGES
						pool = mp.pool.ThreadPool(processes=nthreads)
					if (MPROCESSING or HTHREADING):
						pool.map_async(clean_fits, light_list)
						pool.close()
						pool.join()
					else:
						for image in light_list:
							clean_fits(image) #No Multiprocessing or Hyperthreading


				timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")
				print (timestamp + "Cleaning complete")


def merge_darks_temp(left, right):
	# If the first list is empty, then nothing needs
    # to be merged, and you can return the second list as the result

    if len(left) == 0:

        return right


    # If the second list is empty, then nothing needs
    # to be merged, and you can return the first list as the result

    if len(right) == 0:

        return left


    result = []

    index_left = index_right = 0


    # Now go through both lists until all the elements
    # make it into the resultant list

    while len(result) < len(left) + len(right):

        # The elements need to be sorted to add them to the
        # resultant list, so you need to decide whether to get
        # the next element from the first or the second list

        if get_temp(left[index_left]) <= get_temp(right[index_right]):

            result.append(left[index_left])

            index_left += 1

        else:

            result.append(right[index_right])

            index_right += 1


        # If you reach the end of either list, then you can
        # add the remaining elements from the other list to
        # the result and break the loop

        if index_right == len(right):

            result += left[index_left:]

            break


        if index_left == len(left):

            result += right[index_right:]

            break


    return result


def merge_sort_darks(dark_list):

    # If the input array contains fewer than two elements,
    # then return it as the result of the function

    if len(dark_list) < 2:

        return dark_list


    midpoint = len(dark_list) // 2


    # Sort the list by recursively splitting the input
    # into two equal halves, sorting each half and merging them
    # together into the final result

    return merge_darks_temp(

        left=merge_sort_darks(dark_list[:midpoint]),

        right=merge_sort_darks(dark_list[midpoint:]))


def create_bins(dark_list, standard_width, standard_nbin):
	""" create_bins returns a dynamic distance bin partitioning list of a given dark list. 
        Width seperation is an ascending list of tuples, representing the intervals.
		Bin with insufficient darks for a 'standard' temperature width will become loose in temperature allowance,
		while bins with too many darks for a 'standard' temperature width will become tight in temperature allowance. 
    """ 

# Creating dynamic bins of all collected darks
#-------------------------------------------------------------------
	dark_bins = []
	next_bin = True 
	next_bin_start = 0
	dark_num = 0
	bin_num = 0

	while dark_num < len(dark_list):

		if next_bin:
			dark_bins.append([])
			dark_num = next_bin_start
			bin_lower_temp = get_temp(dark_list[dark_num])
			next_bin = False

		dark_temp = get_temp(dark_list[dark_num])

		if (dark_temp - bin_lower_temp) < standard_width:

			if len(dark_bins[bin_num]) < int(MAX_DARKS):
				dark_bins[bin_num].append(dark_list[dark_num])
			else:
				next_bin_start = dark_num
				bin_num += 1
				next_bin = True

		else:

			next_bin_start = dark_num

			if len(dark_bins[bin_num]) < int(MIN_DARKS):
				dark_bins[bin_num].append(dark_list[dark_num])
			else:
				bin_num += 1
				next_bin = True
		
		dark_num += 1
	
	# If last bin has less darks than MIN_DARKS, we can copy new darks from previous bin(s) 
	# until we have MIN_DARKS. Iterating backwards on the total dark list.
	
	last_bin = len(dark_bins) - 1

	if len(dark_bins) > 1:

		if len(dark_bins[last_bin]) < int(MIN_DARKS):

			for dark in reversed(dark_list):

				if len(dark_bins[last_bin]) < int(MIN_DARKS):

					if dark not in dark_bins[last_bin]:

						dark_bins[last_bin].append(dark)
				else:
					break

	#-------------------------------------------------------------------

	# Print the bin interval distibution 
	bin_intervals = []
	print(f"\nDARK_BINS ({len(dark_bins)}): {dark_bins}")

	for bin in dark_bins:

		if len(bin) != 0:
			bin_temp_low = get_temp(bin[0])
			bin_temp_high = get_temp( bin[len(bin)-1] )

		bin_intervals.append((round(bin_temp_low, 3), round(bin_temp_high, 3)))

	print(f'\nThe bin temperature interval distibution for {len(dark_bins)} bins is as follows in degrees kelvin, tuples: \n{bin_intervals}\n')

	return dark_bins


def get_darks_filter(raster_size, YEAR, DAY, incoming_DAY, incoming_YEAR, dark_span, light_temp):
	'''
		A wrapper filter function for the get_darks function. Prepares by enumerating all total darks
		as well as handles dark list with post-collection binning filter
	'''

	global temp_high, temp_low, standard_bin_width


	#-------------------------------
	# First numerate total darks so we can later stop dark collection early if no longer necessary to continue
	#-------------------------------
	total_darks = 0
	new_darks = []
	max_dark_span = int(MAX_DARK_SPAN)

	print("path for total darks: " + DARKS_DIR + raster_size)
	if os.path.exists(DARKS_DIR + raster_size):
		for dark_year in listdir(DARKS_DIR + raster_size):

			if int(YEAR) <= int(dark_year) <= int(incoming_YEAR):

				for dark_day in listdir(DARKS_DIR + raster_size + "/" + dark_year):

					if dark_day.endswith(".tar.lrz"):
						sys.stdout = terminal
						print(f"\nDecompressing -- dark_library/" + dark_day )
						sys.stdout = cleaner_log
						print("\nDecompressing -- " + DARKS_DIR + raster_size + "/" + dark_year + "/" + dark_day)
						os.chdir(DARKS_DIR + raster_size + "/" + dark_year )
						subprocess.call(["lrzuntar", DARKS_DIR + raster_size + "/" + dark_year + "/" + dark_day])
						os.remove(DARKS_DIR + raster_size + "/" + dark_year + "/" + dark_day)
						dark_day = dark_day.split(".tar.lrz")[0]

					if (YEAR == incoming_YEAR):

						if int(DAY)-max_dark_span <= int(dark_day) <= int(incoming_DAY)+max_dark_span:

							for dark in listdir(DARKS_DIR + raster_size + "/" + dark_year + "/" + dark_day):

								if dark.endswith(".fits"):

									total_darks += 1
							
							if int(DAY)+2 <= int(dark_day) <= int(incoming_DAY): 

								new_darks.append([YEAR, dark_day])

					else:

						if dark_year == YEAR:

							if int(dark_year) % 4 == 0:
								if int(DAY)-max_dark_span <= int(dark_day) <= 366: 
									for dark in listdir(DARKS_DIR + raster_size + "/" + dark_year + "/" + dark_day):
										if dark.endswith(".fits"):
											total_darks += 1

							else:
								if int(DAY)-max_dark_span <= int(dark_day) <= 365:
									for dark in listdir(DARKS_DIR + raster_size + "/" + dark_year + "/" + dark_day):
										if dark.endswith(".fits"):
											total_darks += 1
						
						else:

							if 0 <= int(dark_day) <= int(incoming_DAY):
								for dark in listdir(DARKS_DIR + raster_size + "/" + dark_year + "/" + dark_day):
									if dark.endswith(".fits"):
										total_darks += 1

	#---------------------------------------------------------------------

	filter_parameter = [raster_size, YEAR, DAY, dark_span, max_dark_span, total_darks, new_darks]
	
	dark_list_result = get_darks(filter_parameter)

	dark_list_bins = []

	if len(dark_list_result) >= int(MIN_DARKS):
		dark_list_result = merge_sort_darks(dark_list_result)

		print(f"\nDARK_LIST ({len(dark_list_result)}): {dark_list_result}\n" )

		lower_bound = get_temp(dark_list_result[0])
		print(f"\nBinning lower_bound: {lower_bound} Kelvin")
		upper_bound = get_temp(dark_list_result[len(dark_list_result)-1])
		print(f"Binning upper_bound: {upper_bound} Kelvin")

		binning_range = upper_bound - lower_bound

		standard_nbin = binning_range / standard_bin_width

		# LIMIT THE MAXIMUM ALLOWED TEMP SEPARATION BETWEEN THE BINS
		if standard_bin_width < 0.1:
			standard_bin_width = 0.1

		dark_list_bins = create_bins(dark_list_result, standard_bin_width, standard_nbin)
	
		return dark_list_bins, dark_list_result

	else:
		# Return dark list if insufficient in length anyway

		return dark_list_bins, dark_list_result

def get_new_darks(raster_size, YEAR, DAY, new_darks):
	""" Returns an appended dark list with darks that are considered new for the selected raster size of light.
		Typical dark collection searches for darks in span of one day ahead of light DOY 
		to (light DOY - configured max dark span limit). This accounts for any new detected darks 
		up to the current incoming DOY.
	"""
	global SAA_DARKS

	dark_list = []

	if len(new_darks) != 0:

		for dark in new_darks:

			YEAR = dark[0]
			DAY = dark[1]

			if os.path.exists(DARKS_DIR + raster_size + "/" + YEAR + "/" + DAY):
				if not os.path.exists(DARKS_DIR + raster_size + "/" + YEAR + "/" + DAY + "/PROCESSED_DARKS/"): 
					makedirs(DARKS_DIR + raster_size + "/" + YEAR + "/" + DAY + "/PROCESSED_DARKS/")		

			proc_dark_dir = DARKS_DIR + raster_size + "/" + YEAR + "/" + DAY + "/PROCESSED_DARKS/"
			if os.path.exists(proc_dark_dir):
				print("searching new darks in:" + proc_dark_dir)
				for processed_dark in os.listdir(proc_dark_dir):
					path = proc_dark_dir + processed_dark
					if get_size(path) == raster_size:
						dark_day = int(processed_dark[13:16])
						if int(DAY) == dark_day:
							# Before appending, check if dark is within defined SAA region
							if not is_in_SAA(path):
								dark_list.append(path)
								print(f'\nAdding new dark {path}')
							else:
								SAA_DARKS += 1
								print(f'\nDiscarded new dark {path} in SAA')

				if len(listdir(proc_dark_dir)) == 0:
					print("no new darks")
				else:
					print("found new darks!")


			unproc_dark_dir = DARKS_DIR + raster_size + "/" + YEAR + "/" + DAY
			if os.path.exists(unproc_dark_dir):
				print("\nsearching new darks in:" + unproc_dark_dir)
				for darks in os.listdir(unproc_dark_dir):
					if fnmatch.fnmatch(darks, "*.fits"):
						# Making sure not to repeat any already processed incoming darks
						if not os.path.exists(unproc_dark_dir + "/PROCESSED_DARKS/" + darks.split(".fits")[0] + "_cor.fits"):
							# Before appending, check if dark is within defined SAA region
							if not is_in_SAA(unproc_dark_dir + "/" + darks):
								dark_list.append(unproc_dark_dir + "/" + darks)
								print(f'\nAdding new dark {unproc_dark_dir + "/" + darks}')
							else:
								SAA_DARKS += 1
								print(f'\nDiscarded new dark {unproc_dark_dir + "/" + darks} in SAA')

				if len(listdir(unproc_dark_dir)) == 0:
					print("no new darks")
				else:
					print("found new darks!")

		
		return dark_list
		
	else:
		return dark_list


def get_darks(filter_parameter):
	""" This function will search through dates backwards by a Max Dark Span 
		and forwards in the exact dates of darks available,
		returning a list of all discovered dark frames compatible with the 
		raster size, exposure length, and SSA definition
	"""

	global CURRENT_OBSERVER, INCOMING_DARKS, SAA_DARKS
	
	raster_size = filter_parameter[0]
	YEAR = filter_parameter[1]
	DAY = filter_parameter[2]
	dark_span = filter_parameter[3]
	max_dark_span = filter_parameter[4]
	total_darks = filter_parameter[5]
	
	day = int(DAY)
	srchd_future = False

	print(f"starting get_darks({raster_size}, {YEAR}, {day}): \n")

	dark_list = []
	dark_list = get_new_darks(raster_size, YEAR, DAY, filter_parameter[6])

	print(f"Got {len(dark_list)} new darks; adding darks over span of {max_dark_span} days\n")

	
	while dark_span <= max_dark_span: 

		if day - 1 < 1:
		#Restarting from the previous year
			day = 365
			if (int(YEAR) - 1) % 4 == 0:
				day = 366
			year = str(int(YEAR) - 1)
		else:
			year = YEAR
		
		day_str = str(day)

		if day == int(DAY) - 1 and not srchd_future:
			day += 1
			day_str = str(day+1)
			srchd_future = True

		# Making sure it is in the correct DOY format for finding files 
		if len(day_str) == 1:
			day_str = "00" + day_str
		if len(day_str) == 2:
			day_str = "0" + day_str

		if os.path.exists(DARKS_DIR + raster_size + "/" + year + "/" + day_str):
			if not os.path.exists(DARKS_DIR + raster_size + "/" + year + "/" + day_str + "/PROCESSED_DARKS/"): 
				makedirs(DARKS_DIR + raster_size + "/" + year + "/" + day_str + "/PROCESSED_DARKS/")		

		proc_dark_dir = DARKS_DIR + raster_size + "/" + year + "/" + day_str + "/PROCESSED_DARKS/"
		print("searching in:" + proc_dark_dir) 
		if os.path.exists(proc_dark_dir):
			proc_dark_count = len(listdir(proc_dark_dir))
			if proc_dark_count == 0:
				print("no processed darks but dir exists")
			else:
				print(f"found {proc_dark_count} processed darks!")

			for processed_dark in os.listdir(proc_dark_dir):
				path = proc_dark_dir + processed_dark
				if get_size(path) == raster_size:
					dark_day = int(processed_dark[13:16])
					#print(f"evaluating {processed_dark} - Day ({day}) Dark_Day {dark_day}")
					if day == dark_day or day + 1 == dark_day:
					# Searching for a day ahead as well
					# Before appending, check if dark is within defined SAA region
						if not is_in_SAA(path):
							if path not in dark_list:
								dark_list.append(path)
								print(f'Adding dark {path}')
						else:
							SAA_DARKS += 1
							print(f'Discarded dark {path} in SAA')


		# If all total darks are collected, there is no need to continue the collection process
		if len(dark_list) == (total_darks - SAA_DARKS):
			return dark_list
		

		unproc_dark_dir = DARKS_DIR + raster_size + "/" + year + "/" + day_str
		print("searching in:" + unproc_dark_dir) 
		if os.path.exists(unproc_dark_dir):
			unproc_dark_count = len(listdir(unproc_dark_dir))
			if unproc_dark_count == 0:
				print("no unprocessed darks but dir exists")
			else:
				print(f"found {unproc_dark_count-1} unprocessed darks!")  # Subtract one to account for PROCESSED_DARKS directory

			for darks in os.listdir(unproc_dark_dir):
				if fnmatch.fnmatch(darks, "*.fits"):
					# Making sure not to repeat any already processed incoming darks
					if not os.path.exists(unproc_dark_dir + "/PROCESSED_DARKS/" + darks.split(".fits")[0] + "_cor.fits"):
						# Before appending, check if dark is within defined SAA region
						if not is_in_SAA(unproc_dark_dir + "/" + darks):
							if darks not in dark_list:
								dark_list.append(unproc_dark_dir + "/" + darks)
								print(f'\nAdding dark {unproc_dark_dir + "/" + darks}')
						else:
							SAA_DARKS += 1
							print(f'\nDiscarded dark {unproc_dark_dir + "/" + darks} in SAA')


		# If all total darks are collected, there is no need to continue the collection process
		if len(dark_list) == (total_darks - SAA_DARKS):
			return dark_list

		dark_span += 1
		day -= 1

	# If we are at the end of the while loop, it means all possible searches were conducted
	return dark_list

	
	
def produce_cord(YEAR, DAY, incoming_YEAR):
	print(f"\nstarting produce_cord({YEAR}, {DAY}): \n")
	timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")

	# LOOP THROUGH EACH DAY THAT MAY NEED TO REPROCESS LIGHTS
	proc_start = int(DAY) - int(NODARKS_SPAN)
	proc_end = int(DAY) + 1

	# LOOP THROUGH EACH OF THE OBSERVERS
	for observer in listdir(IMAGES_DIR):

		for day in range(proc_start, proc_end):

			global MASTER_DARK, MASTERDARK_DICTIONARIES, CLEAN_TYPE, CURRENT_OBSERVER, DARKS, nprocesses
			
			MASTER_DARK = None
			dark_list = []
			light_list = []

			if day <= 0: 
				year = str(int(YEAR) - 1)
				if int(YEAR) % 4 == 0:
					if int(YEAR) % 100 == 0:
						if int(YEAR) % 400 == 0:
							day = 366 + day
						else:
							day = 365 + day
					else:
						day = 365 + day
				else:
					day = 365 + day
			else:
				year = YEAR
			
			# Formating DOY to 000 if below 100
			day_str = str(day)
			if len(day_str) == 1:
				day_str = "00" + day_str
			elif len(day_str) == 2:
				day_str = "0" + day_str

			if os.path.isfile(IMAGES_DIR + observer + "/" + year + "/" + day_str + ".tar.lrz"):

				sys.stdout = terminal
				print(f"\nDecompressing -- images/" + observer + "/" + year + "/" + day_str + ".tar.lrz")
				sys.stdout = cleaner_log
				print("\nDecompressing -- " + IMAGES_DIR + observer + "/" + year + "/" + day_str + ".tar.lrz")

				subprocess.call(["lrzuntar", IMAGES_DIR + observer + "/" + year + "/" + day_str + ".tar.lrz"])
				os.remove(IMAGES_DIR + observer + "/" + year + "/" + day_str + ".tar.lrz")

			timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")
			print (timestamp + "looking in " + IMAGES_DIR + observer + "/" + year + "/" + day_str)
			if os.path.exists(IMAGES_DIR + observer + "/" + year + "/" + day_str + "/cord_pending"):
				if len(listdir(IMAGES_DIR + observer + "/" + year + "/" + day_str + "/cord_pending")) != 0:
					print(f"Dark subtracting the images in {year}/{day}/cord_pending")

					# ADD ALL THE IMAGES THAT NEED TO BE CLEANED TO A LIST
					for light_fits in listdir(IMAGES_DIR + observer + "/" + year + "/" + day_str + "/cord_pending"):
						if light_fits.endswith("_cor.fits"):

							# OPEN FITS FILE AND READ HEADER
							light_file = IMAGES_DIR + observer + "/" + year + "/" + day_str + "/cord_pending/" + light_fits

							size = get_size(light_file)
							temp = get_temp(light_file)

							timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")
							print(f"\n{timestamp} CORD_PENDING {light_file} size:", size)
							print(f'Is light in SAA: {is_in_SAA(light_file)}')

							if size not in MASTERDARK_DICTIONARIES:

								dark_list_bins, dark_list = get_darks_filter(size, year, day_str, DAY, incoming_YEAR, 1, temp)

								if len(dark_list) >= int(MIN_DARKS):
									timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")
									print (timestamp + "Master dark is not created, we will create one for " + size)

									sys.stdout = terminal
									print(f"\nProcessing DARKS for {observer}/{year}/{day}")
									sys.stdout = cleaner_log

									bin_temps = []
									MASTERDARK_DICTIONARIES[size] = []											

									for bin in dark_list_bins:
										
										if len(bin) >= int(MIN_DARKS):
											for dark in bin:
												bin_temps.append(get_temp(dark))
										
											sys.stdout = terminal
											print(f" -- {size}, Masterdark for bin #{dark_list_bins.index(bin) + 1}")
											sys.stdout = cleaner_log
											MASTER_DARK = get_master_dark(bin)
											median_temp = median(bin_temps)
											bin_temps.clear()

											MASTERDARK_DICTIONARIES[size].append([median_temp, MASTER_DARK])
										else:
											MASTERDARK_DICTIONARIES[size].append([])
									
									print(f"\nMASTERDARK_DICTIONARIES: \n{MASTERDARK_DICTIONARIES}")

									DARKS[size] = dark_list_bins
									light_list.append(IMAGES_DIR + observer + "/" + year + "/" + day_str + "/cord_pending/" + light_fits)
								else:
									if len(dark_list) < int(MIN_DARKS):
										print("\nThe bare minimum dark(s) required is at least " + MIN_DARKS + ". Insufficient darks detected. Could not create cord file.\n")

									# DIRECTORY TO PASS TO NODARKS LIGHTLIST
									nodarks_light = IMAGES_DIR + observer + "/" + year + "/" + day_str + "/previous/" + light_fits

									# WRITE TO NODARKS LIGHTLIST IF NOT ENOUGH DARKS TO PRODUCE CORD
									nodarks_list = open(DARKS_DIR + "lightlists/nodarks.lightlist", "a+")
									nodarks_list.write(f"{nodarks_light}\t{size}\n")
									nodarks_list.close()

									# REMOVE LIGHT WITH NO DARKS FROM CORD_PENDING
									os.remove(light_file)
							else:
								timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")
								print (timestamp + "Binned master darks already exists for " + size)
								light_list.append(IMAGES_DIR + observer + "/" + year + "/" + day_str + "/cord_pending/" + light_fits)


					CLEAN_TYPE = 1
					CURRENT_OBSERVER = observer
					if len(light_list) > 0:

						sys.stdout = terminal
						print(f"\nProducing _CORD for {observer}/{year}/{day}")
						sys.stdout = cleaner_log

						if (PROGRESS_BAR):

							if (MPROCESSING):
								#USE MULTIPROCESSING TO CLEAN IMAGES
								process_map(clean_fits, light_list, max_workers=nprocesses)
							elif (HTHREADING):
								#USE HYPERTHREADING TO CLEAN IMAGES
								thread_map(clean_fits, light_list, max_workers=nthreads)
							else:
								for image in light_list:
									clean_fits(image) #No Multiprocessing or Hyperthreading
						else:

							if (MPROCESSING):
								#USE MULTIPROCESSING TO CLEAN IMAGES
								pool = mp.Pool(processes=nprocesses)
							if (HTHREADING):
								#USE HYPERTHREADING TO CLEAN IMAGES
								pool = mp.pool.ThreadPool(processes=nthreads)
							if (MPROCESSING or HTHREADING):
								pool.map_async(clean_fits, light_list)
								pool.close()
								pool.join()
							else:
								for image in light_list:
									clean_fits(image) #No Multiprocessing or Hyperthreading
				else:
					for light_fits in listdir(IMAGES_DIR + observer + "/" + year + "/" + day_str):
						if light_fits.endswith("_cor.fits"):
						
							# OPEN FITS FILE AND READ HEADER
							light_file = IMAGES_DIR + observer + "/" + year + "/" + day_str + "/" + light_fits

							size = get_size(light_file)

							# DIRECTORY TO PASS TO NODARKS LIGHTLIST
							nodarks_light = IMAGES_DIR + observer + "/" + year + "/" + day_str + "/previous/" + light_fits

							# WRITE TO NODARKS LIGHTLIST IF NOT ENOUGH DARKS TO PRODUCE CORD
							nodarks_list = open(DARKS_DIR + "lightlists/nodarks.lightlist", "a+")
							nodarks_list.write(f"{nodarks_light}\t{size}\n")
							nodarks_list.close()

							# REMOVE LIGHT WITH NO DARKS FROM CORD_PENDING
							os.remove(light_file)


def queue_cord(incoming_darks, incoming_lights, images_dir):
	'''
		Queues up all "_cor.fits" files in preparation for "_cord.fits" file production
		Reads and writes from relevant "light-lists" for "_cord.fits" re-processing
	'''

	global DAY_RANGE_LIMIT

	# QUEUE UP ALL INCOMING LIGHTS FOR POTENTIAL CORD PRODUCTION
	for light in incoming_lights:

		light_YEAR = light.split("NEOS_SCI_")[1][0:4]
		light_DAY = light.split("NEOS_SCI_")[1][4:7]
		observer = (light.split("images/"))[1].split("/")[0]
		light_name = (light.split(f'{light_DAY}/')[1]).split('.fits')[0] + '_cor.fits'

		if not os.path.exists(f'{images_dir}{observer}/{light_YEAR}/{light_DAY}/cord_pending/{light_name}'):

			if os.path.exists(f'{images_dir}{observer}/{light_YEAR}/{light_DAY}/{light_name}'):

				shutil.move(f'{images_dir}{observer}/{light_YEAR}/{light_DAY}/{light_name}', 
							f'{images_dir}{observer}/{light_YEAR}/{light_DAY}/cord_pending')
			else:

				if os.path.exists(f'{images_dir}{observer}/{light_YEAR}/{light_DAY}/previous/{light_name}'):

					shutil.move(f'{images_dir}{observer}/{light_YEAR}/{light_DAY}/previous/{light_name}', 
								f'{images_dir}{observer}/{light_YEAR}/{light_DAY}/cord_pending')
		

	if len(incoming_darks) != 0:
		sizes = []
		all_lights = []
		pending_lights = {}

		for dark in incoming_darks:
			size_dir = dark.split("dark_library/")[1]
			size = size_dir.split("/")[0]
			dark_YEAR = dark.split("NEOS_SCI_")[1][0:4]
			dark_DAY = dark.split("NEOS_SCI_")[1][4:7]

			if not size in sizes:
				sizes.append(size)

				if os.path.exists(f"{DARKS_DIR}lightlists/{size}.lightlist"):

					# QUEUE UP LIGHTS FOR PROCESSING FROM RELEVANT LIGHTLISTS
					with open(f"{DARKS_DIR}lightlists/{size}.lightlist", "r+") as lightlist:
						for line in lightlist:
							light = line.split()
							if not light[0] in all_lights:
								all_lights.append(light[0])
						lightlist.seek(0)

						for light in all_lights:
							
							lightlist.write(f'{light}\n')

							light_YEAR = light.split("NEOS_SCI_")[1][0:4]
							light_DAY = light.split("NEOS_SCI_")[1][4:7]
							observer = (light.split("images/"))[1].split("/")[0]
							light_name = light.split(f'{light_DAY}/')[1]

							if (light_YEAR == dark_YEAR) and ((int(dark_DAY) - int(light_DAY)) <= int(DAY_RANGE_LIMIT)):

								if not os.path.exists(f'{images_dir}{observer}/{light_YEAR}/{light_DAY}/cord_pending/{light_name}'):

									if os.path.exists(light):
										shutil.move(light, f'{images_dir}{observer}/{light_YEAR}/{light_DAY}/cord_pending')

									elif os.path.exists(f'{images_dir}{observer}/{light_YEAR}/{light_DAY}/previous/{light_name}'):
										shutil.move(f'{images_dir}{observer}/{light_YEAR}/{light_DAY}/previous/{light_name}', 
													f'{images_dir}{observer}/{light_YEAR}/{light_DAY}/cord_pending')

						lightlist.truncate()
					lightlist.close()

				# QUEUE UP LIGHTS FOR PROCESSING FROM PREVIOUSLY UNFINISHED LIGHTS
				if os.path.exists(f"{DARKS_DIR}lightlists/nodarks.lightlist"):
					with open(f"{DARKS_DIR}lightlists/nodarks.lightlist", "r+") as nodarks:
						for line in nodarks:
							(key, val) = line.split()
							if val == size:
								light_path = key.split("previous/")[0]
								light_name = key.split("previous/")[1]
								light_YEAR = light_name.split("NEOS_SCI")[1][0:4]
								light_DAY = light_name.split("NEOS_SCI_")[1][4:7]

								if (light_YEAR == dark_YEAR) and ((int(dark_DAY) - int(light_DAY)) <= int(NODARKS_SPAN)):

									if  not os.path.exists(light_path + "cord_pending/" + light_name):	

										if os.path.exists(key):
											shutil.move(key,  light_path + "cord_pending/" + light_name)
							else:
								pending_lights[key] = val
						if pending_lights == {}:
							os.remove(f"{DARKS_DIR}lightlists/nodarks.lightlist")
						else:
							nodarks.seek(0)
							for light in pending_lights.keys():
								nodarks.write(f"{light}\t{pending_lights.get(light)}\n")
							nodarks.truncate()
					nodarks.close()



def main():
	base_dir = str(sys.argv[2])
	global SAVE_DIR, IMAGES_DIR, DARKS_DIR, LOG_DIR, INCOMING, INCOMING_LIGHTS, INCOMING_DARKS, MPROCESSING, HTHREADING, cleaner_log

	SAVE_DIR = base_dir + SAVE_DIR
	IMAGES_DIR = base_dir + IMAGES_DIR
	DARKS_DIR = base_dir + DARKS_DIR
	LOG_DIR = base_dir + LOG_DIR
	INCOMING = base_dir + INCOMING
	fin_processing = True

	date = datetime.datetime.utcnow()
	timestamp = date.strftime("%Y-%j-%H:%M:%S - ")

	# CLEANER LOG FILE 
	FILE_NAME = LOG_DIR + "NEOS_" + date.strftime("%Y%j%H%M%S") + ".CLEANLOG"
	cleaner_log = open(FILE_NAME, "w+")
	sys.stdout = cleaner_log 

	# VERIFY INPUT FILE IS AVAILABLE
	
	if not os.path.exists(INCOMING):
		print (timestamp + "Run n1sort.py first to inform cleaner which files to process or ensure that files are provided")
		raise Exception('n1clean: No incoming files; ending...')

	# CREATE DICTIONARY FROM INCOMING CLEANLIST
	cleanlist = {}
	with open(INCOMING) as current:
		for line in current:
			(key, val) = line.split()
			cleanlist[key] = val
	current.close()

	#CREATE YEAR, DOY, INCOMING_LIGHTS, INCOMING_DARKS BASED ON INCOMING CLEANLIST
	YEAR_DAY = []
	for incoming in cleanlist.keys():
		year_day = incoming.split("NEOS_SCI_")[1][0:7]
		if not year_day in YEAR_DAY:
			YEAR_DAY.append(year_day)
		if cleanlist.get(incoming) == 'LIGHT':
			INCOMING_LIGHTS.append(incoming)
		else:
			INCOMING_DARKS.append(incoming)

	# ADD INITIAL DATA TO THE LOG FILE
	print (f";FILENAME:\t\t\t {FILE_NAME}")
	print (";FILE_CREATION_TIME:\t\t " + timestamp)
	print (";**************************************************************************")


	# MAKE SURE APPROPRIATE OUTGOING DIRECTORIES ARE CREATED
	for observer in listdir(IMAGES_DIR):
		if not os.path.exists(SAVE_DIR + observer) and not observer == "previous":
			mkdir(SAVE_DIR + observer)

	timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")
	print (timestamp + "Detected new incoming cleanlist; beginning processing\reprocessing \"cor\" and \"cord\"")

	# FINDING THE MOST RECENT INCOMING YEAR IN INCOMING LIST
	#-------------------------------------------------------
	incoming_YEAR = "0"

	for year_day in YEAR_DAY:
		if int(year_day[0:4]) > int(incoming_YEAR):
			incoming_YEAR = year_day[0:4]
	#-------------------------------------------------------

	# TERMINAL EXECUTION INFO
	sys.stdout = terminal					
	print("\nProducing _COR file(s), processing all DARKS, and producing _CORD file(s)...")
	sys.stdout = cleaner_log
	
	# PROCESSING MAIN -- COR PRODUCTION
	print (timestamp + "Processing any \"cor\" first")
	for year_day in YEAR_DAY:
		YEAR = year_day[0:4]
		DAY = year_day[4:7]
		produce_cor(YEAR, DAY)

	# QUEUING CORS FOR CORD PROCESSING
	timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")
	print (timestamp + "Queuing up \"cor\", in preparation for \"cord\" production")
	queue_cord(INCOMING_DARKS, INCOMING_LIGHTS, IMAGES_DIR)

	# PROCESSING MAIN -- CORD PRODUCTION
	timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")
	print (timestamp + "Done processing \"cor\", moving on to \"cord\"")
	for year_day in YEAR_DAY:
		YEAR = year_day[0:4]
		DAY = year_day[4:7]
		produce_cord(YEAR, DAY, incoming_YEAR)

		for observer in listdir(IMAGES_DIR):
			if os.path.exists(f'{IMAGES_DIR}{observer}/{YEAR}/{DAY}/cord_pending'):
				if not len(os.listdir(f'{IMAGES_DIR}{observer}/{YEAR}/{DAY}/cord_pending')) == 0:
					fin_processing = False

			
	timestamp = datetime.datetime.utcnow().strftime("%Y-%j-%H:%M:%S - ")
	print (timestamp + "Completed cleaning")

	if fin_processing:
		timestamp = datetime.datetime.utcnow().strftime("%Y%j%H%M%S")
		os.rename(INCOMING, base_dir + "incoming/" + timestamp + ".cleanlist")
		cleaner_log.close()

if __name__ == "__main__":
	main()
