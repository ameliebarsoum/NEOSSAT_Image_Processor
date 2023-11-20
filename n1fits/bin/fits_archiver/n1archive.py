#------------------------------------------------------------------------------
#
# Name       : n1archive.py
# Written by : Adrian Fagarasanu
# Date       : February 15, 2021
# Description
#
#  This script is designed to perform scheduled intensive archiving/compression 
#  of local library related files that, as pipeline utilization increases, 
#  accumulate to large sizes.
#
#------------------------------------------------------------------------------
import os
from os import listdir, path, remove, chdir
import shutil
import subprocess
import time
import datetime
import configparser
config = configparser.ConfigParser()
config.read(path.join(path.abspath(path.dirname(__file__)), '../fits_cleaner/config/n1clean.config'))

#|========================|
#|-----CONFIGURATION------|
#|========================|

# Change to path of device operating from
DATA_PATH = config['fitsconfig']['DATA_PATH']

# Java convertor directory paths
java_basedir = DATA_PATH + "image/fits_processor/"
gps_dir = java_basedir + "data/gpsextract/"
incoming_prev_dir = java_basedir + "incoming/previous/"

# Stats directory paths
stats_basedir = DATA_PATH + "image/fits_stats/"
stats_comp_dir = stats_basedir + "work/Completed/"

# Cleaning file directory
clean_basedir = DATA_PATH + "image/fits_cleaner/"
clean_darklib = clean_basedir + "data/dark_library/"
clean_imglib = clean_basedir + "data/images/"

# Compression Type
#-------------------------------
# Options affecting compression:
#  -b            bzip2 compression
#  -g            gzip compression using zlib
#  -l            lzo compression (ultra fast)
#  -n            no backend compression - prepare for other compressor
#  -z            zpaq compression (best, extreme compression, extremely slow)

compType = "-z"

# Time (in days) which determines how old files have to be before they are archivable
archiveTime = 30

#|==========================|
#|==========================|

def archive_files(YEAR, Doy):

    if (int(Doy) - archiveTime) > 0:

        year = YEAR
        start_doy = int(Doy) - archiveTime

    else:

        year = str(int(YEAR) - 1)
        start_doy = int(Doy) - archiveTime

        if int(year) % 4 == 0:
            if int(year) % 100 == 0:
                if int(year) % 400 == 0:
                    start_doy = 366 + start_doy
                else:
                    start_doy = 365 + start_doy
            else:
                start_doy = 365 + start_doy
        else:
            start_doy = 365 + start_doy

    # Formating DOY to 000 if below 100
    start_doy = str(start_doy)
    if len(start_doy) == 1:
        start_doy = "00" + start_doy
    elif len(start_doy) == 2:
        start_doy = "0" + start_doy

    # Creating a concat. date
    ref_date = year + start_doy

    #--------------------------
    #--FITPROCESSOR ARCHIVING--
    #--------------------------
    print("\nArchiving fits_processor files...")

    for yr in listdir(gps_dir):

        os.chdir(gps_dir + yr)

        for file in listdir(gps_dir + yr):
            
            if file.endswith(".gpsextract"):
                file_date = file.split("NEOS_")[1].split(".gpsextract")[0][0:7]
                
                if int(file_date) <= int(ref_date):
                    subprocess.call(["lrzip", compType, file])
                    os.remove(file)

    os.chdir(incoming_prev_dir)

    for file in listdir(incoming_prev_dir):
        if file.endswith(".DC1") or file.endswith(".VC1"):
            file_date = file.split("NEOS_")[1].split(".")[0][0:7]
                
            if int(file_date) <= int(ref_date):
                subprocess.call(["lrzip", compType, file])
                os.remove(file)

    #------------------------
    #--FITS STATS ARCHIVING--
    #------------------------
    print("\nArchiving stats files...")

    if os.path.exists(stats_comp_dir):
        os.chdir(stats_comp_dir)

        for file in listdir(stats_comp_dir):

            if file.endswith(".stats"):
                fileYear = file.split("NEOS_")[1].split("-")[0]
                fileDoy = file.split("-")[1].split("_")[0]
                file_date = fileYear + fileDoy

                if int(file_date) <= int(ref_date):
                    subprocess.call(["lrzip", compType, file])
                    os.remove(file)

    #-------------------------------
    #--FITS DARK LIBRARY ARCHIVING--
    #-------------------------------
    print("\nArchiving dark library files...")

    for size in listdir(clean_darklib):
        if size != "lightlists":
            for yr in listdir(clean_darklib + size):
                if int(yr) <= int(year):

                    os.chdir(clean_darklib + size + "/" + yr)

                    for doy in listdir(clean_darklib + size + "/" + yr):

                        if os.path.isdir(doy):
                            if int(yr) == int(year):
                                if int(doy) <= int(Doy):
                                    subprocess.call(["lrztar", compType, doy])
                                    shutil.rmtree(doy, ignore_errors=True)
                            else:
                                subprocess.call(["lrztar", compType, doy])
                                shutil.rmtree(doy, ignore_errors=True)


    #--------------------------------
    #--FITS IMAGE LIBRARY ARCHIVING--
    #--------------------------------
    print("\nArchiving image library files...")
    
    for observer in listdir(clean_imglib):
        for yr in listdir(clean_imglib + observer):
            if int(yr) <= int(year):

                os.chdir(clean_imglib + observer + "/" + yr)

                for doy in listdir(clean_imglib + observer + "/" + yr):

                    if os.path.isdir(doy):
                        if int(yr) == int(year):
                            if int(doy) <= int(Doy):
                                subprocess.call(["lrztar", compType, doy])
                                shutil.rmtree(doy, ignore_errors=True)
                        else:
                            subprocess.call(["lrztar", compType, doy])
                            shutil.rmtree(doy, ignore_errors=True)

def main():

    # Isolate the Year and Doy at time of run, so that we may distingush which files/directories to archive
    YEAR = datetime.datetime.utcnow().strftime("%Y")
    Doy = datetime.datetime.utcnow().strftime("%j")

    print(f"Archiving all relevant uncompressed files 2 weeks from {YEAR}/{Doy} and older...")

    archive_files(YEAR, Doy)


if __name__ == "__main__":
	main()

