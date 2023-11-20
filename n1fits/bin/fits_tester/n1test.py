# n1test.py
#
# REVISION HISTORY:
#  4.00 - A. Fagarasanu - 2021 Mar - Initial version
import os
from os import listdir, makedirs, path
import subprocess
from shutil import copy, move, copyfile
from ftplib import FTP

from subprocess import Popen, call

import configparser
config = configparser.ConfigParser()
config.read(path.join(path.abspath(path.dirname(__file__)), '../fits_cleaner/config/n1clean.config'))

################################
#######  CONFIGURATION  ########
################################

BIN_PATH = config['fitsconfig']['BIN_PATH']
DATA_PATH = config['fitsconfig']['DATA_PATH']

#Java convertor directory
javadir = BIN_PATH + "fits_processor/JneossatTLM.jar"
java_basedir = DATA_PATH + "image/fits_processor/"

#Java data paths
unfinished_tlm = java_basedir + "data/unfinished_tlm/"
published_lists = java_basedir + "data/published_lists/"
incoming = java_basedir + "incoming/"
 
MODE = "0"
SOURCE = None
DESTINATION = DATA_PATH + "image/fits_processor/incoming/"
FTP_DESTINATION = DATA_PATH + "image/fits_processor/outgoing/ASTRO/"

def iterate_FTP_range(YEAR, DOY, direction):

    if direction == "forward":
        DOY = str(int(DOY) + 1)

        if len(DOY) == 1:
            DOY = "00" + DOY
        elif len(DOY) == 2:
            DOY = "0" + DOY
                
        if (int(YEAR) % 400 == 0):
            if int(DOY) > 366:
                DOY = "001"
                YEAR = str(int(YEAR) + 1)
        else:    
            if int(DOY) > 365:
                DOY = "001"
                YEAR = str(int(YEAR) + 1)  

    elif direction == "backward":
        DOY = str(int(DOY) - 1)

        if len(DOY) == 1:
            DOY = "00" + DOY
        elif len(DOY) == 2:
            DOY = "0" + DOY
                
        if (int(YEAR) % 400 == 0):
            if int(DOY) < 1:
                DOY = "366"
                YEAR = str(int(YEAR) - 1)
        else:    
            if int(DOY) < 1:
                DOY = "365"
                YEAR = str(int(YEAR) - 1) 

    return YEAR, DOY  

def main():

    global SOURCE, MODE

    while not MODE in ["1", "1 ", "2", "2 "]:
        MODE = input("Select desired tester mode:\n1 -- FITS Processor Incoming Source\n2 -- FTP Image DOY Range\n\n")
        if not MODE in ["1", "1 ", "2", "2 "]:
            print("\nInvalid input.\n")

    if MODE in ["1", "1 "]:

        SOURCE = input("\nEnter your source path: \n")

        SSCLOGs = []
        NAVSOLs = []
        DC1s_VC1s = []

        for file in os.listdir(SOURCE):

            if file.endswith(".SSCLOG"):
                SSCLOGs.append(file)

            elif file.endswith(".navsol"):
                NAVSOLs.append(file)

            elif file.endswith(".DC1") or file.endswith(".VC1"):
                DC1s_VC1s.append(file)


        for SSCLOG in SSCLOGs:
            copy(SOURCE + SSCLOG, DESTINATION)
        
        for navsol in NAVSOLs:
            copy(SOURCE + navsol, DESTINATION)
        
        DC1s_VC1s.sort()

        #------------------------------------------------------------------
        # Clearing all published lists and unfinished tlm before beginning
        #------------------------------------------------------------------

        pubFiles = []

        for file in DC1s_VC1s:
            pubFiles.append(file + ".LIST")
           
        if os.path.exists(unfinished_tlm):
            if len(listdir(unfinished_tlm)) != 0:
                for file in listdir(unfinished_tlm):
                    os.remove(unfinished_tlm + file)

        if os.path.exists(published_lists):
            if len(listdir(published_lists)) != 0:
                for file in listdir(published_lists):
                    if file in pubFiles:
                        os.remove(published_lists + file)
        #------------------------------------------------------------------

        for file in DC1s_VC1s:
            copy(SOURCE + file, DESTINATION)
            exit_code_pipeline = call("python3 " + BIN_PATH + "launch_fits.py", shell=True)

        
        print("All test files transferred and processed.")

    if MODE in ["2", "2 "]:

        DWN_RANGE_valid = False

        while not DWN_RANGE_valid:
            DWN_RANGE = input("\nEnter desired download range in format: YEAR/DOY-YEAR/DOY\n")
            try:
                if (len(DWN_RANGE) == 17):
                    YEAR_start = DWN_RANGE.split("-")[0].split("/")[0]
                    DOY_start = DWN_RANGE.split("-")[0].split("/")[1]

                    YEAR_end = DWN_RANGE.split("-")[1].split("/")[0]
                    DOY_end = DWN_RANGE.split("-")[1].split("/")[1]

                    if len(YEAR_start) == 4 and len(YEAR_end) == 4 and len(DOY_start) == 3 and len(DOY_end) == 3:
                        if int(DOY_start) <= 365 and int(DOY_end) <= 365:
                            if int(YEAR_start) >= 2017 and int(YEAR_end) >= 2017:
                                DWN_RANGE_valid = True
                                break
                            else:
                                print("\nNo data available before 2017.")
                print("\nInput or formatting is invalid. Example entry: 2019/067-2019/070")
            
            except Exception as e:
                print("\nInput or formatting is invalid. Example entry: 2019/067-2019/070")


        ftp = FTP('ftp.asc-csa.gc.ca')
        ftp.login()

        while True:

            print(f'~/users/OpenData_DonneesOuvertes/pub/NEOSSAT/ASTRO/{YEAR_start}/{DOY_start}')
            
            ftp.cwd(f'~/users/OpenData_DonneesOuvertes/pub/NEOSSAT/ASTRO/{YEAR_start}/{DOY_start}')
            
            filenames = ftp.nlst() # get filenames within the directory

            for filename in filenames:

                print(f"\nDownloading {filename}...")
                local_filename = os.path.join(FTP_DESTINATION, filename)
                file = open(local_filename, 'wb')
                ftp.retrbinary('RETR '+ filename, file.write)

                file.close()
            
            YEAR_start, DOY_start = iterate_FTP_range(YEAR_start, DOY_start, "forward")

            YEAR_prev, DOY_prev = iterate_FTP_range(YEAR_start, DOY_start, "backward")

            if (YEAR_prev == YEAR_end) and (DOY_prev == DOY_end):
                break


        ftp.quit()
        exit_code_pipeline = call("python3 " + BIN_PATH + "launch_fits.py", shell=True)



if __name__ == "__main__":
	main()