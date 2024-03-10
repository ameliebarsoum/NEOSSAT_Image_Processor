# n1publish.py
#
# REVISION HISTORY:
#  3.50 - M. Hamza      - 2019 Aug 19 - Initial release to publisher program that moves files from image/cleaner/outgoing to image/outgoing for FTP transfer 
#  4.00 - A. Fagarasanu - 2021 Jan 25 - Implemented optional fits file compression (fits.gz) before files are transferred.

from os import listdir, path, remove
from shutil import move
import sys
import gzip
import multiprocessing as mp
from multiprocessing import Pool
import configparser
config = configparser.ConfigParser()
config.read(path.join(path.abspath(path.dirname(__file__)), '../fits_cleaner/config/n1clean.config'))

# DIRECTORIES
FROM = "fits_cleaner/outgoing/"
TO   = "outgoing/" 

# CONFIGURATION
COMPRESS = config.getboolean('fitsconfig', 'COMPRESS')
nprocesses = mp.cpu_count() # Gives number of processing cores available

def compress_files(file_dir):
	input = open(file_dir, 'rb')
	s = input.read()
	input.close()

	output = gzip.GzipFile(file_dir + ".gz", 'wb')
	output.write(s)
	output.close()

	remove(file_dir)

def publish_files():
	global FROM, TO
	
	for observer in listdir(FROM):
		for fits in listdir(FROM + observer):

			if (fits.endswith(".fits") or fits.endswith(".fits.gz")):

				move(FROM + observer + "/" + fits, TO + observer + "/" + fits)
	
def main():
	base_dir = str(sys.argv[2])

	global FROM, TO

	FROM = base_dir + FROM
	TO   = base_dir + TO

	if COMPRESS:
		files = []

		for observer in listdir(FROM):
			for fits in listdir(FROM + observer):
				if (fits.endswith(".fits")):
					files.append(FROM + observer + "/" + fits)

		if len(files) != 0:
			print("\nCompressing files...")
		else:
			print("\nn1publish: No incoming files; ending...")

		pool = mp.Pool(processes=nprocesses)
		pool.map(compress_files, files)

	publish_files()

if __name__ == "__main__":
	main()


