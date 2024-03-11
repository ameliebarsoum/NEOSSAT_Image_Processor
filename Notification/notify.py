import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
import numpy as np
import cv2
from scipy.ndimage import shift
from skimage.transform import resize
import random

"""
Algorithm to send email notifications for each detected object that cannot be classified. 
Uses smtplib for Python: https://docs.python.org/3/library/smtplib.html
- Receives FITS file with unclassified object, and unclassified object's coordinates (should be 
contained in flagges folder)
- Returns email containing FITS file, time image was taken, name of FITS file, RA/Dec,
coordinates of detected anomaly (and maybe nearby known objects?)
"""

INPUT_PATH = 'mission_images/'
ALIGNED_PATH = 'aligned/'
GAUSSIANBLUR_PATH = 'gaussian_blur/'
DIFFERENCED_PATH = 'differenced/'
FLAGGED_PATH = 'flagged/'