import numpy as np
from matplotlib import pyplot as plt
from scipy import signal
from scipy import optimize as opt

from astropy import units as u
from astropy.coordinates import SkyCoord

from ctapipe.image import extractor
from ctapipe import utils
from ctapipe.visualization import CameraDisplay
from ctapipe.io import EventSource, DataWriter
from ctapipe.visualization import ArrayDisplay, CameraDisplay
from ctapipe.calib import CameraCalibrator
from ctapipe.calib.camera.gainselection import ThresholdGainSelector
from ctapipe.calib.camera.flatfield import FlasherFlatFieldCalculator
from ctapipe.calib.camera.pedestals import PedestalIntegrator

#Sprime = S*conv2 =S*np.sqrt(conv1)

#np.sqrt(conv1*(Varn(s)**2) + (S**2/conv1)*conv1var**2)

source = EventSource("/mnt/c/Users/ctoen/Documents/processed.h5", max_events=100000)

for i,event in enumerate(source):

  print(event)

