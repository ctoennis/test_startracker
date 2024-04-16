import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

from ctapipe.image import extractor
from ctapipe import utils
from ctapipe.io import EventSource, DataWriter
from ctapipe_io_zfits import ProtozfitsDL0TelescopeEventSource

source = ProtozfitsDL0TelescopeEventSource("/mnt/c/Users/ctoen/Documents/TEL001_SDH3001_20231015T033032_SBID0000000002000000039_OBSID0000000002000000107_CHUNK001.fits.fz", max_events=100000)

#So first i need to check the status of each pixel 
