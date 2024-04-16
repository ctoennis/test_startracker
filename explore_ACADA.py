import numpy as np
from matplotlib import pyplot as plt
from scipy import signal
from scipy import optimize as opt

from astropy import units as u
from astropy.coordinates import SkyCoord

from lstchain.calib.camera.calibration_calculator import LSTCalibrationCalculator
from traitlets.config import Config
from ctapipe.image import extractor
from ctapipe import utils
from ctapipe.visualization import CameraDisplay
from ctapipe.io import EventSource
from ctapipe_io_zfits import ProtozfitsDL0TelescopeEventSource
from ctapipe.visualization import ArrayDisplay, CameraDisplay
from ctapipe.calib import CameraCalibrator
from ctapipe.calib.camera.gainselection import ThresholdGainSelector
from ctapipe.calib.camera.flatfield import FlasherFlatFieldCalculator
from ctapipe.calib.camera.pedestals import PedestalIntegrator

import eventio

import copy

#Sprime = S*conv2 =S*np.sqrt(conv1)

#np.sqrt(conv1*(Varn(s)**2) + (S**2/conv1)*conv1var**2)

def VarN(x):

  avar = np,var(x)
  return avar/np.sqrt(len(x)-1)

def CalculateCorrection(event):

  variance = np.var(event.r0.tel[1].waveform,axis=2)[0]
  return variance

def powerlaw(x, a, b):


  return(a * (np.abs(x)**b))

config_calib = Config({"LSTCalibrationCalculator" : {"sample_size" : 100}})

#source = ProtozfitsDL0TelescopeEventSource("/mnt/c/Users/ctoen/Documents/TEL001_SDH3001_20231015T033032_SBID0000000002000000039_OBSID0000000002000000107_CHUNK001.fits.fz", max_events=100000)

source = EventSource("/mnt/c/Users/ctoen/Documents/interleaved/interleaved_LST-1.Run16643.0000.h5")

crab = SkyCoord.from_name("Crab Nebula")

sub = source.subarray

Charger = extractor.FixedWindowSum(source.subarray, peak_index = 18, window_shift = 6, window_width = 12, apply_integration_correction = False)

geometry = sub.tel[1].camera.geometry

FFC = FlasherFlatFieldCalculator(subarray=source.subarray,tel_id=1,sample_size=100)

#print(FFC.tel_id)

charges   = []
variances = []

calibrator = LSTCalibrationCalculator(subarray=source.subarray,config = config_calib)

calibrator.tel_id = 1

calibrator.pedestal.tel_id = 1
calibrator.flatfield.tel_id = 1

calibrator.pedestal.sample_size = 100
calibrator.flatfield.sample_size = 100

#print(calibrator.tel_id)
for i,event in enumerate(source):

  if event.trigger.event_type.value == 2 or event.trigger.event_type.value == 0:

    event.mon.tel[1].pixel_status.hardware_failing_pixels = np.zeros(event.r1.tel[1].waveform.shape[:2])

    ped, ff = calibrator.process_interleaved(event)

    #print(event)

    #print(event.mon.tel[1].pixel_status.hardware_failing_pixels)

    if ped:

      print("pedestal")
      print(event)

    if ff:

      print("flatfield")
      print(event)


    #print(calibrator.flatfield.charge_medians)

    print("Ped ", calibrator.pedestal.num_events_seen, calibrator.pedestal.sample_size)
    print("FF ", calibrator.flatfield.num_events_seen, calibrator.flatfield.sample_size)

    #print(calibrator.output_interleaved_results(event))

    n_pixels   = len(event.r1.tel[1].waveform[0])
    n_channels = len(event.r1.tel[1].waveform)

    broken_pixels = np.zeros(n_pixels, dtype=bool)

    no_gain_selection = np.zeros((n_channels, n_pixels), dtype=np.int64)

    eventtype = ""

    event.meta.update({'origin':'LST'})

    event.r1.tel[FFC.tel_id].selected_gain_channel = no_gain_selection

    if event.trigger.event_type.value == 2:

      eventtype = "SKY_PEDESTAL"

    else:

      eventtype = "FLATFIELD"

    if event.trigger.event_type.value == 0:     

      image = np.var(event.r1.tel[1].waveform,axis=2)[0]

      charge = Charger(event.r1.tel[1].waveform,1,no_gain_selection,broken_pixels)
   
      plt.plot(np.var(event.r1.tel[1].waveform,axis=2)[1], np.mean(event.r1.tel[1].waveform,axis=2)[1],"k+")

      xvals = np.arange(0,35,0.1)

      #plt.plot(xvals,[powerlaw(x,*popt) for x in xvals],"r-")

      plt.xlabel("Mean")
      plt.ylabel("Variance")

      plt.savefig(str(i) + "_" + eventtype + "_correlation_2.png", format = "png")

    plt.clf()

    """

    plt.figure(figsize=(10, 10))

    
    disp = CameraDisplay(geometry, image=image)
    disp.add_colorbar()
    plt.savefig(str(i) + "_" + eventtype + "_LOW.png", format = "png")
    plt.close()
    """
    if event.trigger.event_type.value == 2:

      image = np.var(event.r1.tel[1].waveform,axis=2)[0]

    elif event.trigger.event_type.value == 0:

      image = np.var(event.r1.tel[1].waveform,axis=2)[0]

      variances.append(image)
      #print(image)

      charge = Charger(event.r1.tel[1].waveform,1,no_gain_selection,broken_pixels)

      charges.append(charge.image[0])

    '''
    plt.figure(figsize=(10, 10))
    disp = CameraDisplay(geometry, image=image)
    disp.add_colorbar()
    plt.savefig(str(i) + "_" + eventtype + "_HIGH.png", format = "png")
    plt.close()
    '''

variances = np.array(variances)
charges   = np.array(charges)

bins = np.where(variances[0] < 1.0,True,False)

for i,check in enumerate(bins):

  if check:
 
    plt.hist2d(charges[:,i],variances[:,i],bins=(np.arange(-2.0,25.0,0.1),np.arange(0.0,3.0,0.02)))
    plt.xlabel("Charge [PE]")
    plt.ylabel("Variance")
    plt.savefig("Histo_" + str(i) + "_histo.png", format = "png")


#source = EventSource("/mnt/c/Users/ctoen/Documents/TEL001_SDH3001_20231015T033032_SBID0000000002000000039_OBSID0000000002000000107_CHUNK001.fits.fz", max_events=100)

