import numpy as np
from matplotlib import pyplot as plt
from scipy import signal
from scipy import optimize as opt

from Pixel_checker import PixelStatusVariance

from astropy import units as u
from astropy.coordinates import SkyCoord

from lstchain.calib.camera.calibration_calculator import LSTCalibrationCalculator
from ctapipe.image import extractor
from ctapipe import utils
from ctapipe.visualization import CameraDisplay
from ctapipe.io import EventSource, DataWriter
from ctapipe_io_zfits import ProtozfitsDL0TelescopeEventSource
from ctapipe.visualization import ArrayDisplay, CameraDisplay
from ctapipe.calib import CameraCalibrator
from ctapipe.calib.camera.gainselection import ManualGainSelector
from ctapipe.calib.camera.flatfield import FlasherFlatFieldCalculator
from ctapipe.calib.camera.pedestals import PedestalIntegrator


import eventio

import copy

#Sprime = S*conv2 =S*np.sqrt(conv1)

#np.sqrt(conv1*(Varn(s)**2) + (S**2/conv1)*conv1var**2)

gain_channel = 1

def VarN(x):

  avar = np,var(x)
  return avar/np.sqrt(len(x)-1)
 
#source = ProtozfitsDL0TelescopeEventSource("/mnt/c/Users/ctoen/Documents/TEL001_SDH3001_20231015T033032_SBID0000000002000000039_OBSID0000000002000000107_CHUNK001.fits.fz", max_events=100000)

source = EventSource("/mnt/c/Users/ctoen/Documents/interleaved/interleaved_LST-1.Run15134.0000.h5")

Stat = PixelStatusVariance(subarray=source.subarray,tel_id=1,sample_size=300,Varcut_Low = 3.0,Varcut_High = 3.0)

#Charger = extractor.FixedWindowSum(source.subarray, peak_index = 18, window_shift = 6, window_width = 12, apply_integration_correction = False)

Charger = extractor.LocalPeakWindowSum(source.subarray, window_shift = 5, window_width = 12, apply_integration_correction = False) 

sub = source.subarray

geometry = sub.tel[1].camera.geometry

#FC = FlasherFlatFieldCalculator(subarray=source.subarray,charge_product="FixedWindowSum",tel_id=1,sample_size=100)


#PIN = PedestalIntegrator(subarray=source.subarray,tel_id=1,sample_size=300)

#calib = CameraCalibrator(subarray=source.subarray)

calibrator = LSTCalibrationCalculator(subarray=source.subarray)

calibrator.tel_id = 1

calibrator.pedestal.tel_id = 1
calibrator.flatfield.tel_id = 1

calibrator.pedestal.sample_size = 400
calibrator.flatfield.sample_size = 410

mon_found  = False
flat_found = False

config = {}

with DataWriter(event_source=source,output_path="/mnt/c/Users/ctoen/Documents/processed.h5",overwrite=True) as write_data:

  for i,event in enumerate(source):

    print("Were at: ",i)

    if event.trigger.event_type.value == 2 or event.trigger.event_type.value == 0:

      n_pixels   = len(event.r1.tel[1].waveform[0])
      n_channels = len(event.r1.tel[1].waveform) 

      event.mon.tel[1].pixel_status.hardware_failing_pixels = np.zeros(event.r1.tel[1].waveform.shape[:2])

      event.meta.update({'origin':'LST'})

      if mon_found:

        event.mon.tel[1].pedestal = mondata

      if flat_found:

        event.mon.tel[1].flatfield = flatdata

      ped, ff = calibrator.process_interleaved(event)

      if ped:

        print("pedestal")
        print(event.mon.tel[1].pedestal)

        mondata = event.mon.tel[1].pedestal
        mon_found = True

      if ff:

        print("flatfield")
        print(event.mon.tel[1].flatfield)

        flatdata = event.mon.tel[1].flatfield
        flat_found = True
        plt.figure(figsize=(10, 10))

        bins = np.arange(40.0,80.0,1.0)
        print("shape: ", calibrator.flatfield.charges.shape)
        meancharges = np.median(calibrator.flatfield.charges,axis=0)
        meancharges = meancharges - event.mon.tel[1].pedestal.charge_median
        meancharges_corrected = np.divide(meancharges,event.mon.tel[1].flatfield.relative_gain_median)
        print("widths: ", np.var(meancharges[gain_channel]), np.var(meancharges_corrected[gain_channel]))
        plt.hist(meancharges[gain_channel],bins=bins,histtype='bar',ec="black",fill=False,label="uncalibrated")
        plt.hist(meancharges_corrected[gain_channel],bins=bins,histtype='bar',ec="red",fill=False,label="calibrated")
        plt.legend()
        plt.savefig(str(i) + "_histo.png", format = "png")
        plt.clf()

        disp = CameraDisplay(geometry, image=event.mon.tel[1].flatfield.relative_gain_mean[1])
        disp.add_colorbar()
        plt.savefig(str(i) + "_gain.png", format = "png")
        plt.close()

    if flat_found and mon_found:

      no_gain_selection = np.zeros((n_channels, n_pixels), dtype=np.int64)

      broken_pixels = np.zeros(n_pixels, dtype=bool)

      uncorrected_charge = Charger(event.r1.tel[1].waveform,1,no_gain_selection,broken_pixels).image[gain_channel]-event.mon.tel[1].pedestal.charge_mean[gain_channel]

      gain_corrected_charge = np.multiply(uncorrected_charge,event.mon.tel[1].flatfield.relative_gain_mean[gain_channel])

      uncorrected_image = np.var(event.r1.tel[1].waveform,axis=2)[gain_channel]
       
      gain_corrected_image = np.multiply(uncorrected_image,np.sqrt(event.mon.tel[1].flatfield.relative_gain_median[gain_channel]))

      plt.figure(figsize=(10, 10))

      disp = CameraDisplay(geometry, image=gain_corrected_image)
      disp.add_colorbar()

      if event.trigger.event_type.value == 2:
      
        plt.savefig(str(i) + "_corrected_sf.png", format = "png")

      else:

        plt.savefig(str(i) + "_corrected_ff.png", format = "png")

      plt.close()

      plt.plot(uncorrected_charge,uncorrected_image,"k+",label="Uncalibrated data")
      #print(gain_corrected_charge.shape,gain_corrected_image.shape)
      plt.plot(gain_corrected_charge,gain_corrected_image,"r+",label="Calibrated data")
      plt.xlabel("Charge [PE]")
      plt.ylabel("Variance of waveform")
      plt.title("Single image")
      plt.legend()
      plt.savefig(str(i) + "_correl.png", format = "png")
      plt.clf()

    write_data(event)
