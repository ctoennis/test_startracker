"""
Tool to check the status of each pixel in a run
"""

from abc import abstractmethod

import numpy as np
from astropy import units as u

from ctapipe.containers import DL1CameraContainer
from ctapipe.core import Component
from ctapipe.core.traits import Int, Float, Unicode
from ctapipe.image.extractor import ImageExtractor

from matplotlib import pyplot as plt

import matplotlib.colors as mcolors

__all__ = ["PixelStatusChecker","PixelStatusVariance"]

class PixelStatusChecker(Component):

    """
    Parent class for the Pixel status checkers.
    Performs check on each pixel and sees what state it is in and creates a record over time. Also calculates 
   
    Parameters
    ----------
    tel_id : int
          id of the telescope (default 0)
    sample_size : int
         number of Flat Field events requested for the statistics
    n_channels : int
         number of waveform channel to be considered
    config : traitlets.loader.Config
        Configuration specified by config file or cmdline arguments.
        Used to set traitlet values.
        Set to None if no configuration to pass.

    kwargs
    """

    tel_id = Int(0, help="id of the telescope to calculate the pedestal values").tag(config=True)
    sample_size = Int(10000, help="sample size").tag(config=True)
    n_channels = Int(2, help="number of channels to be treated").tag(config=True)
    charge_product = Unicode("FixedWindowSum", help="Name of the charge extractor to be used").tag(config=True)

    def __init__(self, subarray, **kwargs):

        """
        Parent class for the Pixel status checkers.
        Performs check on each pixel and sees what state it is in and creates a record over time. Also calculates 
   
        Parameters
        ----------
        tel_id : int
            id of the telescope (default 0)
        sample_size : int
            number of Flat Field events requested for the statistics
        n_channels : int
            number of waveform channel to be considered
        charge_product : str
            Name of the charge extractor to be used
        config : traitlets.loader.Config
            Configuration specified by config file or cmdline arguments.
            Used to set traitlet values.
            Set to None if no configuration to pass.

        kwargs
        """

        super().__init__(**kwargs)

        # load the waveform charge extractor
        self.extractor = ImageExtractor.from_name(
            self.charge_product, parent=self, subarray=subarray
        )
        self.log.info(f"extractor {self.extractor}")

class PixelStatusVariance(PixelStatusChecker):

    """
    Creates a record of th estatus of each pixel based on Flatfield events regarding reduced voltage. 
   
    Parameters
    ----------
    Varcut_Low : Float
        Threshold value applied to the low gain channel to find reduced HV pixels
    Varcut_High : Float
        Threshold value applied to the high gain channel to find reduced HV pixels

    kwargs
    """

    Varcut_Low = Float(0,help="Threshold value applied to the low gain channel to find reduced HV pixels").tag(config=True)
    Varcut_High =  Float(0,help="Threshold value applied to the high  gain channel to find reduced HV pixels").tag(config=True)

    def __init__(self, **kwargs):

        """
        Creates a record of th estatus of each pixel based on Flatfield events regarding reduced voltage.

        Parameters
        ----------
        Varcut_Low : Float
            Threshold value applied to the low gain channel to find reduced HV pixels
        Varcut_High : Float
            Threshold value applied to the high gain channel to find reduced HV pixels

        kwargs
        """

        super().__init__(**kwargs)

        self.log.info("Used events statistics : %d", self.sample_size)

        self.n_events_seen = 0
        self.time_start = None  # trigger time of first event in sample

    def CheckStatus(self,event):

        """
        Check the status of each pixel in the event

        Parameters
        ----------
        event : general event container

        """

        # initialize the np array at each cycle
        waveform = event.r1.tel[self.tel_id].waveform

        # re-initialize counter
        if self.n_events_seen == self.sample_size:
            self.n_events_seen = 0

        # real data
        trigger_time = event.trigger.time

        if self.n_events_seen == 0:
            self.time_start = float(trigger_time.unix_tai)
            self.setup_sample_buffers(waveform, self.sample_size)

        varpic = np.var(event.r1.tel[self.tel_id].waveform,axis=2)

        meanpic = np.mean(event.r1.tel[self.tel_id].waveform,axis=2)

        mean_low  = np.where(meanpic[1] < 1.0,False, True)
        mean_high = np.where(meanpic[0] < 1.0,False, True)
        var_low   = np.where(varpic[1]  < self.Varcut_Low,False, True)
        var_high  = np.where(varpic[0]  < self.Varcut_High,False, True)

        marks_mean = np.logical_and(mean_low,mean_high)
        marks_var  = np.logical_and(var_low,var_high)

        marks = np.logical_or(marks_mean, marks_var)

        self.marks[self.n_events_seen] = marks
        self.times[self.n_events_seen] = float(trigger_time.unix_tai) - self.time_start

        self.n_events_seen += 1

        if self.n_events_seen == self.sample_size:

            self.calculate_status_results("Status_")



    def setup_sample_buffers(self, waveform, sample_size):

        """Initialize sample buffers"""

        n_channels = waveform.shape[0]
        n_pix = waveform.shape[1]

        self.marks = np.zeros((sample_size, n_pix))
        self.times = np.zeros(sample_size)



    def calculate_status_results(self,name):

        """Make a plot of the status over time"""

        colors = list(mcolors.TABLEAU_COLORS.keys())

        for pix in range(self.marks.shape[1]):
            
            plt.plot(self.times,self.marks[:,pix],color=colors[pix%10],label= "status pixel " + str(pix))
            plt.xlabel("Time [s]")
            plt.ylabel("Status")

            if pix%10 == 9:

                plt.legend()
                
                plt.savefig(name+str(pix - (pix%10)))
                plt.clf()
