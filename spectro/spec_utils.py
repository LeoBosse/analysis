#!/usr/bin/python3
# -*-coding:utf-8 -*

import numpy as np
from scipy.interpolate import interpn
from scipy.signal import fftconvolve
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import os
from subprocess import call
import datetime as dt

# import radis
import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun, get_body
from astropy.time import Time
from pybaselines import Baseline, utils
from scipy.ndimage import uniform_filter1d

import spectrapepper as spep


from tkinter import Tk
from tkinter.filedialog import askdirectory
from tkinter import filedialog

data_path = '/home/leob/iasb/analysis/data/skibotn_spectro/' 


class Serie:
    """
    Class containing a spectra kinetic serie. Used to load fits files, clean the data, plot them or isolate emission lines.
    Each spectra is described by a Spectrum object, defined below.
    """
    
    path = data_path

    def __init__(self, header, data, **kwargs):
        """Initialize a time serie from a fits header and a 2d data array, plus optional args."""

        self.header = header

        # Set important parameters
        self.acq_mode = self.header['ACQMODE']
        self.shape = data.shape
        self.nb_spec, self.nb_wavelengths = self.shape

        self.exposure = self.header['EXPOSURE']
        if 'ACT' in self.header.keys():
            self.cycle_time = dt.timedelta(seconds = self.header['ACT']) #exposure + few 1/100 sec
        else:
            self.cycle_time = dt.timedelta(seconds = self.exposure) 

        # Set wavelengths list
        self.wavelengths = self.GetWavelengthList()
        self.spectral_res = np.average(np.diff(self.wavelengths))

        # Set times of each spectra in the serie
        self.start_datetime = dt.datetime.strptime(self.header['DATE'], '%Y-%m-%dT%H:%M:%S')
        self.times = self.start_datetime + np.arange(0, self.nb_spec) * self.cycle_time
        self.timestamps = np.array([t.timestamp() for t in self.times])
        
        self.start = self.times[0]
        self.end = self.times[-1]
        
        self.lines = []
        self.color = 'k'
        self.delay = dt.timedelta(seconds = 0)

        # Create a list containing all Spectrum object
        if 'dont_init_spectra' in kwargs.keys():
            self.spectra = []
        else:
            self.InitSpectraList(data)

        self.GetIntegral()


    def InitSpectraList(self, data):
        """Create a list of Spectrum objects from a 2d data array"""
        self.spectra = [Spectrum(self.header, data[s], wavelengths=self.wavelengths, resolution=self.spectral_res, datetime=self.times[s]) for s in range(self.nb_spec)]

    
    @classmethod
    def FromData(cls, header, data, wavelengths, times):
        """
        Class method returning a Serie object from a 2d data array, the corresponding wavelenghts and times to use as axis.
        Used to create line series easily.
        """
        
        # exposure = times[1] - times[0]
        # header['EXPOSURE'] = exposure.total_seconds()
        # header['ACT'] = exposure.total_seconds()

        serie = Serie(header, data, dont_init_spectra=True)

        serie.wavelengths = wavelengths
        serie.start_datetime = times[0]
        serie.times = times
        serie.timestamps = np.array([t.timestamp() for t in serie.times])
        
        serie.InitSpectraList(data)

        return serie

    @classmethod
    def FromFitsFile(cls, file_name):
        """
        Class method returnong a Serie object from a file name. 
        The file must be a fits file containing a 2d data array with each spectrum as a line.
        """
        if not file_name:
            root = Tk()
            root.withdraw()
            root.update()
            file_name = filedialog.askopenfilename(initialdir = Serie.path, title="Select spectrograph data file (.fits).", filetypes=[("fits files", "*.fits")])
            root.destroy()    

        header, raw_data = Serie.LoadFitsFile(file_name)

        return Serie(header, raw_data)


    @classmethod
    def LoadFitsFile(cls, file_name):
        """
        From a file name (.fits), returns its header and data array
        """
        header = fits.getheader(file_name)
        raw_data = fits.getdata(file_name)
        return header, raw_data


    def SetDark(self, dark_file_30, dark_file_60):
        """
        Set the dark used throughout the Serie, based on the first spectrum of the serie. 
        The same dark will be applied to all spectrum.
        """
        self.dark_current, self.dark_offset, self.dark = self.spectra[0].ComputeDark(dark_file_30, dark_file_60)

        # self.darks = np.zeros_like(self.data)
        # for i, s in enumerate(self.data):
        #     # darks[i] = dark * np.median(s) / np.median(dark)
        #     self.darks[i] = self.dark #* np.median(s) / np.median(dark)
        

    def Subtract(self, to_sub):
        """
        Subtract a spectrum from each spectrum of the serie. used for Dark, Offset, Baseline subtraction.
        to_sub must be a 1d array of the same size as the spectrum of the serie.
        """
        for s in self.spectra:
            s.Subtract(to_sub)
    
    def SubtractDarks(self):
        """
        Subtract the dark of the serie from each spectrum. Must call self.SetDark() first.
        """
        for s in self.spectra:
            s.Subtract(self.dark)
     
     
    def SubtractBaseLines(self):
        """
        Subtract the its baseline from each spectrum. Each spectrum has a different baseline.
        The Baselines must be computed first with self.SetBaseLines().
        """
        for s in self.spectra:
            s.Subtract(s.baseline)

    def SetBaseLines(self):
        """Set the baseline of each spectrum of the serie."""
        for s in self.spectra:
            s.SetBaseLine()


    def GetWavelengthList(self):
        """
        Returns an array containing the wavelengths of the serie. Same is used for all spectra in the serie.
        Computed from the 'CALIB' param of the fits header as a polynomial af degree <= 3
        """
        calibration = [float(c) for c in self.header['CALIB'].split(',')]
        pixels = np.arange(self.header['NAXIS1'])
        wavelengths = pixels**3 * calibration[3] + pixels**2 * calibration[2] + pixels * calibration[1] + calibration[0]
        return wavelengths


    def MaskSunMoon(self, sun_limit=None, moon_limit=None):
        """
        Use astropy to remove all spectra polluted from the Sun and/or the moon.
        Give the maximum elevation of the sun/moon in degrees. If not given, no filters will be applied.
        Use the time of the header computed during initialisation of the Serie and the location of Skibotn observatory.
        """
        if sun_limit is None and moon_limit is None:
            return
        
        #Set the location of the observation and the times of each spectra in the serie.
        skibotn = EarthLocation(lat=69.34822*u.deg, lon=20.36323*u.deg, height=150*u.m)
        # utcoffset = 2*u.hour  # Eastern Daylight Time
        astro_times = Time(self.start_datetime.strftime('%Y-%m-%d %H:%M:%S')) + (self.timestamps - self.timestamps[0])*u.s
        frame_night = AltAz(obstime=astro_times,
                            location=skibotn)
        
        if sun_limit is not None: # If an upper limit on Sun elevation is given
            sun = get_sun(astro_times).transform_to(frame_night) #Get astropy Sun object
            sun_mask = np.where(sun.alt < sun_limit*u.deg) #Get indices of the Spectra when the sun is below the upper elevation limit
            self.ApplyTimeMask(sun_mask) #Remove all spectra where the sun is too high in the sky

        if moon_limit is not None:
            moon = get_body("moon", astro_times).transform_to(frame_night)
            moon_mask = np.where(moon.alt < moon_limit*u.deg)
            self.ApplyTimeMask(moon_mask)


    def ApplyTimeMask(self, mask):
        """
        Given a list of indices, clean the Serie to keep only the spectra at those indices.
        """
        # print(mask)
        self.spectra = np.array(self.spectra)[mask]
        self.times = self.times[mask]
        
        self.start = self.times[0]
        self.end = self.times[-1]

        self.timestamps = self.timestamps[mask]
        self.nb_spec = len(self.spectra)
        self.shape = self.nb_spec, self.nb_wavelengths


    def GetGradient(self, data = None):
        """
        Returns the gradient of the data given in input (or the self.spectra data if not given).
        This is a size 2 list containing the gradient (2d array) with respect to time and to wavelength in that order.
        """
        if data is None:
            data = self.data
        return np.gradient(data, self.cycle_time.total_seconds(), self.wavelengths)


    # def GetCosmics(self):
    #     self.gradient = self.GetGradient()
    #     # self.cosmics  = self.gradient[1] > 100 / self.cycle_time.total_seconds()

    #     # self.cosmics  = self.gradient[1] > 200 / self.spectral_res

    #     conditions = [np.abs(self.gradient[0]) > 8, np.abs(self.gradient[1]) < 10]
        
    #     self.cosmics = np.all(conditions, axis = 0)


    #     cosmic_signal = np.array([[0, 0, 0],
    #                               [0, 1, 0],
    #                               [0, 0, 0]])

    #     cosmic_convolve = fftconvolve(self.data, cosmic_signal, mode='same')

    #     # plt.figure()
    #     # plt.imshow(cosmic_convolve)
    #     # plt.plot()

    # def RemoveCosmicRays(self, method='dd'):

    #     match method:
    #         case 'dd':  clean_data = spep.cosmicdd(self.data, th=100, asy=0.6745, m=10)  # from Whitaker and Hayes
    #         case 'mp':  clean_data = spep.cosmicmp(self.data, alpha=15, avg=10)  # from Barton and Hennelly
    #         case 'med': clean_data = spep.cosmicmed(self.data, sigma=100)  # elimination by similarity
    #         case _: 
    #             print("ERROR: cosmic ray removal method invalid. No correction applied. Correct methods: 'dd', 'med', 'mp'")
    #             clean_data = self.data
        
    #     self.cosmics = clean_data != self.data

    def GetTimeEvolution(self, wl):
        """
        For a given wavelengths, returns an array of its count value at each time.
        Counts are interpolated for wavelengths between 2 wavelengths values of self.wavelengths.
        """
        evol = np.empty(self.nb_spec)
        for i, spec in enumerate(self.spectra):
            counts = np.interp(wl, self.wavelengths, spec.data)
            evol[i] = counts
        return evol

    def AddLine(self, wl, half_width, color = 'k', delay = 0):
        """
        Add a Serie object in the self.lines list.
        See self.GetLine() to know how this Serie is computed.
        """
        self.lines.append(self.GetLine(wl, half_width, color=color, delay = delay))


    def GetLine(self, wl, half_width, color='k', delay = 0):
        """
        Returns a Serie object computed from a subset of this Serie.
        The subset is defined by a center wavelength (wl) and a width (half_width).
        The new Serie contains only the data in this range of wavelengths.
        Usefull to isolate an emisison line.
        color='k': can be used to plot it in a different color than other lines. Black by default
        delay=0: in seconds. desexitation delay of the line. For example for the auroral red line, delay = 100sec
        """

        #Get the wavelengths indices that fall in the correct wavelength range.
        mask = np.where(self.wavelengths < wl + half_width)
        mask = np.where(self.wavelengths[mask] > wl - half_width)
        
        data = np.array([s.data for s in self.spectra])
        
        # Create a Serie object from self.header and self.data cropped to the correct wavelengths range
        line = Serie.FromData(self.header, data[:, mask][:, 0, :], self.wavelengths[mask], self.times)

        # Set optional arguments of the line
        line.color = color
        line.delay = dt.timedelta(seconds = delay)

        return line


    def GetIntegral(self):
        """
        Compute, set and returns the integral of each spectrum over all wavelengths.
        Used mainly for single emission lines Series.

        """
        self.integral = np.array([s.GetIntegral() for s in self.spectra])
        return self.integral

    
    def GetRatio(self, div_line):
        """
        Returns the ratio of self / div_line.
        div_line must be another Serie object containing the same number of spectra.
        If both Series have a relative time delay, will trim the non overlapping spectra, and the results might not have the same size.
        """

        # relative time delay between 2 series.
        line_delay = int((self.delay - div_line.delay).total_seconds() / self.exposure)
        # print(line_delay)
        
        # trim the non overlapping spectra at the start or end of the series.
        if line_delay >= 0: 
            delayed_integral = self.integral[:self.nb_spec - line_delay]
            div_delayed_integral = div_line.integral[line_delay:]
        else:
            line_delay = abs(line_delay)
            delayed_integral = self.integral[line_delay:]
            div_delayed_integral = div_line.integral[:div_line.nb_spec - line_delay ]

        # print(len(delayed_integral), len(div_delayed_integral))
        return delayed_integral / div_delayed_integral

    def GetAllSpectraData(self, spectra = None):
        """
        Returns the data of all spectra in the serie as a 2d numpy array.
        Each spectra is a line of the array.
        """
        if spectra is None:
            spectra = self.spectra

        return np.array([s.data for s in spectra])


    def GetMax(self):
        data = self.GetData()
        index = np.unravel_index(np.argmax(data), data.shape)
        print(index)
        return index, data[index]
    
    def GetPercentile(self, P):
        return np.percentile(self.GetData(), P)
    
    def GetMin(self):
        index = np.argmin(self.GetData())
        return index, self.GetData()[index]

    def GetData(self):
        return np.array([s.data for s in self.spectra])

    def PlotSpectrum(self, ax, spec_id, spectra=None):
        """
        given a matplotlib axes object, will plot the a specific Spectrum of the serie on it.
        """
        if spectra is None:
            spectra = self.spectra

        spectra[spec_id].Plot(ax, label = f"{spec_id} {self.times[spec_id].strftime('%Y-%m-%d %H:%M:%S')}")

        ax.plot(self.wavelengths, self.darks[spec_id])

        for i in np.where(self.cosmics[spec_id]):
            x = [self.wavelengths[i], self.wavelengths[i]]
            y = [0, np.max(spectra[spec_id].data)]
            ax.plot(x, y, 'r')

        ax.legend()


    def MakeSpectrumFigure(self, spec_id, spectra=None):
        """
        Makes a matplotlib figure containing a single spectrum.
        """
        if spectra is None:
            spectra = self.spectra

        spectra[spec_id].MakeFigure()

        # fig, axs = plt.subplots(2, 1, sharex=True)
        # self.PlotSpectrum(axs[0], spec_id, spectra=spectra)
        # # self.PlotSpectrum(axs[0], spec_id+1)
        # # self.PlotSpectrum(axs[0], spec_id-1)

        # # axs[1].plot(self.wavelengths, np.abs(self.gradient[1][spec_id]))
        # axs[1].plot(self.wavelengths, np.abs(self.gradient[0][spec_id]), 'r')

        # axs[1].plot(self.wavelengths, np.abs(self.gradient[0][spec_id-1]), 'k')
        # axs[1].plot(self.wavelengths, np.abs(self.gradient[0][spec_id+1]), 'g')

    def PlotImage(self, spectra=None, wavelengths=None, timestamps=None):
        """
        Makes a plt.figure with the serie as a image.
        """
        if spectra is None:
            spectra = self.spectra 
            wavelengths = self.wavelengths
            timestamps = self.timestamps
        
        fig = plt.figure()
        # W, T = np.meshgrid()
        extent = [0, (timestamps[-1] - timestamps[0])/3600, wavelengths[-1], wavelengths[0]]
        data = self.GetAllSpectraData(spectra)
        plt.imshow(data.T, extent = extent, interpolation='none', aspect="auto")


    def Plot3D(self, spectra=None, wavelengths=None, timestamps=None):
        """
        Makes a plt.figure with the serie as a 3d plot.
        """
        if spectra is None:
            spectra = self.spectra
            wavelengths = self.wavelengths
            timestamps = self.timestamps

        fig = plt.figure()
        ax = plt.axes(projection ='3d')
        W, T = np.meshgrid(wavelengths, timestamps/3600)
        data = self.GetAllSpectraData()
        ax.plot_surface(T, W, data)


    def MakeAnimation(self):
        """
        Makes a plt.figure with an animation of the spectra over time.
        """
        fig = plt.figure() # initialise la figure
        ax = fig.add_subplot(1, 1, 1)
        line, = self.spectra[0].Plot(ax, label = '0 ' + self.times[0].strftime('%Y-%m-%d %H:%M:%S'))
        L=plt.legend(loc=1) #Define legend objects
        ax.set_ylim(top = self.GetPercentile(99.9) * 1.5)


        def animate(i): 
            line.set_data(self.wavelengths, self.spectra[i].data)
            L.get_texts()[0].set_text(str(i) + ' ' + self.times[i].strftime('%Y-%m-%d %H:%M:%S')) #Update label each at frame
            

            return line, L
        
        self.ani = animation.FuncAnimation(fig, animate, frames=self.nb_spec,
                                    interval=100, blit=True, repeat=False)


class Spectrum:
    """
    Class containing a single spectrum. Used to initialize a spectrum, clean and modify it and plot it.
    """
    path = data_path

    def __init__(self, header, data, **kwargs):
        """
        Initialize a spectrum from a header (dictionnary) and a 1d data array.
        """

        self.header, self.raw_data = header, data
        self.data = self.raw_data.copy()

        self.shape = self.raw_data.shape
        self.nb_wavelengths = self.shape[0]

        self.exposure = self.header['EXPOSURE']

        self.wavelengths  = kwargs['wavelengths'] if 'wavelengths' in kwargs.keys() else self.GetWavelengthList()
        self.spectral_res = kwargs['resolution'] if 'resolution' in kwargs.keys()   else np.average(np.diff(self.wavelengths))
        self.datetime     = kwargs['datetime'] if 'datetime' in kwargs.keys()       else dt.datetime.strptime(self.header['DATE'], '%Y-%m-%dT%H:%M:%S')

        self.timestamp = self.datetime.timestamp()

        self.lines = []
        self.integral = None
        self.color = 'k'


    @classmethod
    def FromFitsFile(cls, file_name=None):
        """
        Class method returning a Spectrum object from a file name. 
        The file must be a fits file containing a 2d data array. (with only one line)
        """
        if not file_name:
            root = Tk()
            root.withdraw()
            root.update()
            file_name = filedialog.askopenfilename(initialdir = Serie.path, title="Select spectrograph single spectrum data file (.fits).", filetypes=[("fits files", "*.fits")])
            root.destroy()    

        header, raw_data = Spectrum.LoadFitsFile(Spectrum.path + file_name)
        return Spectrum(header, raw_data)

    @classmethod
    def LoadFitsFile(cls, file_name):
        """
        From a file name (.fits), returns its header and data array.
        The file must be contain a 2d data array with only one line.
        """
        header = fits.getheader(file_name)
        raw_data = fits.getdata(file_name)[0]
        return header, raw_data



    def GetWavelengthList(self):
        """
        Returns an array containing the wavelengths of the spectrum.
        Computed from the 'CALIB' param of the fits header as a polynomial af degree <= 3
        """
        calibration = [float(c) for c in self.header['CALIB'].split(',')]
        pixels = np.arange(self.nb_wavelengths)
        wavelengths = pixels**3 * calibration[3] + pixels**2 * calibration[2] + pixels * calibration[1] + calibration[0]
        return wavelengths


    def ComputeDark(self, dark_file_30, dark_file_60):
        """
        Input: 2 fits files names containing 2 dark accumulation data of different exposure.
        Output: the dark current, offset and total dark of the spectrograph.
        """
        dark_header_30, dark_raw_data_30 = Serie.LoadFitsFile(dark_file_30)
        dark_header_60, dark_raw_data_60 = Serie.LoadFitsFile(dark_file_60)
        dark_30 = dark_raw_data_30[0] / dark_header_30['NUMACC']
        dark_60 = dark_raw_data_60[0] / dark_header_60['NUMACC']

        self.dark_current = (dark_60 - dark_30) / (dark_header_60['EXPOSURE'] - dark_header_30['EXPOSURE'])
        self.dark_offset =  dark_60 - self.dark_current * dark_header_60['EXPOSURE']
        self.dark_offset += dark_30 - self.dark_current * dark_header_30['EXPOSURE']
        self.dark_offset /= 2.
        
        self.dark = self.dark_current * self.exposure + self.dark_offset

        # self.darks = np.zeros_like(self.data)
        # for i, s in enumerate(self.data):
        #     # darks[i] = dark * np.median(s) / np.median(dark)
        #     self.darks[i] = self.dark #* np.median(s) / np.median(dark)

        return self.dark_current, self.dark_offset, self.dark        


    def Subtract(self, to_sub=None):
        """
        Subtract the to_sub 1d data array from the spectrum. If to_sub is not given, use the spectrum dark (must call ComputeDark first.).
        """
        if to_sub is None:
            to_sub = self.dark
        self.data -= to_sub

    def SetBaseLine(self, algo = 'modpoly'):
        """
        Use pybaselines package to compute the baseline of the spectrum.
        Several algorithm are available, default 'modpoly': modpoly, asls, mor, snip
        """

        baseline_fitter = Baseline(x_data=self.wavelengths)
        self.baseline = np.zeros_like(self.wavelengths)
        
        y = uniform_filter1d(self.data, 11)

        match algo:
            case 'modpoly':
                self.baseline = baseline_fitter.modpoly(y, poly_order=3)[0]
            case 'asls':
                self.baseline = baseline_fitter.asls(y, lam=1e7, p=0.02)[0]
            case 'mor':
                self.baseline = baseline_fitter.mor(y, half_window=30)[0]
            case 'snip':
                self.baseline = baseline_fitter.snip(y, max_half_window=40, decreasing=True, smooth_half_window=3)[0]
            case _:
                print('Warning: the baseline algorithm is invalid. Baseline set to 0.')


    def SetIntegral(self):
        """Compute and set the integral of the spectra ovor all wavelengths."""
        self.integral = np.sum(self.data)
    

    def GetIntegral(self):
        """Set and returns the integral of the spectra over all wavelengths."""
        if self.integral is None:
            self.SetIntegral()
        return self.integral


    def Plot(self, ax, label=None, **kwargs):
        """Plot the spectra on a plt.axis abject"""
        if label is None:
            f"{self.datetime.strftime('%Y-%m-%d %H:%M:%S')}"
        im = ax.plot(self.wavelengths, self.data, color=self.color, label = label, **kwargs)
        return im


    def MakeFigure(self, **kwargs):
        """Make a figure with the spectra."""
        f = plt.figure()
        ax = f.add_subplot(1, 1, 1)
        ax.plot(self.wavelengths, self.data, self.color, label = f"{self.datetime.strftime('%Y-%m-%d %H:%M:%S')}", **kwargs)




