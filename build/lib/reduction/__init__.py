# Reduction od DATA Astronomy

"""
Reduction: Python Package to Analysis on Astronomy Data
"""

import numpy as np
from astropy.io import fits
import lmfit
import pandas as pd

from photutils import CircularAperture
from photutils import aperture_photometry
# from photutils import CircularAnnulus
# from photutils import SigmaClip
# from photutils import Background2D, MedianBackground

from photutils import daofind
from astropy.stats import median_absolute_deviation as mad

class psf(object):

    def __init__(self, image_path):
        """

        """
        self.image_path =  image_path
        self.data = np.array(fits.getdata(image_path),dtype='Float64')
        self.header = fits.getheader(image_path)

    def centroid(self, guess_center, delta=10.):
        img = self.data[int(guess_center[1]-delta):int(guess_center[1]+delta),int(guess_center[0]-delta):int(guess_center[0]+delta)]
        center = np.unravel_index(np.argmax(img), img.shape)

        if center[0] < delta:
            new_X = int(guess_center[0] - (delta - center[0]))
        else:
            new_X = int(guess_center[0] + (center[0]-delta))

        if center[1] < delta:
            new_Y = int(guess_center[1] - (delta - center[1]))
        else:
            new_Y = int(guess_center[1] + (center[1]-delta))

        return new_X, new_Y

    def sources_field(self):
        bkg_sigma = 1.48 * mad(self.data)
        sources = daofind(self.data, fwhm=4.0, threshold=3 * bkg_sigma)
        return sources

    def fit(self, center, delta=10., model='gaussian',show=False):
     # PSF Fitting
        '''
        Fitting a PSF model to a column choosing between the Gaussian or pseudo-Voigt profile.
        '''
        counts = self.data[int(center[1]),int(center[0]-delta):int(center[0]+delta)]
        rows = np.arange(0,len(counts),1)
        amp = np.max(counts)/2.
        cen = np.mean(rows)
        wid = 1.

        def gaussian(x, amp, cen, wid):
            '''
            gaussian PSF
            '''
            return amp * np.exp(-(x-cen)**2 /wid)
    
        def pvoigt(x, amp, cen, sigma, W, B):
            '''
            pseudo-Voigt profile
            PSF model based on Kyle's paper
            '''
            return (1-W) * amp * np.exp(-(x - cen)**2/(2.*sigma)) + 1. * W  * (amp*sigma**2)/((x-cen)**2 +sigma**2) + B

        if model == 'gaussian':
            gmodel = lmfit.Model(gaussian)
            result = gmodel.fit(counts, x=rows, amp=amp, cen=cen, wid=wid)
        
        if model == 'pVoigt':
            gmodel = lmfit.Model(pvoigt)
            result = gmodel.fit(counts, x=rows, amp=amp, cen=cen, sigma = wid/2., W = 1., B = 1.)

        if show == True:
            print(result.fit_report())

        return result

    def fwhm(self,sigma):
        return

class photometry(object):


    def __init__(self,image_path):
        self.image_path = image_path
        self.data = np.array(fits.getdata(image_path),dtype='Float64')
        self.header = fits.getheader(image_path)

    def background(self,sky,window=100):
        sky_mean = float(np.median(self.data[int(sky[1]-window):int(sky[1]+window),int(sky[0]-window):int(sky[0]+window)]))
        sky_size = self.data.shape
        return np.random.poisson(sky_mean,sky_size)

    def aperture(self,positions,sky,radius=3.,window=100):

        flux, eflux = np.zeros(len(positions)), np.zeros(len(positions))

        apertures = CircularAperture(positions, r=radius)
        bkg = background(self,sky,window=window)
        phot_table = aperture_photometry(img, apertures, error=bkg)

        for i in range(len(positions)):
            flux[i] = phot_table['aperture_sum'][i]
            flux_err[i] = phot_table['aperture_sum_err'][i]

        frames = [pd.DataFrame(flux).T, pd.DataFrame(eflux).T]
        data_flux = pd.concat(frames)
        data_flux.columns = ['flux','eflux']

        return data_flux

    def airmass(self):
        if 'airmass' in self.header:
            return float(self.header['airmass'])