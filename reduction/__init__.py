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

class psf(object):

    def __init__(self, image_path):
        """

        """
        self.image_path =  image_path
        self.data = np.array(fits.getdata(image_path),dtype='Float64')
        self.header = fits.getheader(image_path)

    def fit(self, center, delta=10., model='gaussian', amp = 5000, cen = 80, wid = 1., W = 1., B=1., show=False):
     # PSF Fitting
    '''
    Fitting a PSF model to a column choosing between the Gaussian or pseudo-Voigt profile.
    '''

    counts = self.data[int(center),int(center-delta):int(center+delta)]
    rows = np.arange(0,len(counts),1)

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