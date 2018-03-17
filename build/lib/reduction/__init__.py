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
# from photutils import SkyAperture
from photutils import CircularAnnulus
# from photutils import SigmaClip
# from photutils import Background2D, MedianBackground

from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from astropy.stats import median_absolute_deviation as mad

class CCD(object):

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



    def sources_field(self,sky,fwhm=4.,threshold_std=3.):
        def background(self, sky, window=100):
            sky_mean = float(np.median(
                self.data[int(sky[1] - window):int(sky[1] + window), int(sky[0] - window):int(sky[0] + window)]))
            sky_size = self.data.shape
            return np.random.poisson(sky_mean, sky_size)

        bkg_sigma = np.std(background(self,sky))
        mean, median, std = sigma_clipped_stats(self.data, sigma=3.0, iters=5)
        daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold_std * std)
        sources = daofind(self.data - median)
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

    def background(self,sky,window=100):
        sky_mean = float(np.median(self.data[int(sky[1]-window):int(sky[1]+window),int(sky[0]-window):int(sky[0]+window)]))
        sky_size = self.data.shape
        return np.random.poisson(sky_mean,sky_size)

    def annulus_background(self,positions,radius=3.):
        #Obtain the sky local background
        annulus_apertures = CircularAnnulus(positions, r_in=radius+2., r_out=radius+4.)
        bkg_table = aperture_photometry(self.data, annulus_apertures)
        bkg_counts = float(bkg_table['aperture_sum'])
        bkg_mean = bkg_counts/annulus_apertures.area()
        bkg = np.random.poisson(bkg_mean,self.data.shape)
        return bkg, bkg_counts, bkg_mean

    def aperture(self,positions,bkg=None,radius=3.):

        flux, eflux = [], []

        apertures = CircularAperture(positions, r=radius)
        
        if bkg is None:
            bkg, bkg_counts, bkg_mean = self.annulus_background(positions,radius=radius)

        phot_table = aperture_photometry(self.data, apertures, error=bkg)

        for i in range(len(positions)):
            flux.append(phot_table['aperture_sum'][i])
            eflux.append(phot_table['aperture_sum_err'][i])

        frames = [pd.DataFrame(flux).T, pd.DataFrame(eflux).T]
        data_flux = pd.concat(frames,axis=1)
        # data_flux.columns = ['flux','eflux']

        return data_flux

    def airmass(self):
        if 'airmass' in self.header:
            return float(self.header['airmass'])

    def snratio(self,guess_center,sky,radius=4.,fwhm=8.,delta=10.,ND=None,NR2=None,exp_time=None,gain=None):
        
        #check the centroid
        X,Y = guess_center
        sources = self.sources_field(sky=sky,fwhm=fwhm)
        sources = sources.to_pandas()
        center = sources[(sources['ycentroid'] < int(Y+delta)) & (sources['ycentroid'] > int(Y-delta)) 
                         & (sources['xcentroid'] < int(X+delta)) & (sources['xcentroid'] > int(X-delta))]
        positions = (float(center['xcentroid']),float(center['ycentroid']))
        
        bkg, bkg_counts, bkg_mean = self.annulus_background(positions,radius=radius)
        counts = self.aperture([positions],bkg,radius=radius)
        counts.columns = ['flux','eflux']
        
        #HEADER Information
        if exp_time is None:
            exp_time = float(self.header['EXPTIME'])
        
        header_list = []
        for i in self.header.keys():
            header_list.append(i)
            
        if gain is None:
            matching = [s for s in header_list if "GAIN" in s]
            gain_list = np.zeros(len(matching))
            for i in range(len(matching)):
                gain_list[i] = self.header[matching[i]]
            gain = np.mean(gain_list)
            
        if ND is None:
            ND = float(self.header['DARKCUR'])

        if NR2 is None:
            if 'RDNOISE' in self.header.keys():
                NR2 =  float(self.header['RDNOISE'])
            else:
                matching = [s for s in header_list if "RDNOISE" in s]
                NR2_partial = np.zeros(len(matching))
                for i in range(len(matching)):
                    NR2_partial[i] = float(self.header[matching[i]])
                NR2 =  np.mean(NR2_partial)**2

        n_pix = np.pi * radius**2
        flux = float(counts['flux'])
        SN = flux/np.sqrt(flux+(n_pix*bkg_mean)+exp_time*(1.*n_pix*ND/gain)+(1.*n_pix*NR2/gain**2))
        return SN

class photometry(object):


    def __init__(self,image_path):
        self.image_path = image_path
        self.data = np.array(fits.getdata(image_path),dtype='Float64')
        self.header = fits.getheader(image_path)

    def background(self,sky,window=100):
        sky_mean = float(np.median(self.data[int(sky[1]-window):int(sky[1]+window),int(sky[0]-window):int(sky[0]+window)]))
        sky_size = self.data.shape
        return np.random.poisson(sky_mean,sky_size)

    def aperture(self,positions,bkg,radius=3.):

        flux, eflux = [], []

        apertures = CircularAperture(positions, r=radius)
        phot_table = aperture_photometry(self.data, apertures, error=bkg)

        for i in range(len(positions)):
            flux.append(phot_table['aperture_sum'][i])
            eflux.append(phot_table['aperture_sum_err'][i])

        frames = [pd.DataFrame(flux).T, pd.DataFrame(eflux).T]
        data_flux = pd.concat(frames,axis=1)
        # data_flux.columns = ['flux','eflux']

        return data_flux

    def airmass(self):
        if 'airmass' in self.header:
            return float(self.header['airmass'])