import pyvo as vo
import csv
import pathlib
import tarfile
from scipy.ndimage import median_filter
from astropy.utils.data import download_file
import matplotlib.pyplot as plt
import shutil
from astroquery.eso import Eso
import pandas as pd
import requests
from astropy.table import Table
from astropy.io import fits
import os
import numpy as np
from scipy.signal import find_peaks
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
import astropy.units as u
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d



def plot_spectra():


    flux = fits.open("lte06500-4.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits")[0].data
    wavelength = fits.open("WAVE_PHOENIX-ACES-AGSS-COND-2011.fits")[0].data
    

    os.makedirs("synthetic_spectra", exist_ok=True)
    plt.plot(wavelength, flux, linewidth=1)
    plt.xlabel("Wavelength")
    plt.ylabel("Flux")
    plt.xlim(4800, 5800)
    plt.show()
    plt.savefig(f"synthetic_spectra/testA.png")
    plt.close()



def normalise_spectra():

    
    flux = fits.open("lte06500-4.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits")[0].data

    print(f"Pre norm {flux.mean()}")
    
    wavelength = fits.open("WAVE_PHOENIX-ACES-AGSS-COND-2011.fits")[0].data

    cut = (wavelength > 4700) & (wavelength < 6900)
    wavelength = wavelength[cut]
    flux = flux[cut]

    
    fitter = fitting.FittingWithOutlierRemoval(fitting.LinearLSQFitter(), sigma_clip, niter=10, sigma_lower=1.5, sigma_upper=3)
    continuum, _ = fitter(models.Chebyshev1D(10, domain=[wavelength.min(), wavelength.max()]), wavelength, flux)

    flux = flux/continuum(wavelength)

    print(f"Post norm {flux.mean()}")

    os.makedirs("synthetic_spectra", exist_ok=True)
    plt.plot(wavelength, flux, linewidth=1)
    plt.xlabel("Wavelength")
    plt.ylabel("Flux")
    plt.ylim(0, 1.4)
    plt.savefig(f"synthetic_spectra/testA_normalised.png")
    plt.show()
    plt.close()




def degrade_spectra():

    target_res = 66320

    flux = fits.open("lte06500-4.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits")[0].data
    wavelength = fits.open("WAVE_PHOENIX-ACES-AGSS-COND-2011.fits")[0].data
    base_res = 500000


##    data = pd.read_csv("bosz2024_mp_t6500_g+4.0_m+0.00_a+0.00_c+0.00_v1_rorig_noresam.txt.gz", sep=r"\s+", header=None, names=["wavelength", "flux", "continuum"])
##    data.to_csv("bosz_t6500_g4.0_m0.00_cut.csv", index=False)
##    bosz_csv = pd.read_csv("bosz_t6500_g4.0_m0.00_cut.csv")
##    flux = bosz_csv["flux"]
##    wavelength = bosz_csv["wavelength"]
##    base_res = 353936

    header=fits.open("lte06500-4.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits")[0].header

    print(header)


    cut = (wavelength > 4700) & (wavelength < 6900)
    wavelength = wavelength[cut]
    flux = flux[cut]

    log_wave_grid = np.linspace(np.log(wavelength.min()), np.log(wavelength.max()), len(wavelength))
    
    flux_log = interp1d(wavelength, flux, kind='linear', fill_value="extrapolate")(np.exp(log_wave_grid))
    
    d_log_lambda = np.mean(np.diff(log_wave_grid))
    sigma_log = (
        np.sqrt(1 / target_res**2 - 1 / base_res**2)
        / (2 * np.sqrt(2 * np.log(2)))
    )
    sigma_pixels = sigma_log / d_log_lambda
    
    degraded_flux_log = gaussian_filter(flux_log, sigma=sigma_pixels)
    
    degraded_flux = interp1d(np.exp(log_wave_grid), degraded_flux_log, kind='linear', fill_value="extrapolate")(wavelength)
    
    flux = degraded_flux
    
    normalised = np.zeros_like(flux)

    fitter = fitting.FittingWithOutlierRemoval(fitting.LinearLSQFitter(), sigma_clip, niter=10, sigma_lower=1.5, sigma_upper=3)
    continuum, _ = fitter(models.Chebyshev1D(10, domain=[wavelength.min(), wavelength.max()]), wavelength, flux)
    normalised = flux / continuum(wavelength)

    flux = normalised

    #flux = 1 - flux
    

    os.makedirs("synthetic_spectra", exist_ok=True)
    plt.plot(wavelength, flux, linewidth=1)
    plt.xlabel("Wavelength")
    plt.ylabel("Flux")
    plt.ylim(0, 1.4)
    plt.savefig(f"synthetic_spectra/testA_normalised.png")
    plt.show()
    plt.close()

  
    flux_hdu = fits.PrimaryHDU(data=flux)
    wavelength_hdu = fits.ImageHDU(data=wavelength, name="WAVELENGTH")

    hdul = fits.HDUList([flux_hdu, wavelength_hdu])
    hdul.writeto("synthetic_spectra/degraded_spectrum.fits", overwrite=True)



#plot_spectra()
#normalise_spectra()
degrade_spectra()
