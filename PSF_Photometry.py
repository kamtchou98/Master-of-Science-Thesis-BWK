# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 14:01:10 2020

@author: Bill Kamtchou
"""

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.nddata import NDData
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.table import Table
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from Background_Reduction import Reduced_Image_Data, Std, R_Band_Files, I_Band_Files, Z_Band_Files, G_Band_Files, U_Band_Files
import csv
import matplotlib.pyplot as plt
from numpy import mean 
import numpy as np
from scipy import spatial
from scipy.stats import zscore
from standard_star_file import std_star_file_processing
import statistics
from photutils import DAOStarFinder
from photutils.psf import IterativelySubtractedPSFPhotometry, DAOGroup
from photutils.psf import extract_stars
from photutils import EPSFBuilder
from photutils import find_peaks




RA, DEC, Position, Umag, Err_Umag, Gmag, Err_Gmag, Rmag, Err_Rmag, Imag, Err_Imag, Zmag, Err_Zmag = std_star_file_processing('Standard_Star_File.csv')



def psf_main(Band, Band_Mag, Band_Mag_Err, image_num, Letter):
       
    PSF_Mags_1, PSF_Mags_Err_1, Flux, Flux_Err, X_0, Y_0 = do_psf(image_num, Letter) 
       
    catalog_matches, c_matches, Cata_RA, Cata_DEC, WX, WY = source_pos_comparison(Band, X_0, X_0, Y_0, RA, DEC, image_num, Letter)
    
    Mag_Cat, Mag_Err_Cat, Mag_WCS, Mag_Err_WCS, WCS_RA_DEC_PSF_, Source_Matching_PSF = cross_correlation(Position, Band_Mag, Band_Mag_Err, PSF_Mags_1, PSF_Mags_Err_1, WX, WY, Cata_RA, Cata_DEC, image_num)
    
    write_out_data(image_num, WX, WY, Flux, Flux_Err, PSF_Mags_1, PSF_Mags_Err_1, Letter)
    
    slope_, intercept_, slope_error_, intercept_error_ = psf_calibration_stats(Mag_WCS, Mag_Err_WCS, Mag_Cat, Mag_Err_Cat, Letter)
    
    
    return slope_, intercept_, slope_error_, intercept_error_, PSF_Mags_1, PSF_Mags_Err_1, WCS_RA_DEC_PSF_



def build_epsf(image_num):
    
    size = 7
    hsize = (size - 1) / 2
    peaks_tbl = find_peaks(Reduced_Image_Data[image_num], threshold = 750.)
    x = peaks_tbl['x_peak']
    y = peaks_tbl['y_peak']  
    mask = ((x > hsize) & (x < (Reduced_Image_Data[image_num].shape[1] -1 - hsize)) & (y > hsize) & (y < (Reduced_Image_Data[image_num].shape[0] -1 - hsize)))  

    stars_tbl = Table()
    stars_tbl['x'] = x[mask]  
    stars_tbl['y'] = y[mask] 
   
    nddata = NDData(data = Reduced_Image_Data[image_num])  

    stars = extract_stars(nddata, stars_tbl, size = 10)

    #nrows = 5
    #ncols = 5
    #fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize = (15, 15), squeeze = True)
    #ax = ax.ravel()

    #for i in range(nrows*ncols):
        #norm = simple_norm(stars[i], 'log', percent = 99.)
        #ax[i].imshow(stars[i], norm = norm, origin = 'lower', cmap = 'gray')

    epsf_builder = EPSFBuilder(oversampling = 8, maxiters = 20, progress_bar = False)  
    epsf, fitted_stars = epsf_builder(stars) 

    #norm = simple_norm(epsf.data, 'log', percent = 99.)
    #plt.imshow(epsf.data, norm = norm, origin = 'lower', cmap = 'viridis')
    #plt.colorbar()
    
    return epsf


def do_psf(image_num, Letter):
        
    epsf = build_epsf(image_num)
    
    daofind = DAOStarFinder(fwhm = 8, threshold = 3. * Std[image_num])
    
    sigma_psf = 2.0
    daogroup = DAOGroup(2.0*sigma_psf*gaussian_sigma_to_fwhm)

   
    photometry = IterativelySubtractedPSFPhotometry(finder = daofind,
                                                    group_maker = daogroup,
                                                    bkg_estimator = None,
                                                    psf_model = epsf,
                                                    fitter = LevMarLSQFitter(),
                                                    aperture_radius = 5,
                                                    niters = 1, fitshape = (11,11))

    result_tab = photometry(image = Reduced_Image_Data[image_num])
    #residual_image = photometry.get_residual_image()
    
    X_0 = result_tab['x_fit']
    Y_0 = result_tab['y_fit']
    Flux = result_tab['flux_fit']
    Flux_Err = np.sqrt(np.array(Flux))
    PSF_Mags_1 = -2.5*np.log10(np.array(Flux)) # Convert aperture sums to magnitudes 
    PSF_Mags_Err_1 = ((-2.5*Flux_Err)/(np.log10(10)*Flux))
    

    return PSF_Mags_1, PSF_Mags_Err_1, Flux, Flux_Err, X_0, Y_0



def source_pos_comparison(List, coord_loop_len, X_coord, Y_coord, RA, DEC, image_num, Letter):
    
    Cata_RA = []
    Cata_DEC = []
        
    header = fits.getheader(f'RAW_DATA/{List[image_num]}')
    w = WCS(header)
        
    WX, WY = w.wcs_pix2world(np.array(X_coord), np.array(Y_coord), 0)

        
    c = SkyCoord(RA*u.deg, DEC*u.deg)
    catalog = SkyCoord(WX*u.deg, WY*u.deg)
    
    max_sep = 1 * u.arcsec
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)
    sep_constraint = d2d < max_sep
    c_matches = c[sep_constraint]
    catalog_matches = catalog[idx[sep_constraint]]

            
    Cata_List = list(catalog_matches)
    
        
    for m in range(len(catalog_matches)):
        
        Cata_Val = str(Cata_List[m])[40 : -3]
        Cata_Array = [float(idx) for idx in Cata_Val.split(', ')]
        Cata_RA_ta = Cata_Array[0]
        Cata_DEC_ta = Cata_Array[1]
    
        Cata_RA.append(Cata_RA_ta)
        Cata_DEC.append(Cata_DEC_ta)
        
            
    
    with open(f'CATALOGS/' + f'PSF_{Letter}_{image_num}_my_catalog.tsv', 'w') as out_file:
            
        tsv_writer = csv.writer(out_file, delimiter = '\t')
            
        for n in range(len(Cata_RA)):
                
            tsv_writer.writerow([Cata_RA[n], Cata_DEC[n]])
                    
                
    return  catalog_matches, c_matches, Cata_RA, Cata_DEC, WX, WY
            
            

def cross_correlation(Position, Band, Band_Err, Mags, Mag_Err, WX, WY, Cata_RA, Cata_DEC, image_num):
    
    Pos_Mag_Idx = []
    WCS_Mag_Idx = []
    
    Mag_Cat = []
    Mag_Err_Cat = []
    
    Mag_WCS = []
    Mag_Err_WCS = []
    
    
    Cata_RA_DEC = [(Cata_RA[i], Cata_DEC[i]) for i in range(0, len(Cata_RA))]
    
    WCS_RA_DEC_PSF = [(WX[i], WY[i]) for i in range(0, len(WX))] 
    

    Tree_Cat = spatial.KDTree(Position)
    
    for i in range(len(Cata_RA_DEC)):
        
        Idx = Tree_Cat.query(Cata_RA_DEC[i])
        Pos_Mag_Idx.append(Idx[1])
        
    
    for j in range(len(Pos_Mag_Idx)):
        
        Mag_Cat.append(Band[Pos_Mag_Idx[j]])
        Mag_Err_Cat.append(Band_Err[Pos_Mag_Idx[j]])
        
        
    Tree_WCS = spatial.KDTree(WCS_RA_DEC_PSF)
    
    for k in range(len(Cata_RA_DEC)):
        
        Idx_WCS = Tree_WCS.query(Cata_RA_DEC[k])
        WCS_Mag_Idx.append(Idx_WCS[1])
    
    
    for l in range(len(WCS_Mag_Idx)):
        
        Mag_WCS.append(Mags[WCS_Mag_Idx[l]]) 
    
    for m in range(len(WCS_Mag_Idx)):
        
        Mag_Err_WCS.append(Mag_Err[WCS_Mag_Idx[m]]) 
        
        
    
    Source_Matching_PSF = [(WCS_Mag_Idx[i], Pos_Mag_Idx[i]) for i in range(0, len(WCS_Mag_Idx))] 
            
    
    return np.array(Mag_Cat), np.array(Mag_Err_Cat), np.array(Mag_WCS), np.array(Mag_Err_WCS), WCS_RA_DEC_PSF, Source_Matching_PSF

    

def psf_calibration_stats(Mag_WCS, Mag_Err_WCS, Mag_Cat, Mag_Err_Cat, Letter):
    
    Mag_Cat = Mag_Cat[np.logical_not(np.isnan(Mag_WCS))]
    Mag_Err_Cat = Mag_Err_Cat[np.logical_not(np.isnan(Mag_WCS))]
    Mag_Err_WCS = Mag_Err_WCS[np.logical_not(np.isnan(np.array(Mag_WCS)))]
    Mag_WCS = Mag_WCS[np.logical_not(np.isnan(np.array(Mag_WCS)))] 

    indices = np.where(np.absolute(zscore(Mag_WCS)) > 5)[0]
    indices_filter = [i for i,n in enumerate(Mag_WCS) if i not in indices]


    PSF_Mags = Mag_WCS[indices_filter]
    PSF_Mags_Err = Mag_Err_WCS[indices_filter]
    Mag_Cat = Mag_Cat[indices_filter]
    Mag_Err_Cat = Mag_Err_Cat[indices_filter]

    plt.figure(1)
    plt.plot(PSF_Mags, Mag_Cat, 'o')
    plt.errorbar(PSF_Mags, Mag_Cat, yerr = np.array(Mag_Err_Cat), xerr = np.array(PSF_Mags_Err), linestyle = 'None')
    linear_coeff, cov = np.polyfit(PSF_Mags, Mag_Cat, 1, w = 1/np.array(Mag_Err_Cat), cov = True)
    slope_ = linear_coeff[0]
    intercept_ = linear_coeff[1]
    fit_error = np.sqrt(np.diagonal(cov))
    slope_error_ = fit_error[0]
    intercept_error_ = fit_error[1]
    x_new = np.linspace(-15, -7.5, 5000)
    y_new = np.poly1d(linear_coeff)
    plt.plot(x_new, y_new(x_new), label = 'Line of Best Fit')
    plt.xlabel('Instrumental Magnitude')
    plt.ylabel('Standard Star File Magnitude')
    plt.title(f'Photometric Calibration Plot Using {Letter} Band Files (PSF)')
    plt.legend(loc = 'best')

    Calibrated_Mag = slope_ * PSF_Mags + intercept_
    Mag_Diff = (Mag_Cat - Calibrated_Mag)*10
    Mean_Diff = mean(Mag_Diff)
    Median_Diff = statistics.median(Mag_Diff)
    Stdev_Diff = statistics.stdev(Mag_Diff, Mean_Diff)

    plt.figure(2)
    plt.hist(Mag_Diff, bins = 'auto', facecolor = 'blue', edgecolor = 'k', alpha = 0.75)
    plt.xlabel('Percentage [%]')
    plt.ylabel('Number')
    plt.title('Histogram of Percentage Differnce in Calibrated Stars vs Standard Star Magnitudes (PSF)')
    plt.text(1, 15, f'$\mu$ = {Mean_Diff} % \n \n $\sigma$ = {Stdev_Diff} % \n \n Median = {Median_Diff} % ')
    plt.axvline(Median_Diff, color = 'k', linestyle = 'dashed', linewidth = 2)
    plt.show() 
    
    return slope_, intercept_, slope_error_, intercept_error_



def write_out_data(image_num, WX, WY, Flux, Flux_Err, Mags, Mag_Err, Letter):
    
    with open(f'FIT_DATA/PSF_{Letter}_Data_{image_num}.csv', 'w') as out_file:
            
        csv_writer = csv.writer(out_file, delimiter = ',')
        csv_writer.writerow(['Label', 'RA', 'DEC', 'Flux', 'Flux_Err', 'Mag', 'Err_Mag'])
            
        for n in range(len(Mags)):
            
            csv_writer.writerow([f'Star{n}', WX[n], WY[n], Flux[n], Flux_Err[n], Mags[n], Mag_Err[n]])
    
    return
            
            

def main_psf_data(Band, SSF_Mag, SSF_Mag_Err, Letter):
    
    PSF_Slope = []
    PSF_Intercept = []
    PSF_Slope_Err = []
    PSF_Intercept_Err = []
    PSF_Mags_1_ = []
    PSF_Mags_Err_1_ = []
    WCS_RA_DEC_PSF_ = []
    
    for i in range(len(Band)):
        
    
        slope_, intercept_, slope_error_, intercept_error_, PSF_Mags_1, PSF_Mags_Err_1, WCS_RA_DEC_PSF = psf_main(Band, SSF_Mag, SSF_Mag_Err, i, Letter)    
    
    
        PSF_Slope.append(slope_)
        PSF_Intercept.append(intercept_)
        PSF_Slope_Err.append(slope_error_)
        PSF_Intercept_Err.append(intercept_error_)
        PSF_Mags_1_.append(PSF_Mags_1)
        PSF_Mags_Err_1_.append(PSF_Mags_Err_1)
        WCS_RA_DEC_PSF_.append(WCS_RA_DEC_PSF)
        
        
    print('PSF Photometry Completed')
    
    return PSF_Slope, PSF_Intercept, PSF_Slope_Err, PSF_Intercept_Err, PSF_Mags_1_, PSF_Mags_Err_1_, WCS_RA_DEC_PSF_



PSF_Slope, PSF_Intercept, PSF_Slope_Err, PSF_Intercept_Err, PSF_Mags_1_, PSF_Mags_Err_1_, WCS_RA_DEC_PSF_ = main_psf_data(Z_Band_Files, Zmag, Err_Zmag, 'Z')


##################################################################
# This script performs PSF photometry on a set of exposure
# fits files by calulating an ePSF for each exposure This sciprt 
# requires the background subtracted exposures, the standard star 
# magnitudes for the band your running. The code is quite automatic, 
# you can vary the aperture radius, threshold and FWHM on the 
# starfinder argument on line 96 at your discretion.
##################################################################

