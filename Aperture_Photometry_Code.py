# -*- coding: utf-8 -*-
"""
Created on Mon May 25 16:35:45 2020

@author: billw
"""

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from Background_Reduction import band_classification, Reduced_Image_Data, Bkg, Bkg_Sigma, Median, Std
import csv
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from numpy import mean 
import numpy as np
from photutils import aperture_photometry, CircularAperture
from photutils import DAOStarFinder
from scipy import spatial
from scipy.stats import zscore
from standard_star_file import std_star_file_processing
import statistics





Header_Data_Dictionary, R_Band_Files, I_Band_Files, Z_Band_Files, G_Band_Files, U_Band_Files =  band_classification()

RA, DEC, Position, Umag, Err_Umag, Gmag, Err_Gmag, Rmag, Err_Rmag, Imag, Err_Imag, Zmag, Err_Zmag = std_star_file_processing('Standard_Star_File.csv')




def aperture_main(Band, Band_Mag, Band_Mag_Err, image_num, Letter):
    
        
    Aperture_Sum, Mags, Mag_Err, X_coord, Y_coord, Phot_Table, DAO_Flux = source_detection(Band, Reduced_Image_Data, Bkg, Bkg_Sigma, Median, Std, image_num)
        
    catalog_matches, c_matches, Cata_RA, Cata_DEC, WX, WY = source_pos_comparison(Band, X_coord, X_coord, Y_coord, RA, DEC, image_num, Letter)
    
    write_out_data(Band, image_num, WX, WY, Aperture_Sum, Mags, Mag_Err, Letter)
    
    Mag_Cat, Mag_Err_Cat, Mag_WCS, Mag_Err_WCS, WCS_RA_DEC, Source_Matching = cross_correlation(Position, Band_Mag, Band_Mag_Err, Mags, Mag_Err, WX, WY, Cata_RA, Cata_DEC, image_num)
    
    slope, intercept, slope_error, intercept_error = calibration_stats(Mag_WCS, Mag_Err_WCS, Mag_Cat, Mag_Err_Cat, Letter)
    
    
    return X_coord, Y_coord, WX, WY, Mags, Mag_Err, WCS_RA_DEC, slope, intercept, slope_error, intercept_error, Source_Matching, Mag_Cat, Mag_Err_Cat, DAO_Flux
    


def source_detection(Band, Reduced_Image_Data, Bkg, Bkg_Sigma, Median, Std, image_num):
    
    
    Find_Source = DAOStarFinder(fwhm = 8, threshold = 3. * Std[image_num]) 
    Sources = Find_Source(Reduced_Image_Data[image_num] - Median[image_num]) # Finding the sources within the defined parameters in DAOStarFinder
        
    Positions = np.transpose((Sources['xcentroid'], Sources['ycentroid']))  
    DAO_Flux = np.array(Sources['flux'])
    Apertures = CircularAperture(Positions, r = 5.)  # Chossing an aperture radius of 5 pixels
    Phot_Table = aperture_photometry(Reduced_Image_Data[image_num], Apertures)
    
    X_coord = np.array(Phot_Table['xcenter'])
    Y_coord = np.array(Phot_Table['ycenter'])
        
    Aperture_Sum = Phot_Table['aperture_sum']
    Flux_Err = 15*0.5*Std[image_num] 
        
     
    Mags = -2.5*np.log10(Aperture_Sum) # Convert aperture sums to magnitudes
    Mag_Err = ((-2.5*Flux_Err)/(np.log10(10)*Aperture_Sum))
    
        
    
    #plt.figure()   # This block of code will plot the identified sources
    #plt.imshow(Reduced_Image_Data[image_num], cmap = 'gray_r', origin = 'lower')
    #Apertures.plot(color = 'blue', lw = 1.5, alpha = 1.5)
    #print(f'Time of Obs : {Header_Data_Dictionary[Band[image_num]][0]}')
    #print(f'Band : {Header_Data_Dictionary[Band[image_num]][2]}')
    
    
    return Aperture_Sum, Mags, Mag_Err, X_coord, Y_coord, Phot_Table, DAO_Flux



def source_pos_comparison(List, coord_loop_len, X_coord, Y_coord, RA, DEC, image_num, Letter):
    
    Cata_RA = []
    Cata_DEC = []
    
        
    header = fits.getheader(f'RAW_DATA/{List[image_num]}')
    w = WCS(header)
        
    WX, WY = w.wcs_pix2world(np.array(X_coord), np.array(Y_coord), 0)

        
    c = SkyCoord(RA*u.deg, DEC*u.deg)
    catalog = SkyCoord(WX*u.deg, WY*u.deg)
    
    max_sep = 1 * u.arcsec # Defining the limit in which to match stars to the catalog.
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
        
            
    
    with open(f'CATALOGS/' + f'{Letter}_{image_num}_my_catalog.tsv', 'w') as out_file:
            
        tsv_writer = csv.writer(out_file, delimiter = '\t')
            
        for n in range(len(Cata_RA)):
                
            tsv_writer.writerow([Cata_RA[n], Cata_DEC[n]])
                    
                
    return  catalog_matches, c_matches, Cata_RA, Cata_DEC, WX, WY



def write_out_data(Band, image_num, WX, WY, Aperture_Sum, Mags, Mag_Err, Letter):
    
    with open(f'FIT_DATA/{Letter}_Data_{image_num}.csv', 'w') as out_file:
            
        csv_writer = csv.writer(out_file, delimiter = ',')
        csv_writer.writerow(['Label', 'RA', 'DEC', 'Flux', 'Mag', 'Err_Mag'])
            
        for n in range(len(Mags)):
            
            csv_writer.writerow([f'Star{n}', WX[n], WY[n], Aperture_Sum[n], Mags[n], Mag_Err[n]])
                   
    return            
            
            
            

def cross_correlation(Position, Band, Band_Err, Mags, Mag_Err, WX, WY, Cata_RA, Cata_DEC, image_num):
    
    Pos_Mag_Idx = []
    WCS_Mag_Idx = []
    
    Mag_Cat = []
    Mag_Err_Cat = []
    
    Mag_WCS = []
    Mag_Err_WCS = []
    
    
    Cata_RA_DEC = [(Cata_RA[i], Cata_DEC[i]) for i in range(0, len(Cata_RA))]
    
    WCS_RA_DEC = [(WX[i], WY[i]) for i in range(0, len(WX))] 
    

    Tree_Cat = spatial.KDTree(Position)
    
    for i in range(len(Cata_RA_DEC)):
        
        Idx = Tree_Cat.query(Cata_RA_DEC[i])
        Pos_Mag_Idx.append(Idx[1])
        
    
    for j in range(len(Pos_Mag_Idx)):
        
        Mag_Cat.append(Band[Pos_Mag_Idx[j]])
        Mag_Err_Cat.append(Band_Err[Pos_Mag_Idx[j]])
        
        
    Tree_WCS = spatial.KDTree(WCS_RA_DEC)
    
    for k in range(len(Cata_RA_DEC)):
        
        Idx_WCS = Tree_WCS.query(Cata_RA_DEC[k])
        WCS_Mag_Idx.append(Idx_WCS[1])
    
    
    for l in range(len(WCS_Mag_Idx)):
        
        Mag_WCS.append(Mags[WCS_Mag_Idx[l]]) 
    
    for m in range(len(WCS_Mag_Idx)):
        
        Mag_Err_WCS.append(Mag_Err[WCS_Mag_Idx[m]]) 
        
        
    
    Source_Matching = [(WCS_Mag_Idx[i], Pos_Mag_Idx[i]) for i in range(0, len(WCS_Mag_Idx))] 
            
    
    return np.array(Mag_Cat), np.array(Mag_Err_Cat), np.array(Mag_WCS), np.array(Mag_Err_WCS), WCS_RA_DEC, Source_Matching
  
 

def calibration_stats(Mag_WCS, Mag_Err_WCS, Mag_Cat, Mag_Err_Cat, Letter):
    
    Mag_Cat = Mag_Cat[np.logical_not(np.isnan(Mag_WCS))]
    Mag_Err_Cat = Mag_Err_Cat[np.logical_not(np.isnan(Mag_WCS))]
    Mag_Err_WCS = Mag_Err_WCS[np.logical_not(np.isnan(Mag_WCS))]
    Mag_WCS = Mag_WCS[np.logical_not(np.isnan(Mag_WCS))] 

    indices = np.where(np.absolute(zscore(Mag_WCS)) > 5)[0]
    indices_filter = [i for i,n in enumerate(Mag_WCS) if i not in indices]


    Mag_WCS = Mag_WCS[indices_filter]
    Mag_Err_WCS = Mag_Err_WCS[indices_filter]
    Mag_Cat = Mag_Cat[indices_filter]
    Mag_Err_Cat = Mag_Err_Cat[indices_filter]

    plt.figure()
    plt.plot(Mag_WCS, Mag_Cat, 'o')
    plt.errorbar(Mag_WCS, Mag_Cat, yerr = np.array(Mag_Err_Cat), xerr = np.array(Mag_Err_WCS), linestyle = 'None')
    linear_coeff, cov = np.polyfit(Mag_WCS, Mag_Cat, 1, w = 1/np.array(Mag_Err_Cat), cov = True)
    slope = linear_coeff[0]
    intercept = linear_coeff[1]
    fit_error = np.sqrt(np.diagonal(cov))
    slope_error = fit_error[0]
    intercept_error = fit_error[1]
    #x_new = np.linspace(-12, -5, 5000)
    x_new = np.linspace(-15.5, -8, 5000)
    y_new = np.poly1d(linear_coeff)
    plt.plot(x_new, y_new(x_new), label = 'Line of Best Fit')
    plt.xlabel('Instrumental Magnitude')
    plt.ylabel('Standard Star File Magnitude')
    plt.title(f'Photometric Calibration Plot Using {Letter} Band Files')
    plt.legend(loc = 'best')

    Calibrated_Mag = slope * Mag_WCS + intercept
    Mag_Diff = (Mag_Cat - Calibrated_Mag)*10
    Mean_Diff = mean(Mag_Diff)
    Median_Diff = statistics.median(Mag_Diff)
    Stdev_Diff = statistics.stdev(Mag_Diff, Mean_Diff)

    plt.figure()
    plt.hist(Mag_Diff, bins = 'auto', facecolor = 'r', edgecolor = 'k', alpha = 0.75)
    plt.xlabel('Percentage [%]')
    plt.ylabel('Number')
    plt.title('Histogram of Percentage Differnce in Calibrated Stars vs Standard Star Magnitudes')
    plt.text(1, 9, f'$\mu$ = {Mean_Diff} % \n \n $\sigma$ = {Stdev_Diff} % \n \n Median = {Median_Diff} % ')
    plt.axvline(Median_Diff, color = 'k', linestyle = 'dashed', linewidth = 2)
    plt.show() 
    
    return slope, intercept, slope_error, intercept_error


#X_coord, Y_coord, WX, WY, Mags, Mag_Err, WCS_RA_DEC, slope, intercept, slope_error, intercept_error, Source_Matching, Mag_Cat, Mag_Err_Cat, DAO_Flux = aperture_main(Z_Band_Files, Zmag, Err_Zmag, 0, 'Z')



def main_aperture_data(Band, SSF_Mag, SSF_Mag_Err, Letter):
    
    
    X_coord_ = []
    Y_coord_ = []
    WX_ = []
    WY_ = []
    Mags_ = []
    Mag_Err_ = []
    WCS_RA_DEC_ = []
    slope_ = []
    intercept_ = [] 
    slope_error_ = [] 
    intercept_error_ = [] 
    Source_Matching_ = []
    Mag_Cat_ = []
    Mag_Err_Cat_ = []
    DAO_Flux_ =  []
        
    
          
    for i in range(len(Band)):
            
        X_coord, Y_coord, WX, WY, Mags, Mag_Err, WCS_RA_DEC, slope, intercept, slope_error,\
        intercept_error, Source_Matching, Mag_Cat, Mag_Err_Cat, DAO_Flux = aperture_main(Band, SSF_Mag, SSF_Mag_Err, i, Letter)
    
                
        X_coord_.append(X_coord)
        Y_coord_.append(Y_coord)
        WX_.append(WX)
        WY_.append(WY)
        Mags_.append(np.array(Mags))
        Mag_Err_.append(np.array(Mag_Err))
        WCS_RA_DEC_.append(np.array(WCS_RA_DEC))
        slope_.append(np.array(slope))
        intercept_.append(np.array(intercept))
        slope_error_.append(np.array(slope_error))
        intercept_error_.append(np.array(intercept_error))
        Source_Matching_.append(Source_Matching)
        Mag_Cat_.append(Mag_Cat) 
        Mag_Err_Cat_.append(Mag_Err_Cat)
        DAO_Flux_.append(DAO_Flux)
        
    print('Aperture Photometry completed and all varibles have been saved')
    
        
    return X_coord_, Y_coord_, WX_, WY_, Mags_, Mag_Err_, WCS_RA_DEC_, slope_, intercept_, slope_error_, intercept_error_, Source_Matching_, Mag_Cat_, Mag_Err_Cat_, DAO_Flux_

  
X_coord_, Y_coord_, WX_, WY_, Mags_, Mag_Err_, WCS_RA_DEC_, slope_, intercept_, slope_error_, intercept_error_, Source_Matching_, Mag_Cat_, Mag_Err_Cat_, DAO_Flux_ = main_aperture_data(Z_Band_Files, Zmag, Err_Zmag, 'Z')
    

##################################################################
# This script performs aperture photometry on a set of exposure
# fits files. This sciprt requires the background subtracted 
# exposures, the standard star magnitudes for the band your 
# running. The code is quite automatic, you can vary the aperture 
# radius, threshold and FWHM at your discretion.
##################################################################