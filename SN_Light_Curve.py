# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 14:29:07 2020

@author: Bill Kamtchou
"""

from Aperture_Photometry_Code import Mags_, Mag_Err_, WCS_RA_DEC_, slope_, intercept_, slope_error_, intercept_error_
from Background_Reduction import band_classification
import csv
import matplotlib.pyplot as plt
import numpy as np
from PSF_Photometry import PSF_Slope, PSF_Intercept, PSF_Slope_Err, PSF_Intercept_Err, PSF_Mags_1_, PSF_Mags_Err_1_, WCS_RA_DEC_PSF_
from scipy import spatial
from standard_star_file import std_star_file_processing



Header_Data_Dictionary, R_Band_Files, I_Band_Files, Z_Band_Files, G_Band_Files, U_Band_Files =  band_classification()

RA, DEC, Position, Umag, Err_Umag, Gmag, Err_Gmag, Rmag, Err_Rmag, Imag, Err_Imag, Zmag, Err_Zmag = std_star_file_processing('Standard_Star_File.csv')


R_Band_Zero_Flux = 2.779E-9
I_Band_Zero_Flux = 1.851E-9 
Z_Band_Zero_Flux = 1.316E-9
G_Band_Zero_Flux = 4.679E-9
U_Band_Zero_Flux = 8.609E-9


def main_lc(Band, Zero_Flux, Letter):
    
    AP_LC_Mag, AP_Calibrated_Mag_Error, PSF_LC_Mag, PSF_Mag_Error, MJD = get_lc_data(Band)
    
    x, y, y_err = plot_lc(AP_LC_Mag, AP_Calibrated_Mag_Error, PSF_LC_Mag, PSF_Mag_Error, MJD, Letter)
    
    AP_True_Mag, PSF_True_Mag = bkg_reduced_lc(AP_LC_Mag, AP_Calibrated_Mag_Error, PSF_LC_Mag, PSF_Mag_Error, MJD, x, y, y_err, Zero_Flux, Letter)
    
    write_out_lc_data(AP_LC_Mag, AP_Calibrated_Mag_Error, PSF_LC_Mag, PSF_Mag_Error, AP_True_Mag, PSF_True_Mag, MJD, Letter)
    
    return AP_LC_Mag, AP_Calibrated_Mag_Error, PSF_LC_Mag, PSF_Mag_Error, AP_True_Mag, PSF_True_Mag, MJD, x, y, y_err
    
    

def get_lc_data(Band):
    
    AP_LC_Mag = []
    AP_Calibrated_Mag_Error = []
    PSF_LC_Mag = []
    PSF_Mag_Error = []
    MJD = []
    SN = [(208.511333, 14.729711)]
        
    
    for i in range(len(Band)):
    
        Tree_SN = spatial.KDTree(WCS_RA_DEC_[i])
        Idx = Tree_SN.query(SN)
        
        AP_Calibrated_Mag = (slope_[i] * Mags_[i][Idx[1][0]]) + intercept_[i]
        AP_Calibrated_Mag_Error_1 = np.sqrt((Mags_[i][Idx[1][0]]*slope_error_[i])**2 + (slope_[i]*Mag_Err_[i][Idx[1][0]])**2 + intercept_error_[i]**2)
        
        AP_Calibrated_Mag_Error.append(AP_Calibrated_Mag_Error_1)
        AP_LC_Mag.append(AP_Calibrated_Mag)
        
        
    for i in range(len(Band)):
        
        Tree_SN = spatial.KDTree(WCS_RA_DEC_PSF_[i])
        Idx = Tree_SN.query(SN)
        
        PSF_Calibrated_Mag = (PSF_Slope[i] * PSF_Mags_1_[i][Idx[1][0]]) + PSF_Intercept[i]
        PSF_Calibrated_Mag_Error_1 = np.sqrt((PSF_Mags_1_[i][Idx[1][0]]*PSF_Slope_Err[i])**2 + (PSF_Slope[i]*PSF_Mags_Err_1_[i][Idx[1][0]])**2 + (PSF_Intercept_Err[i]**2))

        PSF_Mag_Error.append(PSF_Calibrated_Mag_Error_1)
        PSF_LC_Mag.append(PSF_Calibrated_Mag)
        MJD.append(Header_Data_Dictionary[Band[i]][1])
        
        
            
    return AP_LC_Mag, AP_Calibrated_Mag_Error, PSF_LC_Mag, PSF_Mag_Error, MJD



def plot_lc(AP_LC_Mag, AP_Calibrated_Mag_Error, PSF_LC_Mag, PSF_Mag_Error, MJD, Letter):
    
    x, y, y_err = np.loadtxt('lightcurve_ZTF.txt', unpack = True, usecols = (0,1,2))
    
    plt.figure(1)
    plt.scatter(MJD, AP_LC_Mag, label = 'AP Calibrated Mags', color = 'green')
    plt.scatter(MJD, PSF_LC_Mag, label = 'PSF Calibrated Magnitudes', color = 'blue')
    #plt.scatter(x, y, label = 'ZTF Data', color = 'red')
    plt.errorbar(MJD, AP_LC_Mag, yerr = np.array(AP_Calibrated_Mag_Error), linestyle = 'None', color = 'black', elinewidth = 1)
    plt.errorbar(MJD, PSF_LC_Mag, yerr = np.array(PSF_Mag_Error), linestyle = 'None', color = 'purple', elinewidth = 1)
    #plt.errorbar(x, y, yerr = y_err, linestyle = 'None', color = 'orange', elinewidth = 1)
    plt.legend(loc = 'best')
    plt.xlabel('Modified Julian Date [MJD]')
    plt.ylabel('Calibrated Magnitude')
    plt.title(f'Light Curve of SN2019cri Using {Letter} Band Exposures')
    plt.gca().invert_yaxis()
    #plt.ylim(21,17.0)
    plt.show()  
    
    return x, y, y_err


def bkg_reduced_lc(AP_LC_Mag, AP_Calibrated_Mag_Error, PSF_LC_Mag, PSF_Mag_Error, MJD, x, y, y_err, Zero_Flux, Letter):
    
    My_Converted_AP_Flux = 10 ** ((np.array(AP_LC_Mag) - np.array(Zero_Flux)) / -2.5)
    My_Converted_PSF_Flux = 10 ** ((np.array(PSF_LC_Mag) - np.array(Zero_Flux)) / -2.5)
    ZTF_Converted_Flux = 10 ** ((np.array(y) - np.array(Zero_Flux)) / -2.5)
    
    Max_ZTF = len(My_Converted_AP_Flux) + 7
    AP_Flux_Diff = np.nanmean(np.array(ZTF_Converted_Flux[7:Max_ZTF] - My_Converted_AP_Flux))
    PSF_Flux_Diff = np.nanmean(np.array(ZTF_Converted_Flux[7:Max_ZTF] - My_Converted_PSF_Flux))

    AP_Minus_Galaxy_Flux = My_Converted_AP_Flux - AP_Flux_Diff
    PSF_Minus_Galaxy_Flux = My_Converted_PSF_Flux - PSF_Flux_Diff
    AP_True_Mag = -2.5*np.log10(AP_Minus_Galaxy_Flux)
    PSF_True_Mag = -2.5*np.log10(PSF_Minus_Galaxy_Flux)
    
    plt.figure(2)
    plt.scatter(MJD, AP_True_Mag, label = 'AP Bkg Reduced Mags', color = 'green')
    plt.errorbar(MJD, AP_True_Mag, yerr = np.array(AP_Calibrated_Mag_Error), linestyle = 'None', color = 'black', elinewidth = 1)
    plt.scatter(MJD, PSF_True_Mag, label = 'PSF Bkg Reduced Magnitudes', color = 'blue')
    plt.errorbar(MJD, PSF_True_Mag, yerr = np.array(PSF_Mag_Error), linestyle = 'None', color = 'purple', elinewidth = 1)
    #plt.scatter(x,y, label = 'ZTF Data', color = 'red')
    #plt.errorbar(x, y, yerr = y_err, linestyle = 'None', color = 'black', elinewidth = 1)
    plt.legend(loc = 'best')
    plt.xlabel('Modified Julian Date [MJD]')
    plt.ylabel('Calibrated Magnitude')
    plt.title(f'Light Curve of SN2019cri Using Host Galaxy Subtracted {Letter} Band Exposures')
    plt.gca().invert_yaxis()
    #plt.ylim(21,17)
    plt.show()  
    
    return AP_True_Mag, PSF_True_Mag



def write_out_lc_data(AP_LC_Mag, AP_Calibrated_Mag_Error, PSF_LC_Mag, PSF_Mag_Error, AP_True_Mag, PSF_True_Mag, MJD, Letter):
    
    with open(f'LC_DATA/{Letter}_Band_LC.csv', 'w') as out_file:
            
        csv_writer = csv.writer(out_file, delimiter = ',')
        csv_writer.writerow(['Exposure Num', 'AP_LC_Mag', 'AP_Calibrated_Mag_Error', 'PSF_LC_Mag', 'PSF_Mag_Error', 'AP_Bkg_Mag', 'PSF_Bkg_Mag', 'MJD'])
            
        for n in range(len(AP_LC_Mag)):
            
            csv_writer.writerow([f'Exposure{n}', AP_LC_Mag[n], AP_Calibrated_Mag_Error[n], PSF_LC_Mag[n], PSF_Mag_Error[n], AP_True_Mag[n], PSF_True_Mag[n], MJD[n]])
                   
    return  


AP_LC_Mag, AP_Calibrated_Mag_Error, PSF_LC_Mag, PSF_Mag_Error, AP_True_Mag, PSF_True_Mag, MJD, x, y, y_err = main_lc(Z_Band_Files, Z_Band_Zero_Flux, 'Z')


##################################################################
# This script finds the magnitude of the supernova at its given
# coordinates and then calibrates it for each exposure for both
# aperture and PSF photometry methods. The results is then plotted.
# bkg_reduced_lc is only to be used to subtract galaxy background.
# All the relevant data is then written out to a csv file.
##################################################################