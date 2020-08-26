# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 11:14:21 2020

@author: billw
"""

from astropy.cosmology import WMAP5
from astropy.coordinates import Distance
from astropy import units as u
import math
import matplotlib.pyplot as plt
from numpy import mean 
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
import statistics
import pandas as pd 
from pylab import rcParams


def dered_mag(lameff = 0, mag = 0, E = 1e-5):
    
    # Dereddens a magnitude via the extinction law of Cardelli (1989)
	# valid for the range 1250 Angstroms to 30000 Angstroms.
	# This is the ONIR range. Assumes Rv = 3.1
	# >>>dered(3652.0, 17, 0.1)
	# >>>16.5175...
	
    Rv = 3.1
    x = 1.0/lameff
		 
    if 1.1e-4 <= x < 3.3e-4:
        	y = (x*1e4-1.82)
        	a = 1 + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
        	b = 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
    if 0.3e-4 <= x < 1.1e-4:
            a = 0.574*(x*1e4)**1.61
            b = -0.527*(x*1e4)**1.61
        	
    if 3.3e-4 <= x < 8.0e-4:
            if 5.9e-4<x<=8.0e-4:
                Fa = -0.04473*((x*1e4)-5.9)**2 - 0.009779*((x*1e4)-5.9)**3
                Fb = 0.2130*((x*1e4)-5.9)**2 + 0.1207*((x*1e4)-5.9)**3
            if x <= 5.9e-4:
                Fa=0.0
                Fb=0.0
            a = 1.752 - 0.316*(x*1e4) - 0.104*( ((x*1e4)-4.67)**2 + 0.341 )**-1 + Fa
            b = -3.090 + 1.825*(x*1e4) + 1.206*( ((x*1e4)-4.62)**2 + 0.263 ) + Fb
                    
    ratio = a + b/Rv
    dered_mag_val = mag - Rv*E*ratio
    
    return dered_mag_val


def dered(wav,flux,Rv,Ebv):
    
    # This function dereddens spectra..
    
    lam=wav*0.0001
    Av=Ebv*Rv
    x=1/lam
    y=x-1.82
    a=1+(0.17699*y)-(0.50477*y**2)-(0.02427*y**3)+(0.72085*y**4)+(0.01979*y**5)-(0.77530*y**6)+(0.32999*y**7)
    b=(1.41338*y)+(2.28305*y**2)+(1.07233*y**3)-(5.38434*y**4)-(0.62251*y**5)+(5.30260*y**6)-(2.09002*y**7)
    AlAv=a+b/Rv
    Al=AlAv*Av
    F=10**(Al/2.5)
    dered_flux = flux*F
    
    return dered_flux

def spdvega(mag, zero_flux):
    
    # Converts magnitudes to spectral flux density in cgs units.
    
    spd = zero_flux * 10**(-(np.array(mag)/2.5))
    
    return spd



def dl(dist_mod):
    
	# Returns a luminosity distance in Mpc for a given distance modulus  
    
    LD = 10**(1.0 + (dist_mod / 5.0)) / 1.0E6
    
    return LD

 
def lum(flux, DL):
    
    # Returns the luminosity in erg per sec for a given flux value in erg s^-1 cm^-2 A^-1
	# and luminosity distance in Mpc. Use Dl module to convert distance modulus to Dl.
    
    Lum = 4 * np.pi * (DL * 3.086E24)**2 * flux
    
    return Lum 



def create_sed(SN_Data, dist_mod, zero_flux):
    
    Sed = []
    
    for i in range(len(SN_Data)):
        
        Mag_Dis = Distance(z = 0.041, cosmology = WMAP5) / u.Mpc
        Abs_Mag = np.array(SN_Data[i]) - (5 * np.log10(Mag_Dis*10**6)) + 5
        
        Sed.append(Abs_Mag)
        
        
    return Sed


def get_band_data(Letter, Band_Wave):
    
    Band_LC = pd.read_csv(f'LC_DATA/{Letter}_Band_LC.csv')
    Mags_MJD = np.array(Band_LC['MJD']) * (1 / (1 + 0.041))
    
    if Letter == 'G':
        
        Mags = np.array(Band_LC['PSF_Bkg_Mag'])
        Mags_Err = np.array(Band_LC['PSF_Mag_Error'])
        
        Dered_Mags = dered_mag(Band_Wave, Mags, E = 0.02)
        
    if Letter == 'R':
        
        Mags = np.array(Band_LC['PSF_Bkg_Mag'])
        Mags_Err = np.array(Band_LC['PSF_Mag_Error'])
    
        Dered_Mags = dered_mag(Band_Wave, Mags, E = 0.02) 
    
    if Letter == 'U':
        
        Mags = np.array(Band_LC['AP_LC_Mag'])
        Mags_Err = np.array(Band_LC['AP_Calibrated_Mag_Error'])
    
        Dered_Mags = dered_mag(Band_Wave, Mags, E = 0.02)
        
    if Letter == 'I':
        
        Mags = np.array(Band_LC['PSF_LC_Mag'])
        Mags_Err = np.array(Band_LC['PSF_Mag_Error'])
    
        Dered_Mags = dered_mag(Band_Wave, Mags, E = 0.02)
    
    if Letter == 'Z':
        
        Mags = np.array(Band_LC['PSF_LC_Mag'])
        Mags_Err = np.array(Band_LC['PSF_Mag_Error'])
    
        Dered_Mags = dered_mag(Band_Wave, Mags, E = 0.02)
    

    return Dered_Mags, Mags_Err, Mags_MJD


RZTF_MJD, RZTF_Mags, RZTF_Mags_Err = np.loadtxt('lightcurve_ZTF.txt', unpack = True, usecols = (0,3,4))
GZTF_MJD, GZTF_Mags, GZTF_Mags_Err = np.loadtxt('lightcurve_ZTF.txt', unpack = True, usecols = (0,1,2))

UDered_Mags, UMags_Err, UMags_MJD = get_band_data('U', 3652)
GDered_Mags, GMags_Err, GMags_MJD = get_band_data('G', 4825)
RDered_Mags, RMags_Err, RMags_MJD = get_band_data('R', 6410)
IDered_Mags, IMags_Err, IMags_MJD = get_band_data('I', 7980)
ZDered_Mags, ZMags_Err, ZMags_MJD = get_band_data('Z', 9097)
R_ZTFDered_Mags = dered_mag(6410, RZTF_Mags, E = 0.02)
G_ZTFDered_Mags = dered_mag(4825, GZTF_Mags, E = 0.02)

UG_Last = (UDered_Mags - GDered_Mags[0 : len(UDered_Mags)])[2]
U_Add = UG_Last + GDered_Mags[12 : 21]
UMag_Comb = np.append(UDered_Mags, U_Add)
UMag_Comb_Err = np.append(UDered_Mags, GMags_Err[12 : 21])
UMJD_Comb = np.append(UMags_MJD, GMags_MJD[12 : 21])


rcParams['figure.figsize'] = 11, 9



def combine_ztf_data(Mag_ZTF, MJD_ZTF, Err_ZTF, My_Mag, My_Err, My_MJD):
    
    Comb_MJD = np.append(MJD_ZTF * (1 / (1 + 0.041)), My_MJD)
    Comb_Mag = np.append(np.array(Mag_ZTF), My_Mag)
    Comb_Mag_Err = np.append(np.array(Err_ZTF), np.array(My_Err))
    
    b = [(Comb_MJD[i], Comb_Mag[i], Comb_Mag_Err[i]) for i in range(0, len(Comb_MJD))]
    c = sorted(b, key=lambda b: b[0])

    
    True_Comb_MJD = []
    Mag_Comb = []
    Mag_Comb_Err = []
    
    for i in range(len(c)):
        
        True_Comb_MJD.append(c[i][0])
        Mag_Comb.append(c[i][1]) 
        Mag_Comb_Err.append(c[i][2])
        
    
    return Mag_Comb, Mag_Comb_Err, True_Comb_MJD


RMag_Comb, RMag_Comb_Err, R_True_Comb_MJD = combine_ztf_data(R_ZTFDered_Mags, RZTF_MJD, RZTF_Mags_Err, RDered_Mags, RMags_Err, RMags_MJD)
GMag_Comb, GMag_Comb_Err, G_True_Comb_MJD = combine_ztf_data(G_ZTFDered_Mags, GZTF_MJD, GZTF_Mags_Err, GDered_Mags, GMags_Err, GMags_MJD)
 


def add_all_ztf(Band_MJD, Band_Mag, Letter):
    
    xnew = RZTF_MJD * (1 / (1 + 0.041))
    
    if Letter == 'U':
        
        Colour_Corr_Val =  Band_Mag[12] - G_ZTFDered_Mags[1]
        
        Insert_Mags = G_ZTFDered_Mags + Colour_Corr_Val
        
        
    if Letter == 'I':
    
        Colour_Corr_Val = R_ZTFDered_Mags[12] - Band_Mag[1]
    
        Insert_Mags = R_ZTFDered_Mags + Colour_Corr_Val
        
        
    if Letter == 'Z':
    
        Colour_Corr_Val = R_ZTFDered_Mags[12] - Band_Mag[1]
    
        Insert_Mags = R_ZTFDered_Mags + Colour_Corr_Val
        
    return xnew, Insert_Mags



U_Insert_MJDs, U_Insert_Mags = add_all_ztf(UMJD_Comb, UMag_Comb, 'U')
I_Insert_MJDs, I_Insert_Mags = add_all_ztf(IMags_MJD, IDered_Mags, 'I')
Z_Insert_MJDs, Z_Insert_Mags = add_all_ztf(ZMags_MJD, ZDered_Mags, 'Z')



def arange_phot_data(OG_MJDs, OG_Mags, ta_MJDs, ta_Mags):
    
    Pre_sort_MJDs = np.append(OG_MJDs, ta_MJDs)
    Pre_sort_Mags = np.append(OG_Mags, ta_Mags)
    
    b = [(Pre_sort_MJDs[i], Pre_sort_Mags[i]) for i in range(0, len(Pre_sort_Mags))]
    c = sorted(b, key=lambda b: b[0])
    
    All_Comb_MJDs = []
    All_Comb_Mags = []
    
    for i in range(len(c)):
        
        All_Comb_MJDs.append(c[i][0])
        All_Comb_Mags.append(c[i][1]) 


    
    return All_Comb_MJDs, All_Comb_Mags 


All_Comb_UMJDs, All_Comb_UMags = arange_phot_data(UMJD_Comb, UMag_Comb, U_Insert_MJDs, U_Insert_Mags)
All_Comb_IMJDs, All_Comb_IMags = arange_phot_data(IMags_MJD, IDered_Mags, I_Insert_MJDs, I_Insert_Mags)
All_Comb_ZMJDs, All_Comb_ZMags = arange_phot_data(ZMags_MJD, ZDered_Mags, Z_Insert_MJDs, Z_Insert_Mags)




def nan_corr(Band_MJD, Band_Mag):
    
    Zero_Idxs = np.where(Band_Mag == 0.0)
    Find_NaN_Idx = np.where(np.isnan(Band_Mag))
    NaN_Idx = np.append(Find_NaN_Idx[0], Zero_Idxs[0])
    
    MJD_tf = []
    
    for i in range(len(NaN_Idx)):
        
        MJD_tf.append(Band_MJD[NaN_Idx[i]])
        
    Band_MJD_ = np.delete(Band_MJD, [NaN_Idx])
    Band_Mag_ = np.delete(Band_Mag, [NaN_Idx])
    
    
    tck = UnivariateSpline(Band_MJD_, Band_Mag_, k = 1, s = 0)
    
    #plt.figure()
    #plt.plot(Band_MJD_, Band_Mag_, '--o', color = 'blue')
    #plt.plot(Band_MJD_, tck(Band_MJD_), color = 'red')
    #plt.gca().invert_yaxis()

    Insert_Mags = []


    for j in range(len(MJD_tf)):
    
        Insert_Mag = tck(MJD_tf[j])
        Insert_Mags.append(Insert_Mag)
    
    
    Pre_sort_MJDs = np.append(Band_MJD_, MJD_tf[0 : len(MJD_tf) - 1])
    Pre_sort_Mags = np.append(Band_Mag_, Insert_Mags[0 : len(MJD_tf) - 1])
    
    b = [(Pre_sort_MJDs[i], Pre_sort_Mags[i]) for i in range(0, len(Pre_sort_Mags))]
    c = sorted(b, key=lambda b: b[0])
    
    NoNan_MJDs = []
    NoNan_Mags = []
    
    for i in range(len(c)):
        
        NoNan_MJDs.append(c[i][0])
        NoNan_Mags.append(c[i][1]) 
        
    return NoNan_MJDs, NoNan_Mags

NoNan_UMJDs, NoNan_UMags =  nan_corr(All_Comb_UMJDs, All_Comb_UMags) 
NoNan_GMJDs, NoNan_GMags =  nan_corr(G_True_Comb_MJD, GMag_Comb) 
NoNan_RMJDs, NoNan_RMags =  nan_corr(R_True_Comb_MJD, RMag_Comb) 
NoNan_IMJDs, NoNan_IMags =  nan_corr(All_Comb_IMJDs, All_Comb_IMags)
NoNan_ZMJDs, NoNan_ZMags =  nan_corr(All_Comb_ZMJDs, All_Comb_ZMags)
    


U_Flux = spdvega(NoNan_UMags, 8.609E-9) * 1.041
G_Flux = spdvega(NoNan_GMags, 4.679E-9) * 1.041
R_Flux = spdvega(NoNan_RMags, 2.779E-9) * 1.041
I_Flux = spdvega(NoNan_IMags, 1.851E-9) * 1.041
Z_Flux = spdvega(NoNan_ZMags, 1.316E-9) * 1.041



def sed():
    
    U_id = np.linspace(0, 0, len(U_Flux))
    G_id = np.linspace(1, 1, len(G_Flux))
    R_id = np.linspace(2, 2, len(R_Flux))
    I_id = np.linspace(3, 3, len(I_Flux))
    Z_id = np.linspace(4, 4, len(Z_Flux))

    U_Comb = [(np.round(NoNan_UMJDs[i], decimals = 1), U_Flux[i], U_id[i]) for i in range(0, len(U_Flux))]
    G_Comb = [(np.round(NoNan_GMJDs[i], decimals = 1), G_Flux[i], G_id[i]) for i in range(0, len(G_Flux))]
    R_Comb = [(np.round(NoNan_RMJDs[i], decimals = 1), R_Flux[i], R_id[i]) for i in range(0, len(R_Flux))]
    I_Comb = [(np.round(NoNan_IMJDs[i], decimals = 1), I_Flux[i], I_id[i]) for i in range(0, len(I_Flux))]
    Z_Comb = [(np.round(NoNan_ZMJDs[i], decimals = 1), Z_Flux[i], Z_id[i]) for i in range(0, len(Z_Flux))]

    All_Comb_Bands = np.concatenate((U_Comb, G_Comb, R_Comb, I_Comb, Z_Comb))

    d_ = {}
    for elem in All_Comb_Bands:
        if elem[0] not in d_:
            d_[elem[0]] = []
        d_[elem[0]].append(elem[1:])
 
    
    ta = []
    MJD = []
    a = list(d_.keys())
    
    
    for i in range(len(d_)):
    
        if len(d_[a[i]]) == 5:
            ta.append(d_[a[i]])
            MJD.append(a[i])
         
            
    Lum_Vals = []
    
    for i in range(len(ta)):
            
        Eff_Wave = np.array([3000*1.041, 3652, 4825, 6410, 7980, 9097, 10000*1.041]) * (1 / (1 + 0.041))
        Flux_Val = [ta[i][0][0]*0.5, ta[i][0][0], ta[i][1][0], ta[i][2][0], ta[i][3][0], ta[i][4][0], ta[i][4][0]*0.75]
  
        tck = UnivariateSpline(Eff_Wave, Flux_Val, k = 1, s = 0)
        xnew = np.arange(3000, 10000, 1)
        Lum_Val = tck.integral(3000, 10000)
        Lum_Vals.append(Lum_Val)

#        plt.figure()
#        plt.plot(Eff_Wave, Flux_Val, '--*', label = 'SN2019cri', ms = 10, color = 'gold', lw = 1.5)
#        plt.plot(xnew, tck(xnew), color = 'green')
#        plt.xlabel('Wavelength [$\AA$]')
#        plt.ylabel('Flux [$erg\ s^{-1} cm^{-2} \AA^{-1}$]')
#        plt.title('Spectral Energy Distribution of SN2019cri')
#        plt.show()
#        plt.close()
        
    
    return Lum_Vals, MJD, d_

Lum_Vals, MJD, d = sed()


def create_psu_bol(Lum_Data):
    
    
    LD = dl(36.089)
    
    Lum = lum(np.array(Lum_Data), LD)
        
        
    return Lum


def plot_psu_bol():
    
    Psu_Bol = create_psu_bol(Lum_Vals) 
    Psu_Bol = list(Psu_Bol)
    Max_Pos_Idx = Psu_Bol.index(max(Psu_Bol))
    MJD_Corr = MJD[Max_Pos_Idx]

    plt.figure()
    plt.plot(MJD - MJD_Corr, np.log10(Psu_Bol), '--*', color = 'red')
    plt.xlabel('Days since Maximum Luminosity')
    plt.ylabel('Log(L) [$erg\ s^{-1}$]')
    plt.title('Pesudo Bolometric Light Curve of SN2019cri')
    
    Days_Since = MJD - MJD_Corr
    
    return Days_Since, MJD_Corr, Psu_Bol

Days_Since, MJD_Corr, Psu_Bol = plot_psu_bol()


plt.figure()
plt.plot(NoNan_RMJDs, np.array(NoNan_RMags) * 1.2, '--o', label = 'R-Band')
plt.plot(NoNan_GMJDs, np.array(NoNan_GMags) * 1.35, '--o', label = 'G-Band')
plt.plot(NoNan_UMJDs, np.array(NoNan_UMags) * 1.5, '--o', label = 'U-Band')
plt.plot(NoNan_IMJDs, np.array(NoNan_IMags) * 1.9, '--o', label = 'I-Band')
plt.plot(NoNan_ZMJDs, np.array(NoNan_ZMags) * 2.15, '--o', label = 'Z-Band')
plt.gca().invert_yaxis()
plt.xlabel('MJD')
plt.ylabel('Magnitudes + Offset')
plt.title('Plot of SN2019cri Light Curves for Different Bands')
plt.legend(loc = 'best')




def get_rise_time(MJD, Flux, min_day, max_day, Letter):
    
    Corr_MJD = MJD - MJD_Corr 
    
    tck = UnivariateSpline(Corr_MJD[min_day : max_day], Flux[min_day : max_day], k = 3, w = 1/np.array(Flux[min_day : max_day]))
    plt.figure()
    plt.plot(Corr_MJD, Flux, '--o', color = 'blue')
    plt.plot(Corr_MJD[min_day : max_day], tck(Corr_MJD[min_day : max_day]), color = 'red', label = 'Spline Fit')
    plt.xlabel('Days since Maximum Luminosity')
    plt.ylabel(f'{Letter}-Band Log(L) [$erg\ s^{-1}$]')
    plt.title(f'Spline Fit of the {Letter}-Band Light Curve of SN2019cri')
    

    Lumin_Vals = []
    Days = np.arange(Corr_MJD[min_day], Corr_MJD[max_day], 1)

    for i in range(len(Days)):
    
        Vals = tck(Days[i])
        Lumin_Vals.append(Vals)
        
    Max_Lum = max(Lumin_Vals)
    Max_Lum_Idx = Lumin_Vals.index(Max_Lum)
    linear_coeff, cov = np.polyfit(Corr_MJD[min_day : max_day], Flux[min_day : max_day], 1, w = 1/np.array(Flux[min_day : max_day]), cov = True)
    slope = linear_coeff[0]
    intercept = linear_coeff[1]
    x_new = np.linspace(-58, 0, 5000)
    y_new = np.poly1d(linear_coeff)
    plt.plot(x_new, y_new(x_new), label = 'Line of Best Fit') 
    plt.legend(loc = 'best')
    Min_Day = -intercept / slope
    
    Rise_Time = Days[Max_Lum_Idx] - Min_Day
    
    return Rise_Time
    
Rise_Time_1 = get_rise_time(NoNan_RMJDs, R_Flux, 0, 20, 'r')
Rise_Time_2 = get_rise_time(NoNan_RMJDs, R_Flux, 0, 19, 'r')
Rise_Time_3 = get_rise_time(NoNan_RMJDs, R_Flux, 1, 21, 'r')
Rise_Time_4 = get_rise_time(NoNan_RMJDs, R_Flux, 1, 18, 'r')
Rise_Time_5 = get_rise_time(NoNan_RMJDs, R_Flux, 1, 22, 'r')
Rise_Times = [Rise_Time_1, Rise_Time_2, Rise_Time_3, Rise_Time_4, Rise_Time_5]
Rise_Times_Mean = np.mean(np.array(Rise_Times))
Rise_Times_Diff = Rise_Times_Mean - Rise_Times
Rise_Times_Err = (np.sum((Rise_Times_Diff)**2))**0.5
print(f'Range of Rise Time : [{min(Rise_Times)} , {max(Rise_Times)}] +- {Rise_Times_Err}')
    


def fit_lum_peak(min_day, max_day): 
    
    tck = UnivariateSpline(Days_Since[min_day : max_day], (Psu_Bol[min_day : max_day]), k = 2, w = 1/np.array(Psu_Bol[min_day : max_day]))
    plt.figure()
    plt.plot(Days_Since, (Psu_Bol), '--o', color = 'blue')
    plt.plot(Days_Since[min_day : max_day], tck(Days_Since[min_day : max_day]), color = 'red')
    plt.xlabel('Days since Maximum Luminosity')
    plt.ylabel('Log(L) [$erg\ s^{-1}$]')
    plt.title('Spline Fit of the Pseudo Bolometric Light Curve of SN2019cri')
    plt.close()

    Lumin_Vals = []
    Days = np.arange(min_day, max_day, 1)

    for i in range(len(Days)):
    
        Vals = tck(Days[i])
        Lumin_Vals.append(Vals)
        
    Max_Lum = max(Lumin_Vals)
    Max_Lum_Idx = Lumin_Vals.index(Max_Lum)
    
    
    tck = UnivariateSpline(Days_Since[3:28], Psu_Bol[3:28], k = 2, w = 1/np.array(Psu_Bol[3:28]))
    
    plt.figure()
    plt.plot(Days_Since, Psu_Bol, '--o', color = 'blue')
    plt.plot(Days_Since[3:28], tck(Days_Since[3:28]), color = 'red')
    plt.xlabel('Days since Maximum Luminosity')
    plt.ylabel('Log(L) [$erg\ s^{-1}$]')
    plt.title('Spline Fit of the Pseudo Bolometric Light Curve of SN2019cri')
    plt.show()
    plt.close()
    


    Lumin_Vals_2 = []
    Days_2 = np.arange(min(Days_Since), max(Days_Since), 1)

    for i in range(len(Days_2)):
    
        Vals = tck(Days_2[i])
        Lumin_Vals_2.append(Vals)

    Half_Max_Lum_1 = Max_Lum/2
    Half_Max_Lum_2_Rise = min(Lumin_Vals_2[0 : round(len(Lumin_Vals_2)/2)], key=lambda x:abs(x-Half_Max_Lum_1))
    Half_Max_Lum_2_Decay = min(Lumin_Vals_2[round(len(Lumin_Vals_2)/2) : len(Lumin_Vals_2)], key=lambda x:abs(x-Half_Max_Lum_1))
    Half_Max_Lum_Rise_Idx = Lumin_Vals_2.index(Half_Max_Lum_2_Rise)
    Half_Max_Lum_Decay_Idx = Lumin_Vals_2.index(Half_Max_Lum_2_Decay)
    
    
    Half_Rise_Time = Days[Max_Lum_Idx]  - Days_2[Half_Max_Lum_Rise_Idx]
    Decay_Time = Days_2[Half_Max_Lum_Decay_Idx] - Days[Max_Lum_Idx]
    
    LCW = Half_Rise_Time + Decay_Time
    
    
    return Max_Lum, Half_Rise_Time, Decay_Time, LCW

  
    
Max_Lum_1, Half_Rise_Time_1, Decay_Time_1, LCW_1 = fit_lum_peak(0, 18)
Max_Lum_2, Half_Rise_Time_2, Decay_Time_2, LCW_2 = fit_lum_peak(1, 17)
Max_Lum_3, Half_Rise_Time_3, Decay_Time_3, LCW_3 = fit_lum_peak(2, 19)
Max_Lum_4, Half_Rise_Time_4, Decay_Time_4, LCW_4 = fit_lum_peak(0, 19)
Max_Lum_5, Half_Rise_Time_5, Decay_Time_5, LCW_5 = fit_lum_peak(1, 18)



def max_bol_lum(Max_Lum):
    
    Bol_Max_Lum = 1.62337868609 + (0.963813820713 * np.log10(Max_Lum))
    
    Bol_Max_Lum_ = ((Bol_Max_Lum/100) * 1) + Bol_Max_Lum
    
    return Bol_Max_Lum_


Bol_Max_Lum_1 = max_bol_lum(Max_Lum_1)
Bol_Max_Lum_2 = max_bol_lum(Max_Lum_2)
Bol_Max_Lum_3 = max_bol_lum(Max_Lum_3)
Bol_Max_Lum_4 = max_bol_lum(Max_Lum_4)
Bol_Max_Lum_5 = max_bol_lum(Max_Lum_5)
Bol_Max = [Bol_Max_Lum_1, Bol_Max_Lum_2, Bol_Max_Lum_3, Bol_Max_Lum_4, Bol_Max_Lum_5]
Bol_Max_Mean = np.mean(np.array(Bol_Max))
Bol_Max_Diff = Bol_Max_Mean - Bol_Max
Bol_Max_Err = (np.sum((Bol_Max_Diff)**2))**0.5
print(f'Bolometric Peak Luminosity Range : [{min(Bol_Max)} , {max(Bol_Max)}] +- {Bol_Max_Err}')


def calc_ni_mass(Rise_Time, Bol_Max_Lum):
    
    Ni_Mass = 10**Bol_Max_Lum * (1E43**-1) * (6.45 * math.exp(-Rise_Time / 8.8) + 1.45 * math.exp(-Rise_Time / 111.3))**-1

    return Ni_Mass


Ni_Mass_1 = calc_ni_mass(Rise_Time_1, Bol_Max_Lum_1)
Ni_Mass_2 = calc_ni_mass(Rise_Time_2, Bol_Max_Lum_2)
Ni_Mass_3 = calc_ni_mass(Rise_Time_3, Bol_Max_Lum_3)
Ni_Mass_4 = calc_ni_mass(Rise_Time_4, Bol_Max_Lum_4)
Ni_Mass_5 = calc_ni_mass(Rise_Time_5, Bol_Max_Lum_5)
Ni_Mass = [Ni_Mass_1, Ni_Mass_2, Ni_Mass_3, Ni_Mass_4, Ni_Mass_5]
Ni_Mass_Mean = np.mean(np.array(Ni_Mass))
Ni_Mass_Diff = Ni_Mass_Mean - Ni_Mass
Ni_Mass_Err = (np.sum((Ni_Mass_Diff)**2))**0.5
print(f'Ni_56 Mass Range : [{min(Ni_Mass)}, {max(Ni_Mass)}] +- {Ni_Mass_Err}')


def calc_ejecta_mass(Rise_Time):
    
    Rise_Time = Rise_Time * 60 * 60 * 24
    Ejecta_Mass = (1/2 * ((13.7 * 3E10) / 0.07) * Rise_Time**2 * 3.5E8) / 2E33
    
    return Ejecta_Mass


Ejecta_Mass_1 = calc_ejecta_mass(Rise_Time_1)
Ejecta_Mass_2 = calc_ejecta_mass(Rise_Time_2)
Ejecta_Mass_3 = calc_ejecta_mass(Rise_Time_3)
Ejecta_Mass_4 = calc_ejecta_mass(Rise_Time_4)
Ejecta_Mass_5 = calc_ejecta_mass(Rise_Time_5)
Ejecta_Mass = [Ejecta_Mass_1, Ejecta_Mass_2, Ejecta_Mass_3, Ejecta_Mass_4, Ejecta_Mass_5]
Ejecta_Mass_Mean = np.mean(np.array(Ejecta_Mass))
Ejecta_Mass_Diff = Ejecta_Mass_Mean - Ejecta_Mass
Ejecta_Mass_Err = (np.sum((Ejecta_Mass_Diff)**2))**0.5 
print(f'Ejecta_Mass Rnage : [{min(Ejecta_Mass)}, {max(Ejecta_Mass)}] +- {Ejecta_Mass_Err}')


All_SN_Data = pd.read_csv ('ALL_SN_DATA.csv')
SN_Name = list(All_SN_Data['SN'])
Type = list(All_SN_Data['Type'])
Ejecta_Masses = np.array(All_SN_Data['Ejecta_Masses'])           
Ni_Masses = np.array(All_SN_Data['Ni_Mass'])     
Log_Lp = np.array(All_SN_Data['log(L_p)']) 
t_p = np.array(All_SN_Data['t_p'])

Ni_Masses_Err_ = np.array(All_SN_Data['Err_Ni_Mass']) 
Log_Lp_Err_ = np.array(All_SN_Data['Err_log(L_p)'])
t_p_Err_ = np.array(All_SN_Data['Err_t_p'])
Ejecta_Masses_Err_ = (Ejecta_Masses / 100) * 30
    




def histogram_plots():
    
    figure, axes = plt.subplots(nrows = 2, ncols = 2)
    figure.suptitle('Histogram Plots of Supernova Parameters for Different Types of Supernovae')

    axes[0, 0].hist(Log_Lp[0:47], bins = 10, facecolor = 'orange', edgecolor = 'k', alpha = 0.75)
    axes[0, 0].hist(Log_Lp[-1], bins = 5, facecolor = 'pink', edgecolor = 'k', alpha = 0.75)
    axes[0, 0].set_xlabel('Log($L_{p}$) [$erg\ s^{-1}$]')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Log of the Peak Luminosity')
    axes[0, 0].axvline(Log_Lp[-1], color = 'k', linestyle = 'dashed', linewidth = 2)
    axes[0, 0].text(42.51, 12, f'SN2019cri Log($L_p$) = {Log_Lp[-1]}') 
     
    axes[0, 1].hist(t_p[0:47], bins = 'auto', facecolor = 'green', edgecolor = 'k', alpha = 0.75)
    axes[0, 1].hist(t_p[-1], bins = 'auto', facecolor = 'pink', edgecolor = 'k', alpha = 0.75)
    axes[0, 1].set_xlabel('$t_{p}$ [Days]')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Rise Time')
    axes[0, 1].axvline(t_p[-1], color = 'k', linestyle = 'dashed', linewidth = 2)
    axes[0, 1].text(29, 13, f'SN2019cri $t_p$ = {t_p[-1]} ')

    axes[1, 0].hist(Ni_Masses[0:47], bins = 'auto', facecolor = 'purple', edgecolor = 'k', alpha = 0.75)
    axes[1, 0].hist(Ni_Masses[-1], bins = 5, facecolor = 'pink', edgecolor = 'k', alpha = 0.75)
    axes[1, 0].set_xlabel('$Ni^{56}$ Mass [$M_{\odot}$]')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('$Ni^{56}$ Mass')
    axes[1, 0].axvline(Ni_Masses[-1], color = 'k', linestyle = 'dashed', linewidth = 2)
    axes[1, 0].text(0.3, 17, 'SN2019cri $Ni^{56}$ Mass = ' f'{Ni_Masses[-1]}') 
    
    axes[1, 1].hist(Ejecta_Masses[0:47], bins = 15, facecolor = 'brown', edgecolor = 'k', alpha = 0.75)
    axes[1, 1].hist(Ejecta_Masses[-1], bins = 'auto', facecolor = 'pink', edgecolor = 'k', alpha = 0.75)
    axes[1, 1].set_xlabel('Ejecta Mass [$M_{\odot}$]')
    axes[1, 1].set_ylabel('Frequency')
    axes[1, 1].set_title('Ejecta Mass')
    axes[1, 1].axvline(Ejecta_Masses[-1], color = 'k', linestyle = 'dashed', linewidth = 2)
    axes[1, 1].text(4.3, 9.2, f'SN2019cri Ejecta Mass = {Ejecta_Masses[-1]} ') 

    figure.tight_layout()
    figure.subplots_adjust(top = 0.88)
    
    return 

histogram_plots()



def scatter_plots():
    
    def func(x, a, b):
        return (a * x) ** b
    
    Ic_Ni = []
    Ic_Ejecta = []
    Ic_tp = []
    Ic_Lp = []

    Ic_BL_Ni = []
    Ic_BL_Ejecta = []
    Ic_BL_tp = []
    Ic_BL_Lp = []  

    Ib_Ni = []
    Ib_Ejecta = []
    Ib_tp = []
    Ib_Lp = []   
 
    IIb_Ni = []
    IIb_Ejecta = []
    IIb_tp = []
    IIb_Lp = [] 

    GRB_SN_Ni = []
    GRB_SN_Ejecta = []
    GRB_SN_tp = []
    GRB_SN_Lp = [] 

    plt.figure()
    plt.ylabel('Ejecta Mass  [$M_{\odot}$]')
    plt.xlabel('$Ni^{56}$ Mass [$M_{\odot}$]')
    plt.title('Plot of Ni Mass versus Ejecta Mass for Different Types of Supernovae')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(10**-2,10**0.5)
    plt.ylim(10**-0.3,10**1.2)
    plt.minorticks_on()
    
    
    for i in range(len(Ejecta_Masses)-1):
    
        if Type[i] == 'Ic':
            
           Ic_Ni.append(Ni_Masses[i])
           Ic_Ejecta.append(Ejecta_Masses[i])

           Ic = plt.scatter(Ni_Masses[i], Ejecta_Masses[i], s = 80, marker = 'd', color = 'red')
           plt.errorbar(Ni_Masses[i], Ejecta_Masses[i], Ejecta_Masses_Err_[i], Ni_Masses_Err_[i], color = 'red', elinewidth = 1)
    
        if Type[i] == 'Ic-BL':
            
            Ic_BL_Ni.append(Ni_Masses[i])
            Ic_BL_Ejecta.append(Ejecta_Masses[i])
            
            Ic_BL = plt.scatter(Ni_Masses[i], Ejecta_Masses[i], s = 80, marker = '<', color = 'blue')
            plt.errorbar(Ni_Masses[i], Ejecta_Masses[i], Ejecta_Masses_Err_[i], Ni_Masses_Err_[i], color = 'blue', elinewidth = 1)

        if Type[i] == 'Ib':
            
            Ib_Ni.append(Ni_Masses[i])
            Ib_Ejecta.append(Ejecta_Masses[i])

            
            Ib = plt.scatter(Ni_Masses[i], Ejecta_Masses[i], s = 80, marker = 'o', color = 'green')
            plt.errorbar(Ni_Masses[i], Ejecta_Masses[i], Ejecta_Masses_Err_[i], Ni_Masses_Err_[i], color = 'green', elinewidth = 1) 
        
        if Type[i] == 'IIb':
            
            IIb_Ni.append(Ni_Masses[i])
            IIb_Ejecta.append(Ejecta_Masses[i])

            
            IIb = plt.scatter(Ni_Masses[i], Ejecta_Masses[i], s = 80, marker = '>', color = 'orange')
            plt.errorbar(Ni_Masses[i], Ejecta_Masses[i], Ejecta_Masses_Err_[i], Ni_Masses_Err_[i], color = 'orange', elinewidth = 1)

        if Type[i] == 'GRB-SN':
            
            GRB_SN_Ni.append(Ni_Masses[i])
            GRB_SN_Ejecta.append(Ejecta_Masses[i])
            
            GRB_SN = plt.scatter(Ni_Masses[i], Ejecta_Masses[i], s = 100, marker = '*', color = 'palevioletred')
            plt.errorbar(Ni_Masses[i], Ejecta_Masses[i], Ejecta_Masses_Err_[i], Ni_Masses_Err_[i], color = 'palevioletred', elinewidth = 1)
        

    cri = plt.scatter(Ni_Masses[-1], Ejecta_Masses[-1], s = 80, marker = '^', color = 'lawngreen')
    plt.errorbar(Ni_Masses[-1], Ejecta_Masses[-1], Ejecta_Masses_Err_[-1], Ni_Masses_Err_[-1], color = 'lawngreen', elinewidth = 1)
    
    Ni_Masses_1 = Ni_Masses[np.logical_not(np.isnan(Ni_Masses))]
    Ejecta_Masses_1 = Ejecta_Masses[np.logical_not(np.isnan(Ni_Masses))]
    Ejecta_Masses_2 = Ejecta_Masses_1[np.logical_not(np.isnan(Ejecta_Masses_1))]
    Ni_Masses_2 = Ni_Masses_1[np.logical_not(np.isnan(Ejecta_Masses_1))]
    
    plt.legend([Ic, Ic_BL, Ib, IIb, cri, GRB_SN], ['Ic', 'Ic_BL', 'Ib', 'IIb', '2019cri', 'GRB-SN'], loc = 'upper left')
    plt.show()
    plt.close()
    
    
    R_coeff_1 = np.nansum((Ni_Masses - np.nanmean(Ni_Masses)/np.nanstd(Ni_Masses)) * (Ejecta_Masses - np.nanmean(Ejecta_Masses)/np.nanstd(Ejecta_Masses))) / (len(Ejecta_Masses) - 1)

    
    # Second
    plt.figure()
    plt.xlabel('$t_{p}$ [Days]')
    plt.ylabel('Log($L_{p}$) [$erg\ s^{-1}$]')
    plt.xscale('linear') 
    plt.xlim(10**0.8,10**1.88)
    plt.title('Plot of the Log(Peak Luminosity) versus Rise Time for Different Types of Supernovae')
    plt.minorticks_on()
    

    
    
    for i in range(len(Ejecta_Masses)-1):
    
        if Type[i] == 'Ic':

           Ic_tp.append(t_p[i])
           Ic_Lp.append(Log_Lp[i])
        
           Ic = plt.scatter(t_p[i], Log_Lp[i], s = 80, marker = 'd', color = 'red')
           plt.errorbar(t_p[i], Log_Lp[i], Log_Lp_Err_[i], t_p_Err_[i], color = 'red', elinewidth = 1)
          
    
        if Type[i] == 'Ic-BL':
            
            Ic_BL_tp.append(t_p[i])
            Ic_BL_Lp.append(Log_Lp[i])
        
            Ic_BL = plt.scatter(t_p[i], Log_Lp[i],  s = 80, marker = '<', color = 'blue')
            plt.errorbar(t_p[i], Log_Lp[i], Log_Lp_Err_[i], t_p_Err_[i],  color = 'blue', elinewidth = 1)

        if Type[i] == 'Ib':
            
            Ib_tp.append(t_p[i])
            Ib_Lp.append(Log_Lp[i])
        
            Ib = plt.scatter(t_p[i], Log_Lp[i],  s = 80, marker = 'o', color = 'green')
            plt.errorbar(t_p[i], Log_Lp[i], Log_Lp_Err_[i], t_p_Err_[i],  color = 'green', elinewidth = 1) 
        
        if Type[i] == 'IIb':
            
            IIb_tp.append(t_p[i])
            IIb_Lp.append(Log_Lp[i])
        
            IIb = plt.scatter(t_p[i], Log_Lp[i],  s = 80, marker = '>', color = 'orange')
            plt.errorbar(t_p[i], Log_Lp[i], Log_Lp_Err_[i], t_p_Err_[i],  color = 'orange', elinewidth = 1)
            

        if Type[i] == 'GRB-SN':
            
            GRB_SN_tp.append(t_p[i])
            GRB_SN_Lp.append(Log_Lp[i])
        
            GRB_SN = plt.scatter(t_p[i], Log_Lp[i],  s = 100, marker = '*', color = 'palevioletred')
            plt.errorbar(t_p[i], Log_Lp[i], Log_Lp_Err_[i], t_p_Err_[i], color = 'palevioletred', elinewidth = 1)
            
        

    cri = plt.scatter(t_p[-1], Log_Lp[-1],  s = 80, marker = '^', color = 'lawngreen')
    plt.errorbar(t_p[-1], Log_Lp[-1], Log_Lp_Err_[-1], t_p_Err_[-1], color = 'lawngreen', elinewidth = 1)
    
    plt.legend([Ic, Ic_BL, Ib, IIb, cri, GRB_SN], ['Ic', 'Ic_BL', 'Ib', 'IIb', '2019cri', 'GRB_SN'], loc = 'upper right')
    plt.show()
    plt.close()
       
    
    # Third
    plt.figure()
    plt.xlabel('Log($L_{p}$) [$erg\ s^{-1}$]')
    plt.ylabel('Ejecta Mass [$M_{\odot}$]')
    plt.yscale('log')
    plt.ylim(10**-0.3, 10**1.3)
    plt.xlim(41.7,43.8)
    plt.title('Plot of the Log(Peak Luminosity) versus Ejecta Mass for Different Types of Supernovae')
    plt.minorticks_on()
    
    
    for i in range(len(Ejecta_Masses)-1):
    
        if Type[i] == 'Ic':
        
           Ic = plt.scatter(Log_Lp[i], Ejecta_Masses[i], s = 80, marker = 'd', color = 'red')
           plt.errorbar(Log_Lp[i], Ejecta_Masses[i], Ejecta_Masses_Err_[i], Log_Lp_Err_[i], color = 'red', elinewidth = 1)
    
        if Type[i] == 'Ic-BL':
        
            Ic_BL = plt.scatter(Log_Lp[i], Ejecta_Masses[i], s = 80, marker = '<', color = 'blue')
            plt.errorbar(Log_Lp[i], Ejecta_Masses[i], Ejecta_Masses_Err_[i], Log_Lp_Err_[i], color = 'blue', elinewidth = 1)

        if Type[i] == 'Ib':
        
            Ib = plt.scatter(Log_Lp[i], Ejecta_Masses[i], s = 80, marker = 'o', color = 'green')
            plt.errorbar(Log_Lp[i], Ejecta_Masses[i], Ejecta_Masses_Err_[i], Log_Lp_Err_[i], color = 'green', elinewidth = 1) 
        
        if Type[i] == 'IIb':
        
            IIb = plt.scatter(Log_Lp[i], Ejecta_Masses[i], s = 80, marker = '>', color = 'orange')
            plt.errorbar(Log_Lp[i], Ejecta_Masses[i], Ejecta_Masses_Err_[i], Log_Lp_Err_[i], color = 'orange', elinewidth = 1)

        if Type[i] == 'GRB-SN':
        
            GRB_SN = plt.scatter(Log_Lp[i], Ejecta_Masses[i], s = 100, marker = '*', color = 'palevioletred')
            plt.errorbar(Log_Lp[i], Ejecta_Masses[i], Ejecta_Masses_Err_[i], Log_Lp_Err_[i], color = 'palevioletred', elinewidth = 1)
        

    cri = plt.scatter(Log_Lp[-1], Ejecta_Masses[-1], s = 80, marker = '^', color = 'lawngreen')
    plt.errorbar(Log_Lp[-1], Ejecta_Masses[-1], Ejecta_Masses_Err_[-1], Log_Lp_Err_[-1], color = 'lawngreen', elinewidth = 1)
    
    plt.legend([Ic, Ic_BL, Ib, IIb, cri, GRB_SN], ['Ic', 'Ic_BL', 'Ib', 'IIb', '2019cri', 'GRB_SN'], loc = 'upper left')
    plt.show()
    plt.close()
    
    
    R_coeff_2 = np.nansum(((Log_Lp - np.mean(Log_Lp))/np.nanstd(Log_Lp)) * ((Ejecta_Masses - np.nanmean(Ejecta_Masses))/np.nanstd(Ejecta_Masses))) / (len(Ejecta_Masses) - 1)


    
    # Fourth
    plt.figure()
    plt.xlabel('Log($L_{p}$) [$erg\ s^{-1}$]')
    plt.ylabel('$Ni^{56}$ Mass [$M_{\odot}$]')
    plt.yscale('log')
    plt.ylim(10**-2, 10**0.5)
    plt.title('Plot of the Log(Peak Luminosity) versus $Ni^{56}$ Mass for Different Types of Supernovae')
    plt.minorticks_on()
    
    
    for i in range(len(Ejecta_Masses)-1):
    
        if Type[i] == 'Ic':
        
           Ic = plt.scatter(Log_Lp[i], Ni_Masses[i], s = 80, marker = 'd', color = 'red')
           plt.errorbar(Log_Lp[i], Ni_Masses[i], Ni_Masses_Err_[i], Log_Lp_Err_[i], color = 'red', elinewidth = 1)
    
        if Type[i] == 'Ic-BL':
        
            Ic_BL = plt.scatter(Log_Lp[i], Ni_Masses[i], s = 80, marker = '<', color = 'blue')
            plt.errorbar(Log_Lp[i], Ni_Masses[i], Ni_Masses_Err_[i], Log_Lp_Err_[i], color = 'blue', elinewidth = 1)

        if Type[i] == 'Ib':
        
            Ib = plt.scatter(Log_Lp[i], Ni_Masses[i], s = 80, marker = 'o', color = 'green')
            plt.errorbar(Log_Lp[i], Ni_Masses[i], Ni_Masses_Err_[i], Log_Lp_Err_[i], color = 'green', elinewidth = 1) 
        
        if Type[i] == 'IIb':
        
            IIb = plt.scatter(Log_Lp[i], Ni_Masses[i], s = 80, marker = '>', color = 'orange')
            plt.errorbar(Log_Lp[i], Ni_Masses[i], Ni_Masses_Err_[i], Log_Lp_Err_[i], color = 'orange', elinewidth = 1)

        if Type[i] == 'GRB-SN':
        
            GRB_SN = plt.scatter(Log_Lp[i], Ni_Masses[i], s = 100, marker = '*', color = 'palevioletred')
            plt.errorbar(Log_Lp[i], Ni_Masses[i], Ni_Masses_Err_[i], Log_Lp_Err_[i], color = 'palevioletred', elinewidth = 1)
        

    cri = plt.scatter(Log_Lp[-1], Ni_Masses[-1], s = 80, marker = '^', color = 'lawngreen')
    plt.errorbar(Log_Lp[-1], Ni_Masses[-1], Ni_Masses_Err_[-1], Log_Lp_Err_[-1], color = 'lawngreen', elinewidth = 1)
    
    Ni_Masses_1 = Ni_Masses[np.logical_not(np.isnan(Ni_Masses))]
    Log_Lp_1 = Log_Lp[np.logical_not(np.isnan(Ni_Masses))]
    Ni_Masses_Err_1 = Ni_Masses_Err_[np.logical_not(np.isnan(Ni_Masses))]
    Log_Lp_Err_1 = Log_Lp_Err_[np.logical_not(np.isnan(Ni_Masses))]
    
 
    popt, pcov = curve_fit(func, Log_Lp_1, Ni_Masses_1, sigma = 1/Log_Lp_Err_1)
    x_new_1 = np.linspace(41.8,43.5, 10000)
    y_new_1 = (popt[0] * x_new_1) ** popt[1]
    plt.plot(x_new_1, y_new_1, 'r-')
    
    plt.legend([Ic, Ic_BL, Ib, IIb, cri, GRB_SN], ['Ic', 'Ic_BL', 'Ib', 'IIb', '2019cri', 'GRB_SN'], loc = 'upper left')
    plt.show()
    plt.close()
    
    
    Cal_Ni_Mass = (popt[0] * Log_Lp_1) ** popt[1]
    print(f'Fit Equation: Ni_Mass = {popt[0]} * Log(L_p) ** {popt[1]}')
    Ni_Diff = ((Ni_Masses_1 - Cal_Ni_Mass)/Cal_Ni_Mass)*100
    Mean_Diff = mean(Ni_Diff)
    Median_Diff = statistics.median(Ni_Diff)
    Stdev_Diff = statistics.stdev(Ni_Diff, Mean_Diff)

    plt.figure()
    plt.hist(Ni_Diff[0:45], bins = 13, facecolor = 'g', edgecolor = 'k', alpha = 0.75)
    plt.hist(Ni_Diff[-1], bins = 1, facecolor = 'purple', edgecolor = 'k', alpha = 0.75)
    plt.xlabel('Percentage [%]')
    plt.ylabel('Frequency')
    plt.title('Histogram of Percentage Differnce in Calculated $Ni^{56}$ Mass')
    plt.text(30, 9, f'$\mu$ = {Mean_Diff} % \n \n $\sigma$ = {Stdev_Diff} % \n \n Median = {Median_Diff} % ')
    plt.axvline(Median_Diff, color = 'k', linestyle = 'dashed', linewidth = 2)
    plt.text(130, 4, 'SN2019cri ')
    plt.axvline(Ni_Diff[-1], color = 'k', linestyle = 'dashed', linewidth = 0.9)
    plt.minorticks_on()
    
    return 

scatter_plots()

