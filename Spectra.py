# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 12:49:53 2020

@author: Bill Kamtchou
"""

from astropy.cosmology import WMAP5
from astropy.coordinates import Distance
from astropy.io import fits 
from astropy import units as u
import matplotlib.pyplot as plt
import numpy as np
import os 
import pandas as pd 
from pylab import rcParams
from scipy.interpolate import UnivariateSpline
from SN_Analysis import MJD_Corr, dered_mag

All_Spectra_Files = os.listdir('2019criSpectra') 


def get_spectra(spectra_data_index, spectra_date_index):
    
    Wavelength, Flux = np.loadtxt(f'2019criSpectra/{All_Spectra_Files[spectra_data_index]}', unpack = True, usecols = (0,1))


    image_hdu = fits.open(f'2019criSpectra/{All_Spectra_Files[spectra_date_index]}')

    Observation_Date = image_hdu[0].header['DATE-OBS']
    MJD = image_hdu[0].header['MJD']
    Band_Filter = image_hdu[0].header['FILTER1']
    Exposure_Time = image_hdu[0].header['EXPTIME']
    image_hdu.close()
    
    
    Dict_Data = {
        'Wavelength'           : Wavelength,
        'Flux'                 : Flux,
		'Observation_Date'     : Observation_Date,
		'Modified_Julian_Date' : MJD,
        'Filter'               : Band_Filter,
		'Exposure_Duration'    : Exposure_Time  # Exposure duration is given in seconds 

	}

    #print(f'Observation_Date : {Observation_Date}')
    #print(f'Modified Julian Date : {MJD}')
    #print(f'Filter : {Band_Filter}')
    #print(f'Exposure Duration : {Exposure_Time} s')
    #print(f'{All_Spectra_Files[spectra_data_index]} : {All_Spectra_Files[spectra_date_index]}')
    
    #plt.figure()
    #plt.plot(Wavelength, Flux)
    #plt.xlabel('Wavelength [$\AA$]')
    #plt.ylabel('Flux')
    #plt.title(f'SN2019cri Spectra, Date: {Observation_Date}')
    
    return Dict_Data 



def get_spectra_asci(spectra_data_index, spectra_date_index):
    
    Wavelength_, Flux_ = np.loadtxt(f'2019criSpectra/{All_Spectra_Files[spectra_data_index]}', unpack = True, usecols = (0,1))


    image_hdu = fits.open(f'2019criSpectra/{All_Spectra_Files[spectra_date_index]}')

    Observation_Date = image_hdu[0].header['DATE-OBS']
    MJD = image_hdu[0].header['MJD-OBS']
    Exposure_Time = image_hdu[0].header['EXPTIME']
    image_hdu.close()
    
    
    Dict_Data_ = {
        'Wavelength'           : Wavelength_,
        'Flux'                 : Flux_,
		'Observation_Date'     : Observation_Date,
		'Modified_Julian_Date' : MJD,
		'Exposure_Duration'    : Exposure_Time

	}
    
    #print(f'{All_Spectra_Files[spectra_data_index]} : {All_Spectra_Files[spectra_date_index]}')
    
    #plt.figure()
    #plt.plot(Wavelength_, Flux_)
    #plt.xlabel('Wavelength [$\AA$]')
    #plt.ylabel('Magnitude')
    #plt.title(f'SN2019cri Spectra, Date: {Observation_Date}', pad = 20)
    
    return Dict_Data_




def plot_my_spectra():
    
    Dict_Data_420 = get_spectra(1,20)
    Dict_Data_421= get_spectra(2,21)
    Dict_Data_422 = get_spectra(11,22)
    Dict_Data_423 = get_spectra(12,23)
    Dict_Data_424 = get_spectra(13,24)
    Dict_Data_428 = get_spectra(14,25)
    Dict_Data_502 = get_spectra(15,26)
    Dict_Data_506 = get_spectra(16,27)
    Dict_Data_Gr = get_spectra_asci(4,5)
    Dict_Data_514 = get_spectra(17,28)
    Dict_Data_521 = get_spectra(18,29)
    Dict_Data_comb = get_spectra_asci(3,6)
    Dict_Data_702 = get_spectra_asci(7,8)
    Dict_Data_805 = get_spectra_asci(9,10)
    Dict_Data_229 = get_spectra(0,30)
    
    return Dict_Data_420, Dict_Data_421, Dict_Data_422, Dict_Data_423, Dict_Data_424, Dict_Data_428, Dict_Data_502, Dict_Data_506, Dict_Data_Gr, Dict_Data_514, Dict_Data_521, Dict_Data_comb, Dict_Data_702, Dict_Data_805, Dict_Data_229


Dict_Data_420, Dict_Data_421, Dict_Data_422, Dict_Data_423, Dict_Data_424, Dict_Data_428, Dict_Data_502, Dict_Data_506, Dict_Data_Gr, Dict_Data_514, Dict_Data_521, Dict_Data_comb, Dict_Data_702, Dict_Data_805, Dict_Data_229 = plot_my_spectra()







# Above is the spectra data for SN2019cri





def get_compar_sn_data():
    
    SN2007gr_Dict = {
            
            'MW_red'   : 0.055,
            'Host_red' : 0.03,
            
                    }
    
    SN1994I_Dict = {
            
            'MW_red'   : 0.03,
            'Host_red' : 0.3,
            
                    }
    
    SN2007ru_Dict = {
            
            'MW_red'   : 0.27,
            'Host_red' : 0.00000001,
            
                    }
    
    SN2011bm_Dict = {
            
            'MW_red'   : 0.032,
            'Host_red' : 0.032,
            
                    }
    
    SN2004aw_Dict = {
            
            'MW_red'   : 0.021,
            'Host_red' : 0.35,
            
                    }
    
    return SN2007gr_Dict, SN1994I_Dict, SN2007ru_Dict, SN2011bm_Dict, SN2004aw_Dict


SN2007gr_Dict, SN1994I_Dict, SN2007ru_Dict, SN2011bm_Dict, SN2004aw_Dict = get_compar_sn_data() 



def compar_sn(SN, Band, Dict):
    
    LC_Data_df = pd.read_json(f'COMPARISON_SN\{SN}.json')
    

    Mags = []           
    Time_Phot = []
    Init_Wavelength = []
    Wavelength = []
    Wavelength_ = []
    Flux = []
    Flux_ = []
    Time_Spectra = []
    Redshift = float(LC_Data_df.SN['redshift'][0]['value'])


    for i in range(len(LC_Data_df.SN['photometry'])):
    
        test = LC_Data_df.SN['photometry'][i]['band']

        if test == Band:
            Mags.append(float(LC_Data_df.SN['photometry'][i]['magnitude']))
            Time_Phot.append(float(LC_Data_df.SN['photometry'][i]['time']))
            
    
    Distance_Modulus = float(LC_Data_df.SN['maxappmag'][0]['value']) - float(LC_Data_df.SN['maxabsmag'][0]['value']) # Distance modulus is given in pc.  
    
    for j in range(len(LC_Data_df.SN['spectra'][:])):
        
        Init_Wavelength.append(LC_Data_df.SN['spectra'][j]['data'])
        
        
    for k in range(len(Init_Wavelength)):
        
        Time_Spectra.append(float(LC_Data_df.SN['spectra'][k]['time']))
        
        for l in range(len(Init_Wavelength[k])):
            
            Wavelength.append(float(Init_Wavelength[k][l][0]))
            Flux.append(float(Init_Wavelength[k][l][1]))
            
    for m in range(len(Init_Wavelength[:])):
        
        Wavelength_.append(Wavelength[0:len(Init_Wavelength[m])])
        Flux_.append(Flux[0:len(Init_Wavelength[m])])
    
    
    Dict_Data = {
            
		'Mag' : Mags,
		'MJD_Phot' : Time_Phot,
        'Wavelength': Wavelength_,
        'Flux' : Flux_,
        'Time_Spectra' : Time_Spectra,
        'Redshift' : Redshift,
        'Distance_Modulus' : Distance_Modulus,
        'MW_Red' : Dict['MW_red'],
        'Host_Red' : Dict['Host_red']

	            }
    
    #plt.figure()
    #plt.plot(Time_Phot, Mags)
    #plt.gca().invert_yaxis()
    #plt.xlabel('MJD')
    #plt.ylabel('Magnitude')
    #plt.title(f'Light Curve of {SN} in the {Band}-Band')
    
    return Dict_Data



def plot_compar_sn():
    
    SN2007gr_Dict_Data = compar_sn('SN2007gr','R', SN2007gr_Dict)
    SN1994I_Dict_Data = compar_sn('SN1994I','R', SN1994I_Dict)
    SN2007ru_Dict_Data = compar_sn('SN2007ru','R', SN2007ru_Dict)
    SN2011bm_Dict_Data = compar_sn('SN2011bm','R', SN2011bm_Dict)
    SN2004aw_Dict_Data = compar_sn('SN2004aw','R', SN2004aw_Dict)
    
    return SN2007gr_Dict_Data, SN1994I_Dict_Data, SN2007ru_Dict_Data, SN2011bm_Dict_Data, SN2004aw_Dict_Data


SN2007gr_Dict_Data, SN1994I_Dict_Data, SN2007ru_Dict_Data, SN2011bm_Dict_Data, SN2004aw_Dict_Data = plot_compar_sn()



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


def redshift_correction(mw_dered_flux, Dict):
    
    redshift_corr = mw_dered_flux * (1 / (1 + Dict['Redshift']))
    
    return redshift_corr




def my_sn_wavelength_correction(Spectra_Dict_Data):
    
    Flux_Corr = []
    Wave_Corr = []
    
    for i in range(len(Spectra_Dict_Data['Wavelength'])):
            
        mw_dered_flux = dered(Spectra_Dict_Data['Wavelength'][i], Spectra_Dict_Data['Flux'][i], 3.1, 0.0216)
            
        redshift_corr_wavelength = Spectra_Dict_Data['Wavelength'][i] * (1 / (1 + 0.041))
            
        #host_dered_flux = dered(redshift_corr, SN_Dict_Data['Flux'][i][j], 3.1, 0.0216)

            
        Flux_Corr.append(mw_dered_flux)
        Wave_Corr.append(redshift_corr_wavelength)

        
    return Flux_Corr, Wave_Corr 


Flux_Corr_420, Wave_Corr_420  = my_sn_wavelength_correction(Dict_Data_420) 
Flux_Corr_421, Wave_Corr_421  = my_sn_wavelength_correction(Dict_Data_421) 
Flux_Corr_422, Wave_Corr_422  = my_sn_wavelength_correction(Dict_Data_422) 
Flux_Corr_423, Wave_Corr_423  = my_sn_wavelength_correction(Dict_Data_423) 
Flux_Corr_424, Wave_Corr_424  = my_sn_wavelength_correction(Dict_Data_424) 
Flux_Corr_428, Wave_Corr_428  = my_sn_wavelength_correction(Dict_Data_428) 
Flux_Corr_502, Wave_Corr_502  = my_sn_wavelength_correction(Dict_Data_502) 
Flux_Corr_506, Wave_Corr_506  = my_sn_wavelength_correction(Dict_Data_506) 
Flux_Corr_Gr, Wave_Corr_Gr    = my_sn_wavelength_correction(Dict_Data_Gr) 
Flux_Corr_514, Wave_Corr_514  = my_sn_wavelength_correction(Dict_Data_514) 
Flux_Corr_521, Wave_Corr_521  = my_sn_wavelength_correction(Dict_Data_521) 
Flux_Corr_comb, Wave_Corr_comb = my_sn_wavelength_correction(Dict_Data_comb) 
Flux_Corr_702, Wave_Corr_702  = my_sn_wavelength_correction(Dict_Data_702) 
Flux_Corr_805, Wave_Corr_805  = my_sn_wavelength_correction(Dict_Data_805) 
Flux_Corr_229, Wave_Corr_229  = my_sn_wavelength_correction(Dict_Data_229) 


Oxygen = [9263, 8446, 7774, 6158]     
NaI = [5896, 5890]
FeII = [5169, 5018, 4924, 4549, 4515, 4352, 4303]
MgI = [8807, 5528, 5184, 5173, 5167, 4703, 4571, 3838]
SiII = [6371, 6347, 6355, 5958, 5670, 5056, 5041, 4131, 4128, 3856]       
CaII = [8662, 8542, 8498, 3934, 3969]


def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

def get_spec_mjd():
    
    Dicts = [Dict_Data_420, Dict_Data_421, Dict_Data_422, Dict_Data_423, Dict_Data_424, Dict_Data_428, Dict_Data_502, Dict_Data_506, Dict_Data_Gr, Dict_Data_514, Dict_Data_521, Dict_Data_comb]
    
    True_MJD = []
    
    for i in range(len(Dicts)):
    
        MJD = Dicts[i]['Modified_Julian_Date'] * (1 / (1 + 0.041))
    
        True_MJD_ = MJD - MJD_Corr
    
        True_MJD.append(True_MJD_)
    
    return True_MJD

Spec_MJD = get_spec_mjd()


def spec_evol():
    
    Dicts = [Dict_Data_423, Dict_Data_424, Dict_Data_428, Dict_Data_502,  Dict_Data_506,\
             Dict_Data_Gr,  Dict_Data_514, Dict_Data_521, Dict_Data_comb, Dict_Data_702,\
             Dict_Data_805]
    
    True_MJD = []
    
    for i in range(len(Dicts)):
    
        MJD = Dicts[i]['Modified_Julian_Date'] * (1 / (1 + 0.041))
    
        True_MJD_ = MJD - MJD_Corr
    
        True_MJD.append(True_MJD_)
        
    
    Flux_Vals = [Flux_Corr_423, Flux_Corr_424, Flux_Corr_428,Flux_Corr_502, Flux_Corr_506,\
                 Flux_Corr_Gr, Flux_Corr_514, Flux_Corr_521, Flux_Corr_comb,Flux_Corr_702, Flux_Corr_805]
    
    Wave_Vals = [Wave_Corr_423, Wave_Corr_424, Wave_Corr_428,Wave_Corr_502, Wave_Corr_506,\
                 Wave_Corr_Gr,  Wave_Corr_514, Wave_Corr_521, Wave_Corr_comb, Wave_Corr_702, Wave_Corr_805] 
    
    
    plt.figure()
    fig, axs = plt.subplots(11, 1, figsize = (10,25),facecolor = 'w', edgecolor ='k')
    fig.subplots_adjust(hspace = .0)
#    fig.suptitle('Vertically stacked subplots')
    axs = axs.ravel()
    
    for i in range(len(Flux_Vals )):
        

    
        axs[i].plot(Wave_Vals[i], np.array(Flux_Vals[i]), color = 'black')
        axs[i].set_xlabel('Wavelength [$\AA$]')
        axs[i].set_ylabel('Flux')
        axs[i].set_xlim(3500, 9000)
        
        if True_MJD[i] > 0:
            axs[i].legend([f'{np.round(True_MJD[i], decimals = 1)} Days After Maximum Luminosity'])
            
        else:
            axs[i].legend([f'{np.round(-True_MJD[i], decimals = 1)} Days Before Maximum Luminosity'])
        
    
    return

spec_evol()



def get_line_vel(Wave, Flux, Lines, Dict, Vel, Element):
    
    c = 3E5
    wave_obs = []
    
    for i in range(len(Lines)):
    
        wave_obs_ta = Lines[i] * np.sqrt((1 - Vel/c) / (1 + Vel/c))
        wave_obs.append(wave_obs_ta)

    #Obs_Date = Dict['Observation_Date']
    MJD = Dict['Modified_Julian_Date'] * (1 / (1 + 0.041))
    #cmaps = get_cmap(len(Lines))
        
    plt.figure()
    plt.plot(Wave, Flux, color = 'black', lw = 1.5)
    plt.xlabel('Wavelength [$\AA$]')
    plt.ylabel('Flux')
    plt.xlim(3000, 9000)
    
    
    
    if MJD - MJD_Corr < 0:
        
        plt.title(f'SN2019cri Spectra Analysis of {Element} taken {np.round((MJD - MJD_Corr)*-1, decimals = 1)} Days Before Maximum Luminosity')
        
    else:
        
        plt.title(f'SN2019cri Spectra Analysis of {Element} taken {np.round((MJD - MJD_Corr), decimals = 1)} Days After Maximum Luminosity')
        
    for i in range(len(wave_obs)):
        
        plt.axvline(x = wave_obs[i], linestyle = '--', lw = 1)
        
    plt.show()
#    plt.close()
    
    return Vel
        



Vel_420_FeII = get_line_vel(Wave_Corr_420, Flux_Corr_420, FeII, Dict_Data_420, 8000, 'FeII') # 5169 \AA
Vel_421_FeII = get_line_vel(Wave_Corr_421, Flux_Corr_421, FeII, Dict_Data_421, 8000, 'FeII')
Vel_422_FeII = get_line_vel(Wave_Corr_422, Flux_Corr_422, FeII, Dict_Data_422, 8000, 'FeII')
Vel_423_FeII = get_line_vel(Wave_Corr_423, Flux_Corr_423, FeII, Dict_Data_423, 7500, 'FeII')
Vel_424_FeII = get_line_vel(Wave_Corr_424, Flux_Corr_424, FeII, Dict_Data_424, 7000, 'FeII')
Vel_428_FeII = get_line_vel(Wave_Corr_428, Flux_Corr_428, FeII, Dict_Data_428, 7200, 'FeII')
Vel_502_FeII = get_line_vel(Wave_Corr_502, Flux_Corr_502, FeII, Dict_Data_502, 7200, 'FeII')
Vel_506_FeII = get_line_vel(Wave_Corr_506, Flux_Corr_506, FeII, Dict_Data_506, 6800, 'FeII')
Vel_Gr_FeII = get_line_vel(Wave_Corr_Gr, Flux_Corr_Gr, FeII, Dict_Data_Gr, 6500, 'FeII')  
Vel_514_FeII = get_line_vel(Wave_Corr_514, Flux_Corr_514, FeII, Dict_Data_514, 6500, 'FeII')       
Vel_521_FeII = get_line_vel(Wave_Corr_521, Flux_Corr_521, FeII, Dict_Data_521, 5500, 'FeII')  
Vel_comb_FeII = get_line_vel(Wave_Corr_comb, Flux_Corr_comb, FeII, Dict_Data_comb, 5000, 'FeII')  
Vel_FeII = [Vel_420_FeII, Vel_421_FeII, Vel_422_FeII, Vel_423_FeII, Vel_424_FeII, Vel_428_FeII,\
          Vel_502_FeII, Vel_506_FeII, Vel_Gr_FeII, Vel_514_FeII, Vel_521_FeII, Vel_comb_FeII]
Vel_FeII_Err = [500, 500, 700, 700, 800, 500, 800, 800, 700, 500, 700, 700]  
  
Vel_420_O = get_line_vel(Wave_Corr_420, Flux_Corr_420, Oxygen, Dict_Data_420, 8700, 'Oxygen')  # 7774 \AA
Vel_421_O = get_line_vel(Wave_Corr_421, Flux_Corr_421, Oxygen, Dict_Data_421, 8700, 'Oxygen')
Vel_422_O = get_line_vel(Wave_Corr_422, Flux_Corr_422, Oxygen, Dict_Data_422, 10700, 'Oxygen')
Vel_423_O = get_line_vel(Wave_Corr_423, Flux_Corr_423, Oxygen, Dict_Data_423, 9500, 'Oxygen')
Vel_424_O = get_line_vel(Wave_Corr_424, Flux_Corr_424, Oxygen, Dict_Data_424, 8000, 'Oxygen')
Vel_428_O = get_line_vel(Wave_Corr_428, Flux_Corr_428, Oxygen, Dict_Data_428, 8000, 'Oxygen')
Vel_502_O = get_line_vel(Wave_Corr_502, Flux_Corr_502, Oxygen, Dict_Data_502, 9500, 'Oxygen')
Vel_506_O = get_line_vel(Wave_Corr_506, Flux_Corr_506, Oxygen, Dict_Data_506, 9200, 'Oxygen')
Vel_Gr_O = get_line_vel(Wave_Corr_Gr, Flux_Corr_Gr, Oxygen, Dict_Data_Gr, 8500, 'Oxygen')  
Vel_514_O = get_line_vel(Wave_Corr_514, Flux_Corr_514, Oxygen, Dict_Data_514, 8000, 'Oxygen')       
Vel_521_O = get_line_vel(Wave_Corr_521, Flux_Corr_521, Oxygen, Dict_Data_521, 8500, 'Oxygen')  
Vel_comb_O = get_line_vel(Wave_Corr_comb, Flux_Corr_comb, Oxygen, Dict_Data_comb, 8000, 'Oxygen')
Vel_O = [Vel_420_O, Vel_421_O, Vel_422_O, Vel_423_O, Vel_424_O, Vel_428_O,\
          Vel_502_O, Vel_506_O, Vel_Gr_O, Vel_514_O, Vel_521_O, Vel_comb_O]
Vel_O_Err = [1000, 1000, 800, 800, 1000, 1000, 900, 900, 900, 800, 800, 800]


Vel_420_NaI = get_line_vel(Wave_Corr_420, Flux_Corr_420, NaI, Dict_Data_420, 7500, 'NaI')  # 5895 \AA
Vel_421_NaI = get_line_vel(Wave_Corr_421, Flux_Corr_421, NaI, Dict_Data_421, 7500, 'NaI')
Vel_422_NaI = get_line_vel(Wave_Corr_422, Flux_Corr_422, NaI, Dict_Data_422, 7500, 'NaI')   # Weird
Vel_423_NaI = get_line_vel(Wave_Corr_423, Flux_Corr_423, NaI, Dict_Data_423, 7500, 'NaI')   # Weird
Vel_424_NaI = get_line_vel(Wave_Corr_424, Flux_Corr_424, NaI, Dict_Data_424, 6900, 'NaI')   # Could be one on the right?
Vel_428_NaI = get_line_vel(Wave_Corr_428, Flux_Corr_428, NaI, Dict_Data_428, 8000, 'NaI')
Vel_502_NaI = get_line_vel(Wave_Corr_502, Flux_Corr_502, NaI, Dict_Data_502, 6800, 'NaI')
Vel_506_NaI = get_line_vel(Wave_Corr_506, Flux_Corr_506, NaI, Dict_Data_506, 6800, 'NaI')
Vel_Gr_NaI = get_line_vel(Wave_Corr_Gr, Flux_Corr_Gr, NaI, Dict_Data_Gr, 5500, 'NaI')  
Vel_514_NaI = get_line_vel(Wave_Corr_514, Flux_Corr_514, NaI, Dict_Data_514, 5500, 'NaI')       
Vel_521_NaI = get_line_vel(Wave_Corr_521, Flux_Corr_521, NaI, Dict_Data_521, 5500, 'NaI')  
Vel_comb_NaI = get_line_vel(Wave_Corr_comb, Flux_Corr_comb, NaI, Dict_Data_comb, 6000, 'NaI')  
Vel_NaI = [Vel_420_NaI, Vel_421_NaI, Vel_422_NaI, Vel_423_NaI, Vel_424_NaI, Vel_428_NaI,\
          Vel_502_NaI, Vel_506_NaI, Vel_Gr_NaI, Vel_514_NaI, Vel_521_NaI, Vel_comb_NaI]
Vel_NaI_Err = [500, 500, 500, 500, 500, 500, 500, 600, 500, 500, 500, 500]

Vel_420_MgI = get_line_vel(Wave_Corr_420, Flux_Corr_420, MgI, Dict_Data_420, 11000, 'MgI') # Can't tell  #4571 \AA
Vel_421_MgI = get_line_vel(Wave_Corr_421, Flux_Corr_421, MgI, Dict_Data_421, 11000, 'MgI') # Maybe
Vel_422_MgI = get_line_vel(Wave_Corr_422, Flux_Corr_422, MgI, Dict_Data_422, 11000, 'MgI')
Vel_423_MgI = get_line_vel(Wave_Corr_423, Flux_Corr_423, MgI, Dict_Data_423, 11000, 'MgI')
Vel_424_MgI = get_line_vel(Wave_Corr_424, Flux_Corr_424, MgI, Dict_Data_424, 11000, 'MgI')
Vel_428_MgI = get_line_vel(Wave_Corr_428, Flux_Corr_428, MgI, Dict_Data_428, 11000, 'MgI')
Vel_502_MgI = get_line_vel(Wave_Corr_502, Flux_Corr_502, MgI, Dict_Data_502, 11000, 'MgI')
Vel_506_MgI = get_line_vel(Wave_Corr_506, Flux_Corr_506, MgI, Dict_Data_506, 11000, 'MgI')
Vel_Gr_MgI = get_line_vel(Wave_Corr_Gr, Flux_Corr_Gr, MgI, Dict_Data_Gr, 10000, 'MgI')  
Vel_514_MgI = get_line_vel(Wave_Corr_514, Flux_Corr_514, MgI, Dict_Data_514, 10000, 'MgI')   # Can't tell    
Vel_521_MgI = get_line_vel(Wave_Corr_521, Flux_Corr_521, MgI, Dict_Data_521, 10000, 'MgI')   # Can't tell 
Vel_comb_MgI = get_line_vel(Wave_Corr_comb, Flux_Corr_comb, MgI, Dict_Data_comb, 10000, 'MgI')  # Can't tell 
Vel_MgI = np.array([Vel_420_MgI, Vel_421_MgI, Vel_422_MgI, Vel_423_MgI, Vel_424_MgI, Vel_428_MgI,\
          Vel_502_MgI, Vel_506_MgI, Vel_Gr_MgI, Vel_514_MgI, Vel_521_MgI, Vel_comb_MgI])  #0,9,10,11, nan idx
Vel_MgI = Vel_MgI.astype("float")
Vel_MgI[0] = np.nan
Vel_MgI[9] = np.nan
Vel_MgI[10] = np.nan
Vel_MgI[11] = np.nan
Vel_MgI_Err = [np.nan, 800, 500, 600, 600, 500, 700, 600, 800, np.nan, np.nan, np.nan]


Vel_420_SiII = get_line_vel(Wave_Corr_420, Flux_Corr_420, SiII, Dict_Data_420, 4000, 'SiII')  # 6347 \AA # Can't tell
Vel_421_SiII = get_line_vel(Wave_Corr_421, Flux_Corr_421, SiII, Dict_Data_421, 4000, 'SiII')  # Can't tell
Vel_422_SiII = get_line_vel(Wave_Corr_422, Flux_Corr_422, SiII, Dict_Data_422, 5500, 'SiII')
Vel_423_SiII = get_line_vel(Wave_Corr_423, Flux_Corr_423, SiII, Dict_Data_423, 5500, 'SiII')
Vel_424_SiII = get_line_vel(Wave_Corr_424, Flux_Corr_424, SiII, Dict_Data_424, 5500, 'SiII')
Vel_428_SiII = get_line_vel(Wave_Corr_428, Flux_Corr_428, SiII, Dict_Data_428, 4000, 'SiII')
Vel_502_SiII = get_line_vel(Wave_Corr_502, Flux_Corr_502, SiII, Dict_Data_502, 4000, 'SiII')
Vel_506_SiII = get_line_vel(Wave_Corr_506, Flux_Corr_506, SiII, Dict_Data_506, 4000, 'SiII')
Vel_Gr_SiII = get_line_vel(Wave_Corr_Gr, Flux_Corr_Gr, SiII, Dict_Data_Gr, 3500, 'SiII')  
Vel_514_SiII = get_line_vel(Wave_Corr_514, Flux_Corr_514, SiII, Dict_Data_514, 3800, 'SiII')       
Vel_521_SiII = get_line_vel(Wave_Corr_521, Flux_Corr_521, SiII, Dict_Data_521, 2400, 'SiII')  
Vel_comb_SiII = get_line_vel(Wave_Corr_comb, Flux_Corr_comb, SiII, Dict_Data_comb, 1100, 'SiII') 
Vel_SiII = np.array([Vel_420_SiII, Vel_421_SiII, Vel_422_SiII, Vel_423_SiII, Vel_424_SiII, Vel_428_SiII,\
          Vel_502_SiII, Vel_506_SiII, Vel_Gr_SiII, Vel_514_SiII, Vel_521_SiII, Vel_comb_SiII])
Vel_SiII = Vel_SiII.astype("float")
Vel_SiII[0] = np.nan
Vel_SiII[1] = np.nan
Vel_SiII_Err = [1000, 900, 600, 600, 600, 1100, 700, 1100, 700, 600, 700, 800]
  

#Vel_420_CaII = get_line_vel(Wave_Corr_420, Flux_Corr_420, CaII, Dict_Data_420, 10000, 'CaII')  # NIR triplet
#Vel_421_CaII = get_line_vel(Wave_Corr_421, Flux_Corr_421, CaII, Dict_Data_421, 10000, 'CaII')  
#Vel_422_CaII = get_line_vel(Wave_Corr_422, Flux_Corr_422, CaII, Dict_Data_422, 10000, 'CaII')  
#Vel_423_CaII = get_line_vel(Wave_Corr_423, Flux_Corr_423, CaII, Dict_Data_423, 10000, 'CaII')  
#Vel_424_CaII = get_line_vel(Wave_Corr_424, Flux_Corr_424, CaII, Dict_Data_424, 10000, 'CaII')  
#Vel_428_CaII = get_line_vel(Wave_Corr_428, Flux_Corr_428, CaII, Dict_Data_428, 10000, 'CaII')  
#Vel_502_CaII = get_line_vel(Wave_Corr_502, Flux_Corr_502, CaII, Dict_Data_502, 10000, 'CaII')  
#Vel_506_CaII = get_line_vel(Wave_Corr_506, Flux_Corr_506, CaII, Dict_Data_506, 10000, 'CaII')  
Vel_Gr_CaII = get_line_vel(Wave_Corr_Gr, Flux_Corr_Gr, CaII, Dict_Data_Gr, 10000, 'CaII')     
#Vel_514_CaII = get_line_vel(Wave_Corr_514, Flux_Corr_514, CaII, Dict_Data_514, 10000, 'CaII')    
#Vel_521_CaII = get_line_vel(Wave_Corr_521, Flux_Corr_521, CaII, Dict_Data_521, 10000, 'CaII')  
#Vel_comb_CaII = get_line_vel(Wave_Corr_comb, Flux_Corr_comb, CaII, Dict_Data_comb, 10000, 'CaII') 
#Vel_CaII = [Vel_420_CaII, Vel_421_CaII, Vel_422_CaII, Vel_423_CaII, Vel_424_CaII, Vel_428_CaII,\
#          Vel_502_CaII, Vel_506_CaII, Vel_Gr_CaII, Vel_514_CaII, Vel_521_CaII, Vel_comb_CaII]

    
def extract_sn_spec(SN_File, Line_Vel, Line_Vel_Err):
    
    All_Line_Vels = pd.read_csv(SN_File)
    Element_Vel = np.array(All_Line_Vels[Line_Vel]) * 1000
    Element_Vel_Err = np.array(All_Line_Vels[Line_Vel_Err]) * 1000
    Epoch = np.array(All_Line_Vels['epoch']) 
    
    return Epoch, Element_Vel, Element_Vel_Err


SN2011bm_Epoch, SN2011bm_Si_Vels, SN2011bm_Si_Vels_Err = extract_sn_spec('SN2011bm_spec.csv', 'Si', 'dSi') 
SN2011bm_Epoch, SN2011bm_Fe_Vels, SN2011bm_Fe_Vels_Err = extract_sn_spec('SN2011bm_spec.csv', 'Fe', 'dFe') 
SN2011bm_Epoch, SN2011bm_Na_Vels, SN2011bm_Na_Vels_Err = extract_sn_spec('SN2011bm_spec.csv', 'Na', 'dNa') 
SN2011bm_Epoch, SN2011bm_O_Vels,  SN2011bm_O_Vels_Err = extract_sn_spec('SN2011bm_spec.csv', 'O', 'dO')  
 
SN2007gr_Epoch, SN2007gr_Si_Vels, SN2007gr_Si_Vels_Err = extract_sn_spec('SN2007gr_spec.csv', 'Si', 'dSi') 
SN2007gr_Epoch, SN2007gr_Fe_Vels, SN2007gr_Fe_Vels_Err = extract_sn_spec('SN2007gr_spec.csv', 'Fe', 'dFe') 
SN2007gr_Epoch, SN2007gr_Na_Vels, SN2007gr_Na_Vels_Err = extract_sn_spec('SN2007gr_spec.csv', 'Na', 'dNa') 
SN2007gr_Epoch, SN2007gr_O_Vels,  SN2007gr_O_Vels_Err = extract_sn_spec('SN2007gr_spec.csv', 'O', 'dO') 

SN2004aw_Epoch, SN2004aw_Si_Vels, SN2004aw_Si_Vels_Err = extract_sn_spec('SN2004aw_spec.csv', 'Si', 'dSi') 
SN2004aw_Epoch, SN2004aw_Fe_Vels, SN2004aw_Fe_Vels_Err = extract_sn_spec('SN2004aw_spec.csv', 'Fe', 'dFe') 
SN2004aw_Epoch, SN2004aw_Na_Vels, SN2004aw_Na_Vels_Err = extract_sn_spec('SN2004aw_spec.csv', 'Na', 'dNa') 
SN2004aw_Epoch, SN2004aw_O_Vels,  SN2004aw_O_Vels_Err = extract_sn_spec('SN2004aw_spec.csv', 'O', 'dO')  

SN2002ap_Epoch, SN2002ap_Si_Vels, SN2002ap_Si_Vels_Err = extract_sn_spec('SN2002ap_spec.csv', 'Si', 'dSi') 
SN2002ap_Epoch, SN2002ap_Fe_Vels, SN2002ap_Fe_Vels_Err = extract_sn_spec('SN2002ap_spec.csv', 'Fe', 'dFe') 
SN2002ap_Epoch, SN2002ap_Na_Vels, SN2002ap_Na_Vels_Err = extract_sn_spec('SN2002ap_spec.csv', 'Na', 'dNa') 
SN2002ap_Epoch, SN2002ap_O_Vels,  SN2002ap_O_Vels_Err = extract_sn_spec('SN2002ap_spec.csv', 'O', 'dO')  

SN1998bw_Epoch, SN1998bw_Si_Vels, SN1998bw_Si_Vels_Err = extract_sn_spec('SN1998bw_spec.csv', 'Si', 'dSi') 
SN1998bw_Epoch, SN1998bw_Fe_Vels, SN1998bw_Fe_Vels_Err = extract_sn_spec('SN1998bw_spec.csv', 'Fe', 'dFe')
SN1998bw_Epoch, SN1998bw_Na_Vels, SN1998bw_Na_Vels_Err = extract_sn_spec('SN1998bw_spec.csv', 'Na', 'dNa')
SN1998bw_Epoch, SN1998bw_O_Vels, SN1998bw_O_Vels_Err = extract_sn_spec('SN1998bw_spec.csv', 'O', 'dO')

SN1994I_Epoch, SN1994I_Si_Vels, SN1994I_Si_Vels_Err = extract_sn_spec('SN1994I_spec.csv', 'Si', 'dSi') 
SN1994I_Epoch, SN1994I_Fe_Vels, SN1994I_Fe_Vels_Err = extract_sn_spec('SN1994I_spec.csv', 'Fe', 'dFe') 
SN1994I_Epoch, SN1994I_Na_Vels, SN1994I_Na_Vels_Err = extract_sn_spec('SN1994I_spec.csv', 'Na', 'dNa') 
SN1994I_Epoch, SN1994I_O_Vels,  SN1994I_O_Vels_Err = extract_sn_spec('SN1994I_spec.csv', 'O', 'dO')  


Si_Vels = [SN2011bm_Si_Vels, SN2007gr_Si_Vels, SN2004aw_Si_Vels, SN2002ap_Si_Vels, SN1998bw_Si_Vels, SN1994I_Si_Vels]
Fe_Vels = [SN2011bm_Fe_Vels, SN2007gr_Fe_Vels, SN2004aw_Fe_Vels, SN2002ap_Fe_Vels, SN1998bw_Fe_Vels, SN1994I_Fe_Vels]
Na_Vels = [SN2011bm_Na_Vels, SN2007gr_Na_Vels, SN2004aw_Na_Vels, SN2002ap_Na_Vels, SN1998bw_Na_Vels, SN1994I_Na_Vels]
O_Vels  = [SN2011bm_O_Vels, SN2007gr_O_Vels, SN2004aw_O_Vels, SN2002ap_O_Vels, SN1998bw_O_Vels, SN1994I_O_Vels]

Si_Vels_Errs = [SN2011bm_Si_Vels_Err, SN2007gr_Si_Vels_Err, SN2004aw_Si_Vels_Err, SN2002ap_Si_Vels_Err, SN1998bw_Si_Vels_Err, SN1994I_Si_Vels_Err]
Fe_Vels_Errs = [SN2011bm_Fe_Vels_Err, SN2007gr_Fe_Vels_Err, SN2004aw_Fe_Vels_Err, SN2002ap_Fe_Vels_Err, SN1998bw_Fe_Vels_Err, SN1994I_Fe_Vels_Err]
Na_Vels_Errs = [SN2011bm_Na_Vels_Err, SN2007gr_Na_Vels_Err, SN2004aw_Na_Vels_Err, SN2002ap_Na_Vels_Err, SN1998bw_Na_Vels_Err, SN1994I_Na_Vels_Err]
O_Vels_Errs  = [SN2011bm_O_Vels_Err, SN2007gr_O_Vels_Err, SN2004aw_O_Vels_Err, SN2002ap_O_Vels_Err, SN1998bw_O_Vels_Err, SN1994I_O_Vels_Err]
Elements_Epoch = [SN2011bm_Epoch, SN2007gr_Epoch, SN2004aw_Epoch, SN2002ap_Epoch, SN1998bw_Epoch, SN1994I_Epoch]

def plot_line_vel(El_Vel, Element, Vel_Err):
    
    plt.figure()
    plt.plot(Spec_MJD, El_Vel, '--o')
    plt.errorbar(Spec_MJD, El_Vel, yerr = Vel_Err)
    plt.xlabel('Days After Maximum Luminosity')
    plt.ylabel('Velocity [$Km\ s^{-1}$]')
    plt.title(f'Line Velocity Evolution of {Element} in SN2019cri')
    
    plt.show()
    plt.close()
    
    return 
    
#plot_line_vel(Vel_FeII, 'Fe', Vel_FeII_Err)      
#plot_line_vel(Vel_O, 'O', Vel_O_Err) 
#plot_line_vel(Vel_NaI, 'NaI', Vel_NaI_Err) 
#plot_line_vel(Vel_MgI, 'MgI', Vel_MgI_Err) 
#plot_line_vel(Vel_SiII, 'SiII', Vel_SiII_Err) 


def plot_compar_line_vel(My_El_Vel, Element, My_Vel_Err, Compar_El, Compar_El_Err):
    
    labels = ['SN2011bm', 'SN2007gr', 'SN2004aw', 'SN2002ap', 'SN1998bw', 'SN1994I']
    colours = ['brown', 'blue', 'green', 'pink', 'purple', 'orange']
    style = ['s', 'D', '*', 'd', 'x', 'v']
    
    plt.figure()
    plt.plot(Spec_MJD, My_El_Vel, 'o', color = 'red', label = 'SN2019cri', ms = 10)
    plt.errorbar(Spec_MJD, My_El_Vel, yerr = My_Vel_Err, color = 'red', elinewidth = 1)
    plt.xlabel('Days After Maximum Luminosity')
    plt.ylabel('Velocity [$Km\ s^{-1}$]')
    plt.title(f'Line Velocity Evolution of {Element} for different Type Ic')
    
    for i in range(len(Compar_El)):
        
        plt.scatter(Elements_Epoch[i], Compar_El[i], marker = style[i], label = labels[i], color = colours[i], lw = 1)
        plt.errorbar(Elements_Epoch[i], Compar_El[i], yerr = Compar_El_Err[i], color = colours[i], elinewidth = 1)
    
    plt.legend(loc = 'best')
    
    
    return
 
plot_compar_line_vel(Vel_FeII, 'Fe', Vel_FeII_Err, Fe_Vels, Fe_Vels_Errs)
plot_compar_line_vel(Vel_O, 'O', Vel_O_Err, O_Vels, O_Vels_Errs)
plot_compar_line_vel(Vel_NaI, 'NaI', Vel_NaI_Err, Na_Vels, Na_Vels_Errs)
plot_compar_line_vel(Vel_SiII, 'SiII', Vel_SiII_Err, Si_Vels, Si_Vels_Errs)


def compar_sn_corr(SN_Dict_Data):
    
    Flux_corr_1 = []
    Flux_Corr = []
    
    Wave_corr_1 = []
    Wave_Corr = []
    
    for i in range(len(SN_Dict_Data['Wavelength'])):
        
        
        for j in range(len(SN_Dict_Data['Wavelength'][i])):
            
            mw_dered_flux = dered(SN_Dict_Data['Wavelength'][i][j], SN_Dict_Data['Flux'][i][j], 3.1, SN_Dict_Data['MW_Red'])
            
            redshift_corr_wavelength = redshift_correction(SN_Dict_Data['Wavelength'][i][j], SN_Dict_Data)
            
            host_dered_flux = dered(redshift_corr_wavelength, mw_dered_flux, 3.1, SN_Dict_Data['Host_Red'])
            
            #print(mw_dered_flux, redshift_corr_wavelength, host_dered_flux)
            
            Flux_corr_1.append(host_dered_flux)
            Wave_corr_1.append(redshift_corr_wavelength)
    
    
    for k in range(len(SN_Dict_Data['Wavelength'])):
        
       Flux_Corr.append(Flux_corr_1[0:len(SN_Dict_Data['Flux'][k])])
       Wave_Corr.append(Wave_corr_1[0:len(SN_Dict_Data['Wavelength'][k])])
        
    return Flux_Corr, Wave_Corr




#Flux_Corr_SN2007gr, Wave_Corr_SN2007gr = compar_sn_corr(SN2007gr_Dict_Data)
#Flux_Corr_SN2011bm, Wave_Corr_SN2011bm = compar_sn_corr(SN2011bm_Dict_Data)














# Comparing SN Light Curves

       
My_R_Band_LC = pd.read_csv('LC_DATA/R_Band_LC.csv')
RMags = list(My_R_Band_LC['PSF_Bkg_Mag'])
RMags_Err = list(My_R_Band_LC['PSF_Mag_Error'])
RMags_MJD = np.array(My_R_Band_LC['MJD']) * (1 / (1 + 0.041))

ZTF_MJD_, ZTF_Mags_, ZTF_Mags_Err = np.loadtxt('lightcurve_ZTF.txt', unpack = True, usecols = (0,3,4))

RS_ZTF = ZTF_MJD_ * (1 / (1 + 0.041))
RS_SN2007gr = np.array(SN2007gr_Dict_Data['MJD_Phot']) * (1 / (1 + SN2007gr_Dict_Data['Redshift']))
RS_SN2007ru = np.array(SN2007ru_Dict_Data['MJD_Phot']) * (1 / (1 + SN2007ru_Dict_Data['Redshift']))
RS_SN1994I = np.array(SN1994I_Dict_Data['MJD_Phot']) * (1 / (1 + SN1994I_Dict_Data['Redshift']))
RS_SN2004aw = np.array(SN2004aw_Dict_Data['MJD_Phot']) * (1 / (1 + SN2004aw_Dict_Data['Redshift']))
RS_SN2011bm = np.array(SN2011bm_Dict_Data['MJD_Phot']) * (1 / (1 + SN2011bm_Dict_Data['Redshift']))


def mjd_alignment(SN_Mag, SN_MJD):
    
    Min_Pos_Idx = SN_Mag.index(min(SN_Mag))
    MJD_Corr = SN_MJD[Min_Pos_Idx]
    
    return MJD_Corr


def mjd_correction(SN, idx):
    
    My_MJD_Corr = mjd_alignment(RMags, RMags_MJD)
    ZTF_Corr = mjd_alignment(list(ZTF_Mags_), RS_ZTF)
    SN2007gr_MJD_Corr = mjd_alignment(SN2007gr_Dict_Data['Mag'], RS_SN2007gr)
    SN2007ru_MJD_Corr = mjd_alignment(SN2007ru_Dict_Data['Mag'], RS_SN2007ru)
    SN1994I_MJD_Corr = mjd_alignment(SN1994I_Dict_Data['Mag'], RS_SN1994I)
    SN2004aw_MJD_Corr = mjd_alignment(SN2004aw_Dict_Data['Mag'], RS_SN2004aw)
    SN2011bm_MJD_Corr = mjd_alignment(SN2011bm_Dict_Data['Mag'], RS_SN2011bm)

    MJD_Corrs = [ZTF_Corr, SN2007gr_MJD_Corr, SN2007ru_MJD_Corr, SN1994I_MJD_Corr, SN2004aw_MJD_Corr, SN2011bm_MJD_Corr]
    MJD_corr_val = []

    for corr in MJD_Corrs:
    
        MJD_corr_val_ta = My_MJD_Corr - corr
        MJD_corr_val.append(MJD_corr_val_ta)
    
    True_MJD = SN + MJD_corr_val[idx]
    
    return True_MJD


def mag_correction(SN_Mag):
    
    dered_mag_val = []
    
    for i in range(len(SN_Mag)):
        
        dered_mag_val_ta = dered_mag(6410, SN_Mag[i], E = 1e-5)
        dered_mag_val.append(dered_mag_val_ta)
        
    return dered_mag_val



def plot_compar_sn():
    
    Dered_RMag = mag_correction(RMags)
    Dred_ZTF = mag_correction(ZTF_Mags_)
    
    Comb_RMag_App = sorted(np.append(Dred_ZTF, Dered_RMag))
    
    Dered_Mag_SN2007gr = mag_correction(SN2007gr_Dict_Data['Mag'])
    Dered_Mag_SN2007ru = mag_correction(SN2007ru_Dict_Data['Mag'])
    Dered_Mag_SN1994I = mag_correction(SN1994I_Dict_Data['Mag'])
    Dered_Mag_SN2004aw = mag_correction(SN2004aw_Dict_Data['Mag'])
    Dered_Mag_SN2011bm = mag_correction(SN2011bm_Dict_Data['Mag'])

    Abs_Mag_SN2007gr = Dered_Mag_SN2007gr - np.array(SN2007gr_Dict_Data['Distance_Modulus'])
    Abs_Mag_SN2004ru = Dered_Mag_SN2007ru - np.array(SN2007ru_Dict_Data['Distance_Modulus'])
    Abs_Mag_SN1994I = Dered_Mag_SN1994I - np.array(SN1994I_Dict_Data['Distance_Modulus'])
    Abs_Mag_SN2004aw = Dered_Mag_SN2004aw - np.array(SN2004aw_Dict_Data['Distance_Modulus'])
    Abs_Mag_SN2011bm = Dered_Mag_SN2011bm - np.array(SN2011bm_Dict_Data['Distance_Modulus'])


    RMag_Dis = Distance(z = 0.041, cosmology = WMAP5) / u.Mpc
    Abs_RMag = np.array(Dered_RMag) - (5 * np.log10(RMag_Dis*10**6)) + 5

    NaN_Idx = np.where(np.isnan(Abs_RMag))
    Abs_RMag = [x for x in Abs_RMag if ~np.isnan(x)]
    RMags_MJD_ = np.delete(RMags_MJD, [NaN_Idx[0]])
    RMags_Err_ = np.delete(RMags_Err, [NaN_Idx[0]])



    ZTF_Abs_Mag = np.array(Dred_ZTF) - (5 * np.log10(RMag_Dis*10**6)) + 5
    
    NaN_Idx_ZTF = np.where(np.isnan(ZTF_Abs_Mag))
    ZTF_Abs_Mag = [x for x in ZTF_Abs_Mag if ~np.isnan(x)]
    ZTF_Err_ = np.delete(ZTF_Mags_Err, [NaN_Idx_ZTF[0]])
    RS_ZTF_ = np.delete(RS_ZTF, [NaN_Idx_ZTF[0]])
        
        

    True_MJD_SN2007gr = mjd_correction(RS_SN2007gr, 1)
    True_MJD_SN2007ru = mjd_correction(RS_SN2007ru, 2)
    True_MJD_SN1994I = mjd_correction(RS_SN1994I, 3)
    True_MJD_SN2004aw = mjd_correction(RS_SN2004aw, 4)
    True_MJD_SN2011bm = mjd_correction(RS_SN2011bm, 5)
        
    My_MJD_Corr = mjd_alignment(RMags, RMags_MJD)
    
    Comb_RMag_MJD = np.append(RS_ZTF_, RMags_MJD_)
    Comb_RMag = np.append(np.array(ZTF_Abs_Mag), Abs_RMag)
    Comb_RMag_Err = np.append(np.array(ZTF_Err_), np.array(RMags_Err_))
    
    b = [(Comb_RMag_MJD[i], Comb_RMag[i]) for i in range(0, len(Comb_RMag))]
    c = sorted(b)
    
    Mag_Comb_Err = sorted(Comb_RMag_Err)
    
    True_Comb_MJD = []
    RMag_Comb = []
    
    for i in range(len(c)):
        
        True_Comb_MJD.append(c[i][0])
        RMag_Comb.append(c[i][1])

    
    rcParams['figure.figsize'] = 10, 8
    #print(len(Comb_RMag_App), len(Comb_RMag))
    #print(Comb_RMag_App[1] - Abs_RMag[1])
    print(min(RMag_Comb))
    
    
    plt.figure()
    plt.plot(True_MJD_SN2007gr - My_MJD_Corr, Abs_Mag_SN2007gr, '--o', label = 'SN2007gr', ms = 4, lw = 1, mew = 1)
    plt.plot(True_MJD_SN2007ru - My_MJD_Corr, Abs_Mag_SN2004ru, '--o', label = 'SN2007ru', ms = 4, lw = 1, mew = 1)
    plt.plot(True_MJD_SN1994I - My_MJD_Corr,Abs_Mag_SN1994I, '--o', label = 'SN1994I', ms = 4, lw = 1, mew = 1)
    plt.plot(True_MJD_SN2004aw - My_MJD_Corr, Abs_Mag_SN2004aw, '--o', label = 'SN2004aw', ms = 4, lw = 1, mew = 1)
    plt.plot(True_MJD_SN2011bm - My_MJD_Corr, Abs_Mag_SN2011bm , '--o', label = 'SN2011bm', ms = 4, alpha = 0.4, lw = 1, color = 'darkblue', mew = 1)
    plt.plot(True_Comb_MJD - My_MJD_Corr, RMag_Comb, '--*', label = 'SN2019cri', ms = 10, color = 'gold', lw = 1.5)
    plt.errorbar(True_Comb_MJD - My_MJD_Corr, RMag_Comb, yerr = np.array(Mag_Comb_Err), linestyle = 'None', color = 'black')
    plt.minorticks_on()
    plt.xlim(-50,200)
    plt.ylim(-19, -12.5)
    plt.xlabel('Days After r-Band Maximum Light')
    plt.ylabel('Absolute Magnitude')
    plt.title('Plot of Type Ic Supernova Light Curves (r-Band)')
    plt.legend(loc = 'best', prop={'size': 12}, ncol = 2)
    plt.gca().invert_yaxis()
    plt.show()
    
    return True_MJD_SN2007gr, True_MJD_SN2007ru, True_MJD_SN1994I, True_MJD_SN2004aw, True_MJD_SN2011bm, True_Comb_MJD, My_MJD_Corr,\
           Mag_Comb_Err, Comb_RMag_App, Dered_Mag_SN2007gr, Dered_Mag_SN2007ru, Dered_Mag_SN1994I, Dered_Mag_SN2004aw, Dered_Mag_SN2011bm


True_MJD_SN2007gr, True_MJD_SN2007ru, True_MJD_SN1994I, True_MJD_SN2004aw, True_MJD_SN2011bm, True_Comb_MJD,\
My_MJD_Corr, Mag_Comb_Err, Comb_RMag_App, Dered_Mag_SN2007gr, Dered_Mag_SN2007ru, Dered_Mag_SN1994I, Dered_Mag_SN2004aw, Dered_Mag_SN2011bm = plot_compar_sn()

 
  
 
    