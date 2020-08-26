# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 10:26:35 2020

@author: Bill Kamtchou
"""


from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
from astropy.stats import mad_std
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch
import astroscrappy
from fits_file_classification import classifier
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
from photutils import Background2D, MedianBackground







def band_classification():
    
    Header_Data_Dictionary = classifier('RAW_DATA/')

    Bands = ['SDSS-R', 'SDSS-I', 'SDSS-Z', 'SDSS-G', 'SDSS-U']

    R_Band_Files = []
    I_Band_Files = []
    Z_Band_Files = []
    G_Band_Files = []
    U_Band_Files = []
    All_Band_Files = [R_Band_Files, I_Band_Files, Z_Band_Files, G_Band_Files, U_Band_Files]

    for i in range(len(Bands)):

        {All_Band_Files[i].append(file) for (file, data) in Header_Data_Dictionary.items() if Bands[i] in data}
        
    
    return Header_Data_Dictionary, R_Band_Files, I_Band_Files, Z_Band_Files, G_Band_Files, U_Band_Files



Header_Data_Dictionary, R_Band_Files, I_Band_Files, Z_Band_Files, G_Band_Files, U_Band_Files =  band_classification()



def background_reduction(List):
    
    Raw_Image_Data = []
    Reduced_Image_Data = []
    Bkg = []
    Bkg_Sigma = []
    Median = []
    Std = []
    Norm = ImageNormalize(stretch=SqrtStretch())

    for i in range(len(List)):
    
        Image_File = get_pkg_data_filename(f'RAW_DATA/{List[i]}')
        crmask, Raw_Image_Data_1 = astroscrappy.detect_cosmics(fits.getdata(Image_File, ext = 0))
        Raw_Image_Data.append(Raw_Image_Data_1)

        Mask = (Raw_Image_Data_1 == 0)
        Sigma_Clip = SigmaClip(sigma = 3.0)
        Bkg_Estimator = MedianBackground()
        
        Bkg_ta = Background2D(Raw_Image_Data_1, (25, 25), filter_size = (3, 3), 
                   sigma_clip = Sigma_Clip, bkg_estimator = Bkg_Estimator, mask = Mask)
        
        Bkg.append(Bkg_ta)
    
        Bkg_Sigma_ta = mad_std(Raw_Image_Data_1)
        Bkg_Sigma.append(Bkg_Sigma_ta)
        
        Mean, Median_ta, Std_ta = sigma_clipped_stats(Raw_Image_Data_1, sigma = 3.0)
        Median.append(Median_ta)
        Std.append(Std_ta)
        
        Reduced_Image_Data_to_append = Raw_Image_Data_1 - Bkg_ta.background
        Reduced_Image_Data.append(Reduced_Image_Data_to_append)


        #plt.figure()
        #plt.imshow(Raw_Image_Data_1, cmap = 'gray', origin = 'lower', interpolation = 'bicubic', norm = Norm, vmin = 100, vmax = 1000)
        #plt.colorbar()

        #plt.figure()
        #plt.imshow(Bkg_ta.background, origin='lower', cmap='gray', interpolation = 'bicubic')
        #plt.colorbar()
    
        #plt.figure()
        #plt.imshow(Reduced_Image_Data_to_append, norm = Norm, origin = 'lower', cmap = 'gray', interpolation = 'bicubic', vmin = 1, vmax = 1000)
        #plt.colorbar()
        
        #plt.show()
    
  
    return  Raw_Image_Data, Reduced_Image_Data, Bkg, Bkg_Sigma, Median, Std



Raw_Image_Data, Reduced_Image_Data, Bkg, Bkg_Sigma, Median, Std = background_reduction(Z_Band_Files) # Change to band needed



#####################################################################################################
#  The purpose of this scipt is to firstly classify all the exposure fits files
#  into the correct ban d and hence form an array. The second function is to 
#  perform background subtraction using code from Photutils.

#  All this script needs is the loctaion where the fits files are (line 67) and the sorted band 
#  file of exposures, e.g. R_Band_Files. This folder must only contain the exposure fits files. 
#  The scipt will return the indexable raw image data,reduced image data and all the variables 
#  needed for future photometry for all the exposure fits files in the folder. 
#  E.g. Reduced_Image_Data[0] is the reduced image data for the expsoure and so on.
#####################################################################################################