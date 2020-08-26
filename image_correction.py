# -*- coding: utf-8 -*-
"""
Created on Mon May 25 11:20:44 2020

@author: Bill Kamtchou
"""

from matplotlib.colors import LogNorm # Importing functions that are needed from different libraries.             
import matplotlib.pyplot as plt   
import numpy as np
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch
from astropy.utils.data import get_pkg_data_filename
import astroalign as aa




def stack_image(image_name_1, image_name_2, new_file_name):
    
    File_List = [image_name_1, image_name_2]
    
    with fits.open(f'RAW_DATA/{File_List[0]}') as image_hdu_1:
    
        Data_1 = image_hdu_1[0].data
        image_hdu_1.close()
        
        
    with fits.open(f'RAW_DATA/{File_List[1]}') as image_hdu_2:
    
        Data_2 = image_hdu_2[0].data  # Creating a numpy array of the exposure
        image_hdu_2.close()

    
        Image_Data = [Data_1, Data_2]
    
        Stacked_Image_Data = np.median(Image_Data, axis = 0)  # Stacking

        plt.style.use(astropy_mpl_style)
        Norm = LogNorm()

        plt.figure(1)
        plt.imshow(Image_Data[0], cmap = 'gray', interpolation = 'bicubic', norm = Norm, vmin = 10, vmax = 4000)
        plt.colorbar()

        plt.figure(2)
        plt.imshow(Image_Data[1], cmap = 'gray', interpolation = 'bicubic', norm = Norm)
        plt.colorbar()

        plt.figure(3)
        plt.imshow(Stacked_Image_Data, cmap = 'gray', interpolation = 'bicubic', norm = Norm)
        plt.colorbar()
    
    
        image_hdu_2[0].data = Stacked_Image_Data  # Writing to a new fits file
        image_hdu_2.writeto(f'STACKED/{new_file_name}', overwrite = True)
        image_hdu_2.close()

    return 


def align_image(image_name_1, image_name_2, new_file_name):
    
    
    File_List = [image_name_1, image_name_2]
    
    with fits.open(f'RAW_DATA/{File_List[0]}') as image_hdu_1:
    
        Data_1 = image_hdu_1[0].data
        image_hdu_1.close()
        
        
    with fits.open(f'RAW_DATA/{File_List[1]}') as image_hdu_2:
    
        Data_2 = image_hdu_2[0].data  # Creating a numpy array of the exposure
        image_hdu_2.close()

    
    Image_Data = [Data_1, Data_2]
    
    transf, (source_list, target_list) = aa.find_transform(Image_Data[0], Image_Data[1])
    
    Aligned_Image_Data, footprint = aa.apply_transform(transf, Image_Data[0], Image_Data[1]) # Aligning 


    plt.style.use(astropy_mpl_style)
    Norm = ImageNormalize(stretch = SqrtStretch())

    plt.figure(1)
    plt.imshow(Image_Data[0], cmap = 'gray', interpolation = 'bicubic', norm = Norm, vmin = 700, vmax = 800)
    plt.colorbar()

    plt.figure(2)
    plt.imshow(Image_Data[1], cmap = 'gray', interpolation = 'bicubic', norm = Norm)
    plt.colorbar()

    plt.figure(3)
    plt.imshow(Aligned_Image_Data, cmap = 'gray', interpolation = 'bicubic', norm = Norm)
    plt.colorbar()
    
    
    with fits.open(f'RAW_DATA/{File_List[0]}', mode = 'update') as image_hdu_3:
    
        image_hdu_3[0].data = Aligned_Image_Data
        image_hdu_3.flush()
        image_hdu_3.writeto(f'AlIGNED/{new_file_name}', overwrite = True)  # Writing to a new fits file
    

    return 


#################################################################
# The purpose of this script is to align ans stack images that 
# have more than one exposure to get a better flux count and  
# hence have better aperture and PSF photometry. The inputs are
# simply the path to the two fits files and a new file name.
#################################################################



