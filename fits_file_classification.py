# -*- coding: utf-8 -*-
"""
Created on Tue May 25 12:27:57 2020

@author: Bill Kamtchou
"""

import os                                          # Importing all the necessary functions
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename



def classifier(fits_files_path):  # This code classifies fits files header data into a dictionary
    
    All_Fits_Files = os.listdir(fits_files_path)
    Headers = ['DATE-OBS', 'MJD', 'FILTER1', 'EXPTIME']   # The headers can be specified at the users discreton
    Header_Data = []
    n = len(Headers)

    for i in range(len(All_Fits_Files)):
    
        raw_image_file = get_pkg_data_filename(f'RAW_DATA/{All_Fits_Files[i]}')
        image_hdu = fits.open(raw_image_file)
    
        for i in range(len(Headers)):

            Header_Result = image_hdu[0].header[Headers[i]]
            Header_Data.append(Header_Result)

        image_hdu.close()


    Header_Data = [Header_Data[k:k+n] for k in range(0, len(Header_Data), n)]

    Header_Data_Dictionary = dict(zip(All_Fits_Files, Header_Data))
    
    return Header_Data_Dictionary


######################################################################################
#  All the function needs is the path to the fits file and
#  is called as such e.g. Header_Data_Dictionary = classifier('RAW_DATA/')
#
#  Moreover, vlaues of different fits files can be accessed as such e.g.
#  print(Header_Data_Dictionary['name_of_file'][index]). The index (which starts at 0)  
#  corresponds to the position of the header being accessed in Headers. 
#  E.g. index = 0 is 'DATE-OBS'
######################################################################################


