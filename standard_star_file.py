# -*- coding: utf-8 -*-
"""
Created on Wed May 27 12:01:02 2020

@author: Bill Kamtchou

"""

import csv
import numpy as np
import pandas as pd




def std_star_file_processing(path):
    
    Std_Star_Data = pd.read_csv (f'{path}')
    
    Reduced_Std_Star_Data_1 = Std_Star_Data.drop(Std_Star_Data.query(f'cl == 3').index)
    # Removing all class 3 objects as they as they are galaxies and we are interested in stars
    
    Reduced_Std_Star_Data = Reduced_Std_Star_Data_1.drop(Reduced_Std_Star_Data_1.query('rmag > 20').index)
    # Removing any stars higher than 20th order magnitude
    
    RA = np.array(Reduced_Std_Star_Data['_RAJ2000'])
    DEC = np.array(Reduced_Std_Star_Data['_DEJ2000'])
    
    Position = [(RA[i], DEC[i]) for i in range(0, len(RA))] # Creating a coordinate point of RA and DEC
    
    Umag = np.array(Reduced_Std_Star_Data['umag'])
    Err_Umag = np.array(Reduced_Std_Star_Data['e_umag'])
    
    Gmag = np.array(Reduced_Std_Star_Data['gmag'])
    Err_Gmag = np.array(Reduced_Std_Star_Data['e_gmag'])
    
    Rmag = np.array(Reduced_Std_Star_Data['rmag'])
    Err_Rmag = np.array(Reduced_Std_Star_Data['e_rmag'])
    
    Imag = np.array(Reduced_Std_Star_Data['imag'])
    Err_Imag = np.array(Reduced_Std_Star_Data['e_imag'])
    
    Zmag = np.array(Reduced_Std_Star_Data['zmag'])
    Err_Zmag = np.array(Reduced_Std_Star_Data['e_zmag'])
    
    
    
    with open(f'Cleaned_Standard_Star_File.csv', 'w') as out_file:
            
        csv_writer = csv.writer(out_file, delimiter = ',')
        csv_writer.writerow(['Label', 'RA', 'DEC', 'Rmag', 'Err_Rmag', 'Imag', 'Err_Imag', 'Zmag', 'Err_Zmag', 'Gmag', 'Err_Gmag', 'Umag', 'Err_Umag'])
            
        for n in range(len(RA)):
            
            csv_writer.writerow([f'SSFS{n}', RA[n], DEC[n], Rmag[n], Err_Rmag[n], Imag[n], Err_Imag[n], Zmag[n], Err_Zmag[n], Gmag[n], Err_Gmag[n], Umag[n], Err_Umag[n]])
    
    return RA, DEC, Position, Umag, Err_Umag, Gmag, Err_Gmag, Rmag, Err_Rmag, Imag, Err_Imag, Zmag, Err_Zmag
    
    

######################################################################################
#  All the function needs is the path to the standard star file in a csv format.
#  ie.g. Magnitude_Data = classifier('RAW_DATA/standard_star_file.csv')
#
#  Data from the csv file can be specified by the user by acessing the csv file like a
#  dictionary. e.g. RA = np.array(Reduced_Std_Star_Data['_RAJ2000'])
#  This function returns all the magnitudes from the R,U,G,I and Z bands with their 
#  associated error and their corresponding RA and DEC. It will also output a cleaned
#  and labelled standard star file.
######################################################################################
