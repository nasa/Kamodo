# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 11:40:54 2022

@author: rringuet
"""

from glob import glob
import os

# first version: all files had the same name regardless of grid type or
# hemisphere, so had to run this four times to get all the names correct.
'''
new_ending = '_Ndf.txt'  # N hemisphere default grid
file_dir = 'C:/Users/rringuet/Kamodo_Data/SuperDARN/fullday/'
# new_ending = '_Sdf.txt'  # S hemisphere default grid
# file_dir = no data for this
# new_ending = '_Nea.txt'  # N hemisphere equal-area grid
# file_dir = 'C:/Users/rringuet/Kamodo_Data/SuperDARN/fullday-equal/'
# new_ending = '_Sea.txt'  # S hemisphere equal-area grid
# file_dir = 'C:/Users/rringuet/Kamodo_Data/SuperDARN/fullday-equal-sh/'
'''
# second version: all files come in tarballs. The files are named for the
# hemisphere and grid type, but with the wrong endings
# 'eq' instead of 'ea' and 'uni' instead of 'df'. This causes a problem in
# later code because of the unequal endings. Fixing...
file_dir = 'C:/Users/rringuet/Kamodo_Data/SuperDARN/SDconvection/'

import tarfile
tar_files = glob(file_dir+'*.tar')
for file in tar_files:
    tar = tarfile.open(file)
    tar.extractall(file_dir)
    tar.close()

files = glob(file_dir+'*.txt')
for file in files:
    if 'eq' in file:
        new_filename = file[:-6] + 'ea.txt'
    elif 'uni' in file:
        new_filename = file[:-7] + 'df.txt'
    # copy data over to new file
    file_object = open(file)
    data = file_object.readlines()
    file_object.close()
    file_object = open(new_filename, 'w')
    for item in data:
        file_object.write(item)
    file_object.close()
    
    # read in both to compare
    file1 = open(file)
    data1 = file1.readlines()
    file1.close()
    file2 = open(new_filename)
    data2 = file2.readlines()
    file2.close()
    
    # remove old file if identical
    if data1 == data2:
        os.remove(file)
