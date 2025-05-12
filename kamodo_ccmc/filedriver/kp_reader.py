"""
Created on Wed Oct  2 19:34:52 2024

@author: Josh Pettit, GMU/NASA GSFC

Computes KP daily mean for FileDriver code
from a file (kp_2000-2025.txt built from data
pulled from https://kp.gfz.de/en/data
"""

import numpy as np

input_file_dir = ''

with open(input_file_dir+'kp_2000-2025.txt') as f:
    contents = f.readlines()
    
kp_total = np.zeros(len(contents))  
date4 = ["" for x in range(len(contents))]

for i in range(0, len(contents)-1):
    
    kp_temp = str(contents[i])
    kp_total[i]= kp_temp[47:52]
    date1 = kp_temp[0:4]
    date2 = kp_temp[5:7]
    date3 = kp_temp[8:10]
    date4[i] = [date1, date2, date3]
    
number_of_days = int(len(contents)/8)
kp_daily = np.zeros(number_of_days)    
date_daily = ["" for x in range(len(kp_daily))]

for k in range(0, len(kp_daily)):
    
    if k == 0:
        
        j = 0
        m = 8
        kp_temp2 = kp_total[j:m]
        
        kp_daily[k] = np.mean(kp_temp2)
        date_daily[k] = date4[k]
        
        j = 9
        m = 17
                
    if k != 0:     
        
        kp_temp2 = kp_total[j:m]
        kp_daily[k] = np.mean(kp_temp2)
        date_daily[k] = date4[j]
        
        j = j+8
        m = m+8
    
np.savez('KP_Daily.npz', date_daily=date_daily, kp_daily=kp_daily)
