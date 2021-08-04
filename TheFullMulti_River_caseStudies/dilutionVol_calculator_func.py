# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 14:07:59 2021

@author: PradoDomercq
"""
import numpy as np

#Create vector of compartmet volumes per RS for each species
sp_RScomp_vol_m3=[]
def volumesVector(Clist,compartments_prop):
    for sp3 in Clist:
        rivS= sp3[2:-3]
        cmp3=sp3[-3]
        #generate location index from compartment properties
        vol_index = np.where((compartments_prop['riverSection']== int(rivS)) & (compartments_prop['compartment']== int(cmp3)))[0][0]
        
        sp_RScomp_vol_m3.append(compartments_prop.depth_m.loc[vol_index]*compartments_prop.length_m.loc[vol_index]*compartments_prop.width_m.loc[vol_index])
    return sp_RScomp_vol_m3
        
       