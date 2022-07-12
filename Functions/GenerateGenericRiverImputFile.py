# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 19:05:39 2021

@author: PradoDomercq
"""

#Script to generate the Generic River imput file including
#river sections and compartments of constant dimensions (depth_m,length_m,width_m),
#G,T_K,discharge_m3s and SPM_mgL

#Import modules 
import math

#Define river 
compartments = ["surface", "flowingWater", "stagnantWater", "sediment"]
compDepth= [0.1, 4.9, 0.4, 0.02]
widthRS=80
#all alculations are in m!
lengthTotal_km = 1000 #Total modelled length in Km
numRS= 20
G = "10"
T_K = "287.15" #asummed the same for all compartments 
SPM_mgL= ["30", "50", "70", "90"]
vflow_m_s=1.30

indexList= list(range(numRS*4))
#generate list of river sections (RS) with same dimensions and list of their corresponding lengths 
#according to Praetorius et al subdivission.
RS_length_m= int(lengthTotal_km*1000/numRS)
listRS=list(range(numRS))
RSlengths_list= [RS_length_m]*numRS

#Cummulative lengths
RSCumLength=[RS_length_m]
for l in range (1,len(RSlengths_list)):
    RSCumLength.append(RSlengths_list[l]+RSCumLength[l-1])




##CREATE RIVER IMPUT FILE

file_name = "compartmentsGenericRiverSec_prop.txt"

out_file = open(file_name, "w")
out_file.write("riverSection,nameRS,compartment,compType,depth_m,length_m,volume_m3,width_m,G,T_K,vFlow_m_s,SPM_mgL\n")
for rs in range(numRS):
    for comp in range(len(compartments)):
        # if comp == 0:
        #     index = rs*4
        # elif comp== 1:
        #     index = rs*4+1
        # elif comp== 2:
        #      index = rs*4+2
        # else:
        #      index = rs*4+3
        out_file.write(str(rs)+","+ "RS"+str(rs)+ "," + str(comp+1) +","+ compartments[comp]+","+ str(compDepth[comp])+ "," + str(RS_length_m)+ "," + str(compDepth[comp]*RS_length_m*widthRS)+ ","+str(widthRS)+ ","+ G +","+T_K +","+ str(vflow_m_s)+","+ SPM_mgL[0] +"\n")    
out_file.close()

