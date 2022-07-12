# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 10:44:30 2021

@author: PradoDomercq
"""

"""Estimate RATE CONSTANTS"""
##Estimate rate constants for all "species" = combination of RS*riverComp*MPforms*sizeBins (len (Clist))
import numpy as np
import pandas as pd
from Functions.RC_GeneratorRiver_v2_0 import*
from Functions.objectGenerationRiver_func import*

#class to generate environmental
from objects.EnvCompartmentRiver import*
#import file storing required constants
from helpers.GlobalConstants import *
#import classes to generate objects
from objects.Particulates import Particulates #class to generate MP and SPM objects
from objects.ParticulatesBF import ParticulatesBF #class to generate MP and SPM objects
from objects.ParticulatesSPM import ParticulatesSPM #class to generate MP and SPM objects

def RC_estimation_function_v2_0(processList,CombList,Clist,MP_prop,compartments_prop,process_df,numberRS, composition,mode2, mode, date,riverComp,MPforms,sizeBin,river_flows):
    
    RC_df = pd.DataFrame(index=processList , columns= CombList)
        
    if numberRS <= 10:
        print("Less than 10 RS")
    #if len(Clist) <= 10*len(riverComp)*len(MPforms)*len(sizeBin):
        for spc in Clist:
           
            riverSection = spc[2]
            compartment= spc[3]
            aggState = spc[4]
            sizeBinIdx = spc[5]
        
            #Generate particle OBJECTS MPs and SPM for the corresponding size bin
            
            #generate location index from imput tables
            MP_index = np.where((MP_prop['composition']== composition) & (MP_prop['sizeBin']== sizeBinIdx))[0][0]
            SPM_index = np.where((MP_prop['name']== "SPM1") & (MP_prop['sizeBin']== sizeBinIdx))[0][0]
            comp_index = np.where((compartments_prop['compartment']== int(compartment)) & (compartments_prop['riverSection']== int(riverSection)))[0][0]
            rs_indexs=np.where(compartments_prop['riverSection']== int(riverSection))[0]
            #for processes
            idx = np.where((process_df['aggState']==aggState) & (process_df['compartment']== int(compartment)))[0][0]

            
            #Generate compartment objects corresponding to the River Section
            [surface, flowingWater, stagnantWater, sediment]= preProcessLayers(mode, compartments_prop, date, rs_indexs,riverSection)
            #Create a dictionary with the environmental compartments objects 
            comp_dict = {"1": surface, "2": flowingWater, "3": stagnantWater, "4": sediment }

            #Generate particle object:

            #generate SPM object(s)
            #SPM_index = 11 #need to move SPM in own input file
            SPM1 = Particulates(MP_prop, SPM_index)
            SPM1.calc_volume()
            SPM1.calc_numConc(compartments_prop.SPM_mgL[comp_index], 0)

            #Generate particle
            if aggState== "A":
                particle= Particulates(MP_prop, MP_index)
                particle.calc_volume()
            elif aggState == "B":
                MP1= Particulates(MP_prop, MP_index)
                MP1.calc_volume()
                #B: heteroaggregated (MP attached to suspended particulate matter (SPM)) 
                particle = ParticulatesSPM("MP1-SPM", MP1, SPM1) 
                particle.calc_volume(MP1, SPM1)
            
            elif aggState == "C":
                MP1= Particulates(MP_prop, MP_index)
                MP1.calc_volume()
                #C: biofiolm-covered (MP with biofilm (BF) layer on surface)
                particle = ParticulatesBF("MP1-BF", MP1, 1388, 5e-6) 
                particle.calc_volume()
            
            elif aggState == "D":
                MP1= Particulates(MP_prop, MP_index)
                MP1.calc_volume()
                #C: biofiolm-covered (MP with biofilm (BF) layer on surface)
                MP1_BF = ParticulatesBF("MP1-BF", MP1, 1388, 5e-6) 
                MP1_BF.calc_volume()
                #D: biofilm-heteroaggregated (MP with BF layer attached to SPM)
                particle = ParticulatesSPM("MP1-BF-SPM", MP1_BF, SPM1) 
                particle.calc_volume(MP1_BF, SPM1)
        
            #Create a dictionary with the objects 
            #particles_dict = {"A": MP1, "B": MP1_SPM, "C": MP1_BF, "D": MP1_BF_SPM }
                       
        
            RC_df.loc["degradation", spc[2:6]]= degradation(process_df.t_half_d.loc[idx]) 
            RC_df.loc["fragmentation", spc[2:6]]= fragmentation(process_df,idx,particle,sizeBinIdx,aggState)
            RC_df.loc["heteroagg", spc[2:6]]= heteroagg(process_df,idx,particle,SPM1, compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment,aggState)
            RC_df.loc["breakup", spc[2:6]]= breakup(process_df, idx, particle,SPM1,compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment, aggState)
            RC_df.loc["settling", spc[2:6]]= settling(particle, surface.depth_m, "Stokes", compartment)
            RC_df.loc["rising", spc[2:6]]= rising(particle, flowingWater.depth_m, "Stokes", compartment)
            RC_df.loc["advection", spc[2:6]]= advection(compartments_prop,comp_dict, compartment,riverSection,river_flows)
            RC_df.loc["mixing", spc[2:6]]= mixing(flowingWater,compartment, "up", comp_dict) #update creteria for mixing direction
            RC_df.loc["biofilm", spc[2:6]]= biofilm(compartment, process_df, comp_dict, idx, aggState)#add t biofilm growth?? Look at new processess MPLake-Antonia
            RC_df.loc["resusp", spc[2:6]]= resusp(compartment, comp_dict)
            RC_df.loc["burial", spc[2:6]]= burial(compartment, comp_dict)
            RC_df.loc["sedTransport", spc[2:6]]= sedTransport(compartment,comp_dict)
            RC_df.loc["defouling", spc[2:6]]= defouling(compartment, process_df, comp_dict, idx, aggState)
            RC_df.loc["volume_m3", spc[2:6]]= particle.volume_m3
            RC_df.loc["density_kg_m3", spc[2:6]]= particle.density_kg_m3
    
    elif numberRS <= 100:
        print("10 or more RS but less than 100")
        for spc in Clist:
            if spc in Clist[0:10*len(riverComp)*len(MPforms)*len(sizeBin)]:
                riverSection = spc[2]
                compartment= spc[3]
                aggState = spc[4]
                sizeBinIdx = spc[5]
                #spcAll=spc[2:6]
            elif spc in Clist[10*len(riverComp)*len(MPforms)*len(sizeBin):100*len(riverComp)*len(MPforms)*len(sizeBin)]:
                riverSection = spc[2:4]
                compartment= spc[4]
                aggState = spc[5]
                sizeBinIdx = spc[6]
                #spcAll=spc[2:7]
            #Generate particle OBJECTS MPs and SPM for the corresponding size bin
            
            #generate location index from imput tables
            MP_index = np.where((MP_prop['composition']== composition) & (MP_prop['sizeBin']== sizeBinIdx))[0][0]
            SPM_index = np.where((MP_prop['name']== "SPM1") & (MP_prop['sizeBin']== sizeBinIdx))[0][0]
            comp_index = np.where((compartments_prop['compartment']== int(compartment)) & (compartments_prop['riverSection']== int(riverSection)))[0][0]
            rs_indexs=np.where(compartments_prop['riverSection']== int(riverSection))[0]
            #for processes
            idx = np.where((process_df['aggState']==aggState) & (process_df['compartment']== int(compartment)))[0][0]
              
            #Generate compartment objects corresponding to the River Section
            [surface, flowingWater, stagnantWater, sediment]= preProcessLayers(mode, compartments_prop, date, rs_indexs,riverSection)
            #Create a dictionary with the environmental compartments objects 
            comp_dict = {"1": surface, "2": flowingWater, "3": stagnantWater, "4": sediment }

            #Generate particle object:

            #generate SPM object(s)
            #SPM_index = 11 #need to move SPM in own input file
            SPM1 = Particulates(MP_prop, SPM_index)
            SPM1.calc_volume()
            SPM1.calc_numConc(compartments_prop.SPM_mgL[comp_index], 0)

            #Generate particle
            if aggState== "A":
                particle= Particulates(MP_prop, MP_index)
                particle.calc_volume()
            elif aggState == "B":
                MP1= Particulates(MP_prop, MP_index)
                MP1.calc_volume()
                #B: heteroaggregated (MP attached to suspended particulate matter (SPM)) 
                particle = ParticulatesSPM("MP1-SPM", MP1, SPM1) 
                particle.calc_volume(MP1, SPM1)
            
            elif aggState == "C":
                MP1= Particulates(MP_prop, MP_index)
                MP1.calc_volume()
                #C: biofiolm-covered (MP with biofilm (BF) layer on surface)
                particle = ParticulatesBF("MP1-BF", MP1, 1388, 5e-6) 
                particle.calc_volume()
            
            elif aggState == "D":
                MP1= Particulates(MP_prop, MP_index)
                MP1.calc_volume()
                #C: biofiolm-covered (MP with biofilm (BF) layer on surface)
                MP1_BF = ParticulatesBF("MP1-BF", MP1, 1388, 5e-6) 
                MP1_BF.calc_volume()
                #D: biofilm-heteroaggregated (MP with BF layer attached to SPM)
                particle = ParticulatesSPM("MP1-BF-SPM", MP1_BF, SPM1) 
                particle.calc_volume(MP1_BF, SPM1)
        
        
            RC_df.loc["degradation", spc[2:]]= degradation(process_df.t_half_d.loc[idx]) 
            RC_df.loc["fragmentation", spc[2:]]= fragmentation(process_df,idx,particle,sizeBinIdx,aggState)
            RC_df.loc["heteroagg", spc[2:]]= heteroagg(process_df,idx,particle,SPM1,compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment,aggState)
            RC_df.loc["breakup", spc[2:]]= breakup(process_df,idx,particle,SPM1,compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment,aggState)
            RC_df.loc["settling", spc[2:]]= settling(particle, surface.depth_m, "Stokes", compartment)
            RC_df.loc["rising", spc[2:]]= rising(particle, flowingWater.depth_m, "Stokes", compartment)
            RC_df.loc["advection", spc[2:]]= advection(compartments_prop,comp_dict, compartment,riverSection,river_flows)
            RC_df.loc["mixing", spc[2:]]= mixing(flowingWater,compartment, "up", comp_dict) #update creteria for mixing direction
            RC_df.loc["biofilm", spc[2:]]= biofilm(compartment, process_df, comp_dict, idx, aggState)#add t biofilm growth?? Look at new processess MPLake-Antonia
            RC_df.loc["resusp", spc[2:]]= resusp(compartment, comp_dict)
            RC_df.loc["burial", spc[2:]]= burial(compartment, comp_dict)
            RC_df.loc["sedTransport", spc[2:]]= sedTransport(compartment,comp_dict)
            RC_df.loc["defouling", spc[2:]]= defouling(compartment, process_df, comp_dict, idx, aggState)
            RC_df.loc["volume_m3", spc[2:]]= particle.volume_m3
            RC_df.loc["density_kg_m3", spc[2:]]= particle.density_kg_m3
            
    elif numberRS <= 1000:
        print("100 or more RS but less than 1000")
        for spc in Clist:
            if spc in Clist[0:10*len(riverComp)*len(MPforms)*len(sizeBin)]:
                riverSection = spc[2]
                compartment= spc[3]
                aggState = spc[4]
                sizeBinIdx = spc[5]
            elif spc in Clist[10*len(riverComp)*len(MPforms)*len(sizeBin):100*len(riverComp)*len(MPforms)*len(sizeBin)]:
                riverSection = spc[2:4]
                compartment= spc[4]
                aggState = spc[5]
                sizeBinIdx = spc[6]
            elif spc in Clist[100*len(riverComp)*len(MPforms)*len(sizeBin):1000*len(riverComp)*len(MPforms)*len(sizeBin)]:
                riverSection = spc[2:5]
                compartment= spc[5]
                aggState = spc[6]
                sizeBinIdx = spc[7]
            #Generate particle OBJECTS MPs and SPM for the corresponding size bin
            
            #generate location index from imput tables
            MP_index = np.where((MP_prop['composition']== composition) & (MP_prop['sizeBin']== sizeBinIdx))[0][0]
            SPM_index = np.where((MP_prop['name']== "SPM1") & (MP_prop['sizeBin']== sizeBinIdx))[0][0]
            comp_index = np.where((compartments_prop['compartment']== int(compartment)) & (compartments_prop['riverSection']== int(riverSection)))[0][0]
            rs_indexs=np.where(compartments_prop['riverSection']== int(riverSection))[0]
            #for processes
            idx = np.where((process_df['aggState']==aggState) & (process_df['compartment']== int(compartment)))[0][0]
                         
            #Generate compartment objects corresponding to the River Section
            [surface, flowingWater, stagnantWater, sediment]= preProcessLayers(mode, compartments_prop, date, rs_indexs,riverSection)
            #Create a dictionary with the environmental compartments objects 
            comp_dict = {"1": surface, "2": flowingWater, "3": stagnantWater, "4": sediment }
            
            #Generate particle object:

            #generate SPM object(s)
            #SPM_index = 11 #need to move SPM in own input file
            SPM1 = Particulates(MP_prop, SPM_index)
            SPM1.calc_volume()
            SPM1.calc_numConc(compartments_prop.SPM_mgL[comp_index], 0)

            #Generate particle
            if aggState== "A":
                particle= Particulates(MP_prop, MP_index)
                particle.calc_volume()
            elif aggState == "B":
                MP1= Particulates(MP_prop, MP_index)
                MP1.calc_volume()
                #B: heteroaggregated (MP attached to suspended particulate matter (SPM)) 
                particle = ParticulatesSPM("MP1-SPM", MP1, SPM1) 
                particle.calc_volume(MP1, SPM1)
            
            elif aggState == "C":
                MP1= Particulates(MP_prop, MP_index)
                MP1.calc_volume()
                #C: biofiolm-covered (MP with biofilm (BF) layer on surface)
                particle = ParticulatesBF("MP1-BF", MP1, 1388, 5e-6) 
                particle.calc_volume()
            
            elif aggState == "D":
                MP1= Particulates(MP_prop, MP_index)
                MP1.calc_volume()
                #C: biofiolm-covered (MP with biofilm (BF) layer on surface)
                MP1_BF = ParticulatesBF("MP1-BF", MP1, 1388, 5e-6) 
                MP1_BF.calc_volume()
                #D: biofilm-heteroaggregated (MP with BF layer attached to SPM)
                particle = ParticulatesSPM("MP1-BF-SPM", MP1_BF, SPM1) 
                particle.calc_volume(MP1_BF, SPM1)
                   
            RC_df.loc["degradation", spc[2:]]= degradation(process_df.t_half_d.loc[idx]) 
            RC_df.loc["fragmentation", spc[2:]]= fragmentation(process_df,idx,particle,sizeBinIdx,aggState)
            RC_df.loc["heteroagg", spc[2:]]= heteroagg(process_df,idx,particle,SPM1,compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment,aggState)
            RC_df.loc["breakup", spc[2:]]= breakup(process_df,idx,particle,SPM1,compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment,aggState)
            RC_df.loc["settling", spc[2:]]= settling(particle, surface.depth_m, "Stokes", compartment)
            RC_df.loc["rising", spc[2:]]= rising(particle, surface.depth_m, "Stokes", compartment)
            RC_df.loc["advection", spc[2:]]=  advection(compartments_prop,comp_dict, compartment,riverSection,river_flows)
            RC_df.loc["mixing", spc[2:]]= mixing(flowingWater,compartment, "up", comp_dict) #update creteria for mixing direction
            RC_df.loc["biofilm", spc[2:]]= biofilm(compartment, process_df, comp_dict, idx, aggState)#add t biofilm growth?? Look at new processess MPLake-Antonia
            RC_df.loc["resusp", spc[2:]]= resusp(compartment, comp_dict)
            RC_df.loc["burial", spc[2:]]= burial(compartment, comp_dict)
            RC_df.loc["sedTransport", spc[2:]]= sedTransport(compartment,comp_dict)
            RC_df.loc["defouling", spc[2:]]= defouling(compartment, process_df, comp_dict, idx, aggState)
            RC_df.loc["volume_m3", spc[2:]]= particle.volume_m3
            RC_df.loc["density_kg_m3", spc[2:]]= particle.density_kg_m3
    else:
        print("Number of river sections excedes 1000. Implement new code for RCs")    
    #substitute NAN values by 0 in Rate cosntants dataframe 
    for x in RC_df:
        RC_df[x] = RC_df[x].fillna(0)
        

    
    #SPEED PROCESS CONTROL
    #To speed up the solver eliminate those RC for wich proccesses are faster than 30min (1800sec)--> RC<5.56x10-4
    if mode2 == "Timelimit":
        RC_df_templim=pd.DataFrame(index=processList , columns= CombList)
        
        for x in RC_df:
            for y in range(len(RC_df)-2): 
                if type(RC_df[x][y]) is tuple:
                    # when there are fragements formed extract the fragmentation array to multiply kfrag by number of fragments formed
                    frag = RC_df[x][y]
                    fragRC = frag[0]
                    if (fragRC !=0) and (fragRC > 0.000556):
                        #print(RC_df.index[y]+" RC for " +x+ " is fastest than 30min")
                        fragRC=0.000556
                        RC_df_templim[x][y]=(fragRC,frag[1])
                    else:
                        RC_df_templim[x][y]=RC_df[x][y]
                elif type(RC_df[x][y]) == list:
                    lis=RC_df[x][y]
                    RC_df_templim[x][y]=[0.000556 if i> 0.000556 else i for i in lis]
                else:
                    if RC_df[x][y] > 0.000556:
                        RC_df_templim[x][y]= 0.000556
                    else:
                        RC_df_templim[x][y]=RC_df[x][y]
                    
        for x in RC_df_templim:
            if (RC_df_templim.loc["breakup",[x]][0] == 0.000556):
                RC_df_templim.loc["breakup",[x]] = 0.000556/10
        
        #Add a row to store values of volume of the particles and its density
        extraRow=RC_df.drop(processList)
        RC_df=pd.concat([RC_df_templim,extraRow])
    return RC_df

