# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 10:44:30 2021

@author: PradoDomercq
"""

"""Estimate RATE CONSTANTS"""
##Estimate rate constants for all "species" = combination of RS*riverComp*MPforms*sizeBins (len (Clist))
import numpy as np
import pandas as pd
import RC_GeneratorRiver
from objectGenerationRiver_func import*

def RC_estimation_function(processList,CombList,Clist,MP_prop,compartments_prop,process_df,numberRS, composition,mode2, mode, date,riverComp,MPforms,sizeBin):
    
    RC_df = pd.DataFrame(index=processList , columns= CombList)
    extraRow = pd.DataFrame(index= ["volume_m3","density_kg_m3"] , columns= CombList)
    RC_df=pd.concat([RC_df,extraRow])
    
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
            
            
            [MP1,SPM1,MP1_SPM,MP1_BF,MP1_BF_SPM]=preProcessElements(MP_prop, MP_index, SPM_index, compartments_prop, comp_index)#add compartments_prop, comp_index to include corresponding SPM concentration
        
            #Create a dictionary with the objects 
            particles_dict = {"A": MP1, "B": MP1_SPM, "C": MP1_BF, "D": MP1_BF_SPM }
              
            #load compartment objects corresponding to the River Section
            
            [surface, flowingWater, stagnantWater, sediment]= preProcessLayers(mode, compartments_prop, date, comp_index)
            #Create a dictionary with the environmental compartments objects 
            comp_dict = {"1": surface, "2": flowingWater, "3": stagnantWater, "4": sediment }
            
            #import data needed for each spcn(t_half, radius...)
            #for processes
            idx = np.where((process_df['aggState']==aggState) & (process_df['compartment']== int(compartment)))[0][0]
            
        
            RC_df.loc["degradation", spc[2:6]]= RC_GeneratorRiver.degradation(process_df.t_half_d.loc[idx]) 
            RC_df.loc["fragmentation", spc[2:6]]= RC_GeneratorRiver.fragmentation(process_df.t_frag_d.loc[idx], particles_dict[aggState].radius_m,  particles_dict[aggState].volume_m3,  particles_dict[aggState].diameter_um, sizeBinIdx)
            RC_df.loc["heteroagg", spc[2:6]]= RC_GeneratorRiver.heteroagg(process_df.alpha.loc[idx], MP1.radius_m, SPM1.radius_m, MP1.density_kg_m3, SPM1.density_kg_m3, SPM1.concNum_part_m3, compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment)
            RC_df.loc["breakup", spc[2:6]]= RC_GeneratorRiver.breakup(process_df, idx, MP1.radius_m, SPM1.radius_m, MP1.density_kg_m3, SPM1.density_kg_m3, SPM1.concNum_part_m3,compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment, aggState)
            RC_df.loc["settling", spc[2:6]]= RC_GeneratorRiver.settling(particles_dict[aggState].density_kg_m3, particles_dict[aggState].radius_m, surface.depth_m, "Stokes", compartment)
            RC_df.loc["rising", spc[2:6]]= RC_GeneratorRiver.rising(particles_dict[aggState].density_kg_m3, particles_dict[aggState].radius_m, flowingWater.depth_m, "Stokes", compartment)
            RC_df.loc["advection", spc[2:6]]= RC_GeneratorRiver.advection(comp_dict, compartment)
            RC_df.loc["mixing", spc[2:6]]= RC_GeneratorRiver.mixing(flowingWater,compartment, "up", comp_dict) #update creteria for mixing direction
            RC_df.loc["biofilm", spc[2:6]]= RC_GeneratorRiver.biofilm(compartment, process_df, comp_dict, idx, aggState)#add t biofilm growth?? Look at new processess MPLake-Antonia
            RC_df.loc["resusp", spc[2:6]]= RC_GeneratorRiver.resusp(compartment, comp_dict)
            RC_df.loc["burial", spc[2:6]]= RC_GeneratorRiver.burial(compartment, comp_dict)
            RC_df.loc["sedTransport", spc[2:6]]= RC_GeneratorRiver.sedTransport(compartment,comp_dict)
            RC_df.loc["volume_m3", spc[2:6]]= particles_dict[aggState].volume_m3
            RC_df.loc["density_kg_m3", spc[2:6]]= particles_dict[aggState].density_kg_m3
    
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
            #else:
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
            
            
            [MP1,SPM1,MP1_SPM,MP1_BF,MP1_BF_SPM]=preProcessElements(MP_prop, MP_index, SPM_index, compartments_prop, comp_index)#add compartments_prop, comp_index to include corresponding SPM concentration
        
            #Create a dictionary with the objects 
            particles_dict = {"A": MP1, "B": MP1_SPM, "C": MP1_BF, "D": MP1_BF_SPM }
              
            #load compartment objects corresponding to the River Section
            
            [surface, flowingWater, stagnantWater, sediment]= preProcessLayers(mode, compartments_prop, date, comp_index)
            #Create a dictionary with the environmental compartments objects 
            comp_dict = {"1": surface, "2": flowingWater, "3": stagnantWater, "4": sediment }
            
            #import data needed for each spcn(t_half, radius...)
            #for processes
            idx = np.where((process_df['aggState']==aggState) & (process_df['compartment']== int(compartment)))[0][0]
            
        
            RC_df.loc["degradation", spc[2:]]= RC_GeneratorRiver.degradation(process_df.t_half_d.loc[idx]) 
            RC_df.loc["fragmentation", spc[2:]]= RC_GeneratorRiver.fragmentation(process_df.t_frag_d.loc[idx], particles_dict[aggState].radius_m,  particles_dict[aggState].volume_m3,  particles_dict[aggState].diameter_um, sizeBinIdx)
            RC_df.loc["heteroagg", spc[2:]]= RC_GeneratorRiver.heteroagg(process_df.alpha.loc[idx], MP1.radius_m, SPM1.radius_m, MP1.density_kg_m3, SPM1.density_kg_m3, SPM1.concNum_part_m3, compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment)
            RC_df.loc["breakup", spc[2:]]= RC_GeneratorRiver.breakup(process_df, idx, MP1.radius_m, SPM1.radius_m, MP1.density_kg_m3, SPM1.density_kg_m3, SPM1.concNum_part_m3,compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment, aggState)
            RC_df.loc["settling", spc[2:]]= RC_GeneratorRiver.settling(particles_dict[aggState].density_kg_m3, particles_dict[aggState].radius_m, surface.depth_m, "Stokes", compartment)
            RC_df.loc["rising", spc[2:]]= RC_GeneratorRiver.rising(particles_dict[aggState].density_kg_m3, particles_dict[aggState].radius_m, flowingWater.depth_m, "Stokes", compartment)
            RC_df.loc["advection", spc[2:]]= RC_GeneratorRiver.advection(comp_dict, compartment)
            RC_df.loc["mixing", spc[2:]]= RC_GeneratorRiver.mixing(flowingWater,compartment, "up", comp_dict) #update creteria for mixing direction
            RC_df.loc["biofilm", spc[2:]]= RC_GeneratorRiver.biofilm(compartment, process_df, comp_dict, idx, aggState)#add t biofilm growth?? Look at new processess MPLake-Antonia
            RC_df.loc["resusp", spc[2:]]= RC_GeneratorRiver.resusp(compartment, comp_dict)
            RC_df.loc["burial", spc[2:]]= RC_GeneratorRiver.burial(compartment, comp_dict)
            RC_df.loc["sedTransport", spc[2:]]= RC_GeneratorRiver.sedTransport(compartment,comp_dict)
            RC_df.loc["volume_m3", spc[2:]]= particles_dict[aggState].volume_m3
            RC_df.loc["density_kg_m3", spc[2:]]= particles_dict[aggState].density_kg_m3
            
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
            
            
            [MP1,SPM1,MP1_SPM,MP1_BF,MP1_BF_SPM]=preProcessElements(MP_prop, MP_index, SPM_index, compartments_prop, comp_index)#add compartments_prop, comp_index to include corresponding SPM concentration
        
            #Create a dictionary with the objects 
            particles_dict = {"A": MP1, "B": MP1_SPM, "C": MP1_BF, "D": MP1_BF_SPM }
              
            #load compartment objects corresponding to the River Section
            
            [surface, flowingWater, stagnantWater, sediment]= preProcessLayers(mode, compartments_prop, date, comp_index)
            #Create a dictionary with the environmental compartments objects 
            comp_dict = {"1": surface, "2": flowingWater, "3": stagnantWater, "4": sediment }
            
            #import data needed for each spcn(t_half, radius...)
            #for processes
            idx = np.where((process_df['aggState']==aggState) & (process_df['compartment']== int(compartment)))[0][0]
            
        
            RC_df.loc["degradation", spc[2:]]= RC_GeneratorRiver.degradation(process_df.t_half_d.loc[idx]) 
            RC_df.loc["fragmentation", spc[2:]]= RC_GeneratorRiver.fragmentation(process_df.t_frag_d.loc[idx], particles_dict[aggState].radius_m,  particles_dict[aggState].volume_m3,  particles_dict[aggState].diameter_um, sizeBinIdx)
            RC_df.loc["heteroagg", spc[2:]]= RC_GeneratorRiver.heteroagg(process_df.alpha.loc[idx], MP1.radius_m, SPM1.radius_m, MP1.density_kg_m3, SPM1.density_kg_m3, SPM1.concNum_part_m3, compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment)
            RC_df.loc["breakup", spc[2:]]= RC_GeneratorRiver.breakup(process_df, idx, MP1.radius_m, SPM1.radius_m, MP1.density_kg_m3, SPM1.density_kg_m3, SPM1.concNum_part_m3,compartments_prop.G.loc[comp_index], compartments_prop.T_K.loc[comp_index], compartment, aggState)
            RC_df.loc["settling", spc[2:]]= RC_GeneratorRiver.settling(particles_dict[aggState].density_kg_m3, particles_dict[aggState].radius_m, surface.depth_m, "Stokes", compartment)
            RC_df.loc["rising", spc[2:]]= RC_GeneratorRiver.rising(particles_dict[aggState].density_kg_m3, particles_dict[aggState].radius_m, flowingWater.depth_m, "Stokes", compartment)
            RC_df.loc["advection", spc[2:]]= RC_GeneratorRiver.advection(comp_dict, compartment,compartments_prop.discharge_m3s[comp_index])
            RC_df.loc["mixing", spc[2:]]= RC_GeneratorRiver.mixing(flowingWater,compartment, "up", comp_dict) #update creteria for mixing direction
            RC_df.loc["biofilm", spc[2:]]= RC_GeneratorRiver.biofilm(compartment, process_df, comp_dict, idx, aggState)#add t biofilm growth?? Look at new processess MPLake-Antonia
            RC_df.loc["resusp", spc[2:]]= RC_GeneratorRiver.resusp(compartment, comp_dict)
            RC_df.loc["burial", spc[2:]]= RC_GeneratorRiver.burial(compartment, comp_dict)
            RC_df.loc["sedTransport", spc[2:]]= RC_GeneratorRiver.sedTransport(compartment,comp_dict)
            RC_df.loc["volume_m3", spc[2:]]= particles_dict[aggState].volume_m3
            RC_df.loc["density_kg_m3", spc[2:]]= particles_dict[aggState].density_kg_m3
    else:
        print("Number of river sections excedes 1000. Implement new code for RCs")    
    #substitute NAN values by 0 in Rate cosntants dataframe 
    for x in RC_df:
        RC_df[x] = RC_df[x].fillna(0)
    
    #SPEED PROCESS CONTROL
    #To speed up the solver eliminate those RC for wich proccesses are faster than 30min (1800sec)--> RC<5.56x10-4
    if mode2 == "Timelimit":
        RC_df_templim=pd.DataFrame(index=processList , columns= CombList)
        #Add a row to store values of volume of the particles and its density
        extraRow = pd.DataFrame(index= ["volume_m3","density_kg_m3"] , columns= CombList)
        RC_df_templim=pd.concat([RC_df_templim,extraRow])
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
                else:   
                    if (RC_df[x][y]!=0) and (RC_df[x][y] > 0.000556):
                        #print(RC_df.index[y]+" RC for " +x+ " is fastest than 30min")
                        RC_df_templim[x][y]=0.000556
                    else:
                        #print(RC_df[x][y])
                        RC_df_templim[x][y]=RC_df[x][y]
            for y in range(len(RC_df)-2,len(RC_df)):
                RC_df_templim[x][y]= RC_df[x][y]
        for x in RC_df_templim:
            if (RC_df_templim.loc["breakup",[x]][0] == 0.000556):
                RC_df_templim.loc["breakup",[x]] = 0.000556/10
        
        RC_df=RC_df_templim

    return RC_df