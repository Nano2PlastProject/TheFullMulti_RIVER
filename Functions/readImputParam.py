# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 12:42:25 2020

@author: PradoDomercq modifications on Antonia Praetorius code for Lake version
"""
import pandas as pd
import numpy as np
from pathlib import Path
import sys

import datetime

#Script for reading imput files
# fn = "process_paramRiver.txt" 

def sanityCheckData(dataFrame):
    #check if there is any negative value
    #check the dataframe type
    a=dataFrame.dtypes
    for i in range(len(a)):
    #include only the numerical fields
        if a[i]!="object": 
            #accepts only positive values
            if any(dataFrame[a.index[i]]<0):
                print("Error! Negative values in the column :" + a.index[i])
                sys.exit(1)  

def readProcessparam (fileToOpen):
    try:
        process_param = pd.read_csv(fileToOpen, comment='#')
    except IOError:
        print("File not accessible")
        sys.exit(1)
    
    sanityCheckData(process_param)   
    return process_param
    
def microplasticData(fileToOpen):
    try:
        plastics_prop = pd.read_csv(fileToOpen, comment='#')
    except IOError:
        print("File not accessible")
        sys.exit(1)
    
    sanityCheckData(plastics_prop)      #need to make this specific to process info    
    return plastics_prop      

def readCompartmentData (fileToOpen):
    try:
        compartments_prop = pd.read_csv(fileToOpen, comment='#')
    except IOError:
        print("File not accessible")
        sys.exit(1)  

    sanityCheckData(compartments_prop)          
    return compartments_prop   


def generateFinalDataFrame(dayFirst, dayLast, mode, data_folder):
    finalDataframe = pd.DataFrame()
    if mode=="Standard":
        lakeFile = data_folder / "lake.txt"#change for river
        finalDataframe = readCompartmentData(lakeFile)
    elif mode == "Monthly":
        rangeDates = dateRangeGenerator(dayFirst,dayLast,mode)
        for i in range(len(rangeDates)):
            rangeDates[i]=rangeDates[i]
            fileName="lake-"+str(rangeDates[i])+".txt"#Change for the river
            lakeFile = data_folder / fileName
            compartments_prop_month = readCompartmentData(lakeFile)
            compartments_prop_month["date"] = rangeDates[i]
            finalDataframe = finalDataframe.append(compartments_prop_month)
        
    return finalDataframe 

def dateRangeGenerator(dayFirst,dayLast,mode):
    try:
        if mode=="Monthly":
            #dayFirst must be < dayLast
            if dayFirst > dayLast :
               print("Wrong dates order")
               sys.exit(1)
               
            if dayFirst.month == dayLast.month and dayFirst.year == dayLast.year:
                dateRange= pd.Series(pd.date_range(dayFirst, dayFirst ,periods=1))
                dateRange[0] = dateRange[0].to_period('M')
                return dateRange
            else :
                dateRange = pd.Series(pd.date_range(dayFirst, dayLast ,freq='M'))
                for i in range(len(dateRange)):
                    dateRange[i]=dateRange[i].to_period('M')
                return dateRange        
    except:
          #trigger a generic error
          print("Error in date range generation")
          sys.exit(1)
      