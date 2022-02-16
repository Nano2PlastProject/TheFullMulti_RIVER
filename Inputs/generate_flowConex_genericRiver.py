# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 16:15:00 2022

@author: PradoDomercq
"""

#Script to generate the flow_connectivity.csv imput file based on the desired 
#structure of the river

#Create a container (pandas data frame) to populate with data and later convert to csv

import pandas as pd

#Generic River parameters for a representative unit cell (rier section):
v_flow_m_s=1.3
CrossArea_m2=5.42*80

#Create a container (pandas data frame)
flow_connectivity=pd.DataFrame(columns=["Region_I","Region_J","q(m3/h)"])

#Populate dataframe
flow_connectivity["Region_I"]=list(range(1,21))
flow_connectivity["Region_J"]=list(range(2,21))+["nan"]
flow_connectivity["q(m3/h)"]=v_flow_m_s*CrossArea_m2*60*60

title_csv="flow_connectivity.csv"
flow_connectivity.to_csv(title_csv)
