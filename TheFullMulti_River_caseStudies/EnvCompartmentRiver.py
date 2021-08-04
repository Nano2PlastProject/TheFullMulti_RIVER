# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 16:13:51 2019
Modified by Prado Domercq for MPs River Model on Aug 2020

@author: AntoniaPraetorius and Prado Domercq
"""

#In order to generate higher spatial resolution each river section will
# be subdivided in boxes of equal length set as l_box_m im the MAIN script

import math
from GlobalConstants import *


#define environmental compartment class
class EnvCompartment:
    "This is a class to create environmental compartment objects of the river"
    
    #class attribute
    species = "compartment"
        #constructor
    def __init__(self, compartments_prop, comp_index):
        self.name = compartments_prop.compType.loc[comp_index]
        self.RSname = compartments_prop.nameRS.loc[comp_index]
        self.compType = compartments_prop.compType.loc[comp_index]
        self.density_kg_m3 = density_w_21C_kg_m3
        self.mu_w_kg_ms = mu_w_21C_kg_ms
        self.depth_m = compartments_prop.depth_m.loc[comp_index]
        self.length_m = compartments_prop.length_m.loc[comp_index]
        self.width_m = compartments_prop.width_m.loc[comp_index]
        self.surfArea_m2 = compartments_prop.width_m.loc[comp_index]*compartments_prop.length_m.loc[comp_index]
        self.CrossArea_m2 = compartments_prop.width_m.loc[comp_index]*compartments_prop.depth_m.loc[comp_index]
        self.G = compartments_prop.G.loc[comp_index]
        self.T_K = compartments_prop.T_K.loc[comp_index]
        #self.n_boxes = compartments_prop.length_m.loc[comp_index]/l_box_m
        #self.surfAreaBox_m2 = compartments_prop.width_m.loc[comp_index]*compartments_prop.length_m.loc[comp_index]
        if "discharge_m3s" in compartments_prop:
            self.v_riv_flow = compartments_prop.discharge_m3s.loc[comp_index]/(compartments_prop.depthRS_m.loc[comp_index]*compartments_prop.width_m.loc[comp_index])
        else:
            self.v_riv_flow = compartments_prop.vFlow_m_s.loc[comp_index] 
    #methods
    #dimension calculation --> calculates surface area and volume of the individual boxes
    def calc_dimensions(self):
        
        if self.surfArea_m2 == 0:
            self.surfArea_m2 = self.length_m*self.width_m 
        else:
            self.surfArea_m2 = self.surfArea_m2
        #calculates surface area (in m2)
        
        self.volume_m3 = self.depth_m*self.surfArea_m2
        #calculates volume (in m3)

      
