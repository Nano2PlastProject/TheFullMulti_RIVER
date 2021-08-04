# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 16:10:20 2019

@author: AntoniaPraetorius

"""
import math
from GlobalConstants import *

from Particulates import Particulates #class to generate MP and SPM objects


#define ParticulatesBF class (inheriting from Particulates class)
class ParticulatesBF(Particulates):
    "This is a class to create ParticulatesBIOFILM objects" 
    
    #class attribute
    species = "particulate"
    
    #constructor
    def __init__(self, name, parentMP, BF_density_kg_m3, BF_thickness_m):
        self.name = name
        self.parentMP = parentMP
        self.BF_density_kg_m3 = BF_density_kg_m3
        self.BF_thickness_m = BF_thickness_m
        self.radius_m = parentMP.radius_m + BF_thickness_m
        self.diameter_m = self.radius_m*2
        self.diameter_um = self.diameter_m*1e6
        
        if parentMP.length_a_m == 0:
            self.length_a_m = 0
        else:
            self.length_a_m = parentMP.length_a_m + BF_thickness_m*2
            
        if parentMP.length_b_m == 0:
            self.length_b_m = 0
        else:  
            self.length_b_m = parentMP.length_b_m + BF_thickness_m*2
            
        if parentMP.length_c_m == 0:
            self.length_c_m = 0
        else:
            self.length_c_m = parentMP.length_c_m + BF_thickness_m*2
            
        self.composition = parentMP.composition + "-BF"
        self.shape = parentMP.shape #to be updated for biofilm, could argue that shape is retained (unlike for SPM-bound)
        
        self.density_kg_m3 =  (self.parentMP.radius_m**3*self.parentMP.density_kg_m3 + ((self.parentMP.radius_m + self.BF_thickness_m)**3 - self.parentMP.radius_m**3)*self.BF_density_kg_m3)/((self.parentMP.radius_m + self.BF_thickness_m)**3)
        #equation from Kooi et al for density

      # default values
      #BF_thickness_m = 5*10**-6 #assuming a biofilm thickness of 5 um (very roughly estimated from algae volum in Kooi et al, to be refined!!)
      #BF_density_kg_m3 = 1388 #from Kooi et al
            
    
  
    def calc_settling(self):
        
        self.vSet_m_s = 2/9*(self.density_kg_m3-density_w_21C_kg_m3)/mu_w_21C_kg_ms*g_m_s2*(self.radius_m)**2
        
        
        
    def calc_numConc(self, concMass_mg_L, concNum_part_L):
        
        if concNum_part_L == 0:
            self.concNum_part_m3 = concMass_mg_L/1000/self.density_kg_m3/self.volume_m3
            #if mass concentration is given, it is converted to number concentration
        else:
            self.concNum_part_m3 = concNum_part_L*1000
            #if number concentration is given, it is converted from part/L to part/m3
    

    