# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 16:10:20 2019

@author: AntoniaPraetorius

"""
import math
from GlobalConstants import *

from Particulates import Particulates #class to generate MP and SPM objects


#define ParticulatesSPM class (inheriting from Particulates class)
class ParticulatesSPM(Particulates):
    "This is a class to create ParticulatesSPM objects" 
    
    #class attribute
    species = "particulate"
    
    #constructor
    def __init__(self, name, parentMP, parentSPM):
        self.name = name
        self.parentMP = parentMP
        self.parentSPM = parentSPM
        self.density_kg_m3 = parentMP.density_kg_m3*(parentMP.volume_m3/(parentMP.volume_m3+parentSPM.volume_m3))+parentSPM.density_kg_m3*(parentSPM.volume_m3/(parentMP.volume_m3+parentSPM.volume_m3))        
        self.radius_m = (3*(parentMP.volume_m3+parentSPM.volume_m3)/(4*math.pi))**(1/3)  #Note: this is an equivalent radius. MP-SPM most likely not truly spherical
        self.diameter_m = self.radius_m*2
        self.diameter_um = self.diameter_m*1e6
        
        self.composition = parentMP.composition + "-SPM"
        self.shape = parentMP.shape #to be updated for biofilm, could argue that shape is retained (unlike for SPM-bound)
                    

    
    #methods
    
    #volume calculation - currently simple version. 
    #more complexity to be added later:
    #different formulas for different particle shapes.
    #currently defined for spheres, fibres, cylinders, pellets and irregular fragments
    def calc_volume(self, parentMP, parentSPM):
        
        self.volume_m3 = parentMP.volume_m3 + parentSPM.volume_m3

    def calc_settling(self):
        
        self.vSet_m_s = 2/9*(self.density_kg_m3-density_w_21C_kg_m3)/mu_w_21C_kg_ms*g_m_s2*(self.radius_m)**2
        
        
        
    def calc_numConc(self, concMass_mg_L, concNum_part_L):
        
        if concNum_part_L == 0:
            self.concNum_part_m3 = concMass_mg_L/1000/self.density_kg_m3/self.volume_m3
            #if mass concentration is given, it is converted to number concentration
        else:
            self.concNum_part_m3 = concNum_part_L*1000
            #if number concentration is given, it is converted from part/L to part/m3
    
#    #settling calculations
#    """settling can be calculated using different equations (e.g. Stokes, 
#    modified versions of it or others) or be taken from experimental studies
#    !! currently only classical Stokes is implemented (which is probably not 
#    very realistic and will be updated soon !!""" 
#    def calc_settling(self, density_w_21C_kg_m3, mu_w_21C_kg_m3, g_m_s2, settlingMethod):
#        
#        if settlingMethod == "Stokes":
#            self.vSet_Stokes_m_s = 2/9*(self.density_kg_m3-density_w_21C_kg_m3)/mu_w_21C_kg_ms*g_m_s2*(self.radius_m)**2
#        else: 
#            print("Error: cannot calculate settling other than Stokes yet")
#            #print error message settling methods other than Stokes 
#            #(to be removed when other settling calculations are implemented)

    #degradation estimations
    """ relates only to MP & NPs. Full degradation probably extremely slow
    possibly not significant for most simulations. But add anyway for scenario
    analysis or biodegradable polymers. Values currently placeholders
    ! Add a size relation?!"""
    
    
#    def est_degradation(self, t_half_w_1_d, t_half_UV_1_d, t_half_sed_1_d):
#        #degradation half-lifes used as input are in days
#        
#        #degradation rate constant in water
#        self.k_deg_w = math.log(2)/(t_half_w_1_d*24*60*60) 
#        
#        #degradation rate constant by UV radation (relevant for water surface)
#        self.k_deg_UV = math.log(2)/(t_half_UV_1_d*24*60*60) 
#        
#        #degradation rate constant in sediment
#        self.k_deg_sed = math.log(2)/(t_half_sed_1_d*24*60*60)
        

#    #fragmentation estimations
#    """ relates only to MP & NPs. Data scarce. Need to investiage more to get 
#    numbers. Current values/pathways are placeholders for later. Currently 
#    assume fragmentation only in water
#    current function works only for 5 size bin and fragmentation from one
#    size to the next smaller size bin. Need to calculate how many fragments are 
#    created from next bigger size and estimate time that takes
#    """
#    
#    def est_fragmentation(self, k_frag_gen_w1_d, k_frag_gen_ws_d):
#        #inputs k_frag_gen_w1_d & k_frag_gen_ws_d are generic fragmentation rate
#        #constants describing how many MP fragment per day. To be scaled with a
#        #size factor to create size dependence
#        #fragmentation in water surface (ws) assumed to be twice as fast as in 
#        #main water (w1)
#        
#        #estimate fragmentation relation between size bins
#        volume_fragment = 4/3*math.pi*(self.radius_m/10)**3 #!!!only works for bins 10 times smaller!!!
#        self.fragments_formed = self.volume_m3/volume_fragment
#        self.k_frag_w1 = k_frag_gen_w1_d*self.diameter_um/1000/24/60/60
#        self.k_frag_ws = k_frag_gen_ws_d*self.diameter_um/1000/24/60/60
#
    