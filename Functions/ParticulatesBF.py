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
            
    
    #methods
    
#    #volume calculation
#    #different formulas for different particle shapes.
#    #currently defined for spheres, fibres, cylinders, pellets and irregular fragments
#    def calc_volume(self):
#        
#        if self.shape == "sphere":
#            self.volume_m3 = 4/3*math.pi*(self.radius_m)**3
#            #calculates volume (in m3) of spherical particles from MP radius  
#            self.CSF = 1    
#            #calculate corey shape factor (CSF) 
#            #(Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)        
#            #particle number concentration calculation
#         
#        elif self.shape == "fibre" or self.shape == "fiber" or self.shape == "cylinder":
#            self.volume_m3 = math.pi*(self.radius_m)**2*self.length_a_m
#            #calculates volume (in m3) of fibres or cylinders from diameter and  
#            #length assuming cylindrical shape
#            self.CSF = self.radius_m/math.sqrt(self.length_a_m*self.radius_m)    
#            #calculate corey shape factor (CSF) 
#            #(Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)        
#            #particle number concentration calculation
#         
#            
#        elif self.shape == "pellet" or self.shape == "fragment":
#            self.volume_m3 = self.length_a_m*self.length_b_m*self.length_c_m
#            #approximate volume calculation for irregular fragments
#            #approximated as a cuboid using longest, intermediate and shortest length
#            #!! Note: not sure if pellets fits best here or rather as sphere/cylinder
#            #might adjust later!!
#            self.CSF = self.length_c_m/math.sqrt(self.length_a_m*self.length_b_m)    
#            #calculate corey shape factor (CSF) 
#            #(Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)        
#            #particle number concentration calculation
#            
#        else:
#            print("Error: unknown shape")
#            #print error message for shapes other than spheres 
#            #(to be removed when other volume calculations are implemented)

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
    