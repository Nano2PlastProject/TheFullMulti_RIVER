# -*- coding: utf-8 -*-
"""
Created on Thu May 14 15:58:47 2020

@author: A. Praetorius and PradoDomercq
"""


#extension of RS_generator module containing functions to calculate all rate constants 
#Modification of advection and addition of sed transport for rivers

import math

#import file storing required constants
from GlobalConstants import *

    
    
    
def degradation(t_half_d):
    
    #degradation estimations
    """ relates only to MP & NPs. Full degradation probably extremely slow
    possibly not significant for most simulations. But add anyway for scenario
    analysis or biodegradable polymers. Values currently placeholders
    ! Add a size relation?!"""
    #degradation half-life of MPs used as input is in days
        
    #degradation rate constant 
    k_deg = math.log(2)/(t_half_d*24*60*60) 

    return k_deg

    
def fragmentation(t_frag_d, MP_radius_m, MP_volume_m3, MP_diameter_um, sizeBinIdx):   
    #estimate fragmentation relation between size bins (all except smallest size bin)
    
    #since sizeBins should be always organized froom smallest to largest
    
    if sizeBinIdx == "a":
        if t_frag_d == "NAN":
            k_frag = 0
        else:
            #print("Smallest sizeBin, fragments formed will be considered losses")
            k_frag = (1/(float(t_frag_d)*24*60*60))*MP_diameter_um/1000
        
        return (k_frag)
    else:
        if t_frag_d == "NAN":
            k_frag = 0
            fragments_formed=0
        else:
            volume_fragment = 4/3*math.pi*(MP_radius_m/10)**3 #!!!only works for bins 10 times smaller!!!
            fragments_formed = MP_volume_m3/volume_fragment
            k_frag = (1/(float(t_frag_d)*24*60*60))*MP_diameter_um/1000
        
        return (k_frag, fragments_formed)
  
       
    #NOTE: check this again!!
    
    


def settling(MP_density_kg_m3, MP_radius_m, comp_depth_m, settlingMethod, compartment):
    #settling calculations
    """settling can be calculated using different equations (e.g. Stokes, 
    modified versions of it or others) or be taken from experimental studies
    !! currently only classical Stokes is implemented (which is probably not 
    very realistic and will be updated soon !!""" 
    
    #Settling occurs in all compartments but in sediment (compartment 4)
    if compartment == "4":
        k_set = 0
        
    else:
        if settlingMethod == "Stokes":
            vSet_m_s = 2/9*(MP_density_kg_m3-density_w_21C_kg_m3)/mu_w_21C_kg_ms*g_m_s2*(MP_radius_m)**2
        else: 
            print("Error: cannot calculate settling other than Stokes yet")
            #print error message settling methods other than Stokes 
            #(to be removed when other settling calculations are implemented)
            
        #for the water and surface water compartments:
        #settling and rising rate constants for free MP
        if vSet_m_s > 0:
            k_set = vSet_m_s/comp_depth_m
                
        elif vSet_m_s  < 0:
            k_set  = 0
            
        else:
            k_set  = 0
            
    return k_set
        

def rising(MP_density_kg_m3, MP_radius_m, comp_depth_m, settlingMethod, compartment):
    #rising calculations
    """rising is calculated in the same way as settling for particles with negative 
    settling velocitis. 
    It can be calculated using different equations (e.g. Stokes, 
    modified versions of it or others) or be taken from experimental studies
    !! currently only classical Stokes is implemented (which is probably not 
    very realistic and will be updated soon !!""" 
    
    #Rising only occus in the flowing water and stagnant water compartments (2 and 3)
    
    if compartment == "2" or compartment == "3": 
 
        if settlingMethod == "Stokes":
            vSet_m_s = 2/9*(MP_density_kg_m3-density_w_21C_kg_m3)/mu_w_21C_kg_ms*g_m_s2*(MP_radius_m)**2
        else: 
            print("Error: cannot calculate settling other than Stokes yet")
        #print error message settling methods other than Stokes 
        #(to be removed when other settling calculations are implemented)
    else:
         vSet_m_s =0
    #for the water and surface water compartments:
    #settling and rising rate constants for free MP
    if vSet_m_s > 0:
        k_rise = 0
            
    elif vSet_m_s  < 0:
        k_rise = -vSet_m_s/comp_depth_m
        
    else:
        k_rise  = 0
        
    return k_rise

def heteroagg(alpha, MP_radius_m, SPM_radius_m, MP_density_kg_m3, SPM_density_kg_m3, SPM_concNum_part_m3, G, T_K, compartment):

    #Heteroaggregation only occurs for pristine (A) and biofouled MPs (C) and on the water compartments (1, 2 and 3)
    #heteroaggregation for B and D is limited by alpha values given as NA
    if compartment == "4":
        k_hetAgg = 0
    else:
    
        #heteroaggregation rate constants
        """heteroaggregation requires to particles to collide and interact
        favorably for the collision to result in attachment
        the heteroaggregation rate constants is therefore composed of two
        parts, 1) a collision rate constant and 2) and attachement 
        efficiency (alpha) (representing the probability of attachement).
        For heteroaggregation a common simplifaction is the assumption that
        SPM concentration is not signficantly affected by the heteroaggre-
        gation process. Therefore, a pseudo first-order heteroaggregation 
        rate constant is obtained by multiplying collision rate with alpha 
        and with the SPM number concentration"""
        
        #first the different collision mechanisms are calculated
        k_peri = (2*k_B_J_K*T_K)/(3*mu_w_21C_kg_ms)*(MP_radius_m + SPM_radius_m)**2/(MP_radius_m * SPM_radius_m)
        #perikinetic contributions to collision rate constant (Brownian motion)
        
        k_ortho = 4/3*G*(MP_radius_m + SPM_radius_m)**3
        #orthokinetic contributions to collision rate constant (caused by fluid motion)
        
        MP_vSet_m_s = 2/9*(MP_density_kg_m3-density_w_21C_kg_m3)/mu_w_21C_kg_ms*g_m_s2*(MP_radius_m)**2
        SPM_vSet_m_s = 2/9*(SPM_density_kg_m3-density_w_21C_kg_m3)/mu_w_21C_kg_ms*g_m_s2*(SPM_radius_m)**2
        #settling velocity. currently according to classical Stokes law. Need to include other modes and put calculation on its own, so that it can also be accessed for other processes
        
        k_diffSettling = math.pi*(MP_radius_m + SPM_radius_m)**2 * abs(MP_vSet_m_s-SPM_vSet_m_s)
        #differential settling contributions to collision rate constant
    
        k_coll = k_peri + k_ortho + k_diffSettling
        #the collision rate constant
        
        k_hetAgg = alpha*k_coll*SPM_concNum_part_m3
        #the pseudo first-order heteroaggregation rate constant
        
    
    return k_hetAgg


def breakup(process_df, idx, MP_radius_m, SPM_radius_m, MP_density_kg_m3, SPM_density_kg_m3, SPM_concNum_part_m3, G, T_K, compartment, aggState):
    # Breackup doesnt occur in the sediment compartment and only for MP aggregates (B and D), 
    #however Kbreackup is calculated based on Kheter of A and C
    if compartment == "4" or aggState == "A" or aggState =="C":
        k_aggBreakup = 0
    else:
    
            #data is limited on aggregate breakup, but this process is likely
            #more relvant for larger aggregates
            #!! 1/10 of k_hetAgg is just a placeholder,  needs to be refined
            #possibly using a size dependent function !!
    
        #first the different collision mechanisms are calculated
        
        k_peri = (2*k_B_J_K*T_K)/(3*mu_w_21C_kg_ms)*(MP_radius_m + SPM_radius_m)**2/(MP_radius_m * SPM_radius_m)
        #perikinetic contributions to collision rate constant (Brownian motion)
        
        k_ortho = 4/3*G*(MP_radius_m + SPM_radius_m)**3
        #orthokinetic contributions to collision rate constant (caused by fluid motion)
        
        MP_vSet_m_s = 2/9*(MP_density_kg_m3-density_w_21C_kg_m3)/mu_w_21C_kg_ms*g_m_s2*(MP_radius_m)**2
        SPM_vSet_m_s = 2/9*(SPM_density_kg_m3-density_w_21C_kg_m3)/mu_w_21C_kg_ms*g_m_s2*(SPM_radius_m)**2
        #settling velocity. currently according to classical Stokes law. Need to include other modes and put calculation on its own, so that it can also be accessed for other processes
        
        k_diffSettling = math.pi*(MP_radius_m + SPM_radius_m)**2 * abs(MP_vSet_m_s-SPM_vSet_m_s)
        #differential settling contributions to collision rate constant
    
        k_coll = k_peri + k_ortho + k_diffSettling
        #the collision rate constant
        
        k_hetAgg = process_df.alpha.loc[idx-1]*k_coll*SPM_concNum_part_m3
        #the pseudo first-order heteroaggregation rate constant
    
        k_aggBreakup = 1/10*k_hetAgg
    
    return k_aggBreakup


def advection(comp_dict, compartment):
     #advective transport
    
    # Based on Praetorius et al. 2012: Kflow = v_riv_flow*(Aw1/Vw1)
    #Being v_riv_flow the river flow velocity in ms-1, Aw1 is the crossectional 
    #area of the flowing water and Vw1 the volume of the box of moving water.
    #dimensions of the river we estimated resudence time of 26 days in flowing
    #water and 28 min in the surface watercompartment
    
    #RIVER SECTION DEPENDENT WITH VARYING DISCHARGE 
    
    if comp_dict[compartment].name == "flowingWater":
        #k_adv = comp_volume_m3/(29*24*60*60) #Kflow = v_riv_flow = 1.5 m3s-1
        #k_adv = comp_volume_m3/(35*60*60) #Kflow = v_riv_flow = 29.8 m3s-1
        #k_adv = discharge/ comp_dict[compartment].volumeBox_m3 #Kflow = v_riv_flow = 65.8 m3s-1
        k_adv = comp_dict[compartment].v_riv_flow*(comp_dict[compartment].CrossArea_m2/comp_dict[compartment].volume_m3)
    elif comp_dict[compartment].name == "surface":
        #k_adv = comp_volume_m3/(28*60)#Kflow = v_riv_flow = 1.5 m3s-1
        #k_adv = comp_volume_m3/(1.4*60)#Kflow = v_riv_flow = 29.8 m3s-1
        #k_adv = discharge/comp_dict[compartment].volumeBox_m3 #Kflow = v_riv_flow = 65.8 m3s-1
        k_adv = comp_dict[compartment].v_riv_flow*(comp_dict[compartment].CrossArea_m2/comp_dict[compartment].volume_m3)
    else:
        k_adv = 0
        
    return k_adv


def mixing(flowingWater, compartment, updown, comp_dict):
    
    if comp_dict[compartment].name == "flowingWater": #flowingWater= 2
        if updown == "up":
            k_mix = 10e-10 #random number, need to refine!!
        if updown == "down":
            k_mix = 10e-13 #random number, need to refine!!
        
    if comp_dict[compartment].name == "surface":
        k_mix = 10e-10*(flowingWater.volume_m3/comp_dict[compartment].volume_m3)
    
    if comp_dict[compartment].name == "stagnantWater":
        k_mix = 10e-13*(flowingWater.volume_m3/comp_dict[compartment].volume_m3)
          
    if comp_dict[compartment].name == "sediment":
        k_mix = 0
    
    return k_mix



def biofilm(compartment, process_df, comp_dict, idx, aggState):
    #Biofilm formation taken into account for pristin and heteroaggregated MPs (A and C)
    #only takes place in the water compartments ( 1, 2 and 3)
    
    if aggState == "A" or "C":
    
        if comp_dict[compartment].name == "surface":
                if process_df.t_biof_growth_d.loc[idx] == 0:
                    k_biof = 0
                else:
                    k_biof = 1/process_df.t_biof_growth_d.loc[idx]/24/60/60
                
        if comp_dict[compartment].name == "flowingWater":
                if process_df.t_biof_growth_d.loc[idx] == 0:
                    k_biof = 0
                else:
                    k_biof = 1/process_df.t_biof_growth_d.loc[idx]/24/60/60
            
        if comp_dict[compartment].name == "stagnantWater":
                if process_df.t_biof_growth_d.loc[idx] == 0:
                    k_biof = 0
                else:
                    k_biof = 1/process_df.t_biof_growth_d.loc[idx]/24/60/60
            
        if comp_dict[compartment].name == "sediment":
                k_biof = 0
    else:
        k_biof = 0 
            
        #assume it takes x days for biofilm coverage to grow
        #need to update!!    # k_biof = 1/30/24/60/60 #assume it takes 30 days for biofilm coverage to grow
    # #need to update!!
    
    return k_biof

def defouling(compartment, process_df, comp_dict, idx, aggState):
    #Defouling = degradation of Biofilm. for biofouled and heteroaggregated and biofouled particles (B and D)
    #only takes place in the water compartments ( 1, 2 and 3)
    
    if aggState == "B" or "D":
        if comp_dict[compartment].name == "sediment":
            k_defoul = 0
    
        else:
            if type(process_df.t_biof_degrad_d.loc[idx]) == str:
                k_defoul = 0
            else:
                k_defoul = 1/process_df.t_biof_degrad_d.loc[idx]/24/60/60

    else:
        k_defoul = 0 
            
        #assume it takes x days for biofilm coverage to be degraded
        
    return k_defoul

#for the sediment compartment rate constants for resuspension and
            #burial in deep sediment are calculated & degradation rate assigned

def resusp(compartment, comp_dict):
   if comp_dict[compartment].name == "sediment":
       k_resusp = 2.3*10**-7/comp_dict[compartment].depth_m
   else:
       k_resusp=0
    
   return k_resusp



def burial(compartment, comp_dict):
    if comp_dict[compartment].name == "sediment":
    
        k_burial = 5.6*10**-7/comp_dict[compartment].depth_m
    else:
        k_burial = 0
    
    return k_burial

def sedTransport(compartment,comp_dict):
    if comp_dict[compartment].name == "sediment":
        m_sed_kg=(1-sed_porosity)*sed_density*10**3*comp_dict[compartment].volume_m3
        k_sed_trans=v_sed_trans/m_sed_kg
    else:
        k_sed_trans=0
    return k_sed_trans
    

