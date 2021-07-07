# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 16:37:08 2019

@author: Antonia Praetorius
"""

#this files stores all constants needed for running the model (i.e. for 
#particulate & environmental compartment objects & calcualting rate processes)

k_B_J_K = 1.38*10**-23 #Boltzmann constant k_B (in J/K)
g_m_s2 = 9.81 #gravitational acceleration on earth (in m/s2)

density_w_21C_kg_m3 = 998 #density of water at 21 degree C (in kg?m2)

mu_w_21C_mPas = 0.9764 #dynamic viscosit of water at 21 degree C (in mPa*s)
mu_w_21C_kg_ms = mu_w_21C_mPas/1000 #dynamic viscosit of water at 21 degree C (in kg/(m*s))

### Sediment parameters
v_sed_trans=3.0 #in kgs-1 as per Praetorius et al. 2012
sed_porosity=0.85 #as per Praetorius et al. 2012
sed_density=2.5#in g_cm3 as per Praetorius et al. 2012

