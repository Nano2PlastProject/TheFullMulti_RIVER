# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 14:50:56 2020

@author: PradoDomercq
"""

###TRANSPORT MATRIX 

import numpy as np
import pandas as pd
#Script to fillin the initeractions matrix (interactions_df)
  


def fillMatricesOfInteractions(RC_df, Clist,compartments_prop):
    
    
    #Remove volume and density from RC_df
    
    RC_df= RC_df.drop("volume_m3")
    RC_df= RC_df.drop("density_kg_m3")
     
    #remove C from Clist:
    SpeciesList=[]
    for a in Clist:
        SpeciesList.append(a[2:])
    
    #Diagonal of the df corresponds to the losses of each species
    #crate the array of values for the diagonal wich is the sum of all RC corresponding to one species:
    diag_list = []
    
    for sp in SpeciesList:
        if type(RC_df[sp]["fragmentation"]) is tuple:
            # when there are fragements formed extract the fragmentation array to multiply kfrag by number of fragments formed (fro the mas formed)
            #IN the lossess (estimation of mass lost) process we dont need to multiply by the number of fragments!!
            frag = RC_df[sp]["fragmentation"]
            #fragval = frag[0]*frag[1]
            
            RC_df_noFrag = RC_df[sp].drop("fragmentation")
            diag_list.append(-(sum(RC_df_noFrag)+frag[0]))
        
        else:
            diag_list.append(-(sum(RC_df[sp])))
            # diag_array = np.asarray(diag_list)
            # diag_array[np.isnan(diag_array)] = 0 #substitute nan by 0
    
    #Dataframes to indicate interactions between "Species". Specific dataframes for removal, transport and production process and compilation of all in interactions matrix
    allLosses_df= pd.DataFrame(np.diag(diag_list), index=SpeciesList , columns= SpeciesList)
    interactions_df = pd.DataFrame(np.diag(diag_list), index=SpeciesList , columns= SpeciesList)
    #removal_df = pd.DataFrame(np.diag(diag_list), index=SpeciesList , columns= SpeciesList)
    production_df=pd.DataFrame(index=SpeciesList , columns= SpeciesList)
    #transport_df=pd.DataFrame(index=SpeciesList , columns= SpeciesList)
    fragmentation_df=pd.DataFrame(index=SpeciesList , columns= SpeciesList)
    
    #fill in rest of interactions between species 
    for sp1 in interactions_df.index:
       for sp2 in interactions_df.columns:
           #interaction between same species are the losses (already added in diagonal)
           if sp1 == sp2:
               continue
           else:
               #Processess inside the same river section
               if sp1[:-3] == sp2[:-3]:
                   #only different size bins--> fragmentation
                   if (sp1[0:-1] == sp2[0:-1] and sp1[-1] != sp2[-1]):
                       #fragmentation only will occur from bigger to smaller and in consecutive sieBins
                       if (sp2[-1] =="b" and sp1[-1] =="a") or (sp2[-1] =="c" and sp1[-1] =="b") or (sp2[-1] =="d" and sp1[-1] =="c") or (sp2[-1] =="e" and sp1[-1] =="d"):
                           if type(RC_df[sp2]["fragmentation"])is tuple:
                             frag = RC_df[sp2]["fragmentation"]
                             fragval = frag[0]*frag[1]  
                             interactions_df[sp1][sp2] = fragval
                             production_df[sp1][sp2] = fragval
                             fragmentation_df[sp1][sp2] = fragval
                           else:
                             interactions_df[sp1][sp2] = RC_df[sp2]["fragmentation"] 
                             production_df[sp1][sp2] = RC_df[sp2]["fragmentation"] 
                             fragmentation_df[sp1][sp2] = RC_df[sp2]["fragmentation"] 
                       else:
                           interactions_df[sp1][sp2] = 0
                           production_df[sp1][sp2] = 0
                   #only different aggergation states--> heteroagg, biofouling, agg-breackup         
                   elif sp1[0:-2]+sp1[-1] == sp2[0:-2]+sp2[-1]:
                       if (sp2[-2] =="A" and sp1[-2] =="B") or (sp2[-2] =="C" and sp1[-2] =="D"):
                           interactions_df[sp1][sp2] = RC_df[sp2]["heteroagg"]
                           production_df[sp1][sp2] = RC_df[sp2]["heteroagg"]
                       elif (sp2[-2] =="B" and sp1[-2] =="A") or (sp2[-2] =="D" and sp1[-2] =="B"):
                           interactions_df[sp1][sp2] = RC_df[sp2]["breakup"]
                           production_df[sp1][sp2] = RC_df[sp2]["breakup"]
                       elif (sp2[-2] =="A" and sp1[-2] =="C") or (sp2[-2] =="B" and sp1[-2] =="D"):
                           interactions_df[sp1][sp2] = RC_df[sp2]["biofilm"]
                           production_df[sp1][sp2] = RC_df[sp2]["biofilm"]
                   #only different compartments-->settling, rising, mixing, resusp, burial?         
                   elif sp1[:-3]+sp1[-2:]  == sp2[:-3]+sp2[-2:]:
                       if (sp2[-3] =="1" and sp1[-3] =="2") or (sp2[-3] =="2" and sp1[-3] =="3") or (sp2[-3] =="3" and sp1[-3] =="4"):
                           interactions_df[sp1][sp2] = RC_df[sp2]["settling"]+ RC_df[sp2]["mixing"]
                           production_df[sp1][sp2] = RC_df[sp2]["settling"]+ RC_df[sp2]["mixing"]
                       elif (sp2[-3] =="3" and sp1[-3] =="2") or (sp2[-3] =="2" and sp1[-3] =="1"):
                           interactions_df[sp1][sp2] = RC_df[sp2]["rising"] + RC_df[sp2]["mixing"]
                           production_df[sp1][sp2] = RC_df[sp2]["rising"] + RC_df[sp2]["mixing"]
                       elif (sp2[-3] =="4" and sp1[-3] =="3"):
                           interactions_df[sp1][sp2] = RC_df[sp2]["resusp"]
                           production_df[sp1][sp2] = RC_df[sp2]["resusp"]
                
                   
              #If different river sections 
              ####Transport between river sections if consecutive and also                 
               elif int(sp2[:-3])+1 == int(sp1[:-3]):#Consecutive river sections
                       #Transport only in same species, compartment and size bin
                       if sp1[-3:] == sp2[-3:]:
                           if sp1[-3]!="4":
                               interactions_df[sp1][sp2] = RC_df[sp2]["advection"]
                           else:
                               interactions_df[sp1][sp2] = RC_df[sp2]["sedTransport"]
                       else:
                           interactions_df[sp1][sp2] = 0 
               
               
                      
    # for a in transport_df:
    #     transport_df[a] = transport_df[a].fillna(0)
    for b in production_df:
        production_df[b] = production_df[b].fillna(0)
    # for c in removal_df:
    #     removal_df[c] = removal_df[c].fillna(0)
    for d in fragmentation_df:
         fragmentation_df[d] = fragmentation_df[d].fillna(0)
    
               
    return [interactions_df, production_df,allLosses_df,fragmentation_df]        

def fillMatricesOfIrrevLoss(RC_df,Clist,compartments_prop):
    
    #Remove volume and density from RC_df
    
    RC_df= RC_df.drop("volume_m3")
    RC_df= RC_df.drop("density_kg_m3")
     
    #remove C from Clist:
    SpeciesList=[]
    for a in Clist:
        SpeciesList.append(a[2:])
    
    #Diagonal of the df corresponds to the losses of each species
    #crate the array of values for the diagonal wich is the sum of all RC corresponding to one species:
    diag_list = []
    
    for g in SpeciesList:
        if g[:-3]== Clist[-1][2:-3]: #If is the last river section advection is also an irreversible loss process
            if type(RC_df[g]["fragmentation"]) is tuple:
                # when there are fragements formed then fragmentation is not an irreversible loss process
            
                diag_list.append(-(RC_df[g]["degradation"]+RC_df[g]["burial"]+RC_df[g]["advection"]+RC_df[g]["sedTransport"]))
            
            else:
                diag_list.append(-(RC_df[g]["degradation"]+RC_df[g]["burial"]+RC_df[g]["fragmentation"]+RC_df[g]["advection"]+RC_df[g]["sedTransport"]))
                # diag_array = np.asarray(diag_list)
                # diag_array[np.isnan(diag_array)] = 0 #substitute nan by 0
        else:
            if type(RC_df[g]["fragmentation"]) is tuple:
                # when there are fragements formed then fragmentation is not an irreversible loss process
            
                diag_list.append(-(RC_df[g]["degradation"]+RC_df[g]["burial"]))
            
            else:
                diag_list.append(-(RC_df[g]["degradation"]+RC_df[g]["burial"]+RC_df[g]["fragmentation"]))
                # diag_array = np.asarray(diag_list)
                # diag_array[np.isnan(diag_array)] = 0 #substitute nan by 0
    #Dataframes to indicate interactions between "Species". Specific dataframes for removal, transport and production process and compilation of all in interactions matrix
    IrrevLoss_df = pd.DataFrame(np.diag(diag_list), index=SpeciesList , columns= SpeciesList)       

    return IrrevLoss_df 