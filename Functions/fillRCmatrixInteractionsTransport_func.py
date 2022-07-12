# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:35:45 2020

@author: PradoDomercq
"""


###NEW INTERACTIONS MATRIX INCLUDING TRANSPORT

import numpy as np
import pandas as pd
#Script to fillin the initeractions matrix (interactions_df)



def fillMatrixInteractionsALL(RC_df, Clist,compartments_prop):
    #Timer for Interactions matrix
#    startTime = datetime.now()
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
            # when there are fragements formed extract the fragmentation array to multiply kfrag by number of fragments formed
            #IN the lossess (estimation of mass lost) process we dont need to multiply by the number of fragments!!
            frag = RC_df[sp]["fragmentation"]
            #fragval = frag[0]*frag[1]
            
            RC_df_noFrag = RC_df[sp].drop("fragmentation")
            diag_list.append(-(sum(RC_df_noFrag)+frag[0]))
        
        else:
            diag_list.append(-(sum(RC_df[sp])))
            # diag_array = np.asarray(diag_list)
            # diag_array[np.isnan(diag_array)] = 0 #substitute nan by 0
    
    #Dataframe to indicate interactions between "Species" inside one river box
    interactions_df = pd.DataFrame(np.diag(diag_list), index=SpeciesList , columns= SpeciesList)
   
    
    #fill in rest of interactions between species 
    for sp1 in interactions_df.index:
       for sp2 in interactions_df.columns:
           
           #interaction between same species are the losses (already added in diagonal)
           if sp1 == sp2:
               continue
           else:
               #Processess inside the same river section (RS)
               if sp1[:-3] == sp2[:-3]:
                   #only different size bins--> fragmentation
                   if (sp1[:-1] == sp2[:-1] and sp1[-1] != sp2[-1]):
                       #fragmentation only will occur from bigger to smaller and in consecutive sieBins Sizebin = sp[-3]
                       if (sp2[-1] =="b" and sp1[-1] =="a") or (sp2[-1] =="c" and sp1[-1] =="b") or (sp2[-1] =="d" and sp1[-1] =="c") or (sp2[-1] =="e" and sp1[-1] =="d"):
                           if type(RC_df[sp2]["fragmentation"])is tuple:
                             frag = RC_df[sp2]["fragmentation"]
                             fragval = frag[0]*frag[1]  
                             interactions_df[sp1][sp2] = fragval
                           else:
                             interactions_df[sp1][sp2] = RC_df[sp2]["fragmentation"] 
                       else:
                           interactions_df[sp1][sp2] = 0
                   #only different aggergation states--> heteroagg, biofouling, agg-breackup (aggregation state = sp[-4])    
                   elif sp1[:-2]+sp1[-1] == sp2[:-2]+sp2[-1]:
                       if (sp2[-2] =="A" and sp1[-2] =="B") or (sp2[-2] =="C" and sp1[-2] =="D"):
                           interactions_df[sp1][sp2] = RC_df[sp2]["heteroagg"]
                       elif (sp2[-2] =="B" and sp1[-2] =="A") or (sp2[-2] =="D" and sp1[-2] =="B"):
                           interactions_df[sp1][sp2] = RC_df[sp2]["breakup"]
                       elif (sp2[-2] =="A" and sp1[-2] =="C") or (sp2[-2] =="B" and sp1[-2] =="D"):
                           interactions_df[sp1][sp2] = RC_df[sp2]["biofilm"]
                           
                   #only different compartments-->settling, rising, mixing, resusp  (compartment sp[-5])     
                   elif sp1[:-3]+sp1[-2:]  == sp2[:-3]+sp2[-2:]:
                       if (sp2[-3] =="1" and sp1[-3] =="2") or (sp2[-3] =="2" and sp1[-3] =="3") or (sp2[-3] =="3" and sp1[-3] =="4"):
                           interactions_df[sp1][sp2] = RC_df[sp2]["settling"]+ RC_df[sp2]["mixing"]
                       elif (sp2[-3] =="3" and sp1[-3] =="2") or (sp2[-3] =="2" and sp1[-3] =="1"):
                           interactions_df[sp1][sp2] = RC_df[sp2]["rising"] + RC_df[sp2]["mixing"]
                       elif (sp2[-3] =="4" and sp1[-3] =="3"):
                           interactions_df[sp1][sp2] = RC_df[sp2]["resusp"]
                               
                   else:
                      interactions_df[sp1][sp2] = 0 
                   
              #If different river sections:
              ####Transport between river sections only if consecutive river sections
               elif int(sp2[:-3])+1 == int(sp1[:-3]):#Consecutive river sections
                       #And Transport only in same species, compartment and size bin
                       if sp1[-3:] == sp2[-3:]:
                           if sp1[-3]!="4":
                               interactions_df[sp1][sp2] = RC_df[sp2]["advection"]
                           else:
                               interactions_df[sp1][sp2] = RC_df[sp2]["sedTransport"]
                       else:
                           interactions_df[sp1][sp2] = 0 
               else:
                   interactions_df[sp1][sp2] = 0 
               
                      
#    print(datetime.now() - startTime)           
    return interactions_df               