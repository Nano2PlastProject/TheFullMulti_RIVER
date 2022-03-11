# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 11:19:34 2021

@author: PradoDomercq
"""
import numpy as np
import pandas as pd

def fillInteractions_fun (RC_df, Clist,compartments_prop,river_flows):
    
    
    def interact3(sp1,interactions_df,RC_df,river_flows):
        
        def transportProcess (sp1,sp2,RC_df,river_flows):
            
            #If different river sections:
            ####Transport between river sections only if same compartment, species and size fraction
            if sp2[-3:] == sp1[-3:]:#Only different box (but all rest same)
                J=int(sp1[-4])+1
                I=int(sp2[-4])+1
                flowI_df=river_flows[river_flows.Region_I == I]
                if J in flowI_df.Region_J.tolist():
                    if sp1[-3]!="4":
                        if isinstance(RC_df[sp2]["advection"], (int, float)):
                            solution = RC_df[sp2]["advection"]
                        else:
                            idx_ad= np.where(flowI_df.Region_J== J)[0][0]
                            solution = RC_df[sp2]["advection"][idx_ad]
                    else:
                        solution = RC_df[sp2]["sedTransport"] 
                else:
                    solution = 0
            else:
                solution = 0
            return solution
        
        def inboxProcess(sp1,sp2,RC_df):
            if (sp1[:-1] == sp2[:-1] and sp1[-1] != sp2[-1]):
            #fragmentation only will occur from bigger to smaller and in consecutive sieBins Sizebin = sp[-3]
                if (sp2[-1] =="b" and sp1[-1] =="a") or (sp2[-1] =="c" and sp1[-1] =="b") or (sp2[-1] =="d" and sp1[-1] =="c") or (sp2[-1] =="e" and sp1[-1] =="d"):
                    if type(RC_df[sp2]["fragmentation"])is tuple:
                        frag = RC_df[sp2]["fragmentation"]
                        fragval = frag[0]*frag[1]  
                        sol = fragval
                    else:
                        sol = RC_df[sp2]["fragmentation"] 
                else:
                    sol = 0
             #only different aggergation states--> heteroagg, biofouling,defouling and agg-breackup (aggregation state = sp[-4])    
            elif sp1[:-2]+sp1[-1] == sp2[:-2]+sp2[-1]:
                if (sp2[-2] =="A" and sp1[-2] =="B") or (sp2[-2] =="C" and sp1[-2] =="D"):
                    sol = RC_df[sp2]["heteroagg"]
                elif (sp2[-2] =="B" and sp1[-2] =="A") or (sp2[-2] =="D" and sp1[-2] =="B"):
                    sol = RC_df[sp2]["breakup"]
                elif (sp2[-2] =="A" and sp1[-2] =="C") or (sp2[-2] =="B" and sp1[-2] =="D"):
                    sol = RC_df[sp2]["biofilm"]
                elif (sp2[-2] =="C" and sp1[-2] =="A") or (sp2[-2] =="D" and sp1[-2] =="B"):
                    sol = RC_df[sp2]["defouling"]
                else:
                    sol=0
            #only different compartments-->settling, rising, mixing, resusp  (compartment sp[-5])     
            elif sp1[:-3]+sp1[-2:]  == sp2[:-3]+sp2[-2:]:
                if (sp2[-3] =="1" and sp1[-3] =="2") or (sp2[-3] =="2" and sp1[-3] =="3") or (sp2[-3] =="3" and sp1[-3] =="4"):
                   sol = RC_df[sp2]["settling"]+ RC_df[sp2]["mixing"]
                elif (sp2[-3] =="3" and sp1[-3] =="2") or (sp2[-3] =="2" and sp1[-3] =="1"):
                    sol = RC_df[sp2]["rising"] + RC_df[sp2]["mixing"]
                elif (sp2[-3] =="4" and sp1[-3] =="3"):
                    sol = RC_df[sp2]["resusp"]
                else:
                    sol=0
            else:
                sol=0
            return sol    
        
        interact_row=[]
        for sp2 in interactions_df.index.to_list():
        
            if sp1 == sp2:
                interact_row.append(interactions_df[sp2][sp1])
            else:#Processess inside the same river section (RS)
                if sp1[:-3] == sp2[:-3]:
                    interact_row.append(inboxProcess(sp1,sp2,RC_df))
    
              #If different river sections:
              ####Transport between river sections is indicated by the Ocean flow matrix
    
                elif sp2[-3:] == sp1[-3:]:#Only different river section (but all rest same)
                    interact_row.append(transportProcess (sp1,sp2,RC_df,river_flows))
                else:
                    interact_row.append(0)
            
        return interact_row
    
    
    
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
        if isinstance(RC_df[sp]["advection"], (int, float)):
            if type(RC_df[sp]["fragmentation"]) is tuple:
            # when there are fragements formed extract the fragmentation array to multiply kfrag by number of fragments formed (fro the mas formed)
            #IN the lossess (estimation of mass lost) process we dont need to multiply by the number of fragments!!
                frag = RC_df[sp]["fragmentation"]
                RC_df_noFrag = RC_df[sp].drop("fragmentation")
                diag_list.append(-(sum(RC_df_noFrag)+frag[0]))
            else:
                diag_list.append(-(sum(RC_df[sp])))
        else:    
           adv=RC_df[sp]["advection"]
           RC_noAdv=RC_df[sp].drop("advection")
        
           if type(RC_df[sp]["fragmentation"]) is tuple:
               frag = RC_df[sp]["fragmentation"]
            #fragval = frag[0]*frag[1]
               RC_df_noFrag = RC_noAdv.drop("fragmentation")
               diag_list.append(-(sum(RC_df_noFrag)+frag[0]+sum(adv)))
        
           else:
               diag_list.append(-(sum(RC_noAdv)+sum(adv)))
        
                
                
    #Dataframe with lossess rate constants
    interactions_df = pd.DataFrame(np.diag(diag_list), index=SpeciesList , columns= SpeciesList)
    
    list_sp1=interactions_df.index.to_list()

    interactions_df_rows=[]
    for sp1 in interactions_df.index.to_list():
        interactions_df_rows.append(interact3(sp1,interactions_df,RC_df,river_flows))
    #interact3(sp1) for sp1 in interactions_df.index.to_list()]
    array=np.column_stack(interactions_df_rows)
    interactions_df_sol=pd.DataFrame(array,index=list_sp1, columns=list_sp1)
    
    return interactions_df_sol
