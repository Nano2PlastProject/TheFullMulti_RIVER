# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 11:19:34 2021

@author: PradoDomercq
"""
import numpy as np
import pandas as pd

def fillInteractions_fun_v2_0 (RC_df,Clist,river_flows):
    
    def interact3(sp1,interactions_df,RC_df,river_flows):
        sol=[]
        for sp2 in interactions_df.index.to_list():
            
            #Same particle in the same river section and compartment
            if sp1 == sp2:
                sol.append(interactions_df[sp2][sp1])
            
            #Different particle or different river section or compartment
            else:

                #Same river section (RS)--> In box processes
                if sp1[:-3] == sp2[:-3]:
                    sol.append(inboxProcess(sp1,sp2,RC_df))
    
                #Different river section but same particle in same compartment -->Transport (advection or sediment transport determined by flow_connectivity file)
                elif sp2[-3:] == sp1[-3:]:
                    sol.append(transportProcess (sp1,sp2,RC_df,river_flows))
                else:
                    sol.append(0)
            
        return sol
        
    def transportProcess (sp1,sp2,RC_df,river_flows):
        J=int(sp1[:-3])+1
        I=int(sp2[:-3])+1
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
    
        return solution
        
    def inboxProcess(sp1,sp2,RC_df):
        #only different size bins--> fragmentation
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

            #heteroaggregation from A-->B or from C-->D
            if (sp2[-2] =="A" and sp1[-2] =="B") or (sp2[-2] =="C" and sp1[-2] =="D"):
                sol = RC_df[sp2]["heteroagg"]

            #heteroaggregate breackup from B-->A and from D-->C  
            elif (sp2[-2] =="B" and sp1[-2] =="A") or (sp2[-2] =="D" and sp1[-2] =="C"):
                sol = RC_df[sp2]["breakup"]
            
            #Biofouling from A-->C or from B-->D
            elif (sp2[-2] =="A" and sp1[-2] =="C") or (sp2[-2] =="B" and sp1[-2] =="D"):
                sol = RC_df[sp2]["biofilm"]

            #Defouling from C-->A or from D-->B
            elif (sp2[-2] =="C" and sp1[-2] =="A") or (sp2[-2] =="D" and sp1[-2] =="B"):
                sol = RC_df[sp2]["defouling"]

            else:
                sol=0

        #only different compartments-->settling, rising, mixing, resusp  (compartment sp[-5])     
        elif sp1[:-3]+sp1[-2:]  == sp2[:-3]+sp2[-2:]:

            #Downwards movement between water compartments (settling and mixing)
            if (sp2[-3] =="1" and sp1[-3] =="2") or (sp2[-3] =="2" and sp1[-3] =="3"):
                if type(RC_df[sp2]["mixing"]) == tuple:
                    sol = RC_df[sp2]["settling"]+ RC_df[sp2]["mixing"][1]
                else:
                    sol = RC_df[sp2]["settling"]+ RC_df[sp2]["mixing"]

            #Downwards movement with sediment (only sttling)
            elif (sp2[-3] =="3" and sp1[-3] =="4"):
                sol = RC_df[sp2]["settling"]

            #Upwards movement water compartments (rising and mixing)
            elif (sp2[-3] =="3" and sp1[-3] =="2") or (sp2[-3] =="2" and sp1[-3] =="1"):
                if type(RC_df[sp2]["mixing"]) == tuple:
                    sol = RC_df[sp2]["rising"] + RC_df[sp2]["mixing"][0]
                else:
                    sol = RC_df[sp2]["rising"] + RC_df[sp2]["mixing"]
            #Upwards movement from the sediment
            elif (sp2[-3] =="4" and sp1[-3] =="3"):
                sol = RC_df[sp2]["resusp"]
            else:
                sol=0
        else:
            sol=0
        return sol    
    
    #Remove volume and density from RC_df
    RC_df= RC_df.drop("volume_m3")
    RC_df= RC_df.drop("density_kg_m3")
     
    #remove C from Clist:
    SpeciesList=[]
    for a in Clist:
        SpeciesList.append(a[2:])
    
    #Diagonal of the df corresponds to the losses of each species
    #crate the array of values for the diagonal wich is the sum of all RC corresponding to one species:

    #Estimate losses (diagonal)
    diag_list = []
    for sp in SpeciesList:
        RC_sp=RC_df[sp].to_frame()
        losses=[]
        for i in range(len(RC_sp[sp])):
            if type(RC_sp[sp][i])== tuple:
                if (RC_sp.index[i] == "advection") or (RC_sp.index[i]=="fragmentation") :
                    losses.append(RC_sp[sp][i][0])
                elif RC_sp.index[i] == "mixing":
                    losses.append(sum(RC_sp[sp][i]))
                else: 
                    losses.append(RC_sp.index[i])
            else:
                losses.append(RC_sp[sp][i])
        diag_list.append(-(sum(losses)))

    interactions_df = pd.DataFrame(np.diag(diag_list), index=SpeciesList , columns= SpeciesList)

    #Estimate interactions
    interactions_df_rows=[]
    list_sp1=interactions_df.index.to_list()

    for sp1 in list_sp1:
        interactions_df_rows.append(interact3(sp1,interactions_df,RC_df,river_flows))
    
    #interact3(sp1) for sp1 in interactions_df.index.to_list()]
    array=np.column_stack(interactions_df_rows)# vstack it was set as column stack and was wrong!!
    interactions_df_sol=pd.DataFrame(array,index=list_sp1, columns=list_sp1)

    return interactions_df_sol
