# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 11:28:19 2021

@author: PradoDomercq
"""
import numpy as np
def reshape_RC_df(RC_df,CombList):
    ### Reshape table of RC
    RC_dfShort=RC_df[CombList[0:(4*4*5)]]
    RC_dfShort= RC_dfShort.drop("volume_m3")
    RC_dfShort= RC_dfShort.drop("density_kg_m3")
    RC_dfShort
    RC_df_stacked = RC_dfShort.stack(dropna=False)
    RC_df_stacked
    RC_df_long = RC_df_stacked.reset_index()
    RC_df_long['Compartment'] = [x[-3] for x in RC_df_long['level_1']]
    RC_df_long['MP_form'] = [x[-2] for x in RC_df_long['level_1']]
    RC_df_long['SizeFrac'] = [x[-1] for x in RC_df_long['level_1']]
    RC_df_long['k_s-1'] = [x[0] if isinstance(x, tuple) else x for x in RC_df_long[0]]
    RC_df_long['t1/2_h'] = [round(np.log(2)/x/60/60,2) if x!=0 else x for x in RC_df_long['k_s-1'] ]
    RC_df_long['t_h'] = [round(1/x/60/60,2) if x!=0 else x for x in RC_df_long['k_s-1'] ]
    RC_df_long['t1/2_days'] = [round(np.log(2)/x/60/60/24, 2) if x!=0 else x for x in RC_df_long['k_s-1'] ]
    RC_df_long['Process'] = RC_df_long["level_0"]
    RC_df_final=RC_df_long.drop([0 ,'level_1',"level_0"], axis=1)
    RC_df_final.loc[(RC_df_final.Compartment == "1"),'Compartment']='Surface Water'
    RC_df_final.loc[(RC_df_final.Compartment == "2"),'Compartment']='Flowing Water'
    RC_df_final.loc[(RC_df_final.Compartment == "3"),'Compartment']='Stagnant Water'
    RC_df_final.loc[(RC_df_final.Compartment == "4"),'Compartment']='Sediment'
    RC_df_final.loc[(RC_df_final.MP_form == "A"),'MP_form']='Pristine'
    RC_df_final.loc[(RC_df_final.MP_form == "B"),'MP_form']='Heteroaggregate'
    RC_df_final.loc[(RC_df_final.MP_form == "C"),'MP_form']='Biofouled'
    RC_df_final.loc[(RC_df_final.MP_form == "D"),'MP_form']='Biof&Heteroagg'
    RC_df_final.loc[(RC_df_final.SizeFrac == "a"),'SizeFrac']='0.1 um'
    RC_df_final.loc[(RC_df_final.SizeFrac == "b"),'SizeFrac']='1 um'
    RC_df_final.loc[(RC_df_final.SizeFrac == "c"),'SizeFrac']='10 um'
    RC_df_final.loc[(RC_df_final.SizeFrac == "d"),'SizeFrac']='100 um'
    RC_df_final.loc[(RC_df_final.SizeFrac == "e"),'SizeFrac']='1000 um'
    RC_df_final.pivot_table(index=["Compartment", "MP_form", "SizeFrac"],columns='Process', values='k_s-1')
    
    return RC_df_final