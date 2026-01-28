# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 14:51:10 2021

@author: u0139791
Marit Hendrickx

"""

## if multithreading (parallel running)
## via (anaconda) cmd (with correct folder):
## cd C:/Users/u0139791/OneDrive - KU Leuven/Documents/GitHub/Auto-all-fields/Validation 21-22-23
# python Final_mcmc_figure.py 

import os
import time

#### FILL IN YOUR OWN DIRECTORY ####
main_dir='C:\\Users\\u0139791\\OneDrive - KU Leuven\\Documents\\GitHub\\Auto-all-fields\\Validation 21-22-23'
os.chdir(main_dir)

import numpy as np
import pandas as pd
from math import *
from datetime import datetime
import csv

from SWB_model import SWB

def ConvertToSerialDate(date):
    ''' date as datetime(2020,2,6) '''
    temp = datetime(1899, 12, 30)
    delta = date - temp
    return float(delta.days) + (float(delta.seconds) / 86400)

# sla files op met "_simirr" voor irrigatiesimulatie
# simirr=True 

version='v1' #TODO
v_folder='Validatie '+version
if version=='v8': v_folder='Validation v8'
main_dir_v=os.path.join(main_dir, v_folder)
os.chdir(main_dir_v)
listdir = os.listdir()

## PLV 2023
year='2023' #TODO
case = 'C58157' #TODO

## PCG 2022
# year='2022' #TODO
# case = 'C54986' #TODO

# folders_case = [x for x in listdir if case in x and not '.' in x and not 'Timeline' in x and len(x)>23]

if version=='v8': DREAM_log = pd.read_excel('DREAM_auto_log_'+version+'_full_'+year+'.xlsx')
else: DREAM_log = pd.read_excel('DREAM_auto_log_'+version+'_all_'+year+'.xlsx')
folders = DREAM_log.Folder.tolist()
folders = [x for x in folders if isinstance(x, str)]
folders_case = [i for i in folders if i.startswith(case)]

folders_case=folders_case[-2:]

#%% FUNCTION

def final_func(year,df_list,folder,folder2,cal_par_on,cal,forecast,multiplier_bounds,simirr=False,sim_NoMultiplier=False):
    if 'opkomst' in folder: opkomst=True
    else: opkomst=False
    if opkomst: samp='Mean5'
    else: samp='Mean30'
    
    foldersplit = folder.split('_')
    case = foldersplit[0]
    if opkomst: crop_name = foldersplit[1]+'_opkomst'
    else: crop_name = foldersplit[1]
    
    folder_path=os.path.join(main_dir_v, folder)
    if folder2=='':
        save_path=folder_path
    else:
        save_path=os.path.join(folder_path, folder2)
        
    if cal_par_on:
        par_names=['Kcb_ini','Kcb_mid','Kcb_end','L_ini','L_dev','L_mid','fc','log(Ksat)','CN','GWT_max','Zr_max','v_ini','Interc a','Slope b']
    else:
        par_names=['Kcb_ini','Kcb_mid','Kcb_end','L_ini','L_dev','L_mid','fc','log(Ksat)','CN','GWT_max','Zr_max','v_ini']
    
    sim_UncertainMultiplier = False #TODO
    
    n_mult=len(multiplier_bounds)
    if n_mult>0:
        par_names=par_names+['M'+str(x) for x in range(1,n_mult+1)]
    
    if n_mult>0 and sim_NoMultiplier:
        folder_sim_NoMultiplier = 'sim_NoMultiplier'
        save_path=os.path.join(folder_path, folder_sim_NoMultiplier)
        if not os.path.exists(save_path):
            os.mkdir(save_path)
        
        validation_info = pd.read_csv('validation_info.csv')
        last_valdate = validation_info['last_date'][0]
        # mult_lastIndex = np.where(multiplier_bounds['Upper bound date']<=last_valdate)[0][-1] 
        # --> laatste multiplier group volledig binnen calibration period
        mult_lastIndex = np.where(multiplier_bounds['Lower bound date']<=last_valdate)[0][-1] 
        # --> laatste multiplier group volledig OF DEELS binnen calibration period
        print(mult_lastIndex)
        n_pred = n_mult-mult_lastIndex-1
        
    elif n_mult>0 and sim_UncertainMultiplier:
        folder_sim_NoMultiplier = 'sim_UncertainMultiplier'
        save_path=os.path.join(folder_path, folder_sim_NoMultiplier)
        if not os.path.exists(save_path):
            os.mkdir(save_path)
        
        validation_info = pd.read_csv('validation_info.csv')
        last_valdate = validation_info['last_date'][0]
        mult_lastIndex = np.where(multiplier_bounds['Upper bound date']<=last_valdate)[0][-1] 
        # --> laatste multiplier group volledig binnen calibration period
        # mult_lastIndex = np.where(multiplier_bounds['Lower bound date']<=last_valdate)[0][-1] 
        # --> laatste multiplier group volledig OF DEELS binnen calibration period
        print(mult_lastIndex)        
        n_pred = n_mult-mult_lastIndex-1
    else: n_pred=0
    par=len(par_names)
    name=case
    
    os.chdir(folder_path)
    ParSet=np.load('ParSet_'+name+year+crop_name+'.npy')
    ParSetMax=np.load('ParSet_MaxLL_'+name+year+crop_name+'.npy')
    # fx=np.load('fx_'+name+year+crop_name+'.npy')
    
    istart=np.round(ParSet.shape[0]*0.5).astype('int')
    ParSet50 = ParSet[istart:,:-2] # Parameter values excl 50% burn-in
    excl = np.round(ParSet50.shape[0]*0.025).astype('int') # excl upper and lower 2.5% to obtain 95% PI
    sort=ParSet50.argsort(axis=0)
    ParSet50_sort = ParSet50[sort, np.arange(sort.shape[1])] # every column is now sorted
    ParSet_95 = ParSet50_sort[excl:-excl,:] # this is the 95% for each parameter
    ParSet_95_min=np.min(ParSet_95,axis=0) # lower 95% boundary
    ParSet_95_max=np.max(ParSet_95,axis=0) # upper 95% boundary
    
    
    ## Keep only last x parsets
    # Pmult ensemble = 50
    # 50% of all 36000 parsets = 18000
    # 18000/50 = 360 --> equal amount of sim ensemble members
    ParSet50_360 = ParSet50[-360:,:-2] # Last 360 parsets (is 60 parsets of each chain)
    
    # sample 10% random rows (within the last 50%) --> 1800
    ParSet50_10 = pd.DataFrame(ParSet50).sample(frac=0.1, random_state=42)  # Use random_state for reproducibility
    # sample 5% random rows (within the last 50%) --> 900
    ParSet50_5 = pd.DataFrame(ParSet50).sample(frac=0.05, random_state=42)  # Use random_state for reproducibility
    
    
    if n_mult>0 and sim_UncertainMultiplier:
        import warnings
        warnings.filterwarnings("ignore", category=FutureWarning)
        
        N_ens = len(ParSet50) # (10x1800 or 20x900)
        mult_ens = pd.DataFrame(columns=[x for x in range(n_pred)])
        for n in range(N_ens):
            random_values = sample_LaplaceTrunc(n,n_pred)
            mult_ens = mult_ens.append(pd.Series(random_values,index=[x for x in range(n_pred)]), ignore_index=True)

    
    global date_YYYYMMDD
    date_YYYYMMDD='20231002' # TODO
    
    if folder2=='':
        save_path=folder_path
        folder_save = v_folder+'\\'+folder
    else:
        save_path=os.path.join(folder_path, folder2)
        folder_save = v_folder+'\\'+folder+'\\'+folder2

    if not os.path.exists(save_path):
        os.mkdir(save_path)
    
    os.chdir(folder_path)
    if simirr:
        SWC_df = pd.read_csv('SWC_df_simirr_'+date_YYYYMMDD+'.csv')
        sw_df = pd.read_csv('sw_df_simirr_'+date_YYYYMMDD+'.csv')
        sw_vkrit_df = pd.read_csv('sw_krit_df_simirr_'+date_YYYYMMDD+'.csv')
        vkrit_df = pd.read_csv('vkrit_df_simirr_'+date_YYYYMMDD+'.csv')
    else:
        SWC_df = pd.read_csv('SWC_df_'+date_YYYYMMDD+'.csv')
        sw_df = pd.read_csv('sw_df_'+date_YYYYMMDD+'.csv')
        sw_vkrit_df = pd.read_csv('sw_krit_df_'+date_YYYYMMDD+'.csv')
        vkrit_df = pd.read_csv('vkrit_df_'+date_YYYYMMDD+'.csv')
    
    SWC_sort = SWC_df.sort_values(by=SWC_df.columns.to_list()).reset_index(drop=True)
    sw_sort = sw_df.sort_values(by=sw_df.columns.to_list()).reset_index(drop=True)
    # every column is now sorted
    if simirr: sw_sort.to_csv('sw_sorted_df_simirr_'+date_YYYYMMDD+'.csv')
    else: sw_sort.to_csv('sw_sorted_df_'+date_YYYYMMDD+'.csv')
    
    excl = np.round(SWC_df.shape[0]*0.025).astype('int') 
    # excl upper and lower 2.5% to obtain 95% PI
    SWC_95 = SWC_sort[excl:-excl]
    sw_95 = sw_sort[excl:-excl]
    
    ## All (100%) instead of 95%
    # SWC_95 = SWC_sort
    # sw_95 = sw_sort
    ##
    
    SWC_95_min = np.min(SWC_95)
    SWC_95_max = np.max(SWC_95)
    sw_95_min = np.min(sw_95)
    sw_95_max = np.max(sw_95)
    
    CI_SWC = np.array([SWC_95_min,SWC_95_max]).reshape(2,-1)
    CI_sw = np.array([sw_95_min,sw_95_max]).reshape(2,-1)
    CI=np.array([CI_SWC,CI_sw]).reshape(4,-1)
    
    ## Include CI around simulation and enter forecast
    # CI period is limited to the simulation period with available weather data
    # TO DO: enter forecast uncertainty
    if cal_par_on:
        sensor_data_adj=np.load('sensor_data_adj_'+name+year+crop_name+'.npy')
        print("y = ",ParSetMax[7],"x +",ParSetMax[6])
        sensor_cal=ParSetMax[7]*sensor_data_adj[:,1] + ParSetMax[6]
    else:
        sensor_cal=np.empty(0)
    
    ## Calibration time
    sensor = pd.read_csv('sensordata_'+name+year+crop_name+'_'+date_YYYYMMDD+'.csv')
    g_list = pd.read_csv('g_list_'+name+year+crop_name+'_'+date_YYYYMMDD+'.csv')
    meas = pd.read_csv('Measurements_used_DREAM.csv')
    N = meas['N_sensor'][0]
    if N==0: N=len(sensor)
    s_first = sensor['Date'][0]
    s_last = sensor['Date'][N-1]
    s_pred = s_last+10
    if N<len(sensor):
        where = np.where(g_list==s_pred)[0]
        if where.size:
            last_i = where[0]
            CI = CI[:,:last_i+1]
            print(CI.shape[1])
            if CI.shape[1]==len(g_list):
                CI = np.column_stack((CI, [np.nan]*CI.shape[0]))
                print(CI.shape[1])
        else: 
            print(CI.shape[1])
            CI = np.column_stack((CI, [np.nan]*CI.shape[0]))
            print(CI.shape[1])
    
    
    ## interval (CI) komt van hierboven; 
    # plot blauwe lijn = max likelihood, met blauwe band (CI)
    os.chdir(main_dir)
    
    if n_mult>0 and sim_NoMultiplier:
        folder_save_sim_NoMultiplier = folder+'\\sim_NoMultiplier'
        SWC, sw_list, g_list, sensor, covar, df_obs, df_all, df_parameters =SWB(ParSetMax[0],ParSetMax[1],ParSetMax[2], \
             ParSetMax[3],ParSetMax[4],ParSetMax[5],ParSetMax[6],ParSetMax[7],ParSetMax[8],ParSetMax[9],ParSetMax[10],ParSetMax[11], \
             [ParSetMax[x] for x in range(12,12+mult_lastIndex+1)]+[0 for x in range(12+mult_lastIndex+1,12+n_mult)], \
             multiplier_bounds,sensor=True,cal=cal,sensor_cal=sensor_cal,CI=CI,show=[True,folder_save_sim_NoMultiplier], \
             case=name,year=year,forecast=forecast,df_list=df_list)
    else:
        SWC, sw_list, g_list, sensor, covar, df_obs, df_all, df_parameters =SWB(ParSetMax[0],ParSetMax[1],ParSetMax[2], \
                 ParSetMax[3],ParSetMax[4],ParSetMax[5],ParSetMax[6],ParSetMax[7],ParSetMax[8],ParSetMax[9],ParSetMax[10],ParSetMax[11], \
                 [ParSetMax[x] for x in range(12,12+n_mult)],multiplier_bounds, \
                 sensor=True,cal=cal,sensor_cal=sensor_cal,CI=CI,show=[True,folder_save], \
                 case=name,year=year,forecast=forecast,df_list=df_list)
    


#%% run_here: for loop

for folder in folders_case:
    global date_YYYYMMDD
    date_YYYYMMDD='20231002'
    folder2 = 'Extra figures_v0' # 'irr_0.78' # 'no_irr' 'sim_irr' #TODO
    cal_par_on=False
    forecast_on=False #TODO
    multipliers=False #TODO
    sim_NoMultiplier = False # if True, multiplier is zero (factor 1), without uncertainty #TODO
    sim_UncertainMultiplier = False
    simirr=False
    
    os.chdir(main_dir)
    
    ## read csv data
    df1 = pd.read_csv('soildata_212223.csv', encoding= 'unicode_escape')
    df1['year'] = df1['year'].fillna(0).astype(int).astype(str)
    df1['year'] = df1['year'].replace('0', '')
    df1 = df1[df1['year']==year]
    df1.reset_index(inplace=True,drop=True)
    df2 = pd.read_csv('crop_FAO.csv', encoding= 'unicode_escape')
    
    df_ET = pd.read_csv('ETo_212223.csv', encoding= 'unicode_escape')
    df_ET = df_ET[df_ET['year'].astype(str)==year]
    df_ET.reset_index(inplace=True,drop=True)
    df_R = pd.read_csv('Neerslag_ALL_filled_212223_adjusted.csv', encoding= 'unicode_escape')
    df_R = df_R[df_R['year'].astype(str)==year]
    df_R.reset_index(inplace=True,drop=True)
    df_I = pd.read_csv('Irrigatie_212223.csv', encoding= 'unicode_escape')
    df_I = df_I[df_I['year'].astype(str)==year]
    df_I.reset_index(inplace=True,drop=True)
    df_obs=pd.read_csv('Soilobs_212223.csv', encoding= 'unicode_escape')
    df_obs = df_obs[df_obs['year'].astype(str)==year]
    df_obs.reset_index(inplace=True,drop=True)
    
    df_list = [df_R,df_I,df_ET,df1,df_obs,df2]
    
    if multipliers:
        folder_path=os.path.join(main_dir_v, folder)
        os.chdir(folder_path)  
        multiplier_bounds = pd.read_csv('multiplier_bounds.csv', encoding= 'unicode_escape')
    else: multiplier_bounds = []
    
    if forecast_on:
        dir_forecasts='C:\\Users\\u0139791\\OneDrive - KU Leuven\\Documents\\Python Scripts\\ECMWF API\\Fields'+year
        os.chdir(dir_forecasts)
        forecast_R=pd.read_csv(case+'/cf_'+case+'_'+date_YYYYMMDD+'_daily_R.csv', encoding= 'unicode_escape')
        forecast_R.loc[~(forecast_R['Daily precipitation [mm/day]'] > 0), 'Daily precipitation [mm/day]']=0
        forecast_ETo=pd.read_csv(case+'/cf_'+case+'_'+date_YYYYMMDD+'_daily_ETo.csv', encoding= 'unicode_escape')
        forecast = pd.DataFrame({'Date':forecast_R['Date'],'R':forecast_R['Daily precipitation [mm/day]'],'ETo':forecast_ETo['Daily ETo [mm/day]']})
        forecast['Date'] = [int(ConvertToSerialDate(datetime.strptime(date, '%d-%m-%Y'))) \
                               for date in forecast['Date']]
        forecast['I']=0
        print(forecast) # in mm
        os.chdir(main_dir_v)
    else: 
        forecast=np.empty(0)
        forecast_date = np.nan
            
    
    ## Function
    
    final_func(year,df_list,folder,folder2,cal_par_on,'gen',forecast,
               multiplier_bounds,simirr=simirr,sim_NoMultiplier=sim_NoMultiplier)
    
        

