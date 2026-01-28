# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 14:51:10 2021

@author: u0139791
Marit Hendrickx

"""

## if multithreading (parallel running)
## via (anaconda) cmd (with correct folder):
## cd C:/Users/u0139791/OneDrive - KU Leuven/Documents/GitHub/Auto-all-fields/Validation 21-22-23
# python Final_mcmc_enfo.py 

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

def ConvertToDate(serial_date):
    serial_date=int(serial_date)
    return datetime.fromordinal(datetime(1899, 12, 30).toordinal() + serial_date)

# sla files op met "_simirr" voor irrigatiesimulatie
# simirr=True 

### 
global date_YYYYMMDD
## PLV Herent 2023
g0 = 45082
## PCG 2022
# g0 = 44692 #TODO
cal_days = 60 #TODO
g_c = g0+cal_days
date_YYYYMMDD = ConvertToDate(g_c).strftime('%Y%m%d')
# date_YYYYMMDD='20230625' # TODO
year='2023' #TODO

df_list = pd.read_excel('Validatie v1\\DREAM_auto_log_v1_all_'+year+'.xlsx', engine='openpyxl')

name = 'PLV_Herent' # 'PCG_Kruisem'  'PLV_Herent'
case = 'C58157' # PCG: 'C54986' # Herent: 'C58157'
folder = df_list[(df_list['Case']==case) * (df_list['Validatiedag']==int(cal_days))]['Folder'].values[0]

forecast_on=True 
forecast_type = 'CF' # ENS / CF / WC / OBS_0irr #TODO

maxLLH = False # use only maxLLH parameter set #TODO
if maxLLH: forecast_type = 'ENS_maxLLH'

folder2 = forecast_type+'_forecast'
if forecast_type == 'WC': worst_case = True
else: worst_case = False
if forecast_type == 'OBS_0irr': obs_case_0irr = True
else: obs_case_0irr = False


cal_par_on=False
multipliers=False #TODO
sim_NoMultiplier = False # if True, multiplier is zero (factor 1), without uncertainty #TODO
sim_UncertainMultiplier = False

#%%

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
    folder_path=os.path.join(main_dir, folder)
    os.chdir(folder_path)  
    multiplier_bounds = pd.read_csv('multiplier_bounds.csv', encoding= 'unicode_escape')
else: multiplier_bounds = []

simirr=False

#### Weather Forecasts

if forecast_on:
    dir_forecasts='C:\\Users\\u0139791\\OneDrive - KU Leuven\\Documents\\Python Scripts\\ECMWF API'
    os.chdir(dir_forecasts)
    name_folder = 'PLV_2023' # 'PCG_2022' 'PLV_2023' #TODO
    forecast_R_cf=pd.read_csv(name_folder+'/cf_'+name+'_'+date_YYYYMMDD+'_daily_R.csv', encoding= 'unicode_escape')
    forecast_ETo_cf=pd.read_csv(name_folder+'/cf_'+name+'_'+date_YYYYMMDD+'_daily_ETo.csv', encoding= 'unicode_escape')
    forecast_R_pf=pd.read_csv(name_folder+'/'+name+'_'+date_YYYYMMDD+'_daily_R_10days.csv', encoding= 'unicode_escape')
    forecast_ETo_pf=pd.read_csv(name_folder+'/'+name+'_'+date_YYYYMMDD+'_daily_ETo_10days.csv', encoding= 'unicode_escape')
    # all 50 perturbed forecasts + 1 control forecast in one dataframe:
    forecast_R_pf['50']=forecast_R_cf.iloc[:,1]
    forecast_ETo_pf['50']=forecast_ETo_cf.iloc[:,1]

else: 
    raise ValueError("ERROR: This script is made for forecasts.")

global val_v1
val_v1 = True #TODO
    
#%%
def final_func(year,df_list,folder,folder2,cal_par_on,cal,forecast_R,forecast_ETo,multiplier_bounds,simirr=False,sim_NoMultiplier=False):
    if 'opkomst' in folder: opkomst=True
    else: opkomst=False
    if opkomst: samp='Mean5'
    else: samp='Mean30'
    
    foldersplit = folder.split('_')
    case = foldersplit[0]
    if opkomst: crop_name = foldersplit[1]+'_opkomst'
    else: crop_name = foldersplit[1]
    
    folder_path=os.path.join(main_dir, folder)
    if val_v1: folder_path=os.path.join(main_dir,'Validatie v1', folder)
    
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
        folder_sim_NoMultiplier = 'sim_NoMultiplier' #TODO
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
        folder_sim_NoMultiplier = 'sim_UncertainMultiplier' #TODO
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
    else: 
        n_pred=0
        if not os.path.exists(save_path):
            os.mkdir(save_path)
            
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

    
    #%%
    with open('ParSetMax-95_'+name+year+crop_name+'.csv', mode='w',newline='') as file:
        file = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        file.writerow(['Parameter','Max LLH parameter value','Min 95% interval','Max 95% interval']) 
        file.writerows([par_names[i], ParSetMax[i],ParSet_95_min[i],ParSet_95_max[i]] for i in range(par)) 
        
    
    #%% Parallel computation of simulations -> CI (2 minutes)
    os.chdir(main_dir)
    
    #### SUBSAMPLING OF MODEL ENS
    # Stratified Random sampling
    # np.random.seed(40)  # Set seed for reproducibility
    all_resampled_indices = []
    for i in range(0,ParSet50.shape[0],int(ParSet50.shape[0]/6)):
        data = ParSet50[i:i+int(ParSet50.shape[0]/6),:]
        # parameters = data[:, :-2]  # Shape (N, 12)
        # log_likelihoods = data[:, -1]  # Shape (N,)
        ## subsampling #TODO
        # N = 100 # this is 5% of 2000
        N = 200 # this is 10% of 2000
        # N = 500 # this is 25% of 2000
        np.random.seed(40+i)  # Set seed for reproducibility
        resampled_indices = np.random.choice(len(data), size=N, replace=False) # Uniform random sampling
        all_resampled_indices.append(resampled_indices)
    # Combine all resampled indices into a single DataFrame
    df_resampled_indices = np.vstack(all_resampled_indices).flatten()
    
    ## only first stratum or all #TODO
    # df_resampled_indices = np.array(range(0,int(ParSet50.shape[0]/6))) # first stratum only
    # df_resampled_indices = np.array(range(0,int(ParSet50.shape[0]))) # all iterations
    
    #### Only highest likelihood parameter set
    if maxLLH: df_resampled_indices = np.array([np.argmax(ParSet[istart:,-1])])
    
    ## forecast
    forecast_cf = pd.DataFrame({'Date':forecast_R_pf['Date'],'R':forecast_R_pf['50'],
                             'ETo':forecast_ETo_pf['50']})
    forecast_cf['Date'] = [int(ConvertToSerialDate(datetime.strptime(date, '%d-%m-%Y'))) \
                           for date in forecast_cf['Date']]
    forecast_R_pf['Date'] = [int(ConvertToSerialDate(datetime.strptime(date, '%d-%m-%Y'))) \
                           for date in forecast_R_pf['Date']]  
    forecast_cf['I']=0 #TODO
    
    if worst_case:
        df_ET_case = df_ET[['Date',case]]
        df_ET_case = df_ET_case.loc[(df_ET_case['Date']<=np.min(forecast_R_pf['Date']))
                                    * (df_ET_case['Date']>=np.min(forecast_R_pf['Date'])-10)]
        forecast = pd.DataFrame({'Date':forecast_R_pf['Date'],'R':np.zeros(len(forecast_R_pf['Date'])),
                                 'ETo':np.ones(len(forecast_R_pf['Date']))*np.mean(df_ET_case[case])})
        forecast['I']=0 #TODO
        forecast_indices = range(1)
        
    elif obs_case_0irr:
        df_ET_case = df_ET[['Date',case]]
        df_ET_case = df_ET_case.loc[(df_ET_case['Date']<=np.max(forecast_R_pf['Date']))
                                    * (df_ET_case['Date']>=np.min(forecast_R_pf['Date']))]
        df_R_case = df_R[['Date',case]]
        df_R_case = df_R_case.loc[(df_R_case['Date']<=np.max(forecast_R_pf['Date']))
                                    * (df_R_case['Date']>=np.min(forecast_R_pf['Date']))]
        forecast = pd.DataFrame({'Date':df_R_case['Date'].reset_index(drop=True),
                                 'R':df_R_case[case].reset_index(drop=True),
                                 'ETo':df_ET_case[case].reset_index(drop=True)})
        forecast['I']=0 #TODO
        forecast_indices = range(1)
        
    elif forecast_type == 'CF': 
        forecast = forecast_cf
        forecast['I']=0 #TODO
        forecast_indices = range(1)
        
    elif forecast_type == 'ENS' or forecast_type == 'ENS_maxLLH': 
        forecast = pd.DataFrame([])
        forecast_indices = range(51)

        
    #### PARALLEL COMPUTATION subsampled model ens x forecast ens
    start=datetime.now()
    from joblib import Parallel, delayed
    parallel_jobs=int(os.cpu_count()-2)
    parallel_jobs=5 #TODO

    # def iter_process(i):
    def iter_process(i,f,forecast):
        
        if forecast_type == 'ENS' or forecast_type == 'ENS_maxLLH': 
            forecast_R_pf.loc[~(forecast_R_pf[str(f)] > 0), str(f)]=0
            forecast = pd.DataFrame({'Date':forecast_R_pf['Date'],'R':forecast_R_pf[str(f)],
                                 'ETo':forecast_ETo_pf[str(f)]})
            forecast['I']=0 #TODO
            
        print(forecast) # in mm
        
        SWC, sw_list, _, _, _, _, df_all, _ =SWB(ParSet50[i,0],ParSet50[i,1],ParSet50[i,2],ParSet50[i,3], \
                        ParSet50[i,4],ParSet50[i,5],ParSet50[i,6],ParSet50[i,7], \
                        ParSet50[i,8],ParSet50[i,9],ParSet50[i,10],ParSet50[i,11], \
                        [ParSet50[i,x] for x in range(12,12+n_mult)],multiplier_bounds, \
                        sensor=False,cal='',sensor_cal=np.empty(0), \
                        CI=np.empty(0),show=[False,''],case=name,year=year,forecast=forecast,df_list=df_list)
        
        
        if opkomst:
            sw_vkrit_list=df_all['ASW vkrit_5 (mm)']
            vkrit_list=df_all['vkrit_Ze']
        else:
            sw_vkrit_list=df_all['ASW vkrit (mm)']
            vkrit_list=df_all['vkrit']
            
        return SWC, sw_list, sw_vkrit_list, vkrit_list
    
    # SWC, sw_list, sw_vkrit_list, vkrit_list = zip(*Parallel(n_jobs=parallel_jobs)(delayed(iter_process)(i) for i in range(len(ParSet50))))
    SWC, sw_list, sw_vkrit_list, vkrit_list = zip(*Parallel(n_jobs=parallel_jobs)(delayed(iter_process)(i,f,forecast)  \
                                                  for i in df_resampled_indices for f in forecast_indices))

    print(datetime.now()-start)
    
    
    SWC_df = pd.DataFrame(SWC)
    sw_df = pd.DataFrame(sw_list)
    sw_vkrit_df = pd.DataFrame(sw_vkrit_list)
    vkrit_df = pd.DataFrame(vkrit_list)
    os.chdir(save_path)
    if simirr:
        SWC_df.to_csv('SWC_df_simirr_'+date_YYYYMMDD+'.csv',index=False)
        sw_df.to_csv('sw_df_simirr_'+date_YYYYMMDD+'.csv',index=False)
        sw_vkrit_df.to_csv('sw_krit_df_simirr_'+date_YYYYMMDD+'.csv',index=False)
        vkrit_df.to_csv('vkrit_df_simirr_'+date_YYYYMMDD+'.csv',index=False)
    else:
        SWC_df.to_csv('SWC_df_'+date_YYYYMMDD+'.csv',index=False)
        sw_df.to_csv('sw_df_'+date_YYYYMMDD+'.csv',index=False)
        sw_vkrit_df.to_csv('sw_krit_df_'+date_YYYYMMDD+'.csv',index=False)
        vkrit_df.to_csv('vkrit_df_'+date_YYYYMMDD+'.csv',index=False)
    
    SWC_sort = SWC_df.sort_values(by=SWC_df.columns.to_list()).reset_index(drop=True)
    sw_sort = sw_df.sort_values(by=sw_df.columns.to_list()).reset_index(drop=True)
    # every column is now sorted
    if simirr: sw_sort.to_csv('sw_sorted_df_simirr_'+date_YYYYMMDD+'.csv',index=False)
    else: sw_sort.to_csv('sw_sorted_df_'+date_YYYYMMDD+'.csv',index=False)
    
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
    
    os.chdir(main_dir)
    ## Include CI around simulation and enter forecast
    # CI period is limited to the simulation period with available weather data
    # TO DO: enter forecast uncertainty
    if cal_par_on:
        sensor_data_adj=np.load('sensor_data_adj_'+name+year+crop_name+'.npy')
        print("y = ",ParSetMax[7],"x +",ParSetMax[6])
        sensor_cal=ParSetMax[7]*sensor_data_adj[:,1] + ParSetMax[6]
    else:
        sensor_cal=np.empty(0)
    
            
    ## interval (CI) komt van hierboven; 
    # plot blauwe lijn = max likelihood, met blauwe band (CI)
    if folder2=='':
        folder_save = folder
    else:
        folder_save = folder+'\\'+folder2
        
    if val_v1: folder_save='Validatie v1\\'+folder_save #TODO
    
    # if n_mult>0 and sim_NoMultiplier:
    #     folder_save_sim_NoMultiplier = folder+'\\sim_NoMultiplier'
    #     SWC, sw_list, g_list, sensor, covar, df_obs, df_all, df_parameters =SWB(ParSetMax[0],ParSetMax[1],ParSetMax[2], \
    #          ParSetMax[3],ParSetMax[4],ParSetMax[5],ParSetMax[6],ParSetMax[7],ParSetMax[8],ParSetMax[9],ParSetMax[10],ParSetMax[11], \
    #          [ParSetMax[x] for x in range(12,12+mult_lastIndex+1)]+[0 for x in range(12+mult_lastIndex+1,12+n_mult)], \
    #          multiplier_bounds,sensor=True,cal=cal,sensor_cal=sensor_cal,CI=CI,show=[True,folder_save_sim_NoMultiplier], \
    #          case=name,year=year,forecast=forecast_cf,df_list=df_list)
    # else:
    if forecast_type == 'ENS' or forecast_type == 'ENS_maxLLH' or forecast_type == 'CF':
        SWC, sw_list, g_list, sensor, covar, df_obs, df_all, df_parameters =SWB(ParSetMax[0],ParSetMax[1],ParSetMax[2], \
                 ParSetMax[3],ParSetMax[4],ParSetMax[5],ParSetMax[6],ParSetMax[7],ParSetMax[8],ParSetMax[9],ParSetMax[10],ParSetMax[11], \
                 [ParSetMax[x] for x in range(12,12+n_mult)],multiplier_bounds, \
                 sensor=True,cal=cal,sensor_cal=sensor_cal,CI=CI,show=[True,folder_save], \
                 case=name,year=year,forecast=forecast_cf,df_list=df_list)
    elif forecast_type == 'WC' or forecast_type == 'OBS_0irr':    
        SWC, sw_list, g_list, sensor, covar, df_obs, df_all, df_parameters =SWB(ParSetMax[0],ParSetMax[1],ParSetMax[2], \
                 ParSetMax[3],ParSetMax[4],ParSetMax[5],ParSetMax[6],ParSetMax[7],ParSetMax[8],ParSetMax[9],ParSetMax[10],ParSetMax[11], \
                 [ParSetMax[x] for x in range(12,12+n_mult)],multiplier_bounds, \
                 sensor=True,cal=cal,sensor_cal=sensor_cal,CI=CI,show=[True,folder_save], \
                 case=name,year=year,forecast=forecast,df_list=df_list)
    
    sw_vkrit_list=df_all['ASW vkrit (mm)']
    vkrit_list=df_all['vkrit']

    
    os.chdir(save_path)
    pd.DataFrame(covar).to_csv('covar_'+name+year+crop_name+'_'+date_YYYYMMDD+'.csv',index=False)
    df_obs.to_csv('obs_'+name+year+crop_name+'_'+date_YYYYMMDD+'.csv',index=False)
    
    if simirr: 
        df_all.to_csv('df_all_'+name+year+crop_name+'_simirr_'+date_YYYYMMDD+'.csv',index=False)
        df_parameters.to_csv('df_parameters_'+name+year+crop_name+'_simirr_'+date_YYYYMMDD+'.csv',index=False)
        pd.DataFrame(sw_vkrit_list).to_csv('sw_vkrit_list_'+name+year+crop_name+'_simirr_'+date_YYYYMMDD+'.csv',index=False)
        pd.DataFrame(SWC).to_csv('SWC_'+name+year+crop_name+'_simirr_'+date_YYYYMMDD+'.csv',index=False)
        pd.DataFrame(sw_list).to_csv('sw_'+name+year+crop_name+'_simirr_'+date_YYYYMMDD+'.csv',index=False)
        pd.DataFrame(g_list).to_csv('g_list_'+name+year+crop_name+'_simirr_'+date_YYYYMMDD+'.csv',index=False)
        pd.DataFrame(CI).to_csv('CI_'+name+year+crop_name+'_simirr_'+date_YYYYMMDD+'.csv',index=False)
        pd.DataFrame(sensor).to_csv('sensordata_'+name+year+crop_name+'_simirr_'+date_YYYYMMDD+'.csv',index=False)
        if cal_par_on: 
            pd.DataFrame(sensor_cal).to_csv('sensor_calDREAM_'+name+year+crop_name+'_simirr_'+date_YYYYMMDD+'.csv',index=False)
     
    
    else: 
        df_all.to_csv('df_all_'+name+year+crop_name+'_'+date_YYYYMMDD+'.csv',index=False)
        df_parameters.to_csv('df_parameters_'+name+year+crop_name+'_'+date_YYYYMMDD+'.csv',index=False)
        pd.DataFrame(sw_vkrit_list).to_csv('sw_vkrit_list_'+name+year+crop_name+'_'+date_YYYYMMDD+'.csv',index=False)
        pd.DataFrame(SWC).to_csv('SWC_'+name+year+crop_name+'_'+date_YYYYMMDD+'.csv',index=False)
        pd.DataFrame(sw_list).to_csv('sw_'+name+year+crop_name+'_'+date_YYYYMMDD+'.csv',index=False)
        pd.DataFrame(g_list).to_csv('g_list_'+name+year+crop_name+'_'+date_YYYYMMDD+'.csv',index=False)
        pd.DataFrame(CI).to_csv('CI_'+name+year+crop_name+'_'+date_YYYYMMDD+'.csv',index=False)
        pd.DataFrame(sensor).to_csv('sensordata_'+name+year+crop_name+'_'+date_YYYYMMDD+'.csv',index=False)
        if cal_par_on: 
            pd.DataFrame(sensor_cal).to_csv('sensor_calDREAM_'+name+year+crop_name+'_'+date_YYYYMMDD+'.csv',index=False)
       
    
    ### SWC vs sensor, samples
    
    # get SWC simulations of the days with a sensor observation
    g_list=np.array(g_list)
    n=len(sensor['Date'])
    sim=np.empty(n)
    for i in range(0,n):
        index=np.where(g_list==sensor['Date'].iloc[i])
        sim[i]=SWC[index[0][0]]
    
    # compare SWC simulations to mean sensor values
    ME_sensor = np.mean(sensor.Mean-sim)
    MAE_sensor = np.mean(abs(sensor.Mean-sim))
    RMSE_sensor = np.sqrt(np.mean((sensor.Mean-sim)**2))
    
    # compare SWC simulations to soil samples (0-30)
    ME_samp = np.mean(df_obs[samp]-df_obs['Sim'])
    MAE_samp = np.mean(abs(df_obs[samp]-df_obs['Sim']))
    RMSE_samp = np.sqrt(np.mean((df_obs[samp]-df_obs['Sim'])**2))
    n_samp=len(df_obs[samp])
    
    # compare SWC simulations to sensor & soil samples (0-30)
    ME_all = (np.sum(df_obs[samp]-df_obs['Sim'])+np.sum(sensor.Mean-sim))/(n+n_samp)
    MAE_all = (np.sum(abs(df_obs[samp]-df_obs['Sim']))+np.sum(abs(sensor.Mean-sim)))/(n+n_samp)
    RMSE_all = np.sqrt((np.sum((df_obs[samp]-df_obs['Sim'])**2)+np.sum((sensor.Mean-sim)**2))/(n+n_samp))
    
    df_metrics = pd.DataFrame(data={'sensor':[n,ME_sensor,MAE_sensor,RMSE_sensor],
                                            'samples':[n_samp,ME_samp,MAE_samp,RMSE_samp],
                                            'both':[n+n_samp,ME_all,MAE_all,RMSE_all]},
                 index=['n','ME','MAE','RMSE'])
    
    df_metrics.to_csv('metrics_'+name+year+crop_name+'_'+date_YYYYMMDD+'.csv',index=True)
    
    #%% Save for BDB
    
    # headers = ['datum','ruwe sensormeting','modnr','gekalibreerde sensormeting','CImin',
    #             'CImax','veldcapaciteit','stressdrempel','bodemvochtmeting']
    # # modnr=case
    # df_BDB = pd.DataFrame(columns=headers)
    
    # df_BDB['datum']=g_list
    # df_BDB['modnr']=[case]*len(g_list)
    # df_BDB['CImin']=SWC_95_min
    # df_BDB['CImax']=SWC_95_max
    # fc=ParSetMax[6]
    # df_BDB['veldcapaciteit']=[fc]*len(g_list)
    # df_BDB['stressdrempel']=vkrit_list
    
    # list_samp = []
    # for g in g_list:
    #     find = df_obs['Date']==g
    #     if sum(find)==1:
    #         if opkomst: samp30 = df_obs['Mean5'][find]
    #         else: samp30 = df_obs['Mean30'][find]
    #         list_samp.append(samp30)
    #     else: list_samp.append(np.nan)
    # df_BDB['bodemvochtmeting']=list_samp
    
    # df_BDB['ruwe sensormeting']=[np.nan]*len(g_list)
    # df_BDB['gekalibreerde sensormeting']=sensor
    
    # df_BDB.to_csv('df_BDB_'+name+year+crop_name+'_'+date_YYYYMMDD+'.csv',index=False)
    
    
    #%% Risk
    
    os.chdir(save_path)
    sw_df=pd.read_csv('sw_df_'+date_YYYYMMDD+'.csv')
    sw_vkrit_df=pd.read_csv('sw_krit_df_'+date_YYYYMMDD+'.csv')
    g_list=pd.read_csv('g_list_'+name+year+crop_name+'_'+date_YYYYMMDD+'.csv')
    
    risk=30 #30%
    
    df_summary_krit = pd.DataFrame(index=['N(krit)', 'P(krit)','RiskReached('+str(risk)+'%)']) # aantal (N) en kans (P) kritisch
    for h in range(11): # 10 dagen voorspelling + laatste gekende dag
        i=11-h
        krit = np.where(sw_df[sw_df.columns[-i]]<sw_vkrit_df[sw_vkrit_df.columns[-i]])
        krit_n = len(krit[0])
        if krit_n>=0:
            krit_rel = krit_n/sw_df.shape[0]
        if krit_rel>risk/100: below=1
        else: below=0
        df_summary_krit[g_list.iloc[-i][0]]=[krit_n, krit_rel,below]
    print(df_summary_krit)
    
    if simirr: 
        df_summary_krit.to_csv('Critical_simirr_'+date_YYYYMMDD+'.csv',index=True)
    else: 
        df_summary_krit.to_csv('Critical_'+date_YYYYMMDD+'.csv',index=True)
    
    
    ### Critical soil water over the whole period [NEW]
    df_summary_krit_all = pd.DataFrame(index=['N(krit)', 'P(krit)','RiskReached('+str(risk)+'%)']) # aantal (N) en kans (P) kritisch
    for d in range(len(g_list)):
        i=len(g_list)-d
        krit = np.where(sw_df[sw_df.columns[-i]]<sw_vkrit_df[sw_vkrit_df.columns[-i]])
        krit_n = len(krit[0])
        if krit_n>=0:
            krit_rel = krit_n/sw_df.shape[0]
        if krit_rel>risk/100: below=1
        else: below=0
        df_summary_krit_all[g_list.iloc[-i][0]]=[krit_n, krit_rel,below]
    print(df_summary_krit_all)
    
    if simirr: 
        df_summary_krit_all.to_csv('Critical_all_simirr_'+date_YYYYMMDD+'.csv',index=True)
    else: 
        df_summary_krit_all.to_csv('Critical_all_'+date_YYYYMMDD+'.csv',index=True)
    
    # +bekijk ook dikke lijn (most likely) alleen zonder onzekerheid
    
    
    #%% Irrigation advice
    irr_adv_on=False
    if irr_adv_on:
        today = int(ConvertToSerialDate(datetime.today()))
        for day in np.arange(5)+today:
            
            irr_adv, forecast, choice = irr_advice(g_list, sw_list, sw_vkrit_list, 
                                                 swfc_list, R_list, I_list, Ieff, 
                                                 forecast_cf, df_summary_krit, risk=30, 
                                                 day = day)
        
            if choice=='Simuleer': 
                SWC, sw_list, g_list, sensor, covar, df_obs, df_all, df_parameters =SWB(ParSetMax[0],ParSetMax[1],ParSetMax[2], \
                     ParSetMax[3],ParSetMax[4],ParSetMax[5],ParSetMax[6],ParSetMax[7],ParSetMax[8],ParSetMax[9],ParSetMax[10],ParSetMax[11], \
                     [ParSetMax[x] for x in range(12,12+n_mult)],multiplier_bounds, \
                     sensor=True,cal=cal,sensor_cal=sensor_cal,CI=CI,show=[True,''],case=name,year=year,forecast=forecast_cf,df_list=df_list)
    
    
    return df_metrics

#%% Run function

final_func(year,df_list,folder,folder2,cal_par_on,'gen',forecast_R_pf,forecast_ETo_pf,
                multiplier_bounds,simirr=simirr,sim_NoMultiplier=sim_NoMultiplier)


