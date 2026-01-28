# -*- coding: utf-8 -*-
''' 
Soil Water Balance Model based on FAO approaches 

Creator: Marit Hendrickx (KU Leuven)
Last update: September 2025

In this script, there are several #TODO's indicated, where you can change the settings

'''

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)

from cycler import cycler
# Color scheme: color blind & print friendly: 
CB_2_S = ['#d73027', '#fc8d59', '#fee090', '#e0f3f8', '#91bfdb', '#4575b4']
mpl.rcParams['axes.prop_cycle'] = cycler(color=CB_2_S)

from math import *
import numpy as np
from datetime import datetime
# from datetime import date
import pandas as pd

import statistics as stat
import csv
import scipy.stats as stats


#### SPECIFY DIRECTORY ####
import os
main_dir = os.path.dirname(os.path.realpath(__file__))
os.chdir(main_dir)

#### DATE CONVERSION FUNCTIONS ####
def ConvertToSerialDate(date):
    ''' Convert to serial date
    Input: date as datetime(YYYY,MM,DD) '''
    temp = datetime(1899, 12, 30)
    delta = date - temp
    return float(delta.days) + (float(delta.seconds) / 86400)

def ConvertToDate(serial_date):
    ''' Convert from serial date to datetime
    Input: serial date '''
    serial_date=int(serial_date)
    return datetime.fromordinal(datetime(1899, 12, 30).toordinal() + serial_date)


#### SPECIFY SIMULATION YEAR ####
year = '2023'  #TODO


#%% Settings

def settings_func(year):
    drop_samp=np.nan #1 #np.nan #TODO
        # drop a sample: enter location, eg. 1, 2, ..., -1
        # if type(drop_samp)==int: drop starting with that row [drop_samp:]
        # elif type(drop_samp)==list: drop only that row [drop_samp[0]]
        
    validation_days=np.nan #np.nan #20 #10 #TODO
    
    ## Get all sensor-farmer links (year-dependent)
    df_sensor_teler = pd.read_excel(io='Sensors_overview.xlsx', sheet_name='Link sensor-teler')
    df_sensor_teler.sort_values(by=int(year),ignore_index=True,inplace=True)
    df_sensor_teler.dropna(subset=[int(year)], inplace=True)
    df_sensor_teler.reset_index(drop=True, inplace=True) 
    cases = df_sensor_teler['sensor'].tolist()
    
    ## Estimate calibration parameters in DREAM yes or no?
    cal_par_on=False #TODO
    
    ## Estimate sensor error variances in DREAM yes or no?
    lik_sigma_est=False #TODO
    if lik_sigma_est: 
        corr_est=True
    else: 
        corr_est=False
    
    ## Observation data
    sensor=True  #TODO
    DREAM_obsdata='Sensor+stalen'  #'Sensor+stalen' 'Sensor only' 'Samples only' 'Sensor+2samp' #TODO
    if cal_par_on: cal= ''
    else: cal='gen' #'WTLS slope1' or 'WTLS slope_gen'  #TODO
        # WTLS: intc0 or slope1 or slope_gen = general slope possible (or w/o)
        # cal= '' (no cal), 'WTLS', 'gen' (=general)
    
    version = 'DREAM_auto_2025'
    
    ## Weather forecast
    forecast_on=True #TODO
    if not isnan(validation_days):
        forecast_on=False
    
    ## REALTIME APPLICATION?
    realtime=True #TODO
    
    # CHOOSE TODAY AS DATE (realtime application)
    if realtime: date_YYYYMMDD = datetime.strftime(datetime.today(), '%Y%m%d') #TODO
    else:
        # CHOOSE DIFFERENT DATE
        date_YYYYMMDD = year+'1002' #TODO
            # date_YYYYMMDD = '20230625'
            # date_YYYYMMDD='20231002'
    
        ## PLV Herent 2023
        # g0 = 45082
        ## PCG 2022
        # g0 = 44692
        # cal_days = 60 #TODO
        # g_c = g0+cal_days
        # date_YYYYMMDD = ConvertToDate(g_c).strftime('%Y%m%d')
        # print(date_YYYYMMDD)
        
    return drop_samp,cases,df_sensor_teler,cal_par_on,sensor,DREAM_obsdata,lik_sigma_est,corr_est,cal,version,forecast_on,date_YYYYMMDD,validation_days

drop_samp,_,df_sensor_teler,cal_par_on,sensor,_,_,_,cal,_,forecast_on,date_YYYYMMDD,validation_days = settings_func(year)

#%% INITIALIZATION
from rosetta import rosetta, SoilData

def pF_convert_wc_Ksat(x,soil):
    wcr=soil[0][0]
    wcs=soil[0][1]
    alpha=10**soil[0][2]
    n=10**soil[0][3]
    m = 1-1/n
    Ksat=10**soil[0][4]*10 # mm/dag # rosetta returns log10(Ksat[cm/day])
    h=10**x
    wc = wcr+(wcs-wcr)/(1+(alpha*h)**n)**m
    return wc, Ksat

def get_wc_Ksat_from_pF(pF,soiltype):
    if soiltype=='Z':
        sa=90; si=6; cl=4
    elif soiltype=='S':
        sa=77; si=15; cl=8    
    elif soiltype=='P':
        sa=60; si=32; cl=8
    elif soiltype=='L':
        sa=30; si=60; cl=10
    elif soiltype=='A':
        sa=8; si=77; cl=15
    fractions=[sa,si,cl]
    mean, stdev, codes = rosetta(3, SoilData.from_array([fractions]))
    wc, Ksat = pF_convert_wc_Ksat(pF,mean)
    return wc, Ksat


def initial_func(case,df1,df2,forecast_on,show=False):
    try: 
        teler = df_sensor_teler[df_sensor_teler.sensor==case][int(year)].values[0]
        verantw = df_sensor_teler[df_sensor_teler.sensor==case]['Verantwoordelijke'].values[0]
    except: 
        teler = ''
        verantw = ''

    soil=df1.loc[df1['case'] == case]
    p1 = soil['Planting1'].values[0]
    p2 = soil['Planting2'].values[0]
    soil_type=soil['Textuuranalyse_BE'].values[0] # voorlaatste letter
    if isinstance(soil_type, str): # indien geen analysedata: van bodemkaart
        soil_type=soil_type[-2]
    else:
        soil_type=soil['Bodemtype'].values[0][0] # eerste letter
    
    part='0'
    opkomst=False
    if isnan(p2)==False:
        part='2' # 1 or 2 #TODO
        crop_name=soil['Teelt1'].values[0]
        if crop_name=='witloof': part='2' # to skip witloofopkomst #TODO
    
    if part=='0' or part=='1': 
        crop_name=soil['Teelt1'].values[0]
        if crop_name=='witloof' and isnan(p2)==False: part='1'
        elif crop_name=='witloof' and isnan(p2): part='0'
    elif part=='2': crop_name=soil['Teelt2'].values[0]
    
    if case=='C58157' and not year=='2021': part='0' ; opkomst=False
    
    if part=='1' and crop_name=='witloof': 
        opkomst=True
        crop_name = 'witloof_opkomst'
    
    if show:
        print(case)
        print(teler)
        print(soil_type)
        print('Part',part)    
        
    ## Irrigation effiency
    irr_method = soil['Irrigatiemethode'].values[0]
    global Ieff
    if irr_method=='druppel':
        Ieff=0.85
    elif irr_method=='sprinklers' or \
        irr_method=='haspel' or \
            irr_method=='boom':
        Ieff=0.75
    else: Ieff=0.80
    ## Specific irrigation application
    # Ieff=0.95 #TODO
    
    if show:
        print(crop_name)
        print('Irrigation method: '+str(irr_method)+' ('+str(Ieff*100)+'% efficient)')

    ## Crop data
    crop=df2.loc[df2['crop'] == crop_name]
    Zrmax=crop['Zr(max)'].values[0]/100
    fc=soil[['pF2_1','pF2_2','pF2_3']].mean(axis=1,skipna=True,numeric_only=True).values[0] 
    if isnan(fc): 
        try: fc = get_wc_Ksat_from_pF(2,soil_type)[0]
        except: fc = 0.4
        
    K1=crop['Kc ini'].values[0]
    K2=crop['Kcmid2'].values[0]
    K3=crop['Kc end'].values[0]
    Lini=crop['Lini'].values[0]
    if opkomst: 
        if not isnan(p2): Lini=p2-p1
        K1=0
    Ldev=crop['Ldev'].values[0]
    Lmid=crop[['Lmid1','Lmid2','Lmid3']].sum(axis=1,skipna=True).values[0]
    v0=fc*0.8

    if show:
        print('Initial parameters:')
        print('Kc ini = ',K1)
        print('Kc mid = ',K2)
        print('Kc end = ',K3)
        print('Zrmax = ',Zrmax)
        print('L ini = ',Lini)
        print('L dev = ',Ldev)
        print('L mid = ',Lmid)
        print('fc = ',fc)
        print('v0 = ',v0)

    ## Ksat
    try: 
        K_sat_log = log(get_wc_Ksat_from_pF(0,soil_type)[1])
    except: K_sat_log = np.nan()
    if isnan(K_sat_log): # if no data, then use rosetta PFT estimate
        if soil_type=='Z' or soil_type=='S' or soil_type=='P': # Sandy ; Ksat 200-2000
            K_sat_log=log(200)
        elif soil_type=='L' or soil_type=='A': # Loamy ; Ksat 100-750
            K_sat_log=log(100)
        elif soil_type=='E' or soil_type=='U': # Clayey ; Ksat 5-150
            K_sat_log=log(5)

    ## Curve number
    CN = 75
    ## Groundwater depth
    gnj = 150
    
    if show:
        print('log(Ksat) = ',K_sat_log)
        print('Curve number = ',CN)
        print('GWT najaar = ',gnj)
    
    ini=[K1,K2,K3,Lini,Ldev,Lmid,fc,K_sat_log,CN,gnj,Zrmax,v0]

    ## Weather forecasts
    loc_weather = teler
    if forecast_on:
        forecast_R=pd.read_csv(year+'/weather_forecasts/'+case+'/cf_'+case+'_'+date_YYYYMMDD+'_daily_R.csv', encoding= 'unicode_escape')
        forecast_R.loc[~(forecast_R['Daily precipitation [mm/day]'] > 0), 'Daily precipitation [mm/day]']=0
        
        forecast_ETo=pd.read_csv(year+'/weather_forecasts/'+case+'/cf_'+case+'_'+date_YYYYMMDD+'_daily_ETo.csv', encoding= 'unicode_escape')
        
        forecast = pd.DataFrame({'Date':forecast_R['Date'],'R':forecast_R['Daily precipitation [mm/day]'],'ETo':forecast_ETo['Daily ETo [mm/day]']})
        
        if show: print(forecast) # in mm
        forecast['Date'] = [int(ConvertToSerialDate(datetime.strptime(date, '%d-%m-%Y'))) \
                               for date in forecast['Date']]
            
        forecast['I']=0 # set every forecast day irrigation amount to zero
        
        ## ADD IRRIGATION SIMULATION #TODO
            # example:
            # forecast.loc[0,'I']=10
            # forecast.loc[4,'I']=15
        forecast_date = date_YYYYMMDD
        os.chdir(main_dir)
        
    else: 
        forecast=np.empty(0)
        forecast_date = np.nan
        
    if part=='2':
        g0=soil['Planting2'].values[0]
    else: g0=soil['Planting1'].values[0]
    g0=int(g0)
        
    return ini, teler, opkomst, verantw, p1, p2, part, soil_type, crop_name, irr_method, forecast, forecast_date, g0

#%% FUNCTIONS
from Sensordata import sensordata as sd_func

## Convert sensordata (METER equation)
def sensorEQ(mV):
    ''' Equation to calculate the SWC from the mV sensor measurements '''
    SWC=4.824*10**(-10)*mV**3-2.278*10**(-6)*mV**2+3.898*10**(-3)*mV-2.154
    return SWC

## Crop growth
def kcb_func(kcb_ini, kcb_mid, kcb_end, lini, ldev, lmid, llate, ltot):
    kcb = np.empty(ltot + 1)
    kcb[:lini + 1] = kcb_ini
    kcb[lini + 1:lini + ldev + 1] = np.linspace(kcb_ini, kcb_mid, ldev)
    kcb[lini + ldev + 1:lini + ldev + lmid + 1] = kcb_mid
    kcb[lini + ldev + lmid + 1:] = np.linspace(kcb_mid, kcb_end, ltot - lini - ldev - lmid)
    return kcb

## Crop root growth
def zr_func(zrmin, zrmax, lzini, lzdev, lzmid, lzlate, lztot):
    zr = np.empty(lztot + 1)
    zr[:lzini + 1] = zrmin
    zr[lzini + 1:lzini + lzdev + 1] = np.linspace(zrmin, zrmax, lzdev)
    zr[lzini + lzdev + 1:] = zrmax
    return zr

## Kcmax
def kcmax_func(kcmaxini, kcmaxmid, lini, ldev, lmid, llate, ltot):
    kcmax = np.empty(ltot + 1)
    kcmax[:lini + 1] = kcmaxini
    kcmax[lini + 1:lini + ldev + 1] = np.linspace(kcmaxini, kcmaxmid, ldev)
    kcmax[lini + ldev + 1:] = kcmaxmid
    return kcmax

## Aggregations of samples
def agg_func(date_Obs,obs):
    date_obs_grouped=[date_Obs[0]]
    indices0 = [j for j, x in enumerate(date_Obs) if x == date_Obs[0]]
    obs_list0=obs[indices0[0]:indices0[len(indices0)-1]+1]

    if len(obs_list0)==1:
        mean=[obs_list0[0]]
        stdev=[0]
    else:
        mean=[stat.mean(obs_list0)]
        stdev=[stat.stdev(obs_list0)]
        
    for i in range(1,len(date_Obs)):
        if date_Obs[i]!=date_Obs[i-1]:
            date_obs_grouped.append(date_Obs[i])
            indices = [j for j, x in enumerate(date_Obs) if x == date_Obs[i]]
            obs_list=obs[indices[0]:indices[len(indices)-1]+1]
            if len(obs_list)==1:
                mean.append(obs_list[0])
                stdev.append(0)
            else:
                mean.append(stat.mean(obs_list))
                stdev.append(stat.stdev(obs_list))        
        
    return mean, stdev, date_obs_grouped


#%% SWB FUNCTION

I_list_adv=np.array([])

def SWB(K1,K2,K3,Lini,Ldev,Lmid,fc,K_sat_log,CN,gnj,Zrmax,v0,sensor,cal,sensor_cal,CI,show,case,year,forecast,df_list):
    ''' Soil water balance model 
    Input: crop coefficients (3), length growth stages (3),
    field capacity (m³/m³), log Ksat, curve number (runoff),
    maximum groundwater table (cm), max rooting depth, 
    initial soil water content (m³/m³)
    
    Sensordata, with or without calibration,
    simulated SWC and confidence interval that is calculated in DREAM (otherwise put np.empty(0)),
    if show=True: plots are shown,
    name and year to define the case
    
    forecast : columns: Date-R-ETo
    
    Output:
        Simulates soil water at the end of the day
        
    '''
    #runtime start
    start=datetime.now()
    
    date_axis=True # x axis is date instead of growing day
    show_vec = show # [True,folder] or [False, ''] or ..
    if show_vec[1]=='':
        save=False
    else: 
        save=True
        folder=show_vec[1]
    show=show_vec[0] 

    ## dataframes
    [df_R,df_I,df_ET,df1,df_obs,df2] = df_list
    
    ####  Read in climate data (CSV) ------------------------------------------
    today=np.floor(ConvertToSerialDate(datetime.strptime(date_YYYYMMDD,'%Y%m%d')))
    
    df_forecast=forecast
    if not df_forecast.size: ## IF END-OF-CYCLE (no prediction period)
        if year=='2023':
            if case=='C54992': today = np.floor(ConvertToSerialDate(datetime.strptime('20230828','%Y%m%d'))) # selder oogst 28/8
            elif case=='C58157': today = np.floor(ConvertToSerialDate(datetime.strptime('20231030','%Y%m%d'))) # witloof rooi 30/10
            elif case=='C54986': today = np.floor(ConvertToSerialDate(datetime.strptime('20230929','%Y%m%d'))) # prei oogst 29/9
            elif case=='C562C6': today = np.floor(ConvertToSerialDate(datetime.strptime('20231006','%Y%m%d'))) # bataat oogst 6/10
            elif case=='C535C9': today = np.floor(ConvertToSerialDate(datetime.strptime('20230920','%Y%m%d'))) # selder goede data tot 20/9
        elif year=='2022':
            if case=='C58157': today = np.floor(ConvertToSerialDate(datetime.strptime('20220912','%Y%m%d'))) # PLV witloof rond 12/9
            if case=='C54986': today = np.floor(ConvertToSerialDate(datetime.strptime('20220817','%Y%m%d'))) # PCG prei oogst 16/8
        elif year=='2021':
            if case=='C54762': today = np.floor(ConvertToSerialDate(datetime.strptime('20210827','%Y%m%d'))) # einde sensordata rond 27/8

        
    df_R = df_R.loc[df_R['Date']<today]
    if case=='C54BB2': rain = [0 if isnan(x) else x for x in df_R['C58157']]
    else: rain = [0 if isnan(x) else x for x in df_R[case]]
    
    df_I = df_I.loc[df_I['Date']<today]
    irri=np.zeros(len(df_R['Date'])).tolist()
    
    irr_case = case
    # irr_case = 'C54BB2' #TODO
    ## SPECIFIC IRRIGATION TREATMENTS
        ## 2023
        # C562C6_telersGevoel & C562C6_beperkteIrr bataat PCG
        # C54986_telersGevoel prei PCG
        # 'C535C9' PSKW sensorgestuurd
        # 'C54BB2' PLV 0% drip na opkomst
        # 'C58157_0.5naopk' PLV 50% drip na opkomst
        ## 2022
        # 'C54BB2' PLV 0% drip na opkomst

    for i in range(len(df_R['Date'])):
        if df_R['Date'][i] in df_I['Date'].values:
            if isnan(df_I[irr_case][df_I['Date']==df_R['Date'][i]])==False:
                irri[i]=df_I[irr_case][df_I['Date']==df_R['Date'][i]].tolist()[0]

    date_R_I = df_R['Date'].tolist()

    ## Possibility to multiply irrigation by a factor
    fractie=1 #TODO
    irri=[i * fractie for i in irri]
    
    ## If irrigation during 'witloofopkomst' should remain unaffected by the applied factor:
    # na_index = np.where(df_R['Date']==45098)[0][0]
    # for i in range(na_index,len(df_R['Date'])):
    #     irri[i] = irri[i]*fractie
    

    if case=='C54BB2': df_ET = df_ET[['Date','C58157']]
    else: df_ET = df_ET[['Date',case]]
    df_ET = df_ET.loc[df_ET['Date']<today]
    if case=='C54BB2': eto = df_ET['C58157'].tolist()
    else: eto = df_ET[case].tolist()
    date_ETo = df_ET['Date'].tolist()
    if show: print('Last available weather data on', max(date_ETo))
    
    if df_forecast.size:
        if show: print('Forecast:\n',df_forecast)
        forecast_date=min(df_forecast['Date'])-0.5
        date_R_I.extend(df_forecast['Date'].tolist())
        date_ETo.extend(df_forecast['Date'].tolist())
        rain.extend(df_forecast['R'].tolist())
        irri.extend(df_forecast['I'].tolist())
        eto.extend(df_forecast['ETo'].tolist())
        
        
    ####  Get soil and crop data (CSV) ----------------------------------------

    ## SOIL DATA
    soil=df1.loc[df1['case'] == case]#.loc[df1['Diepte'] == '0-30']
    
    # SOIL OBS
    df_obs=df_obs.loc[df_obs['Sensornr'] == case].iloc[:,:4]  ## 'Date', '0_30(grav%)', '30_60(grav%)', '0_5(grav%)'
    df_obs=df_obs.loc[df_obs['Date'] <= today]
    df_obs.reset_index(inplace=True,drop=True)
    
    ## CROP DATA
    _, _, _, _, p1, p2, part, soil_type, crop_name, irr_method, _, _, _ = initial_func(case,df1,df2,forecast_on=False)
    crop=df2.loc[df2['crop'] == crop_name]

    ####  Growing period ------------------------------------------------------

    ## Planting date g0
    
    #as date
    # g0 = datetime(2020,7,11)
    # print(g0)
    # g0=ConvertToSerialDate(g0)
    # print('g0 =',int(g0))

    
    #or as serial number
    if part=='2':
        g0=soil['Planting2'].values[0]
    else: g0=soil['Planting1'].values[0]
    # if case=='C560EE' and part=='2': g0=44357
    g0=int(g0)
    startdate=max(min(date_R_I),min(date_ETo),g0)
    if show:
        print('g0 =',g0)
        print('Planting date =',ConvertToDate(g0))
        print('First available data on',startdate)
    g0=startdate
    

    ####  Crop growth ---------------------------------------------------------
    
    ## Growth phases (calibrated)
    lini=round(Lini)
    ldev=round(Ldev)
    lmid=round(Lmid)
    llate=crop['Lend'].astype(int).values[0]
    ltot_min=int(crop['Lini'].values[0]+crop['Ldev'].values[0]+ \
        crop[['Lmid1','Lmid2','Lmid3']].sum(axis=1,skipna=True).values[0]+llate)
    if crop_name=='witloof_opkomst':
        ltot=lini
        K1=0
        llate=0
    else: 
        ltot=ltot_min
        llate=ltot_min-lmid-ldev-lini
        if llate<0:lmid=ltot_min-ldev-lini

    enddate=g0+ltot+1
    if part=='1' and isnan(p2)==False:
        enddate=min(enddate,p2-1)
    # if today is the last day to simulate
    # enddate=floor(ConvertToSerialDate(datetime.now()))
    if show: print('End growth =',enddate)
    
    ## Function describing evolution in rooting depth (in function of days after planting)
    lzini=lini
    lzdev=ldev+floor(lmid/2)
    lzmid=ceil(lmid/2)
    lzlate=llate
    lztot=ltot
    ## Crop rooting parameters (minimal and maximal rooting depth) (m)
    zrmin=crop['Zr(ini)'].values[0]/100
    if zrmin<0.3: zrmin=0.3
    zrmax=Zrmax
    # zr list from function    
    zr = zr_func(zrmin,zrmax,lzini,lzdev,lzmid,lzlate,lztot)
        
    if show or save:
        plt.figure()
        plt.subplot(211)
        plt.plot(zr)
        plt.ylabel('Rooting depth (m)')
        ax = plt.gca()
        ax.set(xlim=(0, ltot))
    

    ####  K factors -----------------------------------------------------------
    
    kcb_ini=K1
    kcb_mid=K2
    kcb_end=K3
    # kcb list from function    
    kcb = kcb_func(kcb_ini,kcb_mid,kcb_end,lini,ldev,lmid,llate,ltot)    
    
    if show or save:
        plt.subplot(212)
        plt.plot(kcb)
        plt.ylabel('Kcb')
        plt.xlabel('Growing days')
        ax = plt.gca()
        ax.set(xlim=(0, ltot))
        if show and not save: plt.show()
        elif save: 
            plt.savefig(folder+'\\growth.png',dpi=300)
            plt.close()
    
    ## Function describing evolution Kcmax
    # kcmaxini=1 
    kcmaxini=crop['Kcmaxini'].astype(float).values[0]
    # kcmaxmid=1.1
    kcmaxmid=crop['Kcmaxmid'].astype(float).values[0]
    # kcmax list from function  
    kcmax=kcmax_func(kcmaxini,kcmaxmid,lini,ldev,lmid,llate,ltot)
    
    ####  Soil parameters -----------------------------------------------------
    
    # saturation
    sat=soil[['pF0_1','pF0_2','pF0_3']].mean(axis=1,skipna=True,numeric_only=True).values[0] 
 
    # field capacity - is estimated in DREAM
    # fc=soil[['pF2_1','pF2_2','pF2_3']].mean(axis=1,skipna=True).values[0] 
    # if isnan(fc): fc=fc_pf         
    
    pF2=soil[['pF2_1','pF2_2','pF2_3']].mean(axis=1,skipna=True,numeric_only=True).values[0]
    
    # wilting point
    wp=soil[['pF4.2_1','pF4.2_2','pF4.2_3']].mean(axis=1,skipna=True,numeric_only=True).values[0]
    if isnan(sat) and isnan(pF2) and isnan(wp): 
        wc_from_pF = get_wc_Ksat_from_pF(np.array([0,2,4.2]),soil_type)[0]
        sat = wc_from_pF[0]
        pF2 = wc_from_pF[1]
        wp = wc_from_pF[2]
    elif isnan(wp) and not isnan(sat): 
        wp = get_wc_Ksat_from_pF(np.array([4.2]),soil_type)[0][0]

    if fc>=pF2: sat=min(0.5,sat*fc/pF2)
    # bulk density
    bd=soil[['bd_1','bd_2','bd_3']].mean(axis=1,skipna=True,numeric_only=True).values[0]
    if isnan(bd): bd = (1-sat)*2.65 # Estimate BD based on saturation (TPV) and specific weight of soil

    ## initial SWC is estimated in DREAM
    # initial soil moisture at g0
    # v0=soil['v0'].values[0]
    # if isnan(v0): v0=fc*0.8   # dit is random gekozen
    
    opkomst=False
    ## Is there an initial soil sample?
    df_obs=df_obs[df_obs['Date']>=g0].reset_index(drop=True)
    # raw soil sample data (df_obs loaded earlier)
    obs0_5_date=np.array(df_obs['Date'])
    obs0_5=df_obs['0_5(grav%)']
    df_obs=df_obs.drop(['0_5(grav%)'],axis=1)
    # df_obs=df_obs.dropna().reset_index(drop=True)
    df_obs=df_obs.dropna(thresh=2).reset_index(drop=True)
    # date_Obs=[x for x in df_obs['Date']]
    date_Obs=np.array(df_obs['Date'])
    #data was given as percentage gravimetric
    obs0_30=[x*bd/100 for x in df_obs['0_30(grav%)']]
    obs30_60=[x*bd/100 for x in df_obs['30_60(grav%)']]
    # available soil water                      
    obs0_30_asw=[x*300-wp*300 for x in obs0_30]
    obs30_60_asw=[x*300-wp*300 for x in obs30_60]
        
    if df_obs.size:
        ind0_5 = obs0_5[~np.isnan(obs0_5)].index.tolist()
        if len(ind0_5)>0 and crop_name=='witloof_opkomst':
            opkomst=True
            obs0_5 = [obs0_5[ind0_5[i]] for i in range(len(ind0_5))]
            obs0_5=[x*bd/100 for x in obs0_5]
            obs0_5_asw=[x*50-wp*50 for x in obs0_5]
            obs0_5_date=[obs0_5_date[ind0_5[i]] for i in range(len(ind0_5))]
        
        if not opkomst:
            if np.any(date_Obs==g0): 
                ini_i = np.where([date_Obs[i]==g0 for i in range(len(date_Obs))])[0]
                ini_v = [obs0_30[ini_i[i]] for i in range(len(ini_i))]
                # v0 = np.mean(ini_v)
            elif np.any(date_Obs==g0+1): 
                ini_i = np.where([date_Obs[i]==g0+1 for i in range(len(date_Obs))])[0]
                ini_v = [obs0_30[ini_i[i]] for i in range(len(ini_i))]
                # v0 = np.mean(ini_v)
            elif np.any(date_Obs==g0-1): 
                ini_i = np.where([date_Obs[i]==g0-1 for i in range(len(date_Obs))])[0]
                ini_v = [obs0_30[ini_i[i]] for i in range(len(ini_i))]
                # v0 = np.mean(ini_v)
            elif np.any(date_Obs==g0+2): 
                ini_i = np.where([date_Obs[i]==g0+2 for i in range(len(date_Obs))])[0]
                ini_v = [obs0_30[ini_i[i]] for i in range(len(ini_i))]
                # v0 = np.mean(ini_v)
            elif np.any(date_Obs==g0-2): 
                ini_i = np.where([date_Obs[i]==g0-2 for i in range(len(date_Obs))])[0]
                ini_v = [obs0_30[ini_i[i]] for i in range(len(ini_i))]
                # v0 = np.mean(ini_v)
            elif np.any(date_Obs==g0+3): 
                ini_i = np.where([date_Obs[i]==g0+3 for i in range(len(date_Obs))])[0]
                ini_v = [obs0_30[ini_i[i]] for i in range(len(ini_i))]
                # v0 = np.mean(ini_v)
            elif np.any(date_Obs==g0-3): 
                ini_i = np.where([date_Obs[i]==g0-3 for i in range(len(date_Obs))])[0]
                ini_v = [obs0_30[ini_i[i]] for i in range(len(ini_i))]
                # v0 = np.mean(ini_v)
            else: 
                if show: print('No g0 0-30cm sample')
    else: 
        if show: print('No samples')
        
    
    v0_0_5=0
    if opkomst:
        if np.any(obs0_5_date==g0): 
            ini_i = np.where([obs0_5_date[i]==g0 for i in range(len(obs0_5_date))])[0]
            ini_v = [obs0_5[ini_i[i]] for i in range(len(ini_i))]
            v0_0_5 = np.mean(ini_v)
        elif np.any(obs0_5_date==g0+1): 
            ini_i = np.where([obs0_5_date[i]==g0+1 for i in range(len(obs0_5_date))])[0]
            ini_v = [obs0_5[ini_i[i]] for i in range(len(ini_i))]
            v0_0_5 = np.mean(ini_v)
        elif np.any(obs0_5_date==g0+2): 
            ini_i = np.where([obs0_5_date[i]==g0+2 for i in range(len(obs0_5_date))])[0]
            ini_v = [obs0_5[ini_i[i]] for i in range(len(ini_i))]
            v0_0_5 = np.mean(ini_v)
        else: 
            if show: print('No g0 0-5cm sample')
        if v0_0_5>0: v0=v0_0_5
    if show: print('Opkomst: ',opkomst)
    
    ## Critical threshold opkomst
    TAW_opkomst = 0.322 # preliminary value
    vkrit_opkomst = (pF2-wp)*TAW_opkomst + wp
    vkrit_opkomst = (fc-wp)*TAW_opkomst + wp
    
    ## Soil type
    soil_type=soil['Textuuranalyse_BE'].values[0] # voorlaatste letter
    if isinstance(soil_type, str): # indien geen analysedata: van bodemkaart
        soil_type=soil_type[-2]
    else:
        soil_type=soil['Bodemtype'].values[0][0] # eerste letter
        
    ## Define TEW & REW (mm)
    Ze=0.1
    TEW=1000*(fc-0.5*wp)*Ze
    # Relation REW and TEW based on FAO table 19 ~ soil type
    if soil_type=='Z' and TEW<=12: REW=2+(TEW-6)*5/6
    elif (soil_type=='S' or soil_type=='P') and TEW<=20 and TEW>=9: REW=4+(TEW-9)*6/11
    elif (soil_type=='L' or soil_type=='A') and TEW<=26 and TEW>=16: REW=8+(TEW-16)*3/10
    elif (soil_type=='E' or soil_type=='U') and TEW<=29 and TEW>=22: REW=8+(TEW-22)*4/7
    else: REW=2+(TEW-6)*10/23
    
    ## INSERT p (FAO table 22 ~ crop)
    pFAO=crop['pFAO'].values[0]
    
    ## INSERT ksat value
    ksat=e**K_sat_log           # mm/day
    
    if show:
        print('Soil type = ',soil_type)
        print('Ksat = ',ksat)
        print('sat = ',sat)
        print('fc = ',fc)
        print('wp = ',wp)
        print('bd = ',bd)
        print('v0 = ',v0)
        print('TEW = ',TEW)
        print('REW = ',REW)  
        print('pFAO = ',pFAO)  
        
    df_pars = pd.DataFrame({'g0':g0,'Soil type':soil_type,'Ksat':ksat,'sat':sat,
                            'fc':fc,'wp':wp,'bd':bd,'v0':v0,'TEW':TEW,'REW':REW,
                            'pFAO':pFAO,'Irr method':irr_method}, index=[0])
    
    ####  Functions --------------------- # Can these be defined beforehand, outside of SWB function?
                       
    ## Define function that expresses capillary rise
    def cr(x):
        # Parameters depending on soil type
        if soil_type=='Z' or soil_type=='S' or soil_type=='P': # Sandy ; Ksat 200-2000
            a=-0.3112-10**(-5)*ksat
            b=-1.4936+0.2416*log(ksat)
        elif soil_type=='L' or soil_type=='A': # Loamy ; Ksat 100-750
            a=-0.4986+9*10**(-5)*ksat
            b=-2.1320+0.4778*log(ksat)
        elif soil_type=='E' or soil_type=='U': # Clayey ; Ksat 5-150
            a=-0.5677-4*10**(-5)*ksat
            b=-3.7189+0.5922*log(ksat)
        # x = z (m)
        # cr (mm/day)
        return exp((log(x)-b)/a)
    
    ## Define function that expresses ground water table depth in function of day in year
    ## INSERT ground water table data
    def depth(x,gnj):
        ''' x = day of the year (can be calculated: date.strftime(%j))
        depth in cm'''
        gvj=100
        # gvj=soil['gvj'].values[0] # possibility to put gvj in the soil file (location specific)
        # DREAM estimation of gnj (range: 100-200?)
        ampl=(gnj-gvj)/2
        return gvj+ampl-sin((x+15)*pi/180)*ampl
    
    
    ## Function day of the year
    def DOY(g):
        g_date=ConvertToDate(g)
        return int(g_date.strftime('%j')) #day of year
    
    
    ## Precipitation efficiency : runoff (curve number method)
    def runoff(R, CN):
        S = 25400/CN - 254      # Pot. max. retention, depending on soil structure (curve number CN)
        Ia = 0.2*S              # Initial abstraction Ia. Actual retention = R - Ia - Q
        if R>Ia : Q = (R - Ia)**2/(R + 0.8*S)      # Accumulated runoff depth Q
        else: Q=0
        return Q, Ia
    
    ## Function expressing soil water content in the 10 cm evaporative layer (profile)
    def evap_layer(v_Ze,v_10,depth):
        ''' Depth (cm)
            v_Ze (m³/m³): soil moisture content of evap layer (10cm)
            v_10, v_0 (m³/m³): soil moisture content at 10 and 0 cm depth '''
        # v_10 = v  # fc for bare soil
        v_0 = 2*v_Ze-v_10    
        v=v_0+(v_10-v_0)/10*depth
        if v<wp*0.05: v=wp*0.05
        if v>fc: v=fc
        if v_0<wp*0.05: v_0=wp*0.05
        if v_0>fc: v_0=fc
        # v at specific depth
        return v,v_0
    
    ####  Irrigation/precipitation efficiency ---------------------------------
    
    Reff=1    
    ## Rainfall efficiency depending on curve number CN
    ## CN is calibrated
    
    ## Wetting fraction
    fw=1
    # Precipitation & sprinkler = 1
    # Drip irrigation = 0.3-0.4?
    few=fw
    # Drip irrigation: few=(1-(2/3)*cov)*fw?
    
    ####  Water balance loop over growth period -------------------------------

    ## Initial conditions
    sw0=v0*zrmin*1000-wp*zrmin*1000   #sw0 is soil water storage at day 0 in mm
    swfc0=fc*1000*zrmin-wp*1000*zrmin
    v=v0
    sw=sw0
    
    if v0_0_5==0: 
        v0_Ze=v0
    else: 
        v0_5_10=v0_0_5/7+6/7*v0
        v0_Ze=(v0_0_5+v0_5_10)/2    
    
    v0_10=(3*v0-v0_Ze)/2 
    v0_5_layer = np.mean([evap_layer(v0_Ze,v0_10,i)[0] for i in np.arange(0,5.1,0.1)])
    De=fc*1000*Ze-v0_Ze*Ze*1000
    if De<0: De=0    
    
    kr_list=[]
    De_list=[]
    v_Ze_list=[]
    sw_Ze_list=[]
    v_5_list=[]
    sw_5_list=[]
    ke_list=[]
    kcb_list=[]
    kc_list=[]
    kc_a_list=[]
    kcb_a_list=[]
    ETcb_list=[]
    ETcb_a_list=[]
    ETe_list=[]
    etm_list=[]
    etm_cum_list=[]
    etcb_cum_list=[]
    etcb_a_cum_list=[]
    ete_cum_list=[]
    eta_cum_list=[]
    vkrit_list=[]
    vkrit_Ze_list=[]
    sw_vkrit_list=[]
    sw_Ze_vkrit_list=[]
    sw_5_vkrit_list=[]
    # global swfc_list
    swfc_list=[]
    swfc_Ze_list=[]
    swfc_5_list=[]
    Dr_list=[]
    TAW_list=[]
    RAW_list=[]
    eta_list=[]
    Ks_list=[]
    p_list=[]
    cri_list=[]
    v_list=[v0]
    sw_list=[sw0]
    zvirt_list=[zr[0]]

    end=int(min(ltot,max(date_ETo)-g0,max(date_R_I)-g0,enddate-g0))+1
    if show: 
        print('End of sim =',end + int(g0))
        print('End of weather data =',int(min(max(date_ETo)-g0,max(date_R_I)-g0)+1 + int(g0)))
        print('End of growth data =',int(ltot+1 + int(g0)))

    g_list=[int(g0+i) for i in range(0,end)]
    etoi=[eto[date_ETo.index(g_list[i])] for i in range(0,end)]
    rain=np.multiply(rain,Reff)
    runoff_list=[runoff(rain[date_R_I.index(g_list[i])],CN)[0]  for i in range(0,end)]
    R_list=[rain[date_R_I.index(g_list[i])] - runoff_list[i]  for i in range(0,end)]
    I_list=[irri[date_R_I.index(g_list[i])]*Ieff for i in range(0,end)]

    if I_list_adv.size: I_list+=I_list_adv*Ieff
    kcmax_value=[max(kcmax[i],kcb[i]+0.05) for i in range(0,end)]

    ## TEW correction for bare soil
    TEW = np.ones(len(g_list))*TEW
    REW = np.ones(len(g_list))*REW
    for i in range(0,min(lini,end)):
        if etoi[i]<5: TEW[i] = TEW[i]*np.sqrt(etoi[i]/5)
        if REW[i]>TEW[i]: REW[i]=TEW[i]
        
    TEW_swc = [fc - TEW[i]/(1000*Ze) for i in range(0,end)]
    REW_swc = [fc - REW[i]/(1000*Ze) for i in range(0,end)]
    
    for i in range(0,end):
        
        ## Virtual soil depth (compensating)
        if i!=0:
            v_yest = min(v,fc)     # percolatie van vorige dag gebeurt pas vandaag (voor balans)
            ## Virtual depth (m) grows with root growth (bottom of layer)
            if v_yest==wp: X = X + (zr[i]-zr[i-1])*(fc)/(v_yest)
            else: X = X + (zr[i]-zr[i-1])*(fc-wp)/(v_yest-wp)
            ## Virtual depth (m) with precipitation/irrigation (top of layer)
            v_new30 = v_yest + (R_list[i]+I_list[i])/300  #mm
            if v_new30 >= fc: 
                Xfc = X * v_yest/fc   # virtual gaat eerst krimpen om fc te bereiken met originele mm water
                perc = (v_new30-fc)*0.3
                Xperc = Xfc + perc/fc     # virtual gaat dan weer groeien om fc te behouden met extra mm percolatiewater
                X = min(Xperc,zr[i]-0.3)
            else:
                X = X * v_yest/v_new30
            zvirt = 0.3+X
            zvirt_list.append(zvirt)    
        else: X = 0
        

        ## Evaporation (10 cm)
        if De<=REW[i]: kr=1
        elif De>=TEW[i]: kr=0
        else: kr=(TEW[i]-De)/(TEW[i]-REW[i])
        kr_list.append(kr)
        ke=min(kr*(kcmax_value[i]-kcb[i]),few*kcmax_value[i])
        ke_list.append(ke)
        ETcb_list.append(kcb[i]*etoi[i])
        ETe_list.append(ke*etoi[i])
    
        ## Calculate ETm
        kc=min((kcb[i]+ke),kcmax_value[i])
        kc_list.append(kc)
        etm=kc*etoi[i]
        etm_list.append(etm)
        
        vkrit=fc-(pFAO+0.04*(5-etm))*(fc-wp)    # RAW = p*TAW
        vkrit_list.append(vkrit)
        sw_vkrit=vkrit*zr[i]*1000-wp*zr[i]*1000
        sw_vkrit_list.append(sw_vkrit)
        swfc=fc*1000*zr[i]-wp*1000*zr[i]
        swfc_list.append(swfc)
        
        # Emergence threshold (opkomst witloof)
        vkrit_Ze=vkrit_opkomst
        vkrit_Ze_list.append(vkrit_Ze)
        sw_Ze_vkrit=vkrit_Ze*Ze*1000-wp*Ze*1000
        sw_Ze_vkrit_list.append(sw_Ze_vkrit)
        sw_5_vkrit=vkrit_Ze*0.05*1000-wp*0.05*1000
        sw_5_vkrit_list.append(sw_5_vkrit)
        swfc_Ze=fc*1000*Ze-wp*1000*Ze
        swfc_Ze_list.append(swfc_Ze)
        swfc_5=fc*1000*0.05-wp*1000*0.05
        swfc_5_list.append(swfc_5)
        
        ## Calculate percolation
        if v>fc: #when v is above field capacity it drains back to fc
            p=sw-swfc #calculates amount of percolation (mm)
            sw=swfc
        else: p=0
        p_list.append(p)

        ## Calculate eta
        Dr = swfc-sw # root zone depletion Dr
        TAW = swfc
        RAW = swfc-sw_vkrit
        Dr_list.append(Dr)
        TAW_list.append(TAW)
        RAW_list.append(RAW)
        ## Water stress coefficient Ks
        if Dr<=RAW: Ks=1
        elif Dr>=TAW: Ks=0
        else: Ks=(TAW-Dr)/(TAW-RAW)
        
        ## Calculate ETm_adj
        kcb_list.append(kcb[i])
        kc_a=min((Ks*kcb[i]+ke),kcmax_value[i])
        kc_a_list.append(kc_a)
        kcb_a_list.append(Ks*kcb[i])
        eta=kc_a*etoi[i]
        Ks_list.append(Ks) 
        eta_list.append(eta) 
        ETcb_a_list.append(Ks*kcb[i]*etoi[i])
        
        ## Cumulative ET
        if i==0 or i==lini or i==lini+ldev or i == lini+ldev+lmid or i==ltot:
            etm_cum=0
            etcb_cum=0
            etcb_a_cum=0
            ete_cum=0
            eta_cum=0
        else:
            etm_cum+=etm
            etcb_cum+=kcb[i]*etoi[i]
            etcb_a_cum+=Ks*kcb[i]*etoi[i]
            ete_cum+=ke*etoi[i]
            eta_cum+=eta
        etm_cum_list.append(etm_cum)
        etcb_cum_list.append(etcb_cum)
        etcb_a_cum_list.append(etcb_a_cum)
        ete_cum_list.append(ete_cum)
        eta_cum_list.append(eta_cum)
        
        ## Calculates capillary rise, only occuring when v is below fc when above cr is considered 0
        z=depth(DOY(g_list[i]),gnj)/100 #z(m)
        if v<fc:
            cri=cr(z)
        else:cri=0
        cri_list.append(cri)
    
        ## Calculate soil water balance (new sw value)
        if i!=0:
            sw=sw+R_list[i]-eta+cri+I_list[i]   # p is already accounted for (sw = swfc)
            sw=sw+(zr[i]-zr[i-1])*1000*(fc-wp)  # soil underneath root zone is at fc, so extra available water
            if sw<0:
                sw=0
            sw_list.append(sw)
            v=sw/(zvirt*1000)+wp
            v_list.append(v)
            
        
        ## Evaporative 10 cm layer
        # Kr en dus Ke worden berekend adhv de depletie De van de vorige dag (i-1);
        # depletie De van een dag i-1 wordt op zijn beurt berekend adhv de Ke van dag i-1 zelf
        if i!=0:
            if De<0: De=0 # fc  # added (12/6/2023)
            elif De>TEW[i]: De=TEW[i] # added (12/6/2023)
            De+=-(Ze/zr[i]*R_list[i]+Ze/zr[i]*I_list[i]/fw) # Ze/zr[i]* of niet?
            if De<0: De=0 # fc
            elif De>TEW[i]: De=TEW[i]
            De += ke*etoi[i] + Ze/zr[i]*kcb[i-1]*etoi[i] + Ze/zr[i]*p
            # De increases with (part of) ETa and percolation p
        De_list.append(De)
        v_Ze=fc-De/(1000*Ze)
        v_Ze_list.append(v_Ze)
        sw_Ze=v_Ze*Ze*1000-wp*Ze*1000
        sw_Ze_list.append(sw_Ze)
        
        ## 5 cm layer
        v_10=(3*v-v_Ze)/2 
        v_5_layer = np.mean([evap_layer(v_Ze,v_10,i)[0] for i in np.arange(0,5.1,0.1)])
        v_5_list.append(v_5_layer)
        sw_5_layer=v_5_layer*0.05*1000-wp*0.05*1000
        sw_5_list.append(sw_5_layer)
        

    
    #### Soil water content ---------------------------------------------------
    
    # Soil water content in the root zone (assumed to equal swc in upper 30cm; measured by sensors)
    SWC = v_list # is hetzelfde
    
    # SWC_60 = [(zr[i]*1000-300)*SWC[i]/300 + (600 - zr[i]*1000)*fc/300 for i in range(end)]
    # Soil water content in the 30-60 cm soil layer
    # Assuming field capacity under the root zone eg. 30-40 on SWC, 40-60 on fc
        
    #### Soil samples ---------------------------------------------------------

    ## Calculate available soil moisture over entire root zone
    for i in range(0,len(date_Obs)):
        if obs0_30_asw[i]<0:
            obs0_30_asw[i]=0
        if obs30_60_asw[i]<0:
            obs30_60_asw[i]=0

    # obs_asw=[obs0_30_asw[i]+obs30_60_asw[i]*(zr[date_Obs[i]-g0]-0.3)/0.3 for i in range(0,len(date_Obs))] # zrmin is always 0.3
    obs_asw=[0]*len(date_Obs)
    for i in range(0,len(date_Obs)):
        if date_Obs[i]-g0 <= len(zr)-1:
            if zr[date_Obs[i]-g0]>0.3:
                if isnan(obs30_60_asw[i]):
                    obs_asw[i]=obs0_30_asw[i]
                else:
                    obs_asw[i]=obs0_30_asw[i]+obs30_60_asw[i]*(zr[date_Obs[i]-g0]-0.3)/0.3
            else: obs_asw[i]=obs0_30_asw[i]
        else: obs_asw[i]=np.nan
        
    #aggregate and calculate means and stdevs
    if df_obs.size: # if there are soil sample data
        ## ASW
        obs_asw_mean, obs_asw_stdev, date_obs_grouped = agg_func(date_Obs,obs_asw)                
        ## SWC            
        obs_mean=[0]*len(date_obs_grouped)
        obs_stdev=[0]*len(date_obs_grouped)
        for i in range(0,len(date_obs_grouped)):
            if date_obs_grouped[i]-g0 <= len(zr)-1:
                obs_mean[i]=obs_asw_mean[i]/(zr[date_obs_grouped[i]-g0]*1000) + wp 
                obs_stdev[i]=obs_asw_stdev[i]/(zr[date_obs_grouped[i]-g0]*1000)
            else:
                obs_mean[i]=np.nan
                obs_stdev[i]=np.nan
        # print(obs_mean)
        # print(obs_stdev)
        
        ## 0-5cm
        if opkomst:
            obs0_5_mean, obs0_5_stdev, date_obs_grouped_0_5 = agg_func(obs0_5_date,obs0_5)                
            obs0_5_asw_mean, obs0_5_asw_stdev, date_obs_grouped_0_5 = agg_func(obs0_5_date,obs0_5_asw)                

        ## vol SWC of 0-30
        obs0_30_mean, obs0_30_stdev, date_obs_grouped = agg_func(date_Obs,obs0_30)                

        
    else:
        date_obs_grouped=[]
        obs_asw_mean=[]
        obs_asw_stdev=[]
        obs_mean=[]
        obs_stdev=[]
        obs0_30_mean=[]
        obs0_30_stdev=[]
        
    obs_sim=[]
    for i in range(len(date_obs_grouped)):
        if np.any(np.array(g_list)==date_obs_grouped[i]):
            index=np.where(np.array(g_list)==date_obs_grouped[i])[0][0]
            obs_sim.append(SWC[index])
        else: 
            obs_sim.append(np.nan)

    if opkomst:
        obs_sim_5=[]
        for i in range(len(date_obs_grouped_0_5)):
            if np.any(np.array(g_list)==date_obs_grouped_0_5[i]):
                index=np.where(np.array(g_list)==date_obs_grouped_0_5[i])[0][0]
                obs_sim_5.append(v_5_list[index])
            else: 
                obs_sim_5.append(np.nan)  
        df_obs_grouped_0_5 = pd.DataFrame({'Date': date_obs_grouped_0_5, \
                       'Mean5': obs0_5_mean, 'Stdev5': obs0_5_stdev, \
                       'Sim': obs_sim_5})
  
    df_obs_grouped = pd.DataFrame({'Date': date_obs_grouped, 'Mean': obs_mean, 'Stdev': obs_stdev, \
                           'Mean30': obs0_30_mean, 'Stdev30': obs0_30_stdev, \
                           'Sim': obs_sim})
    
    df_obs_grouped.sort_values(by=['Date'],inplace=True,ignore_index=True)
    
    # print(df_obs_grouped)

    #### Sensor data ----------------------------------------------------------

    if sensor:
        # from Sensordata import sensordata
        ## calibration with root zone soil samples
        # sensor_data, covar = sensordata(case,year,g0,date_obs_grouped,obs_mean,obs_stdev,cal,show=show)
        ## or calibration with 0-30cm soil samples
        # sensor_data, covar, df_calpar = sensordata(date_YYYYMMDD,case,year,g0,date_obs_grouped,obs0_30_mean,obs0_30_stdev,cal,show=show)
        if not isnan(p2) and part=='1': enddate_sensor = str(min(int(ConvertToDate(p2).strftime("%Y%m%d")),int(date_YYYYMMDD)))
        else: enddate_sensor = date_YYYYMMDD
        
        if opkomst:
            sensor_data, covar, df_calpar = sd_func(enddate_sensor,case,year,g0,
                       df_obs_grouped_0_5.Date.tolist(),df_obs_grouped_0_5.Mean5.tolist(),
                       df_obs_grouped_0_5.Stdev5.tolist(),cal,show=show_vec,drop_samp=drop_samp)
        else:
            sensor_data, covar, df_calpar = sd_func(enddate_sensor,case,year,g0,
                       df_obs_grouped.Date.tolist(),df_obs_grouped.Mean30.tolist(),
                       df_obs_grouped.Stdev30.tolist(),cal,show=show_vec,drop_samp=drop_samp)

        # print(sensor_data)
        # sensor_data.dropna(inplace=True)
        # sensor_data.reset_index(drop=True,inplace=True)
        if max(sensor_data['Date'])>g0+end: # if more sensordata, longer than growing period
            # sensor_max=sensor_data.loc[(sensor_data['Date']==g0+end)].index.values[0]-1
            low_closest = sensor_data.iloc[(sensor_data['Date'][sensor_data['Date']<=g0+end]-g0-end).abs().argsort()].dropna()
            sensor_max = int(low_closest.index[0])-1
            # print(sensor_max)
        else: sensor_max=len(sensor_data['Date'])-1
        date_sensor=sensor_data['Date'][:sensor_max+1]
        # data given as m³/m³
        sensordata=sensor_data['Mean'][:sensor_max+1]
        sensordata_stdev=sensor_data['Stdev'][:sensor_max+1]
        covar_adj=covar[:sensor_max+1,:sensor_max+1]
        # sensor -> asw
        sensordata_asw=[]
        sensordata_asw_stdev=[]
        for i in range(0,sensor_max+1):
            sensordata_asw.append(sensordata[i]*zr[i+min(date_sensor)-g0]*1000-wp*zr[i+min(date_sensor)-g0]*1000)
            sensordata_asw_stdev.append(sensordata_stdev[i]*zr[i+min(date_sensor)-g0]*1000)
        if opkomst:
            sensordata_asw5=[]
            sensordata_asw5_stdev=[]
            for i in range(0,sensor_max+1):
                sensordata_asw5.append(sensordata[i]*50-wp*50)
                sensordata_asw5_stdev.append(sensordata_stdev[i]*50)

        N=3 # 3 sensors
        CI_sensor=[stats.t.ppf(1-0.05/2, N-1)*x/sqrt(N) for x in sensordata_stdev]
        # print(CI_sensor)
        low_sensor=np.array(sensordata)-np.array(CI_sensor)
        upper_sensor=np.array(sensordata)+np.array(CI_sensor)
        
        CI_sensor_asw=[stats.t.ppf(1-0.05/2, N-1)*x/sqrt(N) for x in sensordata_asw_stdev]
        low_sensor_asw=np.array(sensordata_asw)-np.array(CI_sensor_asw)
        upper_sensor_asw=np.array(sensordata_asw)+np.array(CI_sensor_asw)
        if opkomst:
            CI_sensor_asw5=[stats.t.ppf(1-0.05/2, N-1)*x/sqrt(N) for x in sensordata_asw5_stdev]
            low_sensor_asw5=np.array(sensordata_asw5)-np.array(CI_sensor_asw5)
            upper_sensor_asw5=np.array(sensordata_asw5)+np.array(CI_sensor_asw5)
        
                
        # Array with the sensordata that match soil samples
        sensor_obs=[]
        sensor_stdev_obs=[]
        for x in range(len(date_obs_grouped)):
            if np.any(date_sensor == date_obs_grouped[x]):
                sensor_obs.append(sensordata.loc[(date_sensor == date_obs_grouped[x])].tolist())
                sensor_stdev_obs.append(sensordata_stdev.loc[(date_sensor == date_obs_grouped[x])].tolist())
            else: 
                sensor_obs.append([np.nan])
                sensor_stdev_obs.append([np.nan])
        # print(sensor_obs,sensor_stdev_obs)
        try:
            sensor_obs=np.stack(sensor_obs).reshape((len(date_obs_grouped)))
            sensor_stdev_obs=np.stack(sensor_stdev_obs).reshape((len(date_obs_grouped)))
        except:
            sensor_obs=sensor_obs
            sensor_stdev_obs=sensor_stdev_obs
        # print(sensor_obs,sensor_stdev_obs)
        
        sensor_data_adj = pd.DataFrame({'Date': date_sensor, 'Mean': sensordata, 'Stdev': sensordata_stdev})
        
        df_obs_grouped['Sensor obs']= sensor_obs
        df_obs_grouped['Sensor stdev obs']= sensor_stdev_obs
        print(df_obs_grouped)
        
        if opkomst:
            # Array with the sensordata that match soil samples
            sensor_obs_5=[]
            sensor_stdev_obs_5=[]
            for x in range(len(date_obs_grouped_0_5)):
                if np.any(date_sensor == date_obs_grouped_0_5[x]):
                    sensor_obs_5.append(sensordata.loc[(date_sensor == date_obs_grouped_0_5[x])].tolist())
                    sensor_stdev_obs_5.append(sensordata_stdev.loc[(date_sensor == date_obs_grouped_0_5[x])].tolist())
                else: 
                    sensor_obs_5.append([np.nan])
                    sensor_stdev_obs_5.append([np.nan])
            sensor_obs_5=np.stack(sensor_obs_5).reshape((len(date_obs_grouped_0_5)))
            sensor_stdev_obs_5=np.stack(sensor_stdev_obs_5).reshape((len(date_obs_grouped_0_5)))
            # print(sensor_obs,sensor_stdev_obs)
            df_obs_grouped_0_5['Sensor obs']= sensor_obs_5
            df_obs_grouped_0_5['Sensor stdev obs']= sensor_stdev_obs_5
        
        
        n_samp = 9 # individual samples in a mixed sample
        if opkomst:
            df_obs_std=df_obs_grouped_0_5['Stdev5'].copy()
            for i in range(len(df_obs_grouped_0_5['Stdev5'])):
                if df_obs_grouped_0_5['Stdev5'][i]==0 or isnan(df_obs_grouped_0_5['Stdev5'][i]):
                    # df_obs_std[i]=0.004     # calculated from Kinrooi & PSKW samples    
                    df_obs_std[i]=0.0114/np.sqrt(n_samp)
            n=df_obs_grouped_0_5.count()['Sim']
        else:
            df_obs_std=df_obs_grouped['Stdev30'].copy()
            for i in range(len(df_obs_grouped['Stdev30'])):
                if df_obs_grouped['Stdev30'][i]==0 or isnan(df_obs_grouped['Stdev30'][i]):
                    # df_obs_std[i]=0.004     # calculated from Kinrooi & PSKW samples    
                    df_obs_std[i]=0.0114/np.sqrt(n_samp)
            n=df_obs_grouped.count()['Sim']
        
        N=len(sensordata)
        
        if cal_par_on: double=True # soil samples used twice in DREAM
        else: double=False 
        if double:
            n*=2
            obs_stdev_ = np.concatenate([df_obs_std,df_obs_std])
        else: obs_stdev_=df_obs_std

        # n=0 # if only sensor data, no soil samples # if run here?
        
        ## Add extra rows/cols (samples) to covar matrix with variance of samples on diagonal (covar=0)
        tn=N+n
        meas_covar=np.zeros((tn,tn))
        for i in range(tn):
            for j in range(tn):
                if i<N and j<N:
                    meas_covar[i,j]=covar_adj[i,j]
                else:
                    meas_covar[i,j]=0
                if i>=N and j>=N and j==i:
                    meas_covar[i,j]=obs_stdev_[i-N]**2
        covar_adj=meas_covar
        
        if show:
            print(covar_adj)
            
        def change_diagonal(matrix,factor):
            np.fill_diagonal(matrix,np.diag(matrix)*factor) # heel kleine error term (*1.0001) zorgt voor heel grote correlatie --> aanpassen
            diag = np.diag(matrix)
            diag = diag + (diag==0)*np.mean(diag[diag!=0]) # var waardes van 0 vervangen door estimate
            np.fill_diagonal(matrix,diag)
            return matrix
        
        zero_cov = False #TODO
        if zero_cov:
            # if zero covariances 
            covar_adj_zeros = np.zeros([len(covar_adj),len(covar_adj)])
            np.fill_diagonal(covar_adj_zeros,np.diag(covar_adj))
            covar_adj = covar_adj_zeros
        else:
            # if non-zero, multiplied diagonal is possible
            covar_adj = change_diagonal(covar_adj,factor=1) # possibility to adjust factor #TODO
        
        if show:
            print(covar_adj)        
        ## check inversion
        if show:
            print("Matrix condition number after correction =",np.linalg.cond(covar_adj))
            print("Losing up to",round(log10(np.linalg.cond(covar_adj)),1),"digits of accuracy\ncompared to matrix dimensions",len(covar_adj))
            print("LogDeterminant \n",np.linalg.slogdet(covar_adj)[1])
            # print("Inverse check \n",np.dot(covar_adj,np.linalg.inv(covar_adj)))
            print("Rank \n",np.linalg.matrix_rank(covar_adj))
            

    else: 
        sensor_data_adj=np.empty(0)
        covar_adj=np.empty(0)
        sensor_obs=np.empty(0)
        sensor_stdev_obs=np.empty(0)
    
    covar_adj[np.isinf(covar_adj)] = 0
    
    df_obs_grouped.dropna(inplace=True,subset = ['Sim'])
    df_obs_grouped.reset_index(inplace=True,drop=True)
    
    # print(df_obs_grouped)

    ### DA Simulation --------------------------------------------------------
    if sensor:
        diff=min(date_sensor)-g0
        if sensor_cal.size:
            ASW_sensor_cal=[(sensor_cal[i-diff]-wp)*zr[i]*1000 for i in range(diff,len(sensor_cal)+diff)]
            
    #### Plots ----------------------------------------------------------------
    
    ## Rooting depth and virtual depth
    
    if show or save:
        plt.figure()
        plt.plot(zr,label='Rooting depth')
        plt.plot(zvirt_list,label='Virtual depth')
        plt.ylabel('Depth (m)')
        plt.xlim(0, end-1)
        plt.legend(loc="upper left")
        ax1=plt.gca()
        ax2 = ax1.twinx()
        ax2.set_ylabel('Precipitation & irrigation (mm)')
        ax2.set_ylim(0, 80)
        ax2.bar([x - int(g0) -0.5 for x in g_list],I_list, label='Irrigation',color='#fee090',alpha=0.6)
        ax2.bar([x - int(g0) -0.5 for x in date_R_I],rain, label='Rainfall',color='#91bfdb',alpha=0.6)
        plt.legend(loc="upper right")
        if show and not save: plt.show()
        elif save: 
            plt.savefig(folder+'\\'+case+'_'+year+'_Rooting_virtual_depth.png',dpi=300)
            plt.close()        
    
    if show or save:
        ### Soil water balance (mm)
        plt.figure(case+'_'+year+'_Soil water balance (mm)')
        plt.plot([x - int(g0) for x in g_list],sw_list, label='Soil water storage',color='#4575b4')
        plt.plot([x - int(g0) for x in g_list],swfc_list, label='Field capacity',color='#fc8d59')
        plt.plot([x - int(g0) for x in g_list],sw_vkrit_list, label='Water stress threshold',color='#d73027')
        if CI.size and CI.shape[1] == len(g_list):        
            plt.fill_between([x - int(g0) for x in g_list], CI[2], CI[3], alpha = 0.2,color='#4575b4',label='95% CI')
        plt.bar([x - int(g0) -0.5 for x in g_list],I_list, label='Irrigation',color='#fee090',alpha=0.6)
        plt.bar([x - int(g0) -0.5 for x in date_R_I],rain, label='Rainfall',color='#91bfdb',alpha=0.6)
        # if sensor: 
        #     plt.plot([x - int(g0) for x in date_sensor],sensordata_asw, label='Sensordata TEROS 10',color='k')
        #     plt.fill_between([x - int(g0) for x in date_sensor], low_sensor_asw, upper_sensor_asw, alpha = 0.1,color='k',label='95% CI')
        # if sensor_cal.size:
        #     plt.plot([x - int(g0) for x in date_sensor],ASW_sensor_cal, label='Calibrated sensor data',color='k',linestyle='dashed')
        # plt.errorbar([x - int(g0) for x in date_obs_grouped],obs_asw_mean, obs_asw_stdev, linestyle='None', marker='o', capsize=3, label='Soil water observation',color='k')
        # plt.scatter([x - int(g0) for x in date_Obs],obs0_30_asw, facecolors='none', edgecolors='k',s=6)
        # if opkomst: plt.scatter([x - int(g0) for x in obs0_5_date],obs0_5_asw, marker='d',s=10, facecolors='none', edgecolors='gray', label='Soil sample 5 cm')
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        plt.ylabel('Available soil water storage (mm)')
        plt.xlim([0,end])
        plt.ylim([0,max(swfc_list)+20])
        if date_axis:
            plt.xticks(ticks=range(0,len(g_list),7),
                   minor=True)
            plt.xticks(ticks=range(0,len(g_list),14),
                   labels=[ConvertToDate(g_list[i]).strftime("%d/%m") for i in range(0,len(g_list),14)],
                   ha='center')
            plt.xlabel('Date\n(End of the day)')
        else:
            plt.xlabel('Growing days\n(End of the day)')        
        handles, labels = plt.gca().get_legend_handles_labels()
        if  CI.size and CI.shape[1] == len(g_list):
            comb_handles = [(handles[0],handles[3])]+[handles[i] for i in [1,2,4,5]]
            comb_labels = ['Soil water storage with 95% CI']+[labels[i] for i in [1,2,4,5]]
            plt.legend(comb_handles,comb_labels,loc='upper center',bbox_to_anchor=(0.5,1.1),
                       ncol=2,fancybox=True,shadow=True)
        else: plt.legend(loc='upper center',bbox_to_anchor=(0.5,1.1),ncol=2,fancybox=True,shadow=True)
        plt.tight_layout()  
        plt.subplots_adjust(top=0.9)
        if show and not save: plt.show()
        elif save: 
            plt.savefig(folder+'\\'+case+'_'+year+'_Soil_water_balance_(mm).png',dpi=300)
            plt.close()

        
        ### Soil water content (m³/m³)
        fig, ax1 = plt.subplots(num=case+'_'+year+'_SWC')
        ax1.set_ylabel('Soil water content (m³/m³)')
        ax1.set_ylim(0, 0.5)
        ax1.set_xlim(0,end)
        if date_axis:
            ax1.set_xticks(ticks=range(0,len(g_list),7),
                   minor=True)
            ax1.set_xticks(ticks=range(0,len(g_list),14),
                   labels=[ConvertToDate(g_list[i]).strftime("%d/%m") for i in range(0,len(g_list),14)],
                   ha='center')
            ax1.set_xlabel('Date\n(End of the day)')
        else:
            ax1.set_xlabel('Growing days\n(End of the day)')
        ax1.plot([x - int(g0) for x in g_list],SWC, label='Soil water content',color='#4575b4')
        ax1.plot([x - int(g0) for x in g_list],np.ones(len(g_list))*fc, label='Field capacity',color='#fc8d59')
        ax1.plot([x - int(g0) for x in g_list],vkrit_list, label='Water stress threshold',color='#d73027')
        if sensor:
            ax1.plot([x - int(g0) for x in date_sensor],sensordata, label='Sensordata TEROS 10',color='k')
            ax1.fill_between([x - int(g0) for x in date_sensor], low_sensor, upper_sensor, alpha = 0.1,color='k',label='95% CI')
        if CI.size and CI.shape[1] == len(g_list):        
            ax1.fill_between([x - int(g0) for x in g_list], CI[0], CI[1], alpha = 0.3,color='#4575b4',label='95% CI')
        elif CI.size and CI.shape[1] < len(g_list):        
            ax1.fill_between([x - int(g0) for x in g_list[:CI.shape[1]]], CI[0], CI[1], alpha = 0.3,color='#4575b4',label='95% CI')
        elif CI.size and CI.shape[1] > len(g_list):        
            ax1.fill_between([x - int(g0) for x in g_list], [CI[0][x - int(g0)] for x in g_list],
                                 [CI[1][x - int(g0)] for x in g_list], alpha = 0.3,color='#4575b4',label='95% CI')
        if sensor_cal.size:
            ax1.plot([x - int(g0) for x in date_sensor],sensor_cal, label='Calibrated sensor data',color='k',linestyle='dashed')
        if df_forecast.size: ax1.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        # Cropped
        if CI.size and CI.shape[1] < len(g_list):
            plt.fill_between([g_list[CI.shape[1]-1]+0.5 - int(g0), g_list[-1]+1 - int(g0)], 0, 0.5, color='white',zorder=20)
            plt.fill_between([x - int(g0)-9.5 for x in g_list[CI.shape[1]-1:]], 0, 0.5, color='white', 
                             edgecolor='none',alpha=0.4,zorder=21)
            ax1.axvline(g_list[CI.shape[1]-1]-9.5 - int(g0), color='k', linestyle='dotted',zorder=22)
        # Cropped - end
        elif CI.size and CI.shape[1] > len(g_list):
            if CI.shape[1] == len(g_list)+1:
                n_pred = (CI.shape[1]-2)%10
                plt.fill_between([x - int(g0)+0.5 for x in g_list[CI.shape[1]-n_pred-3:]]+[g_list[-1]+1], 0, 0.5, color='white', 
                                 edgecolor='none',alpha=0.4,zorder=21)
                ax1.axvline(g_list[CI.shape[1]-n_pred-3]+0.5 - int(g0), color='k', linestyle='dotted',zorder=22)
            elif CI.shape[1] > len(g_list)+3: # harvest is earlier (e.g. PCG prei 2022)
                n_pred = (CI.shape[1]-20)%10
                plt.fill_between([x - int(g0)+0.5 for x in g_list[CI.shape[1]-n_pred-3:]]+[g_list[-1]+1], 0, 0.5, color='white', 
                                  edgecolor='none',alpha=0.4,zorder=21)
                ax1.axvline(g_list[CI.shape[1]-n_pred-3]+0.5 - int(g0), color='k', linestyle='dotted',zorder=22)           

        # 2nd axis
        ax2 = ax1.twinx()
        ax2.set_ylabel('Precipitation & irrigation (mm)')
        ax2.set_ylim(0, 80)
        ax2.bar([x - int(g0) -0.5 for x in g_list],I_list, label='Irrigation',color='#fee090',alpha=0.6)
        ax2.bar([x - int(g0) -0.5 for x in date_R_I],rain, label='Rainfall',color='#91bfdb',alpha=0.6)
        ax1.errorbar([x - int(g0) for x in date_obs_grouped],obs0_30_mean,obs0_30_stdev, linestyle='None', marker='o', capsize=3, label='Soil water observation',color='k')
        if opkomst: ax1.scatter([x - int(g0) for x in obs0_5_date],obs0_5,  marker='d',s=10, facecolors='none', edgecolors='gray', label='Soil sample 5 cm')
        
        handles1, labels1 = ax1.get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()
        if  CI.size and sensor:
            comb_handles = [(handles1[0],handles1[5])]+[handles1[i] for i in [1,2]]+ \
                            [plt.Line2D([],[], alpha=0),(handles1[3],handles1[4]),handles1[-1]]+[handles2[i] for i in [0,1]]
            comb_labels = ['Soil water content \nwith 95% CI']+[labels1[i] for i in [1,2]]+ \
                            ['','Sensor data \nwith 95% CI',labels1[-1]]+[labels2[i] for i in [0,1]]
            fig.legend(comb_handles,comb_labels,loc='upper center',bbox_to_anchor=(0.5,1),
                       ncol=2,fancybox=True,shadow=True)
        elif sensor and not CI.size:
            comb_handles = [handles1[i] for i in [0,1,2]]+ \
                            [plt.Line2D([],[], alpha=0),(handles1[3],handles1[4]),handles1[-1]]+[handles2[i] for i in [0,1]]
            comb_labels = [labels1[i] for i in [0,1,2]]+['','Sensor data with 95% CI']+ \
                            [labels1[-1]]+[labels2[i] for i in [0,1]]
            fig.legend(comb_handles,comb_labels,loc='upper center',bbox_to_anchor=(0.5,1),
                       ncol=2,fancybox=True,shadow=True)            
        else: 
            fig.legend(loc='upper center',bbox_to_anchor=(0.5,1),ncol=2,fancybox=True,shadow=True)

        fig.tight_layout()                
        plt.subplots_adjust(top=0.9)
        if show and not save: plt.show()
        elif save: 
            if CI.size and (CI.shape[1] < len(g_list) or CI.shape[1] > len(g_list)):
                plt.savefig(folder+'\\'+case+'_'+year+'_SWC_CROPPED.png',dpi=300)
            else: plt.savefig(folder+'\\'+case+'_'+year+'_SWC.png',dpi=300)
            plt.close()    
            
        
        #### PAPER FIGURES ----------------------------------------------------
        fs=(6,4.5)
        fonts=12
        ### Soil water content (m³/m³) - PAPER FIG
        fig, ax1 = plt.subplots(num=case+'_'+year+'_SWC_PAPER',figsize=fs)
        ax1.set_ylabel('Soil water content (m³/m³)',labelpad=10,fontsize=fonts)
        ax1.set_ylim(0, 0.5)
        ax1.set_xlim(0,end)
        if date_axis:
            ax1.set_xticks(ticks=range(0,len(g_list),7),
                   minor=True)
            ax1.set_xticks(ticks=range(0,len(g_list),14),
                   labels=[ConvertToDate(g_list[i]).strftime("%d/%m") for i in range(0,len(g_list),14)],
                   ha='center')
            ax1.set_xlabel('Date',labelpad=10,fontsize=fonts)
        else:
            ax1.set_xlabel('Growing days\n(End of the day)',labelpad=10,fontsize=fonts)
        ax1.plot([x - int(g0) for x in g_list],SWC, label='Soil water content',color='#4575b4')
        ax1.plot([x - int(g0) for x in g_list],np.ones(len(g_list))*fc, label='Field capacity',color='#fc8d59')
        ax1.plot([x - int(g0) for x in g_list],vkrit_list, label='Water stress threshold',color='#d73027')
        if sensor:
            ax1.plot([x - int(g0) for x in date_sensor],sensordata, label='Sensordata TEROS 10',color='k')
        if CI.size and CI.shape[1] == len(g_list):        
            ax1.fill_between([x - int(g0) for x in g_list], CI[0], CI[1], alpha = 0.3,color='#4575b4',label='95% CI')
        elif CI.size and CI.shape[1] < len(g_list):        
            ax1.fill_between([x - int(g0) for x in g_list[:CI.shape[1]]], CI[0], CI[1], alpha = 0.3,color='#4575b4',label='95% CI')
        elif CI.size and CI.shape[1] > len(g_list):        
            ax1.fill_between([x - int(g0) for x in g_list], [CI[0][x - int(g0)] for x in g_list],
                                 [CI[1][x - int(g0)] for x in g_list], alpha = 0.3,color='#4575b4',label='95% CI')
        if sensor_cal.size:
            ax1.plot([x - int(g0) for x in date_sensor],sensor_cal, label='Calibrated sensor data',color='k',linestyle='dashed')
        if df_forecast.size: ax1.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        # Cropped
        if CI.size and CI.shape[1] < len(g_list):
            plt.fill_between([g_list[CI.shape[1]-1]+0.5 - int(g0), g_list[-1]+1 - int(g0)], 0, 0.5, color='white',zorder=20)
            plt.fill_between([x - int(g0)-9.5 for x in g_list[CI.shape[1]-1:]], 0, 0.5, color='white', 
                             edgecolor='none',alpha=0.4,zorder=21)
            ax1.axvline(g_list[CI.shape[1]-1]-9.5 - int(g0), color='k', linestyle='dotted',zorder=22)
        # Cropped - end
        elif CI.size and CI.shape[1] > len(g_list):
            if CI.shape[1] == len(g_list)+1:
                n_pred = (CI.shape[1]-2)%10
                plt.fill_between([x - int(g0)+0.5 for x in g_list[CI.shape[1]-n_pred-3:]]+[g_list[-1]+1], 0, 0.5, color='white', 
                                 edgecolor='none',alpha=0.4,zorder=21)
                ax1.axvline(g_list[CI.shape[1]-n_pred-3]+0.5 - int(g0), color='k', linestyle='dotted',zorder=22)
            else: # harvest is earlier (e.g. PCG prei 2022)
                n_pred = (CI.shape[1]-40)%10
                plt.fill_between([x - int(g0)+0.5 for x in g_list[CI.shape[1]-n_pred-3:]]+[g_list[-1]+1], 0, 0.5, color='white', 
                                  edgecolor='none',alpha=0.4,zorder=21)
                ax1.axvline(g_list[CI.shape[1]-n_pred-3]+0.5 - int(g0), color='k', linestyle='dotted',zorder=22)           

        # 2nd axis
        ax2 = ax1.twinx()
        ax2.set_ylabel('Precipitation & irrigation (mm)',labelpad=10,fontsize=fonts)
        ax2.set_ylim(0, 80)
        ax2.bar([x - int(g0) -0.5 for x in g_list],I_list, label='Irrigation',color='#fee090',alpha=0.6)
        ax2.bar([x - int(g0) -0.5 for x in date_R_I],rain, label='Rainfall',color='#91bfdb',alpha=0.6)
        ax1.errorbar([x - int(g0) for x in date_obs_grouped],obs0_30_mean,0, linestyle='None', marker='o', capsize=3, label='Soil water observation',color='k')
        if opkomst: ax1.scatter([x - int(g0) for x in obs0_5_date],obs0_5,  marker='d',s=10, facecolors='none', edgecolors='gray', label='Soil sample 5 cm')
        fig.tight_layout()                
        plt.subplots_adjust(top=0.9)
        if show and not save: plt.show()
        elif save: 
            if CI.size and (CI.shape[1] < len(g_list) or CI.shape[1] > len(g_list)):
                plt.savefig(folder+'\\'+case+'_'+year+'_SWC_CROPPED_PAPER.png',dpi=300)
            else: plt.savefig(folder+'\\'+case+'_'+year+'_SWC_PAPER.png',dpi=300)
            plt.close()    
                
        ### Soil water balance (mm) - PAPER FIG
        plt.figure(case+'_'+year+'_Soil water balance (mm) PAPER',figsize=fs)
        plt.plot([x - int(g0) for x in g_list],sw_list, label='Soil water storage',color='#4575b4')
        plt.plot([x - int(g0) for x in g_list],swfc_list, label='Field capacity',color='#fc8d59')
        plt.plot([x - int(g0) for x in g_list],sw_vkrit_list, label='Water stress threshold',color='#d73027')
        if CI.size and CI.shape[1] == len(g_list):        
            plt.fill_between([x - int(g0) for x in g_list], CI[2], CI[3], alpha = 0.2,color='#4575b4',label='95% CI')
        plt.bar([x - int(g0) -0.5 for x in g_list],I_list, label='Irrigation',color='#fee090',alpha=0.6)
        plt.bar([x - int(g0) -0.5 for x in date_R_I],rain, label='Rainfall',color='#91bfdb',alpha=0.6)
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        plt.ylabel('Available soil water storage (mm)''\n''Precipitation & irrigation (mm)',labelpad=10,fontsize=fonts)
        plt.xlim([0,end])
        # plt.ylim([0,max(swfc_list)+20])
        plt.ylim(bottom=0)
        if date_axis:
            plt.xticks(ticks=range(0,len(g_list),7),
                   minor=True)
            plt.xticks(ticks=range(0,len(g_list),14),
                   labels=[ConvertToDate(g_list[i]).strftime("%d/%m") for i in range(0,len(g_list),14)],
                   ha='center')
            plt.xlabel('Date',labelpad=10,fontsize=fonts)
        else:
            plt.xlabel('Growing days\n(End of the day)',labelpad=10, fontsize=fonts)        
        plt.tight_layout()  
        plt.subplots_adjust(top=0.9)
        if show and not save: plt.show()
        elif save: 
            plt.savefig(folder+'\\'+case+'_'+year+'_Soil_water_balance_(mm)_PAPER.png',dpi=300)
            plt.close()
            
        ### Virtual rooting depth
        plt.figure(figsize=fs)
        plt.plot([x - int(g0) for x in g_list],zr[:end],label='Rooting depth')
        plt.plot([x - int(g0) for x in g_list],zvirt_list[:end],label='Virtual depth')
        plt.ylabel('(Virtual) rooting depth (m)',labelpad=10,fontsize=fonts)
        plt.xlim(0, end-1)
        plt.xticks(ticks=range(0,len(g_list),7),
               minor=True)
        plt.xticks(ticks=range(0,len(g_list),14),
               labels=[ConvertToDate(g_list[i]).strftime("%d/%m") for i in range(0,len(g_list),14)],
               ha='center')
        plt.xlabel('Date',labelpad=10,fontsize=fonts)
        plt.legend(loc="upper left")
        ax1=plt.gca()
        ax2 = ax1.twinx()
        ax2.set_ylabel('Precipitation & irrigation (mm)',labelpad=10,fontsize=fonts)
        ax2.set_ylim(0, 80)
        ax2.bar([x - int(g0) -0.5 for x in g_list],I_list,color='#fee090',alpha=0.6)
        ax2.bar([x - int(g0) -0.5 for x in date_R_I],rain,color='#91bfdb',alpha=0.6)
        plt.tight_layout()
        plt.subplots_adjust(top=0.9)
        if show and not save: plt.show()
        elif save: 
            plt.savefig(folder+'\\'+case+'_'+year+'_Rooting_virtual_depth_PAPER.png',dpi=300)
            plt.close()       
            
            
        ### All-plots-plot - PAPER
        fonts_mult=16
        plt.figure(case+'_'+year+'_all-plots',figsize=(16,7))
        plt.subplot(231)
        plt.plot([x - int(g0) for x in g_list],etm_list, label=r'$\mathrm{ET_m}$',color='#4575b4')
        plt.plot([x - int(g0) for x in g_list],eta_list, label=r'$\mathrm{ET_a}$',color='#4575b4',linestyle='dashed')
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        plt.ylabel('ET (mm)',labelpad=8,fontsize=fonts_mult)
        plt.legend(loc='upper right', fontsize=fonts_mult-2)
        plt.xlabel('')
        ax = plt.gca()
        ax.set_ylim(bottom=0)
        ax.set(xlim=(0, end), xticks=range(0,end,10))
        ax.set_xticklabels(ax.get_xticks(), fontsize=fonts_mult-4)
        ax.set_yticklabels([f"{int(y)}" for y in ax.get_yticks()], fontsize=fonts_mult-4)

        
        plt.subplot(233)
        plt.plot([x - int(g0) for x in g_list],Dr_list,label='Depletion root zone',color='#4575b4')
        plt.plot([x - int(g0) for x in g_list],TAW_list, label='TAW',linestyle='dashed',color='#d73027') #TEW
        plt.plot([x - int(g0) for x in g_list],RAW_list, label='RAW',linestyle='dotted',color='#d73027') #REW
        plt.ylabel(r'$D_\mathrm{r}}$ (mm)',labelpad=8,fontsize=fonts_mult)
        plt.xlabel('')
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        plt.legend(loc='upper right', fontsize=fonts_mult-2)
        ax = plt.gca()
        ax.set_ylim(bottom=0)
        ax.set(xlim=(0, end),xticks=range(0,end,10))
        ax.set_xticklabels(ax.get_xticks(), fontsize=fonts_mult-4)
        ax.set_yticklabels([f"{int(y)}" for y in ax.get_yticks()], fontsize=fonts_mult-4)

        
        plt.subplot(236)
        plt.plot([x - int(g0) for x in g_list],De_list,label='Depletion 10 cm',color='#4575b4')
        plt.plot([x - int(g0) for x in g_list],TEW, label='TEW',linestyle='dashed',color='#d73027') #TEW
        plt.plot([x - int(g0) for x in g_list],REW, label='REW',linestyle='dotted',color='#d73027') #REW
        plt.ylabel(r'$D_\mathrm{e}}$ (mm)',labelpad=8,fontsize=fonts_mult)
        plt.xlabel('Growing days',labelpad=8,fontsize=fonts_mult)
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        plt.legend(loc='upper right', fontsize=fonts_mult-2)
        ax = plt.gca()
        ax.set_ylim(bottom=0)
        ax.set(xlim=(0, end),xticks=range(0,end,10))
        ax.set_xticklabels(ax.get_xticks(), fontsize=fonts_mult-4)
        ax.set_yticklabels([f"{int(y)}" for y in ax.get_yticks()], fontsize=fonts_mult-4)

        
        plt.subplot(234)
        plt.plot([x - int(g0) for x in g_list],ETcb_list, label=r'$\mathrm{ET_{cb}}$', color='#fc8d59')
        plt.plot([x - int(g0) for x in g_list],ETcb_a_list, label=r'$\mathrm{ET_{cb,adj}}$', linestyle='dashed', color='#fc8d59')
        plt.plot([x - int(g0) for x in g_list],ETe_list, label=r'$\mathrm{ET_{e}}$', color='#91bfdb')
        plt.ylabel('ET (mm)',labelpad=8,fontsize=fonts_mult)
        plt.xlabel('Growing days',labelpad=8,fontsize=fonts_mult)
        plt.legend(loc='upper right', fontsize=fonts_mult-2)
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        ax = plt.gca()
        ax.set_ylim(bottom=0)
        ax.set(xlim=(0, end),xticks=range(0,end,10))
        ax.set_xticklabels(ax.get_xticks(), fontsize=fonts_mult-4)
        ax.set_yticklabels([f"{int(y)}" for y in ax.get_yticks()], fontsize=fonts_mult-4)

        
        plt.subplot(232)
        plt.plot(kcb, label=r'$K_\mathrm{cb}$', color='#fc8d59')
        plt.plot(kcb_a_list, label=r'$K_\mathrm{cb,adj}$', linestyle='dashed', color='#fc8d59')
        plt.plot([x - int(g0) for x in g_list],ke_list,label=r'$K_\mathrm{e}$', color='#91bfdb')
        plt.plot([x - int(g0) for x in g_list],kc_list,label=r'$K_\mathrm{c}$', color='#4575b4')
        plt.plot([x - int(g0) for x in g_list],kc_a_list,label=r'$K_\mathrm{c,adj}$', linestyle='dashed', color='#4575b4')
        plt.ylabel('Crop coefficient',labelpad=8,fontsize=fonts_mult)
        plt.xlabel('')
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        plt.legend(loc='upper right', fontsize=fonts_mult-2)
        ax = plt.gca()
        ax.set_ylim(bottom=0)
        ax.set(xlim=(0, end),xticks=range(0,end,10))
        ax.set_xticklabels(ax.get_xticks(), fontsize=fonts_mult-4)
        ax.set_yticklabels([f"{round(y,1)}" for y in ax.get_yticks()], fontsize=fonts_mult-4)

        
        plt.subplot(235)
        plt.plot([x - int(g0) for x in g_list],etcb_cum_list, label=r'Cumulative $\mathrm{ET_{cb}}$', color='#fc8d59')
        plt.plot([x - int(g0) for x in g_list],etcb_a_cum_list, label=r'Cumulative $\mathrm{ET_{cb,adj}}$', linestyle='dashed', color='#fc8d59')
        plt.plot([x - int(g0) for x in g_list],ete_cum_list, label=r'Cumulative $\mathrm{ET_{e}}$', color='#91bfdb')
        plt.plot([x - int(g0) for x in g_list],etm_cum_list, label=r'Cumulative $\mathrm{ET_{m}}$', color='#4575b4')
        plt.plot([x - int(g0) for x in g_list],eta_cum_list, label=r'Cumulative $\mathrm{ET_{a}}$', linestyle='dashed', color='#4575b4')
        plt.ylabel('Cumulative ET (mm)',labelpad=8,fontsize=fonts_mult)
        plt.xlabel('Growing days',labelpad=8,fontsize=fonts_mult)
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        plt.legend(loc='upper right', fontsize=fonts_mult-2)
        ax = plt.gca()
        ax.set_ylim(bottom=0)
        ax.set(xlim=(0, end),xticks=range(0,end,10))
        ax.set_xticklabels(ax.get_xticks(), fontsize=fonts_mult-4)
        ax.set_yticklabels([f"{int(y)}" for y in ax.get_yticks()], fontsize=fonts_mult-4)

        plt.subplots_adjust(left=0.04, bottom=0.1, right=0.98, top=0.93, wspace=None, hspace=None)
        plt.tight_layout()
        
        if show and not save: plt.show()
        elif save: 
            plt.savefig(folder+'\\'+case+'_'+year+'_all-plots_PAPER.png',dpi=300)
            plt.close()  
            
       
    if show or save:
        ### All plots
        plt.figure(case+'_'+year+'_all-plots',figsize=(16,8))
        plt.subplot(231)
        plt.plot([x - int(g0) for x in g_list],etm_list, label='ETm',color='#4575b4')
        plt.plot([x - int(g0) for x in g_list],eta_list, label='ETa',color='#4575b4',linestyle='dashed')
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        plt.ylabel('ET (mm)')
        plt.legend(loc='upper left')
        plt.xlabel('Growing days')
        ax = plt.gca()
        ax.set(xlim=(0, end), xticks=range(0,end,10))
        
        plt.subplot(233)
        plt.plot([x - int(g0) for x in g_list],Dr_list,label='Depletion root zone',color='#4575b4')
        plt.plot([x - int(g0) for x in g_list],TAW_list, label='TAW',linestyle='dashed',color='#d73027') #TEW
        plt.plot([x - int(g0) for x in g_list],RAW_list, label='RAW',linestyle='dotted',color='#d73027') #REW
        plt.ylabel('Dr (mm)')
        plt.xlabel('Growing days')
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        plt.legend(loc='upper left')
        ax = plt.gca()
        ax.set(xlim=(0, end),xticks=range(0,end,10))
        
        plt.subplot(236)
        plt.plot([x - int(g0) for x in g_list],De_list,label='Depletion 10 cm',color='#4575b4')
        plt.plot([x - int(g0) for x in g_list],TEW, label='TEW',linestyle='dashed',color='#d73027') #TEW
        plt.plot([x - int(g0) for x in g_list],REW, label='REW',linestyle='dotted',color='#d73027') #REW
        plt.ylabel('De (mm)')
        plt.xlabel('Growing days')
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        plt.legend(loc='upper left')
        ax = plt.gca()
        ax.set(xlim=(0, end),xticks=range(0,end,10))
        
        plt.subplot(234)
        plt.plot([x - int(g0) for x in g_list],ETcb_list, label='ETcb', color='#fc8d59')
        plt.plot([x - int(g0) for x in g_list],ETcb_a_list, label='ETcb_adj', linestyle='dashed', color='#fc8d59')
        plt.plot([x - int(g0) for x in g_list],ETe_list, label='ETe', color='#91bfdb')
        plt.ylabel('ET (mm)')
        plt.xlabel('Growing days')
        plt.legend(loc='upper left')
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        ax = plt.gca()
        ax.set(xlim=(0, end),xticks=range(0,end,10))
        
        plt.subplot(232)
        plt.plot(kcb, label='Kcb', color='#fc8d59')
        plt.plot(kcb_a_list, label='Kcb_adj', linestyle='dashed', color='#fc8d59')
        plt.plot([x - int(g0) for x in g_list],ke_list,label='Ke', color='#91bfdb')
        plt.plot([x - int(g0) for x in g_list],kc_list,label='Kc', color='#4575b4')
        plt.plot([x - int(g0) for x in g_list],kc_a_list,label='Kc_adj', linestyle='dashed', color='#4575b4')
        plt.ylabel('Crop coefficient')
        plt.xlabel('Growing days')
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        plt.legend(loc='upper left')
        ax = plt.gca()
        ax.set(xlim=(0, end),xticks=range(0,end,10))
        
        plt.subplot(235)
        plt.plot([x - int(g0) for x in g_list],etcb_cum_list, label='Cumulative ETcb', color='#fc8d59')
        plt.plot([x - int(g0) for x in g_list],etcb_a_cum_list, label='Cumulative ETcb_adj', linestyle='dashed', color='#fc8d59')
        plt.plot([x - int(g0) for x in g_list],ete_cum_list, label='Cumulative ETe', color='#91bfdb')
        plt.plot([x - int(g0) for x in g_list],etm_cum_list, label='Cumulative ETm', color='#4575b4')
        plt.plot([x - int(g0) for x in g_list],eta_cum_list, label='Cumulative ETa', linestyle='dashed', color='#4575b4')
        plt.ylabel('Cumulative ET (mm)')
        plt.xlabel('Growing days')
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        plt.legend(loc='upper left')
        ax = plt.gca()
        ax.set(xlim=(0, end),xticks=range(0,end,10))
        
        plt.subplots_adjust(left=0.04, bottom=0.1, right=0.98, top=0.93, wspace=None, hspace=None)
        plt.suptitle(case+'_'+year, ha='center',weight='bold')
        
        if show and not save: plt.show()
        elif save: 
            plt.savefig(folder+'\\'+case+'_'+year+'_all-plots.png',dpi=300)
            plt.close()  

    # 5 - 10 - root zone
    if show and opkomst:
        
        ### Soil water balance (mm)
        plt.figure(case+'_'+year+'_Soil water balance (mm) Evaporative layer (5 & 10 cm)')
        plt.plot([x - int(g0) for x in g_list],sw_5_list, label='Soil water storage 5 cm',color='#4575b4')
        plt.plot([x - int(g0) for x in g_list],swfc_5_list, label='Field capacity 5 cm',color='#fc8d59')
        plt.plot([x - int(g0) for x in g_list],sw_5_vkrit_list, label='Emergence threshold 5 cm',color='#d73027')
        plt.plot([x - int(g0) for x in g_list],sw_Ze_list, label='Soil water storage 10 cm',color='#91bfdb')
        plt.plot([x - int(g0) for x in g_list],swfc_Ze_list, label='Field capacity 10 cm', linestyle='dotted',color='#fc8d59',alpha=0.6)
        plt.plot([x - int(g0) for x in g_list],sw_Ze_vkrit_list, label='Emergence threshold 10 cm', linestyle='dotted',color='#d73027',alpha=0.6)
        if CI.size:        
            plt.fill_between([x - int(g0) for x in g_list], CI[2], CI[3], alpha = 0.2,color='#4575b4',label='95% CI')
        plt.bar([x - int(g0) -0.5 for x in g_list],I_list, label='Irrigation',color='#fee090',alpha=0.6)
        plt.bar([x - int(g0) -0.5 for x in date_R_I],rain, label='Rainfall',color='#91bfdb',alpha=0.6)
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        plt.ylabel('Available soil water storage (mm)')
        plt.xlim([0,end])
        plt.ylim([0,max(swfc_Ze_list)+20])
        if date_axis:
            plt.xticks(ticks=range(0,len(g_list),7),
                   minor=True)
            plt.xticks(ticks=range(0,len(g_list),7),
                   labels=[ConvertToDate(g_list[i]).strftime("%d/%m") for i in range(0,len(g_list),7)],
                   ha='center')
            plt.xlabel('Date\n(End of the day)')
        else:
            plt.xlabel('Growing days\n(End of the day)')        
        handles, labels = plt.gca().get_legend_handles_labels()
        if  CI.size:
            comb_handles = [(handles[0],handles[6])]+[handles[i] for i in [1,2,7,3,4,5,8]]
            comb_labels = ['Soil water storage 5cm with 95% CI']+[labels[i] for i in [1,2,7,3,4,5,8]]
            plt.legend(comb_handles,comb_labels,loc='upper center',bbox_to_anchor=(0.5,1.1),
                        ncol=2,fancybox=True,shadow=True)
        else: plt.legend(loc='upper center',bbox_to_anchor=(0.5,1.1),ncol=2,fancybox=True,shadow=True)
        # plt.legend(loc='upper center',bbox_to_anchor=(0.5,1.1),ncol=2,fancybox=True,shadow=True)
        plt.tight_layout()  
        plt.subplots_adjust(top=0.9)
        if show and not save: plt.show()
        elif save: 
            plt.savefig(folder+'\\'+case+'_'+year+'_Soil_water_balance_Ze_5_(mm).png',dpi=300)
            plt.close()

        ### Soil water balance 5cm (mm)
        plt.figure(case+'_'+year+'_Soil water balance 5cm (mm)')
        plt.plot([x - int(g0) for x in g_list],sw_5_list, label='Soil water storage 5 cm',color='#4575b4')
        plt.plot([x - int(g0) for x in g_list],swfc_5_list, label='Field capacity 5 cm',color='#fc8d59')
        plt.plot([x - int(g0) for x in g_list],sw_5_vkrit_list, label='Emergence threshold 5 cm',color='#d73027')
        if CI.size:        
            plt.fill_between([x - int(g0) for x in g_list], CI[2], CI[3], alpha = 0.2,color='#4575b4',label='95% CI')
        plt.bar([x - int(g0) -0.5 for x in g_list],I_list, label='Irrigation',color='#fee090',alpha=0.6)
        plt.bar([x - int(g0) -0.5 for x in date_R_I],rain, label='Rainfall',color='#91bfdb',alpha=0.6)
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        plt.ylabel('Available soil water storage (mm)')
        plt.xlim([0,end])
        plt.ylim([0,max(swfc_5_list)+20])
        if date_axis:
            plt.xticks(ticks=range(0,len(g_list),7),
                   minor=True)
            plt.xticks(ticks=range(0,len(g_list),7),
                   labels=[ConvertToDate(g_list[i]).strftime("%d/%m") for i in range(0,len(g_list),7)],
                   ha='center')
            plt.xlabel('Date\n(End of the day)')
        else:
            plt.xlabel('Growing days\n(End of the day)')        
        handles, labels = plt.gca().get_legend_handles_labels()
        if  CI.size:
            comb_handles = [(handles[0],handles[3])]+[handles[i] for i in [1,2,4,5]]
            comb_labels = ['Soil water storage 5cm with 95% CI']+[labels[i] for i in [1,2,4,5]]
            plt.legend(comb_handles,comb_labels,loc='upper center',bbox_to_anchor=(0.5,1.1),
                        ncol=2,fancybox=True,shadow=True)
        else: plt.legend(loc='upper center',bbox_to_anchor=(0.5,1.1),ncol=2,fancybox=True,shadow=True)
        # plt.legend(loc='upper center',bbox_to_anchor=(0.5,1.1),ncol=2,fancybox=True,shadow=True)
        plt.tight_layout()  
        plt.subplots_adjust(top=0.9)
        if show and not save: plt.show()
        elif save: 
            plt.savefig(folder+'\\'+case+'_'+year+'_Soil_water_balance_5_(mm).png',dpi=300)
            plt.close()        
        
        ### Soil water content (m³/m³)
        fig, ax1 = plt.subplots(num=case+'_'+year+'_SWC 5 cm')
        ax1.set_ylabel('Soil water content (m³/m³)')
        ax1.set_ylim(0, 0.5)
        ax1.set_xlim(0,end)
        if date_axis:
            ax1.set_xticks(ticks=range(0,len(g_list),7),
                   minor=True)
            ax1.set_xticks(ticks=range(0,len(g_list),7),
                   labels=[ConvertToDate(g_list[i]).strftime("%d/%m") for i in range(0,len(g_list),7)],
                   ha='center')
            ax1.set_xlabel('Date\n(End of the day)')
        else:
            ax1.set_xlabel('Growing days\n(End of the day)')
        # ax1.plot([x - int(g0) for x in g_list],v_Ze_list, label='SWC 10 cm',color='#91bfdb')
        ax1.plot([x - int(g0) for x in g_list],v_5_list, label='SWC 5 cm',color='#4575b4')
        ax1.plot([x - int(g0) for x in g_list],np.ones(len(g_list))*fc, label='Field capacity',color='#fc8d59')
        ax1.plot([x - int(g0) for x in g_list],vkrit_Ze_list, label='Emergence threshold',color='#d73027')
        if sensor:
            ax1.plot([x - int(g0) for x in date_sensor],sensordata, label='Sensordata 5cm',color='k')
            ax1.fill_between([x - int(g0) for x in date_sensor], low_sensor, upper_sensor, alpha = 0.1,color='k',label='95% CI')
        ax1.errorbar([x - int(g0) for x in date_obs_grouped_0_5],obs0_5_mean,obs0_5_stdev, linestyle='None', marker='o', capsize=3, label='Soil sample 5 cm',color='k')
        # if opkomst: ax1.scatter([x - int(g0) for x in obs0_5_date],obs0_5,  marker='d',s=12, facecolors='none', edgecolors='gray', label='Soil sample 5 cm')
        if CI.size:        
            ax1.fill_between([x - int(g0) for x in g_list], CI[0], CI[1], alpha = 0.3,color='#4575b4',label='95% CI')
        if sensor_cal.size:
            ax1.plot([x - int(g0) for x in date_sensor],sensor_cal, label='Calibrated sensor data',color='k',linestyle='dashed')
        if df_forecast.size: ax1.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        ax2 = ax1.twinx()
        ax2.set_ylabel('Precipitation & irrigation (mm)')
        ax2.set_ylim(0, 80)
        ax2.bar([x - int(g0) -0.5 for x in g_list],I_list, label='Irrigation',color='#fee090',alpha=0.6)
        ax2.bar([x - int(g0) -0.5 for x in date_R_I],rain, label='Rainfall',color='#91bfdb',alpha=0.6)

        handles1, labels1 = ax1.get_legend_handles_labels()
        handles2, labels2 = ax2.get_legend_handles_labels()
        if  CI.size and sensor:
            comb_handles = [(handles1[0],handles1[5])]+[handles1[i] for i in [1,2]]+ \
                            [plt.Line2D([],[], alpha=0),(handles1[3],handles1[4]),handles1[-1]]+[handles2[i] for i in [0,1]]
            comb_labels = ['Soil water content 5cm\nwith 95% CI']+[labels1[i] for i in [1,2]]+ \
                            ['','Sensor data 5cm\nwith 95% CI',labels1[-1]]+[labels2[i] for i in [0,1]]
            fig.legend(comb_handles,comb_labels,loc='upper center',bbox_to_anchor=(0.5,1),
                        ncol=2,fancybox=True,shadow=True)
        elif sensor and not CI.size:
            comb_handles = [handles1[i] for i in [0,1,2]]+ \
                            [(handles1[3],handles1[4]),handles1[-1]]+[handles2[i] for i in [0,1]]
            comb_labels = [labels1[i] for i in [0,1,2]]+['Sensor data 5cm with 95% CI']+ \
                            [labels1[-1]]+[labels2[i] for i in [0,1]]
            fig.legend(comb_handles,comb_labels,loc='upper center',bbox_to_anchor=(0.5,1),
                        ncol=2,fancybox=True,shadow=True)            
        else: 
            fig.legend(loc='upper center',bbox_to_anchor=(0.5,1),ncol=2,fancybox=True,shadow=True)
        
        fig.tight_layout()                
        plt.subplots_adjust(top=0.9)
        if show and not save: plt.show()
        elif save: 
            plt.savefig(folder+'\\'+case+'_'+year+'_SWC_5.png',dpi=300)
            plt.close()            
        
        
        ### --------------------
        
        plt.figure(case+'_'+year+'_5-10-root-zone',figsize=(12,8))
        plt.subplot(221)
        plt.gca().set_title('Evaporative layer (5 & 10 cm)',y=1.1)
        plt.plot([x - int(g0) for x in g_list],sw_Ze_list, label='Soil water storage 10 cm',color='#91bfdb')
        plt.plot([x - int(g0) for x in g_list],sw_5_list, label='Soil water storage 5 cm',color='#4575b4')
        plt.plot([x - int(g0) for x in g_list],swfc_Ze_list, label='Field capacity 10 cm', linestyle='dotted',color='#fc8d59',alpha=0.5)
        plt.plot([x - int(g0) for x in g_list],swfc_5_list, label='Field capacity 5 cm', linestyle='dotted',color='#fc8d59')
        plt.plot([x - int(g0) for x in g_list],sw_Ze_vkrit_list, label='Emergence threshold 10 cm (pF 3)',color='#d73027',alpha=0.5)
        plt.plot([x - int(g0) for x in g_list],sw_5_vkrit_list, label='Emergence threshold 5 cm (pF 3)',color='#d73027')
        plt.bar([x - int(g0) for x in g_list],I_list, label='Irrigation',color='#fee090',alpha=0.4)
        plt.bar([x - int(g0) for x in date_R_I],rain, label='Rainfall',color='#91bfdb',alpha=0.4)
        if opkomst: 
            plt.scatter([x - int(g0) for x in obs0_5_date],obs0_5_asw,  marker='d',s=10, facecolors='none', edgecolors='gray', label='Soil sample 5 cm')
            plt.plot([x - int(g0) for x in date_sensor],sensordata_asw5, label='Sensordata 5cm',color='darkgray')
            plt.fill_between([x - int(g0) for x in date_sensor], low_sensor_asw5, upper_sensor_asw5, alpha = 0.1,color='darkgray',label='95% CI')
        # plt.errorbar([x - int(g0) for x in date_obs_grouped],obs_asw_mean, obs_asw_stdev, linestyle='None', marker='o', capsize=3, label='Soil water observation',color='k')
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        plt.ylabel('Available soil water (mm)\nPrecipitation & irrigation (mm)')
        plt.xlabel('Growing days')
        plt.legend(loc='upper center',bbox_to_anchor=(0.5,1.1),ncol=2,fontsize='small',fancybox=True,shadow=True)
        ax = plt.gca()
        ax.set(xlim=(0, end),xticks=range(0,end,10),ylim=(0,60))
        
        plt.subplot(222)
        plt.gca().set_title('Root zone layer',y=1.1)
        plt.plot([x - int(g0) for x in g_list],sw_list, label='Soil water storage',color='#4575b4')
        plt.plot([x - int(g0) for x in g_list],swfc_list, label='Field capacity',color='#fc8d59')
        plt.plot([x - int(g0) for x in g_list],sw_vkrit_list, label='Water stress threshold',color='#d73027')
        plt.errorbar([x - int(g0) for x in date_obs_grouped],obs_asw_mean, obs_asw_stdev, linestyle='None', marker='o', capsize=3, label='Soil water observation',color='k')
        plt.bar([x - int(g0) for x in g_list],I_list, label='Irrigation',color='#fee090',alpha=0.6)
        plt.bar([x - int(g0) for x in date_R_I],rain, label='Rainfall',color='#91bfdb',alpha=0.6)
        # plt.scatter([x - int(g0) for x in date_obs_grouped],obs_asw_mean, marker='o', label='Soil water observation',color='k')
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        plt.ylabel('Available soil water (mm)\nPrecipitation & irrigation (mm)')
        plt.xlabel('Growing days')
        plt.legend(loc='upper center',bbox_to_anchor=(0.5,1.1),ncol=3,fontsize='small',fancybox=True,shadow=True)
        ax = plt.gca()
        ax.set(xlim=(0, end),xticks=range(0,end,10))
        
        plt.subplot(223)
        ax1 = plt.gca()
        ax1.set_xlabel('Growing days')
        ax1.set_ylabel('Soil water content (m³/m³)')
        ax1.set_ylim(0, 0.5)
        ax1.plot([x - int(g0) for x in g_list],v_Ze_list, label='SWC 10 cm',color='#91bfdb')
        ax1.plot([x - int(g0) for x in g_list],v_5_list, label='SWC 5 cm',color='#4575b4')
        ax1.plot([x - int(g0) for x in g_list],np.ones(len(g_list))*fc, label='Field capacity',color='#fc8d59')
        ax1.plot([x - int(g0) for x in g_list],vkrit_Ze_list, label='Emergence threshold (pF 3)',color='#d73027')
        if sensor and opkomst:
            ax1.plot([x - int(g0) for x in date_sensor],sensordata, label='Sensordata 5cm',color='darkgray')
            ax1.fill_between([x - int(g0) for x in date_sensor], low_sensor, upper_sensor, alpha = 0.1,color='darkgray',label='95% CI')
        # ax1.plot([x - int(g0) for x in g_list],TEW_swc, label='Evaporation threshold (TEW)',linestyle='dashed',color='#d73027') #TEW
        # ax1.plot([x - int(g0) for x in g_list],REW_swc, label='REW',linestyle='dotted',color='#d73027') #REW
        if opkomst: ax1.scatter([x - int(g0) for x in obs0_5_date],obs0_5,  marker='d',s=10, facecolors='none', edgecolors='gray', label='Soil sample 5 cm')
        ax2 = ax1.twinx()
        ax2.set_ylabel('Precipitation & irrigation (mm)')
        ax2.set_ylim(0, 80)
        ax2.bar([x - int(g0) for x in g_list],I_list, label='Irrigation',color='#fee090',alpha=0.4)
        ax2.bar([x - int(g0) for x in date_R_I],rain, label='Rainfall',color='#91bfdb',alpha=0.4)
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        ax1.set(xlim=(0, end),xticks=range(0,end,10))
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        plt.legend(lines + lines2, labels + labels2,loc='upper center',bbox_to_anchor=(0.5,1.1),ncol=3,fontsize='small',fancybox=True,shadow=True)
        
        
        plt.subplot(224)
        ax1 = plt.gca()
        ax1.set_xlabel('Growing days')
        ax1.set_ylabel('Soil water content (m³/m³)')
        ax1.set_ylim(0, 0.5)
        ax1.plot([x - int(g0) for x in g_list],v_list, label='SWC',color='#4575b4')
        ax1.plot([x - int(g0) for x in g_list],np.ones(len(g_list))*fc, label='Field capacity',color='#fc8d59')
        ax1.plot([x - int(g0) for x in g_list],vkrit_list, label='Water stress threshold',color='#d73027')
        if sensor:
            if opkomst==False:
                ax1.plot([x - int(g0) for x in date_sensor],sensordata, label='Sensordata',color='darkgray')
                ax1.fill_between([x - int(g0) for x in date_sensor], low_sensor, upper_sensor, alpha = 0.3,color='darkgray',label='95% CI')
        ax2 = ax1.twinx()
        ax2.set_ylabel('Precipitation & irrigation (mm)')
        ax2.set_ylim(0, 80)
        ax2.bar([x - int(g0) for x in g_list],I_list, label='Irrigation',color='#fee090',alpha=0.6)
        ax2.bar([x - int(g0) for x in date_R_I],rain, label='Rainfall',color='#91bfdb',alpha=0.6)
        # plt.scatter([x - int(g0) for x in date_obs_grouped],obs_asw_mean, marker='o', label='Soil water observation',color='k')
        ax1.errorbar([x - int(g0) for x in date_obs_grouped],obs_mean, obs_stdev, linestyle='None', marker='o', capsize=3, label='Soil water observation',color='k')
        if df_forecast.size: plt.axvline(forecast_date - int(g0), color='lightgray', linestyle='dotted')
        ax1.set(xlim=(0, end),xticks=range(0,end,10))
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        plt.legend(lines + lines2, labels + labels2,loc='upper center',bbox_to_anchor=(0.5,1.1),ncol=3,fontsize='small',fancybox=True,shadow=True)
        plt.tight_layout()
        plt.subplots_adjust(top=0.9)
        plt.suptitle(case+'_'+year,ha='center',weight='bold')
        if show and not save: plt.show()
        elif save: 
            plt.savefig(folder+'\\'+case+'_'+year+'_5-10-root-zone.png',dpi=300)
            plt.close()
            
            
    ## Remove NaN within time series (sensor data & covar matrix)
    if sensor: 
        i_notnan=np.argwhere(~np.isnan(np.array(sensor_data_adj['Mean'])))[:,0]
        i_nan=np.argwhere(np.isnan(np.array(sensor_data_adj['Mean'])))[:,0]
        sensor_data_adj = sensor_data_adj.iloc[i_notnan]
        sensor_data_adj.reset_index(drop=True,inplace=True)
        covar_adj = np.delete(covar_adj,i_nan, axis=0)
        covar_adj = np.delete(covar_adj,i_nan, axis=1)
        # print(sensor_data_adj)
        
    if part=='1': df_obs_grouped=df_obs_grouped[df_obs_grouped['Date']<=enddate]
    # print(df_obs_grouped)
    
    ## summary dataframe all data
    if show:
        df_K_ET = pd.DataFrame({'Kcmax':kcmax[:end],'Kcb': kcb_list,'Kcb_act':kcb_a_list,'Ks':Ks_list,
                                'Ke':ke_list,'Kr':kr_list,'Kc':kc_list,'Kc_act':kc_a_list,
                                'Zr (m)':zr[:end],'Zr virtual (m)':zvirt_list[:end],'Dr (mm)':Dr_list,'De (mm)':De_list,
                                'ETcb (mm)':ETcb_list,'ETcb_act (mm)':ETcb_a_list,
                                'ETe (mm)':ETe_list,'ETm (mm)':etm_list,'ET_act (mm)':eta_list})
        df_krit_fc = pd.DataFrame({'vkrit':vkrit_list,'vkrit_Ze': vkrit_Ze_list,'ASW vkrit (mm)':sw_vkrit_list,
                                'ASW vkrit_Ze (mm)':sw_Ze_vkrit_list,'ASW vkrit_5 (mm)':sw_5_vkrit_list,
                                'FC (mm)':swfc_list,'FC_Ze (mm)':swfc_Ze_list,'FC_5 (mm)':swfc_5_list})
        df_balance = pd.DataFrame({'SWC':v_list,'SWC_Ze':v_Ze_list,'SWC_5':v_5_list,
                                   'ASW (mm)':sw_list,'ASW_Ze (mm)':sw_Ze_list,'ASW_5 (mm)':sw_5_list,
                                   'Percolation (mm)':p_list,'CR (mm)':cri_list,
                                   'Runoff (mm)':runoff_list,'Precipitation (mm)':R_list,'Irrigation (mm)':I_list,
                                   'TAW (mm)':TAW_list,'RAW (mm)':RAW_list})
        
        df_all = pd.concat([df_K_ET, df_krit_fc, df_balance], axis=1, join='inner')
    
    else:
        df_all = pd.DataFrame({'vkrit':vkrit_list,'vkrit_Ze': vkrit_Ze_list,
                               'ASW vkrit (mm)':sw_vkrit_list, 'ASW vkrit_5 (mm)':sw_5_vkrit_list})
    
    if sensor:
        df_parameters = pd.concat([df_pars,df_calpar], axis=1, join='inner')
    else:
        df_parameters = df_pars.copy()
    
    if opkomst:
        SWC=v_5_list
        sw_list=sw_5_list
        df_obs_grouped=df_obs_grouped_0_5
    
    #runtime
    print(datetime.now()-start)
    
    return SWC, sw_list, g_list, sensor_data_adj, covar_adj, df_obs_grouped, df_all, df_parameters

