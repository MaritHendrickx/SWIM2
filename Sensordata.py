# -*- coding: utf-8 -*-
''' 
Sensor data to feed the SWIM² framework

Creator: Marit Hendrickx (KU Leuven)
Last update: September 2025

'''

import matplotlib.pyplot as plt
from math import *
import numpy as np

from datetime import datetime
import pandas as pd
import scipy.stats as stats

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import os

#### CONVERT RAW SENSOR DATA (METER equation) ####
def sensorEQ(mV):
    ''' Equation to calculate the SWC from the mV sensor measurements '''
    SWC=4.824*10**(-10)*mV**3-2.278*10**(-6)*mV**2+3.898*10**(-3)*mV-2.154
    return SWC

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


def sensordata(date_YYYYMMDD,case,year,g0,date_obs_grouped,obs_mean,obs_stdev,cal,show,drop_samp=np.nan):
    ''' 
    Input:
        Date of the (virtual) 'today', case reference number, year, planting date g0
        Soil samples: dates, mean and stdev
    '''
    main_dir = os.path.dirname(os.path.realpath(__file__))
    os.chdir(main_dir)

    show_vec = show
    if show_vec[1]=='':
        save=False
    else: 
        save=True
        folder=show_vec[1]
    show=show_vec[0] 

    df = pd.read_csv(year+'\\sensordata\\dataViaAPI_'+case+'.csv')
    df.sort_values(by='Datetime',inplace=True)
    df.reset_index(drop=True,inplace=True)

    df_sensor=pd.DataFrame([datetime.strptime(df['Datetime'][i], '%Y-%m-%d %H:%M:%S.%f').strftime('%Y-%m-%d') for i in range(len(df))],columns=['Date'])
    df_sensor['Serial date']=[ConvertToSerialDate(datetime.strptime(df['Datetime'][i], '%Y-%m-%d %H:%M:%S.%f')) for i in range(len(df))]
    df_sensor['Datetime']=[datetime.strptime(df['Datetime'][i], '%Y-%m-%d %H:%M:%S.%f') for i in range(len(df))]
    
    df_sensor['Sensor 0 (m³/m³)']=sensorEQ(np.array(df['Adc0 (mV)']))
    df_sensor['Sensor 1 (m³/m³)']=sensorEQ(np.array(df['Adc1 (mV)']))
    df_sensor['Sensor 2 (m³/m³)']=sensorEQ(np.array(df['Adc2 (mV)']))

    df_sensor.loc[~(df_sensor['Sensor 0 (m³/m³)'] > 0.01), 'Sensor 0 (m³/m³)']=np.nan
    df_sensor.loc[(df_sensor['Sensor 0 (m³/m³)'] > 1), 'Sensor 0 (m³/m³)']=np.nan
    df_sensor.loc[~(df_sensor['Sensor 1 (m³/m³)'] > 0.01), 'Sensor 1 (m³/m³)']=np.nan
    df_sensor.loc[(df_sensor['Sensor 1 (m³/m³)'] > 1), 'Sensor 1 (m³/m³)']=np.nan
    df_sensor.loc[~(df_sensor['Sensor 2 (m³/m³)'] > 0.01), 'Sensor 2 (m³/m³)']=np.nan
    df_sensor.loc[(df_sensor['Sensor 2 (m³/m³)'] > 1), 'Sensor 2 (m³/m³)']=np.nan

    df_sensor=df_sensor[df_sensor['Serial date']>=g0] 
    today=np.floor(ConvertToSerialDate(datetime.strptime(date_YYYYMMDD,'%Y%m%d')))
    df_sensor=df_sensor[df_sensor['Serial date']<=today]
    df_sensor.reset_index(drop=True,inplace=True)
    # print(df_sensor)
    
    #### Daily sensor values        
    dates = df_sensor.Date.unique()
    daily = []
    median_daily = []
    start_daily = []
    end_daily = []
    for i in range(len(dates)):
        a = df_sensor.loc[(df_sensor['Date'] == dates[i])]
        mean = np.mean(a,axis=0)[-3:].tolist()
        median = np.median(a.iloc[:,-3:],axis=0).tolist()
        start = a.iloc[0,-3:]
        end = a.iloc[-1,-3:]
        daily.append(mean)
        median_daily.append(median)
        start_daily.append(start)
        end_daily.append(end)
    daily = np.array(daily)
    median_daily = np.array(median_daily)
    start_daily = np.array(start_daily)
    end_daily = np.array(end_daily)
    
    dates = [ConvertToSerialDate(datetime.strptime(dates[i], '%Y-%m-%d')) for i in range(len(dates))]
    
    df_daily = pd.DataFrame({'Date': dates,
                             'Mean 0 (m³/m³)': daily[:,0].tolist(),
                             'Mean 1 (m³/m³)': daily[:,1].tolist(), 
                             'Mean 2 (m³/m³)': daily[:,2].tolist(),
                             'Median 0 (m³/m³)': median_daily[:,0].tolist(),
                             'Median 1 (m³/m³)': median_daily[:,1].tolist(), 
                             'Median 2 (m³/m³)': median_daily[:,2].tolist(),
                             'Start 0 (m³/m³)': start_daily[:,0].tolist(),
                             'Start 1 (m³/m³)': start_daily[:,1].tolist(), 
                             'Start 2 (m³/m³)': start_daily[:,2].tolist(),
                             'End 0 (m³/m³)': end_daily[:,0].tolist(),
                             'End 1 (m³/m³)': end_daily[:,1].tolist(), 
                             'End 2 (m³/m³)': end_daily[:,2].tolist(),
                             })
    
    
    daily='end'
    #Based on daily median or start of the day (is end day)
    if daily=='mean': 
        i=1
        j=4
    if daily=='median': 
        i=4
        j=7
    if daily=='start':
        i=7 
        j=10
    if daily=='end':
        i=10
        j=13
        
        
    # ----------------------------------------------------------    
    #### Soil samples - matching

    # Array with the daily means that match soil samples
    daily_obs=[]
    for x in range(len(date_obs_grouped)):
        if np.any(df_daily['Date'] == date_obs_grouped[x]):
            daily_obs.append(df_daily.loc[(df_daily['Date'] == date_obs_grouped[x])].iloc[0,i:j].tolist())
        else: daily_obs.append(np.full(3, np.nan).tolist())
    try:
        daily_obs=np.stack(daily_obs) #.reshape((len(date_obs),4))
    except:
        daily_obs=[daily_obs]

    if show or save:
        plt.figure('Original measurements')
        plt.plot(df_sensor.iloc[:,1],df_sensor.iloc[:,-3:])
        plt.scatter(date_obs_grouped,obs_mean,color='k',s=15)
        # plt.errorbar(date_obs_grouped,obs_mean,obs_stdev,color='k')
        # for x in range(3): plt.scatter(date_obs_grouped,daily_obs[:,x])
        plt.ylabel('Soil water content (m³/m³)')
        plt.xlim(left=g0-1)
        plt.legend(['Sensor 0','Sensor 1','Sensor 2','Soil sample'])
        if show and not save: plt.show()
        elif save: 
            plt.savefig(folder+'\\Original measurements.png',dpi=300)
            plt.close()    
            
    # ----------------------------------------------------------
    ### Mean sensor value & stdev

    sensor_mean = np.array(df_daily.iloc[:,i:j].mean(axis=1,numeric_only=True))
    sensor_stdev = np.array(df_daily.iloc[:,i:j].std(axis=1,numeric_only=True))
    sensor_dates = np.array([int(v) for v in df_daily['Date']])
        
            
    # print(daily_obs)
    # I expect to see Warnings in this block
    with warnings.catch_warnings():
        warnings.filterwarnings(action='ignore', message='Mean of empty slice')
        warnings.filterwarnings(action='ignore', message='Degrees of freedom <= 0 for slice')
        x_mean=np.nanmean(daily_obs,axis=1) # sensor
        x_var=np.nanvar(daily_obs,axis=1)   # sensor var
    y_mean=np.array(obs_mean)           # samples #np.tile(obs_5,3)
    y_var=np.array(obs_stdev)**2        # np.ones(len(y))*0.0001 or =y_var
    if np.all((np.array(y_var) == 0)): y_var=x_var  # assumption of equal uncertainty
    elif np.any((np.array(y_var) == 0)): 
        zero = np.where(y_var==0)[0].tolist()
        for v in range(len(zero)): y_var[zero[v]]=x_var[zero[v]]
        
    y=np.delete(y_mean,np.argwhere(np.isnan(x_mean)+np.isnan(y_mean)))
    x=np.delete(x_mean,np.argwhere(np.isnan(x_mean)+np.isnan(y_mean)))
    x_var=np.delete(x_var,np.argwhere(np.isnan(x_mean)+np.isnan(y_mean)))
    y_var=np.delete(y_var,np.argwhere(np.isnan(x_mean)+np.isnan(y_mean)))
    
    ## drop samples
    if not isnan(drop_samp):
        if type(drop_samp)==int:
            if len(x)>=drop_samp:
                x = x[:drop_samp]
                y = y[:drop_samp]
                x_var = x_var[:drop_samp]
                y_var = y_var[:drop_samp]
            else: print(f"Warning: Not enough samples.")
        elif type(drop_samp)==list:
            if len(x)>=drop_samp[0]+1:
                x = np.delete(x, drop_samp[0])
                y = np.delete(y, drop_samp[0])
                x_var = np.delete(x_var, drop_samp[0])
                y_var = np.delete(y_var, drop_samp[0])
            else: print(f"Warning: Not enough samples.")
            
    #### Sensor data calibration
    if not 'slope1' in cal and not 'slope_gen' in cal and not 'intc0' in cal:
        if isnan(drop_samp) and (len(y)<=1 or len(x)<=1): cal='gen'
    else:
        if isnan(drop_samp) and (len(y)==0 or len(x)==0): cal='gen'    

    if cal == 'gen': # general calibration
        ## Pooled sensor calibration parameters (Hendrickx et al., 2025)
        a = -0.006
        b = 1.260
        
        gen_cal = a+b*df_daily.iloc[:,i:j]
        sensor_mean_gen = np.array(gen_cal.mean(axis=1,numeric_only=True))
        sensor_stdev_gen = np.array(gen_cal.std(axis=1,numeric_only=True))
        
        if show or save:        
            plt.figure('Mean daily ('+daily+' of the day) sensor data',figsize=(12,6))
            plt.plot(sensor_dates,sensor_mean, label='Mean sensor')
            plt.fill_between(sensor_dates,sensor_mean-sensor_stdev,sensor_mean+sensor_stdev,alpha=0.3, label='Standard deviation')
            plt.plot(sensor_dates,sensor_mean_gen, label='Calibrated sensor')
            plt.fill_between(sensor_dates,sensor_mean_gen-sensor_stdev_gen,sensor_mean_gen+sensor_stdev_gen,alpha=0.3, label='Cal. standard deviation')
            plt.scatter(date_obs_grouped,obs_mean,color='k',s=15,label='Soil sample')
            plt.xlabel('Day')
            plt.ylabel('Soil water content (m³/m³)')
            plt.xlim(left=g0-1)
            plt.ylim(0,0.45)
            plt.legend()
            if show and not save: plt.show()
            elif save: 
                plt.savefig(folder+'\\Mean daily ('+daily+' of the day) sensor data.png',dpi=300)
                plt.close()
                
        sensor_mean=sensor_mean_gen
        sensor_stdev=sensor_stdev_gen
        
    elif cal == '':
        a=0
        b=1      
        
    else:
        rho = np.cov(x,y)[0][1]/np.sqrt(np.cov(x,y)[0][0]*np.cov(x,y)[1][1])
        if show:
            print('Correlation coefficient =', rho)
            print('R² =', rho**2)
        
        if 'WTLS' in cal:
            ''' Weighted total least squares (Deming)''' # y = a + b*x
            add_zero=True
            if 'slope1' in cal or 'slope_gen' in cal or 'intc0' in cal: add_zero=False
            if add_zero:
                # Add zero
                if show: print('(0,0) point added')
                # add intercept:
                x=np.append(x,0)
                y=np.append(y,0)
                x_var=np.append(x_var,0)
                y_var=np.append(y_var,0)
                rho = np.cov(x,y)[0][1]/np.sqrt(np.cov(x,y)[0][0]*np.cov(x,y)[1][1])
                if show: 
                    print('Correlation coefficient =', rho)
                    print('R² =', rho**2)
                    print(x,y)
                     
            # Calculation
            delta = np.array(1)     # Defines angle of uncertainy (delta = 1: orthogonal)
            cov = np.cov(x.reshape(-1),y.reshape(-1))
            mean_x = np.nanmean(x)
            mean_y = np.nanmean(y)
            s_xx = cov[0,0]
            s_yy = cov[1,1]
            s_xy = cov[0,1]
    
            if 'slope1' in cal: # estimate intercept
                b=1
                a = mean_y - b* mean_x 
            elif 'slope_gen' in cal: # estimate intercept with the slope of the general calibration
                b = 1.260
                a = mean_y - b* mean_x 
            elif 'intc0' in cal: # estimate slope
                a=0
                b = mean_y/mean_x 
            else:
                b = (s_yy  - delta * s_xx + np.sqrt((s_yy - delta * s_xx) ** 2 + 4 * delta * s_xy ** 2)) / (2 * s_xy)
                a = mean_y - b* mean_x 
    
            ## Calculation Fuller (1987)
            # Regression model variance sigma2
            N = len(x)
            sigma2 = np.sum((y-mean_y-b*(x-mean_x))**2)/(N-2)
            sigma_xx = (np.sqrt((s_yy-delta*s_xx)**2 + 4*delta*s_xy**2) - (s_yy - delta*s_xx)) / (2*delta)
            sigma_uu = (s_yy + delta*s_xx - np.sqrt((s_yy-delta*s_xx)**2 + 4*delta*s_xy**2)) / (2*delta)
            if 'slope1' in cal or 'slope_gen' in cal: # estimate intercept
                var_b = 0
                var_a = sigma2/N
            elif 'intc0' in cal: # estimate slope
                var_a = 0
                var_b = abs(-(sigma2/N)/mean_x**2)
            else:
                var_b = (s_xx*sigma2 - b**2*sigma_uu**2)/((N-1)*sigma_xx**2)
                var_a = sigma2/N + mean_x**2*var_b
                
            covar_ab = -mean_x*var_b
            corr_ab = covar_ab/np.sqrt(var_a*var_b)
            if show:
                print('Standard deviation \n a: ',round(np.sqrt(var_a),6),
                  '\n b: ',round(np.sqrt(var_b),6))
                print('Correlation cal. param. = ',round(corr_ab,4))
                print('y =', np.round(a,3),'+', np.round(b,3),'x')
    
            t = stats.t.ppf(1-0.05/2,N-2)
            x_range=y_range=np.arange(0,0.51,0.01)
            if 'slope1' in cal or 'slope_gen' in cal: # estimate intercept
                se2_ygivenxrange = var_a
                se2_ygivenxrange_pred = var_a*(1+N) # = var_a+sigma2
            elif 'intc0' in cal: # estimate slope
                se2_ygivenxrange = var_b*(x_range)**2
                se2_ygivenxrange_pred = se2_ygivenxrange+sigma2
            else:
                se2_ygivenxrange = sigma2/N*(1+(x_range-mean_x)**2/s_xx)
                se2_ygivenxrange_pred = sigma2/N*(1+N+(x_range-mean_x)**2/s_xx)

            if show or save:
                # for x=sensor
                plt.figure('Calibrated sensordata: a = '+str(round(a,2))+', b = '+str(round(b,2)))
                plt.plot(df_sensor.iloc[:,1],a+b*df_sensor.iloc[:,3:6])
                plt.scatter(date_obs_grouped,obs_mean,color='k',s=15)
                for v in range(3): plt.scatter(date_obs_grouped,a+b*daily_obs[:,v])
                plt.xlabel('Day')
                plt.ylabel('Soil water content (m³/m³)')
                plt.xlim(left=g0-1)
                plt.ylim(bottom=0)
                plt.legend(['Sensor 0','Sensor 1','Sensor 2','Soil sample'])
                if show and not save: plt.show()
                elif save: 
                    plt.savefig(folder+'\\Calibrated sensordata_a_=_'+str(round(a,2))+'_b_=_'+str(round(b,2))+'.png',dpi=300)
                    plt.close() 
                        
                plt.figure('Calibration')
                plt.plot(x_range,a+b*x_range)
                # confidence interval
                plt.fill_between(x_range,a+b*x_range-t*np.sqrt(se2_ygivenxrange),a+b*x_range+t*np.sqrt(se2_ygivenxrange),color='gray',alpha=0.3, label='95% CI')
                # prediction interval
                plt.fill_between(x_range,a+b*x_range-t*np.sqrt(se2_ygivenxrange_pred),a+b*x_range+t*np.sqrt(se2_ygivenxrange_pred),color='lightgray',alpha=0.3, label='95% PI')
                plt.scatter(x,y,color='k',s=15)
                plt.ylabel('Soil water content (m³/m³)')
                plt.ylim((0,0.45))
                plt.xlabel('Sensor measurement (m³/m³)')
                plt.xlim((0,0.45))
                plt.title('Calibration: \n y = '+str(round(a,2))+' + '+str(round(b,2))+'x')
                if show and not save: plt.show()
                elif save: 
                    plt.savefig(folder+'\\Calibration.png',dpi=300)
                    plt.close()    

            wtls_cal = a+b*df_daily.iloc[:,i:j]
            sensor_mean_wtls = np.array(wtls_cal.mean(axis=1,numeric_only=True))
            sensor_stdev_wtls = np.array(wtls_cal.std(axis=1,numeric_only=True))            
        
            if show or save: 
                plt.figure('Mean daily ('+daily+' of the day) calibrated sensor data',figsize=(12,6))
                plt.plot(sensor_dates,sensor_mean, label='Mean sensor')
                plt.fill_between(sensor_dates,sensor_mean-sensor_stdev,sensor_mean+sensor_stdev,alpha=0.3, label='Standard deviation')
                plt.plot(sensor_dates,sensor_mean_wtls, label='Calibrated sensor')
                plt.fill_between(sensor_dates,sensor_mean_wtls-sensor_stdev_wtls,sensor_mean_wtls+sensor_stdev_wtls,alpha=0.3, label='Cal. standard deviation')
                plt.scatter(date_obs_grouped,obs_mean,color='k',s=15,label='Soil sample')
                plt.xlabel('Day')
                plt.ylabel('Soil water content (m³/m³)')
                plt.xlim(left=g0-1)
                plt.ylim(0,0.45)
                plt.title('Mean daily ('+daily+' of the day) sensor data')
                plt.legend()
                if show and not save: plt.show()
                elif save: 
                    plt.savefig(folder+'\\Mean daily ('+daily+' of the day) calibrated sensor data.png',dpi=300)
                    plt.close()
            
            sensor_mean=sensor_mean_wtls
            sensor_stdev=sensor_stdev_wtls
    
    sensor_data = pd.DataFrame({'Date':sensor_dates,'Mean':sensor_mean,'Stdev':sensor_stdev})
    df_calpar = pd.DataFrame({'intercept a':a,'slope b':b}, index=[0])


    #### COVARIANCE MATRIX ----------------------------------------------------

    ## Covariance matrix based on pooled error model
    # Paper by Hendrickx et al. (2025)
    rand_error_var = 0.000998
    syst_error_covar = 0.001070
    
    covar = np.ones((len(sensor_mean),len(sensor_mean)))*syst_error_covar
    np.fill_diagonal(covar, np.diag(covar) + rand_error_var)
        
    
    return sensor_data, covar, df_calpar


