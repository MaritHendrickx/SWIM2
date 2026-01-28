# SWIM2

_The SWIM² framework has been developed by Marit G. A. Hendrickx (KU Leuven), under guidance of Jan Diels, Jan Vanderborght and Pieter Janssens, and in collaboration with Eric Laloy, Sander Bombeke, Evi Matthyssen and Anne Waverijn._

_Reference: Hendrickx, M. G. A., Vanderborght, J., Janssens, P., Laloy, E., Bombeke, S., Matthyssen, E., et al. (2026). Field‐scale soil moisture predictions in real time using in situ sensor measurements in an inverse modeling framework: SWIM2. Water
Resources Research, 62, e2025WR041324. https://doi.org/10.1029/2025WR041324_

All scripts should be in the same directory. The data folders should also be in this directory. If this is not the case, mind that you adjust the paths.

How to run it:
- Open your prompt, activate your environment.
- Go to the correct directory, example:
``cd C:/Users/u0139791/OneDrive - KU Leuven/Documents/SWIM2``
- Run the DREAM script:
``python DREAM_auto.py``

## Model script (SWB_model.py)
In this script, there are several #TODO indicated, where you can change the settings.
#### drop_samp
- You can drop a sample by entering the index, e.g., 1, [2], …, -1
- if type(drop_samp)==int: drop starting with that row ``[drop_samp:]``
- if type(drop_samp)==list: drop only that row ``[drop_samp[0]]``

#### validation_days
``np.nan`` if using it in realtime with weather forecasts
#### cal_par_on
Estimate calibration parameters in DREAM yes or no?
#### lik_sigma_est
Estimate sensor error (co)variances in DREAM yes or no?
#### DREAM_obsdata
Options: 
-	'Sensor+stalen' : use sensor data and samples
-	'Sensor only' : use sensor data only
-	'Samples only' : use samples only
-	'Sensor+2samp' : use sensor data and first two samples
#### cal
Define which sensor calibration is used.
- 'gen' = general / pooled sensor calibration (Hendrickx et al., 2025)
- 'WTLS' = Weighted total least squares (Deming) with ‘intc0’ (zero intercept) or ‘slope1’ (slope of 1) or ‘slope_gen’ (slope of general calibration)
- '' = no sensor calibration

These different approaches are programmed in Sensordata.py
#### zero_cov
Set error covariances to zero, yes or no?

## MCMC script (mcmc.py)
Nothing needs changing here. Everything is selected automatically based on choices made in SWB_model.py and DREAM_auto.py. CaseStudy==1 is current framework. Different priors can be defined here. 

## MCMC functions (mcmc_func.py)
Nothing needs changing here. 
soil_water_model  is a function at the end of this script where the SWB model is run with the parameter set of one iteration. The output (fx) contains the model output that is comparable with the observation data (sensor data and/or samples).

## Sensor data script (Sensordata.py)
Nothing needs changing here. 
Here, you can adjust the file name (the csv with the data), or make adjustments depending on the format of your data.
sensordata  is a function which reads the sensor data, calibrates the sensor data, and creates the sensor error covariance matrix.

## DREAM automation script (DREAM_auto.py)
In this script, there are several #TODO indicated, where you can change the settings.
