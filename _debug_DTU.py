#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json 
import scipy.integrate as spi
from scipy.interpolate import interp1d  

from algorithmos_DTU import Hansen_Algorithm
%load_ext autoreload
%autoreload 2
       
#%% 
blade_geom_file_2 = "blade_geom_DTU.json"
hansen_DTU = Hansen_Algorithm(
    blade_geom_DTU=blade_geom_file_2,
    B=3,
    air_density=1.225,
    csv_data_file='csv_data_file_DTU.csv'
)




# ΔΙΑΓΡΑΜΜΑ Power Coefficient Cp - Tip Speed Ratio λ for DTU geometry
wind_speed_V0=10
rotation_speed_values = np.linspace(0, 1.3, 50)


V0=10
rotation_speed =0.5
pitch=hansen_DTU.pitch[0]
chordn= hansen_DTU.chords[0]
tc_ratio = hansen_DTU.tc_ratios[0]
res_dict = hansen_DTU.segment_calculation(
    wind_speed_V0=V0,
    omega_rad_sec=rotation_speed,
        r=2.8, 
        chord=chordn,
        pitch_angle_deg=pitch,
        twist_deg= 0,
        tc_ratio=tc_ratio,
        f=0.3
        )

for key, value in res_dict.items():
    print(f"{key:20s}: {value}")
# %%
hansen_DTU.pitch

# %%
