#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json 
import scipy.integrate as spi
from scipy.interpolate import interp1d  

from algorithmos_DTU import Hansen_Algorithm

       
 
#%%
if __name__ == "__main__":
    blade_geom_file_2 = "blade_geom_DTU.json"
    hansen_DTU = Hansen_Algorithm(
        blade_geom_DTU=blade_geom_file_2,
        B=3,
        air_density=1.225,
        csv_data_file='csv_data_file_DTU.csv'
    )


    
#%%    
    # ΔΙΑΓΡΑΜΜΑ Power Coefficient Cp - Tip Speed Ratio λ for DTU geometry
    wind_speed_V0=10
    rotation_speed_values = np.linspace(0, 1.3, 50)
    lambda_values = []
    cp_values = []
    
    for rotational_speed in rotation_speed_values:
        λ = (rotational_speed * hansen_DTU.R) / wind_speed_V0
        lambda_values.append(λ)
        
        results_for_Naca_geometry, total_power, total_torque, total_thrust = hansen_DTU.DTU_blade_calculation(wind_speed_V0=wind_speed_V0, rotation_speed=rotational_speed)
        cp = hansen_DTU.calculation_of_coefficient_of_power_cp_for_DTU(total_power, wind_speed_V0=wind_speed_V0)
        cp_values.append(cp)
        
    plt.figure(figsize=(10, 8))
    plt.plot(lambda_values, cp_values, 'o-', label="$C_p$ vs $λ$")
    plt.xlabel("$λ$ (Tip speed ratio)", fontsize=15)
    plt.ylabel("$C_p$ (Power coefficient)", fontsize=15)
    plt.title("Coefficient of Power $C_p$ as a function of Tip Speed Ratio $λ$ for DTU geometry")
    plt.legend()
    plt.grid()
    # plt.show()
#%%    
    # ΔΙΑΓΡΑΜΜΑ ταχύτητας n(rpm) - ισχύος P(kWatt) for DTU geometry
    rpm_values = np.linspace(0, 25, 50) # Τιμές RPM
    wind_speed_values = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15] # Ταχύτητες ανέμου

    plt.figure(figsize=(10, 6))

    for V0 in wind_speed_values:
        power_values = []
        for n in rpm_values:
            ω = (2 * np.pi * n) / 60 # Γωνιακή ταχύτητα σε rad/s
            λ = (ω * hansen_DTU.R) / V0 # Tip Speed Ratio
            
            results_for_DTU_geometry, total_power, total_torque, total_thrust = hansen_DTU.DTU_blade_calculation(wind_speed_V0=V0, rotation_speed=ω)
            power_values.append(total_power*1e-3) # Ισχύς σε kW
            
        plt.plot(rpm_values, power_values,'o', label=f"{V0} m/s", )

    plt.xlabel("Rotational Speed (RPM)", fontsize=12)
    plt.ylabel("Power (kW)", fontsize=12)
    plt.title("Power vs Rotational Speed for Different Wind Speeds for DTU geometry", fontsize=14)
    plt.legend(title="Wind Speed [m/s]")
    plt.grid()
    plt.show()
    
    df_DTU_results = pd.DataFrame(results_for_DTU_geometry)
    print(df_DTU_results)
    
    print(f"Η συνολική ισχύς της ανεμογεννήτριας είναι {total_power:.2f} Watt")
    print(f"H συνολική ροπή της ανεμογεννήτριας είναι {total_torque:.2f} Nm")
    print(f"H συνολική ώση της ανεμογεννήτριας είναι {total_thrust:.2f} N")
    # df_DTU_results.to_csv("output_tab.csv", sep="\t", index=False)
#%%
if __name__== "__main__":

    V0=10
    rotation_speed =0.5
    pitch=hansen_DTU.pitch[0]
    chordn= hansen_DTU.chords[0]
    tc_ratio = hansen_DTU.tc_ratios[0]
    hansen_DTU.segment_calculation(
        wind_speed_V0=V0,
        omega_rad_sec=rotation_speed,
           r=2.8, 
            chord=chordn,
            pitch_angle_deg=pitch,
            twist_deg= 0,
            tc_ratio=tc_ratio,
            f=0.3
            )

    # %%
    hansen_DTU.pitch

# %%
