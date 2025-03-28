#%%
import numpy as np
import pandas as pd 
import json
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from algorithmos import Hansen_Algorithm

#%%
class Hansen_Algorithm_for_DTU_geometry(Hansen_Algorithm): # Νέα κλάση για τη γεωμετρία του πτερυγίου (DTU airfoil)
    def __init__(self, wind_speed_V0, rotation_speed, blade_geom_file_2, B=3, air_density=1.225, csv_data_file_DTU='csv_data_file_DTU.csv'):
        with open(blade_geom_file_2, 'r') as f:
            blade_geom_file_2 = json.load(f)
        self.r_is = blade_geom_file_2["r_is"] # πίνακας με την ακτινική θέση του κάθε τμήματος του πτερυγίου
        self.chords = blade_geom_file_2["chords"] # πίνακας με το μήκος χορδής του κάθε τμήματος του πτερυγίου
        self.pitch = blade_geom_file_2["pitch"] # πίνακας με τη γωνία βήματος του κάθε τμήματος του πτερυγίου
        self.tc_ratios = blade_geom_file_2["tc_ratios"] # πίνακας με το λόγο πάχους / χορδή αεροτομής (t/c ratio) του κάθε τμήματος του πτερυγίου
        self.no_sections = blade_geom_file_2["no_sections"] # αριθμός τμημάτων του πτερυγίου
        R = blade_geom_file_2["R"] # ακτίνα του ρότορα
        super().__init__(wind_speed_V0, R, rotation_speed, B, air_density, airfoil_type="DTU", csv_data_file=csv_data_file_DTU)
    
    def DTU_blade_calculation(self):
        results_list_for_DTU_airfoil = [] # η λίστα που θα αποθηκεύει τα αποτελέσματα για το κάθε τμήμα του πτερυγίου
        total_power = 0 # αρχικά η συνολική ισχύς είναι 0
        total_torque = 0 # αρχικά η συνολική ροπή είναι 0
        total_thrust = 0 # αρχικά η συνολική ώση είναι 0
        for i in range(self.no_sections):
            r = self.r_is[i]
            chord = self.chords[i]
            pitch_angle = self.pitch[i]
            tc_ratio = self.tc_ratios[i]
            twist = 0 
            try:
                results_for_DTU_airfoil = self.segment_calculation(r=r, chord=chord, pitch_angle_deg=pitch_angle, twist_deg=twist, tc_ratio=tc_ratio)
                # Η εφαπτομενική συνιστώσα είναι αυτή που παράγει τη ροπή
                dr = (self.r_is[i+1] - self.r_is[i]) if i < self.no_sections - 1 else (self.R - r)
                flow_angle, a, a_p = results_for_DTU_airfoil["flow_angle (degrees)"], results_for_DTU_airfoil["a"], results_for_DTU_airfoil["a_p"]
                Cn = results_for_DTU_airfoil["Cn"]
                Ct = results_for_DTU_airfoil["Ct"]
                dM = (
                    0.5 * self.air_density * self.B *
                    ((self.wind_speed_V0 * (1 - a) * self.rotation_speed * r * (1 + a_p)) /
                     (np.sin(flow_angle) * np.cos(flow_angle))) *
                    chord * Ct * r * dr
                )
                dT = ( 
                    0.5 * self.air_density * self.B * 
                ((self.wind_speed_V0**2 * (1 - a)**2) / np.sin(flow_angle)**2) * chord * Cn * dr 
                )
                # dM = r * self.B * results["pt"] * dr
                # dM = 4 * np.pi * self.air_density * self.wind_speed_V0 * self.rotation_speed * (1 - a) * a_p * r**3 * dr
                power = self.rotation_speed * dM
                total_power += power # Συνολική ισχύς όλου του ρότορα 
                total_torque += dM # Συνολική ροπή του ρότορα
                total_thrust += dT # Συνολική ώση του ρότορα
                # dP = 0.5 * self.air_density * 2 * np.pi * self.wind_speed_V0**3 * Ct * r * dr
                # total_power += dP
                results_for_DTU_airfoil["t/c ratio"] = tc_ratio
                results_for_DTU_airfoil["dT (Ν)"] = (dT/3) # Διαιρώ δια 3 καθώς η συγκεκριμένη τιμή αφορά και τα τρία πτερύγια
                results_for_DTU_airfoil["dM (Nm)"] = (dM/3) # Διαιρώ δια 3 καθώς η συγκεκριμένη τιμή αφορά και τα τρία πτερύγια
                results_for_DTU_airfoil["Power (Watt)"] = ((power)/3) # Διαιρώ δια 3 καθώς η συγκεκριμένη τιμή αφορά και τα τρία πτερύγια
                results_list_for_DTU_airfoil.append(results_for_DTU_airfoil)
            except Exception as e:
                print(f"Section {i} at radius {r}: {e}")
        return results_list_for_DTU_airfoil, total_power, total_torque, total_thrust
    
    def calculation_of_coefficient_of_power_cp_for_DTU(self, total_power):
        swept_area = np.pi * self.R**2 # επιφάνεια σάρωσης
        wind_power = 0.5 * self.air_density * swept_area * self.wind_speed_V0**3
        DTU_cp = total_power / wind_power
        return DTU_cp
    
    def calculation_of_coefficient_of_thrust_CT_for_DTU(self, total_thrust):
        swept_area = np.pi * self.R**2 # επιφάνεια σάρωσης
        wind_force = 0.5 * self.air_density * swept_area * self.wind_speed_V0**2
        DTU_CT = total_thrust / wind_force
        return DTU_CT
    
#%%
if __name__ == "__main__":
    blade_geom_file_2 = "blade_geom_file_2.json"
    hansen_DTU = Hansen_Algorithm_for_DTU_geometry(
        wind_speed_V0=10,
        rotation_speed=0,
        blade_geom_file_2=blade_geom_file_2,
        B=3,
        air_density=1.225
    )
    
    # ΔΙΑΓΡΑΜΜΑ Power Coefficient Cp - Tip Speed Ratio λ for DTU geometry
    rotation_speed_values_for_DTU = np.linspace(0.01, 0.56, 50)
    lambda_values_for_DTU = []
    cp_values = []
    
    for rotational_speed in rotation_speed_values_for_DTU:
        hansen_DTU.rotation_speed = rotational_speed
        λ = (rotational_speed * hansen_DTU.R) / hansen_DTU.wind_speed_V0
        lambda_values_for_DTU.append(λ)
        
        results_for_Naca_geometry, total_power, total_torque, total_thrust = hansen_DTU.DTU_blade_calculation()
        cp = hansen_DTU.calculation_of_coefficient_of_power_cp_for_DTU(total_power)
        cp_values.append(cp)
        
    plt.figure(figsize=(10, 8))
    plt.plot(lambda_values_for_DTU, cp_values, 'o-', label="$C_p$ vs $λ$")
    plt.xlabel("$λ$ (Tip speed ratio)", fontsize=15)
    plt.ylabel("$C_p$ (Power coefficient)", fontsize=15)
    plt.title("Coefficient of Power $C_p$ as a function of Tip Speed Ratio $λ$ for DTU geometry")
    plt.legend()
    plt.grid()
    plt.show()
    
    # ΔΙΑΓΡΑΜΜΑ ταχύτητας n(rpm) - ισχύος P(kWatt) for DTU geometry
    rpm_values_for_DTU = np.linspace(0, 6, 50) # Τιμές RPM
    wind_speed_values_for_DTU = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15] # Ταχύτητες ανέμου

    plt.figure(figsize=(10, 6))

    for V0 in wind_speed_values_for_DTU:
        power_values_for_DTU = []
        for n in rpm_values_for_DTU:
            ω = (2 * np.pi * n) / 60 # Γωνιακή ταχύτητα σε rad/s
            λ = (ω * hansen_DTU.R) / V0 # Tip Speed Ratio
            hansen_DTU.wind_speed_V0 = V0
            hansen_DTU.rotation_speed = ω
            results_for_DTU_geometry, total_power, total_torque, total_thrust = hansen_DTU.DTU_blade_calculation()
            power_values_for_DTU.append(total_power*1e-3) # Ισχύς σε kW
            
        plt.plot(rpm_values_for_DTU, power_values_for_DTU,'o', label=f"{V0} m/s", )

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
# %%
