#%%
import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from algorithmos import Hansen_Algorithm

#%%
class Hansen_Algorithm_for_Naca_geometry(Hansen_Algorithm): # Νέα κλάση για τη γεωμετρία του πτερυγίου (Naca airfoil)
    def __init__(self, wind_speed_V0, rotation_speed, blade_geom_file, B=3, air_density=1.225, csv_data_file_Naca='csv_data_file_Naca.csv'):
        with open(blade_geom_file, 'r') as f:
            blade_geom_file = json.load(f)
        self.r_is = blade_geom_file["r_is"] # πίνακας με την ακτινική θέση του κάθε τμήματος του πτερυγίου
        self.chords = blade_geom_file["chords"] # πίνακας με το μήκος χορδής του κάθε τμήματος του πτερυγίου
        self.pitch = blade_geom_file["pitch"] # πίνακας με τη γωνία βήματος του κάθε τμήματος του πτερυγίου
        self.lambda0 = blade_geom_file["lambda0"] # ο λόγος ταχυτήτων ακροπτερυγίου
        self.no_sections = blade_geom_file["no_sections"] # αριθμός τμημάτων του πτερυγίου
        R = blade_geom_file["R"] # ακτίνα του ρότορα
        super().__init__(wind_speed_V0, R, rotation_speed, B, air_density, airfoil_type="NACA", csv_data_file=csv_data_file_Naca)
    
    def Naca_blade_calculation(self):
        results_list_for_Naca_airfoil = [] # η λίστα που θα αποθηκεύει τα αποτελέσματα για το κάθε τμήμα του πτερυγίου
        total_power = 0 # αρχικά η συνολική ισχύς είναι 0
        total_torque = 0 # αρχικά η συνολική ροπή είναι 0
        total_thrust = 0 # αρχικά η συνολική ώση είναι 0
        for i in range(self.no_sections):
            r = self.r_is[i]
            chord = self.chords[i]
            pitch_angle = self.pitch[i]
            twist = 0 
            try:
                results_for_Naca_airfoil = self.segment_calculation(r=r, chord=chord, pitch_angle_deg=pitch_angle, twist_deg=twist, tc_ratio=None)
                # Η εφαπτομενική συνιστώσα είναι αυτή που παράγει τη ροπή
                dr = (self.r_is[i+1] - self.r_is[i]) if i < self.no_sections - 1 else (self.R - r)
                flow_angle, a, a_p = results_for_Naca_airfoil["flow_angle (rads)"], results_for_Naca_airfoil["a"], results_for_Naca_airfoil["a_p"]
                Cn = results_for_Naca_airfoil["Cn"]
                Ct = results_for_Naca_airfoil["Ct"]
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
                power = self.rotation_speed * dM
                total_power += power # Συνολική ισχύς όλου του ρότορα 
                total_torque += dM # Συνολική ροπή του ρότορα
                total_thrust += dT # Συνολική ώση του ρότορα
                results_for_Naca_airfoil["dT (N)"] = dT/3 # Διαιρώ δια 3 καθώς η συγκεκριμένη τιμή αφορά και τα τρία πτερύγια
                results_for_Naca_airfoil["dM (Nm)"] = dM/3 # Διαιρώ δια 3 καθώς η συγκεκριμένη τιμή αφορά και τα τρία πτερύγια
                results_for_Naca_airfoil["Power (Watt)"] = power/3 # Διαιρώ δια 3 καθώς η συγκεκριμένη τιμή αφορά και τα τρία πτερύγια
                results_list_for_Naca_airfoil.append(results_for_Naca_airfoil)
            except Exception as e:
                print(f"Section {i} at radius {r}: {e}")
        return results_list_for_Naca_airfoil, total_power, total_torque, total_thrust
    
    def calculation_of_coefficient_of_power_cp_for_Naca(self, total_power):
        swept_area = np.pi * self.R**2 # επιφάνεια σάρωσης
        wind_power = 0.5 * self.air_density * swept_area * self.wind_speed_V0**3
        return total_power / wind_power
    
    def calculation_of_coefficient_of_thrust_CT_for_Naca(self, total_thrust):
        swept_area = np.pi * self.R**2 # επιφάνεια σάρωσης
        wind_force = 0.5 * self.air_density * swept_area * self.wind_speed_V0**2
        return total_thrust / wind_force
        
#%%
# Εκτέλεση του αλγορίθμου
if __name__ == "__main__":
    blade_geom_file = "blade_geom_file.json"
    hansen_Naca = Hansen_Algorithm_for_Naca_geometry(
        wind_speed_V0=10,
        rotation_speed=0,
        blade_geom_file=blade_geom_file,
        B=3,
        air_density=1.225
    )
    
    # ΔΙΑΓΡΑΜΜΑ Power Coefficient Cp - Tip Speed Ratio λ for Naca geometry
    rotation_speed_values = np.linspace(1, 100, 50)
    lambda_values = []
    cp_values = []
    
    for rotational_speed in rotation_speed_values:
        hansen_Naca.rotation_speed = rotational_speed
        λ = (rotational_speed * hansen_Naca.R) / hansen_Naca.wind_speed_V0
        lambda_values.append(λ)
        
        results_for_Naca_geometry, total_power, total_torque, total_thrust = hansen_Naca.Naca_blade_calculation()
        cp = hansen_Naca.calculation_of_coefficient_of_power_cp_for_Naca(total_power)
        cp_values.append(cp)
        
    plt.figure(figsize=(10, 8))
    plt.plot(lambda_values, cp_values, 'o-', label="$C_p$ vs $λ$")
    plt.xlabel("$λ$ (Tip speed ratio)", fontsize=15)
    plt.ylabel("$C_p$ (Power coefficient)", fontsize=15)
    plt.title("Coefficient of Power $C_p$ as a function of Tip Speed Ratio $λ$ for Naca geometry")
    plt.legend()
    plt.grid()
    plt.show()
    
    # ΔΙΑΓΡΑΜΜΑ ταχύτητας n(rpm) - ισχύος P(kWatt) for Naca geometry
    rpm_values = np.linspace(0, 1000, 50) # Τιμές RPM
    wind_speed_values = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15] # Ταχύτητες ανέμου

    plt.figure(figsize=(10, 6))

    for V0 in wind_speed_values:
        power_values = []
        for n in rpm_values:
            ω = (2 * np.pi * n) / 60 # Γωνιακή ταχύτητα σε rad/s
            λ = (ω * hansen_Naca.R) / V0 # Tip Speed Ratio
            hansen_Naca.wind_speed_V0 = V0
            hansen_Naca.rotation_speed = ω
            results_for_Naca_geometry, total_power, total_torque, total_thrust = hansen_Naca.Naca_blade_calculation()
            power_values.append(total_power*1e-3) # Ισχύς σε kW

        plt.plot(rpm_values, power_values,'o', label=f"{V0} m/s", )

    plt.xlabel("Rotational Speed (RPM)", fontsize=12)
    plt.ylabel("Power (kW)", fontsize=12)
    plt.title("Power vs Rotational Speed for Different Wind Speeds for Naca geometry", fontsize=14)
    plt.legend(title="Wind Speed [m/s]")
    plt.grid()
    plt.show()
    
    df_Naca_results = pd.DataFrame(results_for_Naca_geometry)
    print(df_Naca_results)
    
    print(f"Η συνολική ισχύς του πτερυγίου είναι {(total_power/3):.2f} Watt")
    print(f"Η συνολική ροπή του πτερυγίου είναι {(total_torque/3):.2f} Nm")
    print(f"Η συνολική ώση του πτερυγίου είναι {(total_thrust/3):.2f} N")
    print(f"Η συνολική ισχύς της ανεμογεννήτριας είναι {total_power:.2f} Watt")
    print(f"H συνολική ροπή της ανεμογεννήτριας είναι {total_torque:.2f} Nm")
    print(f"H συνολική ώση της ανεμογεννήτριας είναι {total_thrust:.2f} N")