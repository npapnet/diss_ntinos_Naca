#%%
import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from algorithmos import Hansen_Algorithm

#%%
def new_blade_geometry(r_is, chords, pitch, tc_ratios,
                            r_first, r_last, num_sections=10):
    """
    Δημιουργεί νέο διάνυσμα ακτινικών θέσεων r_new με μήκος num_sections
    (δηλ. num_sections σημεία) γραμμικά κατανεμημένα από r_first μέχρι r_last.
    Στη συνέχεια υπολογίζει με γραμμική παρεμβολή τις τιμές chord, pitch, tc_ratio
    στα r_new.
    """
    # ορίζουμε τα νέα ακτινικά σημεία
    r_new = np.linspace(r_first, r_last, num_sections)

    # φτιάχνουμε τις συναρτήσεις παρεμβολής από τα αρχικά δεδομένα
    chord_interp = interp1d(r_is, chords, kind='linear')
    pitch_interp = interp1d(r_is, pitch, kind='linear')
    tc_interp = interp1d(r_is, tc_ratios, kind='linear')
    
    # υπολογίζουμε τις νέες τιμές
    chords_new = chord_interp(r_new)
    pitch_new = pitch_interp(r_new)
    tc_new = tc_interp(r_new)

    return r_new, chords_new, pitch_new, tc_new

#%%
class Hansen_Algorithm_for_DTU_geometry(Hansen_Algorithm): # Νέα κλάση για τη γεωμετρία του πτερυγίου (DTU airfoil)
    def __init__(self, blade_geom_DTU, B=3, air_density=1.225, csv_data_file_DTU='csv_data_file_DTU.csv'):
        with open(blade_geom_DTU, 'r') as f:
            blade_geom_DTU = json.load(f)
        r_is_original = blade_geom_DTU["r_is"]
        chords_original = blade_geom_DTU["chords"]
        pitch_original = blade_geom_DTU["pitch"]
        tc_ratios_original = blade_geom_DTU["tc_ratios"]
        no_sections = blade_geom_DTU["no_sections"]

        r_first = r_is_original[0] 
        r_last = np.max(r_is_original)             

        # Δημιουργούμε 10 σημεία με γραμμική παρεμβολή
        r_new, chords_new, pitch_new, tc_new = new_blade_geometry(
            r_is_original,
            chords_original,
            pitch_original,
            tc_ratios_original,
            r_first,
            r_last,
            num_sections = 10
        )

        self.r_is = r_new
        self.chords = chords_new
        self.pitch = pitch_new
        self.tc_ratios = tc_new
        self.no_sections = len(r_new)
        R_new = np.max(r_new)
        
        super().__init__(R_new,
                         B, air_density,
                         airfoil_type="DTU",
                         csv_data_file=csv_data_file_DTU)
        
    def DTU_blade_calculation(self, wind_speed_V0,
        rotation_speed):
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
                results_for_DTU_airfoil = self.segment_calculation(
                    wind_speed_V0=wind_speed_V0, 
                    omega_rad_sec=rotation_speed,
                    r=r, chord=chord, pitch_angle_deg=pitch_angle, twist_deg=twist, tc_ratio=tc_ratio)
                # Η εφαπτομενική συνιστώσα είναι αυτή που παράγει τη ροπή
                dr = (self.r_is[i+1] - self.r_is[i]) if i < self.no_sections - 1 else (self.R - r)
                flow_angle_rad, a, a_p = results_for_DTU_airfoil["flow_angle (rads)"], results_for_DTU_airfoil["a"], results_for_DTU_airfoil["a_p"]
                Cn = results_for_DTU_airfoil["Cn"]
                Ct = results_for_DTU_airfoil["Ct"]

                dM = (
                    0.5 * self.air_density * self.B *
                    ((wind_speed_V0 * (1 - a) * rotation_speed * r * (1 + a_p)) /
                     (np.sin(flow_angle_rad) * np.cos(flow_angle_rad))) *
                    chord * Ct * r * dr
                )
                dT = ( 
                    0.5 * self.air_density * self.B * 
                ((wind_speed_V0**2 * (1 - a)**2) / (np.sin(flow_angle_rad)**2)) * chord * Cn * dr 
                )
                # dM = r * self.B * results["pt"] * dr
                # dM = 4 * np.pi * self.air_density * self.wind_speed_V0 * self.rotation_speed * (1 - a) * a_p * r**3 * dr
                power = rotation_speed * dM
                total_power += power # Συνολική ισχύς όλου του ρότορα 
                total_torque += dM # Συνολική ροπή του ρότορα
                total_thrust += dT # Συνολική ώση του ρότορα
                # dP = 0.5 * self.air_density * 2 * np.pi * self.wind_speed_V0**3 * Ct * r * dr
                results_for_DTU_airfoil["t/c ratio"] = tc_ratio
                results_for_DTU_airfoil["dT (Ν)"] = (dT/3) # Διαιρώ δια 3 καθώς η συγκεκριμένη τιμή αφορά και τα τρία πτερύγια
                results_for_DTU_airfoil["dM (Nm)"] = (dM/3) # Διαιρώ δια 3 καθώς η συγκεκριμένη τιμή αφορά και τα τρία πτερύγια
                results_for_DTU_airfoil["Power (Watt)"] = ((power)/3) # Διαιρώ δια 3 καθώς η συγκεκριμένη τιμή αφορά και τα τρία πτερύγια
                results_list_for_DTU_airfoil.append(results_for_DTU_airfoil)
            except Exception as e:
                print(f"Section {i} at radius {r}: {e}")
        return results_list_for_DTU_airfoil, total_power, total_torque, total_thrust
    
    def calculation_of_coefficient_of_power_cp_for_DTU(self, total_power, wind_speed_V0):
        swept_area = np.pi * self.R**2 # επιφάνεια σάρωσης
        wind_power = 0.5 * self.air_density * swept_area * wind_speed_V0**3
        return total_power / wind_power
    
    def calculation_of_coefficient_of_thrust_CT_for_DTU(self, total_thrust,wind_speed_V0):
        swept_area = np.pi * self.R**2 # επιφάνεια σάρωσης
        wind_force = 0.5 * self.air_density * swept_area * wind_speed_V0**2
        return total_thrust / wind_force

#%%
if __name__ == "__main__":
    blade_geom_file_2 = "blade_geom_DTU.json"
    hansen_DTU = Hansen_Algorithm_for_DTU_geometry(

        blade_geom_DTU=blade_geom_file_2,
        B=3,
        air_density=1.225
    )
    
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
