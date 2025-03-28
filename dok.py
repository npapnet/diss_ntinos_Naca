#%%
import numpy as np
from naca4415_2 import Naca4415
import pandas as pd
import json
import matplotlib.pyplot as plt

class HansenAlgorithm:
    tolerance = 1e-4 # Σταθερά για τον έλεγχο σύγκλισης των τιμών των συντελεστών a και a'
    max_iter = 1000 # Μέγιστος αριθμός επαναλήψεων
    
    def __init__(self, wind_speed_V0, R, rotation_speed, B=3, air_density=1.225, csv_data_file='csv_data_file.csv'):
        self.wind_speed_V0 = wind_speed_V0 # ταχύτητα του ανέμου (σε m/sec)
        self.R = R # ακτίνα του ρότορα
        self.rotation_speed = rotation_speed # ταχύτητα περιστροφής του ρότορα (σε rad/sec)
        self.B = B # αριθμός πτερυγίων 
        self.air_density = air_density # πυκνότητα του αέρα (σε kg/m^3)
        self.naca4415 = Naca4415(csv_data_file)
        
    def calculation_of_flow_angle(self, a, a_p, r): # μέθοδος για τον υπολογισμό της γωνίας ροής φ (ΒΗΜΑ 2 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        # a_p = a'
        v_axial = self.wind_speed_V0 * (1 - a) # αξονική ταχύτητα
        v_tangential = self.rotation_speed * r * (1 + a_p) # εφαπτομενική ταχύτητα 
        return np.arctan(v_axial / v_tangential)
    
    def calculation_of_local_angle_of_attack(self, phi, theta_p, twist): # μέθοδος για τον υπολογισμό της τοπικής γωνίας προσβολής (ΒΗΜΑ 3 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        # phi = η γωνία ροής φ
        # theta = θ
        # theta_p = η γωνία κλίσης θp
        # twist = η συστροφή του πτερυγίου β
        # theta = theta_p + twist
        return phi - (theta_p + twist)
    
    def calculation_of_Cl_and_Cd(self, alpha): # πίνακας για τους συντελεστές άνωσης και οπισθέλκουσας Cl and Cd βάσει του alpha (γωνία προσβολής) (ΒΗΜΑ 4 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        Cl = self.naca4415.cl(alpha)
        Cd = self.naca4415.cd(alpha)
        return Cl, Cd
    
    def calculation_of_Cn_and_Ct(self, Cl, Cd, phi): # μέθοδος για τον υπολογισμό των συντελεστών Cn(normal) και Ct(tangential) (ΒΗΜΑ 5 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        Cn = (Cl * np.cos(phi)) + (Cd * np.sin(phi))
        Ct = (Cl * np.sin(phi)) - (Cd * np.cos(phi))
        return Cn, Ct
    
    def updated_induction_factors(self, Cn, Ct, r, chord, phi): # μέθοδος για τον υπολογισμό των νέων τιμών των συντελεστών a και a΄ (ΒΗΜΑ 6 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        sigma = (self.B * chord) / (2 * np.pi * r) # solidity factor (στερεότητα)
        a_new = 1 / ((4 * np.sin(phi)**2) / (sigma * Cn) + 1)
        a_p_new = 1 / ((4 * np.sin(phi)*np.cos(phi)) / (sigma * Ct) - 1)
        return a_new, a_p_new

    def calculation_of_local_forces(self, r, a, a_p, chord, phi, Cl, Cd): # υπολογισμός των τοπικών φορτίων στα τμήματα του πτερυγίου (ΒΗΜΑ 8 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        v_axial = self.wind_speed_V0 * (1 - a)
        v_tangential = self.rotation_speed * r * (1 + a_p)
        Vrel = np.sqrt((v_axial)**2 + (v_tangential)**2)
        L = 0.5 * self.air_density * (Vrel)**2 * Cl * chord
        D = 0.5 * self.air_density * (Vrel)**2 * Cd * chord
        pn = L * np.cos(phi) + D * np.sin(phi)
        pt = L * np.sin(phi) - D * np.cos(phi)
        return L, D, pn, pt
     
    def segment_calcultion(self, r, theta_p, twist, chord, f=0.3): # εκτέλεση του αλγόριθμου
        a, a_p = 0, 0  # αρχικοποίηση των συντελεστών επαγωγής a και a' σε 0
        exit_flag = False
        counter = 0 
        
        while not exit_flag:
            phi = self.calculation_of_flow_angle(a, a_p, r)
            alpha = self.calculation_of_local_angle_of_attack(phi, theta_p, twist)
            Cl, Cd = self.calculation_of_Cl_and_Cd(alpha)
            Cn, Ct = self.calculation_of_Cn_and_Ct(Cl, Cd, phi)
            a_new, a_p_new = self.updated_induction_factors(Cn, Ct, r, chord, phi)
            
            if abs(a - a_new) < self.tolerance and abs(a_p - a_p_new) < self.tolerance: # έχω σύγκλιση
                exit_flag = True
            else: # συνεχίζω τον αλγόριθμο
                a = a * (1 - f) + f * a_new
                a_p = a_p * (1 - f) + f * a_p_new  
            counter += 1
            if counter > self.max_iter:
                raise Exception("Δεν έχω σύγκλιση")
            L, D, pn, pt = self.calculation_of_local_forces(r, a, a_p, chord, phi, Cl, Cd)
            
        
        return {
            "r_i": r,
            "chord": chord,
            "theta_p": theta_p,
            "twist": twist,
            "a": a,
            "a_p": a_p,
            "phi": phi,
            "alpha": alpha,
            "Cl": Cl,
            "Cd": Cd,
            "Cn": Cn,
            "Ct": Ct,
            "a_new": a_new,
            "a_p_new": a_p_new,
            "Lift": L,
            "Drag": D,
            "pn": pn,
            "pt": pt,
            "counter": counter
        }
#%%
class HansenAlgorithmWithBladeGeom(HansenAlgorithm): # Νέα κλάση για τη γεωμετρία του πτερυγίου
    def __init__(self, wind_speed_V0, rotation_speed, blade_geom_file, B=3, air_density=1.225):
        with open(blade_geom_file, 'r') as f:
            blade_geom_file = json.load(f)
        self.r_is = blade_geom_file["r_is"]
        self.chords = blade_geom_file["chords"]
        self.pitch = blade_geom_file["pitch"]
        self.lambda0 = blade_geom_file["lambda0"]
        self.no_sections = blade_geom_file["no_sections"]
        R = blade_geom_file["R"]
        super().__init__(wind_speed_V0, R, rotation_speed, B, air_density)
    
    def blade_calculation(self):
        results_list = []
        total_power = 0
        for i in range(self.no_sections):
            r = self.r_is[i]
            chord = self.chords[i]
            theta_p = self.pitch[i]
            twist = 0 
            try:
                results = self.segment_calcultion(r=r, theta_p=theta_p, twist=twist, chord=chord)
                # Ft = results["Lift"] * np.sin(results["phi"]) - results["Drag"] * np.cos(results["phi"]) # Η εφαπτομενική συνιστώσα είναι αυτή που παράγει την ροπή
                # torque = Ft * r # η ροπή σε κάθε τμήμα του φτερού
                # power = torque * self.rotation_speed
                # results["Power"] = power
                # results["Ft"] = Ft
                # total_power += power
                # results_list.append(results)
                dr = (self.r_is[i+1] - self.r_is[i]) if i < self.no_sections - 1 else (self.R - r)
                phi, a, a_p = results["phi"], results["a"], results["a_p"]
                Ct = results["Ct"]
                dM = (
                    0.5 * self.air_density * self.B *
                    ((self.wind_speed_V0 * (1 - a) * self.rotation_speed * r * (1 + a_p)) /
                     (np.sin(phi) * np.cos(phi))) *
                    chord * Ct * r * dr
                )
                # dM = r * self.B * results["pt"] * dr
                power = self.rotation_speed * dM
                total_power += power
                # results["r"] = r
                results["chord"] = chord
                results["power"] = power
                results_list.append(results)
            except Exception as e:
                print(f"Section {i} at radius {r}: {e}")
        return results_list, total_power
    
    def calculate_coefficient_of_power_cp(self, total_power):
        swept_area = np.pi * self.R**2 # επιφάνεια σάρωσης
        wind_power = 0.5 * self.air_density * swept_area * self.wind_speed_V0**3
        cp = total_power / wind_power
        return cp
#%%
# Εκτέλεση του αλγορίθμου
if __name__ == "__main__":
    blade_geom_file = "blade_geom_file.json"
    hansen = HansenAlgorithmWithBladeGeom(
        wind_speed_V0=10,
        rotation_speed=0,
        blade_geom_file=blade_geom_file,
        B=3,
        air_density=1.225
    )
    
    # ΔΙΑΓΡΑΜΜΑ Power Coefficient Cp - λ
    rotation_speed_values = np.linspace(1, 100, 50)
    lambda_values = []
    cp_values = []
    
    for rotational_speed in rotation_speed_values:
        hansen.rotation_speed = rotational_speed
        λ = (rotational_speed * hansen.R) / hansen.wind_speed_V0
        lambda_values.append(λ)
        
        results, total_power = hansen.blade_calculation()
        cp = hansen.calculate_coefficient_of_power_cp(total_power)
        cp_values.append(cp)
        total_power_of_wind_turbine = hansen.B * total_power
        
    plt.figure(figsize=(10, 8))
    plt.plot(lambda_values, cp_values, 'o-', label="$C_p$ vs $λ$")
    plt.xlabel("$λ$ (Tip speed ratio)", fontsize=15)
    plt.ylabel("$C_p$ (Power coefficient)", fontsize=15)
    plt.title("Coefficient of Power $C_p$ as a function of Tip Speed Ratio $λ$")
    plt.legend()
    plt.grid()
    plt.show()
    
    # ΔΙΑΓΡΑΜΜΑ ταχύτητας n(rpm) - ισχύος P(kWatt) 
    rpm_values = np.linspace(0, 1000, 50) # Τιμές RPM
    wind_speed_values = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15] # Ταχύτητες ανέμου

    plt.figure(figsize=(10, 6))

    for V0 in wind_speed_values:
        power_values = []
        for n in rpm_values:
            ω = (2 * np.pi * n) / 60 # Γωνιακή ταχύτητα σε rad/s
            λ = (ω * hansen.R) / V0 # Tip Speed Ratio
            hansen.wind_speed_V0 = V0
            hansen.rotation_speed = ω
            results, total_power = hansen.blade_calculation()
            power_values.append(total_power*1e-3)  # Ισχύς σε kW

        plt.plot(rpm_values, power_values,'o', label=f"{V0} m/s", )

    plt.xlabel("Rotational Speed (RPM)", fontsize=12)
    plt.ylabel("Power (kW)", fontsize=12)
    plt.title("Power vs Rotational Speed for Different Wind Speeds", fontsize=14)
    plt.legend(title="Wind Speed [m/s]")
    plt.grid(True)
    plt.show()
   
    df_results = pd.DataFrame(results)
    print(df_results)
    print(f"Συνολική Ισχύς του φτερού: {total_power:.2f} Watt")
    print(f"Συντελεστής Απόδοσης (Cp): {cp:.4f}")
    print(f"Συνολική Ισχύς της Ανεμογεννήτριας: {total_power_of_wind_turbine:.2f} Watt")
# %%
