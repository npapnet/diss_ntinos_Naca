
#%%
import numpy as np
from fakelos.naca4415 import Naca4415
import pandas as pd
import json

class HansenAlgorithm:
    tolerance = 1e-6 # Σταθερά για τον έλεγχο σύγκλισης των τιμών των συντελεστών a και a'
    max_iter = 10000 # Μέγιστος αριθμός επαναλήψεων
    
    def __init__(self, wind_speed_V0, R, angular_speed, B=3, air_density=1.225, csv_data_file='csv_data_file.csv'):
        self.wind_speed_V0 = wind_speed_V0 # ταχύτητα του ανέμου (σε m/sec)
        self.R = R
        self.angular_speed = angular_speed # γωνιακή ταχύτητα (σε rad/sec)
        self.B = B # αριθμός πτερυγίων 
        self.air_density = air_density # πυκνότητα του αέρα (σε kg/m^3)
        self.naca4415 = Naca4415(csv_data_file)
         
    def calculation_of_flow_angle(self, a, a_accent, r): # μέθοδος για τον υπολογισμό της γωνίας ροής φ (ΒΗΜΑ 2 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        # a_accent = a'
        v_axial = self.wind_speed_V0 * (1 - a) # αξονική ταχύτητα
        v_tangential = self.angular_speed * r * (1 + a_accent) # εφαπτομενική ταχύτητα 
        return np.arctan(v_axial / v_tangential)
    
    def calculation_of_local_angle_of_attack(self, phi, theta_p, twist): # μέθοδος για τον υπολογισμό της τοπικής γωνίας προσβολής (ΒΗΜΑ 3 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        # phi = η γωνία ροής φ
        # theta = θ
        # theta_p = η γωνία κλίσης θp
        # twist = η συστροφή του πτερυγίου β
        return phi - (theta_p + twist)
    
    def cl_and_cd(self, alpha): # πίνακας για τους συντελεστές άνωσης και οπισθέλκουσας Cl and Cd βάσει του alpha (γωνία προσβολής) (ΒΗΜΑ 4 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        # Θα πρέπει να εισάγω κάποιον πίνακα για τους συντελεστές αυτούς
        Cl = self.naca4415.cl(alpha)
        Cd = self.naca4415.cd(alpha)
        return Cl, Cd
    
    def calculation_of_Cn_and_Ct(self, Cl, Cd, phi): # μέθοδος για τον υπολογισμό των συντελεστών Cn(normal) και Ct(tangential) (ΒΗΜΑ 5 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
        Ct = Cl * np.sin(phi) - Cd * np.cos(phi)
        return Cn, Ct
    
    def updated_induction_factors(self, Cn, Ct, r, chord, phi): # μέθοδος για τον υπολογισμό των νέων τιμών των συντελεστών a και a΄ (ΒΗΜΑ 6 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        sigma = (self.B * chord) / (2 * np.pi * r)
        a_new = 1 / ((4 * np.sin(phi)**2) / (sigma * Cn) + 1)
        a_accent_new = 1 / ((4 * np.sin(phi)*np.cos(phi)) / (sigma * Ct) - 1)
        return a_new, a_accent_new

    def calculation_of_local_forces(self, r, a, a_accent, chord, Cl, Cd): # υπολογισμός των τοπικών δυνάμεων στα τμήματα του πτερυγίου (ΒΗΜΑ 8 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        v_axial = self.wind_speed_V0 * (1 - a)
        v_tangential = self.angular_speed * r * (1 + a_accent)
        Vrel = np.sqrt((v_axial)**2 + (v_tangential)**2)
        L = 0.5 * self.air_density * (Vrel)**2 * Cl * chord
        D = 0.5 * self.air_density * (Vrel)**2 * Cd * chord
        return L, D
     
    def run_the_algorithm(self, r, theta_p, twist, chord, f = 0.3): # εκτέλεση του αλγόριθμου
        a, a_accent = 0 , 0 # αρχικοποίηση των συντελεστών a και a' σε 0
        exit_flag = False
        counter = 0 
        
        while not exit_flag:
            phi = self.calculation_of_flow_angle(a, a_accent, r)
            alpha = self.calculation_of_local_angle_of_attack(phi, theta_p, twist)
            Cl, Cd = self.cl_and_cd(alpha)
            Cn, Ct = self.calculation_of_Cn_and_Ct(Cl, Cd, phi)
            a_new, a_accent_new = self.updated_induction_factors(Cn, Ct, r, chord, phi)
            
            if abs(a - a_new) < self.tolerance and abs(a_accent - a_accent_new) < self.tolerance: # έχω σύγκλιση
                exit_flag = True
            else: # συνεχίζω τον αλγόριθμο
                a = a * (1 - f) + f * a_new
                a_accent = a_accent * (1 - f) + f * a_accent_new  
            counter += 1
            if counter > self.max_iter:
                raise Exception("Δεν έχω σύγκλιση")
            L, D = self.calculation_of_local_forces(r, a, a_accent, chord, Cl, Cd)
        
        return {
            "r_i": r,
            "theta_p": theta_p,
            "twist": twist,
            "chord": chord,
            "a": a,
            "a_accent": a_accent ,
            "phi": phi,
            "alpha": alpha,
            "Cl": Cl,
            "Cd": Cd,
            "Cn": Cn,
            "Ct": Ct,
            "a_new": a_new,
            "a_accent_new": a_accent_new,
            "Lift": L,
            "Drag": D,
            "counter": counter
        }
#%%
dic_list = []
hansen = HansenAlgorithm(wind_speed_V0=8, R=50, angular_speed=2, B=3, air_density=1.225)
for k in np.arange(0,4,step=0.1):
    try:
        results = hansen.run_the_algorithm(r=25, theta_p=2, twist=k, chord=14, f=0.3)
        dic_list.append(results)
    except:
        print(k)
        break
    
df_res = pd.DataFrame(dic_list)
#%%
hansen = HansenAlgorithm(wind_speed_V0=8, R=50, angular_speed=2, B=3, air_density=1.225)
results = hansen.run_the_algorithm(r=25, theta_p=2, twist=2.4, chord=14, f=0.3)
#%%
# hansen = HansenAlgorithm(wind_speed_V0=8, R=50, angular_speed=2, B=3, air_density=1.225)
# results = hansen.run_the_algorithm(r=25, theta_p=2, twist=0.05, chord=4)
for key, value in results.items():
    print(f"{key:15s}: {value:.3g}")
#%%
if __name__ == "__main__" :
    hansen = HansenAlgorithm(wind_speed_V0=8, R=50, angular_speed=2, B=3, air_density=1.225)
  
print(df_res)