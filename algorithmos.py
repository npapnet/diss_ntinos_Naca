#%% 
import json
import pandas as pd
import numpy as np
from Naca_table import Naca_calc
from Dtu_table import DTU_calc
import matplotlib.pyplot as plt

#%%
class Hansen_Algorithm:
    tolerance = 1e-4 # Σταθερά για τον έλεγχο σύγκλισης των τιμών των συντελεστών αξονικής και εφαπτομενικής επαγωγής a και a'
    max_iter = 1000 # Μέγιστος αριθμός επαναλήψεων
    
    def __init__(self, R, B=3, air_density=1.225, airfoil_type=None, csv_data_file=None):
        # self.wind_speed_V0 = wind_speed_V0 # ταχύτητα του ανέμου (σε m/sec)
        self.R = R # ακτίνα του ρότορα
        # self.rotation_speed = rotation_speed # ταχύτητα περιστροφής του ρότορα (σε rad/sec)
        self.B = B # αριθμός πτερυγίων 
        self.air_density = air_density # πυκνότητα του αέρα (σε kg/m^3)
        
        if airfoil_type == "NACA":
            self.airfoil_calc = Naca_calc(csv_data_file) # Χρήση NACA δεδομένων
        elif airfoil_type == "DTU":
            self.airfoil_calc = DTU_calc(csv_data_file) # Χρήση DTU δεδομένων
        else:
            self.airfoil_calc = None
        
    def calculation_of_flow_angle_rad(self, a, a_p, r, v0, w_rps): # μέθοδος για τον υπολογισμό της γωνίας ροής φ (ΒΗΜΑ 2 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        # a_p = a'
        
        v_axial = v0 * (1 - a) # συνιστώσα αξονικής ταχύτητας
        v_tangential = w_rps* r * (1 + a_p) # συνιστώσα εφαπτομενικής ταχύτητας
        
        if a_p == -1:
            raise ValueError("Διαίρεση με το 0")
        else:
            return np.arctan(v_axial / v_tangential)
    
    def calculation_of_local_angle_of_attack_rad(self, flow_angle_rad:float, pitch_angle_deg:float, twist_deg:float): # μέθοδος για τον υπολογισμό
        """ της τοπικής γωνίας προσβολής (ΒΗΜΑ 3 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)

        Args:
            flow_angle_rad (float): η γωνία της ροής φ σε rad
            pitch_angle_deg (float): η γωνία βήματος του πτερυγίου
            twist_deg (float):η γωνία συστροφής του πτερυγίου β

        Returns:
             (float): angle of attack in rad
        """
        return (flow_angle_rad - (np.radians(pitch_angle_deg) + twist_deg))
    
    def calculation_of_Cl_and_Cd(self, angle_of_attack_deg, tc_ratio=None): # πίνακας για τους συντελεστές άνωσης και οπισθέλκουσας Cl and Cd
        # βάσει του angle_of_attack (γωνία προσβολής), αλλά και του λόγου t/c (thickness/chord ratio) (ΒΗΜΑ 4 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        # Cl = self.DTU_calc.cl(angle_of_attack, tc_ratio)
        # Cd = self.DTU_calc.cd(angle_of_attack, tc_ratio)
        # return Cl, Cd
        if isinstance(self.airfoil_calc, Naca_calc):
            Cl = self.airfoil_calc.cl(angle_of_attack_deg) # Υπολογισμός Cl μόνο με βάση τη γωνία προσβολής
            Cd = self.airfoil_calc.cd(angle_of_attack_deg) # Υπολογισμός Cd μόνο με βάση τη γωνία προσβολής
        elif isinstance(self.airfoil_calc, DTU_calc):
            Cl = self.airfoil_calc.cl(angle_of_attack_deg, tc_ratio) # Υπολογισμός Cl με βάση τη γωνία προσβολής και t/c
            Cd = self.airfoil_calc.cd(angle_of_attack_deg, tc_ratio) # Υπολογισμός Cd με βάση τη γωνία προσβολής και t/c
        else:
            raise ValueError("Άγνωστος τύπος αεροτομής")
        return Cl, Cd
    
    def calculation_of_Cn_and_Ct(self, Cl, Cd, flow_angle_rad): # μέθοδος για τον υπολογισμό των συντελεστών 
        # Cn(normal) και Ct(tangential) (ΒΗΜΑ 5 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        Cn = (Cl * np.cos(flow_angle_rad)) + (Cd * np.sin(flow_angle_rad))
        Ct = (Cl * np.sin(flow_angle_rad)) - (Cd * np.cos(flow_angle_rad))
        return Cn, Ct
    
    def calculation_of_updated_induction_factors(self, Cn, Ct, r, chord, flow_angle_rad): # μέθοδος για τον υπολογισμό των νέων τιμών 
        # των συντελεστών επαγωγής a και a΄ (ΒΗΜΑ 6 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        solidity_factor = (self.B * chord) / (2 * np.pi * r) # solidity factor (συντελεστής στερεότητας)
        a_new = 1 / ((4 * np.sin(flow_angle_rad)**2) / (solidity_factor * Cn) + 1)
        a_p_new = 1 / ((4 * np.sin(flow_angle_rad)*np.cos(flow_angle_rad)) / (solidity_factor * Ct) - 1)
        return a_new, a_p_new

    def calculation_of_local_loads(self, r, a, a_p, chord, flow_angle, Cl, Cd): # υπολογισμός των τοπικών φορτίων
        # στα τμήματα του πτερυγίου (ΒΗΜΑ 8 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        v_axial = self.wind_speed_V0 * (1 - a)
        v_tangential = self.rotation_speed * r * (1 + a_p)
        Vrel = np.sqrt((v_axial)**2 + (v_tangential)**2)
        L = 0.5 * self.air_density * (Vrel)**2 * Cl * chord
        D = 0.5 * self.air_density * (Vrel)**2 * Cd * chord
        pn = (L * np.cos(flow_angle)) + (D * np.sin(flow_angle))
        pt = (L * np.sin(flow_angle)) - (D * np.cos(flow_angle))
        return L, D, pn, pt
     
    def segment_calculation(self, wind_speed_V0, omega_rad_sec, r, chord, pitch_angle_deg, twist_deg, tc_ratio, f=0.3): # εκτέλεση του αλγορίθμου
        v0=wind_speed_V0
        w_rps = omega_rad_sec
        # για κάθε τμήμα του πτερυγίου
        a, a_p = 0, 0 # αρχικοποίηση των συντελεστών επαγωγής a και a' σε 0
        exit_flag = False
        counter = 0 # αρχικά ο μετρητής έχει την τιμή 0
        
        while not exit_flag:
            flow_angle_rad = self.calculation_of_flow_angle_rad(a, a_p, r)
            angle_of_attack_rad = self.calculation_of_local_angle_of_attack_rad(flow_angle_rad=flow_angle_rad, pitch_angle_deg=pitch_angle_deg, twist_deg=twist_deg)
            Cl, Cd = self.calculation_of_Cl_and_Cd(angle_of_attack_deg=np.degrees(angle_of_attack_rad), tc_ratio=tc_ratio)
            Cn, Ct = self.calculation_of_Cn_and_Ct(Cl, Cd, flow_angle_rad)
            a_new, a_p_new = self.calculation_of_updated_induction_factors(Cn, Ct, r, chord, flow_angle_rad)
            
            if abs(a - a_new) < self.tolerance and abs(a_p - a_p_new) < self.tolerance: # έχω σύγκλιση των τιμών
                exit_flag = True
            else: # συνεχίζω τον αλγόριθμο
                a = a * (1 - f) + f * a_new
                a_p = a_p * (1 - f) + f * a_p_new  
            counter += 1 # αύξηση της τιμής του μετρητή κατά 1 
            if counter > self.max_iter:
                raise Exception("Δεν έχω σύγκλιση")
            L, D, pn, pt = self.calculation_of_local_loads(r, a, a_p, chord, flow_angle_rad, Cl, Cd)
            
        return {
            # "r_i (m)": r,
            # "chord (m)": chord,
            # "pitch_angle (degrees)": pitch_angle,
            # "twist (degrees)": twist,
            "a": a,
            "a_p": a_p,
            "flow_angle (rads)": flow_angle_rad,
            "flow angle (degrees)": np.degrees(flow_angle_rad),
            "angle_of_attack (rads)": angle_of_attack_rad,
            "angle of attack": np.degrees(angle_of_attack_rad),
            "Cl": Cl,
            "Cd": Cd,
            "Cn": Cn,
            "Ct": Ct,
            "a_new": a_new,
            "a_p_new": a_p_new,
            "Lift (N/m)": L,
            "Drag (N/m)": D,
            "pn (N/m)": pn,
            "pt (N/m)": pt,
            "counter": counter
        }
# #%% 
# # Εκτέλεση του αλγορίθμου και για τα δύο είδη αεροτομών
# if __name__ == "__main__": 
        
#     # print(f"Συνολική Ισχύς του φτερού: {total_power:.2f} Watt")
#     # print(f"Συντελεστής Απόδοσης (Cp): {cp:.4f}")
#     # print(f"Συνολική Ισχύς της Ανεμογεννήτριας: {3*total_power:.2f} Watt")

