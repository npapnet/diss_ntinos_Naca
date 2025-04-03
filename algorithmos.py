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
        """ 
        μέθοδος αρχικοποίησης των βασικών μεταβλητών 

        Args:
            R (float): η ακτίνα του ρότορα σε m
            B (int, optional): ο αριθμός των πτερυγίων. Defaults to 3.
            air_density (float, optional): η πυκνότητα του αέρα σε kg/m^3. Defaults to 1.225.
            airfoil_type (string, optional): ο τύπος της αεροτομής που χρησιμοποιείται. Defaults to None.
            csv_data_file (_type_, optional): το αρχείο csv που περιέχει τα δεδομένα για τους συντελεστές Cl και Cd. Defaults to None.
        """
        # self.wind_speed_V0 = wind_speed_V0 # ταχύτητα του ανέμου (σε m/sec)
        self.R = R 
        # self.rotation_speed = rotation_speed # ταχύτητα περιστροφής του ρότορα (σε rad/sec)
        self.B = B 
        self.air_density = air_density 
        if airfoil_type == "NACA":
            self.airfoil_calc = Naca_calc(csv_data_file) # Χρήση NACA δεδομένων
        elif airfoil_type == "DTU":
            self.airfoil_calc = DTU_calc(csv_data_file) # Χρήση DTU δεδομένων
        else:
            self.airfoil_calc = None
        
    def calculation_of_flow_angle_rad(self, a, a_p, r:float, v0:int, w_rps:float):
        """
        μέθοδος για τον υπολογισμό της γωνίας ροής φ (ΒΗΜΑ 2 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)

        Args:
            a (float): ο συντελεστής αξονικής επαγωγής α 
            a_p (float): ο συντελεστής εφαπτομενικής επαγωγής α΄
            r (float): η ακτίνα του κάθε τμήματος του πτερυγίου σε m
            v0 (int): η ταχύτητα του ανέμου σε m/sec
            w_rps (float): η ταχύτητα περιστροφής του ρότορα σε rad/sec

        Raises:
            ValueError: Εάν η τιμή του συντελεστή εφαπτομενικής επαγωγής α΄ έχει την τιμή -1
             τότε η εφαπτομενική ταχύτητα που βρίσκεται στον παρονομαστή γίνεται μηδέν και δεν 
             είναι δυνατόν να γίνει διαίρεση κλάσματος με παρονομαστή το μηδέν

        Returns:
            (float): η γωνία της ροής φ σε rad
        """
        v_axial = v0 * (1 - a) # συνιστώσα αξονικής ταχύτητας
        v_tangential = w_rps* r * (1 + a_p) # συνιστώσα εφαπτομενικής ταχύτητας

        if a_p == -1:
            raise ValueError("Διαίρεση με το 0")
        else:
            return np.arctan(v_axial / v_tangential)
    
    def calculation_of_local_angle_of_attack_rad(self, flow_angle_rad:float, pitch_angle_deg:float, twist_deg:float): 
        """ 
        μέθοδος για τον υπολογισμό της τοπικής γωνίας προσβολής (ΒΗΜΑ 3 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)

        Args:
            flow_angle_rad (float): η γωνία της ροής φ σε rad
            pitch_angle_deg (float): η γωνία βήματος του πτερυγίου σε μοίρες
            twist_deg (float):η γωνία συστροφής του πτερυγίου β σε μοίρες

        Returns:
             (float): η γωνία προσβολής σε rad
        """
        return (flow_angle_rad - (np.radians(pitch_angle_deg) + np.radians(twist_deg)))
    
    def calculation_of_Cl_and_Cd(self, angle_of_attack_deg:float, tc_ratio=None):
        """ 
        πίνακας για τους συντελεστές άνωσης και οπισθέλκουσας Cl and Cd
        βάσει του angle_of_attack (γωνία προσβολής), αλλά και του λόγου t/c (thickness/chord ratio) (ΒΗΜΑ 4 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        
        Args:
            angle_of_attack_deg (float): η γωνία προσβολής α σε μοίρες
            tc_ratio (float): o λόγος πάχους αεροτομής / χορδή αεροτομής
            
        Returns:
             (float): αεροδυναμικούς συντελεστές άνωσης και οπισθέλκουσας Cl και Cd
        """
        if isinstance(self.airfoil_calc, Naca_calc):
            Cl = self.airfoil_calc.cl(angle_of_attack_deg) # Υπολογισμός Cl μόνο με βάση τη γωνία προσβολής
            Cd = self.airfoil_calc.cd(angle_of_attack_deg) # Υπολογισμός Cd μόνο με βάση τη γωνία προσβολής
        elif isinstance(self.airfoil_calc, DTU_calc):
            Cl = self.airfoil_calc.cl(angle_of_attack_deg, tc_ratio) # Υπολογισμός Cl με βάση τη γωνία προσβολής και t/c
            Cd = self.airfoil_calc.cd(angle_of_attack_deg, tc_ratio) # Υπολογισμός Cd με βάση τη γωνία προσβολής και t/c
        else:
            raise ValueError("Άγνωστος τύπος αεροτομής")
        return Cl, Cd
    
    def calculation_of_Cn_and_Ct(self, Cl:float, Cd:float, flow_angle_rad:float):
        """ 
        μέθοδος για τον υπολογισμό των συντελεστών 
        Cn(normal) και Ct(tangential) (ΒΗΜΑ 5 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)

        Args:
            Cl (float): ο συντελεστής άνωσης Cl
            Cd (float): o συντελεστής οπισθέλκουσας Cd
            flow_angle_rad (float): η γωνία της ροής φ σε rad

        Returns:
             (float): συντελεστές Cn(normal) και Ct(tangential)
        """
        
        Cn = (Cl * np.cos(flow_angle_rad)) + (Cd * np.sin(flow_angle_rad))
        Ct = (Cl * np.sin(flow_angle_rad)) - (Cd * np.cos(flow_angle_rad))
        return Cn, Ct
    
    def calculation_of_updated_induction_factors(self, Cn:float, Ct:float, r:float, chord:float, flow_angle_rad:float):
        """
        μέθοδος για τον υπολογισμό των νέων τιμών 
        των συντελεστών επαγωγής a και a΄ (ΒΗΜΑ 6 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)

        Args:
            Cn (float): o συντελεστής Cn
            Ct (float): o συντελεστής Ct
            r (float): η ακτίνα του κάθε τμήματος του πτερυγίου σε m
            chord (float): το μήκος χορδής του κάθε τμήματος του πτερυγίου σε m
            flow_angle_rad (float): η γωνία της ροής φ σε rad

        Returns:
             (float): οι νέοι συντελεστές αξονικής και εφαπτομενικής επαγωγής α και α΄
        """
        
        solidity_factor = (self.B * chord) / (2 * np.pi * r) # solidity factor (συντελεστής στερεότητας)
        a_new = 1 / ((4 * np.sin(flow_angle_rad)**2) / (solidity_factor * Cn) + 1)
        a_p_new = 1 / ((4 * np.sin(flow_angle_rad)*np.cos(flow_angle_rad)) / (solidity_factor * Ct) - 1)
        return a_new, a_p_new

    def calculation_of_local_loads(self, r:float, a:float, a_p:float, v0:int, w_rps:float,
                                   chord:float, flow_angle_rad:float, Cl:float, Cd:float):
        """
        μέθοδος για τον υπολογισμό των τοπικών φορτίων
        στα τμήματα του πτερυγίου (ΒΗΜΑ 8 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)

        Args:
            r (float): η ακτίνα του κάθε τμήματος του πτερυγίου σε m
            a (float): ο συντελεστής αξονικής επαγωγής α 
            a_p (float): ο συντελεστής εφαπτομενικής επαγωγής α΄
            chord (float): το μήκος χορδής του κάθε τμήματος του πτερυγίου σε m
            flow_angle_rad (float): η γωνία της ροής φ σε rad
            Cl (float): ο συντελεστής άνωσης Cl
            Cd (float): o συντελεστής οπισθέλκουσας Cd

        Returns:
             (float): οι δυνάμεις άνωσης και οπισθέλουσας L και D καθώς επίσης και οι συνιστώσες pn και pt
        """
        
        v_axial = v0 * (1 - a)
        v_tangential = w_rps * r * (1 + a_p)
        Vrel = np.sqrt((v_axial)**2 + (v_tangential)**2)
        L = 0.5 * self.air_density * (Vrel)**2 * Cl * chord
        D = 0.5 * self.air_density * (Vrel)**2 * Cd * chord
        pn = (L * np.cos(flow_angle_rad)) + (D * np.sin(flow_angle_rad))
        pt = (L * np.sin(flow_angle_rad)) - (D * np.cos(flow_angle_rad))
        return L, D, pn, pt
     
    def segment_calculation(self, wind_speed_V0, omega_rad_sec, r, chord, pitch_angle_deg, twist_deg, tc_ratio, f=0.3):
        """ εκτέλεση του αλγορίθμου για κάθε τμήμα του πτερυγίου """
        v0 = wind_speed_V0
        w_rps = omega_rad_sec
        
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
            "r_i (m)": r,
            "chord (m)": chord,
            "pitch_angle (degrees)": pitch_angle_deg,
            "twist (degrees)": twist_deg,
            "a": a,
            "a_p": a_p,
            "flow_angle (rads)": flow_angle_rad,
            "flow angle (degrees)": np.degrees(flow_angle_rad),
            "angle_of_attack (rads)": angle_of_attack_rad,
            "angle of attack (degrees)": np.degrees(angle_of_attack_rad),
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

