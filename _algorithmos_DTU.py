#%% 
import json
import pandas as pd
import numpy as np
from Dtu_table import DTU_calc
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

#%%
def blade_geometry_seperation(r_is, chords, pitch, tc_ratios,
                            r_first, r_last, num_sections=10):
    """
    Μέθοδος που θα διαβάσει το αρχείο json με τα δεδομένα της αεροτομής και δημιουργεί ένα νέο πίνακα ακτινικών θέσεων r_new με μήκος num_sections
    (δηλ. 10 σημεία που μελετάμε στην προκειμένη περίπτωση) γραμμικά κατανεμημένα από r_first μέχρι r_last, σε ίση απόσταση περίπου 10 m μεταξύ τους.
    Στη συνέχεια υπολογίζει με γραμμική παρεμβολή τις νέες τιμές chord, pitch, tc_ratio στα r_new. Με τη συγκεκριμένη μέθοδο είναι εφικτό να γίνει ο 
    διαχωρισμός του φτερού σε όσα τμήματα επιθυμούμε.
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
class Hansen_Algorithm:
    """
    Κλάση που εμπεριέχει όλα τα βήματα του αλγορίθμου της αεροδυναμικής θεωρίας Blade Element Momentum (θεωρία στοιχείων πτερύγωσης - ορμής)
    για την ανάλυση της αεροδυναμικής συμπεριφοράς πτερυγίων ανεμογεννητριών.
    """
    tolerance = 1e-4 # Σταθερά για τον έλεγχο σύγκλισης των τιμών των συντελεστών αξονικής και εφαπτομενικής επαγωγής a και a'
    max_iter = 100 # Μέγιστος αριθμός επαναλήψεων κατά τη διαδικασία σύγκλισης
    
    def __init__(self, blade_geom_DTU, B=3, air_density=1.225, airfoil_type=None, csv_data_file='csv_data_file_DTU.csv'):
        """ 
        μέθοδος αρχικοποίησης των βασικών μεταβλητών 

        Args:
            R (float): η ακτίνα του ρότορα σε m
            B (int, optional): ο αριθμός των πτερυγίων. Defaults to 3.
            air_density (float, optional): η πυκνότητα του αέρα σε kg/m^3. Defaults to 1.225.
            airfoil_type (string, optional): ο τύπος της αεροτομής που χρησιμοποιείται. Defaults to None.
            csv_data_file (csv file, optional): το αρχείο csv που περιέχει τα δεδομένα για τους αεροδυναμικούς συντελεστές Cl και Cd. Defaults to None.
        """
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
        r_new, chords_new, pitch_new, tc_new = blade_geometry_seperation(
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

        # self.wind_speed_V0 = wind_speed_V0 # ταχύτητα του ανέμου (σε m/sec)
        self.R = r_last
        # self.rotation_speed = rotation_speed # ταχύτητα περιστροφής του ρότορα (σε rad/sec)
        self.B = B 
        self.air_density = air_density 
        self.airfoil_calc = DTU_calc(csv_data_file) # Χρήση DTU δεδομένων

    def calculation_of_flow_angle_rad(self, a, a_p, r:float, v0:int, w_rps:float):
        """
        μέθοδος για τον υπολογισμό της γωνίας ροής φ σε rad (ΒΗΜΑ 2 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)

        Args:
            a (float): ο συντελεστής αξονικής επαγωγής α 
            a_p (float): ο συντελεστής εφαπτομενικής επαγωγής α΄
            r (float): η ακτίνα του κάθε τμήματος του πτερυγίου σε m
            v0 (int): η ταχύτητα του ανέμου σε m/sec
            w_rps (float): η ταχύτητα περιστροφής του ρότορα σε rad/sec

        Raises:
            ValueError: Εάν η τιμή του συντελεστή εφαπτομενικής επαγωγής α΄ έχει την τιμή -1
            τότε η εφαπτομενική ταχύτητα γίνεται μηδέν και δεν είναι δυνατόν να γίνει διαίρεση 
            κλάσματος με παρονομαστή το μηδέν

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
        μέθοδος για τον υπολογισμό της τοπικής γωνίας προσβολής σε rad (ΒΗΜΑ 3 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)

        Args:
            flow_angle_rad (float): η γωνία της ροής φ σε rad
            pitch_angle_deg (float): η γωνία βήματος του πτερυγίου σε μοίρες
            twist_deg (float): η γωνία συστροφής του πτερυγίου β σε μοίρες

        Returns:
             (float): η γωνία προσβολής σε rad
        """
        return (flow_angle_rad - (np.radians(pitch_angle_deg) + np.radians(twist_deg)))
    
    def calculation_of_Cl_and_Cd(self, angle_of_attack_deg:float, tc_ratio=None):
        """ 
        πίνακας για τους αεροδυναμικούς συντελεστές άνωσης και οπισθέλκουσας Cl and Cd
        βάσει του angle_of_attack (γωνία προσβολής), αλλά και του λόγου t/c (thickness/chord ratio) (ΒΗΜΑ 4 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        
        Args:
            angle_of_attack_deg (float): η γωνία προσβολής α σε μοίρες
            tc_ratio (float): o λόγος πάχους αεροτομής / χορδή αεροτομής
            
        Returns:
             (float): αεροδυναμικούς συντελεστές άνωσης και οπισθέλκουσας Cl και Cd
        """
        Cl = self.airfoil_calc.cl(angle_of_attack_deg, tc_ratio) # Υπολογισμός Cl με βάση τη γωνία προσβολής και t/c
        Cd = self.airfoil_calc.cd(angle_of_attack_deg, tc_ratio) # Υπολογισμός Cd με βάση τη γωνία προσβολής και t/c
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
            chord (float): το μήκος της χορδής του κάθε τμήματος του πτερυγίου σε m
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
            chord (float): το μήκος της χορδής του κάθε τμήματος του πτερυγίου σε m
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
     
    def segment_calculation(self, wind_speed_V0, omega_rad_sec, r, chord, pitch_angle_deg, twist_deg, tc_ratio, 
                            f=0.3, debug_mode=False):
        """ εκτέλεση του αλγορίθμου για κάθε τμήμα του πτερυγίου """
        v0 = wind_speed_V0
        w_rps = omega_rad_sec
        
        _debug_lst = []
        a, a_p = 0, 0 # αρχικοποίηση των συντελεστών επαγωγής a και a' σε 0
        exit_flag = False
        counter = 0 # αρχικά ο μετρητής έχει την τιμή 0
        
        
        while not exit_flag:
            flow_angle_rad = self.calculation_of_flow_angle_rad(a=a, a_p=a_p, r=r, v0=wind_speed_V0, w_rps=omega_rad_sec)
            angle_of_attack_rad = self.calculation_of_local_angle_of_attack_rad(flow_angle_rad=flow_angle_rad, pitch_angle_deg=pitch_angle_deg, twist_deg=twist_deg)
            Cl, Cd = self.calculation_of_Cl_and_Cd(angle_of_attack_deg=np.degrees(angle_of_attack_rad), tc_ratio=tc_ratio)
            Cn, Ct = self.calculation_of_Cn_and_Ct(Cl=Cl, Cd=Cd, flow_angle_rad=flow_angle_rad)
            a_new, a_p_new = self.calculation_of_updated_induction_factors(Cn=Cn, Ct=Ct, r=r, chord=chord, flow_angle_rad=flow_angle_rad)
            
            if abs(a - a_new) < self.tolerance and abs(a_p - a_p_new) < self.tolerance: # έχω σύγκλιση των τιμών
                exit_flag = True
            else: # συνεχίζω τον αλγόριθμο
                a = a * (1 - f) + f * a_new
                a_p = a_p * (1 - f) + f * a_p_new  
            counter += 1 # αύξηση της τιμής του μετρητή κατά 1 
            if counter > self.max_iter:
                break
                raise Exception("Δεν έχω σύγκλιση")
            L, D, pn, pt = self.calculation_of_local_loads(r=r, a=a, a_p=a_p, v0=v0, w_rps=omega_rad_sec, chord=chord, flow_angle_rad=flow_angle_rad, Cl=Cl, Cd=Cd)
            temp_dict = {
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
                "counter":counter
            }
            _debug_lst.append(temp_dict)
            
        if debug_mode:
            df = pd.DataFrame(_debug_lst)
            df.to_excel("save_res.xlsx")
              
        res_dict =  {
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
        return res_dict

    def DTU_blade_calculation(self, wind_speed_V0, rotation_speed):
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
                power = rotation_speed * dM
                total_power += power # Συνολική ισχύς όλου του ρότορα 
                total_torque += dM # Συνολική ροπή του ρότορα
                total_thrust += dT # Συνολική ώση του ρότορα
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
    
    def calculation_of_coefficient_of_thrust_CT_for_DTU(self, total_thrust, wind_speed_V0):
        swept_area = np.pi * self.R**2 # επιφάνεια σάρωσης
        wind_force = 0.5 * self.air_density * swept_area * wind_speed_V0**2
        return total_thrust / wind_force

# #%% 
# if __name__ == "__main__":  
#     # print(f"Συνολική Ισχύς του φτερού: {total_power:.2f} Watt")
#     # print(f"Συντελεστής Απόδοσης (Cp): {cp:.4f}")
#     # print(f"Συνολική Ισχύς της Ανεμογεννήτριας: {3*total_power:.2f} Watt")
#%%
if __name__ == "__main__":
    blade_geom_DTU = "blade_geom_DTU.json"
    hansen = Hansen_Algorithm(
        blade_geom_DTU=blade_geom_DTU,
        B=3,
        air_density=1.225,
        csv_data_file="csv_data_file_DTU.csv" 
    )

    # Καλούμε τη μέθοδο DTU_blade_calculation για μια συγκεκριμένη ταχύτητα ανέμου & ταχύτητα περιστροφής
    results_for_DTU_geometry, total_power, total_torque, total_thrust = hansen.DTU_blade_calculation(
        wind_speed_V0=10, 
        rotation_speed=1.3
    )

    df_DTU_results = pd.DataFrame(results_for_DTU_geometry)
    print(df_DTU_results)

    # print(df_DTU_results.head())
    # print(df_DTU_results.columns)
    


