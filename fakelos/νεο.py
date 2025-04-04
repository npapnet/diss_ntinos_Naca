import numpy as np
from fakelos.naca4415 import Naca4415

class HansenAlgorithm:
    tolerance = 1e-6
    max_iter = 10000
    
    def __init__(self, wind_speed_V0, R, angular_speed, B=3, air_density=1.225, csv_data_file='csv_data_file.csv'):
        self.wind_speed_V0 = wind_speed_V0
        self.R = R
        self.angular_speed = angular_speed
        self.B = B
        self.air_density = air_density
        self.naca4415 = Naca4415(csv_data_file)

    def calculation_of_flow_angle(self, a, a_accent, r):
        v_axial = self.wind_speed_V0 * (1 - a)
        v_tangential = self.angular_speed * r * (1 + a_accent)
        return np.arctan(v_axial / v_tangential)
    
    def calculation_of_local_angle_of_attack(self, phi, theta_p, twist):
        return phi - (theta_p + twist)
    
    def cl_and_cd(self, alpha):
        Cl = self.naca4415.cl(alpha)
        Cd = self.naca4415.cd(alpha)
        return Cl, Cd
    
    def calculation_of_Cn_and_Ct(self, Cl, Cd, phi):
        Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
        Ct = Cl * np.sin(phi) - Cd * np.cos(phi)
        return Cn, Ct
    
    def updated_induction_factors(self, Cn, Ct, r, chord, phi):
        sigma = (self.B * chord) / (2 * np.pi * r)
        a_new = 1 / ((4 * np.sin(phi)**2) / (sigma * Cn) + 1)
        a_accent_new = 1 / ((4 * np.sin(phi) * np.cos(phi)) / (sigma * Ct) - 1)
        return a_new, a_accent_new

    def calculation_of_local_forces(self, r, a, a_accent, chord, Cl, Cd):
        v_axial = self.wind_speed_V0 * (1 - a)
        v_tangential = self.angular_speed * r * (1 + a_accent)
        Vrel = np.sqrt(v_axial**2 + v_tangential**2)
        L = 0.5 * self.air_density * Vrel**2 * Cl * chord
        D = 0.5 * self.air_density * Vrel**2 * Cd * chord
        return L, D

    def run_algorithm_for_section(self, r, theta_p, twist, chord):
        a, a_accent = 0, 0
        exit_flag = False
        counter = 0

        while not exit_flag:
            phi = self.calculation_of_flow_angle(a, a_accent, r)
            alpha = self.calculation_of_local_angle_of_attack(phi, theta_p, twist)
            Cl, Cd = self.cl_and_cd(alpha)
            Cn, Ct = self.calculation_of_Cn_and_Ct(Cl, Cd, phi)
            a_new, a_accent_new = self.updated_induction_factors(Cn, Ct, r, chord, phi)
            
            if abs(a - a_new) < self.tolerance and abs(a_accent - a_accent_new) < self.tolerance:
                exit_flag = True
            else:
                a, a_accent = a_new, a_accent_new
            counter += 1
            if counter > self.max_iter:
                raise Exception("Δεν έχω σύγκλιση στο τμήμα του φτερού")

        L, D = self.calculation_of_local_forces(r, a, a_accent, chord, Cl, Cd)
        
        return {
            "a": a,
            "a_accent": a_accent,
            "phi": phi,
            "alpha": alpha,
            "Cl": Cl,
            "Cd": Cd,
            "Cn": Cn,
            "Ct": Ct,
            "Lift": L,
            "Drag": D
        }

    def run_algorithm_for_blade(self, R_start, R_end, chord_distribution, twist_distribution, theta_p):
        sections = 10
        section_results = []
        r_values = np.linspace(R_start, R_end, sections)

        for i, r in enumerate(r_values):
            chord = chord_distribution[i]  # Χορδή για το τμήμα i
            twist = twist_distribution[i]  # Συστροφή για το τμήμα i
            result = self.run_algorithm_for_section(r, theta_p, twist, chord)
            result['r'] = r
            section_results.append(result)
            print(f"Αποτελέσματα για το τμήμα {i + 1} (r = {r:.2f} m):")
            for key, value in result.items():
                print(f"  {key}: {value}")
            print("\n")
        
        return section_results

# Παράδειγμα χρήσης
# Ορισμός κατανομής χορδής και συστροφής για κάθε τμήμα του φτερού
R_start, R_end = 3, 50  # Ξεκινώντας και τελική ακτίνα του φτερού
chord_distribution = np.linspace(3, 1, 10)  # Κατανομή χορδής από τη βάση ως την άκρη
twist_distribution = np.linspace(0.2, 0, 10)  # Κατανομή συστροφής από τη βάση ως την άκρη
theta_p = 2  # Γωνία κλίσης του φτερού

# Εκτέλεση του αλγορίθμου για ολόκληρο το φτερό
hansen = HansenAlgorithm(wind_speed_V0=8, R=50, angular_speed=2, B=3, air_density=1.225, csv_data_file='csv_data_file.csv')
blade_results = hansen.run_algorithm_for_blade(R_start, R_end, chord_distribution, twist_distribution, theta_p)
