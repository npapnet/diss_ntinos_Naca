import numpy as np

class HansenAlgorithm:
    def __init__(self, wind_speed_V0, R, angular_speed, B, air_density):
        self.wind_speed_V0 = wind_speed_V0 # ταχύτητα του ανέμου (σε m/sec)
        self.R = R
        self.angular_speed = angular_speed # γωνιακή ταχύτητα (σε rad/sec)
        self.B = B # αριθμός πτερυγίων 
        self.air_density = air_density # πυκνότητα του αέρα (σε kg/m^3)
        
    def initialize_induction_factors(self, a, a_accent): # αρχικοποίηση των συντελεστών a και a' στο 0 (ΒΗΜΑ 1 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        self.a = 0.0
        self.a_accent = 0.0
        
    def calculation_of_flow_angle(self, a, a_accent, r): # μέθοδος για τον υπολογισμό της γωνίας ροής φ (ΒΗΜΑ 2 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        # a_accent = a'
        v_axial = self.wind_speed_V0 * (1 - a) # αξονική ταχύτητα
        v_tangential = self.angular_speed * r * (1 + a_accent) # εφαπτομενική ταχύτητα 
        return np.arctan(v_axial / v_tangential)
    
    def calculation_of_local_angle_of_attack(self, phi, theta_p, beta): # μέθοδος για τον υπολογισμό της τοπικής γωνίας προσβολής (ΒΗΜΑ 3 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        # phi = η γωνία ροής φ
        # theta = θ
        # theta_p = η γωνία κλίσης θp
        # beta = η συστροφή του πτερυγίου β
        return phi - (theta_p + beta)
    
    def cl_and_cd(self, alpha): # πίνακας για τους συντελεστές άνωσης και οπισθέλκουσας Cl and Cd βάσει του alpha (γωνία προσβολής) (ΒΗΜΑ 4 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        # Θα πρέπει να εισάγω κάποιον πίνακα για τους συντελεστές αυτούς
        Cl = 1 # ΥΠΟΘΕΣΗ
        Cd = 1 # ΥΠΟΘΕΣΗ
        return Cl, Cd
    
    def calculation_of_Cn_and_Ct(self, Cl, Cd, phi): # μέθοδος για τον υπολογισμό των συντελεστών Cn(normal) και Ct(tangential) (ΒΗΜΑ 5 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
        Ct = Cl * np.sin(phi) - Cd * np.cos(phi)
        return Cn, Ct
    
    def updated_induction_factors(self, Cn, Ct, chord, r, phi): # μέθοδος για τον υπολογισμό των νέων τιμών των συντελεστών a και a΄ (ΒΗΜΑ 6 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        sigma = (self.B * chord) / (2 * np.pi * r)
        a_new = 1 / ((4 * np.sin(phi)**2) / (sigma * Cn) + 1)
        a_accent_new = 1 / ((4 * np.sin(phi)*np.cos(phi)) / (sigma * Ct) - 1)
        return a_new, a_accent_new
    
    def check_the_convergence(self, a, a_new, a_accent, a_accent_new, tolerance=1): # έλεγχος σύγκλισης των τιμών των συντελεστών a και a΄ (ΒΗΜΑ 7 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        return abs(a - a_new) < tolerance and abs(a_accent - a_accent_new) < tolerance
    
    def calculation_of_local_forces(self, Vrel, Cl, Cd, chord): # υπολογισμός των τοπικών δυνάμεων στα τμήματα του πτερυγίου (ΒΗΜΑ 8 ΤΟΥ ΑΛΓΟΡΙΘΜΟΥ)
        L = 0.5 * self.air_density * (Vrel)**2 * Cl * chord
        D = 0.5 * self.air_density * (Vrel)**2 * Cd * chord
        return L, D
    
    def calculation_of_relative_spped(self, a, a_accent, r): # υπολογισμός της σχετικής ταχύτητας Vrel
         v_axial = self.wind_speed_V0 * (1 - a)
         v_tangential = (1 + a_accent) * self.angular_speed * r
         return np.sqrt((v_axial)**2 + (v_tangential)**2)
     
    def simulation_of_algorithm(self, r, chord, theta_p, beta): # εκτέλεση του αλγόριθμου
        a, a_accent = self._initialize_induction_factors()
        while True:
            phi = self.calculation_of_flow_angle(a, a_accent, r)
            alpha = self.calculation_of_local_angle_of_attack(phi, theta_p, beta)
            Cl, Cd = self.cl_and_cd(alpha)
            Cn, Ct = self.calculation_of_Cn_and_Ct(Cl, Cd, phi)
            a_new, a_accent_new = self.updated_induction_factors(Cn, Ct, chord, r, phi)
            
            if self.check_the_convergence(a, a_new, a_accent, a_accent_new):
                break
            else:
                a, a_accent = a_new, a_accent_new  
        Vrel = self.calculation_of_relative_spped(a, a_accent, r)
        L, D = self.calculation_of_local_forces(Vrel, Cl, Cd, r, phi, chord)
        
        return {
            "a" : a,
            "a_accent" : a_accent,
            "phi" : phi,
            "alpha" : alpha,
            "Cl" : Cl,
            "Cd" : Cd,
            "Cn" : Cn,
            "Ct" : Ct,
            "Vrel" : Vrel,
            "Lift" : L,
            "Drag" : D  
        }
    
