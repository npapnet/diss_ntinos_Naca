import numpy as np
import json
import pandas as pd
import logging
from hansen_algorithm import HansenAlgorithm  # Εισαγωγή της κλάσης του αλγορίθμου Hansen

class BladeGeometry:
    _bg_version: str = "1.0"
    R: float = None
    bl_Ri: np.ndarray = None
    bl_c_i: np.ndarray = None
    _no_sections: int = None
    pitch: np.ndarray = None  # γωνία κλίσης (pitch angle)
    lambda0: float = None  # lambda στην άκρη του φτερού

    def __init__(self, R: float, r_is: np.ndarray, lambda0: float, chords: np.ndarray, pitch: np.ndarray, no_sections: int = None):
        self.R = R
        self.bl_c_i = chords
        self.bl_Ri = r_is
        self._no_sections = len(chords)
        
        # Διασφάλιση ότι η γωνία κλίσης έχει τον σωστό τύπο και μέγεθος
        if pitch is None:
            pitch = np.zeros(self._no_sections)
        elif isinstance(pitch, (int, float)):
            pitch = np.ones(self._no_sections) * pitch
        elif isinstance(pitch, np.ndarray) and len(pitch) == self._no_sections:
            self.pitch = pitch
        else:
            raise ValueError("Η pitch πρέπει να είναι float, int ή numpy array με το ίδιο μέγεθος με τα chords.")
        self.lambda0 = lambda0

    @property  
    def no_sections(self):
        """Αριθμός τμημάτων του φτερού"""
        return self._no_sections

    def run_hansen_algorithm(self, wind_speed_V0, angular_speed, air_density=1.225, theta_p=0):
        """Εκτέλεση αλγορίθμου Hansen για κάθε τμήμα του φτερού"""
        results = []
        for i in range(self._no_sections):
            r = self.bl_Ri[i]
            chord = self.bl_c_i[i]
            twist = self.pitch[i]
            hansen = HansenAlgorithm(wind_speed_V0, self.R, angular_speed, air_density=air_density)
            result = hansen.run_algorithm_for_section(r, theta_p, twist, chord)
            results.append(result)
        
        return results

    # Υπόλοιπες μέθοδοι της BladeGeometry (π.χ., to_dict, to_df, to_json) παραμένουν οι ίδιες

# Παράδειγμα χρήσης
R_start, R_end = 3, 50  # Ακτίνες για τις διατομές
chord_distribution = np.linspace(3, 1, 10)  # Κατανομή χορδής
twist_distribution = np.linspace(0.2, 0, 10)  # Κατανομή γωνίας στρέψης
theta_p = 2  # Γωνία κλίσης φτερού

# Δημιουργία αντικειμένου BladeGeometry και εκτέλεση του αλγορίθμου
bg = BladeGeometry(R=50, r_is=np.linspace(R_start, R_end, 10), lambda0=8, chords=chord_distribution, pitch=twist_distribution)
blade_results = bg.run_hansen_algorithm(wind_speed_V0=8, angular_speed=2)

# Εκτύπωση αποτελεσμάτων
for i, result in enumerate(blade_results):
    print(f"Αποτελέσματα για το τμήμα {i + 1}:")
    for key, value in result.items():
        print(f"  {key}: {value}")
    print("\n")
