#%%
"""
Κλάση για τον υπολογισμό των αεροδυναμικών συντελεστών άνωσης και οπισθέλκουσας
για την αεροτομή DTU με βάση τη γωνία προσβολής αλλά και το λόγο t/c (thickness / chord ratio)
"""

import csv
import numpy as np

class DTU_calc:
    def __init__(self, csv_data_file_DTU):
        """ 
        Δημιουργούμε ένα λεξικό που θα αποθηκεύει τα δεδομένα της αεροτομής. Στη συνέχεια
        χρησιμοποιούμε ένα σύνολο (set) για να αποθηκεύουμε τις μοναδικές τιμές του λόγου t/c που υπάρχουν στο αρχείο CSV. 
        Τέλος καλούμε τη μέθοδο load_data, η οποία θα διαβάσει το αρχείο CSV και θα γεμίσει τo self.data.
        """
        self.data = {} 
        self.tc_values = set() 
        self.load_data(csv_data_file_DTU) 

    def load_data(self, csv_data_file_DTU):
        with open(csv_data_file_DTU, mode='r') as file:
            reader = csv.reader(file, delimiter=';')
            next(reader) # Αγνοούμε την πρώτη γραμμή του αρχείου που είναι οι επικεφαλίδες από κάθε στήλη
            for row in reader:
                angle_of_attack = round(float(row[0]), 3)
                cl_value = float(row[1])
                cd_value = float(row[2])
                cm_value = float(row[3])
                tc_ratio = float(row[4]) # Ο λόγος t/c
                
                """ 
                Αν το tc_ratio δεν υπάρχει ήδη στο self.data, δημιουργούμε ένα νέο υπολεξικό για αυτό.
                Αποθηκεύουμε στο λεξικό self.data τα δεδομένα των συντελεστών Cl, Cd, Cm με βάση τη γωνία προσβολής.
                Προσθέτουμε την τιμή του tc_ratio στο σύνολο self.tc_values
                """
                if tc_ratio not in self.data: 
                    self.data[tc_ratio] = {}
                self.data[tc_ratio][angle_of_attack] = {'Cl': cl_value, 'Cd': cd_value, 'Cm': cm_value} 
                self.tc_values.add(tc_ratio) 
        
        """ 
        Ταξινόμηση των λόγων t/c, δηλαδή ταξινομούμε τη λίστα self.tc_values, ώστε οι τιμές t/c να είναι με αύξουσα σειρά.
        Ταξινόμηση των γωνιών προσβολής για κάθε t/c, δηλαδή για κάθε tc_ratio ταξινομούμε τις γωνίες προσβολής (angle_of_attack) μέσα στο λεξικό.
        """
        self.tc_values = sorted(self.tc_values) 
        for tc in self.data:
            self.data[tc] = dict(sorted(self.data[tc].items())) 

    def interpolate(self, x, x1, x2, y1, y2):
        return y1 + (y2 - y1) * ((x - x1) / (x2 - x1)) 
    """
    από τη σχέση της γραμμικής παρεμβολής y = y1 + [(x-x1)/(x2-x1)]*(y2-y1
    """
    
    def get_nearest_value(self, values, desired_value): 
        """ μέθοδος που βρίσκει τα δύο κοντινότερα διαθέσιμα σημεία στο dataset που περιβάλλουν την επιθυμητή τιμή (desired_value) """
        for i in range(len(values) - 1):
            if values[i] <= desired_value <= values[i + 1]:
                return values[i], values[i + 1]
        raise ValueError(f'Η τιμή {desired_value} είναι εκτός του εύρους δεδομένων που περιλαμβάνει το αρχείο csv.')
    
    def get_interpolated_value(self, angle_of_attack, tc_ratio, coef):
        """ 
        Βρίσκουμε τα δύο πλησιέστερα t/c γύρω από το tc_ratio και στη συνέχεια
        βρίσκουμε τις δύο πλησιέστερες γωνίες προσβολής γύρω από το angle_of_attack
        """
        tc1, tc2 = self.get_nearest_value(self.tc_values, tc_ratio) 
        angle1, angle2 = self.get_nearest_value(list(self.data[tc1].keys()), angle_of_attack) 
        
        # Παρεμβολή πρώτα ως προς τη γωνία προσβολής
        coef1 = self.interpolate(angle_of_attack, angle1, angle2, 
                                 self.data[tc1][angle1][coef], self.data[tc1][angle2][coef])
        coef2 = self.interpolate(angle_of_attack, angle1, angle2, 
                                 self.data[tc2][angle1][coef], self.data[tc2][angle2][coef])
        
        # Παρεμβολή ως προς το λόγο t/c
        return self.interpolate(tc_ratio, tc1, tc2, coef1, coef2)

    def cl(self, angle_of_attack, tc_ratio): 
        """ μέθοδος που θα μας επιστρέψει την τιμή του συντελεστή άνωσης Cl """
        return self.get_interpolated_value(angle_of_attack, tc_ratio, 'Cl')

    def cd(self, angle_of_attack, tc_ratio): 
        """ μέθοδος που θα μας επιστρέψει την τιμή του συντελεστή αντίστασης (ή οπισθέλκουσας) Cd """
        return self.get_interpolated_value(angle_of_attack, tc_ratio, 'Cd') 
    
    def cm(self, angle_of_attack, tc_ratio): 
        """ μέθοδος που θα μας επιστρέψει την τιμή του συντελεστή ροπής Cm """
        return self.get_interpolated_value(angle_of_attack, tc_ratio, 'Cm')

# # Παράδειγμα χρήσης
# naca4415 = DTU_calc('csv_data_file_DTU.csv')
# a_o_a = 0.119217
# tc = 24.1

# cl_value = naca4415.cl(a_o_a, tc)
# cd_value = naca4415.cd(a_o_a, tc)
# cm_value = naca4415.cm(a_o_a, tc)

# print(f"Για γωνία προσβολής {a_o_a}° και t/c ratio {tc} οι συντελεστές είναι:")
# print(f"Cl: {cl_value}, Cd: {cd_value}, Cm: {cm_value}")



# %%
