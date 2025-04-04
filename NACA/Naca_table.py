#%%
"""
Κλάση για τον υπολογισμό των αεροδυναμικών συντελεστών άνωσης και οπισθέλκουσας
για την αεροτομή Naca4415 με βάση τη γωνία προσβολής, σε αριθμό Reynolds = 200.000
"""

import csv
import numpy as np

class Naca_calc:
    def __init__(self, csv_data_file_Naca):
        self.data = {}
        self.load_data(csv_data_file_Naca)

    def load_data(self, csv_data_file_Naca):
        with open(csv_data_file_Naca, mode = 'r') as file:
            reader = csv.reader(file)
            next(reader) # Αγνοούμε την πρώτη γραμμή του αρχείου που είναι οι επικεφαλίδες από κάθε στήλη
            for row in reader:
                angle_of_attack = round(float(row[0]), 3)
                cl_value = float(row[1])
                cd_value = float(row[2])
                self.data[angle_of_attack] = {'Cl': cl_value, 'Cd': cd_value}
                
        self.sorted_angles = sorted(self.data.keys()) 
        """ ταξινόμηση των γωνιών προσβολής """

    def interpolate(self, angle_of_attack, val1, val2, angle1, angle2):
        return val1 + (val2 - val1) * (angle_of_attack - angle1) / (angle2 - angle1) 
    """
    από τη σχέση της γραμμικής παρεμβολής y = y1 + [(x-x1)/(x2-x1)]*(y2-y1)
    όπου val1 και val2 οι τιμές για τους αεροδυναμικούς συντελεστές cl και cd αντίστοιχα 
    """
    
    def get_nearest_angles(self, angle_of_attack): 
        """ μέθοδος που βρίσκει τις δύο πλησιέστερες γωνίες προσβολής """ 
        for i in range(len(self.sorted_angles) - 1):
            if self.sorted_angles[i] <= angle_of_attack <= self.sorted_angles[i + 1]:
                return self.sorted_angles[i], self.sorted_angles[i + 1]
        raise ValueError(f'Η γωνία προσβολής {angle_of_attack} είναι εκτός του εύρους δεδομένων που περιλαμβάνει το αρχείο csv.')

    def cl(self, angle_of_attack): 
        """ μέθοδος που θα μας επιστρέψει την τιμή του συντελεστή άνωσης cl """
        angle_of_attack = round(angle_of_attack, 3) # Στρογγυλοποίηση στα 3 δεκαδικά ψηφία
        if angle_of_attack in self.data:
            return self.data[angle_of_attack]['Cl']
        else:
            # Βρες τις δύο πλησιέστερες γωνίες προσβολής
            angle1, angle2 = self.get_nearest_angles(angle_of_attack)
            # Γραμμική παρεμβολή για το συντελεστή Cl
            cl1 = self.data[angle1]['Cl']
            cl2 = self.data[angle2]['Cl']
            return self.interpolate(angle_of_attack, cl1, cl2, angle1, angle2)

    def cd(self, angle_of_attack): 
        """ μέθοδος που θα μας επιστρέψει την τιμή του συντελεστή αντίστασης (ή οπισθέλκουσας) cd """
        angle_of_attack = round(angle_of_attack, 3) # Στρογγυλοποίηση στα 3 δεκαδικά ψηφία 
        if angle_of_attack in self.data:
            return self.data[angle_of_attack]['Cd']
        else:
            # Βρες τις δύο πλησιέστερες γωνίες προσβολής
            angle1, angle2 = self.get_nearest_angles(angle_of_attack)
            # Γραμμική παρεμβολή για το συντελεστή Cd
            cd1 = self.data[angle1]['Cd']
            cd2 = self.data[angle2]['Cd']
            return self.interpolate(angle_of_attack, cd1, cd2, angle1, angle2)

# naca4415 = Naca_calc('csv_data_file_Naca.csv')

# cl_value = naca4415.cl(8.55)
# cd_value = naca4415.cd(8.55)
# print(f"Για την συγκεκριμένη γωνία προσβολής οι συντελεστές είναι: Cl: {cl_value}, Cd: {cd_value}")

