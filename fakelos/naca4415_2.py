#Κλάση για την αεροτομή Naca

import csv

class Naca4415:
    def __init__(self, csv_data_file):
        self.data = {}
        self.load_data(csv_data_file)

    def load_data(self, csv_data_file):
        with open(csv_data_file, mode='r') as file:
            reader = csv.reader(file, delimiter=';')
            next(reader) # Αγνοούμε την πρώτη γραμμή του αρχείου που είναι οι επικεφαλίδες από κάθε στήλη
            for row in reader:
                angle_of_attack = round(float(row[0]), 3)
                cl_value = float(row[1])
                cd_value = float(row[2])
                cm_value = float(row[3]) 
                tc_ratio = float(row[4]) # Ο λόγος t/c
                self.data[angle_of_attack] = {'Cl': cl_value, 'Cd': cd_value, 'Cm': cm_value, 't/c ratio': tc_ratio}

        self.sorted_angles = sorted(self.data.keys()) # Ταξινόμηση των γωνιών προσβολής

    def interpolate(self, angle_of_attack, val1, val2, angle1, angle2):
        # Γραμμική παρεμβολή για να υπολογίσουμε τις τιμές μεταξύ των γωνιών προσβολής
        return val1 + (val2 - val1) * (angle_of_attack - angle1) / (angle2 - angle1) # από τη σχέση της γραμμικής παρεμβολής y = y1 + [(x-x1)/(x2-x1)]*(y2-y1)
    
    def get_nearest_angles(self, angle_of_attack): # μέθοδος που βρίσκει τις δύο πλησιέστερες γωνίες προσβολής 
        for i in range(len(self.sorted_angles) - 1):
            if self.sorted_angles[i] <= angle_of_attack <= self.sorted_angles[i + 1]:
                return self.sorted_angles[i], self.sorted_angles[i + 1]
        raise ValueError(f'Η γωνία προσβολής {angle_of_attack} είναι εκτός του εύρους.')

    def cl(self, angle_of_attack): # μέθοδος που θα μας επιστρέψει την τιμή του συντελεστή άνωσης Cl
        angle_of_attack = round(angle_of_attack, 3) # Στρογγυλοποίηση στα 3 δεκαδικά ψηφία
        if angle_of_attack in self.data:
            return self.data[angle_of_attack]['Cl']
        else:
            # Βρες τις δύο πλησιέστερες γωνίες προσβολής
            angle1, angle2 = self.get_nearest_angles(angle_of_attack)
            cl1 = self.data[angle1]['Cl']
            cl2 = self.data[angle2]['Cl']
            return self.interpolate(angle_of_attack, cl1, cl2, angle1, angle2)

    def cd(self, angle_of_attack): # μέθοδος που θα μας επιστρέψει την τιμή του συντελεστή αντίστασης (ή οπισθέλκουσας) Cd
        angle_of_attack = round(angle_of_attack, 3) # Στρογγυλοποίηση στα 3 δεκαδικά ψηφία
        if angle_of_attack in self.data:
            return self.data[angle_of_attack]['Cd']
        else:
            # Βρες τις δύο πλησιέστερες γωνίες προσβολής
            angle1, angle2 = self.get_nearest_angles(angle_of_attack)
            cd1 = self.data[angle1]['Cd']
            cd2 = self.data[angle2]['Cd']
            return self.interpolate(angle_of_attack, cd1, cd2, angle1, angle2)

    def cm(self, angle_of_attack): # μέθοδος που θα μας επιστρέψει την τιμή του συντελεστή ροπής Cm
        angle_of_attack = round(angle_of_attack, 3) # Στρογγυλοποίηση στα 3 δεκαδικά ψηφία
        if angle_of_attack in self.data:
            return self.data[angle_of_attack]['Cm']
        else:
            # Βρες τις δύο πλησιέστερες γωνίες προσβολής
            angle1, angle2 = self.get_nearest_angles(angle_of_attack)
            cm1 = self.data[angle1]['Cm']
            cm2 = self.data[angle2]['Cm']
            return self.interpolate(angle_of_attack, cm1, cm2, angle1, angle2)

    def tc_ratio(self, angle_of_attack): # μέθοδος που θα μας επιστρέψει το λόγο t/c ratio
        angle_of_attack = round(angle_of_attack, 3) # Στρογγυλοποίηση στα 3 δεκαδικά ψηφία
        if angle_of_attack in self.data:
            return self.data[angle_of_attack]['t/c ratio']
        else:
            # Βρες τις δύο πλησιέστερες γωνίες προσβολής
            angle1, angle2 = self.get_nearest_angles(angle_of_attack)
            tc1 = self.data[angle1]['t/c ratio']
            tc2 = self.data[angle2]['t/c ratio']
            return self.interpolate(angle_of_attack, tc1, tc2, angle1, angle2)

# Φόρτωση το αρχείο CSV
naca4415 = Naca4415('csv_data_file.csv')

# Παράδειγμα κλήσης για Cl, Cd, Cm και t/c ratio για μια δεδομένη γωνία προσβολής
a_o_a = 27

cl_value = naca4415.cl(a_o_a)
cd_value = naca4415.cd(a_o_a)
cm_value = naca4415.cm(a_o_a)
tc_value = naca4415.tc_ratio(a_o_a)

# print(f"Για τη γωνία προσβολής {a_o_a}° οι συντελεστές είναι:")
# print(f"Cl: {cl_value}, Cd: {cd_value}, Cm: {cm_value}, t/c ratio: {tc_value}")
