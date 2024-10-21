#Κλάση για την αεροτομή Naca4415
#Reynolds = 200.000

import csv

class Naca4415:
    def __init__(self, csv_data_file):
        self.data = {}
        self.load_data(csv_data_file)

    def load_data(self, csv_data_file):
        with open(csv_data_file, mode = 'r') as file:
            reader = csv.reader(file)
            next(reader) # Αγνοούμε την πρώτη γραμμή του αρχείου που είναι οι επικεφαλίδες από κάθε στήλη
            for row in reader:
                angle_of_attack = round(float(row[0]), 3)
                cl_value = float(row[1])
                cd_value = float(row[2])
                self.data[angle_of_attack] = {'Cl': cl_value, 'Cd': cd_value}
                
        self.sorted_angles = sorted(self.data.keys()) # ταξινόμηση

    def interpolate(self, angle_of_attack, val1, val2, angle1, angle2):
        return val1 + (val2 - val1) * (angle_of_attack - angle1) / (angle2 - angle1) # από τη σχέση της γραμμικής παρεμβολής y = y1 + [(x-x1)/(x2-x1)]*(y2-y1)
    #val1 και val2 οι τιμές για τους συντελεστές cl και cd αντίστοιχα 
    
    def get_nearest_angles(self, angle_of_attack): # μέθοδος που βρίσκει τις δύο πλησιέστερες γωνίες προσβολής 
        for i in range(len(self.sorted_angles) - 1):
            if self.sorted_angles[i] <= angle_of_attack <= self.sorted_angles[i + 1]:
                return self.sorted_angles[i], self.sorted_angles[i + 1]
        raise ValueError(f'Η γωνία προσβολής {angle_of_attack} είναι εκτός του εύρους.')

    def cl(self, angle_of_attack): # μέθοδος που θα μας επιστρέψει την τιμή του συντελεστή άνωσης cl
        angle_of_attack = round(angle_of_attack, 3) # Στρογγυλοποίηση στα 3 δεκαδικά ψηφία
        if angle_of_attack in self.data:
            return self.data[angle_of_attack]['Cl']
        else:
            # Βρες τις δύο πλησιέστερες γωνίες προσβολής
            angle1, angle2 = self.get_nearest_angles(angle_of_attack)
            # Γραμμική παρεμβολή για τον Cl
            cl1 = self.data[angle1]['Cl']
            cl2 = self.data[angle2]['Cl']
            return self.interpolate(angle_of_attack, cl1, cl2, angle1, angle2)

    def cd(self, angle_of_attack): # μέθοδος που θα μας επιστρέψει την τιμή του συντελεστή αντίστασης (ή οπισθέλκουσας) cd
        angle_of_attack = round(angle_of_attack, 3) # Στρογγυλοποίηση στα 3 δεκαδικά ψηφία 
        if angle_of_attack in self.data:
            return self.data[angle_of_attack]['Cd']
        else:
            # Βρες τις δύο πλησιέστερες γωνίες προσβολής
            angle1, angle2 = self.get_nearest_angles(angle_of_attack)
            # Γραμμική παρεμβολή για τον Cd
            cd1 = self.data[angle1]['Cd']
            cd2 = self.data[angle2]['Cd']
            return self.interpolate(angle_of_attack, cd1, cd2, angle1, angle2)

naca4415 = Naca4415('csv_data_file.csv')

cl_value = naca4415.cl(8.55)
cd_value = naca4415.cd(8.55)
print(f"Για την συγκεκριμένη γωνία προσβολής οι συντελεστές είναι: Cl: {cl_value}, Cd: {cd_value}")

    
        