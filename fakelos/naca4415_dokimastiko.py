import csv

class Naca4415:
    def __init__(self, csv_data_file):
        self.data = {}
        self.tc_values = set()
        self.load_data(csv_data_file)

    def load_data(self, csv_data_file):
        with open(csv_data_file, mode='r') as file:
            reader = csv.reader(file, delimiter=';')
            next(reader)  # Αγνοούμε την πρώτη γραμμή του αρχείου
            for row in reader:
                angle_of_attack = round(float(row[0]), 3)
                cl_value = float(row[1])
                cd_value = float(row[2])
                cm_value = float(row[3])
                tc_ratio = float(row[4])
                self.tc_values.add(tc_ratio)
                
                if angle_of_attack not in self.data:
                    self.data[angle_of_attack] = {}
                self.data[angle_of_attack][tc_ratio] = {'Cl': cl_value, 'Cd': cd_value, 'Cm': cm_value}

        self.sorted_angles = sorted(self.data.keys())
        self.sorted_tc_values = sorted(self.tc_values)

    def interpolate(self, x, x1, x2, y1, y2):
        return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
    
    def get_nearest(self, value, sorted_list):
        for i in range(len(sorted_list) - 1):
            if sorted_list[i] <= value <= sorted_list[i + 1]:
                return sorted_list[i], sorted_list[i + 1]
        raise ValueError(f'Η τιμή {value} είναι εκτός εύρους.')
    
    def get_coefficients(self, angle_of_attack, tc_ratio):
        angle_of_attack = round(angle_of_attack, 3)
        tc_ratio = round(tc_ratio, 3)
        
        if angle_of_attack in self.data and tc_ratio in self.data[angle_of_attack]:
            return self.data[angle_of_attack][tc_ratio]
        
        angle1, angle2 = self.get_nearest(angle_of_attack, self.sorted_angles)
        tc1, tc2 = self.get_nearest(tc_ratio, self.sorted_tc_values)
        
        values = {}
        for coef in ['Cl', 'Cd', 'Cm']:
            f11 = self.data[angle1][tc1][coef] if tc1 in self.data[angle1] else self.data[angle1][tc2][coef]
            f12 = self.data[angle1][tc2][coef] if tc2 in self.data[angle1] else self.data[angle1][tc1][coef]
            f21 = self.data[angle2][tc1][coef] if tc1 in self.data[angle2] else self.data[angle2][tc2][coef]
            f22 = self.data[angle2][tc2][coef] if tc2 in self.data[angle2] else self.data[angle2][tc1][coef]
            
            f1 = self.interpolate(tc_ratio, tc1, tc2, f11, f12)
            f2 = self.interpolate(tc_ratio, tc1, tc2, f21, f22)
            values[coef] = self.interpolate(angle_of_attack, angle1, angle2, f1, f2)
        
        return values

# Παράδειγμα χρήσης
naca4415 = Naca4415('csv_data_file.csv')
angle_of_attack = 30
chosen_tc_ratio = 24.1  # Τυχαία επιλεγμένο t/c ratio από τα διαθέσιμα δεδομένα
coefficients = naca4415.get_coefficients(angle_of_attack, chosen_tc_ratio)

print(f"Για γωνία προσβολής {angle_of_attack}° και t/c ratio {chosen_tc_ratio}:")
print(f"Cl: {coefficients['Cl']}, Cd: {coefficients['Cd']}, Cm: {coefficients['Cm']}")
