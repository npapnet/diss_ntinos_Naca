import csv

class Naca4415:
    def __init__(self, csv_data_file):
        self.data = {}
        self.load_data(csv_data_file)
        
    def load_data(self, csv_data_file):
        with open(csv_data_file, mode='r') as file:
            reader = csv.reader(file)
            next(reader) 
            for row in reader:
                angle_of_attack = float(row[0])
                cl_value = float(row[1])
                cd_value = float(row[2])
                self.data[angle_of_attack] = {'Cl': cl_value, 'Cd': cd_value}

    def cl(self, angle_of_attack): # μέθοδος που θα μας επιστρέψει την τιμή του συντελεστή άνωσης cl
        if angle_of_attack in self.data:
            return self.data[angle_of_attack]['Cl']
        else:
            raise ValueError("Η δοθείσα γωνία προσβολής δεν βρέθηκε")

    def cd(self, angle_of_attack): # μέθοδος που θα μας επιστρέψει την τιμή του συντελεστή αντίστασης (ή οπισθέλκουσας) cd
        if angle_of_attack in self.data:
            return self.data[angle_of_attack]['Cd']
        else:
            raise ValueError("Η δοθείσα γωνία προσβολής δεν βρέθηκε")

naca4415 = Naca4415('csv_data_file.csv')

cl_value = naca4415.cl(-5)
cd_value = naca4415.cd(-5)
print(f"Για την συγκεκριμένη γωνία προσβολής οι συντελεστές είναι: Cl: {cl_value}, Cd: {cd_value}")
        
        