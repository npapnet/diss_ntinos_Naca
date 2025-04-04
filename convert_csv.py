import csv
import math

def convert_csv_degs_to_radians_naca(input_file, output_file):
    """
    Διαβάζει το αρχείο CSV όπου η 1η στήλη είναι γωνίες σε μοίρες,
    με διαχωριστικό κόμμα ','.
    Μετατρέπει την 1η στήλη σε ακτίνια
    και γράφει το αποτέλεσμα σε output_file.
    """
    with open(input_file, mode='r', newline='', encoding='utf-8') as fin, \
         open(output_file, mode='w', newline='', encoding='utf-8') as fout:
        
        reader = csv.reader(fin, delimiter=',')    # Naca => comma
        writer = csv.writer(fout, delimiter=',')
        
        # Διαβάζουμε την 1η γραμμή (header). Το γράφουμε όπως είναι.
        header = next(reader)
        writer.writerow(header)
        
        for row in reader:
            # row[0] = angle σε μοίρες (string)
            angle_deg = float(row[0])
            angle_rad = math.radians(angle_deg)
            # Δημιουργούμε νέα γραμμή, 1η στήλη = angle σε ακτίνια
            new_row = [f"{angle_rad:.6f}"] + row[1:]
            writer.writerow(new_row)

def convert_csv_degs_to_radians_dtu(input_file, output_file):
    """
    Διαβάζει το αρχείο CSV όπου η 1η στήλη είναι γωνίες σε μοίρες,
    με διαχωριστικό ';'.
    Μετατρέπει την 1η στήλη σε ακτίνια
    και γράφει το αποτέλεσμα σε output_file.
    """
    with open(input_file, mode='r', newline='', encoding='utf-8') as fin, \
         open(output_file, mode='w', newline='', encoding='utf-8') as fout:
        
        reader = csv.reader(fin, delimiter=';')    # DTU => semicolon
        writer = csv.writer(fout, delimiter=';')
        
        header = next(reader)
        writer.writerow(header)
        
        for row in reader:
            angle_deg = float(row[0])
            angle_rad = math.radians(angle_deg)
            new_row = [f"{angle_rad:.6f}"] + row[1:]
            writer.writerow(new_row)

if __name__ == "__main__":
    # Μετατροπή για Naca
    convert_csv_degs_to_radians_naca(
        "csv_data_file_Naca.csv",
        "csv_data_file_Naca_in_radians.csv"
    )
    
    # Μετατροπή για DTU
    convert_csv_degs_to_radians_dtu(
        "csv_data_file_DTU.csv",
        "csv_data_file_DTU_in_radians.csv"
    )

    print("Η μετατροπή ολοκληρώθηκε με επιτυχία!")
