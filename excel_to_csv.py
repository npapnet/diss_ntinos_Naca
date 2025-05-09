import pandas as pd

# Φρόντισε να εγκαταστήσεις "pandas" και "openpyxl" πριν τρέξεις το script:
#   pip install pandas openpyxl

excel_file = "FFA-W3-CL_CD_CM_long.xlsx"  # Το όνομα του Excel αρχείου σου
output_csv = "FFA-W3-CL_CD_CM_long.csv"   # Το CSV αρχείο που θέλεις να δημιουργηθεί


df = pd.read_excel(excel_file, sheet_name=0)

# Αποθήκευση σε CSV με delimiter=";"
df.to_csv(output_csv, index=False, sep=';', encoding='utf-8')

print(f"Το αρχείο {output_csv} δημιουργήθηκε με επιτυχία!")
