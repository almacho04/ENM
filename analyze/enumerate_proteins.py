import os
import csv

folder_path = os.path.expanduser("~/Desktop/soldier_girl/everything/processed_MD/universes/")
output_csv_path = os.path.expanduser("~/Desktop/soldier_girl/everything/protein_names.csv")

pickle_files = [file for file in os.listdir(folder_path) if file.endswith(".pkl")]



protein_names = {}

for idx, file_name in enumerate(pickle_files):
    abstracted_name = f"protein_{str(idx+1).zfill(3)}"
    protein_names[abstracted_name] = file_name

sizes = [os.path.getsize(os.path.join(folder_path, file)) for file in pickle_files]

with open(output_csv_path, "w") as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow([ "file_name","abstracted_name", "size"])
    for abstracted_name, file_name, size in zip(protein_names.keys(), protein_names.values(), sizes):
        ## only name without extension
        writer.writerow([file_name[:-4], abstracted_name, size])
        

print(f"CSV file '{output_csv_path}' created with {len(protein_names)} entries.")
