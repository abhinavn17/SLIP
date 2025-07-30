import sys
import os 
import subprocess

fits_files = sys.argv[1:]

with open("sofia.par", "r") as f:
    par_lines = f.readlines()

for fits_file in fits_files:
    # Update the input.data line
    new_lines = [
        f"input.data = {fits_file}\n" if line.startswith("input.data") else line
        for line in par_lines
    ]
    # Write to a temporary par file
    with open("sofia.par", "w") as temp_f:
        temp_f.writelines(new_lines)
    # Run SOFIA (replace 'sofia' with the actual command if needed)
    subprocess.run(["sofia", "sofia.par"])