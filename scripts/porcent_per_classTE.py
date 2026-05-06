#!/usr/bin/env python

'''                             Description
        Quantification of the total number of base pairs annotated for 
        each major TE subclass in each genome. 
        Generates the file 'subclass_TEs-all.tsv' to subsequently create 
        the treemap with the script 'treemap_TEclassi.R' included in this 
        repository.
'''

import pandas as pd
import subprocess
import os


# 1. Configure file paths and species names
# Add or adjust paths according to your directory structure and BED output files from EarlGrey
files_map = {
    "Rice": "/EarlGrey/rice_summaryFiles/rice.filteredRepeats.bed",
    "Sorghum": "/EarlGrey/sorghum_summaryFiles/sorghum.filteredRepeats.bed",
    "Maize": "/EarlGrey/maize_summaryFiles/maize.filteredRepeats.bed",
    "CC-01-1940": "/EarlGrey/CC-01-1940/sugarcane_summaryFiles/sugarcane.filteredRepeats.bed",
    "R570": "/EarlGrey/R570/sugarcane_summaryFiles/sugarcane.filteredRepeats.bed",
    "KK3": "/EarlGrey/KK3/sugarcane_summaryFiles/sugarcane.filteredRepeats.bed",
    "SP80-3280": "/EarlGrey/SP80-3280/sugarcane_summaryFiles/sugarcane.filteredRepeats.bed",
    "AP85-441": "/EarlGrey/AP85-441/sugarcane_summaryFiles/sugarcane.filteredRepeats.bed",
    "Np-X": "/EarlGrey/Np-X/sugarcane_summaryFiles/sugarcane.filteredRepeats.bed",
    "GXS87-16": "/EarlGrey/GXS87-16/sugarcane_summaryFiles/sugarcane.filteredRepeats.bed",
    "LA-Purple": "/EarlGrey/LA-Purple/sugarcane_summaryFiles/sugarcane.filteredRepeats.bed"
}

# 2. TE subclasses of interest (as they appear at the beginning of column 4 in the BED file)
te_classes = ["DNA", "RC", "LTR", "LINE", "SINE", "Unknown"]

def process_bed(file_path, te_class):
    """Filter, sort, merge, and sum sequence lengths."""
    if not os.path.exists(file_path):
        print(f"Warning: File not found: {file_path}")
        return 0

    try:
        # Filtering logic for Unknown and Unclassified
        if te_class == "Unknown":
            filter_cmd = "awk '$4 ~ /^Unknown/ || $4 ~ /^Unclassified/' " + file_path
        else:
            filter_cmd = "awk '$4 ~ /^" + te_class + "/' " + file_path

        # Pipeline: Filter -> Sort -> Bedtools Merge
        cmd = "{} | sort -k1,1 -k2,2n | bedtools merge -i stdin".format(filter_cmd)
        
        # Compatible version with Python < 3.7 (using stdout=subprocess.PIPE)
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()
        
        if stderr:
            # If there is an error in the bedtools command, notify the user
            pass 

        total_bases = 0
        for line in stdout.strip().split('\n'):
            if line:
                cols = line.split('\t')
                # Sum (End - Start)
                total_bases += int(cols[2]) - int(cols[1]) + 1
        
        return total_bases
    except Exception as e:
        print("Error processing {} in {}: {}".format(te_class, file_path, e))
        return 0

# 3. Main execution
data = {"Class_TE": te_classes}

for species, path in files_map.items():
    print("Processing: {}...".format(species))
    counts = []
    for te_class in te_classes:
        total = process_bed(path, te_class)
        counts.append(total)
    data[species] = counts

# 4. Create DataFrame and save
df = pd.DataFrame(data)

# Organize columns to maintain dictionary order
cols_order = ["Subclasse_TE"] + list(files_map.keys())
df = df[cols_order]

df.to_csv("subclass_TEs-all.tsv", sep="\t", index=False)
print("\nDone! Table saved as 'subclass_TEs-all.tsv'.")
