'''                             Description
        Separates cases in which the TE intersects the exon from those in which
        the TE is only nearby, and summarizes this information by genome and for
        the combined set of species in tables.
        Output files: "confidence_intervals_TEs-TOTALexons.antisense.tsv"
                      "confidence_intervals_TEs-TOTALexons.sense.tsv"
'''

import pandas as pd
import numpy as np
import glob
import math

# Z value for a 95% confidence level
Z = 1.96

# Search for all BED files
files = glob.glob("*-closest_TEs-TOTALexons.antisense-FILTRADO.bed")      # *-closest_TEs-TOTALexons.sense-FILTRADO.bed

# Global lists to store data from all combined species
all_overlaps = []
all_distances = []
results = []

def calculate_statistics(data, group, category):
    """Calculate mean, standard deviation, margin of error, and 95% CI."""
    n = len(data)
    if n < 2:  # Standard deviation requires at least 2 observations
        return None
    
    mean = np.mean(data)
    std_dev = np.std(data, ddof=1)  # ddof=1 for the correct sample standard deviation
    std_error = std_dev / math.sqrt(n)
    margin_error = Z * std_error
    
    lower_limit = mean - margin_error
    upper_limit = mean + margin_error
    
    return {
        'Group': group,
        'Category': category,
        'Sample_N': n,
        'Mean': round(mean, 2),
        'Standard_Deviation': round(std_dev, 2),
        'Margin_Error_95%': round(margin_error, 2),
        'Lower_Limit': round(lower_limit, 2),
        'Upper_Limit': round(upper_limit, 2)
    }

for file in files:
    cultivar = file.split('-')[0]
    df = pd.read_csv(file, sep=r'\s+', header=None)
    
    # ---------------------------------------------------------
    # 1. TEs that INTERSECT exons (last column == 0)
    #    Calculation: Extract the actual overlap distance
    # ---------------------------------------------------------
    df_intersect = df[df[14] == 0].copy()
    if not df_intersect.empty:
        # Overlap logic: min(ends) - max(starts)
        min_ends = df_intersect[[2, 10]].min(axis=1)
        max_starts = df_intersect[[1, 9]].max(axis=1)
        overlap = (min_ends - max_starts).values
        
        # Filter only valid overlaps (> 0)
        valid_overlap = overlap[overlap > 0]
        
        all_overlaps.extend(valid_overlap)
        stats = calculate_statistics(valid_overlap, cultivar, 'Overlap (Intersects)')
        if stats:
            results.append(stats)
        
    # ---------------------------------------------------------
    # 2. TEs that DO NOT INTERSECT exons (last column > 0)
    #    Calculation: Directly use the reported distance value
    # ---------------------------------------------------------
    df_not_intersect = df[df[14] > 0].copy()
    if not df_not_intersect.empty:
        # The distance in base pairs is already provided by bedtools closest in the last column
        distances = df_not_intersect[14].values
        
        all_distances.extend(distances)
        stats = calculate_statistics(distances, cultivar, 'Distance (Does Not Intersect)')
        if stats:
            results.append(stats)

# ---------------------------------------------------------
# 3. Calculate metrics for the complete set (i.e., all species)
# ---------------------------------------------------------
if all_overlaps:
    stats_all_overlaps = calculate_statistics(all_overlaps, 'COMPLETE_SET', 'Overlap (Intersects)')
    if stats_all_overlaps:
        results.append(stats_all_overlaps)

if all_distances:
    stats_all_distances = calculate_statistics(all_distances, 'COMPLETE_SET', 'Distance (Does Not Intersect)')
    if stats_all_distances:
        results.append(stats_all_distances)

# ---------------------------------------------------------
# 4. Generate and export the table
# ---------------------------------------------------------
if results:
    df_results = pd.DataFrame(results)
    
    # Sort so that COMPLETE_SET appears at the end and cultivars are grouped by category
    df_results = df_results.sort_values(by=['Category', 'Group'])
    
    print(df_results.to_string(index=False))
    
    output_name = "confidence_intervals_TEs-TOTALexons.antisense.tsv"     # "confidence_intervals_TEs-TOTALexons.sense.tsv"
    df_results.to_csv(output_name, sep='\t', index=False)
    print(f"\n[+] Analysis completed. Table saved as: {output_name}")
else:
    print("Not enough data were found to calculate statistics.")
