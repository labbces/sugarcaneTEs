#!/usr/bin/env python3

'''                             Description
        Script used for the integration step of the reclassification results 
        of the TEdistill library 'distilledTE.flTE.fa'.
        The classification files generated independently by TEsorter, DeepTE, 
        and RepeatClassifier are required.
        Generates a tabulated report of the final reclassification decisions: 
        'summary_decisions.tsv', available in this repository.
'''

import sys

def get_truncated(full_string):
    """
    Returns the normalized name up to the first slash in the classification.
    Example: 'TE_00009237#LINE/L1' -> 'TE_00009237#LINE'
    """
    if full_string == "NA" or "#" not in full_string:
        return full_string
    
    parts = full_string.split("#")
    id_part = parts[0]
    class_part = parts[1].split("/")[0]
    
    return f"{id_part}#{class_part}"


def main():
    # Input files
    file_sorterlib = "TEsorter_TEdistill_lib/Lib_TEdistll.TEsorter.cls.lib"
    file_sorterdom = "TEsorter_TEdistill_lib/Lib_TEdistll.TEsorter.dom.tsv"
    file_deepte = "DeepTE_TEdistill_lib/opt_DeepTE.txt"
    file_classified = "RepeatClassifier_TEdistill_lib/distilledTE.flTE.fa.classified"
    file_output = "summary_decisions.tsv"

    sorter_ids = set()
    deepte_map = {}
    classified_map = {}

    # 1. Load TEsorter IDs
    try:
        with open(file_sorterdom, 'r') as f:
            for line in f:
                parts = line.split('\t')[0].split('|')
                id_part = parts[0].split('#')[0]
                sorter_ids.add(id_part)
    except FileNotFoundError: pass

    # 2. Load DeepTE
    try:
        with open(file_deepte, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2: continue
                id_p = parts[0].split('#')[0]
                cls = parts[1].replace("ClassI_", "").replace("ClassII_", "").replace("_", "/", 1)
                deepte_map[id_p] = f"{id_p}#{cls}"
    except FileNotFoundError: pass

    # 3. Load classified FASTA
    try:
        with open(file_classified, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    full_name = line.strip().lstrip("> ").split()[0]
                    id_p = full_name.split("#")[0]
                    classified_map[id_p] = full_name
    except FileNotFoundError: pass

    # 4. Process main file and apply rules
    try:
        with open(file_sorterlib, 'r') as f_in, open(file_output, 'w') as f_out:
            for line in f_in:
                if not line.startswith("> "): continue
                
                col1 = line.strip().lstrip("> ").split()[0]
                id_seq = col1.split("#")[0]
                
                col2 = "YES" if id_seq in sorter_ids else "NO"
                col3 = deepte_map.get(id_seq, "NA")
                col4 = classified_map.get(id_seq, "NA")

                # Truncated names for comparison
                t1 = get_truncated(col1)
                t3 = get_truncated(col3)
                t4 = get_truncated(col4)

                # --- Compute Agreement (column 5) and Decision (column 6) ---
                agreement = 0
                decision = ""

                # Rule 1: All three are equal
                if t1 == t3 == t4 and t1 != "NA":
                    agreement = 3
                    decision = col1 if col2 == "YES" else col4
                
                # Rule 2: Columns 1 and 3 are equal
                elif t1 == t3 and t1 != "NA":
                    agreement = 2
                    decision = col1
                
                # Rule 3: Columns 1 and 4 are equal
                elif t1 == t4 and t1 != "NA":
                    agreement = 2
                    decision = col1 if col2 == "YES" else col4
                
                # Rule 4: Columns 3 and 4 are equal
                elif t3 == t4 and t3 != "NA":
                    agreement = 2
                    decision = col4
                
                # Rule 5: All different (new criterion)
                else:
                    if col2 == "YES":
                        agreement = 1
                        decision = col1
                    else:
                        agreement = 0
                        decision = f"{id_seq}#Unknown"

                # Write final line with TAB separator
                f_out.write(f"{col1}\t{col2}\t{col3}\t{col4}\t{agreement}\t{decision}\n")

        print(f"Report successfully generated: {file_output}")
    except FileNotFoundError:
        print("Error: Files not found.")

if __name__ == "__main__":
    main()
