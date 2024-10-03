#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 15:39:04 2024

@author: f99731hc
"""
import pandas as pd
import subprocess 
import glob, os 
import concurrent.futures
import time
__name__ == '__main__'



inter_sample_results = []
for i in range(0, 20):
    

    blast_wildcard = "".join(["your/filepath/barcode_",str(i),"/*.clean"])
    blast_files = glob.glob(blast_wildcard)
    
    
    barcode_id = blast_files[0].split("/")[-1].split(".")[0]
    
    blast_summary = []
    
    
    for blast_filename in blast_files:
    
        fasta_folder = "your/filepath/cluster/"
        raw_filename = blast_filename.split("/")[-1].split(".")
        barcode_id = raw_filename[0]
        cluster_id = raw_filename[1]
        # Input filenames
        fasta_filename = "".join([fasta_folder, barcode_id, "/", barcode_id, ".", cluster_id])
        # Read input file
        blast_data = pd.read_csv(blast_filename, delimiter='\t', header=None)
        blast_data = blast_data.sort_values(by=[2], ascending=False)
        # saccver qcovs bitscore nident mismatch gaps pident evalue
        # Genus level ID using species name
        blast_data[8] = [x.split(" ")[0] for x in blast_data[0]]
        # Identify maximum bit score
        max_bit_score = blast_data.iloc[0,2]
        
        # Check whether the read can be classified or not
        class_level = "Unknown fungal classification"
        within_range_species = set(blast_data[0][blast_data[2] > max_bit_score - 40])
        within_range_genus = set(blast_data[8][blast_data[2] > max_bit_score - 40])
        
        # Calculate the top and second species
        top_species = blast_data.iloc[0,0]
        top_genus = blast_data.iloc[0,8]
        species_level = "Unknown fungal classification"
        genus_level = "Unknown fungal classification"
        # Is this at species level
        if len(within_range_species) == 1:
            class_level = "Species"
            second_species = blast_data.iloc[0,0]
            second_genus = top_genus
            bit_difference = 0
            species_level = top_species
            genus_level = top_genus
        else:
            reduced_blast_data=blast_data[blast_data[0] != top_species]
            second_species = reduced_blast_data.iloc[0,0]
            second_genus = reduced_blast_data.iloc[0,8]
            bit_difference = max_bit_score - reduced_blast_data.iloc[0, 2]
            # Is this at genus level?
            if (len(within_range_genus) == 1):
                class_level = "Genus"
                reduced_blast_data=blast_data[blast_data[0] != top_species]
                second_species = reduced_blast_data.iloc[0,0]
                second_genus = reduced_blast_data.iloc[0,8]
                bit_difference = max_bit_score - reduced_blast_data.iloc[0, 2]
                genus_level = second_genus
            else:
                reduced_blast_data=blast_data[blast_data[8] != top_genus]
                second_species = reduced_blast_data.iloc[0,0]
                second_genus = reduced_blast_data.iloc[0,8]
                bit_difference = max_bit_score - reduced_blast_data.iloc[0, 2]
    
        read_count = subprocess.run(['grep', '-c', '>', fasta_filename], 
                                    capture_output=True, text=True).stdout.rstrip('\n')
        read_count = int(read_count)
        # Store result
        blast_summary += [[barcode_id, cluster_id, top_species, top_genus, second_species, second_genus, species_level, genus_level, bit_difference, read_count, class_level]]
        # Top species, second species, bit difference, read count
    
        
    calculated_taxa = pd.DataFrame(blast_summary)
    calculated_taxa.columns =['BarcodeID', 'ClusterID', 'TopSpecies', 'TopGenus', 'SecondSpecies', 'SecondGenus', 'SpeciesLevel', 'GenusLevel', 'BitDiff', 'ReadCount', 'ClassLevel']
    
    
    output = "".join(["your/filepath/data/taxa_calc/", barcode_id, ".csv"])
    calculated_taxa.to_csv(output, header = True, index = False)
