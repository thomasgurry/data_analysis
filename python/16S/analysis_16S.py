"""

OVERVIEW:  

Miscellanous analysis tools and scripts for 16S genomic data analysis.

"""

import pandas as pd
import numpy as np

def collapse_OTU_table(OTU_table, taxonomic_level):
    # Collapse OTU table to a given taxonomic level
    # Input: OTU table dataframe with OTUs in a column called 'OTU_ID'
    # Output: OTU table collapsed to a given taxonomic level
    taxonomic_levels = {'phylum': {'RDP_string': 'p__', 'position': 1},
                        'class': {'RDP_string': 'c__', 'position': 2},
                        'order': {'RDP_string': 'o__', 'position': 3},
                        'family': {'RDP_string': 'f__', 'position': 4},
                        'genus': {'RDP_string': 'g__', 'position': 5}}

    # Loop through each OTU and add it to a given collapsed taxon
    taxon_dict = {}
    RDP_position = taxonomic_levels[taxonomic_level]['position']
    for i in OTU_table.index:
        OTU = OTU_table['OTU_ID'].ix[i]
        taxon = ';'.join(OTU.split(';')[:RDP_position+1])
        if taxon not in taxon_dict:
            taxon_dict[taxon] = []
        taxon_dict[taxon].append(i)

    # Create new taxonomic level table
    column_names = OTU_table.columns.tolist()
    sampleIDs = column_names[1:]
    collapsed_OTU_table = pd.DataFrame(columns = [taxonomic_level] + sampleIDs)
    for taxon in taxon_dict:
        OTU_indices = taxon_dict[taxon]
        taxon_abundances = np.array(len(sampleIDs)*[0.0])
        for i in OTU_indices:
            taxon_abundances = taxon_abundances + OTU_table[sampleIDs].ix[i]
        tmp = [taxon] + list(taxon_abundances)
        newentry = pd.DataFrame([tmp], columns = [taxonomic_level] + sampleIDs)
        collapsed_OTU_table = pd.concat([collapsed_OTU_table, newentry], axis=0, ignore_index=True)
    
    return collapsed_OTU_table
