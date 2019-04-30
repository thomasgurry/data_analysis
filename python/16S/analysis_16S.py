"""

OVERVIEW:  

Miscellanous analysis tools and scripts for 16S genomic data analysis.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns


# Plotting parameters
colorListHex = ["#F6CC70","#6DABD4","#79CBBF","#EC80AD","#B092C0","#B7CC93","#AB6A94","#c59d94","#f9d8ae","#C4F598","#E2347C", "#FF7F50", "#FF4500", "#FFA500", "#FFFACD", "#FFE4B5", "#228B22", "#8FBC8F", "#FFA07A", "#8B0000"];
colorsRGB = []
for hex in colorListHex:
    colorsRGB.append(colors.colorConverter.to_rgb(hex))
plt.rc("figure", figsize=(12, 8))
sns.set_style('white')



def normalise_OTU_table(OTU_table):
    # Normalise OTU table by dividing by total read count for each sample
    OTU_table_norm = OTU_table.copy()
    sampleIDs = OTU_table_norm.columns.tolist()[1:]
    for sid in sampleIDs:
        OTU_table_norm[sid] = np.divide(OTU_table_norm[sid], float(np.sum(OTU_table_norm[sid])))
    return OTU_table_norm


def collapse_OTU_table_RDP(OTU_table, taxonomic_level):
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


def collapse_OTU_table_SILVA(OTU_table, taxonomic_level):
    # Collapse OTU table to a given taxonomic level
    # Input: OTU table dataframe with OTUs in a column called 'OTU_ID'
    # Output: OTU table collapsed to a given taxonomic level
    taxonomic_levels = {'phylum': 1,
                        'class': 2,
                        'order': 3,
                        'family': 4,
                        'genus': 5,
                        'species': 6}

    # Loop through each OTU and add it to a given collapsed taxon
    taxon_dict = {}
    SILVA_string_position = taxonomic_levels[taxonomic_level]
    for i in OTU_table.index:
        OTU = OTU_table['OTU_ID'].ix[i]
        taxon = ';'.join(OTU.split(';')[:SILVA_string_position+1])
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


def stacked_barplot_RDP_table(OTU_table, taxonomic_level, sampleIDs=0, sample_labels=0):
    # Returns a stacked barplot of the taxa in an OTU table (at a specified taxonomic level) 
    # for each sample ID, which can be provided in a specific order.  If none are provided, all 
    # sample IDs in the OTU table will be plotted in the order they appear and labeled as in the table.
    
    # Collapse taxonomies
    otu_table_collapsed = collapse_OTU_table_RDP(OTU_table, taxonomic_level)
    fig_size = (12, 8)

    if sampleIDs == 0:
        sampleIDs = otu_table_collapsed.columns.tolist()[1:]
        sample_labels = sampleIDs
    
    # If phylum, display all, else display top 20
    if taxonomic_level == 'phylum':
        taxa = otu_table_collapsed[taxonomic_level]
    else: 
        all_taxa = np.array(otu_table_collapsed[taxonomic_level].tolist())
        taxa_median_abundances = []
        for i in range(len(all_taxa)):
            taxa_median_abundances.append(np.median(otu_table_collapsed[sampleIDs].ix[i]))
        sorted_inds = np.argsort(taxa_median_abundances)[::-1]
        otu_table_collapsed = otu_table_collapsed.ix[sorted_inds[:20]].copy()
        taxa = all_taxa[sorted_inds[:20]]

    print(taxa)

    stacked_barplot_dict = {sid: otu_table_collapsed[sid].tolist() for sid in sampleIDs}


    stacked_barplot_dict['taxa'] = taxa
    f, ax = plt.subplots(figsize=fig_size)
    inds = np.arange(len(sampleIDs))
    last_barplot_tuple = ()
    for sid in sampleIDs:
        last_barplot_tuple = last_barplot_tuple + (0.0, )

    for i in range(len(taxa)):
        taxon = taxa[i]
        barplot_tuple = ()
        for sid in sampleIDs:
            barplot_tuple = barplot_tuple + (stacked_barplot_dict[sid][i], )
        ax.bar(inds, barplot_tuple, bottom=last_barplot_tuple, color=colorsRGB[i%len(colorsRGB)])
        last_barplot_tuple = [sum(x) for x in zip(last_barplot_tuple, barplot_tuple)]

    taxa_names = [taxon.split(';')[-1] for taxon in taxa]
    ax.legend(taxa_names, bbox_to_anchor=(1.35, 1.0))    
    ax.set_xticks([a+0.5 for a in range(len(sample_labels))])
    ax.set_xticklabels(sample_labels, fontsize=12, rotation='vertical')    
    ax.set_ylim([0, 1.0])
    ax.set_ylabel('Relative abundance', fontsize=18)
    return ax


def stacked_barplot_SILVA_table(OTU_table, taxonomic_level, sampleIDs=0, sample_labels=0):
    # Returns a stacked barplot of the taxa in an OTU table (at a specified taxonomic level) 
    # for each sample ID, which can be provided in a specific order.  If none are provided, all 
    # sample IDs in the OTU table will be plotted in the order they appear.
    
    # Collapse taxonomies
    otu_table_collapsed = collapse_OTU_table_SILVA(OTU_table, taxonomic_level)
    fig_size = (12, 8)

    if sampleIDs == 0:
        sampleIDs = otu_table_collapsed.columns.tolist()[1:]
        sample_labels = sampleIDs

    # If phylum, display all, else display top 20
    if taxonomic_level == 'phylum':
        taxa = otu_table_collapsed[taxonomic_level]
    else: 
        all_taxa = np.array(otu_table_collapsed[taxonomic_level].tolist())
        taxa_median_abundances = []
        for i in range(len(all_taxa)):
            taxa_median_abundances.append(np.median(otu_table_collapsed[sampleIDs].ix[i]))
        sorted_inds = np.argsort(taxa_median_abundances)[::-1]
        taxa = all_taxa[sorted_inds[:20]]
        otu_table_collapsed = otu_table_collapsed.ix[sorted_inds[:20]].copy()

    stacked_barplot_dict = {sid: otu_table_collapsed[sid].tolist() for sid in sampleIDs}


    stacked_barplot_dict['taxa'] = taxa
    f, ax = plt.subplots(figsize=fig_size)
    inds = np.arange(len(sampleIDs))
    last_barplot_tuple = ()
    for sid in sampleIDs:
        last_barplot_tuple = last_barplot_tuple + (0.0, )

    for i in range(len(taxa)):
        taxon = taxa[i]
        barplot_tuple = ()
        for sid in sampleIDs:
            barplot_tuple = barplot_tuple + (stacked_barplot_dict[sid][i], )
        ax.bar(inds, barplot_tuple, bottom=last_barplot_tuple, color=colorsRGB[i%len(colorsRGB)])
        last_barplot_tuple = [sum(x) for x in zip(last_barplot_tuple, barplot_tuple)]

    taxa_names = [taxon.split(';')[-1] for taxon in taxa]
    ax.legend(taxa_names, bbox_to_anchor=(1.35, 1.0))    
    ax.set_xticks([a+0.5 for a in range(len(sample_labels))])
    ax.set_xticklabels(sample_labels, fontsize=12, rotation='vertical')    
    ax.set_ylim([0, 1.0])
    ax.set_ylabel('Relative abundance', fontsize=18)
    return ax
