#! /usr/bin/python

import pandas as pd
import numpy as np
from itertools import repeat
import matplotlib.pyplot as plt

barcode_list = ['01','02','03','04','05','06','08','09','10','11']

# calculate the average 5mC level in CpG sites and number of CpG sites that have at least one positive methylation calls

methylation_data_barcodes = pd.DataFrame()

for barcode in barcode_list:

    # import data
    methylation_data = pd.read_csv('barcode{}_methylation_frequency.tsv'.format(barcode),sep='\t')
    
    # total number of CpG islands 
    total_CpG_count = np.sum(methylation_data['num_motifs_in_group'])

    methylation_frequency_list = []
    for i in range(len(methylation_data)):
        methylation_frequency_list.extend(repeat(methylation_data['methylated_frequency'][i],
                                          methylation_data['num_motifs_in_group'][i]))
    
    # average methylation frequency at all CpG islands
    average_methylation_level = np.average(methylation_frequency_list)

    # number/proportion of CpG islands with non-zero methylation call
    threshold = 0
    n_methylation =  [i for i in methylation_frequency_list if i > threshold]
    p_methylation = len(n_methylation)/total_CpG_count

    # total number of cytosines in Pfalciparum reference genome 
    n_cytosines = sum(x.count('C') for x in open('Plasmodium_falciparum.ASM276v2.fasta','r') if x[0] != '>')

    # percentage of methylated 5mC in CpG islands 
    p_5mC_CpG_overC = len(n_methylation)/n_cytosines

    data_list = [average_methylation_level, len(n_methylation), total_CpG_count, p_methylation, p_5mC_CpG_overC]
    
    # save the data into file
    methylation_data_barcodes['barcode_{}'.format(barcode)] = data_list

methylation_data_barcodes.index = ['average_5mC_CpG_level', 'number of CpG sites methylated', 'total number of CpG sites', 
                                   'proportion of CpG sites methylated', 'proportion of cytosines methylation (CpG)']
methylation_data_barcodes.to_csv('methylation_data_barcodes.csv')  

    





