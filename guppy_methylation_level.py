#! /usr/bin/python

import pandas as pd
import numpy as np
from itertools import repeat
import matplotlib.pyplot as plt

barcode_list = ['01','02','03','04','05','06','08','09','10','11']
modification_type_list = ['5mC','5hmC']
sites_type_list = ['cpg','chh','chg']
strand_list = ['+','-']

# summarise methylation informarion for all types of modification in different sites and on different strands
for strand in strand_list:
    for modification in modification_type_list: 
        for site in sites_type_list:
            methylation_data_barcodes = pd.DataFrame()

            for barcode in barcode_list:
                #load file 
                f1 = open("barcode{}/guppy.{}.{}.bam".format(barcode,modification,site),'r')

                #creates container 
                chromosome_list = list()
                start_position_list = list()
                end_position_list = list()
                modified_base_percentage_list = list()
                unmodified_base_count_list = list()
                modified_base_count_list = list()

                lines = f1.readlines()
                

                for line in lines:

                    #split the line into a list by whitespace
                    splitLine = line.rstrip().split()
                    
                    if splitLine[5] == strand:
                        
                        chromosome = float(splitLine[0])
                        start_position = float(splitLine[1])
                        end_position = float(splitLine[2])
                        modified_base_percentage = float(splitLine[10])
                        unmodified_base_count = float(splitLine[11])
                        modified_base_count = float(splitLine[12])

                        chromosome_list.append(chromosome)
                        start_position_list.append(start_position)
                        end_position_list.append(end_position)
                        modified_base_percentage_list.append(modified_base_percentage)
                        unmodified_base_count_list.append(unmodified_base_count)
                        modified_base_count_list.append(modified_base_count)

                f1.close()

                # number of CpG sites in the genome
                num_CpG = len(modified_base_percentage_list)

                # average value of methylation over CpG sites
                average_methylation_level = np.nanmean(modified_base_percentage_list)

                # number and percentage of sites methylated
                threshold = 0
                n_methylation_sites =  [i for i in modified_base_count_list if i > threshold]
                p_methylation_sites = len(n_methylation_sites)/len(modified_base_percentage_list)

                # total number of cytosines in Pfalciparum reference genome 
                n_cytosines = sum(x.count('C') for x in open('Plasmodium_falciparum.ASM276v2.fasta','r') if x[0] != '>')

                # percentage of methylated 5mC in CpG islands 
                p_5mC_CpG_overC = len(n_methylation_sites)/n_cytosines

                data_list = [average_methylation_level, len(n_methylation_sites), num_CpG, p_methylation_sites, p_5mC_CpG_overC]
                methylation_data_barcodes['barcode_{}'.format(barcode)] = data_list

                methylation_data_barcodes.index = ['average_5mC_CpG_level', 'number of CpG sites methylated', 'total number of CpG sites', 
                                   'proportion of CpG sites methylated', 'proportion of cytosines methylation (CpG)']

            methylation_data_barcodes.to_csv('methylation_data_{}_{}_{}.csv'.format(strand,modification,site))  
