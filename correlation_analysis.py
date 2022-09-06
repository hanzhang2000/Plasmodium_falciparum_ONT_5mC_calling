#! /usr/bin/python

import pandas as pd
import numpy as np
from itertools import repeat
import matplotlib.pyplot as plt
import seaborn as sns

# calculate correlation between methylation output of replicates

# combine chromosome and bp position to coordinates
def mgc1(chromo, stt):
    return chromo.astype(str) + ':' + stt.astype(str)

# Nanopolish
# load data 
df = {}
barcode_list = ['01','02','03','04','05','06','08','09','10','11']
for barcode_number in barcode_list:
    df[barcode_number] = pd.read_csv('barcode{}_methylation_frequency.tsv'.format(barcode_number),sep='\t')
    df[barcode_number]['coordinates'] = mgc1(df[barcode_number]['chromosome'],df[barcode_number]['start']).rename('new')
    df[barcode_number] = df[barcode_number].drop(columns=['chromosome','start','end','num_motifs_in_group','called_sites_methylated','called_sites','group_sequence'])
    df[barcode_number] = df[barcode_number].rename(columns={'methylated_frequency':'barcode{}'.format(barcode_number)})

merge1 = pd.merge(df['01'],df['03'],on='coordinates')
merge2 = pd.merge(merge1,df['05'],on='coordinates')
merge3 = pd.merge(merge2,df['08'],on='coordinates')
final_methylation_list = pd.merge(merge3,df['10'],on='coordinates')

common_30hrs_methylation_list = pd.DataFrame()
j = ['BR1','BR2','BR3','BR4','BR5']
for count,i in enumerate(['08','10','01','03','05']):
    common_30hrs_methylation_list[j[count]] = final_methylation_list['barcode{}'.format(i)]

# merge1 = pd.merge(met_barcode02,met_barcode04,on='coordinates')
# merge2 = pd.merge(merge1,met_barcode06,on='coordinates')
# merge3 = pd.merge(merge2,met_barcode09,on='coordinates')
# final_methylation_list = pd.merge(merge3,met_barcode11,on='coordinates')

# common_36hrs_methylation_list = pd.DataFrame()
# j = ['BR1','BR2','BR3','BR4','BR5']
# for count,i in enumerate(['09','11','02','04','06']):
#     common_36hrs_methylation_list[j[count]] = final_methylation_list['barcode{}'.format(i)]

# correlation between different variables
corr = common_methylation_list.corr()

# Set up the matplotlib plot configuration
f, ax = plt.subplots(figsize=(12, 10))

# # Configure a custom diverging colormap
cmap = sns.diverging_palette(230, 20, as_cmap=True)
sns.set(font_scale = 2)

# Draw the heatmap and save the figure
sns_plot = sns.heatmap(corr, annot=True, cmap=cmap,vmin=0.5, vmax=1, annot_kws={'fontsize': 12, 'color':'black'})
sns_plot.set_title('30hpi Samples (Nanopolish)',fontsize = 20)
fig = sns_plot.get_figure()
fig.savefig("30hrs_correlation.png",dpi=500)

# Guppy

# import guppy methylation results iteratvely
df = {}
barcode_list = ['01','02','03','04','05','06','08','09','10','11']
for barcode_number in barcode_list:
    df[barcode_number] = pd.DataFrame()
    f1 = open("barcode{}/guppy.5mC.cpg.bam".format(barcode_number),'r')
    chromosome_list = list()
    start_position_list = list()
    end_position_list = list()
    modified_base_percentage_list = list()
    unmodified_base_count_list = list()
    modified_base_count_list = list()

    lines = f1.readlines()
    # strand = ['+','-']

    for line in lines:

        #split the line into a list by whitespace
        splitLine = line.rstrip().split()

        if splitLine[5] == '+':

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
    
    df[barcode_number]['chromosome'] = chromosome_list
    df[barcode_number] = df[barcode_number].astype({'chromosome':'int'})
    df[barcode_number]['start'] = start_position_list
    df[barcode_number] = df[barcode_number].astype({'start':'int'})
    df[barcode_number]['end'] = end_position_list
    df[barcode_number]['modified_base_percentage'] = modified_base_percentage_list
    df[barcode_number]['unmodified_base_count'] = unmodified_base_count_list
    df[barcode_number]['modified_base_count'] = modified_base_count_list
    df[barcode_number]['coordinates'] = mgc1(df[barcode_number]['chromosome'],df[barcode_number]['start']).rename('new')
    df[barcode_number] = df[barcode_number].rename(columns={'modified_base_percentage':'barcode{}_methylated_frequency'.format(barcode_number)})

merge1 = pd.merge(df['01'],df['03'],on='coordinates')
merge2 = pd.merge(merge1,df['05'],on='coordinates')
merge3 = pd.merge(merge2,df['08'],on='coordinates')
final_methylation_list = pd.merge(merge3,df['10'],on='coordinates')

common_30hrs_methylation_list = pd.DataFrame()
j = ['BR1','BR2','BR3','BR4','BR5']
for count,i in enumerate(['08','10','01','03','05']):
    common_30hrs_methylation_list[j[count]] = final_methylation_list['barcode{}_methylated_frequency'.format(i)]

corr = common_methylation_list.corr()

# Set up the matplotlib plot configuration
f, ax = plt.subplots(figsize=(12, 10))

# # Generate a mask for upper traingle
# mask = np.triu(np.ones_like(corr, dtype=bool))

# # Configure a custom diverging colormap
cmap = sns.diverging_palette(230, 20, as_cmap=True)
sns.set(font_scale = 2)

# Draw the heatmap and save the figure
sns_plot = sns.heatmap(corr, annot=True, cmap=cmap,vmin=0.5, vmax=1, annot_kws={'fontsize': 12, 'color':'black'})
sns_plot.set_title('30hpi Samples (Guppy)',fontsize = 20)
fig = sns_plot.get_figure()
fig.savefig("36hrs_correlation_guppy.png",dpi=500)