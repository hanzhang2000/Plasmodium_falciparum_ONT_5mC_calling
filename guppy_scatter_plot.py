#! /usr/bin/python

import pandas as pd
import numpy as np
from itertools import repeat
import matplotlib.pyplot as plt
import seaborn as sns
import itertools

# load guppy data
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

        # change the condition here for data in either positive or negative strand
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

# scatter plot of region in chromosome 7
SMALL_SIZE = 20
MEDIUM_SIZE = 20
BIGGER_SIZE = 30

sns.set_style('darkgrid')
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=40)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=50)  # fontsize of the figure title

plt.figure()
plt.figure(figsize=(35, 15), dpi=80)

barcode_list = ['08','09','10','11','01','02','03','04','05','06']
legend_list = ['30hpi_BR1','36hpi_BR1','30hpi_BR2','36hpi_BR2','30hpi_BR3','36hpi_BR3','30hpi_BR4','36hpi_BR4',
               '30hpi_BR5','36hpi_BR5']
for count, barcode_number in enumerate(barcode_list):
    data =  df[barcode_number][ ( df[barcode_number]['chromosome']==7) & ( df[barcode_number]['start']>=500000) & (df[barcode_number]['end']<=600000)]
    if barcode_number in ['01','03','05','08','10']:
        plt.plot(data['start'],data['modified_base_percentage'],'o',label='{}'.format(legend_list[count]),color='blue')
    
    if barcode_number in ['02','04','06','09','11']:
        plt.plot(data['start'],data['modified_base_percentage'],'o',label='{}'.format(legend_list[count]),color='green')
    

plt.xlabel("Chr7:500000-600000",labelpad=15,fontsize=40)
plt.ylabel("modified base percentage(%)",labelpad=15,fontsize=40)
plt.ylim(0,100)
plt.title("Positive Strand",fontsize=60)
plt.legend()
file_name = 'chr7:500000_600000_positive_strand.png'
plt.savefig(file_name,dpi=500, bbox_inches='tight')
plt.show()