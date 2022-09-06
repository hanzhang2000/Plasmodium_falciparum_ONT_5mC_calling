#! /usr/bin/python

import pandas as pd
import numpy as np
from itertools import repeat
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interp1d
import itertools

f1 = open("TSS.bed",'r')

chromosome_TSS_list = list()
start_position_TSS_list = list()
end_position_TSS_list = list()
strand_TSS_list = list()

lines = f1.readlines()
# strand = ['+','-']

for line in lines:

    #split the line into a list by whitespace
    splitLine = line.rstrip().split()


    chromosome = float(splitLine[0])
    start_position = float(splitLine[1])
    end_position = float(splitLine[2])
    strand = splitLine[3]

    chromosome_TSS_list.append(chromosome)
    start_position_TSS_list.append(start_position)
    end_position_TSS_list.append(end_position)
    strand_TSS_list.append(strand)

f1.close()

TSS = pd.DataFrame()
TSS['chromosome'] = chromosome_TSS_list
TSS['start'] = start_position_TSS_list
TSS['end'] = end_position_TSS_list
TSS['strand'] = strand_TSS_list

TSS_plus_strand = TSS[TSS['strand'] == '+']
TSS_minus_strand = TSS[TSS['strand'] == '-']

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

df2={}
barcode_list = ['01','02','03','04','05','06','08','09','10','11']
for barcode_number in barcode_list:
    
    df2[barcode_number] = pd.DataFrame()
    for i in range(len(TSS_minus_strand)-1):

        l = np.arange(int(TSS_minus_strand[i:i+1]['end'])-3000, int(TSS_minus_strand[i:i+1]['end'])+3000 , 50)

        range_value = []

        for j in range(len(l)-1):
            range_list = df[barcode_number][(df[barcode_number]['chromosome'] == int(TSS_minus_strand[i:i+1]['chromosome'])) 
                                                                     & (df[barcode_number]['start'] > l[j]) & (df[barcode_number]['start'] < l[j+1])]
            if len(range_list) != 0:
                proportion = len(range_list[range_list['modified_base_count'] != 0])/len(range_list)
            else:
                proportion = 1

            range_value.append(proportion)

        range_value = pd.DataFrame({'TSS_{}'.format(i): range_value})

        df2[barcode_number] = pd.concat([df2[barcode_number], range_value], axis=1)

# plot of methylation level near TSS

SMALL_SIZE = 20
MEDIUM_SIZE = 20
BIGGER_SIZE = 30

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

plt.figure()
plt.figure(figsize=(35, 15), dpi=80)
    
# barcode_list = ['01','03','05','08','10']
barcode_list = ['08','09','10','11','01','02','03','04','05','06']
legend_list = ['30hpi_BR1','36hpi_BR1','30hpi_BR2','36hpi_BR2','30hpi_BR3','36hpi_BR3','30hpi_BR4','36hpi_BR4',
               '30hpi_BR5','36hpi_BR5']
for count, barcode_number in enumerate(barcode_list):
    df2[barcode_number]['mean_rows'] = df2[barcode_number].mean(axis = 1)
    x = np.arange(-2975,2975,50)
    y = df2[barcode_number]['mean_rows']
#     plt.plot(x,y,'o-',label='Barcode{}'.format(barcode_number))
    
    if barcode_number in ['01','03','05','08','10']:
        plt.plot(x,y,'o-',color='blue', label='{}'.format(legend_list[count]))
    
    if barcode_number in ['02','04','06','09','11']:
        plt.plot(x,y,'o-',color='red', label='{}'.format(legend_list[count]))
    
plt.xlabel("Binned Distance to TSS (bp)",fontsize=40)
plt.ylabel("% of CpG sites methylated",fontsize=40)
plt.title("Negative Strand",fontsize=60)
# plt.legend()
file_name = 'TSS_negative_strand_BR.png'
plt.savefig(file_name,dpi=300, bbox_inches='tight')
plt.show()

