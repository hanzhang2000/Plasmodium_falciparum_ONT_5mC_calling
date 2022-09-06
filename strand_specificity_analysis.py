#! /usr/bin/python

import pandas as pd
import numpy as np
from itertools import repeat
import matplotlib.pyplot as plt
import itertools
import seaborn as sns
from scipy.stats import mannwhitneyu

# Box Plots and Mann Whitney U tests between strands

barcode_list = ['01','02','03','04','05','06','08','09','10','11']
modification_type_list = ['5mC','5hmC']
sites_type_list = ['cpg']
strand_list = ['+','-']
j = ['BR1','BR2','BR3','BR4','BR5']

modified_base_percentage_all_list = []
sample_list = []

for count,barcode in enumerate(['08','10','01','03','05']):
    for strand in strand_list:
        #load file 
        f1 = open("barcode{}/guppy.5mC.cpg.bam".format(barcode),'r')

        #creates container 
        chromosome_list = list()
        start_position_list = list()
        end_position_list = list()
#         globals()[f"modified_base_percentage_list_{barcode}"] = list()
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
#                 globals()[f"modified_base_percentage_list_{barcode}"].append(modified_base_percentage)
                modified_base_percentage_list.append(modified_base_percentage)
                unmodified_base_count_list.append(unmodified_base_count)
                modified_base_count_list.append(modified_base_count)

        f1.close()

        modified_base_percentage_all_list = modified_base_percentage_all_list + modified_base_percentage_list
        it = itertools.repeat('{}_strand{}'.format(j[count],strand),len(modified_base_percentage_list))
        sample_list = sample_list + list(it)
        
df = pd.DataFrame()
df['methylation frequency'] = modified_base_percentage_all_list
df['sample_strand'] = sample_list

# Mann Whitney U test
barcode_list = ['01','02','03','04','05','06','08','09','10','11']
for barcode_number in barcode_list:
    U1, p = mannwhitneyu(df[df['sample_strand']=='barcode{}_strand+'.format(barcode_number)]['methylation frequency'].dropna(), 
                     df[df['sample_strand']=='barcode{}_strand-'.format(barcode_number)]['methylation frequency'].dropna(), method="auto")
    print(p)

# plot Box plots
from matplotlib import rcParams
df['methylation frequency'] = np.arcsinh(df['methylation frequency'])
sns.set(rc={'figure.figsize':(25,10)})
sns.set_style('whitegrid')
sns.set(font_scale = 2)
rcParams['axes.titlepad'] = 30

my_pal = {"BR1_strand+": "r", "BR1_strand-": "b","BR2_strand+": "r", "BR2_strand-": "b",
         "BR3_strand+": "r", "BR3_strand-": "b","BR4_strand+": "r", "BR4_strand-": "b",
         "BR5_strand+": "r", "BR5_strand-": "b"}
sns_plot = sns.boxplot(x=df['sample_strand'],y=df['methylation frequency'],width=0.3, palette=my_pal);
sns_plot.set_xlabel("Sample", fontsize = 25,labelpad=15)
sns_plot.set_ylabel("Methylation Frequency (arcsinh transformed)", fontsize = 25,labelpad=15)
sns_plot.set_title("The Mann-Whitney U test between + and - strands (30hpi) (5mC)", fontsize = 30)


x1, x2 = 0, 1   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
y, h, col = df['methylation frequency'].max() + 0.1, 0.1, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, "***", ha='center', va='bottom', color=col,fontsize = 20)

x3, x4 = 2, 3  # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
y, h, col = df['methylation frequency'].max() + 0.1, 0.1, 'k'
plt.plot([x3, x3, x4, x4], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x3+x4)*.5, y+h, "***", ha='center', va='bottom', color=col,fontsize = 20)

x5, x6 = 4, 5   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
y, h, col = df['methylation frequency'].max() + 0.1, 0.1, 'k'
plt.plot([x5, x5, x6, x6], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x5+x6)*.5, y+h, "***", ha='center', va='bottom', color=col,fontsize = 20)

x1, x2 = 6, 7   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
y, h, col = df['methylation frequency'].max() + 0.1, 0.1, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, "***", ha='center', va='bottom', color=col,fontsize = 20)

x1, x2 = 8, 9   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
y, h, col = df['methylation frequency'].max() + 0.1, 0.1, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, "***", ha='center', va='bottom', color=col,fontsize = 20)

plt.show()
fig = sns_plot.get_figure()
file_name = 'guppy_strand_Mann_Whitney_30hpi_5mC.png'
fig.savefig(file_name,dpi=500, bbox_inches='tight')
fig.clf()



