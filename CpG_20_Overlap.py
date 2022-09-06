#! /usr/bin/python

import pandas as pd
import numpy as np
from itertools import repeat
import matplotlib.pyplot as plt
import seaborn as sns

# calculate overlap of CpG sites with methylation frequency > 20% between replicates

# methylation sites in all samples with frquency >20%

# Nanopolish
df = {}
barcode_list = ['01','02','03','04','05','06','08','09','10','11']
for barcode_number in barcode_list:
    df[barcode_number] = pd.read_csv('barcode{}_methylation_frequency.tsv'.format(barcode_number),sep='\t')
    df[barcode_number] = df[barcode_number][df[barcode_number]['methylated_frequency']>=0.20]
    df[barcode_number]['new'] = mgc1(df[barcode_number]['chromosome'],df[barcode_number]['start']).rename('new')


CpG_sites_overlap = {
    "BR1": {i for i in df['08']['new'].values.tolist()},
    "BR2": {i for i in df['10']['new'].values.tolist()},
    "BR3": {i for i in df['01']['new'].values.tolist()},
    "BR4": {i for i in df['03']['new'].values.tolist()},
    "BR5": {i for i in df['05']['new'].values.tolist()}
}

# plot Venn diagram
from venn import venn
venn(CpG_sites_overlap)
plt.savefig('venn_30hpi.png',dpi=500)

# Guppy
# compare guppy and nanopolish results
barcode_list = ['01','02','03','04','05','06','08','09','10','11']
df = {}
for barcode_number in barcode_list:
    df[barcode_number] = pd.DataFrame()
    df[barcode_number]['modified_base_percentage'] = df[barcode_number][df[barcode_number]['modified_base_percentage']>=20]['modified_base_percentage']
    df[barcode_number]['start'] = df[barcode_number][df[barcode_number]['modified_base_percentage']>=20]['start']
    df[barcode_number]['chromosome'] = df[barcode_number][df[barcode_number]['modified_base_percentage']>=20]['chromosome']
    df[barcode_number]['new'] = mgc1(df[barcode_number][df[barcode_number]['modified_base_percentage']>20]['chromosome'],
                                                                  df[barcode_number][df[barcode_number]['modified_base_percentage']>20]['start']).rename('new')

CpG_sites_overlap = {
    "BR1": {i for i in df['08']['new'].values.tolist()},
    "BR2": {i for i in df['10']['new'].values.tolist()},
    "BR3": {i for i in df['01']['new'].values.tolist()},
    "BR4": {i for i in df['03']['new'].values.tolist()},
    "BR5": {i for i in df['05']['new'].values.tolist()}
}


from venn import venn
venn(CpG_sites_overlap)
plt.savefig('overlap_30hpi_guppy.png',dpi=500)
