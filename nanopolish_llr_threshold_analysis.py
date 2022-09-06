#! /usr/bin/python

import pandas as pd
import numpy as np
from itertools import repeat
import seaborn as sns
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

# Analyse Nanopolish output with different stringency - different log-likelihood ratio threshold

# define functions
# replace first line of BED files
def replace_first_line( src_filename, target_filename, replacement_line):
    f = open(src_filename)
    first_line, remainder = f.readline(), f.read()
    t = open(target_filename,"w")
    t.write(replacement_line + "\n")
    t.write(remainder)
    t.close()

# combine chromosome and bp position to coordinates
def mgc1(chromo, stt):
    return chromo.astype(str) + ':' + stt.astype(str)

# def mgc2(chromo, stt,end):
#     return chromo.astype(str) + ':' + stt.astype(str) + '-' + end.astype(str)

barcode_number_list = ['01','02','03','04','05','06','08','09','10','11']

# set strand and llr_threshold
for strand in ['minus']:

    for llr_threshold in [10]:
        data_sum = []

        for barcode_number in barcode_number_list:
            
            # load data
            nanopolish_calls = pd.read_csv('/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_22072022/barcode{}/barcode{}_methylation_calls.tsv'.format(barcode_number,barcode_number),sep='\t')
            
            if strand == 'plus':
                nanopolish_calls_strand = nanopolish_calls[nanopolish_calls['strand']=='+']
            
            if strand == 'minus':
                nanopolish_calls_strand = nanopolish_calls[nanopolish_calls['strand']=='-']
            
            # filter reads based on the llr threshold
            nanopolish_calls_strand_llr = nanopolish_calls_strand[nanopolish_calls_strand['log_lik_ratio']>=llr_threshold]

            # create a column of coodinates for both original dataframe and filtered dataframe
            nanopolish_calls_strand['coordinates'] = mgc1(nanopolish_calls_strand['chromosome'],nanopolish_calls_strand['start']).rename('coordinates')
            nanopolish_calls_strand_llr['coordinates'] = mgc1(nanopolish_calls_strand_llr['chromosome'],
                                                                        nanopolish_calls_strand_llr['start']).rename('coordinates')
            
            # group methylation calls from the same CpG site
            nanopolish_calls_strand_agg = nanopolish_calls_strand.groupby('coordinates').agg(list)
            nanopolish_calls_strand_agg['coordinates'] = nanopolish_calls_strand_agg.index

            # summarise the number of calls and proportion of calls over llr threshold for each CpG site
            data = []
            for location in set(nanopolish_calls_strand_llr['coordinates']):
                calls_list = nanopolish_calls_strand_agg[nanopolish_calls_strand_agg['coordinates']==location]
                chromosome = calls_list['chromosome'].tolist()[0][0]
                start = calls_list['start'].tolist()[0][0]
                end = calls_list['start'].tolist()[0][0]
                count = len(nanopolish_calls_strand_llr[nanopolish_calls_strand_llr['coordinates']==location])
                count_sum = len(calls_list['chromosome'].tolist()[0])
                count_proportion = len(nanopolish_calls_strand_llr[nanopolish_calls_strand_llr['coordinates']==location])/count_sum
                data.append([chromosome, start, end, location, count,count_sum,count_proportion])

            llr_count = pd.DataFrame(data,columns = ['chromosome','start','end','coordinates','count','count_sum','proportion'])

            # make BED file (methylation map based on methylation information)

            # below 0.02 - light orange 255 229 204
            # 0.02-0.05 - 255 153 153
            # 0.05-0.15 - 255 120 120
            # 0.15-0.30 - 255 0 0
            # >0.30 - 204 0 0 
            RgB_list = []
            for i in llr_count['proportion']:
                if i <=0.02:
                    RgB_list.append('255,229,204')
                if 0.02 < i <= 0.05:
                    RgB_list.append('255,153,51')
                if 0.05 < i <= 0.15:
                    RgB_list.append('255,128,0')
                if 0.15 < i <= 0.30:
                    RgB_list.append('255,102,102')
                if 0.30 < i:
                    RgB_list.append('255,0,0')

            llr_count = llr_count.drop(columns = ['coordinates','count','count_sum','proportion'])

            llr_count['NA1'] = '.'
            llr_count['NA2'] = '0'
            llr_count['NA3'] = '.'
            llr_count['repeat1'] = llr_count['start']
            llr_count['repeat2'] = llr_count['end']
            llr_count['RgB'] = RgB_list

            llr_count.to_csv('/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_22072022/new_bed_files/llr{}_strand{}_barcode{}.bed'.format(llr_threshold,strand,barcode_number), sep='\t', index=False)

            replace_first_line('/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_22072022/new_bed_files/llr{}_strand{}_barcode{}.bed'.format(llr_threshold,strand,barcode_number),'/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_22072022/new_bed_files/llr{}_strand{}_barcode{}.bed'.format(llr_threshold,strand,barcode_number), "track name=llr{}_strand{}_barcode{} itemRgb=On".format(llr_threshold,strand,barcode_number))
    
   

