# Identifying 5mC (5-Methylcytosine) and 5hmC (5-Hydroxymethylcytosine) in *Plasmodium falciparum*

 [Nanopolish](https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html) and [Guppy](https://labs.epi2me.io/gm24385-5mc/) pipelines were used to identify cytosine methylation during intra-erythrocytic stage of *Plasmodium falciparum*. 

The steps are shown below:
![alt text](https://github.com/Hann-Zhang/Plasmodium_falciparum_ONT_5mC_calling/blob/main/figures/Pipelines.png)

The examples bash scipts for requred softwares in the pipelines are in the the folder 'bash_scripts'. 

The code files are for parsing Nanopolish output (methylation_calls.csv) and Guppy output (csv files from modbam2bed), calculating methylation levels and making plots. 



