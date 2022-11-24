import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import dinuq as dn  # This package provides a range of metrics for quantifying dinucleotide in genetic sequences. (pip install dinuq)
from scipy.stats import f_oneway
import random
from scipy.stats import tukey_hsd # Este análisis no fue posible correrlo desde la terminal, para correrlo debe hacerse en la consola de python inactivar sys.argv y activar line 46 to 48
import scipy.stats as stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from random import choice

# Create random sequence with a length of 6156701 to compare with pathogen and No pathogen sequences

print(f'Create random sequence')
length = 6156701  # Define length of sequence

bases = ["A", "G", "C", "T"]  # supply list of nucleotides

# Generate the random sequence



# create empty sequence
random_sequence = ""

# randomly select base and add to sequence(repeat 100 times for length=6156701)
for i in range(length):
    base = choice(bases)
    random_sequence += base

# define function to create a fasta file with the random sequence

f = open("random_seq.fasta", "w")
f.write(f'>random \n{random_sequence}')

# print(sequence) # Just to verify sequence

# # Activate to run code from terminal and use: python iRutina_P4.py Pseudomonas_aeruginosa.fasta Pseudomonas_putida.fasta random_seq.fasta
# Pathogen = sys.argv[1]
# No_pathogen = sys.argv[2]
# Other = sys.argv[3]

# Activate to run code python console
Pathogen = 'Pseudomonas_aeruginosa.fasta'
No_pathogen = 'Pseudomonas_putida.fasta'
Other = 'random_seq.fasta'

##############################---Relative dinucleotide abundance---#####################################
# Calculate relative dinucleotide abundance(RDA) in Pseudomonas Aeruginosa

print(f'Creating list of dinucleotides')

dinucl_list = ['TpA', 'TpC', 'TpG', 'TpT', 'ApA',
               'ApC', 'ApG', 'ApT', 'CpA', 'CpC',
               'CpG', 'CpT', 'GpA', 'GpC', 'GpG', 'GpT']  # list of dinucleotides

print(f'Calculate RDAs from fasta file')

# RDA 1 - pseudomonas aeruginosa

rda_pathogen = dn.RDA(Pathogen, dinucl_list)
print(rda_pathogen)

RDA1_df = pd.DataFrame.from_dict(rda_pathogen).reset_index()  # Create data frame with RDAs values
RDA1_df.rename(columns={'index': 'Dinucleotide'}, inplace=True)
print(RDA1_df)

# RDA 2 - Psedomonas putida

rda_nopathogen = dn.RDA(No_pathogen, dinucl_list)
print(rda_nopathogen)

RDA2_df = pd.DataFrame.from_dict(rda_nopathogen).reset_index()  # Create data frame with RDAs values
RDA2_df.rename(columns={'index': 'Dinucleotide'}, inplace=True)
print(RDA2_df)

# RDA 3 - Random sequence

rda_other = dn.RDA(Other, dinucl_list)
print(rda_other)

RDA3_df = pd.DataFrame.from_dict(rda_other).reset_index()  # Create data frame with RDAs values
RDA3_df.rename(columns={'index': 'Dinucleotide'}, inplace=True)
print(RDA3_df)

################################--Data Frame ---########################################################

print(f'Joining RDAs tables')

# Join RDA's Tables
df_conc = RDA1_df.join(RDA2_df[list(rda_nopathogen.keys())[0]])  # Concatenate RDA1 and RDA2 data frame
print(df_conc)

df_conc_2 = df_conc.join(RDA3_df[list(rda_other.keys())[0]])  # Concatenate df_conc and RDA3 data frame
print(df_conc_2)

###############################################################################
print(f'editing concatening data frame')
# Add names to columns in the fd_conc_2
dfm = df_conc_2.melt('Dinucleotide', var_name='Accession', value_name='RDA')
print(dfm)
# dfm.to_csv('Results_RDA.csv') # Use if you need a table with RDAs values
####################################--GRAFICS--###########################################
# Generate plot with RDAs Results


print(f'plots')

# SCATTERPLOT
e = sns.scatterplot(data=dfm,
                    x="Dinucleotide",
                    y="RDA",
                    hue="Accession",
                    size="RDA",
                    sizes=(20, 200),
                    palette="blend:#7AB,#EDA",
                    legend="brief")
plt.show()
plt.savefig('scatterplot.png')

# Relplot
k = sns.relplot(data=dfm,
                x="Dinucleotide",
                y="RDA",
                hue="Accession",
                size="RDA",
                palette=("flare"),
                sizes=(10, 200))
plt.show()
plt.savefig('relplot.png')

# BARPLOT
f = sns.barplot(data=dfm,
                x="Dinucleotide",
                y="RDA",
                palette="Set2",
                hue="Accession")
plt.show()
plt.savefig('barplot.png')
# LINE-PLOT
g = sns.lineplot(data=dfm,
                 x="Dinucleotide",
                 y="RDA",
                 hue="Accession",
                 style="Accession",
                 palette="pastel",
                 markers=True,
                 dashes=False)
plt.show()
plt.savefig('lineplot.png')


#############################-ANÁLISIS ESTADÍSTICO-###############################

# ANOVA

print(f'calculating fvalue and pvalue using ANOVA')

fvalue, pvalue = stats.f_oneway(df_conc_2['NC_002516.2'], df_conc_2['NC_021505.1'], df_conc_2['random'])
print(fvalue, pvalue)

# BOX-PLOT
# Box-plot - each dinucleotide
h = sns.boxplot(x='Dinucleotide',
                y='RDA',
                palette="ch:s=.25,rot=-.25",
                data=dfm,
                )
plt.show()
plt.savefig('boxplot_dinucl.png')
# Box-plot - all dinucleotide
h = sns.boxplot(x='Accession',
                y='RDA',
                hue='Accession',
                data=dfm,
                palette='Set2')
plt.show()
plt.savefig('boxplot.png')

# TUKEY_HSD (Solo puede correrse desde la consola de python, no fue posible correrlo desde la terminal)

print(f'Tukey_HSD')

Tukey_HSD = tukey_hsd(df_conc_2['NC_002516.2'], df_conc_2['NC_021505.1'], df_conc_2['random'])

print(Tukey_HSD)

