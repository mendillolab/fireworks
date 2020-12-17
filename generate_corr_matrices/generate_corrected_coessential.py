import pandas as pd
import numpy as np
import os

df_wide = pd.read_csv('20q4_cleaned.csv', index_col=0) #wide format, cell lines index, gene symbols as columns
df1 = pd.read_csv('genome_annotations_biomart_updated.csv')
dgd = pd.read_csv('dgd_Hsa_all_v71.tsv', sep='\t', index_col='Name')

rank_threshold = 20 #on each side
dup_genes = dgd.index.tolist()

a = pd.DataFrame()
for gene in df_wide.columns.tolist():
	if gene in dup_genes:
		a[gene] = df_wide[gene] #just copy the original essentiality score if a duplicate gene
	else:
		chromosome = df1[df1['Gene'] == gene].Chromosome.values[0]
		df2 = df1[df1['Chromosome'] == chromosome].sort_values(by='End').reset_index()
		pos = df2[df2['Gene'] == gene].index[0]
		locus_genes = df2.loc[pos-rank_threshold:pos+rank_threshold, :].Gene.tolist()
		locus_genes.remove(gene)

		a[gene] = df_wide[gene] - (df_wide[[i for i in locus_genes if i in df_wide]].median(axis=1)/2)

b = a.corr()
b.to_hdf('corrected_coessentiality_matrix.h5', key='stage', mode='w')
print('Done')