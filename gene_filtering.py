#!/bin/sh

#  gene_test.py
#  
#
#  Created by Ruby Mattingley on 14/10/2024.
#  

import pandas as pd
import numpy as np
import os.path
import sys
import argparse

def genes_of_interest_jvar(file_loc, genes_of_int_file, new_file_name, sample_file):
  """
  Function to format a jvar file by sample ID, and only including variants in
  specific genes of interest.

  file_loc: Filtered file created by format_jvar_w_filtering python script.
  genes_of_int_file: tsv file containing genes of interest, with each gene being on a
    new line, as seen in the ANN[0].GENE column of the original jvar file
  sample_file: a tsv file containing all samples from the study, with each sample ID
    being on a new line. If this is provided, all samples will be included in resulting
    file, even if they don't contain a variant in the genes of interest.
  new_file_name: string for the name of resulting csv file

  Returns a csv file with the first column being the sample IDs and columns 2
  onwards are genes of interest, with the specific HGVSc listed, separated by ;
  """
  
  # Getting filtered file and genes of interest
  data=pd.read_csv(file_loc)
  genes_df = pd.read_csv(genes_of_int_file, sep='\t', header=None)
  genes_of_int = list(genes_df[0].values)

  # Creating empty dataframe with Sample ID column and genes of interest
  df = pd.DataFrame(columns = ['Sample ID'] + genes_of_int)

  for gene_ind in range(len(genes_of_int)):
    # Make smaller dataframe for specific gene
    gene_df = data.loc[data['ANN[0].GENE']==genes_of_int[gene_ind]]

    # For each variant
    for ind in range(len(gene_df)):
      hgvs_c = gene_df.iloc[ind]['ANN[0].HGVS_C']
      sample = gene_df.iloc[ind]['SAMPLE']

      #if already in dataframe
      if sample in df['Sample ID'].values:
        sam_ind = df.loc[df['Sample ID']==sample].index[0]
        # If variants exist for current gene:
        if df.at[sam_ind, genes_of_int[gene_ind]] != '':
          df.at[sam_ind, genes_of_int[gene_ind]] += ' ; ' + hgvs_c
        # No variants exist for current gene:
        else:
          df.at[sam_ind, genes_of_int[gene_ind]] += hgvs_c

      # not in dataframe
      else:
        vals = [''] * (len(genes_of_int) + 1)
        vals[0] = sample
        vals[gene_ind+1] = hgvs_c
        df.loc[len(df.index)] = vals
        
  # If a sample file is provided:
  if sample_file != None:
      # Get sample values, and samples with variants in dataframe
      samples_df = pd.read_csv(sample_file, sep='\t', header=None)
      sample_list = list(samples_df[0].values)
      variant_samples = df['Sample ID'].values
      
      # For each sample:
      for samp in sample_list:
        # If sample not in the dataframe:
        if samp not in variant_samples:
          # Add empty row with that sample ID
          vals = [''] * (len(genes_of_int) + 1)
          vals[0] = sample
          df.loc[len(df.index)] = vals

  # Sort values by Sample ID and create new file
  df.sort_values(by='Sample ID', ignore_index=True, inplace=True)
  df.to_csv(new_file_name + '.csv', index=False)


def gene_counts(file_loc, new_file_name, genes_of_int_file, all_genes):
  """
  Function to produce a new file that counts the number of variants per
  gene of interest, with the option to include all genes.

  file_loc: Filtered file created by format_jvar_w_filtering python script.
  new_file_name: string for the name of resulting csv file
  genes_of_int_file: tsv file containing genes of interest, with each gene being on a
    new line, as seen in the ANN[0].GENE column of the original jvar file
  all_genes: option to include all genes in the resulting file rather than just
    the genes of interest.

  Returns a csv file with the first column being the genes with the second column
  being the total number of variants for the gene.
  """
  # Create dataframe for filtered file
  df = pd.read_csv(file_loc)
  
  # If all_genes option is false:
  if all_genes==False:
    # Create a dictionary with all genes in filtered file and create file.
    gene_dict = {'GENES':list(df['ANN[0].GENE'].value_counts().index), 'COUNT':list(df['ANN[0].GENE'].value_counts().values)}
    count_df = pd.DataFrame(gene_dict)
    count_df.to_csv(new_file_name + '.csv', index=False)
    
  # If all_genes is true:
  else:
    # Get genes of interest information, and create reduced dataframe with only genes of interest
    genes_df = pd.read_csv(genes_of_int_file, sep='\t', header=None)
    genes_of_int = list(genes_df[0].values)
    mask = df['ANN[0].GENE'].isin(genes_of_int)
    filtered_df = df[mask]
    
    # Create a dictionary with gene counts and create file
    gene_dict = {'GENES':list(filtered_df['ANN[0].GENE'].value_counts().index), 'COUNT':list(filtered_df['ANN[0].GENE'].value_counts().values)}
    count_df = pd.DataFrame(gene_dict)
    count_df.to_csv(new_file_name + '.csv', index=False)


# Arguments

parser = argparse.ArgumentParser(prog="gene_filtering", description="Program takes a csv file that has already been formatted/filtered by format_jvar_w_filtering and crates a new csv that is organised by SAMPLE, with other columns being genes of interest (provided by user), with each entry being the HGVS_c values. Additional arguments to include every sample (even if they don't have a variant in the genes of interest) and to output a second csv file for total counts per gene (each sample and unique HGVS_c).")

parser.add_argument("filtered_file", type=str, help="Location of the filtered csv file (created from format_jvar_w_filtering script), can just be file name if the script is in the same working directory")
parser.add_argument("genes_of_int_file", type=str, help="A tsv file that contains a list of genes we want to filter for, with each gene on a new line")
parser.add_argument("new_file_by_gene", type=str, help="Name of resulting csv file to be saved, can use whole file path if desired, otherwise new file will be saved in current working directory")
parser.add_argument("-s", "--sample_file", type=str, help="A tsv file containing all samples, if you want resulting csv file to contain every sample even if they don't have a variant in the genes of interest", default=None)
parser.add_argument("-g", "--new_file_gene_counts", type=str, help="Name of file that has the counts of all genes (number of times that gene has a unique variant in a sample), with the option to filter based on genes_of_int_file", default=None)
parser.add_argument("-a", "--all_genes_in_count", help="Name of file that has the counts of all genes (number of times that gene has a unique variant in a sample), with the option to filter based on genes_of_int_file", action="store_true")

args = parser.parse_args()


# Testing

if os.path.exists(args.filtered_file)==False:
    sys.exit("Input filtered file path does not exist")
if (args.filtered_file.endswith('.csv')==False):
    sys.exit("Input filtered file must be a csv file")
df = pd.read_csv(args.filtered_file)
col_names = ['SAMPLE','CHROM','POS','REF','ALT','ANN[0].GENE','ANN[0].HGVS_C','ANN[0].HGVS_P','DP','AF']
for col in col_names:
  if col not in df.columns:
    sys.exit("File does not contain required column: "+ str(col))

if os.path.exists(args.genes_of_int_file)==False:
    sys.exit("Genes of interest file does not exist")
if (args.genes_of_int_file.endswith('.tsv')==False):
    sys.exit("Genes of interest file must be a tsv file")

if args.sample_file != None:
    if os.path.exists(args.sample_file)==False:
        sys.exit("Sample file does not exist")
    if (args.sample_file.endswith('.tsv')==False):
        sys.exit("Sample file must be a tsv file")

if args.all_genes_in_count == True:
    if args.new_file_gene_counts == None:
        sys.exit("Cannot have all_genes in count active if there isn't a gene count file created")
    

def main():
    genes_of_interest_jvar(args.filtered_file, args.genes_of_int_file, args.new_file_by_gene, args.sample_file)
    if args.new_file_gene_counts != None:
        gene_counts(args.filtered_file, args.new_file_gene_counts, args.genes_of_int_file, args.all_genes_in_count)

if __name__ == "__main__":
    main()
