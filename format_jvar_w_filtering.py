# filtering_ap_df

import pandas as pd
import numpy as np
import os.path
import sys
import argparse

def create_df_variants_dp_af(file_loc, additional_cols_file):
  """
  Function to create a dataframe organised by samples and including the columns:
  ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'ANN[0].GENE', 'ANN[0].HGVS_C', 'ANN[0].HGVS_P', 'DP', 'AF']

  file_loc: string location of jvar file to be formatted.
  additional_cols: a list containing strings of additional columns to include, default is no additional columns
  """
  # Creating the dataframe from the jvar file
  ext_file = pd.read_csv(file_loc,sep='\t')
  
  # Getting additional column names to include in resulting csv file
  if additional_cols_file != None:
    new_file = pd.read_csv(additional_cols_file, header=None)
    additional_cols = list(new_file[0].values)
  else:
    additional_cols = []


  dfl = []
  # This is where you could change the resulting columns that are in the resulting csv file
  cols = ['CHROM', 'POS', 'REF', 'ALT', 'ANN[0].GENE', 'ANN[0].HGVS_C', 'ANN[0].HGVS_P'] + additional_cols
  
  # For each variant
  for row_ind in range(len(ext_file)):
    
    # Get values from specified columns
    row = ext_file.iloc[row_ind]
    var_info = list(row[cols].values)

    # If variant is homozygous:
    if type(row['sample.HOM_VAR']) is str:
      # Extract the sample, AF and DP values
      samples = row['sample.HOM_VAR'].split(',')
      dps = row['DP.HOM_VAR'].split(',')
      afs = row['AF.HOM_VAR'].split(',')

      # Add to dataframe
      for i in range(len(samples)):
        dfl += [[samples[i]] + var_info + [int(dps[i])] + [float(afs[i])]]

    # If variant is heterozygous:
    if type(row['sample.HET']) is str:
      # Extract the sample, AF and DP values
      samples = row['sample.HET'].split(',')
      dps = row['DP.HET'].split(',')
      afs = row['AF.HET'].split(',')
      
      # Add to dataframe
      for i in range(len(samples)):
        dfl += [[samples[i]] + var_info + [int(dps[i])] + [float(afs[i])]]

  # Update column names and create dataframe
  col_names = ['SAMPLE'] + cols + ['DP', 'AF']
  df = pd.DataFrame(dfl, columns = col_names)
  return df


def filter_dp_af(df, new_file_name, dp_min, dp_max, af_min, af_max, var_filt):
  """
  Function to filter a file based on DP and AF thresholds.

  df: dataframe in format created by create_df_variants_dp_af function
  new_file_name: string for the name of resulting csv file
  dp_min: minimum depth
  dp_max: maximum depth
  af_min: minimum allele frequency
  af_max: maximum allele frequency
  """
  # Copy of dataframe from function above
  new_df = df.copy()
  
  # For each row and sample:
  for ind in range(len(new_df)):
    # If AF or DP value does not fall within the specified range, remove row from dataframe
    if (dp_min <= new_df['DP'][ind] <= dp_max and af_min <= new_df['AF'][ind] <= af_max) == False:
      new_df.drop(ind, inplace=True)

  # Create list for indices to drop based on frequency filtering
  ind_to_drop = []
  # Create dictionary based on counts of variants
  var_dict = dict(new_df[['ANN[0].GENE','ANN[0].HGVS_C']].value_counts())
  # For each variant:
  for key_val in var_dict.keys():
    # If the count of that variant is above or equal to the filtering frequency
    if var_dict[key_val] >= var_filt:
      # Get indices if too many variants
      ind_to_drop += list(new_df.loc[(df['ANN[0].GENE']==key_val[0]) & (new_df['ANN[0].HGVS_C']==key_val[1])].index)
  new_df.drop(ind_to_drop, inplace=True)
  
  # Create new file
  new_df.to_csv(new_file_name + '.csv', index=False)


# Arguments

parser = argparse.ArgumentParser(prog="format_jvar_w_filtering", description="Program takes a jvar file organised by variants, and outputs a csv file organised by samples with the following information: chromosome, position, ref, alt, gene, HGVS_c, HGVS_p, depth and allele frequency. Additional arguments to include additional columns from the jvar file or to filter based on AF or DP values.")

parser.add_argument("jvar_file_loc", help="Location of the jvar file to be formatted, can just be file name if script is in the same working directory", type=str)
parser.add_argument("-a", "--additional_columns", type=str, help="File containing additional columns to include in resulting output, with rows being additional columns to include", default=None)
parser.add_argument("new_file_name", type=str, help="Name of resulting csv file to be saved, can use whole file path if desired, otherwise new file will be saved in current working directory")
parser.add_argument("-minaf", "--min_allele_frequency", type=float, help="Minimum allele frequency to filter variants (value between 0 and 1)", default=0)
parser.add_argument("-maxaf", "--max_allele_frequency", type=float, help="Maximum allele frequency to filter variants (value between 0 and 1 and must be larger than min allele frequency)", default=1)
parser.add_argument("-mindp", "--min_depth", type=int, help="Minimum depth value to filter variants (must be an integer)", default=0)
parser.add_argument("-maxdp", "--max_depth", type=int, help="Maximum depth value to filter variants (must be an integer and must be larger than min depth)", default=1000000)
parser.add_argument("-s", "--study_freq_filter", type=int, help="Filter used to remove variants if the count of that variant is more or equal to frequency filter", default=100000000)

args = parser.parse_args()

# Testing

if os.path.exists(args.jvar_file_loc)==False:
    sys.exit("jvar file path does not exist")
    
if (args.jvar_file_loc.endswith('.tsv')==False):
    sys.exit("Input jvar file must be tsv file")

if args.additional_columns != None:
    if os.path.exists(args.additional_columns) == False:
        sys.exit("csv file for additional columns does not exist")
    if (args.additional_columns.endswith('.csv')==False):
        sys.exit("Additional columns file must be a csv file")

if (args.min_depth < 0):
    sys.exit("Minimum depth values must be an integer greater than 0")
    
if (args.max_depth <= args.min_depth):
    sys.exit("Maximum depth must be larger than the minimum depth.")

if (args.min_allele_frequency < 0) or (args.min_allele_frequency > 1):
    sys.exit("Minimum allele frequency must be a decimal between 0 and 1.")

if (args.max_allele_frequency < 0) or (args.max_allele_frequency > 1):
    sys.exit("Maximum allele frequency must be a decimal between 0 and 1.")
    
if (args.max_allele_frequency <= args.min_allele_frequency):
    sys.exit("Maximum allele frequency must be larger than the minimum allele frequency.")


def main():
    df = create_df_variants_dp_af(args.jvar_file_loc, args.additional_columns)
    filter_dp_af(df, args.new_file_name, args.min_depth, args.max_depth, args.min_allele_frequency, args.max_allele_frequency, args.study_freq_filter)

if __name__ == "__main__":
    main()
