# Copyright (C) 2021 Sebastian Hupfauf

import argparse
import sys
import time
import pandas as pd
from datetime import date

start = time.time()

##Data import

#Setup of argparse
parser = argparse.ArgumentParser(description='This tool can be used to filter a virus mutation table.', add_help=False, epilog='Thank you for using this filtering tool!\n\n(C) 2021\nSebastian Hupfauf')
requiredNamed = parser.add_argument_group('Required arguments')
optionalNamed = parser.add_argument_group('Optional arguments')
requiredNamed.add_argument('-d', '--data', action="store", type=str, required=True, help='Name of the input data file.')
requiredNamed.add_argument('-m', '--metadata', action="store", type=str, required=True, help='Name of the input metadata file.')
optionalNamed.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Show this help message and exit.')
optionalNamed.add_argument('-v', '--version', action='version', version='%(prog)s 1.0', help="Show program's version number and exit.")
optionalNamed.add_argument('-o', '--output_dir', action='store', type=str, help="Ourput directory where all result files are stored. If no output directory is provided, files are stored in the current working directory!")
optionalNamed.add_argument('-n', '--name', action='store', type=str, help="Name of the analysis, it will be used for labelling all result files. If no name is provided, data files will be labelled with 'DATA.txt' and 'META.txt'!")
optionalNamed.add_argument('-s', '--start', action="store", default='1', type=int, help='First epidemiological week to be included in the analysis. [Default: 1]')
optionalNamed.add_argument('-e', '--end', action="store", default='-1', type=int, help='Last epidemiological week to be included in the analysis. Use -1 for the last week. [Default: -1]')
optionalNamed.add_argument('-r', '--region', action="store", default='all', type=str, help='Regions to be included in the analysis. In case of more than one region, seperate each region by a comma sign. Use "all" to include all regions. [Default: all]')
optionalNamed.add_argument('-l', '--location', action="store", default='all', type=str, help='Locations to be included in the analysis. In case of more than one location, seperate each location by a comma sign. Use "all" to include all locations. [Default: all]')
optionalNamed.add_argument('-rd', '--read_depth', action="store", default='-1', type=int, help='Threshold for read depth, lower values are interpreted as 0. Use -1 if you do not want to apply a threshold for read depth. [Default: -1]')
optionalNamed.add_argument('-af', '--allele_frequency', action="store", default='-1', type=float, help='Threshold for allele frequency, lower values are interpreted as 0. Use -1 if you do not want to apply a threshold for allele frequency. Input must be a number between 0 and 1, e.g. 0.05 or 0.01. [Default: -1]')
optionalNamed.add_argument('-ms', '--min_samples', action="store", default='-1', type=int, help='Minimum number of samples in which a mutation must be present in order to be retained. Use -1 if you do not want to apply this filter. [Default: -1]')
optionalNamed.add_argument('-um', '--undetected_mutations', action="store", default='-1', type=float, help='Maximum percentage of undetected mutations in a sample to be retained. Use -1 if you do not want to apply this filter. Input must be a number between 0 and 1, e.g. 0.05 or 0.01. [Default: -1]')
optionalNamed.add_argument('-cv', '--coefficient_variant', action="store", default='-1', type=float, help='Maximum coefficient of variant of a sample to be retained. Use -1 if you do not want to apply this filter. Input must be a number between 0 and 1, e.g. 0.7 or 1.3. [Default: -1]')
args = parser.parse_args()
data = args.data
meta = args.metadata
#Creating dataframes
D = pd.read_csv(data, delimiter="\t", dtype=str)
M = pd.read_csv(meta, delimiter="\t", dtype=str, index_col=0)
#Creating an ID column and setting it as index for the DF
D["ID"] = D["REF"] + D["POS"] + D["ALT"]
D.set_index("ID", inplace = True)

###Cleaning of the metadata file

sample_count = len(M.index)
#Remove failed-tagged samples from the DF
M = M[M["failed"] != "1"]
failed = sample_count - len(M.index)
#Remove samples without an information on the location, sample date, and facility
M = M[M["sample_source_location_state"].notnull()]
M = M[M["sample_date"].notnull()]
M = M[M["LocationName_coronA"].notnull()]
missing_information = sample_count - failed - len(M.index)
#Output
print("\nINITIAL DATA CLEAN-UP:")
print("Original number of samples:\t\t", sample_count)
print("Removed because fail-tagged:\t\t", failed)
print("Removed because of missing information:\t", missing_information)
#print("Remaining samples:\t\t\t", len(M.index))
samples = list(M.index)

###Cleaning of the data file

#Removing unneccesary columns including the quality (PQ) to be more memory efficient
remaining_columns = [k for k in D.columns if "PQ" not in k]
D[remaining_columns].to_csv("datasub.txt", sep="\t")
#Removing the unneccesary columns from the DF:
remaining_columns = remaining_columns[18:]
#Removing samples that were dropped from the metadata file
remaining_columns = [k for k in remaining_columns if k.split("_")[0] + "_" + k.split("_")[1] in samples]
D = D[remaining_columns]

###Removing samples from metadata that are not included in the data

before = len(samples)
remaining_samples = []
for x in remaining_columns:
    if "AF" in x:
        remaining_samples.append(x.split("_")[0] + "_" + x.split("_")[1])
samples = [k for k in samples if k in remaining_samples]
samples_after_cleanup = len(samples)
M = M.loc[samples, :]
print("Removed because missing in the data:\t", before - samples_after_cleanup)
print("Remaining samples:\t\t\t", samples_after_cleanup)

###Making a subset based on epidemiological week

#Create metadata columns for epidemiological day and week
M["ep_day"] = M["sample_date"].apply(lambda x: (date(*[int(k) for k in x.split("-")]) - date(2020, 2, 26)).days)
M["ep_week"] = ((M["ep_day"] - 1) // 7) + 1
#Application
if args.end == -1:
    args.end = M["ep_week"].max()
if not (args.start <= 1 and args.end >= M["ep_week"].max()):
    beg = args.start
    ex = args.end
    M = M[(M["ep_week"] >= beg) & (M["ep_week"] <= ex)]
    for col in D.columns:
        if not (col.split("_")[0] + "_" + col.split("_")[1]) in M.index:
            del D[col]

###Making a subset based on regions

if args.region != "all":
    if "," in args.region:
        sel = args.region.split(",")
        M = M[M["sample_source_location_state"].isin(sel)]
    else:
        M = M[M["sample_source_location_state"] == args.region]
    for col in D.columns:
        if not (col.split("_")[0] + "_" + col.split("_")[1]) in M.index:
            del D[col]

###Making a subset based on locations

if args.location != "all":
    if "," in args.location:
        sel = args.location.split(",")
        M = M[M["LocationName_coronA"].isin(sel)]
    else:
        M = M[M["LocationName_coronA"] == args.location]
    for col in D.columns:
        if not (col.split("_")[0] + "_" + col.split("_")[1]) in M.index:
            del D[col]

###Report at the beginning of the analysis

print("\nBEGINNING OF THE ANALYSIS:")
print("Number of samples: ", int(len(D.columns) / 2))
print("Number of mutations: ", len(D.index))

###Remove missing values (".") and replace them with 0, and convert strings to float

D.replace('^\.','0', regex=True,inplace=True)
D = D.astype(float)

###Applying a threshold for the read depth

if args.read_depth != -1:
    thresh = args.read_depth
    for col in D.columns:
        if ":DP" in col:
            D[col] = D[col].apply(lambda x: 0.0 if x < thresh else x)
            D[col.split(":")[0] + ":AF"].loc[D.loc[D[col] == 0].index] = 0

###Applying a threshold for the allele frequency

if args.allele_frequency != -1:
    thresh = args.allele_frequency
    for col in D.columns:
        if ":AF" in col:
            D[col] = D[col].apply(lambda x: 0.0 if x < thresh else x)

###Delete columns containing only 0 and adjust metadata

D = D.loc[:, (D != 0).any(axis=0)]
af = [x.strip(":AF") for x in D.columns if ":AF" in x]
dp = [x.strip(":DP") for x in D.columns if ":DP" in x]
overlap = sorted(list(set(af).intersection(dp)))
remaining_samples = [x + ":AF" for x in overlap] + [x + ":DP" for x in overlap]
D = D[remaining_samples]
M = M.loc[[x.split("_")[0] + "_" + x.split("_")[1] for x in overlap], :]

###Delete rows containing only 0

#AF values
cols = [x for x in D.columns if ":AF" in x]
AF = D[cols]
AF = AF[(AF.T != 0).any()]
D = D.loc[AF.index, :]
#DP values
cols = [x for x in D.columns if ":DP" in x]
DP = D[cols]
DP = DP[(DP.T != 0).any()]
D = D.loc[DP.index, :]

###Threshold for DP: Minimum number of samples in which a mutation must be present

if args.min_samples > -1:
    thresh = args.min_samples
    D["temp"] = (D[[x for x in D.columns if ":DP" in x]] == 0).astype(int).sum(axis=1)
    D = D[D["temp"] <= len(D.columns) / 2 - thresh]
    del D["temp"]

###Threshold for DP: Percentage of detected mutations in a sample to be retained

if args.undetected_mutations != -1:
    thresh = args.undetected_mutations
    for col in D[[x for x in D.columns if ":DP" in x]].columns:
        if (D[col] == 0).sum() / len(D.index) > thresh:
            del D[col]
            del D[col.strip(":DP") + ":AF"]
            M = M.drop(index=(col.split("_")[0] + "_" + col.split("_")[1]))

###Threshold for DP: Coefficient of variation per sample

if args.coefficient_variant > -1:
    thresh = args.coefficient_variant
    for col in D[[x for x in D.columns if ":DP" in x]].columns:
        if (D[col].std() / D[col].mean()) > thresh:
            del D[col]
            del D[col.strip(":DP") + ":AF"]
            M = M.drop(index=(col.split("_")[0] + "_" + col.split("_")[1]))

###End of the program

#Report at the end of the analysis
print("\nEND OF THE ANALYSIS:")
print("Number of samples: ", int(len(D.columns) / 2))
print("Number of mutations: ", len(D.index))
#Output of filtered data file and metadata file
if args.name == None:
    name = ""
else:
    name = args.name + "_"
if args.output_dir == None:
    D.to_csv(name + "DATA.txt", sep="\t")
    M.to_csv(name + "META.txt", sep="\t")
else:
    args.output_dir = "/" + args.output_dir.strip("/") + "/"
    D.to_csv(args.output_dir + name + "DATA.txt", sep="\t")
    M.to_csv(args.output_dir + name + "META.txt", sep="\t")
#Print duration of the run
print("\nDuration: ", int(round(time.time() - start, 0)), " s")
