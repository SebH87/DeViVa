# Copyright (C) 2021 Sebastian Hupfauf

import argparse
import sys
import pandas as pd
import zenipy as zp
from datetime import date

##Data import

#Setup of argparse
#parser = argparse.ArgumentParser(description='This tool can be used to filter a virus mutation table.')
#requiredNamed = parser.add_argument_group('required arguments')
#requiredNamed.add_argument('-d', '--data', action="store", type=str, required=True, help='Name of the input data file.')
#requiredNamed.add_argument('-m', '--metadata', action="store", type=str, required=True, help='Name of the input metadata file.')
##parser.add_argument('-k', '--key', action="store", type=str, help='specified key determining the output')
##parser.add_argument('-v', '--verbose', action="store_true", help='specified key determining the output')
#args = parser.parse_args()
#data = args.data
#meta = args.metadata

#Selecting input files with GUI
data = zp.file_selection(title="Please select the data file in tab-delimited format:")
meta = zp.file_selection(title="Please select the metadata file in tab-delimited format:")
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
#User input and application
q = zp.question(title="Mutation Filtering", text="Do you want to select a subset of epidemiological weeks? Click 'No' in order to analyze all timepoints!")
if q == True:
    beg = int(zp.scale(title="Mutation Filtering", text="Ranging from epidemiological week:", value=M["ep_week"].min(), min=M["ep_week"].min(), max=M["ep_week"].max(), step=1))
    ex = int(zp.scale(title="Mutation Filtering", text="Until epidemiological week:", value=beg, min=beg, max=M["ep_week"].max(), step=1))
    M = M[(M["ep_week"] >= beg) & (M["ep_week"] <= ex)]
    for col in D.columns:
        if not (col.split("_")[0] + "_" + col.split("_")[1]) in M.index:
            del D[col]

###Making a subset based on regions

#User input and application
q = zp.question(title="Mutation Filtering", text="Do you want to select a subset of regions? Click 'No' in order to analyze all regions!")
if q == True:
    regs = M["sample_source_location_state"].unique()
    regions = [str(a + 1) + " " + b for (a, b) in zip(range(len(regs)), regs)]
    selection = zp.entry(title="Mutation Filtering", text="Please enter the NUMBER of the regions you want to include! e.g. 1,2,4\n\n" + "\n".join(regions))
    sel = [regs[i] for i in [int(j) - 1 for j in selection.split(",")]]
    M = M[M["sample_source_location_state"].isin(sel)]
    for col in D.columns:
        if not (col.split("_")[0] + "_" + col.split("_")[1]) in M.index:
            del D[col]

###Making a subset based on locations

#User input and application
q = zp.question(title="Mutation Filtering", text="Do you want to select a subset of locations? Click 'No' in order to analyze all locations!")
if q == True:
    locs = M["LocationName_coronA"].unique()
    locations = [str(a + 1) + " " + b for (a, b) in zip(range(len(locs)), locs)]
    selection = zp.entry(title="Mutation Filtering", text="Please enter the NUMBER of the locations you want to include! e.g. 1,2,4\n\n" + "\n".join(locations))
    sel = [locs[i] for i in [int(j) - 1 for j in selection.split(",")]]
    M = M[M["LocationName_coronA"].isin(sel)]
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

#User input and input check
q = zp.question(title="Mutation Filtering", text="Do you want to apply a threshold for read depth? Values below the threshold will be treated as absent.")
if q == True:
    thresh = zp.entry(title="Mutation Filtering", text="Threshold:")
    try:
        thresh = float(int(thresh))
    except:
        zp.error(title="Mutation Filtering", text="Invalid input! Input must be an integer number!")
        sys.exit(1)
#Application
    for col in D.columns:
        if ":DP" in col:
            D[col] = D[col].apply(lambda x: 0.0 if x < thresh else x)
            D[col.split(":")[0] + ":AF"].loc[D.loc[D[col] == 0].index] = 0

###Applying a threshold for the allele frequency

#User input and input check
q = zp.question(title="Mutation Filtering", text="Do you want to apply a threshold for allele frequency? Values below the threshold will be treated as absent.")
if q == True:
    thresh = zp.entry(title="Mutation Filtering", text="Threshold: [0.0 - 1.0]")
    try:
        thresh = float(thresh)
    except:
        zp.error(title="Mutation Filtering", text="Invalid input! Input must be a number!")
        sys.exit(1)
#Application
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

#User input and input check
q = zp.question(title="Mutation Filtering", text="Do you want to set a minimum number of samples in which a mutation must be present in order to be retained?")
if q == True:
    thresh = zp.entry(title="Mutation Filtering", text="Threshold: [0 - %i]"%(len(D.columns) / 2))
    try:
        thresh = int(thresh)
    except:
        zp.error(title="Mutation Filtering", text="Invalid input! Input must be an integer number!")
        sys.exit(1)
#Application
    D["temp"] = (D[[x for x in D.columns if ":DP" in x]] == 0).astype(int).sum(axis=1)
    D = D[D["temp"] <= len(D.columns) / 2 - thresh]
    del D["temp"]

###Threshold for DP: Percentage of detected mutations in a sample to be retained

#User input and input check
q = zp.question(title="Mutation Filtering", text="Do you want to filter out samples with a high percentage of undetected mutations? Samples with a higher percentage are excluded from the analysis!")
if q == True:
    thresh = zp.entry(title="Mutation Filtering", text="Threshold: [0.0 - 1.0]")
    try:
        thresh = float(thresh)
    except:
        zp.error(title="Mutation Filtering", text="Invalid input! Input must be a number!")
        sys.exit(1)
#Application
    for col in D[[x for x in D.columns if ":DP" in x]].columns:
        if (D[col] == 0).sum() / len(D.index) > thresh:
            del D[col]
            del D[col.strip(":DP") + ":AF"]
            M = M.drop(index=(col.split("_")[0] + "_" + col.split("_")[1]))

###Threshold for DP: Coefficient of variation per sample

#User input and input check
q = zp.question(title="Mutation Filtering", text="Do you want to filter out samples with a high coefficient of variation (CV)? Samples with a higher CV are excluded from the analysis!")
if q == True:
    thresh = zp.entry(title="Mutation Filtering", text="Threshold: e.g. 0.01 or 0.1")
    try:
        thresh = float(thresh)
    except:
        zp.error(title="Mutation Filtering", text="Invalid input! Input must be a number!")
        sys.exit(1)
#Application
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
D.to_csv("DATA.txt", sep="\t")
M.to_csv("META.txt", sep="\t")
#Final message
zp.message(title="Mutation Filtering", text="Thank you for using this filtering tool!\n\n(C) 2021\nSebastian Hupfauf")
