# Copyright (C) 2022 Sebastian Hupfauf

import argparse
import time
import sys
import os
import collections
import shutil
import matplotlib
import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from Bio import SeqIO, Entrez
from matplotlib import cm
from matplotlib.lines import Line2D
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, to_tree
from scipy.spatial.distance import pdist
from skbio.diversity import beta_diversity
from skbio.stats.composition import clr, clr_inv, ilr, ilr_inv
from skbio.stats.ordination import pcoa
from skbio.stats.distance import anosim, permanova
from sklearn.metrics import silhouette_score

start = time.time()

###Definition of functions

def var_call(L, ref):
    mut_list = []
    for mut in L:
        mut_list.append(re.split('(\d+)', mut))
    mut_list = [[a, int(b), c] for a, b, c in mut_list]
    mut_list.sort(key = lambda x: x[1])
    reff = open(ref)
    reff.readline()
    ref_str = reff.read().replace("\n", "")
    reff.close()
    new_str = ""
    ctrl = 0
    extra = 0
    for mut in mut_list:
        if mut[1] == ctrl:
            print("ATTENTION! There are multiple mutations at position ", ctrl, "!")
        new_str += ref_str[(ctrl + extra) : (mut[1] - 1)]
        new_str += mut[2]
        if len(mut[0]) > 1:
            extra = len(mut[0]) - 1
        else:
            extra = 0
        ctrl = mut[1]
    new_str += ref_str[ctrl:]
    return new_str

def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

##Data import

#Setup of argparse
parser = argparse.ArgumentParser(description='DeViVa can be used to analyze a virus mutation table.', add_help=False, epilog='Thank you for using DeViVa!\n\n(C) 2021\nSebastian Hupfauf')
requiredNamed = parser.add_argument_group('Required arguments')
optionalNamed = parser.add_argument_group('Optional arguments')
requiredNamed.add_argument('-d', '--data', action="store", type=str, required=True, help='Name of the input data file.')
optionalNamed.add_argument('-m', '--metadata', action="store", type=str, help='Name of the input metadata file. A metadata file is only needed when PCoA plots should colored based on metadata!')
optionalNamed.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Show this help message and exit.')
optionalNamed.add_argument('-v', '--version', action='version', version='%(prog)s 1.0', help="Show program's version number and exit.")
optionalNamed.add_argument('-ds', '--data_separator', action='store', default='TAB', type=str, help="Separator for the data input file, e.g. ',' or ';'. ATTENTION: You need to use quotation marks here! For tabulator-separated data use TAB. [Default: TAB]")
optionalNamed.add_argument('-ms', '--metadata_separator', action='store', default='TAB', type=str, help="Separator for the metadata input file, e.g. ',' or ';'. ATTENTION: You need to use quotation marks here! For tabulator-separated data use TAB. [Default: TAB]")
optionalNamed.add_argument('-os', '--output_separator', action='store', default='TAB', type=str, help="Separator for all output files, e.g. ',' or ';'. ATTENTION: You need to use quotation marks here! For tabulator-separated data use TAB. [Default: TAB]")
optionalNamed.add_argument('-c1mtr', '--c1metric', action='store', default='euclidean', type=str, help="Metric used for the first cluster analysis, available options are: euclidean, minkowski, cityblock, seuclidean, sqeuclidean, cosine, correlation, hamming, jaccard, chebyshev, canberra, braycurtis, mahalanobis, yule, matching, dice, kulsinski, rogerstanimoto, russellrao, sokalmichener, sokalsneath, and wminkowski. If you want to use multiple metrics, separate them by a ',', e.g. 'euclidean,braycurtis'. ATTENTION: Only the last metric/method combination will be used as starting point for the second cluster analysis step! [Default: euclidean]")
optionalNamed.add_argument('-c1mth', '--c1method', action='store', default='ward', type=str, help="Method used for the first cluster analysis, available options are: single, complete, average, weighted, centroid, median, and ward. If you want to use multiple methods, separate them by a ',', e.g. 'complete,average'. ATTENTION: Only the last metric/method combination will be used as starting point for the second cluster analysis step! ATTENTION: Methods centroid, median, and ward are correctly defined only for euclidean distance and hence, no other metric will be accepted in these cases! [Default: ward]")
optionalNamed.add_argument('-c2mtr', '--c2metric', action='store', default='seuclidean', type=str, help="Metric used for the second cluster analysis, available options are: euclidean, minkowski, cityblock, seuclidean, sqeuclidean, cosine, correlation, hamming, jaccard, chebyshev, canberra, braycurtis, mahalanobis, yule, matching, dice, kulsinski, rogerstanimoto, russellrao, sokalmichener, sokalsneath, and wminkowski. [Default: seuclidean]")
optionalNamed.add_argument('-c2mth', '--c2method', action='store', default='complete', type=str, help="Method used for the second cluster analysis, available options are: single, complete, average, weighted, centroid, median, and ward. ATTENTION: Methods centroid, median, and ward are correctly defined only for euclidean distance and hence, no other metric will be accepted in these cases! [Default: complete]")
optionalNamed.add_argument('-c1', '--clusters1', action='store', default=3, type=int, help="Number of Clusters in the first cluster analysis. This is important to correctly distinguish between mutations of interest (MOIs) and background noise [Default: 3]")
optionalNamed.add_argument('-c2', '--clusters2', action='store', default=-1, type=int, help="Number of Clusters in the second cluster analysis (= assumed virus variants). Use -1 to use the optimum cluster count according to the silhouette analysis! [Default: -1]")
optionalNamed.add_argument('-a', '--assignment', action='store', default='pangolin', type=str, help="Strategy to assign mutations of a cluster to a SARS-CoV-2 variant. Available options: pangolin (Phylogenetic Assignment of Named Global Outbreak Lineages), individual (each individual mutation is assigned to a SARS-CoV-2 variant and the most common hit is presented as most likely variant for the cluster). [Default: pangolin]")
optionalNamed.add_argument('-rf', '--reference_file', action='store', default='auto', type=str, help="Reference file for the construction of the FASTA file when using the 'pangolin' mode from --assignment. Provide 'auto' to automatically download the suggested SARS-CoV-2 file from NCBI. [Default: auto]")
optionalNamed.add_argument('-mf', '--mutation_file', action='store', default='mutations_list_20210519.csv', type=str, help="Reference file for variant assignment when using the 'individual' mode from --assignment. [Default: mutations_list_20210519.csv]")
optionalNamed.add_argument('-mfs', '--mutation_file_separator', action='store', default=',', type=str, help="Separator for the mutation file(s), e.g. ',' or ';'. ATTENTION: You need to use quotation marks here! For tabulator-separated data use TAB. ATTENTION: This option is only needed when using the 'individual' mode from --assignment! [Default: ',']")
optionalNamed.add_argument('-u', '--usher', action='store_true', help="Use UShER mode instead of default pangoLEARN when assigning clusters with pangolin. ATTENTION: This option is only needed when --assignment is set to 'pangolin'!")
optionalNamed.add_argument('-ft', '--frequency_threshold', action='store', default=0.05, type=float, help="Minimum mean frequency (= relative abundance) of all mutations of a cluster after the first cluster analysis. Mutations of clusters fulfilling this criterion are added to the list of MOIs, while the others are categorized as bias. ATTENTION: Input must be a number between 0 and 1! [Default: 0.05]")
optionalNamed.add_argument('-st', '--silhouette_threshold', action='store', default=20, type=int, help="Number of clusters to be covered in the silhouette analysis of the second cluster analysis in order to find the best number of clusters. [Default: 20]")
optionalNamed.add_argument('-omtr', '--ometric', action='store', type=str, help="Metric used for ordination analysis, available options are: euclidean, minkowski, cityblock, seuclidean, sqeuclidean, cosine, correlation, hamming, jaccard, chebyshev, canberra, braycurtis, mahalanobis, yule, matching, dice, kulsinski, rogerstanimoto, russellrao, sokalmichener, sokalsneath, and wminkowski. If you want to use multiple metrics, separate them by comma, e.g., 'euclidean,braycurtis'. Leave this parameter out if you do not want to compute an ordination analysis!")
optionalNamed.add_argument('-mv', '--meta_variable', action='store', type=str, default='None', help="Metadata variable that is used to colour the data points in the PCoA plot(s) and to create bar charts (VOC composition) specifically for this parameter. If you want to use multiple variables, separate them by comma, e.g., 'sample_source_location_state,ep_week'. Provide 'None' if you do not want to colour the data points and create specific charts. [default: None]")
optionalNamed.add_argument('-p', '--permutations', action='store', default=999, type=int, help="Number of permutations performed during the ANOSIM and PERMANOVA analysis when doing ordination analysis with metadata variables. [Default: 999]")
optionalNamed.add_argument('-c', '--colors', action='store', default='standard', type=str, help="Define colours for all plots (dendrogram, bar charts). You can choose between all common Python colormaps (https://matplotlib.org/stable/tutorials/colors/colormaps.html), e.g., 'viridis'. Use 'standard' for the default colour style of the clustering tool! [Default: standard]")
optionalNamed.add_argument('-f', '--format', action='store', default='jpeg', type=str, help="Fileformat for all output graphics. You can choose between eps, jpeg, pdf, png, ps, raw, rgba, svg, svgz, and tiff. You can provide multiple file formats by separating them with a comma sign, e.g. 'jpeg,tiff'. [Default: jpeg]")
optionalNamed.add_argument('-dpi', '--dpi', action='store', default=100, type=int, help="Pixel density for all plots; this is only used if a raster file format is selected at --format! [Default: 100]")
optionalNamed.add_argument('-l', '--label', action='store', default='nucleotide', type=str, help="Label mutations in the second HCA either based on the nucleotide or amino acid level. ATTENTION: In case of amino acid based labels, a reference file ('aa_label.txt') must be provided! For creating the reference file, please use the 'ProLab.py' script. Available options: 'nucleotide', 'aa' [Default: nucleotide]")
optionalNamed.add_argument('-t', '--transformation', action='store', default='original', type=str, help="Data transformation of relative mutation counts per sample between the first and the second HCA. You can choose between the following options: clr (centre log ration) and inv_clr (inverse centre log ration). Use 'original' if you do not want to apply any data transformation. [Default: original]") ###, ilr (isometric log ration), and ilr_inv (inverse isometric log ratio)
args = parser.parse_args()
data = args.data
#Check of input
if args.label not in ["aa", "nucleotide"]:
    print("ATTENTION: Incorrect labelling, -l must be either nucleotide or aa! Process terminated!")
    sys.exit(0)
#Creating dataframes
if args.meta_variable != "None":
    try:
        meta = args.metadata
        if args.metadata_separator == "TAB":
            args.metadata_separator = "\t"
        m = pd.read_csv(meta, delimiter=args.metadata_separator, dtype=str, index_col=0)
    except:
        print("ATTENTION: You need to provide a compatible metadata file!")
        sys.exit(1)
if args.data_separator == "TAB":
    args.data_separator = "\t"
if args.output_separator == "TAB":
    args.output_separator = "\t"
d = pd.read_csv(data, delimiter=args.data_separator, dtype=str, index_col=0)
#Convert values to float
d=d.astype(float)
#Keep only rows containing the allele frequency
d = d[[x for x in d.columns if not "DP" in x]]
#Bring sample names to the format CoV_xxxx
d.columns = [x.split("_")[0] + "_" + x.split("_")[1] for x in d.columns]
#Remove all-zero rows
d = d[(d.T != 0.0).any()]
#Check if data was imported correctly
if len(d.columns) == 0:
    print("ATTENTION: No data detected, please check your input file or consider using another data separator!")
    sys.exit(0)
if len(d.index) == 0:
    print("ATTENTION: No data detected, please check your input file or consider using another data separator!")
    sys.exit(0)

###Calculate first HCA

if args.c1method != None and args.c1metric != None:
    matplotlib.use("Agg")
    for meth in args.c1method.split(","):
        for metr in args.c1metric.split(","):
            plt.figure(figsize=(20, 10))
            lmatrix = linkage(d, method=meth, metric=metr)
            x = dendrogram(lmatrix, leaf_rotation=90., labels=d.index, above_threshold_color='black', truncate_mode='lastp', p=args.clusters1)
            plt.tight_layout()
            for ff in args.format.split(","):
                plt.savefig("hca1_" + meth + "_" + metr + "." + ff, dpi=args.dpi)
                print("\nGraphic created: hca1_" + meth + "_" + metr + "." + ff)
            plt.clf()

###Selection of MOIs

cluster = fcluster(lmatrix, args.clusters1, criterion='maxclust')
cluster_index = list(zip(d.index, cluster))
mois = []
bias = []
for n in range(args.clusters1):
    cl = [k for k, j in cluster_index if j == (n + 1)]
    df = d.loc[cl]
    if df.mean().mean() > args.frequency_threshold:
        mois += cl
    else:
        bias += cl
biasd = d.loc[bias]
biasd.to_csv("bias.txt", sep=args.output_separator)
print("\nFile created: bias.txt")
d = d.loc[mois]
d.to_csv("mois.txt", sep=args.output_separator)
print("\nFile created: mois.txt\n")

###Data transformation
if args.transformation != "original":
    col_labels = d.columns
    row_labels = d.index
    d = d.replace(0, 0.0001)
if args.transformation == "clr":
    d = pd.DataFrame(clr(d), index=row_labels)
    d.columns = col_labels
if args.transformation == "clr_inv":
    d = pd.DataFrame(clr_inv(d), index=row_labels)
    d.columns = col_labels
#if args.transformation == "ilr":
#    d = pd.DataFrame(ilr(d), index=row_labels)
#    d.columns = col_labels
#if args.transformation == "ilr_inv":
#    d = pd.DataFrame(ilr_inv(d), index=row_labels)
#    d.columns = col_labels

###Second HCA

lmatrix2 = linkage(d, method=args.c2method, metric=args.c2metric)
#Silhouette Analysis
X = []
Y = []
for k in range(args.silhouette_threshold - 1):
    cluster = fcluster(lmatrix2, k + 2, criterion='maxclust')
    silhouette_avg = silhouette_score(d, cluster)
    X.append(k + 2)
    Y.append(silhouette_avg)
plt.bar(X, Y)
plt.xticks(X, [str(x) for x in X])
plt.xlabel("# of clusters")
plt.ylabel("Mean silhouette score")
plt.tight_layout()
for ff in args.format.split(","):
    plt.savefig("silhouette." + ff, dpi=args.dpi)
    print("Graphic created: silhouette." + ff + "\n")
plt.clf()
if args.clusters2 < 0:
    args.clusters2 = Y.index(max(Y)) + 2
#Ausgabe der gefundenen Varianten
cluster2 = fcluster(lmatrix2, args.clusters2, criterion='maxclust')
cluster2_index = list(zip(d.index, cluster2))
try:
    os.mkdir("variants")
except:
    shutil.rmtree("variants", ignore_errors=True)
    os.mkdir("variants")
fasta_outf = open("variants.fasta", "w")
if args.assignment == "individual":
    if args.mutation_file_separator == "TAB":
        args.mutation_file_separator = "\t"
    ref_df = pd.read_csv(args.mutation_file, sep=args.mutation_file_separator, index_col=0)
    var_outfile = open("variants.txt", "w")
    var_outfile.write("ID" + args.output_separator + "Variant" + args.output_separator + "Frquency\n")
if args.reference_file == "auto":
    Entrez.email = "sebastian.hupfauf@uibk.ac.at"
    hdl = Entrez.efetch(db="nucleotide", id="NC_045512.2", rettype="gb")
    recs = list(SeqIO.parse(hdl, "gb"))
    SeqIO.write(recs, "pangolin_reference.fasta", "fasta")
    print("File created: pangolin_reference.fasta\n")
    args.reference_file = "pangolin_reference.fasta"
labels = []
if args.colors != "standard":
    if args.clusters2 < 10:
        colors = list(cm.get_cmap(args.colors)(np.linspace(0, 1, args.clusters2)))
    else:
        colors = list(cm.get_cmap(args.colors)(np.linspace(0, 1, 10)))
    cols = (colors * ((args.clusters2 // 10) + 1))[:args.clusters2][::-1]
else:
    cols = (sns.color_palette() * ((args.clusters2 // 10) + 1))[:args.clusters2][::-1]
colorlist = []
var_df = pd.DataFrame(index=d.columns)
for n in range(args.clusters2):
    cl2 = [k for k, j in cluster2_index if j == (n + 1)]
    df2 = d.loc[cl2]
    #Creating colorlist for second HCA
    if len(cl2) > 1:
        colorlist.append(cols.pop())
    #Calculating variant counts for each sample
    if len(cl2) > 1:
        df2_temp = df2.replace(0.0, np.NaN)
        var_df["var_" + str(n + 1)] = df2_temp.median(axis = 0).replace(np.NaN, 0.0)
    #Berechnung der wahrscheinlichsten Variante "individual"
    if args.assignment == "individual":
        if len(cl2) > 1:
            assigned_var = []
            for mut in df2.index:
                if list(ref_df[ref_df["AA"] == mut].index) != []:
                    assigned_var += [k.replace(";", ",") for k in list(ref_df[ref_df["AA"] == mut].index)]
                else:
                    assigned_var.append("")
            df2["assigned_var"] = assigned_var
    #Ausgabe der Dataframes für jeden Cluster
    if len(cl2) > 1:
        df2.to_csv("variants/variant_" + str(n + 1) + ".txt", sep=args.output_separator)
    #Ausgabe der FASTA Files
    if len(cl2) > 1:
        fasta_outf.write(">Variant_" + str(n + 1) + "\n")
        fasta_outf.write(var_call(df2.index, args.reference_file) + "\n")
    #Ausgabe der wahrscheinlichsten Variante "individual"
    if args.assignment == "individual":
        if len(cl2) > 1:
            var_flatten = [m for sublist in [x.split(",") for x in assigned_var] for m in sublist]
            var_outfile.write("Variant_" + str(n + 1) + args.output_separator)
            if len(collections.Counter(var_flatten).most_common(10)) == 10:
                if collections.Counter(var_flatten).most_common(10)[0][1] == collections.Counter(var_flatten).most_common(10)[9][1]:
                    var_outfile.write("Basal mutations (>10 lineages)" + args.output_separator)
                    labels.append("Basal mutations (>10 lineages)")
                elif collections.Counter(var_flatten).most_common(10)[0][0] != "":
                    var_outfile.write(collections.Counter(var_flatten).most_common(10)[0][0] + args.output_separator)
                    labels.append(collections.Counter(var_flatten).most_common(10)[0][0] + " (" + str(round(int(collections.Counter(var_flatten).most_common(1)[0][1]) / len(df2.index), 2)) + ")")
                else:
                    var_outfile.write("?" + args.output_separator)
                    labels.append("?")
            else:
                if collections.Counter(var_flatten).most_common(1)[0][0] != "":
                    var_outfile.write(collections.Counter(var_flatten).most_common(1)[0][0] + args.output_separator)
                    labels.append(collections.Counter(var_flatten).most_common(1)[0][0] + " (" + str(round(int(collections.Counter(var_flatten).most_common(1)[0][1]) / len(df2.index), 2)) + ")")
                else:
                    var_outfile.write("?" + args.output_separator)
                    labels.append("?")
            var_outfile.write(str(int(collections.Counter(var_flatten).most_common(1)[0][1]) / len(df2.index)) + "\n")
print("\nFile created: variants.fasta")
print("\nFile created: variants.txt")
print("\nFile(s) created: variants/variant_x.txt")
#Ausgabe der variant counts
var_df.to_csv("variants_per_sample.txt" , sep=args.output_separator)
print("\nFile created: variants_per_sample.txt")
#Ausgabe der wahrscheinlichsten Variante "pangolin"
if args.assignment == "pangolin":
    if args.usher == True:
        os.system("pangolin --usher variants.fasta --outfile pangolin_lineages.csv")
    else:
        os.system("pangolin variants.fasta --outfile pangolin_lineages.csv")
    pangolin_df = pd.read_csv("pangolin_lineages.csv", sep=",")
    labels = list(pangolin_df["lineage"])
if args.assignment == "individual":
    var_outfile.close()
fasta_outf.close()
#Barplots for variants per sample
fig = var_df.plot(kind='bar', stacked=True, color=colorlist)
handles, lbls = fig.get_legend_handles_labels()
lgd = plt.legend(handles[::-1], labels[::-1], loc=2, bbox_to_anchor=(1.05, 1), frameon=False)
for ff in args.format.split(","):
    plt.savefig("variants_per_sample." + ff, dpi=args.dpi, bbox_extra_artists=(lgd,), bbox_inches='tight')
    print("\nFile created: variants_per_sample." + ff)
plt.clf()
#Adjusting color threshold
y = dendrogram(lmatrix2, leaf_rotation=90., labels=d.index, above_threshold_color='black', truncate_mode='lastp', p=args.clusters2)
nodes = []
for node in y["dcoord"]:
    nodes.append(node[1])
plt.clf()
#Constructing Dendrogram & label with most likely variant
lab = pd.read_csv("aa_label.txt", sep="\t", index_col=0)
plt.figure(figsize=(20, 10))
if args.colors != "standard":
    hierarchy.set_link_color_palette([mpl.colors.rgb2hex(rgb[:3]) for rgb in colors])
else:
    hierarchy.set_link_color_palette([mpl.colors.rgb2hex(rgb[:3]) for rgb in sns.color_palette()])
if args.label == "aa":
    y = dendrogram(lmatrix2, leaf_rotation=90., labels=lab.loc[d.index]["label"], above_threshold_color='black', color_threshold=min(nodes))
if args.label == "nucleotide":
    y = dendrogram(lmatrix2, leaf_rotation=90., labels=d.index, above_threshold_color='black', color_threshold=min(nodes))
if labels != []:
    legend_elements = []
    for n in range(len(labels)):
        legend_elements.append(Line2D([0], [0], color=colorlist[n], lw=4, label=labels[n]))
    plt.legend(handles=legend_elements, prop={"size":12}, frameon=False)
for ff in args.format.split(","):
    plt.savefig("hca2_" + args.c2method + "_" + args.c2metric + "." + ff, dpi=args.dpi, bbox_inches='tight')
    print("\nGraphic created: hca2_" + args.c2method + "_" + args.c2metric + "." + ff)
plt.clf()
#Export dendrogram as NEWICK tree
tr = to_tree(lmatrix2, False)
if args.label == "aa":
    n_tr = getNewick(tr, "", tr.dist, [x.replace(" (", "_").replace(")", "").replace(":", "/") for x in list(lab.loc[d.index]["label"])])
#print([x.replace(" (", "_").replace(")", "") for x in list(lab.loc[d.index]["label"])])
if args.label == "nucleotide":
    n_tr = getNewick(tr, "", tr.dist, d.index)
with open("tree.tr", "w") as tf:
    tf.write(n_tr)
print("\nFile created: tree.tre")

###Calculate ordination analysis

if args.ometric != None:
    print("\n")
    matplotlib.use("Agg")
    for mtr in args.ometric.split(","):
        if args.meta_variable == "None":
            try:
                div = beta_diversity(mtr, d.T, d.columns)
            except:
                print("ATTENTION: No distance matrix could be constructed with this ordination metric: ", mtr, "!")
                sys.exit(1)
            div_df = div.to_data_frame()
            div_df.to_csv("distancematrix_" + mtr + ".txt", sep=args.output_separator)
            print("\nFile created: distancematrix_" + mtr + ".txt")
            pc = pcoa(div)
            warn = 0
            for ev in pc.eigvals:
                if ev < 0:
                    warn += 1
            if warn > 0:
                print("ATTENTION: At least 1 of your eigenvalues is negative when using the metric: ", mtr, ", potentially leading to problems! You may want to choose another metric for distance calculation or apply data transformation on the distance matrix (e.g. square root) to get rid of this problem.")
            eig_dm = pd.DataFrame(pc.eigvals, columns=["Eigenvalue"])
            eig_dm["Explained"] = pc.proportion_explained
            eig_dm["Summed_explanation"] = pc.proportion_explained.cumsum()
            eig_dm.to_csv("eigenvalues_" + mtr + ".txt", sep=args.output_separator)
            print("\nFile created: eigenvalues_" + mtr + ".txt")
            plt.figure(figsize=[12, 8])
            fig = sns.scatterplot(data=pc.samples, x="PC1", y="PC2")
            plt.xlabel("PC 1 (" + str(round(pc.proportion_explained[0] * 100, 2)) + "%)")
            plt.ylabel("PC 2 (" + str(round(pc.proportion_explained[1] * 100, 2)) + "%)")
            handles, lbls = fig.get_legend_handles_labels()
            lgd = plt.legend(handles[::-1], lbls[::-1], loc=2, bbox_to_anchor=(1.05, 1), frameon=False)
            for ff in args.format.split(","):
                plt.savefig("PCoA_" + mtr + "." + ff, dpi=args.dpi, bbox_extra_artists=(lgd,), bbox_inches='tight')
                print("\nGraphic created: PCoA_" + mtr + "." + ff)
            plt.clf()
        else:
            try:
                m = m.apply(pd.to_numeric, errors='ignore')
            except:
                print("ATTENTION: You need to provide a compatible metadata file!")
                sys.exit(1)
            try:
                div = beta_diversity(mtr, d.T, d.columns)
            except:
                print("ATTENTION: No distance matrix could be constructed with this ordination metric: ", mtr, "!")
                sys.exit(1)
            div_df = div.to_data_frame()
            div_df.to_csv("distancematrix_" + mtr + ".txt", sep=args.output_separator)
            pc = pcoa(div)
            pc_meta = pd.concat([pc.samples, m], axis=1)
            pc_meta = pd.concat([pc_meta, var_df], axis=1)
            pc_meta_freq = pd.concat([pc_meta, d.T], axis=1)
            print("\nFile created: distancematrix_" + mtr + ".txt")
            warn = 0
            for ev in pc.eigvals:
                if ev < 0:
                    warn += 1
            if warn > 0:
                print("ATTENTION: At least 1 of your eigenvalues is negative when using the metric: ", mtr, ", potentially leading to problems! You may want to choose another metric for distance calculation or apply data transformation on the distance matrix (e.g. square root) to get rid of this problem.")
            eig_dm = pd.DataFrame(pc.eigvals, columns=["Eigenvalue"])
            eig_dm["Explained"] = pc.proportion_explained
            eig_dm["Summed_explanation"] = pc.proportion_explained.cumsum()
            eig_dm.to_csv("eigenvalues_" + mtr + ".txt", sep=args.output_separator)
            print("\nFile created: eigenvalues_" + mtr + ".txt")
            for mv in args.meta_variable.split(","):
                #PCoA plot
                plt.figure(figsize=[12, 8])
                fig = sns.scatterplot(data=pc_meta_freq, x="PC1", y="PC2", hue=mv)
                plt.xlabel("PC 1 (" + str(round(pc.proportion_explained[0] * 100, 2)) + "%)")
                plt.ylabel("PC 2 (" + str(round(pc.proportion_explained[1] * 100, 2)) + "%)")
                handles, lbls = fig.get_legend_handles_labels()
                lgd = plt.legend(handles[::-1], lbls[::-1], loc=2, bbox_to_anchor=(1.05, 1), frameon=False)
                for ff in args.format.split(","):
                    plt.savefig("PCoA_" + mtr + "_" + mv + "." + ff, dpi=args.dpi, bbox_extra_artists=(lgd,), bbox_inches='tight')
                    print("\nGraphic created: PCoA_" + mtr + "_" + mv + "." + ff)
                plt.clf()
                #Statistics
                anos = anosim(div, pc_meta_freq, column=mv, permutations=args.permutations)
                perm = permanova(div, pc_meta_freq, column=mv, permutations=args.permutations)
                with open("statistics_" + mtr + "_" + mv + ".txt", "w") as st:
                    st.write("ANOSIM" + args.output_separator + "Permutations: " + str(args.permutations) + "\n\n")
                    st.write("R" + args.output_separator + str(anos["test statistic"]) + "\n")
                    st.write("p-value" + args.output_separator + str(anos["p-value"]) + "\n\n")
                    st.write("PERMANOVA" + args.output_separator + "Permutations: " + str(args.permutations) + "\n\n")
                    st.write("F" + args.output_separator + str(perm["test statistic"]) + "\n")
                    st.write("p-value" + args.output_separator + str(perm["p-value"]) + "\n\n")
                print("\nFile created: statistics_" + mtr + "_" + mv + ".txt")

###Variants per metadata variable:
if args.meta_variable != "None":
    var_meta_df = pd.concat([var_df, m], axis=1)
    for mv in args.meta_variable.split(","):
        var_per_mv = var_meta_df.groupby(mv).median()[var_df.columns]
        var_per_mv.to_csv("variants_" + mv + ".txt", sep=args.output_separator)
        print("\nFile created: variants_" + mv + ".txt")
        #Barplots for variants per MV
        fig = var_per_mv.plot(kind='bar', stacked=True, color=colorlist)
        handles, lbls = fig.get_legend_handles_labels()
        lgd = plt.legend(handles[::-1], labels[::-1], loc=2, bbox_to_anchor=(1.05, 1), frameon=False)
        for ff in args.format.split(","):
            plt.savefig("variants_" + mv + "." + ff, dpi=args.dpi, bbox_extra_artists=(lgd,), bbox_inches='tight')
            print("\nGraphic created: variants_" + mv + "." + ff)
        plt.clf()

###Print duration of the run

print("\nRuntime: ", int(round(time.time() - start, 0)), " s\n")
