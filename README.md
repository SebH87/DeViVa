# DeViVa (Deconvolution of Virus Variants)

DeViVa is a tool designed for the _de-novo_ detection of novel emerging virus variants (e.g., SARS-CoV-2). Associated mutations are conflated by their corresponding frequency pattern. With DeViVa, it is possible to designate different mutation constellations based on their observed frequencies across different samples in time and space by a hierarchical, unsupervised two-step clustering approach. The produced constellations of correlated mutations can then be characterised by either incorporating them into the SARS-CoV-2 reference genome with subsequent variant typing using the pangolin software (https://cov-lineages.org/resources/pangolin.html), or by assigning them to a specific reference database for virus variants.

DeViVa was originally designed for data from an amplicon-based sequencing approach relying on a modified version of the ARTIC primer set and amplicons of around 400 bases, but can easily be used also for other sequencing data. However, data need to be formatted in the right way (more information will follow soon)!

__USAGE:__

DeViVa is a Python script that can be run on any operating system. Please ensure that you are using Python 3.x! To start a DeViVa analysis, go to a terminal and type in:

_python DeViVa.py -d \<data file\>_

You can add additional arguments if needed, but -d is always needed!

__DEPENDENCIES:__

DeViVa needs some additional Python packages to be installed. You can install them manually or use the DeViVa.install file when working on Linux. Here is a full list of all required Python packages:

* Bio

* matplotlib

* numpy

* pandas

* scikit-bio

* scikit-learn

* scipy

* seaborn

__OPTIONS:__

With the following options, you can adjust your DeViVa run. All parameters are optional except for -d, which needs to be always provided! Here is a list of all available options:

__-a, --assignment__ &emsp; Strategy to assign mutations of a cluster to a SARS-CoV-2 variant. Available options: pangolin (Phylogenetic Assignment of Named Global Outbreak Lineages) or individual (each individual mutation is assigned to a SARS-CoV-2 variant and the most common hit is presented as most likely variant for the cluster). [Default: pangolin]

__-c, --colors__ &emsp; Define colours for all plots (dendrogram, bar charts). You can choose between all common Python colormaps (https://matplotlib.org/stable/tutorials/colors/colormaps.html), e.g., 'viridis'. Use 'standard' for the default colour style of the clustering tool! [Default: standard]

__-c1, --clusters1__ &emsp; Number of Clusters in the first cluster analysis. This is important to correctly distinguish between mutations of interest (MOIs) and background noise [Default: 3]

__-c1mth, --c1method__ &emsp; Method used for the first cluster analysis, available options are: single, complete, average, weighted, centroid, median, and ward. If you want to use multiple methods, separate them by a ',', e.g. 'complete,average'. ATTENTION: Only the last metric/method combination will be used as starting point for the second cluster analysis step! ATTENTION: Methods centroid, median, and ward are correctly defined only for Euclidean distance and hence, no other metric will be accepted in these cases! [Default: ward]

__-c1mtr, --c1metric__ &emsp; Metric used for the first cluster analysis, available options are: euclidean, minkowski, cityblock, seuclidean, sqeuclidean, cosine, correlation, hamming, jaccard, chebyshev, canberra, braycurtis, mahalanobis, yule, matching, dice, kulsinski, rogerstanimoto, russellrao, sokalmichener, sokalsneath, and wminkowski. If you want to use multiple metrics, separate them by a ',', e.g. 'euclidean,braycurtis'. ATTENTION: Only the last metric/method combination will be used as starting point for the second cluster analysis step! [Default: euclidean]

__-c2, --clusters2__ &emsp; Number of Clusters in the second cluster analysis (= assumed virus variants). Use -1 to use the optimum cluster count according to the silhouette analysis! [Default: -1]

__-c2mth, --c2method__ &emsp; Method used for the second cluster analysis, available options are: single, complete, average, weighted, centroid, median, and ward. ATTENTION: Methods centroid, median, and ward are correctly defined only for Euclidean distance and hence, no other metric will be accepted in these cases! [Default: complete]

__-c2mtr, --c2metric__ &emsp; Metric used for the second cluster analysis, available options are: euclidean, minkowski, cityblock, seuclidean, sqeuclidean, cosine, correlation, hamming, jaccard, chebyshev, canberra, braycurtis, mahalanobis, yule, matching, dice, kulsinski, rogerstanimoto, russellrao, sokalmichener, sokalsneath, and wminkowski. [Default: seuclidean]

__-d, --data__ &emsp; Name of the input data file.

__-dpi, --dpi__ &emsp; Pixel density for all plots; this is only used if a raster file format is selected at --format! [Default: 100]

__-ds, --data_separator__ &emsp; Separator for the data input file, e.g. ',' or ';'. ATTENTION: You need to use quotation marks here! For tabulator-separated data use TAB. [Default: TAB]

__-f, --format__ &emsp; Fileformat for all output graphics. You can choose between eps, jpeg, pdf, png, ps, raw, rgba, svg, svgz, and tiff. You can provide multiple file formats by separating them with a comma sign, e.g., 'jpeg,tiff'. [Default: jpeg]

__-ft, --frequency_threshold__ &emsp; Minimum mean frequency (= relative abundance) of all mutations of a cluster after the first cluster analysis. Mutations of clusters fulfilling this criterion are added to the list of mutations of interest (MOIs), while the others are categorized as bias. ATTENTION: Input must be a number between 0 and 1! [Default: 0.05]

__-h, --help__ &emsp; Show help message and exit.

__-l, --label__ &emsp; Label mutations in the second HCA either based on the nucleotide or amino acid level. ATTENTION: In case of amino acid based labels, a reference file ('aa_label.txt') must be provided! For creating the reference file, please use the 'ProLab.py' script from the DeViVa GitHub page. Available options: 'nucleotide', 'aa' [Default: nucleotide]

__-m, --metadata__ &emsp; Name of the input metadata file. A metadata file is only needed when PCoA plots should be coloured based on metadata!

__-mf, --mutation_file__ &emsp; Reference file for variant assignment when using the 'individual' mode from --assignment. Keep in mind that the mutation file must be stored in the same directory as the DeViVa Python script! You can download the default file from the DeViVa GitHub repository. [Default: mutations_list_20210519.csv]

__-mfs, --mutation_file_separator__ &emsp; Separator for the mutation file(s), e.g. ',' or ';'. ATTENTION: You need to use quotation marks here! For tabulator-separated data, use TAB. ATTENTION: This option is only needed when using the 'individual' mode from --assignment! [Default: ',']

__-ms, --metadata_separator__ &emsp; Separator for the metadata input file, e.g. ',' or ';'. ATTENTION: You need to use quotation marks here! For tabulator-separated data use TAB. [Default: TAB]

__-mv, --meta_variable__ &emsp; Metadata variable that is used to colour the data points in the PCoA plot(s) and to create bar charts (VOC composition) specifically for this parameter. If you want to use multiple variables, separate them by comma, e.g. 'sample_source_location_state,ep_week'. Provide 'None' if you do not want to colour the data points and create specific charts. [default: None]

__-omtr, --ometric__ &emsp; Metric used for ordination analysis, available options are: euclidean, minkowski, cityblock, seuclidean, sqeuclidean, cosine, correlation, hamming, jaccard, chebyshev, canberra, braycurtis, mahalanobis, yule, matching, dice, kulsinski, rogerstanimoto, russellrao, sokalmichener, sokalsneath, and wminkowski. If you want to use multiple metrics, separate them by comma, e.g., 'euclidean,braycurtis'. Leave this parameter out if you do not want to compute an ordination analysis!

__-os, --output_separator__ &emsp; Separator for all output files, e.g. ',' or ';'. ATTENTION: You need to use quotation marks here! For tabulator-separated data use TAB. [Default: TAB]

__-p, --permutations__ &emsp; Number of permutations performed during the ANOSIM and PERMANOVA analysis when doing ordination analysis with metadata variables. [Default: 999]

__-rf, --reference_file__ &emsp; Reference file for the construction of the FASTA file when using the 'pangolin' mode from --assignment. Provide 'auto' to automatically download the suggested SARS-CoV-2 file from NCBI. [Default: auto]

__-st, --silhouette_threshold__ &emsp; Number of clusters to be covered in the silhouette analysis of the second cluster analysis in order to find the best number of clusters. [Default: 20]

__-t, --transformation__ &emsp; Data transformation of relative mutation counts per sample between the first and the second HCA. You can choose between the following options: clr (centre log ration) and inv_clr (inverse centre log ration). Use 'original' if you do not want to apply any data transformation. [Default: original]

__-u, --usher__ &emsp; Use UShER mode instead of default pangoLEARN when assigning clusters with pangolin. ATTENTION: This option is only needed when --assignment is set to 'pangolin'!

__-v, --version__ &emsp; Show program's version number and exit.
