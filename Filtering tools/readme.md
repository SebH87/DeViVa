# DeViVa (Deconvolution of Virus Variants)

DeViVa is a tool designed for the _de-novo_ detection of novel emerging virus variants (e.g., SARS-CoV-2). Associated mutations are conflated by their corresponding frequency pattern. With DeViVa, it is possible to designate different mutation constellations based on their observed frequencies across different samples in time and space by a hierarchical, unsupervised two-step clustering approach. The produced constellations of correlated mutations can then be characterised by either incorporating them into the SARS-CoV-2 reference genome with subsequent variant typing using the pangolin software (https://cov-lineages.org/resources/pangolin.html), or by assigning them to a specific reference database for virus variants.

DeViVa was originally designed for data from an amplicon-based sequencing approach relying on a modified version of the ARTIC primer set and amplicons of around 400 bases, but can easily be used also for other sequencing data. However, data need to be formatted in the right way (more information will follow soon)!

__USAGE:__

DeViVa is a Python script that can be run on any operating system. Please ensure that you are using Python 3.x! To start a DeViVa analysis, go to a terminal and type in:

_python DeViVa.py -d_

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
