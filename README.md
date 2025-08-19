# MSA-and-tree-metrics-exploration
MSA and Tree metrics exploration and MSA simulation of empirical phylogenomic datasets

This document describes the workflow of the MSA-and-tree-metrics-exploration system, detailing how the Python feature extraction tool and R Shiny application work together as an integrated bioinformatics pipeline. The system implements a two-tier architecture where batch processing and interactive analysis are separated into distinct but complementary components. The Python tool handles computationally intensive feature extraction, while the R Shiny application provides user-friendly data exploration and simulation planning capabilities. The two main components communicate through CSV files that serve as the data interchange format. The Python tool produces structured CSV output that the R application can directly consume.

## Feature Extraction Tool 
The Feature Extraction Tool is a Python command-line application that processes multiple sequence alignments (MSAs) and phylogenetic tree files to extract comprehensive statistical and evolutionary features for phylogenomic analysis. This tool serves as the batch processing component of the system, generating standardized feature datasets that feed into the interactive analysis platform.

Required Python libraries are numpy, tqdm, and Biopython, which can be installed simply by:

pip install numpy tqdm biopython

The script extracts two main categories of features: MSA alignments statistics and tree-based branch and distance metrics. MSA features include character ('ACGT' for nucleotides or 'ARNDCQEGHILKMFPSTWYV' for amino acids) frequencies, and 16 other features, such as number of sequences, alignment length, unique site patterns, parsimony-informative sites, gap proportion, average site entropy, proportion of invariant sites, relative composition frequency variation, average pairwise identity, Bollback multinomial etc. Fourteen tree-based features include bralch-related metrics, average bipartition support, Distance Variance from Molecular Clock (DVMC), evolutionary rate, treeness, and normalised Robinson-Foulds distance, and so on.

The extract_msa_tree_features.py script provides a unified command-line interface for extracting features from MSA and tree files. The tool supports three primary modes of operation: MSA-only processing, tree-only processing, and combined MSA+tree processing with file matching by prefix. At least one of --msa_dir or --tree_dir must be provided. It requires consistent file prefixes between MSA and tree files.

usage: extract_msa_tree_features.py

                                    [-h] [--msa_dir MSA_DIR] [--tree_dir TREE_DIR] [--msa_type {AA,DNA}] [-o OUTPUT]
                                    [-t THREADS] [--skip_freq_statistics] [--pseudo_tree_metrics] [--fasttree FASTTREE]
                                    [--outgroup_list OUTGROUP_LIST] [--ref_tree REF_TREE]


options:

  -h, --help&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;show this help message and exit
  
  --msa_dir MSA_DIR&emsp;&emsp;&emsp;&emsp;Directory containing MSA files (FASTA)

  --tree_dir TREE_DIR&emsp;&emsp;&emsp;&emsp;Directory containing tree files (Newick)

  --msa_type {AA,DNA}&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;MSA type: AA (amino acid) or DNA (nucleotide). Required when --msa_dir is provided

  -o OUTPUT, --output OUTPUT&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Output CSV file
  
  -t THREADS, --threads THREADS&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Number of threads for parallel processing
  
  --skip_freq_statistics&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Skip frequency statistics calculation
  
  --pseudo_tree_metrics&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Calculate pseudo tree metrics using faster FastTree execution (-noml -boot 500)
  
  --fasttree FASTTREE&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Path to FastTree executable
  
  --outgroup_list OUTGROUP_LIST&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;File containing outgroup species for DVMC calculation
                        
  --ref_tree REF_TREE&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Reference tree for RF distance calculation



#### Examples:

  \# MSA features only
  
  python extract_msa_tree_features.py --msa_dir /path/to/MSAs --msa_type DNA -o results.csv

  \# Tree features only 
  
  python extract_msa_tree_features.py --tree_dir /path/to/trees -o results.csv

  \# Both MSA and tree features
  
  python extract_msa_tree_features.py --msa_dir /path/to/MSAs --tree_dir /path/to/trees --msa_type DNA -o results.csv

  \# Skip frequency statistics
  
  python extract_msa_tree_features.py --msa_dir /path/to/MSAs --tree_dir /path/to/trees --msa_type DNA --skip_freq_statistics -o results.csv

  \# With pseudo tree metrics
  
  python extract_msa_tree_features.py --msa_dir /path/to/MSAs --msa_type DNA --pseudo_tree_metrics --fasttree /usr/bin/FastTree -o results.csv

  \# With outgroups and reference tree
  
  python extract_msa_tree_features.py --msa_dir /path/to/MSAs --tree_dir /path/to/trees --msa_type DNA --outgroup_list outgroups.txt --ref_tree reference.tree -o results.csv



## Interactive Analysis Platform

The Interactive Analysis Platform is a comprehensive R Shiny web application that provides researchers with an intuitive interface for exploring phylogenomic datasets, performing statistical analysis, and generating MSA simulation commands. This document covers the overall architecture, core functionality, and technical implementation of the MSA_Tree_metrics_exploration.R application. The platform provides three main analytical capabilities through a tabbed interface, each with specialized tools and visualizations. For detailed information about specific components, see User Interface Components, Data Loading and Filtering, Statistical Analysis and Visualization, and MSA Simulation Tools. For information about preparing input data, see Feature Extraction Tool. The built-in DNA and AA empirical datasets collected relevant metrics for 733,430 and 232,115 MSAs, respectively. Most features have been generated from the above extract_msa_tree_features.py, except for Gamma shape, which could be mined from tree inference tools, e.g. FastTree or IQ-TREE.

The Shiny application requires these R packages, which can be installed by:

&emsp;&emsp;install.packages(c("readr", "ggplot2", "dplyr", "shiny", "bslib", "corrplot"))

Launch the Shiny Application:

&emsp;&emsp;\# Local activation in R console or RStudio

&emsp;&emsp;shiny::runApp("MSA_Tree_metrics_exploration.R")

&emsp;&emsp;\# Or Online Shiny server

&emsp;&emsp;http://124.222.255.151:3838/apps/phylogenomics/


#### The Distribution Analysis tab enables detailed exploration of individual variables with advanced filtering and visualization options:

The distribution analysis tab provides histogram visualization with density overlays for exploring the statistical properties of individual variables. The main visualization generates histograms with density overlays. The plotting process handles filtered data and applies user-specified axis ranges. Users can adjust visualization parameters through UI controls in the sidebar.

&emsp;&emsp;Variable categorization: Groups variables into sequence-related, model-related, and gene tree-related categories.

&emsp;&emsp;Tukey's Fences filtering: Implements function for statistical outlier detection.

&emsp;&emsp;Range-based exclusion: Allows filtering samples based on specified value ranges.

&emsp;&emsp;Interactive histograms: Generated by with customizable binning and axis ranges.


<img width="2472" height="1166" alt="image" src="https://github.com/user-attachments/assets/23229ab0-9bab-4551-87ce-d8fd91a62807" />


#### Correlation Analysis Tab enables multi-variable correlation analysis with interactive heatmap visualization:

This function is enabled by the 'Enable Correlation Analysis' checkbox in the main sidebar. The correlation analysis uses Spearman's rank correlation and provides hierarchical clustering visualization through the *corrplot* package. 

<img width="2310" height="1094" alt="image" src="https://github.com/user-attachments/assets/f40b8034-bc66-4ff7-b163-6328ef9cd0de" />

#### MSA Simulation Tab generates IQ-TREE simulation commands using various parameter sampling strategies:

Complete empirical match: ensures parameter consistency by selecting all model parameters from the same empirical observation. This approach maintains natural covariances between parameters.

Mixed empirical sampling: samples each parameter independently from empirical distributions, allowing exploration of parameter combinations not present in the original data.

PDF-based estimation: creates continuous distributions from empirical data for selected parameters. 

The simulation tools allow researchers to generate IQ-TREE commands for phylogenomic data simulation using empirical parameter distributions extracted from real datasets. It implements three distinct parameter sampling strategies for generating realistic phylogenomic simulations. Each strategy represents a different approach to leveraging empirical data for simulation parameter selection. The simulation interface is organized into collapsible cards within the MSA Simulation tab, enabled by the 'Enable MSA Simulation' checkbox in the main sidebar. The interface adapts model component options based on the sequence type. For DNA sequences using the GTR model, the system implements flexible GTR rate parameter sourcing. PDF Comparison Plots generate side-by-side comparisons of empirical and simulated distributions.

<img width="2502" height="1080" alt="image" src="https://github.com/user-attachments/assets/0ee6d968-ff07-4b08-a178-23bbafbc1f3f" />

<img width="1562" height="614" alt="image" src="https://github.com/user-attachments/assets/0331004e-aaf8-4ce5-9c5f-2041eab70c7a" />




