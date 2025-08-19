# MSA-and-tree-metrics-exploration
MSA and Tree metrics exploration and MSA simulation of empirical phylogenomic datasets

This document describes the workflow of the MSA-and-tree-metrics-exploration system, detailing how the Python feature extraction tool and R Shiny application work together as an integrated bioinformatics pipeline. The system implements a two-tier architecture where batch processing and interactive analysis are separated into distinct but complementary components. The Python tool handles computationally intensive feature extraction, while the R Shiny application provides user-friendly data exploration and simulation planning capabilities. The two main components communicate through CSV files that serve as the data interchange format. The Python tool produces structured CSV output that the R application can directly consume.

## Feature Extraction Tool 
The Feature Extraction Tool is a Python command-line application that processes multiple sequence alignments (MSAs) and phylogenetic tree files to extract comprehensive statistical and evolutionary features for phylogenomic analysis. This tool serves as the batch processing component of the system, generating standardized feature datasets that feed into the interactive analysis platform. Required Python libraries are numpy, tqdm, and Biopython, which can be installed simply by 'pip install numpy tqdm biopython'.

The script extracts two main categories of features: MSA alignments statistics and tree-based branch and distance metrics.  MSA features include character ('ACGT' for nucleotides or 'ARNDCQEGHILKMFPSTWYV' for amino acids) frequencies, and 16 other features, such as number of sequences, alignment length, unique site patterns, parsimony-informative sites, gap proportion, average site entropy, proportion of invariant sites, relative composition frequency variation, average pairwise identity, Bollback multinomial etc. Fourteen tree-based features include bralch-related metrics, average bipartition support, Distance Variance from Molecular Clock (DVMC), evolutionary rate, treeness, and normalised Robinson-Foulds distance, and so on.

The extract_msa_tree_features.py script provides a unified command-line interface for extracting features from biological data files. The tool supports three primary modes of operation: MSA-only processing, tree-only processing, and combined MSA+tree processing with file matching by prefix. At least one of --msa_dir or --tree_dir must be provided. 

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
  
  python extract_msa_tree.py --msa_dir /path/to/MSAs --msa_type DNA -o results.csv

  \# Tree features only 
  
  python extract_msa_tree.py --tree_dir /path/to/trees -o results.csv

  \# Both MSA and tree features
  
  python extract_msa_tree.py --msa_dir /path/to/MSAs --tree_dir /path/to/trees --msa_type DNA -o results.csv

  \# Skip frequency statistics
  
  python extract_msa_tree.py --msa_dir /path/to/MSAs --tree_dir /path/to/trees --msa_type DNA --skip_freq_statistics -o results.csv

  \# With pseudo tree metrics
  
  python extract_msa_tree.py --msa_dir /path/to/MSAs --msa_type DNA --pseudo_tree_metrics --fasttree /usr/bin/FastTree -o results.csv

  \# With outgroups and reference tree
  
  python extract_msa_tree.py --msa_dir /path/to/MSAs --tree_dir /path/to/trees --msa_type DNA --outgroup_list outgroups.txt --ref_tree reference.tree -o results.csv



  ## Interactive Analysis Platform

  
