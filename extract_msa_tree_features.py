#!/usr/bin/env python3
"""
Unified MSA and Tree Feature Extractor
Combines functionality for MSA frequency analysis, MSA feature extraction, and tree metrics calculation
"""

import argparse
import os
import csv
import sys
import math
import statistics
import warnings
import time
import logging
import threading
from pathlib import Path
from datetime import datetime
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from tempfile import NamedTemporaryFile
import subprocess
from io import StringIO
import copy

import numpy as np
from tqdm import tqdm
from Bio import AlignIO, SeqIO, Phylo
from Bio.Align import MultipleSeqAlignment

# Define valid characters
DNA_BASES = 'ACGT'
AMINO_ACIDS = 'ARNDCQEGHILKMFPSTWYV'
GAP_CHARS = '-'

def setup_logger():
    """Setup logger for the application"""
    log_file = f"features_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    logger = logging.getLogger('UnifiedFeatures')
    logger.setLevel(logging.INFO)
    
    # File handler for all messages
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    
    # Console handler for key messages only
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    class InfoFilter(logging.Filter):
        def filter(self, record):
            if record.levelno >= logging.WARNING:
                return True
            if record.levelno == logging.INFO:
                main_keywords = [
                    "Starting feature extraction",
                    "Command line arguments",
                    "Found",
                    "File validation completed",
                    "Processing",
                    "Features saved to",
                    "Processing Summary",
                    "Total files",
                    "Successfully processed",
                    "Failed",
                    "Total runtime"
                ]
                return any(keyword in record.msg for keyword in main_keywords)
            return False
    
    ch.addFilter(InfoFilter())
    
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger

def format_float(value):
    """Format float value to 6 decimal places"""
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return f"{float(value):.6f}"
    return value

def get_file_prefix(filename):
    """Extract file prefix (name without extension)"""
    return Path(filename).stem

def validate_files(msa_dir, tree_dir, msa_type, logger):
    """Validate MSA and tree files for consistency"""
    msa_files = []
    tree_files = []
    
    # Get MSA files if directory provided
    if msa_dir:
        msa_path = Path(msa_dir)
        if not msa_path.exists():
            raise FileNotFoundError(f"MSA directory {msa_dir} does not exist")
        
        extensions = ['.fas', '.fasta', '.fa', '.fna', '.faa']
        for ext in extensions:
            msa_files.extend(msa_path.glob(f'*{ext}'))
            msa_files.extend(msa_path.glob(f'*{ext.upper()}'))
        
        if not msa_files:
            raise ValueError(f"No FASTA files found in {msa_dir}")
        
        # Validate FASTA format
        for msa_file in msa_files[:5]:  # Check first 5 files
            try:
                with open(msa_file) as f:
                    first_line = f.readline().strip()
                    if not first_line.startswith('>'):
                        raise ValueError(f"File {msa_file} does not appear to be in FASTA format")
            except Exception as e:
                raise ValueError(f"Error reading {msa_file}: {e}")
    
    # Get tree files if directory provided
    if tree_dir:
        tree_path = Path(tree_dir)
        if not tree_path.exists():
            raise FileNotFoundError(f"Tree directory {tree_dir} does not exist")
        
        tree_extensions = ['.tre', '.tree', '.nwk', '.newick', '.ph', '.treefile', '.bestTree']
        for ext in tree_extensions:
            tree_files.extend(tree_path.glob(f'*{ext}'))
            tree_files.extend(tree_path.glob(f'*{ext.upper()}'))
        
        if not tree_files:
            raise ValueError(f"No tree files found in {tree_dir}")
        
        # Validate Newick format
        for tree_file in tree_files[:5]:  # Check first 5 files
            try:
                Phylo.read(tree_file, 'newick')
            except Exception as e:
                raise ValueError(f"File {tree_file} does not appear to be in valid Newick format: {e}")
    
    # Check file prefix consistency if both directories provided
    if msa_dir and tree_dir:
        msa_prefixes = set(get_file_prefix(f.name) for f in msa_files)
        tree_prefixes = set(get_file_prefix(f.name) for f in tree_files)
        
        if msa_prefixes != tree_prefixes:
            missing_in_msa = tree_prefixes - msa_prefixes
            missing_in_tree = msa_prefixes - tree_prefixes
            
            error_msg = "File prefix mismatch between MSA and tree directories:\n"
            if missing_in_msa:
                error_msg += f"Missing in MSA directory: {sorted(missing_in_msa)}\n"
            if missing_in_tree:
                error_msg += f"Missing in tree directory: {sorted(missing_in_tree)}\n"
            
            raise ValueError(error_msg)
        
        if len(msa_files) != len(tree_files):
            raise ValueError(f"File count mismatch: {len(msa_files)} MSA files vs {len(tree_files)} tree files")
    
    logger.info(f"Found {len(msa_files)} MSA files and {len(tree_files)} tree files")
    return msa_files, tree_files

def check_fasttree():
    """Check if FastTree is available"""
    try:
        result = subprocess.run(['which', 'FastTree'], capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.strip()
    except Exception:
        pass
    return None

def standardize_sequence(sequence, msa_type):
    """Standardize sequence by replacing invalid characters with gaps"""
    sequence = sequence.upper()
    
    if msa_type == 'DNA':
        # Replace U with T first
        sequence = sequence.replace('U', 'T')
        valid_chars = DNA_BASES
    else:  # AA
        valid_chars = AMINO_ACIDS
    
    # Replace any character not in valid_chars with gap
    standardized = ''.join(c if c in valid_chars else '-' for c in sequence)
    return standardized

def calculate_frequencies(msa_file, msa_type):
    """Calculate character frequencies in MSA"""
    try:
        # Read all sequences and concatenate
        all_sequence = ""
        with open(msa_file) as f:
            for line in f:
                if not line.startswith('>'):
                    all_sequence += line.strip()
        
        # Standardize sequence
        standardized_seq = standardize_sequence(all_sequence, msa_type)
        
        # Get valid characters
        if msa_type == 'DNA':
            valid_chars = DNA_BASES
        else:
            valid_chars = AMINO_ACIDS
        
        # Only count valid characters (exclude gaps)
        valid_sequence = ''.join([c for c in standardized_seq if c in valid_chars])
        
        if not valid_sequence:
            return {f'freq{char}': format_float(0.0) for char in valid_chars}
        
        # Calculate frequencies
        counter = Counter(valid_sequence)
        total = len(valid_sequence)
        
        frequencies = {}
        for char in valid_chars:
            frequencies[f'freq{char}'] = format_float(counter.get(char, 0) / total)
        
        return frequencies
    except Exception as e:
        print(f"Error calculating frequencies for {msa_file}: {e}")
        return {f'freq{char}': format_float(0.0) for char in (DNA_BASES if msa_type == 'DNA' else AMINO_ACIDS)}

def calculate_msa_features(msa_file, msa_type, fasttree_path=None, pseudo_tree_metrics=False):
    """Calculate MSA features"""
    try:
        # Read MSA
        msa = AlignIO.read(msa_file, "fasta")
        
        # Standardize sequences
        for record in msa:
            standardized_seq = standardize_sequence(str(record.seq), msa_type)
            record.seq = record.seq.__class__(standardized_seq)
        
        # Basic metrics
        ntaxa = len(msa)
        nsites = msa.get_alignment_length()
        
        # Calculate patterns and other statistics
        sites = []
        gaps_count = 0
        invariant_count = 0
        parsimony_informative_count = 0
        
        valid_chars = DNA_BASES if msa_type == 'DNA' else AMINO_ACIDS
        
        for i in range(nsites):
            column = msa[:, i]
            sites.append(str(column))
            
            # Count gaps
            gaps_count += sum(1 for c in column if c == '-')
            
            # Count invariant sites
            unique_chars = set(c for c in column if c in valid_chars)
            if len(unique_chars) <= 1:
                invariant_count += 1
            
            # Count parsimony-informative sites
            char_counts = Counter(c for c in column if c in valid_chars)
            if len(char_counts) >= 2 and sum(count >= 2 for count in char_counts.values()) >= 2:
                parsimony_informative_count += 1
        
        patterns = len(set(sites))
        gaps_prop = gaps_count / (nsites * ntaxa)
        invariant_prop = invariant_count / nsites
        
        # Calculate entropy
        entropies = []
        for i in range(nsites):
            column = ''.join(c for c in msa[:, i] if c in valid_chars)
            if column:
                entropy = 0
                for char in valid_chars:
                    count = column.count(char)
                    if count > 0:
                        prob = count / len(column)
                        entropy -= prob * math.log2(prob)
                entropies.append(entropy)
        
        avg_entropy = statistics.mean(entropies) if entropies else 0
        
        # Pattern entropy
        site_counts = Counter(sites)
        pattern_entropy = sum(N_i * math.log(N_i) for N_i in site_counts.values())
        
        # Bollback multinomial
        bollback = pattern_entropy - nsites * math.log(nsites)
        
        # RCFV calculation
        rcfv = 0
        for state in valid_chars:
            state_freqs = []
            for seq in msa:
                seq_str = str(seq.seq)
                total = sum(seq_str.count(s) for s in valid_chars)
                if total > 0:
                    freq = seq_str.count(state) / total
                    state_freqs.append(freq)
            
            if state_freqs:
                mean_freq = sum(state_freqs) / len(state_freqs)
                rcfv += sum(abs(freq - mean_freq) for freq in state_freqs) / ntaxa
        
        # Average pairwise identity
        total_identity = 0
        pairs = 0
        for i in range(ntaxa):
            for j in range(i + 1, ntaxa):
                seq1 = str(msa[i].seq)
                seq2 = str(msa[j].seq)
                identical = sum(1 for a, b in zip(seq1, seq2) if a == b)
                total_identity += identical / len(seq1)
                pairs += 1
        
        avg_pairwise_identity = total_identity / pairs if pairs > 0 else 0
        
        features = {
            "num_taxa": ntaxa,
            "num_sites": nsites,
            "num_patterns": patterns,
            "num_parsimony_sites": parsimony_informative_count,
            "num_sites/num_taxa": format_float(nsites / ntaxa),
            "num_patterns/num_taxa": format_float(patterns / ntaxa),
            "num_parsimony_sites/num_taxa": format_float(parsimony_informative_count / ntaxa),
            "num_patterns/num_sites": format_float(patterns / nsites),
            "num_parsimony_sites/num_sites": format_float(parsimony_informative_count / nsites),
            "proportion_gaps": format_float(gaps_prop),
            "proportion_invariant": format_float(invariant_prop),
            "entropy": format_float(avg_entropy),
            "bollback": format_float(bollback),
            "pattern_entropy": format_float(pattern_entropy),
            "rcfv": format_float(rcfv),
            "average_pairwise_identity": format_float(avg_pairwise_identity)
        }
        
        # Add pseudo tree metrics if requested
        if pseudo_tree_metrics and fasttree_path:
            tree_metrics = calculate_pseudo_tree_metrics(msa_file, msa_type, fasttree_path)
            if tree_metrics:
                features.update(tree_metrics)
        
        return features
    except Exception as e:
        print(f"Error calculating MSA features for {msa_file}: {e}")
        return None

def calculate_pseudo_tree_metrics(msa_file, msa_type, fasttree_path):
    """Calculate tree metrics using FastTree"""
    try:
        # Prepare FastTree command
        if msa_type == 'DNA':
            cmd = [fasttree_path, '-nt', '-gtr', '-noml', '-boot', '500']
        else:
            cmd = [fasttree_path, '-lg', '-noml', '-boot', '500']
        
        # Run FastTree
        with open(msa_file) as f:
            msa_content = f.read()
        
        process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate(input=msa_content)
        
        if process.returncode != 0:
            print(f"FastTree failed for {msa_file}: {stderr}")
            return None
        
        # Parse tree and calculate metrics
        tree = Phylo.read(StringIO(stdout), 'newick')
        
        bootstrap_values = []
        internal_branch_lengths = []
        terminal_branch_lengths = []
        total_length = 0
        patristic_distances = []
        
        for clade in tree.find_clades():
            if clade.branch_length is not None:
                total_length += clade.branch_length
                if clade.is_terminal():
                    terminal_branch_lengths.append(clade.branch_length)
                else:
                    internal_branch_lengths.append(clade.branch_length)
                    if hasattr(clade, 'confidence') and clade.confidence is not None:
                        bootstrap_values.append(clade.confidence)
        
        # Calculate patristic distances
        terminals = tree.get_terminals()
        for i, term1 in enumerate(terminals):
            for term2 in terminals[i+1:]:
                distance = tree.distance(term1, term2)
                patristic_distances.append(distance)
        
        return {
            'average_BS_FT': format_float(np.mean(bootstrap_values)) if bootstrap_values else "NA",
            'sd_BS_FT': format_float(np.std(bootstrap_values)) if bootstrap_values else "NA",
            'total_tree_length_FT': format_float(total_length),
            'average_internal_branch_length_FT': format_float(np.mean(internal_branch_lengths)),
            'sd_internal_branch_length_FT': format_float(np.std(internal_branch_lengths)),
            'average_terminal_branch_length_FT': format_float(np.mean(terminal_branch_lengths)),
            'sd_terminal_branch_length_FT': format_float(np.std(terminal_branch_lengths)),
            'tree_diameter_FT': format_float(max(patristic_distances)),
            'average_patristic_distance_FT': format_float(np.mean(patristic_distances)),
            'sd_patristic_distance_FT': format_float(np.std(patristic_distances))
        }
    except Exception as e:
        print(f"Error calculating pseudo tree metrics for {msa_file}: {e}")
        return None

def calculate_dvmc(tree, outgroup_list=None):
    """Calculate DVMC after removing outgroups if specified"""
    try:
        tree_copy = copy.deepcopy(tree)
        
        # Remove outgroups if specified
        if outgroup_list:
            with open(outgroup_list) as f:
                outgroups = set(line.strip() for line in f)
            terminal_names = set(term.name for term in tree_copy.get_terminals())
            outgroups_in_tree = outgroups.intersection(terminal_names)
            if outgroups_in_tree:
                for tip in list(tree_copy.get_terminals()):
                    if tip.name in outgroups_in_tree:
                        tree_copy.prune(tip)
        
        # Calculate DVMC
        num_spp = tree_copy.count_terminals()
        sum_dist = 0
        sumi2N = 0
        
        for term in tree_copy.get_terminals():
            dist = tree_copy.distance(term)
            sum_dist += dist
            sumi2N += dist ** 2
        
        avg_dist = sum_dist / num_spp
        squared_diff_sum = sumi2N - num_spp * (avg_dist ** 2)
        
        return math.sqrt(squared_diff_sum / (num_spp - 1))
    except Exception as e:
        print(f"Error calculating DVMC: {e}")
        return None

def calculate_rf_distance(tree, ref_tree):
    """Calculate Robinson-Foulds distance"""
    try:
        tree_zero = copy.deepcopy(tree)
        tree_one = copy.deepcopy(ref_tree)
        
        # Get shared tips
        tree_zero_tips = [term.name for term in tree_zero.get_terminals()]
        tree_one_tips = [term.name for term in tree_one.get_terminals()]
        shared_tips = set(tree_zero_tips) & set(tree_one_tips)
        
        # Prune tips not in shared set
        tips_to_prune_zero = set(tree_zero_tips) - shared_tips
        tips_to_prune_one = set(tree_one_tips) - shared_tips
        
        for tip_name in tips_to_prune_zero:
            for tip in tree_zero.get_terminals():
                if tip.name == tip_name:
                    tree_zero.prune(tip)
                    break
        
        for tip_name in tips_to_prune_one:
            for tip in tree_one.get_terminals():
                if tip.name == tip_name:
                    tree_one.prune(tip)
                    break
        
        # Root both trees with the same tip
        for term in tree_zero.get_terminals():
            tip_for_rooting = term.name
            break
        tree_zero.root_with_outgroup(tip_for_rooting)
        tree_one.root_with_outgroup(tip_for_rooting)
        
        def compare_trees(plain_rf, tree_a, tree_b):
            for clade_a in tree_a.get_nonterminals()[1:]:
                tips_a = set(term.name for term in clade_a.get_terminals())
                clade_b = tree_b.common_ancestor([name for name in tips_a])
                tips_b = set(term.name for term in clade_b.get_terminals())
                if tips_a != tips_b:
                    plain_rf += 1
            return plain_rf
        
        plain_rf = 0
        plain_rf = compare_trees(plain_rf, tree_zero, tree_one)
        plain_rf = compare_trees(plain_rf, tree_one, tree_zero)
        
        tip_count = tree_zero.count_terminals()
        normalized_rf = plain_rf / (2 * (tip_count - 3))
        
        return normalized_rf
    except Exception as e:
        print(f"Error calculating RF distance: {e}")
        return None

def calculate_tree_features(tree_file, outgroup_list=None, ref_tree=None):
    """Calculate tree features"""
    try:
        tree = Phylo.read(tree_file, 'newick')
        
        # Basic tree metrics
        bootstrap_values = []
        internal_branch_lengths = []
        terminal_branch_lengths = []
        total_length = 0
        patristic_distances = []
        
        for clade in tree.find_clades():
            if clade.branch_length is not None:
                total_length += clade.branch_length
                if clade.is_terminal():
                    terminal_branch_lengths.append(clade.branch_length)
                else:
                    internal_branch_lengths.append(clade.branch_length)
                    if hasattr(clade, 'confidence') and clade.confidence is not None:
                        bootstrap_values.append(clade.confidence)
        
        # Calculate patristic distances
        terminals = tree.get_terminals()
        for i, term1 in enumerate(terminals):
            for term2 in terminals[i+1:]:
                distance = tree.distance(term1, term2)
                patristic_distances.append(distance)
        
        # Calculate metrics
        avg_bs = np.mean(bootstrap_values) if bootstrap_values else None
        sd_bs = np.std(bootstrap_values) if bootstrap_values else None
        avg_internal_branch_length = np.mean(internal_branch_lengths)
        sd_internal_branch_length = np.std(internal_branch_lengths)
        avg_terminal_branch_length = np.mean(terminal_branch_lengths)
        sd_terminal_branch_length = np.std(terminal_branch_lengths)
        tree_diameter = max(patristic_distances)
        avg_patristic_distance = np.mean(patristic_distances)
        sd_patristic_distance = np.std(patristic_distances)
        
        # Additional metrics
        evo_rate = total_length / len(terminals)
        treeness = sum(internal_branch_lengths) / total_length if total_length > 0 else 0
        
        # DVMC
        dvmc = calculate_dvmc(tree, outgroup_list)
        
        # RF distance
        rf_dist = None
        if ref_tree:
            rf_dist = calculate_rf_distance(tree, ref_tree)
        
        features = {
            'average_BS': format_float(avg_bs) if avg_bs is not None else "NA",
            'sd_BS': format_float(sd_bs) if sd_bs is not None else "NA",
            'total_tree_length': format_float(total_length),
            'average_internal_branch_length': format_float(avg_internal_branch_length),
            'sd_internal_branch_length': format_float(sd_internal_branch_length),
            'average_terminal_branch_length': format_float(avg_terminal_branch_length),
            'sd_terminal_branch_length': format_float(sd_terminal_branch_length),
            'tree_diameter': format_float(tree_diameter),
            'average_patristic_distance': format_float(avg_patristic_distance),
            'sd_patristic_distance': format_float(sd_patristic_distance),
            'evo_rate': format_float(evo_rate),
            'treeness': format_float(treeness),
            'dvmc': format_float(dvmc) if dvmc is not None else "NA"
        }
        
        if ref_tree:
            features['RF_distance'] = format_float(rf_dist) if rf_dist is not None else "NA"
        
        return features
    except Exception as e:
        print(f"Error calculating tree features for {tree_file}: {e}")
        return None

def process_file_pair(args):
    """Process a single MSA-tree file pair"""
    (loci, msa_file, tree_file, msa_type, skip_freq, pseudo_tree_metrics, 
     fasttree_path, outgroup_list, ref_tree) = args
    
    try:
        result = {"loci": loci}
        
        # Add DataType only if MSA processing is involved
        if msa_file:
            result["DataType"] = msa_type
        
        # Calculate frequency statistics
        if msa_file and not skip_freq:
            freq_features = calculate_frequencies(msa_file, msa_type)
            result.update(freq_features)
        
        # Calculate MSA features
        if msa_file:
            msa_features = calculate_msa_features(msa_file, msa_type, fasttree_path, pseudo_tree_metrics)
            if msa_features:
                result.update(msa_features)
        
        # Calculate tree features
        if tree_file:
            tree_features = calculate_tree_features(tree_file, outgroup_list, ref_tree)
            if tree_features:
                result.update(tree_features)
        
        return result
    except Exception as e:
        print(f"Error processing {loci}: {e}")
        return None

def create_parser():
    """Create argument parser with examples"""
    parser = argparse.ArgumentParser(
        description="Unified MSA and Tree Feature Extractor",
        epilog="""
Examples:
  # MSA features only
  python extract_msa_tree_features.py --msa_dir /path/to/MSAs --msa_type DNA -o results.csv

  # Tree features only  
  python extract_msa_tree_features.py --tree_dir /path/to/trees -o results.csv

  # Both MSA and tree features
  python extract_msa_tree_features.py --msa_dir /path/to/MSAs --tree_dir /path/to/trees --msa_type DNA -o results.csv

  # Skip frequency statistics
  python extract_msa_tree_features.py --msa_dir /path/to/MSAs --tree_dir /path/to/trees --msa_type DNA --skip_freq_statistics -o results.csv

  # With pseudo tree metrics
  python extract_msa_tree_features.py --msa_dir /path/to/MSAs --msa_type DNA --pseudo_tree_metrics --fasttree /usr/bin/FastTree -o results.csv

  # With outgroups and reference tree
  python extract_msa_tree_features.py --msa_dir /path/to/MSAs --tree_dir /path/to/trees --msa_type DNA --outgroup_list outgroups.txt --ref_tree reference.tree -o results.csv
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("--msa_dir", help="Directory containing MSA files (FASTA)")
    parser.add_argument("--tree_dir", help="Directory containing tree files (Newick)")
    parser.add_argument("--msa_type", choices=['AA', 'DNA'], 
                       help="MSA type: AA (amino acid) or DNA (nucleotide). Required when --msa_dir is provided")
    parser.add_argument("-o", "--output", help="Output CSV file")
    parser.add_argument("-t", "--threads", type=int, default=1, 
                       help="Number of threads for parallel processing")
    parser.add_argument("--skip_freq_statistics", action="store_true",
                       help="Skip frequency statistics calculation")
    parser.add_argument("--pseudo_tree_metrics", action="store_true",
                       help="Calculate pseudo tree metrics using faster FastTree execution (-noml -boot 500)")
    parser.add_argument("--fasttree", help="Path to FastTree executable")
    parser.add_argument("--outgroup_list", help="File containing outgroup species for DVMC calculation")
    parser.add_argument("--ref_tree", help="Reference tree for RF distance calculation")
    
    return parser

def main():
    start_time = time.time()
    
    parser = create_parser()
    args = parser.parse_args()
    
    # Validate arguments
    if not args.msa_dir and not args.tree_dir:
        parser.error("At least one of --msa_dir or --tree_dir must be provided")
    
    if args.msa_dir and not args.msa_type:
        parser.error("--msa_type is required when --msa_dir is provided")
    
    logger = setup_logger()
    logger.info("Starting feature extraction")
    logger.info(f"Command line arguments: {vars(args)}")
    
    # Validate files
    try:
        msa_files, tree_files = validate_files(args.msa_dir, args.tree_dir, args.msa_type, logger)
        logger.info("File validation completed successfully")
    except Exception as e:
        logger.error(f"File validation failed: {e}")
        return
    
    # Check FastTree if needed
    fasttree_path = None
    if args.pseudo_tree_metrics:
        fasttree_path = args.fasttree or check_fasttree()
        if not fasttree_path:
            logger.error("FastTree not found. Please provide --fasttree path or install FastTree")
            return
        logger.info(f"FastTree found at: {fasttree_path}")
    
    # Read reference tree if provided
    ref_tree = None
    if args.ref_tree:
        try:
            ref_tree = Phylo.read(args.ref_tree, 'newick')
            logger.info(f"Reference tree loaded: {args.ref_tree}")
        except Exception as e:
            logger.error(f"Error reading reference tree: {e}")
            return
    
    # Create file pairs for processing
    msa_dict = {get_file_prefix(f.name): f for f in msa_files} if msa_files else {}
    tree_dict = {get_file_prefix(f.name): f for f in tree_files} if tree_files else {}
    
    all_loci = set(msa_dict.keys()) | set(tree_dict.keys())
    
    process_args = []
    for loci in sorted(all_loci):
        msa_file = msa_dict.get(loci)
        tree_file = tree_dict.get(loci)
        
        process_args.append((
            loci, msa_file, tree_file, args.msa_type, args.skip_freq_statistics,
            args.pseudo_tree_metrics, fasttree_path, args.outgroup_list, ref_tree
        ))
    
    logger.info(f"Processing {len(process_args)} loci using {args.threads} threads")
    
    # Process files with progress bar
    all_features = []
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        results = list(tqdm(
            executor.map(process_file_pair, process_args),
            total=len(process_args),
            desc="Processing loci",
            unit="loci"
        ))
        all_features.extend([r for r in results if r is not None])
    
    if not all_features:
        logger.error("No features were successfully extracted")
        return
    
    # Write results to CSV
    output = args.output or "features.csv"
    
    # Determine field order
    fieldnames = ["loci"]
    
    # Add DataType if MSA processing was involved
    if args.msa_dir:
        fieldnames.append("DataType")
    
    # Add frequency fields if not skipped
    if not args.skip_freq_statistics and args.msa_dir:
        if args.msa_type == 'DNA':
            fieldnames.extend([f'freq{base}' for base in DNA_BASES])
        else:
            fieldnames.extend([f'freq{aa}' for aa in AMINO_ACIDS])
    
    # Add other fields from first result
    if all_features:
        for key in all_features[0].keys():
            if key not in fieldnames:
                fieldnames.append(key)
    
    with open(output, "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_features)
    
    logger.info(f"Features saved to: {output}")
    
    # Print summary
    end_time = time.time()
    total_time = end_time - start_time
    hours = int(total_time // 3600)
    minutes = int((total_time % 3600) // 60)
    seconds = int(total_time % 60)
    
    logger.info(f"\nProcessing Summary:")
    logger.info(f"Total files found: {len(process_args)}")
    logger.info(f"Successfully processed: {len(all_features)}")
    logger.info(f"Failed: {len(process_args) - len(all_features)}")
    logger.info(f"Total runtime: {hours:02d}:{minutes:02d}:{seconds:02d}")
    print(f"\nIdentical sequences are removed during the analysis, so the number of species may not be consistent with the number of sequences in the MSA alignments.")

if __name__ == "__main__":
    main()
