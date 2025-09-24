#!/usr/bin/env python3
"""
Python implementation of BAliBASE scoring algorithm for MSA evaluation.

This script implements the same logic as the original BAliBASE R9 bali_score.c
but accepts FASTA format inputs and handles sequence order differences between
reference and test alignments.

Usage:
    # Single-file mode
    python bali_score_single.py ref_alignment.fasta test_alignment.fasta [-v]

    # Batch mode (auto-discover references by filename)
    python bali_score_batch.py --data-path ./data/input \
                                --pred-path ./outputs \
                                --ref-path ./data/v2/aln \
                                [-v]
                                
    python ./../bali_score_batch.py --data-path ./v2/no_aln/ --pred-path ./base/ltf2_noflags/ --ref-path ./v2/aln/
    python ./../bali_score_batch.py --data-path ./v2/no_aln/ --pred-path ./base/ltf2_res/ --ref-path ./v2/aln/
    python ./../bali_score_batch.py --data-path ./v2/no_aln/ --pred-path ./base/ltf2_lazy/ --ref-path ./v2/aln/
    python ./../bali_score_batch.py --data-path ./v2/no_aln/ --pred-path ./base/ltf2_lazy_res/ --ref-path ./v2/aln/

    python ./../bali_score_batch.py --data-path ./v2/no_aln/ --pred-path ./no0.1threshold/ltf2_noflags/ --ref-path ./v2/aln/
    python ./../bali_score_batch.py --data-path ./v2/no_aln/ --pred-path ./no0.1threshold/ltf2_res/ --ref-path ./v2/aln/
    python ./../bali_score_batch.py --data-path ./v2/no_aln/ --pred-path ./no0.1threshold/ltf2_lazy/ --ref-path ./v2/aln/
    python ./../bali_score_batch.py --data-path ./v2/no_aln/ --pred-path ./no0.1threshold/ltf2_lazy_res/ --ref-path ./v2/aln/


Scores:
    SP score: Sum-of-Pairs score (0.0 to 1.0)
    TC score: Total Column score (0.0 to 1.0)
"""

import sys
import argparse
import os
from typing import Dict, List, Tuple, Optional
from collections import defaultdict


class Sequence:
    """Represents a sequence with name and aligned data."""
    
    def __init__(self, name: str, data: str):
        norm_name = name.strip()
        if " (k_id:" in norm_name:
            norm_name = norm_name.split(" (k_id:", 1)[0].rstrip()
        self.name = norm_name
        self.data = data.upper()
        self.length = len(data)


class Alignment:
    """Represents a multiple sequence alignment."""
    
    def __init__(self):
        self.sequences: List[Sequence] = []
        self.nseqs = 0
        self.maxlen = 0
    
    def add_sequence(self, name: str, data: str):
        """Add a sequence to the alignment."""
        seq = Sequence(name, data)
        self.sequences.append(seq)
        self.nseqs += 1
        self.maxlen = max(self.maxlen, seq.length)
    
    def get_sequence_by_name(self, name: str) -> Optional[Sequence]:
        """Find sequence by name (case-insensitive)."""
        for seq in self.sequences:
            if seq.name.lower() == name.lower():
                return seq
        return None


def read_fasta(filename: str) -> Alignment:
    """Read sequences from a FASTA file.
    Behavior:
      - If the file contains a line starting with the marker '[Multi Sequence Alignment]',
        parse FASTA only after that marker and stop if another '['-started marker appears.
      - Otherwise, parse the entire file as standard FASTA.
    """
    alignment = Alignment()
    current_name = None
    current_seq = []
    
    try:
        with open(filename, 'r') as f:
            lines = [ln.rstrip('\n') for ln in f]

        # Find canonical marker only
        marker_idx = -1
        for idx, ln in enumerate(lines):
            s = ln.strip()
            if s.startswith('[Multi Sequence Alignment]'):
                marker_idx = idx
                break

        # Decide range to parse
        parse_lines = lines[marker_idx + 1:] if marker_idx >= 0 else lines

        # Parse until another marker line (starting with '[')
        for raw_line in parse_lines:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith('['):
                break
            if line.startswith('>'):
                if current_name is not None:
                    alignment.add_sequence(current_name, ''.join(current_seq))
                current_name = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)

        if current_name is not None:
            alignment.add_sequence(current_name, ''.join(current_seq))

        # Fallback: only if no canonical marker was present and nothing parsed
        if alignment.nseqs == 0 and marker_idx == -1:
            current_name = None
            current_seq = []
            for raw_line in lines:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_name is not None:
                        alignment.add_sequence(current_name, ''.join(current_seq))
                    current_name = line[1:].strip()
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_name is not None:
                alignment.add_sequence(current_name, ''.join(current_seq))

    except FileNotFoundError:
        print(f"Error: Cannot open file '{filename}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file '{filename}': {e}", file=sys.stderr)
        sys.exit(1)
    
    return alignment


def create_sequence_mapping(ref_aln: Alignment, test_aln: Alignment) -> Dict[int, int]:
    """
    Create mapping from test sequence indices to reference sequence indices.
    Returns dict: test_seq_index -> ref_seq_index
    """
    seq_mapping = {}
    
    for test_idx, test_seq in enumerate(test_aln.sequences):
        ref_seq = ref_aln.get_sequence_by_name(test_seq.name)
        if ref_seq is None:
            print(f"Error: sequence '{test_seq.name}' not found in reference alignment", 
                  file=sys.stderr)
            sys.exit(1)
        
        # Find reference sequence index
        ref_idx = -1
        for i, seq in enumerate(ref_aln.sequences):
            if seq.name.lower() == test_seq.name.lower():
                ref_idx = i
                break
        
        seq_mapping[test_idx] = ref_idx
    
    return seq_mapping


def find_valid_columns(ref_aln: Alignment, gap_cutoff_percent: float = 20.0) -> List[bool]:
    """
    Identify columns in reference alignment with fewer gaps than cutoff.
    Returns list of booleans indicating valid columns.
    """
    valid_cols = []
    gap_cutoff = max(1, int(ref_aln.nseqs * gap_cutoff_percent / 100.0))
    
    for col_idx in range(ref_aln.maxlen):
        gap_count = 0
        
        for seq in ref_aln.sequences:
            if col_idx < seq.length and seq.data[col_idx] == '-':
                gap_count += 1
        
        # Column is valid if gap count is below cutoff
        valid_cols.append(gap_count < gap_cutoff)
    
    return valid_cols


def code_reference_alignment(ref_aln: Alignment, valid_cols: List[bool]) -> List[List[int]]:
    """
    Code reference alignment: assign column numbers to each residue position.
    Gap positions are coded as 0.
    """
    ref_codes = []
    
    for seq in ref_aln.sequences:
        seq_codes = []
        for col_idx in range(ref_aln.maxlen):
            if col_idx >= seq.length or seq.data[col_idx] == '-':
                seq_codes.append(0)  # Gap
            elif valid_cols[col_idx]:
                seq_codes.append(col_idx + 1)  # Column number (1-based)
            else:
                seq_codes.append(0)  # Invalid column
        ref_codes.append(seq_codes)
    
    return ref_codes


def code_test_alignment(ref_aln: Alignment, test_aln: Alignment, 
                       seq_mapping: Dict[int, int], ref_codes: List[List[int]],
                       valid_cols: List[bool]) -> List[List[int]]:
    """
    Code test alignment by mapping each residue to reference column numbers.
    """
    test_codes = []
    
    for test_idx, test_seq in enumerate(test_aln.sequences):
        ref_idx = seq_mapping[test_idx]
        seq_codes = [0] * test_aln.maxlen
        
        # Track position in reference sequence
        ref_pos = 0
        
        for test_pos in range(test_seq.length):
            if test_seq.data[test_pos] == '-':
                seq_codes[test_pos] = 0  # Gap in test
            else:
                # Find next non-gap position in reference
                while (ref_pos < ref_aln.maxlen and 
                       (ref_pos >= ref_aln.sequences[ref_idx].length or 
                        ref_aln.sequences[ref_idx].data[ref_pos] == '-')):
                    ref_pos += 1
                
                if ref_pos < ref_aln.maxlen and valid_cols[ref_pos]:
                    seq_codes[test_pos] = ref_codes[ref_idx][ref_pos]
                else:
                    seq_codes[test_pos] = 0
                
                ref_pos += 1
        
        test_codes.append(seq_codes)
    
    return test_codes


def calculate_max_score(valid_cols: List[bool], ref_aln: Alignment) -> float:
    """Calculate maximum possible SP score for the reference alignment."""
    max_score = 0.0
    
    for col_idx, is_valid in enumerate(valid_cols):
        if is_valid:
            # Count non-gap residues in this column
            non_gap_count = 0
            for seq in ref_aln.sequences:
                if col_idx < seq.length and seq.data[col_idx] != '-':
                    non_gap_count += 1
            
            if non_gap_count > 1:
                # Add pairs: n*(n-1)/2
                max_score += non_gap_count * (non_gap_count - 1) / 2.0
    
    return max_score


def calculate_scores(test_aln: Alignment, test_codes: List[List[int]], 
                    valid_cols: List[bool], ref_aln: Alignment, 
                    verbose: bool = False) -> Tuple[float, float, float]:
    """
    Calculate SP and TC scores for the test alignment.
    Returns (sp_score, tc_score).
    """
    sp_score = 0.0
    tc_score = 0
    total_columns = 0
    predicted_pairs_total = 0.0
    
    if verbose:
        print("\nDetailed scoring:")
        print("Column analysis (1=correctly aligned, 0=incorrect, .=gap):")
    
    # Process each column in test alignment
    for col_idx in range(test_aln.maxlen):
        # Group sequences by their reference column assignments
        col_groups = defaultdict(list)
        
        for seq_idx in range(test_aln.nseqs):
            if col_idx < len(test_codes[seq_idx]):
                ref_col = test_codes[seq_idx][col_idx]
                if ref_col > 0:  # Non-gap, valid reference column
                    col_groups[ref_col].append(seq_idx)
        
        # Calculate SP score for this column
        col_sp_score = 0.0
        for ref_col, seq_indices in col_groups.items():
            if len(seq_indices) > 1:
                # Limit group size to reference column size
                ref_col_idx = ref_col - 1
                max_group_size = 0
                if ref_col_idx < len(valid_cols) and valid_cols[ref_col_idx]:
                    # Count non-gaps in reference column
                    for seq in ref_aln.sequences:
                        if ref_col_idx < seq.length and seq.data[ref_col_idx] != '-':
                            max_group_size += 1
                
                actual_group_size = min(len(seq_indices), max_group_size)
                if actual_group_size > 1:
                    col_sp_score += actual_group_size * (actual_group_size - 1) / 2.0
                # Count predicted pairs (denominator for precision) using the actual predicted group size
                predicted_pairs_total += len(seq_indices) * (len(seq_indices) - 1) / 2.0
        
        sp_score += col_sp_score
        
        # Calculate TC score for this column
        col_tc_score = 0
        if col_groups:
            # Check if any group is complete (matches reference column size)
            for ref_col, seq_indices in col_groups.items():
                ref_col_idx = ref_col - 1
                if ref_col_idx < len(valid_cols) and valid_cols[ref_col_idx]:
                    expected_size = 0
                    for seq in ref_aln.sequences:
                        if ref_col_idx < seq.length and seq.data[ref_col_idx] != '-':
                            expected_size += 1
                    
                    if len(seq_indices) >= expected_size:
                        col_tc_score = 1
                        break
            
            total_columns += 1
        
        tc_score += col_tc_score
        
        # Verbose output
        if verbose and col_idx < 50:  # Limit verbose output
            if col_groups:
                print(f"{col_tc_score}", end="")
            else:
                print(".", end="")
            
            if (col_idx + 1) % 50 == 0:
                print()
    
    if verbose:
        print("\n")
    
    # Convert TC score to percentage
    if total_columns > 0:
        tc_score = tc_score / total_columns
    else:
        tc_score = 0.0
    
    return sp_score, tc_score, predicted_pairs_total


def main():
    """Main function to run BAliBASE scoring."""
    parser = argparse.ArgumentParser(
        description='Python implementation of BAliBASE scoring for MSA evaluation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python bali_score_python.py reference.fasta test.fasta
    python bali_score_python.py reference.fasta test.fasta -v
    python bali_score_python.py --data-path ./data/input --pred-path ./outputs --ref-path ./data/v2/aln
        """
    )
    
    # Single-file mode
    parser.add_argument('ref_alignment', nargs='?',
                       help='Reference alignment file in FASTA format')
    parser.add_argument('test_alignment', nargs='?',
                       help='Test alignment file in FASTA format')

    # Batch mode
    parser.add_argument('-d', '--data-path', type=str, default='',
                       help='Input folder (used for filename basis)')
    parser.add_argument('-p', '--pred-path', type=str, default='',
                       help='Predicted outputs folder')
    parser.add_argument('--ref-path', type=str,
                       default='/nfs/hpc/share/naukarkr/LiangHuang/LTF2/LinearX/eval/rnastralign/data/v2/aln',
                       help='Folder containing golden reference alignments')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output showing column-by-column scoring')
    
    args = parser.parse_args()

    def score_pair(ref_file: str, test_file: str) -> Optional[Tuple[float, float, float, float]]:
        if args.verbose:
            print(f"\nComparing test alignment in {test_file}")
            print(f"with reference alignment in {ref_file}")
        ref_aln = read_fasta(ref_file)
        test_aln = read_fasta(test_file)
        if ref_aln.nseqs != test_aln.nseqs:
            if args.verbose:
                print(f"Error: {ref_aln.nseqs} sequences in {ref_file} and {test_aln.nseqs} in {test_file}", file=sys.stderr)
            return None
        if args.verbose:
            print(f"\nFound {ref_aln.nseqs} sequences")
        seq_mapping = create_sequence_mapping(ref_aln, test_aln)
        valid_cols = find_valid_columns(ref_aln)
        valid_count = sum(valid_cols)
        if args.verbose:
            print(f"Using {valid_count} columns (gap cutoff: 20%)")
        ref_codes = code_reference_alignment(ref_aln, valid_cols)
        test_codes = code_test_alignment(ref_aln, test_aln, seq_mapping, ref_codes, valid_cols)
        max_score = calculate_max_score(valid_cols, ref_aln)
        if max_score <= 0:
            if args.verbose:
                print("Error: Invalid reference alignment (max score = 0)", file=sys.stderr)
            return None
        sp_raw, tc_score, predicted_pairs_total = calculate_scores(test_aln, test_codes, valid_cols, ref_aln, args.verbose)
        sp_score = sp_raw / max_score if max_score > 0 else 0.0
        # Precision approximation: correct predicted pairs / total predicted pairs
        precision_score = (sp_raw / predicted_pairs_total) if predicted_pairs_total > 0 else 0.0
        recall_score = sp_score  # SP equals recall against reference-coded columns
        f1_score = (2 * precision_score * recall_score / (precision_score + recall_score)) if (precision_score + recall_score) > 0 else 0.0
        if args.verbose:
            print(f"\n\tRecall (SP)= {recall_score:.3f}")
            print(f"\n\tPrecision= {precision_score:.3f}")
            print(f"\n\tF1= {f1_score:.3f}")
            print(f"\n\tTC score= {tc_score:.3f}")
            print(f"auto {test_file} recall={recall_score:.3f} precision={precision_score:.3f} f1={f1_score:.3f} tc={tc_score:.3f}")
        return recall_score, precision_score, f1_score, tc_score

    # Batch mode if pred-path provided (and optional data-path)
    if args.pred_path:
        ref_dir = args.ref_path
        pred_dir = args.pred_path
        data_dir = args.data_path or ''

        try:
            pred_files = sorted([f for f in os.listdir(pred_dir) if (f.endswith('.fasta') or f.endswith('.txt'))])
        except FileNotFoundError:
            print(f"Error: pred-path not found: {pred_dir}", file=sys.stderr)
            sys.exit(1)

        results = []
        families = set()
        family_precision = defaultdict(lambda: list())
        family_recall = defaultdict(lambda: list())
        family_f1 = defaultdict(lambda: list())
        family_tc = defaultdict(lambda: list())
        for fname in pred_files:
            ref_candidate = os.path.join(ref_dir, fname)
            if not os.path.isfile(ref_candidate) and fname.endswith('.fasta'):
                alt = fname[:-6] + '.aln.fasta'  # replace .fasta with .aln.fasta
                ref_candidate = os.path.join(ref_dir, alt)
            if not os.path.isfile(ref_candidate) and fname.endswith('.txt'):
                alt = fname[:-4] + '.aln.fasta'
                ref_candidate = os.path.join(ref_dir, alt)

            if not os.path.isfile(ref_candidate):
                print(f"Skipping {fname}: reference not found in {ref_dir}", file=sys.stderr)
                continue

            test_file = os.path.join(pred_dir, fname)
            scored = score_pair(ref_candidate, test_file)
            if scored is not None:
                results.append(scored)
                # Derive family name similar to eval_perf (prefix before first dot)
                family = fname.split('.')[0]
                families.add(family)
                recall_score, precision_score, f1_score, tc_score = scored
                family_precision[family].append(precision_score)
                family_recall[family].append(recall_score)
                family_f1[family].append(f1_score)
                family_tc[family].append(tc_score)

        if results:
            # Print summary table similar to eval_perf.py
            print()
            header = "{:<8}\t{:>10}\t{:>10}\t{:>10}\t{:>20}".format(
                "Family",
                "Precision",
                "Sensitivity",
                "F1-score (min, avg, max)",
                "TC score (min, avg, max)",
            )
            print(header)

            total_f1_avg = 0.0
            total_tc_avg = 0.0
            for family in sorted(families):
                prec_list = family_precision[family]
                rec_list = family_recall[family]
                f1_list = family_f1[family]
                tc_list = family_tc[family]

                prec_avg = sum(prec_list) / len(prec_list) if prec_list else 0.0
                rec_avg = sum(rec_list) / len(rec_list) if rec_list else 0.0
                f1_min = min(f1_list) if f1_list else 0.0
                f1_avg = sum(f1_list) / len(f1_list) if f1_list else 0.0
                f1_max = max(f1_list) if f1_list else 0.0
                tc_min = min(tc_list) if tc_list else 0.0
                tc_avg = sum(tc_list) / len(tc_list) if tc_list else 0.0
                tc_max = max(tc_list) if tc_list else 0.0

                print(
                    "{:<8}\t{:>6.2f}\t\t{:>6.2f}\t\t({:>4.2f}  {:>6.2f}  {:>6.2f})\t\t({:>6.2f}  {:>6.2f}  {:>6.2f})".format(
                        family,
                        prec_avg * 100,
                        rec_avg * 100,
                        f1_min * 100,
                        f1_avg * 100,
                        f1_max * 100,
                        tc_min,
                        tc_avg,
                        tc_max,
                    )
                )

                total_f1_avg += f1_avg
                total_tc_avg += tc_avg

            # Averages across families
            total_f1_avg /= len(families) if families else 1.0
            total_tc_avg /= len(families) if families else 1.0

            print("Average F1-score: %0.2f" % (total_f1_avg * 100))
            print("Average TC score: %0.2f\n" % (total_tc_avg))
        else:
            print("\nNo successful comparisons.\n")
        return

    # Single-file mode
    if not args.ref_alignment or not args.test_alignment:
        parser.error('Either provide two FASTA files (single-file mode) or --pred-path (batch mode).')

    scored = score_pair(args.ref_alignment, args.test_alignment)
    if scored is None:
        sys.exit(1)


if __name__ == "__main__":
    main()
