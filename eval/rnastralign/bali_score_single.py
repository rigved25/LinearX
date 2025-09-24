#!/usr/bin/env python3
"""
Python implementation of BAliBASE scoring algorithm for MSA evaluation.

This script implements the same logic as the original BAliBASE R9 bali_score.c
but accepts FASTA format inputs and handles sequence order differences between
reference and test alignments.

Usage:
    python bali_score_python.py ref_alignment.fasta test_alignment.fasta [-v]

    python bali_score_python.py ./data/v2/aln/5S.RNAStrAlign.k_30_1.aln.fasta ./ltf2_noflags/5S.RNAStrAlign.k_30_1.fasta -v

    python bali_score_python.py 1ref.fasta 2test.fasta -v
    python bali_score_python.py 2ref.fasta 3test.fasta -v
    python ./../bali_score_python.py ./../2ref.fasta ./../3test.fasta -v

Scores:
    SP score: Sum-of-Pairs score (0.0 to 1.0)
    TC score: Total Column score (0.0 to 1.0)
"""

import sys
import argparse
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

        # Find the canonical marker only
        marker_idx = -1
        for idx, ln in enumerate(lines):
            s = ln.strip()
            if s.startswith('[Multi Sequence Alignment]'):
                marker_idx = idx
                break

        # Decide range to parse
        parse_lines = lines[marker_idx + 1:] if marker_idx >= 0 else lines

        # Parse until another marker line (starting with '[') appears
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

        # Fallback: if no canonical marker was present and nothing parsed, read as plain FASTA
        if alignment.nseqs == 0 and marker_idx == -1:
            alignment = Alignment()
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
                    verbose: bool = False) -> Tuple[float, float]:
    """
    Calculate SP and TC scores for the test alignment.
    Returns (sp_score, tc_score).
    """
    sp_score = 0.0
    tc_score = 0
    total_columns = 0
    
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
    
    return sp_score, tc_score


def main():
    """Main function to run BAliBASE scoring."""
    parser = argparse.ArgumentParser(
        description='Python implementation of BAliBASE scoring for MSA evaluation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python bali_score_python.py reference.fasta test.fasta
    python bali_score_python.py reference.fasta test.fasta -v
        """
    )
    
    parser.add_argument('ref_alignment', 
                       help='Reference alignment file in FASTA format')
    parser.add_argument('test_alignment', 
                       help='Test alignment file in FASTA format')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output showing column-by-column scoring')
    
    args = parser.parse_args()
    
    # Read alignments
    print(f"\nComparing test alignment in {args.test_alignment}")
    print(f"with reference alignment in {args.ref_alignment}")
    
    ref_aln = read_fasta(args.ref_alignment)
    test_aln = read_fasta(args.test_alignment)
    
    # Check sequence counts match
    if ref_aln.nseqs != test_aln.nseqs:
        print(f"Error: {ref_aln.nseqs} sequences in {args.ref_alignment} "
              f"and {test_aln.nseqs} in {args.test_alignment}", file=sys.stderr)
        sys.exit(1)
    
    print(f"\nFound {ref_aln.nseqs} sequences")
    
    # Create sequence mapping
    seq_mapping = create_sequence_mapping(ref_aln, test_aln)
    
    # Find valid columns (with gap cutoff)
    valid_cols = find_valid_columns(ref_aln)
    valid_count = sum(valid_cols)
    print(f"Using {valid_count} columns (gap cutoff: 20%)")
    
    # Code alignments
    ref_codes = code_reference_alignment(ref_aln, valid_cols)
    test_codes = code_test_alignment(ref_aln, test_aln, seq_mapping, ref_codes, valid_cols)
    
    # Calculate maximum possible score
    max_score = calculate_max_score(valid_cols, ref_aln)
    if max_score <= 0:
        print("Error: Invalid reference alignment (max score = 0)", file=sys.stderr)
        sys.exit(1)
    
    # Calculate scores
    sp_raw, tc_score = calculate_scores(test_aln, test_codes, valid_cols, ref_aln, args.verbose)
    sp_score = sp_raw / max_score if max_score > 0 else 0.0
    
    # Output results
    print(f"\n\tSP score= {sp_score:.3f}")
    print(f"\n\tTC score= {tc_score:.3f}")
    print(f"auto {args.test_alignment} {sp_score:.3f} {tc_score:.3f}")


if __name__ == "__main__":
    main()
