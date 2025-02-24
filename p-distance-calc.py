import sys
from Bio import SeqIO, Align

def calculate_p_distance(seq1, seq2):
    # Create a PairwiseAligner object with appropriate scoring
    aligner = Align.PairwiseAligner(scoring="blastn")
    aligner.mode = "global"

    # Perform the alignment 
    alignments = aligner.align(seq1, seq2)

    alignment = next(alignments)

    # Extract aligned sequences as strings
    aligned_seq1 = str(alignment[0])
    aligned_seq2 = str(alignment[1])

    # Count mismatches and valid positions (ignoring gaps)
    mismatches = 0
    total_positions = 0
    
    for a, b in zip(aligned_seq1, aligned_seq2):
        if a != '-' or b != '-':  # Count all positions where at least one is valid
            total_positions += 1
            if a != b and a != '-' and b != '-':
                mismatches += 1
    
    # Calculate p-distance
    if total_positions == 0:
        return float('nan')  # Avoid division by zero if no valid positions
    p_distance = mismatches / total_positions
    return p_distance

def main(input_fasta, reference_fasta):
    input_seq = next(SeqIO.parse(input_fasta, "fasta")).seq
    input_id = next(SeqIO.parse(input_fasta, "fasta")).id  # Get the ID of the input sequence
    min_distance = float('inf')
    min_id = None

    # Read all reference sequences first to determine the total count
    reference_records = list(SeqIO.parse(reference_fasta, "fasta"))
    total_references = len(reference_records)
    
    print('Comparing input sequence to reference sequences...')

    for index, record in enumerate(reference_records):
        ref_seq = record.seq
        distance = calculate_p_distance(str(input_seq), str(ref_seq))  
        
        # Check if this distance is the smallest found so far
        if distance < min_distance:
            min_distance = distance
            min_id = record.id
            print(f"\nNew closest reference sequence: {min_id}, p-distance = {min_distance:.4f}")

        # Update progress after each comparison
        progress_percentage = (index + 1) / total_references * 100
        print(f"\rProgress: {progress_percentage:.2f}% ({index + 1}/{total_references})", end='')

    print()  # Print a newline after progress update

    # Print the comparison with the lowest p-distance at the end
    if min_id is not None:
        print(f"\033[92mClosest reference sequence is {min_id}, p-distance = {min_distance:.4f}\033[0m")
    
    if min_distance <= 0.1833:
        print("Distance is within subtype cutoff range. Proceeding to ML subtype inference...")
    elif 0.1833 < min_distance <= 0.2145:
        print("\033[93mP-distance to closest reference is greater than subtype cutoff. This sequence may belong to a novel subtype. Proceeding to ML subtype inference...\033[0m")
    else:
        print("\033[93mP-distance to closest reference is greater than clade cutoff. Please make sure this sequence belongs to rat hepatitis E. Proceeding to ML subtype inference...\033[0m")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python p_distance.py <input_fasta> <reference_fasta>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    reference_fasta = sys.argv[2]
    
    main(input_fasta, reference_fasta)
