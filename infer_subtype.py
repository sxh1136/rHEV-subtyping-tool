import subprocess
import sys
import os
from Bio import SeqIO

def extract_fasta_header(input_fasta):
    # Extract the FASTA header name
    with open(input_fasta, 'r') as file:
        for record in SeqIO.parse(file, "fasta"):
            return record.id

def calculate_p_distance(input_fasta, reference_fasta):
    # Run the first script to calculate p-distance
    command = f"python3 p-distance-calc.py {input_fasta} {reference_fasta}"
    subprocess.run(command, shell=True)

def infer_new_tree(existing_alignment, new_sequence, existing_tree, output_dir):
    # Run the second script to infer a new ML tree
    output_alignment = os.path.join(output_dir, f"{os.path.basename(existing_alignment).split('.')[0]}_updated.fasta")
    output_tree = os.path.join(output_dir, f"{os.path.basename(existing_tree).split('.')[0]}_reoptimised")
    
    command = f"python3 infer_new_ML_tree.py {existing_alignment} {new_sequence} {existing_tree} {output_alignment} {output_tree}"
    subprocess.run(command, shell=True)
    
    return output_alignment, output_tree

def infer_subtype(input_newick, predefined_label, csv_file):
    # Run the third script to infer subtype
    command = f"python3 ML_patristic-dist_calc.py {input_newick} {predefined_label} {csv_file}"
    subprocess.run(command, shell=True)

def main(input_fasta):
    # Define fixed file names
    reference_fasta = "reference_genomes.fa"
    existing_alignment = "reference_alignment.fa"
    existing_tree = "reference_tree.tree"
    csv_file = "reference_subtypes.csv"

    # Create output directory
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    # Calculate p-distance
    print("\nCalculating p-distance...")
    calculate_p_distance(input_fasta, reference_fasta)

    # Infer new ML tree
    print("\nInferring new ML tree...")
    output_alignment, output_tree = infer_new_tree(existing_alignment, input_fasta, existing_tree, output_dir)

    # Infer subtype
    print("\nInferring subtype...")
    input_newick = f"{output_tree}.treefile"
    predefined_label = extract_fasta_header(input_fasta)  # Use the FASTA header name as the taxon label
    infer_subtype(input_newick, predefined_label, csv_file)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python infer_subtype.py <input_fasta>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    
    main(input_fasta)
