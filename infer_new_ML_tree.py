import subprocess
import sys

def add_sequence_to_msa(existing_alignment, new_sequence, output_alignment):
    mafft_command = f"mafft --quiet --add {new_sequence} --keeplength {existing_alignment} > {output_alignment}"
    try:
        subprocess.run(mafft_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to add sequence to alignment: {e}")
        sys.exit(1)

def run_phylogenetic_placement(output_alignment, existing_tree):
    # Run IQ-TREE with the guide tree
    iqtree_command = f"iqtree2 -redo --quiet -s {output_alignment} -g {existing_tree} -pre {output_alignment}_pp -m GTR+F+G4"
    try:
        subprocess.run(iqtree_command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to perform phylogenetic placement: {e}")
        sys.exit(1)


def infer_global_optimization_tree(output_alignment, output_tree):
    # Run IQ-TREE with the constraint tree for optimization
    iqtree_command2 = f"iqtree2 -redo --quiet -s {output_alignment} -t {output_alignment}_pp.treefile -pre {output_tree} -m GTR+F+G4"
    try:
        subprocess.run(iqtree_command2, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to infer optimized tree: {e}")
        sys.exit(1)

def main():
    if len(sys.argv) != 4 and len(sys.argv) != 6:
        print("Usage: python script_name.py existing_alignment.fasta new_sequence.fasta existing_tree.treefile [output_alignment.fasta output_tree_prefix]")
        sys.exit(1)
    
    existing_alignment = sys.argv[1]
    new_sequence = sys.argv[2]
    existing_tree = sys.argv[3]
    
    if len(sys.argv) == 6:
        output_alignment = sys.argv[4]
        output_tree = sys.argv[5]
    else:
        output_alignment = f"{existing_alignment.split('.')[0]}_updated.fasta"
        output_tree = f"{existing_tree.split('.')[0]}_reoptimised"
    
    print(f"\nOutput alignment will be saved as: {output_alignment}")
    print(f"Output tree will be saved as: {output_tree}.treefile")
    
    # Add new sequence to existing alignment
    add_sequence_to_msa(existing_alignment, new_sequence, output_alignment)
    
    # Run IQ-TREE phylogenetic placement first
    run_phylogenetic_placement(output_alignment, existing_tree)
    
    # Run IQ-TREE with the constraint tree for optimization
    infer_global_optimization_tree(output_alignment, output_tree)

    print("\nOperation completed successfully.")

if __name__ == "__main__":
    main()