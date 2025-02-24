import subprocess
import sys

def add_sequence_to_msa(existing_alignment, new_sequence, output_alignment):
    mafft_command = f"mafft --quiet --add {new_sequence} --keeplength {existing_alignment} > {output_alignment}"
    subprocess.run(mafft_command, shell=True)

def add_sequence_to_tree(output_alignment, existing_tree, output_tree):
    iqtree_command = f"iqtree2 -redo --quiet -s {output_alignment} -g {existing_tree} -pre {output_tree} -m GTR+F+G4"
    subprocess.run(iqtree_command, shell=True)

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
    
    # Add new sequence to existing tree
    add_sequence_to_tree(output_alignment, existing_tree, output_tree)

    print("\nOperation completed successfully.")

if __name__ == "__main__":
    main()
