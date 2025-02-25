import streamlit as st
import subprocess
import sys
import os
import json
from Bio import SeqIO
import time

def extract_fasta_header(input_fasta):
    # Extract the FASTA header name
    with open(input_fasta, 'r') as file:
        for record in SeqIO.parse(file, "fasta"):
            return record.id

def calculate_p_distance(input_fasta, reference_fasta):
    try:
        command = f"python3 p-distance-calc.py {input_fasta} {reference_fasta}"
        subprocess.run(command, shell=True, check=True)
        
        # Read the output from the file
        with open("output/p_distance_output.json", "r") as file:
            result = json.load(file)
        return result
    except subprocess.CalledProcessError as e:
        st.error(f"Error calculating p-distance: {e}")
        return None

def infer_new_tree(existing_alignment, new_sequence, existing_tree, output_dir):
    try:
        output_alignment = os.path.join(output_dir, f"{os.path.basename(new_sequence).split('.')[0]}_updated.fasta")
        output_tree = os.path.join(output_dir, f"{os.path.basename(new_sequence).split('.')[0]}_reoptimised")
        command = f"python3 infer_new_ML_tree.py {existing_alignment} {new_sequence} {existing_tree} {output_alignment} {output_tree}"
        subprocess.run(command, shell=True, check=True)
        
        # Read the output from the file
        with open("output/ml_tree_output.json", "r") as file:
            result = json.load(file)
        return result, output_alignment, output_tree
    except subprocess.CalledProcessError as e:
        st.error(f"Error inferring new ML tree: {e}")
        return None, None, None

def infer_subtype(input_newick, predefined_label, csv_file):
    try:
        command = f"python3 ML_patristic-dist_calc.py {input_newick} {predefined_label} {csv_file}"
        subprocess.run(command, shell=True, check=True)
        
        # Read the output from the file
        with open("output/subtype_output.json", "r") as file:
            result = json.load(file)
        return result
    except subprocess.CalledProcessError as e:
        st.error(f"Error inferring subtype: {e}")
        return None

def main():
    st.title("Rat Hepatitis E Subtyping Tool v1.0")
    st.header("Sridhar Group")
    
    # File upload with specific types
    input_fasta = st.file_uploader("Upload FASTA file", type=["fasta", "fas", "fa"])
    
    if input_fasta is not None:
        st.write("\n**Summary Statistics:**")
        with open(input_fasta.name, 'r') as file:
            for record in SeqIO.parse(file, "fasta"):
                st.write(f"**Query ID:** {record.id}")
                st.write(f"**Query Length:** {len(record.seq)}")
                break  # Assuming there's only one sequence in the file
        
        # Create output directory
        output_dir = "output"
        os.makedirs(output_dir, exist_ok=True)
        
        # Define fixed file names
        reference_fasta = "reference_genomes.fa"
        existing_alignment = "reference_alignment.fa"
        existing_tree = "reference_tree.tree"
        csv_file = "reference_subtypes.csv"
        
        # Progress bar
        progress_bar = st.progress(0)
        
        # Calculate p-distance
        st.write("\nCalculating p-distance...")
        for i in range(30):  # Move progress bar from 0 to 30%
            progress_bar.progress(i / 100)
            time.sleep(0.1)  # Wait for 0.1 seconds
        p_distance_output = calculate_p_distance(input_fasta.name, reference_fasta)
        if p_distance_output:
            st.success("P-distance calculation completed.")
        else:
            st.error("Failed to calculate p-distance.")
            return
        
        # Infer new ML tree
        st.write("\nInferring new ML tree...")
        for i in range(69):  # Move progress bar from 20 to 99%
            progress_value = 0.3 + (i / 69) * 0.69  # Calculate progress 20 to 99%
            progress_bar.progress(progress_value)
            time.sleep(0.1)  # Wait for 0.1 seconds
        tree_output, output_alignment, output_tree = infer_new_tree(existing_alignment, input_fasta.name, existing_tree, output_dir)
        if tree_output:
            st.success("New ML tree inference completed.")
        else:
            st.error("Failed to infer new ML tree.")
            return
        
        # Infer subtype
        st.write("\nInferring subtype...")
        input_newick = f"{output_tree}.treefile"
        predefined_label = extract_fasta_header(input_fasta.name)  # Use the FASTA header name as the taxon label
        subtype_output = infer_subtype(input_newick, predefined_label, csv_file)
        if subtype_output:
            st.success("Subtype inference completed.")
            progress_bar.progress(100)
        else:
            st.error("Failed to infer subtype.")
            return
        
        # Wait for a moment to ensure files are generated
        time.sleep(1)  # Wait for 1 second
        
        # Display results
        try:
            with open("output/p_distance_output.json", "r") as file:
                p_distance_result = json.load(file)
            with open("output/ml_tree_output.json", "r") as file:
                tree_result = json.load(file)
            with open("output/subtype_output.json", "r") as file:
                subtype_result = json.load(file)
            
            st.write("\nResults:")
            st.write(f"1. **Closest reference by p-distance**: {p_distance_result['closest_reference']} ({p_distance_result['p_distance']:.4f})")
            st.write(f"2. **p-distance below cutoff**: {p_distance_result['below_cutoff']}")
            st.write(f"3. **Closest reference by ML patristic distance**: {subtype_result['closest_reference_ml']} ({subtype_result['ml_distance']:.4f})")
            st.write(f"4. **Conflicts in subtype assignment**: {subtype_result['conflicts']}")
            st.write(f"5. **Subtype assignment using ML patristic distance**: {subtype_result['subtype_assignment']}")
            
        except FileNotFoundError as e:
            st.error(f"Error loading results: File not found - {e}")
        except json.JSONDecodeError as e:
            st.error(f"Error loading results: JSON decode error - {e}")
        except KeyError as e:
            st.error(f"Error loading results: Key error - {e}")
        except Exception as e:
            st.error(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
