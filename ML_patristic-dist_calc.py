import sys
import dendropy
import os
import csv
from collections import Counter

def main(input_newick, predefined_label, csv_file):
    # Check if the input files exist
    if not os.path.isfile(input_newick):
        print(f"File '{input_newick}' does not exist.")
        return
    if not os.path.isfile(csv_file):
        print(f"File '{csv_file}' does not exist.")
        return

    try:
        # Load the tree from the Newick file
        tree = dendropy.Tree.get(
            path=input_newick,
            schema="newick"
        )
    except Exception as e:
        print(f"Failed to load tree from '{input_newick}': {e}")
        return

    # Create a phylogenetic distance matrix
    pdc = tree.phylogenetic_distance_matrix()

    # Find the predefined taxon
    predefined_taxon = None
    for taxon in tree.taxon_namespace:
        if taxon.label == predefined_label:
            predefined_taxon = taxon
            break

    if predefined_taxon is None:
        print(f"Taxon '{predefined_label}' not found in the tree.")
        return

    # List to store distances and corresponding taxa
    distance_list = []

    # Calculate distances from the predefined taxon to all other taxa
    for taxon in tree.taxon_namespace:
        if taxon != predefined_taxon:  # Skip the predefined taxon itself
            distance = pdc(predefined_taxon, taxon)
            distance_list.append((distance, taxon.label))

    # Sort the list by distance (first element of tuple)
    distance_list.sort(key=lambda x: x[0])

    # ML nucleotide threshold
    threshold = 0.5585

    # Filter distances under the threshold
    filtered_distances = [(dist, label) for dist, label in distance_list if dist < threshold]

    # Read the CSV file
    csv_data = {}
    with open(csv_file, mode='r', encoding='utf-8-sig') as file:
        reader = csv.DictReader(file, delimiter=',')
        for row in reader:
            csv_data[row['Name']] = {'Clade': row['Clade'], 'Subtype': row['Subtype']}

    # Report distances under the threshold
    print(f"Distances under subtyping threshold ({threshold})")
    for dist, label in filtered_distances:
        if label in csv_data:
            print(f"{label}: Distance={dist:.4f}, Clade={csv_data[label]['Clade']}, Subtype={csv_data[label]['Subtype']}")
        else:
            print(f"{label}: Distance={dist:.4f}, Not found in the CSV file.")

    # Check for conflicting clades and subtypes
    clades = []
    subtypes = []
    for dist, label in filtered_distances:
        if label in csv_data:
            clades.append(csv_data[label]['Clade'])
            subtypes.append(csv_data[label]['Subtype'])

    if clades and subtypes:
        if len(set(clades)) > 1 or len(set(subtypes)) > 1:
            print(f"\033[93mConflicting clade/subtype inferred. Cutoff thresholds may need to be updated. Please contact the authors.\033[0m")
        else:
            most_common_clade = Counter(clades).most_common(1)[0][0]
            most_common_subtype = Counter(subtypes).most_common(1)[0][0]
            print(f"\033[92mInferred clade and subtype for '{predefined_label}': Clade={most_common_clade}, Subtype={most_common_subtype}\033[0m")
    else:
        print(f"No clade and subtype could be inferred for '{predefined_label}'. This query may belong to a novel subtype.")

    # Extract closest label after processing filtered distances
    closest_label = None
    if filtered_distances:  # Only try to get the label if filtered_distances is not empty
        closest_label = filtered_distances[0][1]  # First element's label
    # Now print the closest label
    if closest_label:
        print(f"\033[92mClosest reference genome by ML patristic distance is {closest_label}\033[0m")
    else:
        print(f"\033[91mNo reference genome under the distance threshold found.\033[0m")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_newick> <predefined_label> <csv_file>")
        sys.exit(1)

    input_newick = sys.argv[1]
    predefined_label = sys.argv[2]
    csv_file = sys.argv[3]
    
    main(input_newick, predefined_label, csv_file)
