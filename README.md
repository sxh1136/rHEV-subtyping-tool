# Tool for subtyping rHEV sequences using p-distances and maximum likelihood patristic distances. 

Accompanies manuscript under review.

Note: Currently only supports individual sequences. If multifasta is passed to script, only the first sequence will be analysed. 

Tool workflow:
1. Calculates p-distance of input sequence to reference genomes.
2. Returns reference with the lowest p-distance.
3. Adds input sequence into an existing reference MSA with mafft --add.
4. Performs phylogenetic placement using existing maximum likelihood constraint tree built from reference genomes.
5. Calculates patrisitic distance of input sequence to reference genomes.
6. Reports references with patristic distance to input less than cutoff threshold and checks for conflicts in subtype assignment.
7. Returns inferred clade and subtype of input sequence.

If p-distance and patristic distance of input sequence to reference genomes are higher than cutoff thresholds, this will be noted by tool. This may occur in the follow scenarios:
  1. Input sequence is highly divergent from reference genomes and belongs to a novel subtype.
  2. Input sequence is highly divergent from reference genomes and does not belong to rHEV.

Hypothetically, conflicts in subtype assignment can occur and will be reported by the tool (step 6). For example, imagine a scenario where 10 references with patristic distance less than the cutoff threshold to the input sequence are reported. It is possible for 9 of the references to belong to one subtype, with the one remaining reference belonging to another. If this occurs, an update of the cutoff thresholds using the new input sequence as a reference is required. 

### Usage
```
infer_subtype.py <input.fasta>
```
