# Identify deletions with homeology
This resource deposites the scripts used to identify junction-associated break ends, or SV deletions with homeology, where a pair of genomic sequences have high similarity. 

# Algorithm:
The prediction approach is similar to Setton et al. 2023 paper, where they applied a sliding bin approach to search for the similar sequences around each end break pair. 
In this approach, we search for the similar sequences using BLASTn algorithm. More detailed, we extracted 200 bp sequences flanking the junction sites, extending 100 bp upstream and downstream. These sequences were then aligned using BLASTn to identify the most similar regions. Alignments were considered homeologous if they exhibited at least 80% similarity, a minimum alignment length of 30 bp, and were located on the same strand. When the best alignments have the same similarity, the longest alignment is selected.

# Usage:
## Input:
Deletions identified by structural variant predictors, at least include sample name, chromosome, positions of two break ends.

## Output:
Deletions with similar sequences around the break ends meeting the criteria mentioned in Algorithm.