# Identify structural deletions with homeology
This resource deposites the scripts used to identify junction-associated break ends, or SV deletions with homeology, where a pair of genomic sequences with high similarity can be used for single-strand anealing. 

## Algorithm
The prediction approach is similar to Setton et al. 2023 paper, where they applied a sliding bin approach to search for the similar sequences around each end break pair. 
In this method, we search for the similar sequences using BLASTn algorithm. More detailed, we extracted 200 bp sequences flanking the junction sites, extending 100 bp upstream and downstream. These sequences were then aligned using BLASTn to identify the most similar regions. Alignments were considered homeologous if they exhibited at least 80% similarity, a minimum alignment length of 30 bp, and were located on the same strand. When the best alignments have the same similarity, the longest alignment is selected.

## Installation and System Requirements
For detailed installation instructions and requirements, please see [INSTALLATION.md](INSTALLATION.md).

## Directory Structure Setup
The script expects a specific directory structure. Make sure you have:
```
.
├── example/
│   ├── input/
│   │   └── structural_deletions_example.tsv
│   └── output/ (will be created by the script)
├── lib/ (contains library scripts)
└── R/ (contains main R scripts)
```
## Input Data Format
The script expects deletion data in a tab-separated format (TSV) with the following columns:
- SAMPLE.TUMOR: Sample identifier
- CHROM: Chromosome
- start_position: Start position of the deletion
- end_position: End position of the deletion
## Running the Analysis
To run the example analysis:
```bash
# Make sure your conda environment is active
conda activate svdelhomeo
# Edit the parameters (deletions input file, output directory)
vi R/predict_svdelhomeo_example.R
# Run the analysis
Rscript R/predict_svdelhomeo_example.R
```
## Output
The script will generate several output files in the `example/output` directory:
- `ssa_events_candidates.tsv`: Detailed information about deletions with homeology
- `ssa_events_samples.tsv`: Summary of samples with homeology events, the proportion of deletions with homeology is calculated.

The demo run for the 50 deletions takes about 5-10 mins.

## License
This project is covered under the **MIT License**.