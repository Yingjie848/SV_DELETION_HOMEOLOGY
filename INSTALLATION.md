# SV Deletion Homeology Analysis - Requirements
This document outlines the requirements and installation steps needed to run the SV (Structural Variant) Deletion Homeology analysis scripts.
## System Requirements
- **Operating System**: Linux (tested on Linux)
- **R**: Version 4.0.0 or higher
- **BLAST+**: NCBI BLAST+ suite (specifically blastn and makeblastdb)
- **Reference Genome**: Human reference genome GRCh37/b37
## Installation Steps
### 1. Set up Conda Environment for R and Dependencies
Using conda makes it easier to manage dependencies and isolate the environment for this project:
```bash
# Install Miniconda if you don't have it already
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# Create and activate a new conda environment for the project
conda create -n svdelhomeo r-base=4.0 r-essentials
conda activate svdelhomeo
# Install required R packages within the conda environment
conda install -c conda-forge r-data.table r-magrittr r-dplyr r-ggplot2 r-ggforce r-ggpubr r-tidyr
conda install -c bioconda bioconductor-biostrings
# Verify R installation
R --version
```
Alternatively, you can create the environment from a YAML file:
```yaml
name: svdelhomeo
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - r-base=4.0
  - r-essentials
  - r-data.table
  - r-magrittr
  - r-dplyr
  - r-ggplot2
  - r-ggforce
  - r-ggpubr
  - r-tidyr
  - bioconductor-biostrings
```
Save this YAML to `environment.yml` and run:
```bash
conda env create -f environment.yml
conda activate svdelhomeo
```
### 2. Install BLAST+ from Source
Install BLAST+ from the official NCBI source:
```bash
# Download the latest BLAST+ source
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.15.0+-src.tar.gz
# Extract the archive
tar -xzvf ncbi-blast-2.15.0+-src.tar.gz
cd ncbi-blast-2.15.0+-src/c++
# Configure and build
./configure
make
# Install (requires admin privileges)
sudo make install
# If you don't have admin privileges, you can install to a custom location
# ./configure --prefix=$HOME/blast
# make
# make install
# Then add $HOME/blast/bin to your PATH
# Verify installation
blastn -version
makeblastdb -version
```
### 3. Reference Genome Setup
The script requires the GRCh37/b37 reference genome. You need modify the path in `lib/lib_blast.R` to point to your reference genome file.