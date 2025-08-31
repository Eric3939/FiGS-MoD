# linear_motif

**Automating Linear Motif Predictions to Map Human Signaling Networks**
Yitao (Eric) Sun, Yu (Brandon) Xia, Jasmin Coulombe-Huntington,
Department of Bioengineering, McGill University

## Table of Contents
- [Abstract](#abstract)
- [Running the algorithm](#running-the-algorithm)
- [Predicted Motif Database](#predicted-motif-database)
- [Paper code repository](#paper-code-repository)
- [Citation](#citation)

## Abstract
Short linear motifs (SLiMs) are short sequence patterns that mediate transient protein-protein interactions, often within disordered regions of proteins. SLiMs play central roles in signaling, trafficking, and post-translational regulation, but their short length and low complexity make them difficult to identify both experimentally and computationally. Since the release of motif discovery tools like MEME Suite, the availability of protein-protein interaction data (e.g., BioGRID) has increased by more than five-fold, and recent advances in machine learning offer new opportunities for large-scale, high-resolution motif discovery. Here, we present a new Gibbs sampling-based SLiM discovery method that introduces two key innovations: First, we replace the traditional position-specific scoring matrix (PSSM) with a Hidden Markov Model (HMM) to better accommodate insertions and deletions common in disordered regions; Second, we introduce biased sampling guided by pre-annotated residue-level features derived from Protein Language Models (PLMs), AlphaFold2-derived predictions (disorder, solvent accessibility), and evolutionary conservation. We evaluate our approach using the ELM database and show improved recovery of known motif instances compared to existing tools. We presented three case studies at the end to showcase potential applications of our algorithm. 

<img width="1200" height="800" alt="image" src="https://github.com/user-attachments/assets/debcdcef-579a-43c2-a7ad-ec1cf1b0d472" />


## Running the algorithm
Before you run:<br>
  We call a hub protein (a protein that interacts with several other proteins) a linear motif binding domain protein (**LMBD protein**). The algorithm (method 1 and 3) finds all the interacting proteins of a LMBD protein and search for a motif within them. We call these interacting proteins a **LMBD protein network**. The protein-protein interactions are obtained from BioGRID, filtered as described in our paper (restricted to low-throughput data or interactions supported by ≥2 sources).


### Step 1 - Clone repository and install dependencies
Clone this repository to your local machine:
```sh
git clone https://github.com/Eric3939/linear_motif.git
cd linear_motif
```
Install the required packages. 
```sh
pip install -r requirements.txt
```
**Note**: `pomegranate==0.15.0` is a strict requirement, as the latest version 1.0.0 does not support HMM hidden states. Earlier versions may work but have not been tested.

### Step 2 - Download data
Download these files from Google Drive and place them in the `data/` directory:
https://drive.google.com/drive/folders/1472iWG8U6g5XaJBz2bdI_UI-kOFbpU-n?usp=sharing
- `protein_database_1.pickle`: contains all human proteins used in our study, with pre-annotated feature scores (PLM, disorder, solvent accessibility, conservation). 
- `biogrid_net.gpickle`: contains all human PPIs from BioGRID, annotated with number of citations and throughput.

### Step 3 - Running models
The algorithm can be run at three different scales:
1. a single LMBD protein network in BioGRID
2. a user-defined set of proteins
3. all LMBD protein networks in BioGRID

**1) Single LMBD protein network (BioGRID proteins)**
```sh
cd script/
python search.py [LMBD_protein] [output_directory]
```
In this mode, the user provides only a single LMBD protein. The algorithm retrieves all interacting partners of the LMBD protein from BioGRID and searches for motifs within them (a LMBD protein network). All proteins are pre-annotated in our database. This is the simplest way to run the algorithm, but it does not allow customization of the input set.

**2) User defined proteins**
```sh
cd script/
```
Prepare a text file with one UniProt ID per line.
```sh
python run_users_proteins.py [proteins_file_path]
```
In this mode, the user specifies a custom set of proteins (via UniProt IDs) that are hypothesized to share a common motif. The algorithm retrieves these proteins from our pre-annotated database.
Note: Our protein database contains 20419 human proteins from UniProt. Proteins not included in our database are not currently supported.

**3) Full proteome (all LMBD networks)**
```sh
cd script/
python run_proteome.py [biogrid_path] [output_folder]
```
This mode reproduces the results described in our paper by running motif discovery across all LMBD protein networks in BioGRID. Parameters can be adjusted by modifying `search.py`.

**(Optional) High-Performance Computing with SLURM**
If user is using a High Performance Computing (HPC) system with SLURM scheduler, we provide a helper script to submit motif-search jobs automatically (adaptation to local pipelines may be required):
```sh
cd script/
python submit_slurm.py [biogrid_path] [results_folder] [dataframe_path]
```

### Step 4 - Results
Each run produces one pickle file per LMBD protein network, saved in the user-specified output folder. For a full proteome run (method 3), the algorithm will therefore generate one pickle file for each LMBD protein network.  
To aggregate results, use the script `result_table.py` (in the `script/` folder) to convert all pickle files into a single table listing the discovered motifs.


## Predicted Motif Database
We applied our algorithm to the entire human proteome as described in the paper, generating 221,840 predicted motifs. The results are provided in `results_table.csv`.

Example records from the database:

```sh
protein  motif    start  end  length  info_content        plm     disorder  solvent_acc  conservation
A0AV96   PPPFQGR  510    516  7       3.251699685122221   1.584   1.128     1.165        1.548
A0AV96   SDSAAGS  12     18   7       2.2353015793897852  0.544   1.38      0.4985       -1.421
A0AV96   STAAMSS  6      12   7       2.415746185344202   1.232   1.299     0.624        0.1235
A0AV96   AAAVIP   500    505  6       2.1140273005497536  0.01503 1.332     0.6396       -1.903
```


## Paper code repository
Additional scripts used in the paper are provided in the `repostory/` folder, including those used to compute the four feature scores (PLM, disorder, solvent accessibility, conservation).


## Citation
If you use this repository in your research, please cite our paper:
