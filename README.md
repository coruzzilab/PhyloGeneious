# PhyloGeneious
 High Performance Computing optimized orthology inference

PhyloGeneious is an improved version of the OrthologID pipeline, optimized for HPC clusters (Works on Slurm and PBS job schedulers).

### Dependencies:
The pipeline uses various programs to complete the analysis. Please make sure to install the following programs and add them to your $PATH variable:
1.  BLAST+
2.  MCL
3.  MAFFT
4.  TNT # Before running the pipeline for the first time, type "tnt", and agree to the terms and conditions
6.  Perl > v.5.20


### Run setup:
1. Create a project folder
2. Copy the "procfiles.txt" and "config" files from the provided testdata folder to your project folder.
   - The "config" file contains important instructions to run the pipeline (see instructions below).
   - The "procfiles.txt" contains instructions for building the gene family trees with TNT.

3. Create a "blastdb" subfolder in your project folder. This folder should include all the protein sequences (a fasta file) for each of the included species (see instructions below)


