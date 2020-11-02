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

#### Input protein fasta files:
- Give a simple and short name for each of the included species. We usually take the first 3 letters from the genus and species names (e.g. "Aratha" for "Arabidopsis thaliana").
- The input fasta file name should be the species short name + ".faa", e.g. Aratha.faa.
- Each sequence header in the faa file should start with the species short name + "#" + sequence id, e.g. ">Aratha#AT1G50030"
- Make sure there are no special characters in the sequence headers (e.g. "@","-") or in the peptid sequences (e.g. "*" stop codons).

#### config file:
1. Required arguments:
- INGROUP= 	List of species short names for the ingroup taxa, separated by spaces.
- OUTGROUP= 	At least one species short name for an outgroup taxon (multiple outgroup taxa should be space-separated)
- HPC=  		S (slurm) or P (pbs) - which job scheduling system is supported by the HPC system.
2. Not required but recommended:
- NCPU=  		Number of cpus to ask for starting a BLAST or TNT type job.
- MAXQS= 		Number of process to run simultaneously - nominally BLAST and TNT jobs.
- TNTA=  		Lowest family size to be pooled in a group job running TNT.
- TNTB=  		Family size requiring a stand-alone process to run TNT.
- BLSTMIN=  	Number of minutes for an average BLAST job to run.

#### Environment variables:
- Use the "$OID_HOME" and "$OID_USER_DIR" environment variables to indicate the path to the PhyloGeneious folder (code) and to your project working directory, respectively (You have to do that before initiating a run):
export OID_HOME=/path/to/PhyloGeneious/folder
export OID_USER_DIR=/path/to/project/folder
- For large analyses (usually above 30 taxa), you can increase the amount of requested memory per CPU (default 2GB per CPU), using the $OID_MCL variable, e.g.:
export OID_MCL=6 


### Disclamer:
- This pipeline is under development and only been tested and used on New York University's HPC cluster (Slurm and PBS job schedulers). While it is likely to work on PBS or Slurm HPC systems, it is not guaranteed.

