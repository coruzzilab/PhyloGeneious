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




Required arguments:
        INGROUP= 	List of taxa labels (faa files without the suffix) for the ingroup taxa separated by spaces.
        OUTGROUP= 	At least one faa file (without suffix) for an outgroup taxon.
        HPC=  		S (slurm) or P (pbs) - which job scheduling system is supported by the HPC system.
Not required but recommended:
        NCPU=  		Number of cpus to ask for starting a BLAST or TNT type job.
        MAXQS= 		Number of process to run simultaneously - nominally BLAST and TNT jobs.
        TNTA=  		Lowest family size to be pooled in a group job running TNT.
        TNTB=  		Family size requiring a stand-alone process to run TNT.
        BLSTMIN=  	Number of minutes for an average BLAST job to run.
