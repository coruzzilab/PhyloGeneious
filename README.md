# PhyloGeneious
 High Performance Computing optimized parsimony-based orthology inference.

PhyloGeneious is an improved version of the OrthologID pipeline, optimized for HPC clusters (Slurm and PBS job schedulers).

![Pipeline](Pipeline_steps.png)

By using this pipeline, you also agree to the [TNT Personal Use License](https://www.lillo.org.ar/phylogeny/tnt/files/LicenseAgreement_1.5.htm).

## System requirements
This pipeline is under development and has only been tested and used on New York University's Prince and Greene HPC clusters (Slurm and PBS job schedulers). The basic requirement is to be able to start batch jobs from running batch jobs, with instant feedback of the job ID created. This way
various programs can not only submit jobs but also monitor them for
completion.

### Dependencies:

The Phylogeneious pipeline uses various programs to complete the analysis. Please make sure to install the following programs and add the path to their executables to your $PATH variable:

1.  BLAST+ 
2.  MCL
3.  MAFFT
4.  TNT # Before running the pipeline for the first time, type "tnt", and agree to the terms and conditions
5.  Perl > v.5.20 # Needs the DB_File.pm module
6.  RAxML (raxmlHPC-PTHREADS-AVX) # for the post-analysis total evidence tree search
7.  Python >= 3.6 # for post-analysis scripts
8.  R # for post-analysis GO enrichment

## Run setup:

To run the pipeline, the user can create a project folder by running  
```sh
setup_oidrundir.sh -d <Path to directory containing named sequence files> 
Optional arguments:
-a      Run directory name [default: OIDrun]
-s      Tab-delimited file with species IDs, names
-o      List of outgroup species IDs [eg. "Species1 Species2 Species3"]
```
or build the project folder manually. This folder should contain the following:
1. blastdb - a folder with the input protein sequence fasta files, where the file names should include a short species name label (e.g. Aratha for Arabidopsis thaliana) and an `.faa` suffix (e.g. `Aratha.faa`). The sequence headers in these fasta files should start with the species label follow by a "#" character and the sequence id (e.g. >Aratha#AT5G47240).
2. config - a configuration file telling basic runtime parameters:

   Required arguments:

   - `INGROUP`= List of taxa labels (faa files without the suffix) for the ingroup taxa separated by spaces.
   - `OUTGROUP`= At least one faa file (without suffix) for an outgroup taxon.
   - `HPC`= S (slurm) or P (pbs) - which job scheduling system is supported by the HPC system.

   Not required but recommended:

   - `NCPU`= Number of cpus to ask for starting a BLAST or TNT type job. `default=1`
   - `MAXQS`= Number of process to run simultaneously - nominally BLAST and TNT jobs. `default=10`
   - `TNTA`= Lowest family size to be pooled in a group job running TNT. `default=200`
   - `TNTB`= Family size requiring a stand-alone process to run TNT. `default=500`
   - `BLSTMIN`= Number of minutes for an average BLAST job to run. `default=60`
   - `SEARCHTYPE`= Program to use for sequence similarity search. B=BLASTP `default`, D=DIAMOND, M=mmseqs2
   - `TREEPROGRAM`= Program to use for gene tree building. TNT=TNT (protein-based) `default`, OBLONG=Oblong+TNT (codon-based) [IN BETA, SLURM only]

3. procfiles.txt - The set of TNT commands for searching the gene family trees (copy from distribution).
4. species.txt - a two column tab-delimited table of species short labels (e.g Athal) and species full names (e.g. Arabidopsis thaliana).
5. run.sh - a shell script to set the `OID_HOME` and `OID_USER_DIR` environment variables to indicate the paths to the PhyloGeneious folder (code) and to your project working directory, respectively (required before initiating a run):
```sh
export OID_HOME=/path/to/PhyloGeneious/folder
export OID_USER_DIR=/path/to/project/folder
```
   - This file is optional (variables can be defined manually), however, it can be helpful if the pipeline needs to be restarted.
   - For large analyses (usually above 30 taxa), you can increase the amount of requested memory per CPU (default 2GB per CPU), using the `OID_MCL` variable, e.g.:
```sh
export OID_MCL=6
```

#### Running command:
Once you have set all the above, including the environment variables, run the following commands while in your project folder:

```sh
. ./setoid.sh

# (You can also add PBS or SLURM args to this command, for instance "--mem 32GB" )
nohup sh $OID_HOME/bin/topshell.sh > run.log 2>&1 &
```

#### Test run:
To do a test run of the pipeline, copy the testdata directory provided here to your desired location, modify the paths in setoid.sh, and follow the above instructions. Run time will vary depending on resource availability, but a successful test run should complete in less than twelve hours, with TNT jobs running for less than an hour.

#### Some notes:
- A log of the main job is saved in the "toplog" subfolder
- The pipeline is set to cancels it self every 20 hr, and restart automatically. This is because often jobs on HPC have a walltime limit. In the furue we will allow to costumize this variable
- It creates ".done" files, so that if a job crashes, it picks up where it left (just restart the run again, don't forget the environment variables...)
- If TNT jobs are running for long periods without completion or log output, double-check that the TNT license has been agreed to.

### Outputs:
The pipeline will identify ortholog groups and create concatenated matrices and partition files:
1. Matrix.nex - Concatenated matrix in a nexus format for searching a total evidence (species) tree using PAUP `default: all ortholog groups`
2. Matrix.tnt - Concatenated matrix in a TNT expected format for searching a total evidence (species) tree using TNT. It includes only parsimony informative characters to save space/memory
3. Matrix.phy - Concatenated matrix in a PHYLIP format for searching a total evidence (species) tree using RAxML
4. partitions.txt - Partition (ortholog group) coordinates on the Matrix.phy for RAxML `default: all ortholog groups`
5. partitionMembers.txt - Partition sequence identifiers (one per member species) and the familyID from where the partition/ortholog group was identified.
6. partitionMembers_paralogs.txt - Same as partitionMembers.txt, but all paralog sequence identifiers are included.
7. MatrixRecording.log - Log file indicating the parsimony-informative sites identified in each partition (ortholog group) alignment, if any. Non-informative partitions are not included in the TNT matrix.
8. blast/ - Directory with sequence similarity search results.
9. data/ - Directory with sequence alignments and gene trees for each family cluster. Contains:
  - Gene family subfolder (numbered from 1, which is the largest family, till family X).
    - "FAMILY" file: fasta file containing all the gene family sequences
    - "FAMILY.aligned" file: fasta file containing the MAFFT alignment for all gene family sequences
    - "oid.tre" file: gene family tree file in a newick format
    - "orthologs" file: each line containing the sequence ids for an identifed ortholog group
    - all other files are less important...
  - Gene families with less than 4 sequences (those will start with "S" followed by a number) - these will not be further processed (alignment, tree search, and orthlog calling), and therefore not be included in the concatenated matrix, as they have too few sequences for tree search.
  - a "singlets" file, containing all the sequences that were not clustered into gene families using MCL.

Matrix and partition files can be recreated if some modifications are required using the PhyloGeneious/bin/orth2matrix.pl script (note that it usually requires a lot of memory, so run it on a high memory node. Usually ≥128GB, or even 250GB, depends of the size of the data set), e.g.:
```sh
export OID_HOME=/path/to/PhyloGeneious/folder

export OID_USER_DIR=/path/to/project/folder

# "-m" for setting a minmum cutoff for number of taxa per ortholog group, for that ortholog group to be included in concatenated matrix:

perl $OID_HOME/bin/orth2matrix.pl -m 10

# "-x" to exclude a given taxa from the concatenated matrix (for tree search). Multiple taxa can be indicated (comma-seperated):

perl $OID_HOME/bin/orth2matrix.pl -x Aratha

# "-O" flag will remove all partitions that do not contain a sequence from the outgroup:

perl $OID_HOME/bin/orth2matrix.pl -O
```

## Troubleshooting:

A. Data entry issues: Issues with data correctness and completeness can cause the pipeline to fail in certain cases. Here are some common issues:

1. Protein sequence files have to be in FASTA format and need to have a name that ends in .faa The pipeline does not recognize files that have a different name such as .fas OR .fasta

2. The species short name (e.g., SPECIES1) has to be identical in all three places: 1. config file where INGROUP and OUTGROUP are specified. 2. File name for fasta file (SPECIES1.faa) and 3. sequence headers in Fasta file (e.g., >SPECIES1#xxxxx) 

B. Pipeline customization: Compute clusters come in all sizes and configurations, with key differences in job management software, walltime limits, memory limits etc. It's impossible to have the pipeline auto-adjust to all these settings. So, please make sure the pipeline is customized for your cluster by doing the following:

1. Provide a template job submission script in $OID_HOME/bin/ For clusters using SLURM job management, name this file pipe.slu and if cluster uses PBS, name the file pipe.pbs In either case, this file should include all the arguments typically provided to the job manager. Example scripts are provided in the distribution and should be edited to fit your cluster.

2. In the file $OID_HOME/bin/run_pipeline.sh the very first non-comment line defines the name of the queue available to you for job submissions. Please edit this line to provide the correct queue name.

C. Failure modes: The PhyloGeneious pipeline uses an array of software tools in multiple stages. With certain input data sets, especially those that are very large, some of the steps will fail due to resource limitations or other unavoidable issues. 

1. After the gene families are reconstructed some of the largest gene families may be very large (> 1000 proteins). These large families need an extraordinary amount of compute resources, i.e., memory and walltime, to be phylogenetically resolved. If the ortholog resolution in a gene family keeps failing for any reason, you can skip the family by creating an empty oid.tre file in that family folder (e.g., for family 25 `touch $OID_USER_DIR/data/25/oid.tre`)

