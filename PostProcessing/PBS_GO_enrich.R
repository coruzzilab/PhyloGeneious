# Gil Eshel, Sep 14, 2017
# Script to run GO overrepresentation analysis for the high pbs partitions for each node on the tree
#
# # synopsis : Rscript PBS_GO_enrich.R </path/to/parse_partitionMembers.txt file> </path/to/pbs.csv file> </path/to/gene2go_with_parents.txt> </path/to/GoDesc.txt> <node_number_to_label.txt> <pbs absolute cutoff>
#
# Outputs:
# partition_info.txt (all the information per partition), partition_annotation_with_pbs.txt and partition_annotation_with_pbs_filtered.txt - tables with the clade (number and label), partition, pbs, model species identifer and GO annotation information (for all partitions or partitions with abs(PBS) >= PBS_cutoff) 
# 
# Commant1 - The user needs to input the gene2go_parent.txt tab-delim file with two columns: Gene_id and GO_id list (comma-separated). For model organisms, the user can download gene2go mapping files (e.g for plants: http://plantregmap.cbi.pku.edu.cn/download.php) or make one out of the Gene Ontology Consortium GAF files (http://www.geneontology.org/page/download-annotations and see http://geneontology.org/page/go-annotation-file-gaf-format-21). The user needs to run the script parse_obo_add_parentGO.py to retrieve parent GO_ids (input gene2go.txt, output gene2go_with_parents.txt), and input the gene2go_with_parents file to this script to properly account for the GO hierarchy structure in the enrichment analysis. The user have to make sure that the gene_ids in the gene2go.txt file correspont exactly to the input sequence identifers (e.g.if a sequence header is Arath_At4g25480.1, the gene_id in gene2go needs to be At4g25480.1, not At4g25480 nor AT4G25480.1) - (* try to change so it will at least be case-insensitive that with ignore.case=TRUE
# Commant2 - This code needs to hendle more than one specified model species (e.g. the user can specify Arabidopsis and maize) for GO annotation - We can suggest that the gene2go file will include the GO annotations for all the model organisms, and so we can loop at each partition (if ";" exist) and will get all the GO_ids for that partition (based on more than one species GO annotation), and keep only unique set
# 

### Get PhyloGeneious environment variables ###
#oid_home_path = Sys.getenv("OID_HOME") 
#oid_user_dir_path = Sys.getenv("OID_USER_DIR")

### Inputs ###

args = commandArgs(TRUE)

# read the input files and the arguments:
partition_info = read.delim(args[1], sep = "\t", header=TRUE) # Input the parse_partitionMembers.txt file (with the Model_species_identifier column)
pbs = read.csv(args[2], header=TRUE) #  Input the pbs.csv file
gene2go = read.delim(args[3], sep = "\t", header=FALSE) # Input the gene2go_with_parents.txt file for a model organism(s) that will be used for the PBS GO enrichment analysis
goDesc = read.delim(args[4], sep = "|", header=FALSE, strip.white=TRUE) # Input the 'GO_id | GO_term' information file - use it for the GO enrichment tables
node_number_to_label = read.delim(args[5], sep = "\t", header=TRUE) # Input the node_number_to_label.txt for the node labels

# Check for the pbs_cutt argument, if not specify, set it to "4":
if (length(args)>5){
  if (!is.na(as.numeric(args[6]))) {
    pbs_cutt = as.numeric(args[6]) # Set the pbs_cutt according to the user specification
  } 
  else {
    pbs_cutt = 4 # Set the default value (4) 
  }
}



### Parse the partition and gene2go information ###

# Currently, the python script that generate parse_partitionMembers.txt allow the user to specify more than one species as a model (i.e. the Model_species_identifier column may have more than one identifier if more than one model was specified - ";" seperated)
partition_info$Model_species_identifier = as.character(partition_info$Model_species_identifier)
for (i in 1:length(partition_info$Model_species_identifier)) { # For partitions without model species identifier, put the PartitionID instead of "NA" so that they will still be represented in PhyloBrowse
  if (is.na(partition_info$Model_species_identifier[i])) {
    partition_info$Model_species_identifier[i]=paste("Partition",partition_info$PartitionID[i], sep="_")
  }
}  

# loop each object in the vector partition_info$Model_species_identifier to get the GOids for all the identifiers of that partition
GO_id = c() # store the GO_ids in this vector (for each partition a sting with a list of GO_ids or "NA")
for (i in 1:length(partition_info$Model_species_identifier)) {
	model_id = unlist(strsplit(partition_info$Model_species_identifier[i],split=";")) # split if more than one identifier exist
	goid = c()
	for (i in 1:length(model_id)){ # can handle more than one identifier
    	goid = c(goid, unlist(strsplit(as.character(gene2go$V2[match(model_id[i], gene2go$V1)]), ",", fixed = TRUE))) # find GO_ids for all identifiers and concatinate to one vector
		goid = unique(goid) # remove redundant GO_ids
		goid = goid[!is.na(goid)]# get rid of "NA" and put it back if there are any annotated GO_ids (if for one identifer there are GO_ids and for the other identifier(s) there are none...)
	}
	if (length(goid) == 0){
      goid = c("NA")
	}
	goid = paste(goid,collapse=",")
	GO_id = c(GO_id, goid)
}

partition_info$GO_id = GO_id # add the geneToGo info to the partition_info table

partition_info[is.na(partition_info)] = "NA" # change NA to "NA" for the partitions that miss either the model organism ortholog or its GO annotation
write.table(partition_info, sep="\t", file="partition_info.txt", quote = FALSE, row.names = FALSE, col.names = TRUE) # save partition_info

### Add partition and node information to the PBS table ###
pbs$node_id = node_number_to_label$label[match(pbs$clade,node_number_to_label$number)] # add a node label column (i.e. the "N#"), based on the get_node_labels.R script output
pbs$taxa_occupancy = partition_info$Partition_taxa_occupancy[match(pbs$partition,partition_info$PartitionID)] # add taxa occupancy
pbs$model_gene_id = partition_info$Model_species_identifier[match(pbs$partition,partition_info$PartitionID)] # add model species gene ids for each partition
pbs$GO_id = partition_info$GO_id[match(pbs$partition,partition_info$PartitionID)] # add GO_id to each partition
pbs[is.na(pbs)] = "NA" # put "NA" for missing values

# Save the pbs data frame into a file that will be used to generate the bpPBS.json file for PhyloBrowse using python. PhyloBrowse will need the clade, model_gene_id and the pbs, plus the GO enrichment results...
write.table(pbs, sep="\t", file="partition_annotation_with_all_pbs.txt", quote = FALSE, row.names = FALSE, col.names = TRUE) # save

# Filter the partitions based on the pbs threshold (default -4,4 if not specified):
pbs_filtered = subset(pbs, pbs >= pbs_cutt | pbs <= -pbs_cutt , select=c(clade, partition, pbs, node_id, taxa_occupancy, model_gene_id, GO_id))
write.table(pbs_filtered, sep="\t", file="partition_annotation_with_pbs_filtered.txt", quote = FALSE, row.names = FALSE, col.names = TRUE) # Save the pbs_filtered if will want to use it

### Run GO enrichment on the filtered PBS genes for each node ###
# generating a unique list of the node_ids to run GO enrichment on each node seperatly  
node = unique(pbs_filtered$node_id)

# Define the universe/background gene list and count gene frequency per GO_id for the background:
universe_Gene_list = partition_info$Model_species_identifier[partition_info$GO_id!="NA"] # generate a list of "geneIDs" - the background for the GO enrichment - of the model organism(s), for all partition that were included in the phylogenetic matrix, and have GO annotation
universe_GO_list = partition_info$GO_id[partition_info$GO_id!="NA"] # generate a list of "GOID" lists - the background for the GO enrichment - for counting genes per GO_id
universe_GO_list = unlist(strsplit(paste(universe_GO_list, collapse = ','), split=",")) # generate a vector with all GO_ids (a redundant set - each occurrence = a different gene)
universe_GO_freq = as.data.frame(table(universe_GO_list)) # dataframe with GO frequencies in the universe/background set of genes 

# extract the filtered PBS genes for each node and run GO enrichment:
for(i in node){
  node_PBS_model_gene_list = pbs_filtered$model_gene_id[pbs_filtered$node_id==i & pbs_filtered$GO_id != "NA"]  # list of partitions with absPBS > pbs_cutoff for node i & gene model ids and GO annotation - will be included in the GO enrichment analysis
  node_PBS_model_GO_list = pbs_filtered$GO_id[pbs_filtered$node_id==i & pbs_filtered$GO_id != "NA"]  # list of the GO_ids of the node_PBS_model_gene_list partitions
  node_PBS_model_GO_list = unlist(strsplit(paste(node_PBS_model_GO_list, collapse = ','), split=",")) # generate a vector with all node_PBS_model_gene_list GO_ids (a redundant set - each occurrence = a different gene)
  node_PBS_model_GO_freq = as.data.frame(table(node_PBS_model_GO_list)) # dataframe with GO frequencies in the node_PBS_model_gene_list set of genes 
  
  fisher_table = data.frame(node_PBS_model_GO_freq$node_PBS_model_GO_list) # generate a dataframe to store all the values for fisher.test, and results
  colnames(fisher_table) = c("GOID")
  fisher_table$GOTERM = goDesc$V2[match(fisher_table$GOID,goDesc$V1)]  ###******!!!!!!
  fisher_table$Filtered_PBS_with_GOID = node_PBS_model_GO_freq$Freq[match(fisher_table$GOID,node_PBS_model_GO_freq$node_PBS_model_GO_list)]
  fisher_table$Filtered_PBS_total = rep(length(node_PBS_model_gene_list),length(fisher_table$GOID))
  fisher_table$BG_with_GOID = universe_GO_freq$Freq[match(fisher_table$GOID,universe_GO_freq$universe_GO_list)]
  fisher_table$BG_total = rep(length(universe_Gene_list),length(fisher_table$GOID))  
  
  Filtered_PBS_without_GOID = fisher_table$Filtered_PBS_total - fisher_table$Filtered_PBS_with_GOID
  BG_without_GOID = fisher_table$BG_total - fisher_table$BG_with_GOID
  
  pvalue = c()
  for(j in 1:length(fisher_table$GOID)){
    pvalue = c(pvalue, fisher.test(matrix(c(fisher_table$Filtered_PBS_with_GOID[j],Filtered_PBS_without_GOID[j],fisher_table$BG_with_GOID[j],BG_without_GOID[j]),nrow=2,ncol=2),alternative="greater")$p.value) # maybe change to alternative="two.sided" ?
  }
  p_adj = round(p.adjust(pvalue, method = "BH"),3) # FDR p-value correction
  fisher_table$pvalue = pvalue
  fisher_table$p_adj_BH = p_adj
  fisher_table = fisher_table[order(fisher_table$p_adj_BH),] # sort by p_adj values (ascending)
  fisher_table_filtered = fisher_table[fisher_table$pvalue<0.05,] # save the enriched GO terms before applying FDR cutoff so that the user can apply his own on phylobrowse
  write.csv(fisher_table_filtered,paste(paste(i,sep="_"), "_all_enriched_terms_table.csv", sep=""), row.names = FALSE)
}