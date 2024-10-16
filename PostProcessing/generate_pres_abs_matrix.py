# Gil Eshel, Oct 23, 2018
# Script to generate presence absence matrix for heatmap visualization

# # synopsis : generate_pres_abs_matrix.py parse_partitionMembers.txt species.txt <ortho_presence_absence_matrix.txt output name>

import sys

#sys.argv = ['generate_pres_abs_matrix.py','parse_partitionMembers.txt', 'species.txt' ,'ortho_presence_absence_matrix.txt']

### Parse input files ###

lab_to_num = {}	# dictionary to contain the species name per species label (e.g. Aratha:'Arabidopsis thaliana') from the species.txt file
all_taxa_labels = []	# list of all the taxa labels (needed for identifying if they are missing)
all_taxa_names = []	# list of all the taxa names for the table column headers
with open(sys.argv[2], 'r') as species_file:
	for line in species_file:
		lab_to_num[line.strip().split('\t')[0]] = line.strip().split('\t')[1]
		all_taxa_labels.append(line.strip().split('\t')[0])
		all_taxa_names.append(line.strip().split('\t')[1])

with open(sys.argv[1], 'r') as partMem_file:
	partition_info = partMem_file.read().splitlines()

#partition_taxa_dict = {} # store the list of taxa for each partition (use this later to assess the taxa coverage within and outside the clades)
partition_presence_absence_dict = {} # store for each partition the presence ("1") or absence ("0") of each taxon
all_partition_ids = [] # store all partition ids for later
for i in partition_info[1:]:
	id = i.split('\t')[0]
	all_partition_ids.append(id)
	partition_taxa_list = ' '.join(i.split('\t')[3].split(';')).split()
	#partition_taxa_dict[id] = partition_taxa_list
	taxa_presence_absence = {}
	for t in all_taxa_labels:
		if t in partition_taxa_list:
			taxa_presence_absence[lab_to_num[t]] = "1"
		else:
			taxa_presence_absence[lab_to_num[t]] = "0"
	partition_presence_absence_dict[id] = taxa_presence_absence

### Generate a presence/absence table ###

presence_absence_all_partitions = ["PartitionID\t" + "\t".join(all_taxa_names)] # start a presence/absence table of all partitions (this is the header line)
for p in all_partition_ids: # loop through all the partitions
	get_taxa_presence_absence = partition_presence_absence_dict[p]
	taxa_pres_abs_list = [p] # start a list with the partition id
	for t in all_taxa_names: # for each taxa append the presence/absence status to the taxa_pres_abs_list
		taxa_pres_abs_list.append(get_taxa_presence_absence[t])
	presence_absence_all_partitions.append("\t".join(taxa_pres_abs_list))

### Save the presence/absence tables ###
f_all_table = open(sys.argv[3], 'w')
f_all_table.write("\n".join(presence_absence_all_partitions))
f_all_table.close()

