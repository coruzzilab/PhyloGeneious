""" 
Gil Eshel May 18 2017

Parse the partitionMembers.txt output file into a tabular (tab separated) file containing the following columns:
(1) Partition (gene) ID
(2) Family ID (gene family from which the partition/ortholog group was inferred)
(3) Partition_taxa_occupancy - the number of represented taxa
(4) Represented taxa list (same name as in the config file) - This column will be useful for pruning the TE tree for each partition, to include only the taxa that were aligned / have representative sequence...
(5) Ortholog sequence identifiers of each species
(6) Model species identifer - (the user can select one or more taxa as model species) all the identifiers are separated by ";", if none put "NA"  

-----------------------------------------------------------------------------------------------------
# input: partitionMembers.txt
# synopsis : parse_partitionMembers.py partitionMembers.txt <list of model species, same as in the config file separated with ",">

Patch notes:
Veronica Sondervan November 2022 - edit to accomodate paralogs in ortholog groups
"""
import sys

tabular_table = ["PartitionID","\t","FamilyID","\t","Partition_taxa_occupancy","\t","Taxa","\t","Ortholog_sequence_identifiers","\t","Model_species_identifier","\n"] # start a table to feed the information
if len(sys.argv) > 2: # Check if model species was/were suggested
	models = sys.argv[2].split(",") # retrieve the string of a list of model species from the command argument, and split it to a models list
else:
	models = 0 # a flag to check...

with open(sys.argv[1], 'r') as f:
    for line in f:
        row = [] # open an array to store the row information
        li=line.strip() # remove whitespace characters from the begining or end of the line, if exist...
        partition = li.split("=")[0].replace("p","").replace(" ","") # get the number of partition for that line
        row.append(partition) # add partition number
        
        FamilyID = li.split("\t")[1].replace("Family","") # get FamilyID (and remove the "Family" string to have only the family number)
        row.append(FamilyID) # add the FamilyID
        
        partition_identifiers = li.split(" = ")[1].split("\t")[0].replace(" ", ";") # isolate a string with all the identifiers and seperate them with ";"
        taxa = []
        for i in range(len(partition_identifiers.split(";"))):
        	taxa.append(partition_identifiers.split(";")[i].split("#")[0])
        #partition_taxa_occupancy = partition_identifiers.count('#') # This assumes that each species is represented with only one sequence
        taxa = reduce(lambda l, x: l if x in l else l+[x], taxa, [])
        partition_taxa_occupancy = len(taxa)-1
        if (partition_taxa_occupancy == 0):
            partition_taxa_occupancy = "NA"
        row.append(str(partition_taxa_occupancy))
        row.append(str(";".join((taxa))))
        row.append(str(partition_identifiers)) 
        
        #used_sequence_identifiers = [] # a list which will store the first identifier per taxa for this partition
        #for i in xrange(len(taxa)): #
        #	first_id = [s for s in all_species_identifier if taxa[i] in s][0]
        #	used_sequence_identifiers.append(first_id)       	
        #row.append(';'.join(sorted(set(used_sequence_identifiers), key=lambda x: used_sequence_identifiers.index(x))))

	      
        
        model_species_identifier = [] # open list to record the identifier of the model species - deal with the case where there is 0,1 or more models specifed 
        id = list()

        if models == 0:
            model_species_identifier = 'no_model_species_was_indicated'
        elif (len(models)==1):
            id = [s for s in partition_identifiers.split(";") if models[0] in s]
            if len(id) == 0:
                model_species_identifier = "NA"
            else:
                if len(id) !=0:
                    for j in range(len(id)):
                        ii = id[j].split("#")[1]
                        model_species_identifier.append(ii)
                if len(model_species_identifier) != 0:
                    model_species_identifier=';'.join(model_species_identifier)
                else:
                    model_species_identifier = "NA"
                
        else:
            for i in range(len(models)):
                id = ([s for s in partition_identifiers.split(";") if models[i] in s])
                if len(id) !=0:
                    for j in range(len(id)):
                        ii = id[j].split("#")[1]
                        
                        model_species_identifier.append(ii)
            if len(model_species_identifier) != 0:
                model_species_identifier=';'.join(model_species_identifier)
            else:
                model_species_identifier = "NA"

        row.append(str(model_species_identifier)) # add only the identifiers of the model species  
        tabular_table.append("\t".join((row)))         
        tabular_table.append("\n") # add end of line        
fout = open('parse_' + sys.argv[1], 'w')    ### open output file for writting
fout.write("".join((tabular_table)))   ### save 
fout.close()    
