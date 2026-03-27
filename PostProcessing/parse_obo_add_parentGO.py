#!/usr/bin/env python
#
# By Gil Eshel, June 22, 2020
#
# A script to include parent GO terms (is_a and part_of) to a provided gene2go annotation tab-delimited gene_id<\t>GO_id_list_comma_separated) - for Fisher Exact test
# Takes as input the gene2go.txt (tab-delimited gene_id<\t>GO_id_list_comma_separated), downloads the latest go-basic.obo file (http://geneontology.org/page/download-ontology) with the GO parent / alternative information, parse it and add parent GO terms to the gene2go annotation file (keep a unique set of GO_ids for each gene_id...) - this output will be use fot the GO enrichment analysis.
# I adopted the obo file parsing code from Damian Kao (http://blog.nextgenetics.net/?e=6) - I modified it to include both is_a and part_of as parent GO_ids. Thank you Damian Kao for a great post!
#
# synopsis : parse_obo_add_parentGO.py <gene2go.txt> <gene2go_with_parents.txt>

import sys, subprocess, itertools, os # sys is for retrieving the script arguments, subprocess used for executing bash commands from within the python script, itertools is for grouping the line of each GO term in the go-basic.obo file together (separated by a blank line), os is for checking if the go-basic.obo file exists already in the working DIR

### get most recent go-basic.obo file
if not os.path.exists('go-basic.obo'):
	try:
		subprocess.call(" ".join(["wget", "http://purl.obolibrary.org/obo/go/go-basic.obo"]), shell=True)
	except:
		print("go-basic.obo could not be download")
	
	print('downloaded obo file - yeee')

### parse the go-basic.obo to retrieve parent GO terms:
# few parsing functions:
def getTerm(stream):
	block = []
	for line in stream:
		if line.strip() == "[Term]" or line.strip() == "[Typedef]":
			break
		else:
			if line.strip() != "":
				block.append(line.strip())
	return(block)

print('loaded getTerm function')

def parseTagValue(term):
	data = {}
	for line in term:
		if line.startswith('relationship: part_of'):
			tag = line.split(' ')[1] # tag = 'part_of'
			value = line.split(' ')[2] # collect the GO_id
		elif line.startswith('is_a'):
			tag = line.split(': ',1)[0] # tag = 'is_a'
			value = line.split(': ',1)[1].split()[0] # collect the GO_id
		else:
			tag = line.split(': ',1)[0] # other tags like 'id:', 'name:', 'namespace:' etc. are collected
			value = line.split(': ',1)[1]
		if tag not in data:
			data[tag] = []
		data[tag].append(value)
	return(data)

print('loaded parseTagValue function')

# A recursive function that take a GO term id and look for all parents (full path until the root) in the terms dictionary (created when we parsed the go-basic.obo file).
# Return a list of parent ids for the input GO id:
def find_parents(go_id, terms_dict, go_term_set=[], ret=False):
	for term2 in terms_dict[go_id]['p']:
		#print term2
		go_term_set.append(term2)
		# Recurse on term to find all parents
		if term2 == 'GO:0003674' or term2 == 'GO:0008150' or term2 == 'GO:0005575':
			ret=True
		else:
			find_parents(term2, terms_dict, go_term_set)          
	return(go_term_set)

print('loaded find_parents function')


oboFile = open('go-basic.obo','r')
print('opened oboFile - yeee!')

terms = {} # declare a blank dictionary - keys are the GOids and the values will be their parent and Child GO_ids
getTerm(oboFile) #skip the file header lines
go_names_dict = {}	# to store the 'GO_id<\t>GO_term' info
go_cat_dict = {}	# to store the GO category 'GO_id<\t>GO_category' info (GO_categories: biological_process (BP), molecular_function (MF), cellular_component (CC))
while 1: #infinite loop to go through the obo file. Breaks when the term returned is empty, indicating end of file
	term = parseTagValue(getTerm(oboFile)) #get the term using the two parsing functions
	if len(term) != 0:
		termID = term['id'][0]
		multi_termID = [termID]
		if 'alt_id' in term:
			multi_termID.extend(term['alt_id'])
		termParents = []
		for i in list(term):	# loop through the term fields
			if i == 'is_a' or i == 'part_of':
				termParents.extend(term[i])
		for t in multi_termID:
			if t not in terms:
				terms[t] = {'p':[],'c':[]} # each GO_id will have two arrays of parents and children
			terms[t]['p'] = termParents # append parents of the current term
			for j in termParents: #for every parent term, add this current term as children
				if j not in terms:
					terms[j] = {'p':[],'c':[]}
				terms[j]['c'].append(t)
		#elif term.has_key('name'):
		for g in multi_termID:
			if g.startswith('GO:'):
				if g not in go_names_dict:
					go_names_dict[g] = ''.join(g) + '\t' + ''.join(term['name'])	# store the 'GO_id<\t>GO_term' info, e.g. 'GO:0000008<\t>obsolete thioredoxin' 
					go_cat_dict[g] = ''.join(g) + '\t' + ''.join(term['namespace'])	# store the 'GO_id<\t>GO_category' info, e.g. 'GO:0000008<\t>molecular_function' 
	else:
		break

oboFile.close()
print('parsed obo file - yeee')

# For each GO term, get a list of all parent GO terms (until the root):
all_parent_terms = {}
for g in terms:
	if len(terms[g]['p']) == 0:
		all_parent_terms[g] = ''
	elif g == 'negatively_regulates' or g == 'positively_regulates':
		continue
	else:
		all_parent_terms[g] = list(set(find_parents(g,terms,go_term_set=[], ret=False)))

### Parse the user provided gene2go.txt file to add the parent GO_ids:
gene_to_go = [] # a list to store the gene<\t>GO_all<\n> strings
with open(sys.argv[1], 'r') as goFile:
	for line in goFile:
		li=line.strip()
		gene = li.split('\t')[0]
		GO = li.split('\t')[1].split(',')
		GO_parents = []
		for go in GO:
			if go in all_parent_terms:
				GO_parents.extend(all_parent_terms[go])
		GO_all = ','.join(list(set(','.join(GO).split(',') + GO_parents))) # merge the GO list (original) with the GO_parents list to a unique list of GO_ids - join them to a comma-separated string
		gene_to_go.append(gene + '\t' + GO_all)		

print('parsed gene2go file - yeee!')

### open an output file to print into the new annotation file with the parent GO_ids included
with open(sys.argv[2], 'w') as output:
    output.write('\n'.join(gene_to_go) + '\n')

### open an output file to print into the GO_id to GO description file
og_desc = ''
for d in go_names_dict:
	og_desc = og_desc + go_names_dict[d] + '\n'

with open(sys.argv[2].rsplit('.',1)[0] + '_GoDesc.txt', 'w') as output:
    output.write(og_desc)

### open an output file to print into the GO_id to GO category file
og_cat = ''
for c in go_cat_dict:
	og_cat = og_cat + go_cat_dict[c] + '\n'

with open(sys.argv[2].rsplit('.',1)[0] + '_GOCategory.txt', 'w') as output:
    output.write(og_cat)
