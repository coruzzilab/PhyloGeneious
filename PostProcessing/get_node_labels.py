# Gil Eshel, Sep 30, 2017
# Script to add node labels ("N#") to the newick tree and to retreive node number to node/tip label table to be used when parsing the pbs.csv file
# # synopsis : get_node_labels.py <mpt_fixed.tre file> <mpt_nel_fixed.tre file>
# Output 1: node_number_to_label.txt - a table of tip/node number to tip/node label assostiation (for parsing the pbs.csv results)
# Output 2: mpt_nel_fixed_labeled.tre - the tree with the tip and node labels (short species names and the "N#" notation on the tree for PhyloBrowse) - use the labels_to_full_names.py script to replace species labels with full names
# Output 3: patristic_distances_all.csv - the patristic distances matrix between any node or tip in the tree (this can be used for normalizing the PBS values of more higher level nodes)

import sys, csv, re

### Input ###

with open(sys.argv[1], 'r') as infile1:
    mpt = infile1.read().replace('\n', '')
with open(sys.argv[2], 'r') as infile2:
    mpt_nel = infile2.read().replace('\n', '')

out_name = sys.argv[1].rsplit('.', 1)[0] + '_labeled.' + sys.argv[1].rsplit(
    '.', 1)[1]

### Get node labels ###

num = mpt.replace("(", "").replace(")", "").replace(";", "").split(
    ",")  # get a list of the terminal node species numbers
labels = [label.split(":")[0] for label in re.findall(r'\w+:', mpt_nel)]

# for later, to label the tree we will need a dictionary with num_to_labels
num_to_labels_dict = {}
for i in range(len(num)):
    num_to_labels_dict[num[i]] = labels[i]

# get the clades and add them to num and labels (N1, N2, ...)
brack_starts_positions = [pos for pos, char in enumerate(mpt)
                          if char == "("]  # get list of positions of "("
brack_ends_positions = [pos for pos, char in enumerate(mpt)
                        if char == ")"]  # get list of positions of ")"
any_brack = brack_starts_positions + brack_ends_positions  # include all positions and sort by position
any_brack_sort = sorted(any_brack)

any_brack_sort_binary = list(
    any_brack_sort
)  # get an corresponding list of "1" or "-1" for the position sorted list based on if it is "(" or ")", respectively
for i in range(len(any_brack_sort_binary)):
    if any_brack_sort_binary[i] in brack_starts_positions:
        any_brack_sort_binary[i] = 1
    else:
        any_brack_sort_binary[i] = -1

clades = []  # start a list to add the identified clades
clades_with_structure = [
]  # save also the clades but with structure for later extracting internal branch lengths
for i in range(
        len(any_brack_sort_binary)
):  # go through the 1 (=="(") and -1 (==")") list to get the positions where the accumulative sum goes back to 0 (track the positions in the any_brack_sort list to get the end position index in the tree string for that clade)
    numList = any_brack_sort_binary[
        i]  # if numList equal to 1 (= a starting position for a clade) if it is equal to "-1", then skip...
    #print numList
    if numList != 1:
        continue
    else:
        for j in range(
                len(any_brack_sort_binary[i + 1:])
        ):  # accumulate the sum of the numbers from that index forward, until accumulative sum equal "0"
            numList = numList + any_brack_sort_binary[i + j + 1]
            if numList == 0:
                end_position = i + j + 1
                break  # once reached numList = 0, quite, add the clade string (the mpt tree between "i" and end_position), and move to the next "i"
            else:
                continue
    clades.append(mpt[any_brack_sort[i]:any_brack_sort[end_position]].replace(
        "(", "").replace(")", "").replace(",", " "))
    clades_with_structure.append(
        mpt[any_brack_sort[i]:any_brack_sort[end_position]] +
        ")")  # Use this to find the internal branch lengths

clade_labels = []  # start a list to add the corresponding clades (N1, N2, ...)
for i in range(len(clades)):
    clade_labels.append("N" + str(i + 1))

# for later, to label the tree we will need a dictionary with clades_with_structure and their node labels ("N#")
clade_struct_label_dict = {}
for i in range(len(clades_with_structure)):
    clade_struct_label_dict[clades_with_structure[i]] = clade_labels[i]

all_num = ["number"] + num + clades
all_labels = ["label"] + labels + clade_labels

### Save the node_number_to_label file ###
with open('node_number_to_label.txt', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(all_num, all_labels))

### label the tree with node labels ###
terminal_branch_lengths = {}  # get a dictionary with terminal branch lengths
mpt_nel_no_terminal_branch_len = str(
    mpt_nel
)  # remove terminal branch lengths from the tree to retrieve internal branch lengths
# get terminal branch lengths
for i in labels:  # get the list of branch lengths for terminal nodes
    branch_length = mpt_nel.split(i + ":", 1)[1].split(",",
                                                       1)[0].split(")", 1)[0]
    terminal_branch_lengths[i] = branch_length
    mpt_nel_no_terminal_branch_len = mpt_nel_no_terminal_branch_len[
        0:mpt_nel_no_terminal_branch_len.
        index(":" + branch_length)] + mpt_nel_no_terminal_branch_len[
            mpt_nel_no_terminal_branch_len.index(":" + branch_length) +
            len(":" + branch_length):]

# replace label with number in the tree to compare with the list of the internal nodes (in numbers...)
for i in range(len(labels)):
    mpt_nel_no_terminal_branch_len = mpt_nel_no_terminal_branch_len.replace(
        labels[i], num[i])

# get the branch lengths of the internal nodes
internal_branch_lengths = {}  # get a dictionary with internal branch lengths
clades_with_structure_sorted = list(clades_with_structure)
#clades_with_structure_sorted.sort(lambda x,y: cmp(len(x), len(y))) # sort the clades_with_structure list from the smallest (in string length) to the largest
clades_with_structure_sorted.sort(
    key=len, reverse=True
)  # sort the clades_with_structure list from the smallest (in string length) to the largest
to_do_list = list(clades_with_structure_sorted)
mpt_nel_no_terminal_branch_len_tree = mpt_nel_no_terminal_branch_len
should_restart = True  # restart the loop until the to_do_list is empty...
while should_restart:
    should_restart = False  # only activate if list is not empty
    for i in range(len(clades_with_structure_sorted)):
        found_node = mpt_nel_no_terminal_branch_len_tree.find(
            clades_with_structure_sorted[i] + ":")
        if found_node != -1:
            int_branch_len = mpt_nel_no_terminal_branch_len_tree[
                found_node:].split(":", 1)[1].split(",", 1)[0].split(")", 1)[0]
            internal_branch_lengths[clades_with_structure_sorted[
                i]] = int_branch_len  # add branch length to the internal_branch_lengths
            mpt_nel_no_terminal_branch_len_tree = mpt_nel_no_terminal_branch_len_tree[
                0:found_node +
                len(clades_with_structure_sorted[i]
                    )] + mpt_nel_no_terminal_branch_len_tree[
                        found_node + len(clades_with_structure_sorted[i]) + 1 +
                        len(int_branch_len):]
            del to_do_list[to_do_list.index(clades_with_structure_sorted[i])]
            if not to_do_list:
                break
            else:
                should_restart = True
        else:
            if not to_do_list:
                should_restart = False

# generate a tree with node labels and branch lengths
labeled_tree = str(mpt)  # label using the mpt tree template (numbers)
order_of_clades_to_add_info = list(
    clades_with_structure
)  # list that will dictate the order of adding the node label and branch length information
order_of_clades_to_add_info.sort(
    key=len, reverse=True)  # this time we what the larger clade to be first
#order_of_clades_to_add_info.sort(lambda x,y: cmp(len(y), len(x))) # this time we what the larger clade to be first

# label internal nodes:
for i in order_of_clades_to_add_info:
    if i in list(internal_branch_lengths.keys()):
        labeled_node = i + clade_struct_label_dict[
            i] + ":" + internal_branch_lengths[i]
        labeled_tree = labeled_tree.replace(i, labeled_node)
    else:  # in the case of the largest clade where the is no branch length info...
        labeled_node = i + clade_struct_label_dict[i]
        labeled_tree = labeled_tree.replace(i, labeled_node)

# label terminal nodes:
for i in num:
    label_of_taxa = num_to_labels_dict[i]  # get the label
    labeling = label_of_taxa + ":" + terminal_branch_lengths[
        label_of_taxa]  # generate the label with the branch length
    if "(" + i + "," in labeled_tree:
        labeled_tree = labeled_tree.replace("(" + i + ",",
                                            "(" + labeling + ",")
    elif "," + i + ")" in labeled_tree:
        labeled_tree = labeled_tree.replace("," + i + ")",
                                            "," + labeling + ")")

### Save the labeled tree ###
fout = open(out_name, 'w')  ### open output file for writting
fout.write(labeled_tree + '\n')  ### save
fout.close()

### calculate the patristic distances between any pair of nodes (both internal and terminal nodes)###

# convert the num to labels for the list of clades:
clades_with_species_labels = list(clades_with_structure)
for i in range(len(clades_with_species_labels)):
    for j in num_to_labels_dict:
        if "(" + j + "," in clades_with_species_labels[i]:
            clades_with_species_labels[i] = clades_with_species_labels[
                i].replace("(" + j + ",", "(" + num_to_labels_dict[j] + ",")
        elif "," + j + ")" in clades_with_species_labels[i]:
            clades_with_species_labels[i] = clades_with_species_labels[
                i].replace("," + j + ")", "," + num_to_labels_dict[j] + ")")

all_branch_lengths = internal_branch_lengths.copy(
)  # copy the internal_branch_lengths to change the numbers to labels
for i in list(all_branch_lengths):
    try:
        ind = clades_with_structure.index(i)
        all_branch_lengths[
            clades_with_species_labels[ind]] = all_branch_lengths[i]
        del all_branch_lengths[i]
    except:
        continue
all_branch_lengths.update(
    terminal_branch_lengths)  # add the terminal_branch_lengths
all_branch_lengths.update(
    {clades_with_species_labels[0]: '0'}
)  # add the largest node (that wasn't included and left in the to_do_list) to all_branch_lengths with '' as branch length

# calculate the patristic distances:
all_nodes = labels + clades_with_species_labels

# create a patristic distances matrix
w, h = len(all_nodes), len(all_nodes)
# set dimensions
all_nodes_patristic_matrix = [[0 for x in range(w)] for y in range(h)]
for i in all_nodes:
    for j in all_nodes:
        i_nodes = [x for x in all_branch_lengths
                   if i in x]  # all nodes that include node i
        j_nodes = [x for x in all_branch_lengths
                   if j in x]  # all nodes that include node j
        list_for_calculating_dist = list(
            set(i_nodes).symmetric_difference(j_nodes)
        )  # a list that contain only one of the two species (= without the common ancestors) for calculating the patristic distance (sum of the branch lengths)
        patristic_dist = sum(
            [int(all_branch_lengths[x]) for x in list_for_calculating_dist])
        all_nodes_patristic_matrix[all_nodes.index(i)][all_nodes.index(
            j)] = patristic_dist  # assign the patristic distance values

# add row labels to the matrix:
all_labels = labels + clade_labels
for i in range(len(all_nodes_patristic_matrix)):
    all_nodes_patristic_matrix[i].insert(0, all_labels[i])

# add column labels to the matrix:
all_nodes_patristic_matrix.insert(0, all_labels)  # insert the column headers
all_nodes_patristic_matrix[0].insert(0, '')  # insert the column headers

### Save patristic distances matrix ###
with open("patristic_distances_all.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerows(all_nodes_patristic_matrix)
