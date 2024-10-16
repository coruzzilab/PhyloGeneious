"""
Gil Eshel October 4 2017

This script will parse the jac_fixed.tre file with multiple trees, and calculates the node support values for the mpt_nel_fixed_labeled.tre tree
# synopsis : parse_jac_trees.py <jac_fixed.tre input> <mpt_nel_fixed_labeled.tre  input> <mpt_nel_labeled_with_jac_support_values.tre output>

"""

import sys

# sys.argv = ['script.py', 'jac_expressed_in_leaf_ovule_fixed.tre' , 'mpt_expressed_in_leaf_ovule_fixed_labeled.tre']


### Functions ###
def get_node_list(
    tree=None,
    node_label_list=None
):  # a function to get a list of node strings (e.g. ['(Grobu,Solyc)', '(Glyma,Teich)']) arguments should be the newick string and a node_label_list ("N#") - if the tree is not labeled - keep it empty...
    brack_starts_positions = [
        pos for pos, char in enumerate(tree) if char == "("
    ]  # get list of positions of "("
    brack_ends_positions = [
        pos for pos, char in enumerate(tree) if char == ")"
    ]  # get list of positions of ")"
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
    node_strings = []
    for i in range(
            len(any_brack_sort_binary)
    ):  # go through the 1 (=="(") and -1 (==")") list to get the positions where the accumulative sum goes back to 0 (track the positions in the any_brack_sort list to get the end position index in the tree string for that clade)
        numList = any_brack_sort_binary[
            i]  # if numList equal to 1 (= a starting position for a clade) if it is equal to "-1", then skip...
        if numList != 1:
            continue
        else:
            for j in range(
                    len(any_brack_sort_binary[i + 1:])
            ):  # accumulate the sum of the numbers from that index forward, until accumulative sum equal "0"
                numList = numList + any_brack_sort_binary[i + j + 1]
                if numList == 0:
                    end_position = i + j + 1
                    break  # once reached numList = 0, quite, add the clade string (the tree string between "i" and end_position), and move to the next "i"
                else:
                    continue
        node_strings.append(
            tree[any_brack_sort[i]:any_brack_sort[end_position]] +
            ")")  # Use this to find the internal branch lengths
    if node_label_list is not None:
        node_label_list_only_number_sorted = sorted(
            [int(s.replace('N', '')) for s in node_label_list],
            key=int,
            reverse=True)
        node_label_list_sorted = [
            'N' + str(s) for s in node_label_list_only_number_sorted
        ]
        node_strings_cleaned = []
        for s in node_strings:
            node_cleaned = str(s)
            for j in node_label_list_sorted:
                if ')' + j in node_cleaned:
                    node_cleaned = node_cleaned[:node_cleaned.find(')' + j) +
                                                1] + node_cleaned[
                                                    node_cleaned.find(')' +
                                                                      j) +
                                                    len(')' + j):]
            node_strings_cleaned.append(node_cleaned)
        return node_strings_cleaned
    else:
        return node_strings


def node_taxa_lists(
    node_string
):  # a function that take a node string and break it to two lists of the decedent clades - e.g.'(Grobu,Solyc)' to: [[Grobu],[Solyc]]
    cut = 0
    for i in range(len(list(node_string))):  # go position at a time, count
        if node_string[i] == "(":
            cut = cut + 1
        elif node_string[i] == ")":
            cut = cut - 1
        elif node_string[
                i] == ",":  # if the character is ",", and cut == 1, then cut the node into two parts/clades and store the taxa in a list of the two strings - order the taxa before converting back to a string
            if cut == 1:
                clade_1 = ','.join(
                    sorted(node_string[:i].replace("(",
                                                   "").replace(")",
                                                               "").split(",")))
                clade_2 = ','.join(
                    sorted(node_string[i + 1:].replace("(", "").replace(
                        ")", "").split(",")))
                node_descendants = [clade_1, clade_2]
                return node_descendants


### Input ###

with open(sys.argv[1], 'r') as infile1:
    jac = infile1.read().splitlines()

with open(sys.argv[2], 'r') as infile2:
    mpt_nel = infile2.read().replace('\n', '')

### Parse the mpt_nel_fixed_labeled.tre ###

# to (1) get the node labels and branch lengths stored for later, and (2) get the species lists of each side of each node
branch_lengths = mpt_nel.replace(";", "").replace(")", ",").replace(
    "(", ""
).split(
    ","
)  # currently a list of "label:branch_length" - for both internal and terminal nodes
branch_lengths_dict = {}
mpt_nel_no_branch_lengths = str(
    mpt_nel)  # get the mpt_nel_labeled without branch lengths
for i in branch_lengths:
    if ":" in i:
        label = i.split(":")[0]
        length = i.split(":")[1]
        mpt_nel_no_branch_lengths = mpt_nel_no_branch_lengths[
            0:mpt_nel_no_branch_lengths.find(label + ":") +
            len(label)] + mpt_nel_no_branch_lengths[
                mpt_nel_no_branch_lengths.find(label + ":") +
                len(label + ":") + len(length):]
        branch_lengths_dict[label] = length
    else:
        branch_lengths_dict[i] = ''

# get internal nodes dictionary (node ID:a_two_string_list of taxa from both sides) - need to match the N# id with its correct clade string - therefore I need to store the N# id from end to start...:
internal_node_labels_list = []  # get internal node labels
internal_node_labels_order_dict = {
}  # Important for later, need to keep track of the order (from the end) since I will need to sort the ids later, and then match back the clade strings to the correct node ids...
end_of_clade_positions = sorted(
    [i for i, ltr in enumerate(mpt_nel_no_branch_lengths) if ltr == ')'],
    key=int,
    reverse=True
)  # list of the positions of ')' in the tree (just before the N# label) - I sort from end to start...
posible_after_node_id = {
    ",", ":", ";", ")"
}  # to identify the position where N# ends - used below
count = 0
for i in end_of_clade_positions:  # loop the tree and extract the N# ids from end to start:
    end = next((i for i, ch in enumerate(mpt_nel_no_branch_lengths[i + 1:])
                if ch in posible_after_node_id), None)
    lab = mpt_nel_no_branch_lengths[i + 1:i + 1 + end]
    internal_node_labels_list.append(mpt_nel_no_branch_lengths[i + 1:i + 1 +
                                                               end])
    internal_node_labels_order_dict[lab] = count
    count = count + 1

# sort the list of node ids from largest to smallest numbers:
internal_node_labels_only_number_sorted = sorted(
    [int(s.replace('N', '')) for s in internal_node_labels_list],
    key=int,
    reverse=True)
internal_node_labels_sorted = [
    'N' + str(s) for s in internal_node_labels_only_number_sorted
]
internal_node_strings_list = get_node_list(mpt_nel_no_branch_lengths,
                                           internal_node_labels_sorted)

node_label_string_dict = {
}  # for each internal node store a list of two strings - of taxa from both sides of the fork - use this to compare the jac trees nodes to
for i in range(len(internal_node_strings_list)):
    node = 'N' + str(i + 1)
    node_label_string_dict[node] = node_taxa_lists(
        internal_node_strings_list[i])

### Parse the jac trees to get support values ###

total = 0
supported_nodes = [
]  # start a list, to add the node ids that have jac support - then count the occurrence of each node id
for i in jac[:
             -1]:  # ignore the last tree, which is the 50% majority rule tree from the jackknife run that is automatically added by TNT to the end of jac.tre file
    total = total + 1
    jac_internal_node_list = get_node_list(i)  # get internal node strings
    for j in jac_internal_node_list:
        jac_string_list = node_taxa_lists(j)
        for k in internal_node_labels_sorted:
            if set(jac_string_list) == set(node_label_string_dict[k]):
                supported_nodes.append(k)

node_support_percent = {
}  # used to be (but took to much time...) = {x:float("{0:.2f}".format(supported_nodes.count(x)*100./total)) for x in supported_nodes}
for x in internal_node_labels_sorted:
    node_support_percent[x] = "%.2f" % (supported_nodes.count(x) * 100. /
                                        total)

# If a node don't exist in the jac trees (wasn't recorded in the supported_nodes list), add a '0' as support value
for i in internal_node_labels_sorted:
    if i not in node_support_percent.keys():
        node_support_percent[i] = float(0)

### Build back the tree with node labels, jackknife support values and branch lengths ###

mpt_nel_with_support = str(mpt_nel)
for i in internal_node_labels_sorted:
    possible_matches = [
        ")" + i + ":", ")" + i + ";", ")" + i + ",", ")" + i + ")"
    ]
    for p in possible_matches:
        if p in mpt_nel_with_support:
            match = p
            start_insert = mpt_nel_with_support.find(
                match
            )  # get the position where to insert the ":" + support value
    add_support_str = match[:-1] + ":" + str(node_support_percent[i])
    mpt_nel_with_support = mpt_nel_with_support[:
                                                start_insert] + add_support_str + mpt_nel_with_support[
                                                    start_insert + len(match) -
                                                    1:]

### Save the tree with node labels, node supports and branch lengths ###

fout = open(sys.argv[3], 'w')  ### open output file for writting
fout.write(mpt_nel_with_support)  ### save
fout.close()

### Save the internal node labels and support values (tab-delimited, one node per line). For PhyloBrowse
with open('node_support_values.txt', 'w') as f:
    for key, value in sorted(node_support_percent.items()):
        f.write(str(key) + '\t' + str(value) + '\n')
