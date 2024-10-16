"""
Gil Eshel June 2 2017

This script will replace, for any input text file, the short species labels with their full names based on a label2fullname tab-seperated file (e.g. species.txt)
# synopsis : labels_to_full_names.py <label2fullname_file> <input_file> <output_file>

"""
import sys


# a function to replace the short with the long name
def replace_all(text, dic):
    #for i, j in dic.iteritems():
    for i, j in dic.items():
        text = text.replace(i, j)
    return text


# parse the label2fullname_file into a dictionary
f = open(sys.argv[1], 'r')  # open the label2fullname_file for reading
nameDict = {
}  # create a dictionary to store the short and long names to be replaced
for line in f:
    short = line.split("\t")[0].strip()
    full = line.split("\t")[1].strip()
    nameDict[short] = full
f.close()  # close the label2fullname_file

# Replace short with full taxa names
f = open(sys.argv[2], 'r')  # open input_file
text = f.read()  # read the file into a string
f.close()  # close the input_file
text_with_full_names = replace_all(
    text, nameDict)  # apply the function and convert the names

# Save the text_with_full_names into the output file
fout = open(sys.argv[3],
            'w')  ### open output file for writing the text_with_full_names
fout.write(text_with_full_names)  ### save
fout.close()
