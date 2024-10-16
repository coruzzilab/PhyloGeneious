"""
Gil Eshel June 5 2017

This script parse a TNT "tre" file output tree into a proper newick format (the ".nel" tree is a special case - use the fix_nel_tree.py for that file)
# input: any TNT.tre
# output: fixed TNT.tre file name
# synopsis : fix_tree.py TNT.tre TNT_fixed.tre

"""
import sys

# Parse the "tre" file
f = open(sys.argv[1], 'r')  # open TNT tree file for reading
tre_lines = f.readlines(
)  # read the lines into a list (essentially three lines...)
f.close()  # close the TNT tree file
tre_lines.pop(
    0
)  # remove the first line in the TNT file (something like: "tread 'tree(s) from TNT, for data in matrix.tnt'")
tre_lines = tre_lines[:-1]  # remove the last line in the TNT file ("proc-;")

# each of the "jac.tre", "mpt.tre" and "mpt.nel" output files are a bit different, therefore need to deal with each case

if len(tre_lines) == 1:  # if only one tree - it is either mpt.nel or mpt.tre
    if '=' in tre_lines[
            0]:  # if contains a "=" character, it means that it is a mpt.nel file
        fixed_tre = tre_lines[0].replace(
            " =", ":"
        ).replace(" )=", "):").replace(" ", ",").replace(",)", ")").replace(
            "@OUTGROUP_", ""
        ).replace(
            "@INGROUP_", ""
        )  # parse the ".nel" file tree into a proper newick format (comma seperation, and replace "=" with ":")
    else:  # it is a mpt.tre file or the branch lengths are already fixed to have ":" and not "="
        fixed_tre = tre_lines[0].replace(" ", ",").replace(",)", ")").replace(
            ")(",
            "),(")  # parse the tre to proper newick format (comma seperation)
else:  # it has more than one tree - it is a jac.tre file
    fixed_tre = []
    for tre in tre_lines:
        fixed_tre.append(
            tre.replace(" ",
                        ",").replace(",)", ")").replace(")(", "),(").replace(
                            "*", ";").replace("@OUTGROUP_",
                                              "").replace("@INGROUP_", "")
        )  # parse the tre to proper newick format (comma seperation)
    fixed_tre = ''.join(fixed_tre)

fout = open(sys.argv[2], 'w')  ### open output file for writing the fixed tree
fout.write(fixed_tre)  ### save
fout.close()
