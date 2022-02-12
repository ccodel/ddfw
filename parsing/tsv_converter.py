# @file tsv_converter.py
# @brief Converts output from smart_converter.py and starexec.py into groups
#        of TSV files for reading by gnuplot.
#
# @usage python tsv_converter.py <file> <x> <y> <z> [split cols]
#
# @author Cayden Codel
# @bug No known bugs.
#
# @todos validation logic, e.g. file existing, etc.
#

import sys
import os.path

def print_data_helper(d, s):
    for k, v in d.items():
        if isinstance(v, list):
            if len(s) == 0:
                f = open(k + ".tsv", "w+")
            else:
                f = open(s + "_{0}.tsv".format(k), "w+")
            for (x, y, z) in v:
                f.write("{0}\t{1}\t{2}\n".format(str(x), str(y), str(z)))
        else:
            if len(s) == 0:
                print_data_helper(v, k)
            else:
                print_data_helper(v, s + "_{0}".format(k))

def print_data(d):
    print_data_helper(d, "")

if len(sys.argv) < 5:
    print("Usage is tsv_converter.py <file> <x> <y> <z> [split cols]")
    exit(1)

file_name = sys.argv[1]
x_col = sys.argv[2]
y_col = sys.argv[3]
z_col = sys.argv[4]

# Open the file
f = open(file_name, "r")
lines = f.readlines()
f.close()

# Extract the header, find out which cols x, y, z correspond to
header = lines[0].strip().split(",")
if x_col not in header or y_col not in header or z_col not in header:
    print("Column not in header")
    exit(1)

x = header.index(x_col)
y = header.index(y_col)
z = header.index(z_col)

# Find out the split columns, too
split = []
if len(sys.argv) > 5:
    split = sys.argv[5:]
    for s in split:
        if s not in header:
            print("Split not in header")
            exit(1)
    
    split = list(map(lambda x: header.index(x), split))

# For use in split
data = {}

lines = list(map(lambda x: x.strip().split(","), lines[1:]))
for line in lines:
    xd = line[x]
    yd = line[y]
    zd = line[z]

    if len(split) > 0:
        d = data
        for i in range(0, len(split) - 1):
            s = split[i]
            sd = line[s]
            if d.get(sd) is None:
                d[sd] = {}
            d = d[sd]

        sd = line[split[-1]]
        if d.get(sd) is None:
            d[sd] = [(xd, yd, zd)]
        else:
            d[sd].append((xd, yd, zd))
    else:
        print("{0}\t{1}\t{2}".format(str(xd), str(yd), str(zd)))

if len(split) > 0:
    # For each array of elems, create a new file and print list
    print_data(data)
