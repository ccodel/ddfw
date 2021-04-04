# @file smart_converter.py
# @brief Converts DDFW into CSV output with many options for data aggregation.
#
# @usage python starexec_converter.py [options]
#   -a
#   -c <column list> [...]
#   -f <input_file>
#   -m <method>
#   -o <output_file>
#   -p <list_file>
#   -r <dir>
#
# @author Cayden Codel (ccodel@andrew.cmu.edu)
#
# @bug No known bugs.

from __future__ import print_function
import sys
import string
from string import ascii_letters
from functools import reduce
import os.path
from os import walk
from optparse import OptionParser

class DDFWParser:
    all_normal_cols = ["instance", "ddfw_version", "a", "A", "c", "C", "w", 
                "g", "G", "o", "cpu_timeout", "flip_timeout", "method", "seed",
                "solve", "flips", "best", "best_step", "time"]
    normal_input_cols = ["instance", "ddfw_version", "a", "A", "c", "C", "w",
            "g", "G", "o", "cpu_timeout", "flip_timeout", "method", "seed"]
    normal_output_cols = ["solve", "flips", "best", "best_step", "time"]

    all_weight_cols = ["instance", "ddfw_version", "a", "A", "c", "C", "w",
            "g", "G", "o", "cpu_timeout", "flip_timeout", "method", "seed", 
            "flips", "transfers", "unsat_weight", "unsat_count",
            "best_unsat_count", "avg_transfer", "min_unsat_weight", 
            "max_unsat_weight", "min_sat_weight", "max_sat_weight"]
    weight_input_cols = ["instance", "ddfw_version", "a", "A", "c", "C", "w",
            "g", "G", "o", "cpu_timeout", "flip_timeout", "method", "seed", 
            "flips"]
    weight_output_cols = ["transfers", "unsat_weight", "unsat_count",
            "best_unsat_count", "avg_transfer", "min_unsat_weight", 
            "max_unsat_weight", "min_sat_weight", "max_sat_weight"]

    all_cols = []
    input_cols = []
    output_cols = []

    normal_method = True
    fixed_flip = None
    config = {}
    data_store = None

    def __init__(self, method="normal", aggregate=False, columns=None, fixed_flip=None):
        if method == "normal":
            self.normal_method = True
            self.all_cols = self.all_normal_cols
            self.input_cols = self.normal_input_cols
            self.output_cols = self.normal_output_cols
        elif method == "weight":
            self.normal_method = False
            if fixed_flip is not None:
                self.fixed_flip = int(fixed_flip)
            self.all_cols = self.all_weight_cols
            self.input_cols = self.weight_input_cols
            self.output_cols = self.weight_output_cols
        else:
            print("Unrecognized method\n", sys.stderr)
            exit()

        self.aggregate = aggregate

        if columns is not None:
            ic = self.input_cols
            oc = self.output_cols

            self.input_cols = [s for s in self.input_cols if s in columns]
            self.output_cols = [s for s in self.output_cols if s in columns]

            if len(self.input_cols) == 0:
                self.input_cols = ic

            if len(self.output_cols) == 0:
                self.output_cols = oc

            self.all_cols = self.input_cols + self.output_cols
        
        if self.aggregate:
            self.all_cols.append("trials")
            self.data_store = {}

    def parse_line(self, line, outfile, path=None):
        # Disregard empty lines or those that aren't comment lines
        if len(line) == 0:
            return

        if "bce" in path:
            return

        if "Alarm" in line:
            return

        # Chop off the first 'c ' part of the string, split into words
        words = line[line.find("c") + 2:].strip().split(" ")
        if len(words) <= 1:
            return

        # Extract configuration and run statistics
        if "Version" == words[0]:
            self.config["ddfw_version"] = int(words[1])
        elif "-f" == words[0]:
            if "theBenchmark" in words[1]:
                self.config["instance"] = path.split("/")[-2]
            else:
                self.config["instance"] = words[1]
        elif "-m" == words[0]:
            self.config["method"] = words[1]
        elif "-w" == words[0]:
            self.config["w"] = float(words[1])
        elif "-a" == words[0]:
            self.config["a"] = float(words[1])
        elif "-A" == words[0]:
            self.config["A"] = float(words[1])
        elif "-c" == words[0]:
            self.config["c"] = float(words[1])
        elif "-C" == words[0]:
            self.config["C"] = float(words[1])
        elif "-g" == words[0]:
            self.config["g"] = words[1]
        elif "-G" == words[0]:
            self.config["G"] = words[1]
        elif "-o" == words[0]:
            self.config["o"] = words[1]
        elif "-s" == words[0]:
            self.config["seed"] = int(words[1])
        elif "-t" == words[0]:
            if not "infty" == words[1]:
                self.config["cpu_timeout"] = int(words[1])
            else:
                self.config["cpu_timeout"] = -1
        elif '-T' == words[0]:
            if not "infty" == words[1]:
                self.config["flip_timeout"] = int(words[1])
            else:
                self.config["flip_timeout"] = -1
        elif "Stats:" == words[0] and self.normal_method:
            stats = line.replace(" ", "").replace(":", "").strip().strip(ascii_letters).split("|")
            if len(stats) < 6:
                return
            
            solve = int(stats[1])
            flips = int(stats[2])
            best = int(stats[3])
            best_step = int(stats[4])
            time = float(stats[5])

            if self.aggregate:
                # Store based on those input columns specified
                l = self.data_store
                for c in self.input_cols[:-1]:
                    if l.get(self.config[c]) is None:
                        l[self.config[c]] = {}
                    l = l[self.config[c]]

                c = self.input_cols[-1]
                if l.get(self.config[c]) is None:
                    l[self.config[c]] = []
                l = l[self.config[c]]

                # Now go through each output column to construct a list to add
                to_add = []
                
                if "solve" in self.output_cols:
                    to_add.append(solve)
                if "flips" in self.output_cols:
                    to_add.append(flips)
                if "best" in self.output_cols:
                    to_add.append(best)
                if "best_step" in self.output_cols:
                    to_add.append(best_step)
                if "time" in self.output_cols:
                    to_add.append(time)

                l.append(to_add)
            else:
                string_comps = []
                for column in self.input_cols:
                    if self.config.get(column) is not None:
                        string_comps.append(str(self.config[column]))
                
                if "solve" in self.output_cols:
                    string_comps.append(str(solve))
                if "flips" in self.output_cols:
                    string_comps.append(str(flips))
                if "best" in self.output_cols:
                    string_comps.append(str(best))
                if "best_step" in self.output_cols:
                    string_comps.append(str(best_step))
                if "time" in self.output_cols:
                    string_comps.append(str(time))

                string = ",".join(string_comps)
                string += "\n"
                outfile.write(string)
        elif "In-stats:" == words[0] and not self.normal_method:
            stats = line.replace(" ", "").replace(":", "").strip().strip(ascii_letters).split("|")
            if len(stats) < 11:
                return
          
            self.config["flips"] = int(stats[1])
            transfers = int(stats[2])
            unsat_count = int(stats[3])
            best_unsat_count = int(stats[4])
            unsat_weight = float(stats[5])
            avg_transfer = float(stats[6])
            min_unsat_weight = float(stats[7])
            max_unsat_weight = float(stats[8])
            min_sat_weight = float(stats[9])
            if len(stats[10]) == 0:
                return
            max_sat_weight = float(stats[10])

            if self.fixed_flip is not None and self.fixed_flip != int(stats[1]):
                return

            if self.aggregate:
                # Store based on those input columns specified
                l = self.data_store
                for c in self.input_cols[:-1]:
                    if l.get(self.config[c]) is None:
                        l[self.config[c]] = {}
                    l = l[self.config[c]]
                
                c = self.input_cols[-1]
                if l.get(self.config[c]) is None:
                    l[self.config[c]] = []
                l = l[self.config[c]]

                # Now go through each output column to construct a list to add
                to_add = []
                if "transfers" in self.output_cols:
                    to_add.append(transfers)
                if "unsat_weight" in self.output_cols:
                    to_add.append(unsat_weight)
                if "unsat_count" in self.output_cols:
                    to_add.append(unsat_count)
                if "best_unsat_count" in self.output_cols:
                    to_add.append(best_unsat_count)
                if "avg_transfer" in self.output_cols:
                    to_add.append(avg_transfer)
                if "min_unsat_weight" in self.output_cols:
                    to_add.append(min_unsat_weight)
                if "max_unsat_weight" in self.output_cols:
                    to_add.append(max_unsat_weight)
                if "min_sat_weight" in self.output_cols:
                    to_add.append(min_sat_weight)
                if "max_sat_weight" in self.output_cols:
                    to_add.append(max_sat_weight)

                l.append(to_add)
            else:
                string_comps = []
                for key, value in self.config.items():
                    if key in self.input_cols:
                        string_comps.append(str(value))

                # Emit only those columns specified
                if "transfers" in self.output_cols:
                    string_comps.append(str(transfers))
                if "unsat_weight" in self.output_cols:
                    string_comps.append(str(unsat_weight))
                if "unsat_count" in self.output_cols:
                    string_comps.append(str(unsat_count))
                if "best_unsat_count" in self.output_cols:
                    string_comps.append(str(best_unsat_count))
                if "avg_transfer" in self.output_cols:
                    string_comps.append(str(avg_transfer))
                if "min_unsat_weight" in self.output_cols:
                    string_comps.append(str(min_unsat_weight))
                if "max_unsat_weight" in self.output_cols:
                    string_comps.append(str(max_unsat_weight))
                if "min_sat_weight" in self.output_cols:
                    string_comps.append(str(min_sat_weight))
                if "max_sat_weight" in self.output_cols:
                    string_comps.append(str(max_sat_weight))

                string = ",".join(string_comps)
                string += "\n"
                outfile.write(string)


    # Recursive function, build string through levels of dictionary
    def aggregate_helper(self, outfile, s, d):
        if isinstance(d, list):
            # Average out all entries
            length = len(d)
            summed = reduce(lambda x, y: [a + b for a, b in zip(x, y)], d)
            averaged = map(lambda x: float(x) / float(length), summed)
            # averaged = map(lambda x: float(x) / float(length), summed)
            stringed = map(lambda x: str(x), averaged)
            outfile.write(s + ",".join(stringed) + "," + str(length) + "\n")
        else:
            items = sorted(d.items())
            for key, val in items:
                self.aggregate_helper(outfile, s + str(key) + ",", val)

    def print_aggregate(self, outfile):
        if not self.aggregate:
            return
        self.aggregate_helper(outfile, "", self.data_store)

###############################################################################
# Main execution
###############################################################################

# Recursive function for calculating input files
def calc_files(dir_str):
    l = []
    if os.path.isdir(dir_str):
        _, dir_names, files = next(walk(dir_str))
        for f in files:
            if "DS_Store" not in f:
                l.append(dir_str + "/" + f)
        for d in dir_names:
            res = calc_files(dir_str + "/" + d)
            l.extend(res)
    return l

# Set up command line options
clparser = OptionParser()
clparser.add_option("-a", "--aggregate", dest="aggregate", action="store_true",
                    default=False, help="Parse aggregate averaging data")
clparser.add_option("-c", "--columns", dest="columns", action="append",
                    help="Columns from the DDFW output to take")
clparser.add_option("-f", "--infile", dest="infile",
                    help="Specify an input file. -f and -r not at same time.")
clparser.add_option("-F", "--fixed_flip", dest="fixed_flip",
                    help="Set a fixed flip to extract")
clparser.add_option("-l", "--list_columns", dest="list_columns", 
                    help="List the possible columns.", action="store_true")
clparser.add_option("-m", "--method", dest="method", default="normal",
                    help="Method of parsing, depending on data to extract")
clparser.add_option("-o", "--outfile", dest="outfile",
                    help="Specify a CSV outfile. stdout is default")
clparser.add_option("-p", "--parsed_list", dest="parsed_list",
                    help="List of parsed files, blacklist")
clparser.add_option("-r", "--directory", dest="directory",
                    help="Directory from which to pull DDFW output files")

(options, args) = clparser.parse_args()

# Parsing of CLE args
input_files = None
parsed_list_file = None
outfile = None

# Create a parser depending on the provided parsing method
if options.method is not None:
    parser = DDFWParser(aggregate=options.aggregate,
                        method=options.method,
                        columns=options.columns,
                        fixed_flip=options.fixed_flip)
else:
    print("No method supplied", file=sys.stderr)
    exit()

# If the user would like to have the columns listed, list them here and exit
if options.list_columns:
    print("\nColumns for normal method are:\n" + str(parser.all_normal_cols))
    print("\nColumns for weight method are:\n" + str(parser.all_weight_cols))
    print("\n")
    exit()

# Take an infile or a directory
if options.infile is not None and options.directory is not None:
    print("Both infile and directory options cannot be supplied together",
          file=sys.stderr)
    exit()
elif options.infile is None and options.directory is None:
    print("Must supply an infile with -f or -r", file=sys.stderr)
    exit()
elif options.infile is not None:
    input_files = [options.infile]
else:
    input_files = calc_files(options.directory)

# Map every file onto its absolute path
input_files = sorted(list(map(lambda x: os.path.abspath(x), input_files)))

# If parsed list provided, pull out lines and remove from input files
# TODO The parsed list cannot contain relative directories, change later
#   e.g. using os.path.abspath()
if options.parsed_list is not None:
    if os.path.isfile(options.parsed_list):
        parsed_list_file = open(options.parsed_list, "r")
        parsed_list_lines = parsed_list_file.readlines()
        parsed_list_lines = list(map(lambda x: x.strip(), parsed_list_lines))
        input_files = [x for x in input_files if x not in parsed_list_lines]
        parsed_list_file.close()
    
    parsed_list_file = open(options.parsed_list, "a+")

# Check existence of CSV output file
if options.outfile is not None:
    if not os.path.isfile(options.outfile) and os.path.exists(options.outfile):
        print("File exists but is not a text file", file=sys.stderr)
        exit()
    else:
        outfile = open(options.outfile, "w+")
        outfile.write(",".join(parser.all_cols) + "\n")
else:
    outfile = sys.stdout
    outfile.write(",".join(parser.all_cols) + "\n")

# Loop through each input file and parse
for f in input_files:
    infile = open(f, 'r')
    inlines = infile.readlines()
    if (len(inlines) < 4):
        print(f)
        continue

    infile_absolute_path = os.path.abspath(f)
    infile.close()

    # Pass each line to the parser
    for line in inlines:
        parser.parse_line(line, outfile, path=infile_absolute_path)

    # Record the file in the blacklist
    if parsed_list_file is not None:
        parsed_list_file.write(infile_absolute_path + "\n")

# Flush the aggregate in the parser
parser.print_aggregate(outfile)

if parsed_list_file is not None:
    parsed_list_file.close()

if outfile != sys.stdout:
    outfile.close()
