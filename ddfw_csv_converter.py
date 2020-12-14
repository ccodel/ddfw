# @file ddfw_to_csv.py
# @brief Converts DDFW output into a CSV file.
#
# @usage python ddfw_to_csv.py <ddfw_output> <csv_output> [parsed_file_list]
#
# <ddfw_output> is the output file from DDFW
# <csv_output> is the file to write the CSV output to. If there is no file,
#   a new CSV file is made with appropriate headers. If a file does exist,
#   then the data is appended to the end
# [parsed_file_list] is an optional argument that will cause the name of
#   the file <ddfw_output> to be appended to the end. At the start of
#   parsing, this file is checked to see if <ddfw_output> should be parsed 
#   again.
#
# @author Cayden Codel (ccodel@andrew.cmu.edu)
#
# @bug No known bugs.

import sys
import string
import os.path

# Verify number of passed-in arguments is correct
if len(sys.argv) != 3 and len(sys.argv) != 4:
    print('Incorrect number of arguments. Usage:')
    print(sys.argv[0] + ' <ddfw_output> <csv_output> [parsed_file_list]')
    sys.exit()

# First, check if optional blacklist file is provided
if len(sys.argv) == 4:
    if os.path.exists(sys.argv[3]):
        parsed_list_file = open(sys.argv[3], 'r')
        parsed_list_lines = parsed_list_file.readlines()

        for line in parsed_list_lines:
            if sys.argv[1] in line:
                sys.exit()

        parsed_list_file.close()

# Next, open the DDFW output file to read in solver data
ddfw_output_file = open(sys.argv[1], 'r')
ddfw_output_lines = ddfw_output_file.readlines()
ddfw_output_file.close()

# Check if the CSV output file exists
csv_file = None
if not os.path.exists(sys.argv[2]):
    csv_file = open(sys.argv[2], 'w+')
    csv_file.write('instance,ddfw_version,a,A,c,C,w,cpu_timeout,flip_timeout,seed,solve,flips,best,best_step,time\n')
else:
    csv_file = open(sys.argv[2], 'a')

# Allocate the output columns
instance_str = sys.argv[1]
version = -1
cpu_timeout = -1
flip_timeout = -1
seed = 0
a = 0.0
A = 0.0
c = 0.0
C = 0.0
w = 0.0
runs = 0
solve = []
flips = []
best = []
best_step = []
time = []

# Next, read in the data from the output file
for line in ddfw_output_lines:
    # Disregard empty lines or those that aren't comment lines
    if len(line) == 0 or line[0] != 'c':
        continue

    # Chop off the first 'c ' part of the comment string
    words = line[2:].strip().split(' ')

    if len(words) <= 1:
        continue

    # Extract configuration and run statistics
    if 'Version' == words[0]:
        version = int(words[1])
    elif '-w' == words[0]:
        w = float(words[1])
    elif '-a' == words[0]:
        a = float(words[1])
    elif '-A' == words[0]:
        A = float(words[1])
    elif '-c' == words[0]:
        c = float(words[1])
    elif '-C' == words[0]:
        C = float(words[1])
    elif '-r' == words[0]:
        runs = int(words[1])
    elif '-s' == words[0]:
        seed = int(words[1])
    elif '-t' == words[0]:
        if not 'infty' == words[1]:
            cpu_timeout = int(words[1])
    elif '-T' == words[0]:
        if not 'infty' == words[1]:
            flip_timeout = int(words[1])
    elif 'Stats:' == words[0]:
        stats = line.replace(' ', '').replace(':', '').strip().strip(string.ascii_letters).split('|')

        solve.append(int(stats[1]))
        flips.append(int(stats[2]))
        best.append(int(stats[3]))
        best_step.append(int(stats[4]))
        time.append(float(stats[5]))

# Check if the number of runs equals the claimed number
if (len(solve) != runs):
    print(sys.argv[1] + ' claims ' + str(runs) + ' runs, but ' + str(len(solve)) + ' were found')
    sys.exit()

print(len(solve))

# For each run, print its data into the CSV file
for i in range(0, len(solve)):
    print(i)
    csv_file.write('{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n'.format(
        instance_str, version, a, A, c, C, w, cpu_timeout, flip_timeout,
        seed, solve[i], flips[i], best[i], best_step[i], time[i]))
csv_file.close()

# Write the name of the instance to the list of parsed files
if len(sys.argv) == 4:
    parsed_list_file = open(sys.argv[3], 'a+')
    parsed_list_file.write(sys.argv[1] + '\n')
    parsed_list_file.close()
