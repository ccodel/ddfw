# DDFW - Divide and Distribute Fixed Weights
## Implemented by Cayden Codel, advised by Marijn Heule

This repository contains an implementation of the "Divide and Distribute Fixed Weights" algorithm presented in the paper ["Neighbourhood Clause Weight Redistribution in Local Search for SAT"](https://www.researchgate.net/publication/29453919_Neighbourhood_Clause_Weight_Redistribution_in_Local_Search_for_SAT "DDFW paper") by Ishtaiwi, Thornton, Sattar, and Pham. The implementation in this repository differs from the original paper in that the clause weights are floating point numbers instead of integers, which allows for greater expression in the kinds of weight distribution policies available to the algorithm.

## Compilation
To compile, change your directory to `src/` and type
```bash
make
```
and an executable `ddfw` will be generated at the root of the project. Typing
```bash
make clean
```
will delete the executable, `.o` files, and other files used in compilation. To compile a version of the algorithm that runs verification of flips, weight distribution, and other checks, type
```bash
make verify
```
to generate an executable `ddfw_verify`. The executable options for `ddfw_verify` are the same as for `ddfw`, but additional code is included.

## Execution options
For an up-to-date listing of the executable options, type `./ddfw -h` for a help message. Below are the compilation options and some documentation surrounding each.

```bash
./ddfw -f <cnf_filename> [-options]

  -a <double>
  -A <double>
  -c <double>
  -C <double>
  -d
  -f <filename>
  -h
  -l <flips>
  -m <method>
  -q
  -r <runs>
  -s <seed>
  -t <timeout>
  -T <flips>
  -v
  -w <double>
```

#### The -a, -A, -c, and -C options
The `-a`, `-A`, `-c`, and `-C` options each allow the constants used in distributing weights to be specified. This implementation of DDFW updates weights of neighboring clauses according to the following formula:
```C
// Transfer weight from Cf to Ct
if (Cf.weight < W_init) {
  Cf.weight = A * Cf.weight + C;
} else {
  Cf.weight = a * Cf.weight + c;
}

Ct.weight += drop_in_Cf_weight
```
In other words, if the current weight of the satisfied clause is greater than the initial weight given to each clause, then its weight is updated by multiplying by a constant and adding a constant. If the weight is below the initial weight, then a different pair of constants are used.

The original DDFW paper took a weight of 2 from `Cf` if the weight was greater than the original weight, and 1 otherwise. Thus, the settings of the original DDFW paper in this implementation are `-a 1.0 -A 1.0 -c -2.0 -C -1.0`. To use the original settings of the DDFW paper, use the `-d` flag.

If none of the four above options are specified, then their default values are 1.0 for `-a` and `-A`, and -2.0 for `-c` and `-C`.

#### The -d option
As mentioned above, the `-d` option uses the default settings in the original DDFW paper. When using `-d`, do not specify any of the `-aAcC` or `-w` options, as those are taken care of with the `-d` flag.

#### The -f option
Specify a path to a CNF file. The file does not need to end in `.cnf`, but does need to conform to the DIMACS CNF format. The `-f` option is required in order to run DDFW.

#### The -h option
Including the `-h` option displays a help message and ends execution, regardless of however many other options are specified.

#### The -l option
In-run statistics regarding weight transfer can be logged. Supplying a `-l`
option will cause those statistics to be logged every `<flips>` flips. If the
`-l` option is not provided, then no in-run statistics are logged.

#### The -m option
By default, method is WEIGHTED.

#### The -q option
By default, `ddfw` runs at a "NORMAL" level of verbosity. Specifying the `-q` option makes `ddfw` run in "QUIET" mode, and so only necessary output is printed, such as which settings are specified for the current run and the solution, if found.

#### The -Q option
The `-q` option, but with solution output suppressed.

#### The -r option
By default, `ddfw` runs only a single time. To run the algorithm multiple times on the same CNF file, use the `-r` option to specify the number of restarts. The algorithm will restart when it times out, depending on the method.

#### The -s option
Flipping variables and choosing neighboring clauses depends on a pseudo-random number generator. The generator is initialized with a default seed value. A different seed value can be specified with the `-s` option. The seed must fit within a 32-bit signed integer value, but can be negative.

#### The -t option
The number of seconds until `ddfw` times out can be specified with the `-t` option. By default, `ddfw` runs for 100 seconds.

#### The -T option
The number of flips until `ddfw`  times out. The algorithm is assumed to run until a certain number of seconds has passed, but specifying the `-T` option will instead make the algorithm time out after a certain number of flips. Do not specify both the `-t` and the `-T` options together, as behavior is then undefined.

#### The -v option
By default, `ddfw` runs at a "NORMAL" level of verbosity. Specifying the `-v` option makes `ddfw` run in "VERBOSE" mode, and so much more debugging output is printed.

#### The -w option
After reading in the CNF file, all clauses are given the same initial weight. Specify what initial weight is given with the `-w` option. The original DDFW paper used `8.0`. By default, `ddfw` gives a weight of 100.0 to each clause.
