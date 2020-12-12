/** @file main.c
 *  @brief Contains the main() function to run DDFW. Parses command-line args.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

#include "cnf_parser.h"
#include "clause.h"
#include "ddfw.h"
#include "logger.h"

/** @brief Main function. Processes command-line arguments and kicks off DDFW.
 *
 *  @param argc  The number of arguments given on the command line.
 *  @param argv  The array of string arguments given on the command line.
 */
int main(int argc, char *argv[]) {
  if (argc == 1) {
    print_usage(argv[0]);
    return 0;
  }

  int ao, Ao, co, Co, d, q, Q, v, w;
  ao = Ao = co = Co = d = q = Q = v = w = 0;
  int seed = DEFAULT_RAND_SEED;
  char *filename = NULL;
  extern char *optarg;
  char opt;
  while ((opt = getopt(argc, argv, "dhvqQa:A:c:C:f:m:r:s:t:T:w:")) != -1) {
    switch (opt) {
      case 'a': ao = 1; mult_a = atof(optarg); break;
      case 'A': Ao = 1; mult_A = atof(optarg); break;
      case 'c': co = 1; add_c = atof(optarg);  break;
      case 'C': Co = 1; add_C = atof(optarg);  break;
      case 'd': d = 1;                         break;
      case 'f': filename = optarg;             break;
      case 'h': print_help(argv[0]);           return 0;
      case 'm':
        switch (optarg[0]) {
          case 'U':
            selection_method = UNIFORM;
            break;
          case 'W':
            selection_method = WEIGHTED;
            break;
          case 'B':
            selection_method = BEST;
            break;
          default:
            fprintf(stderr, "Unrecognized variable selection method\n");
            return 0;
        }
        break;
      case 'q': q = 1; break;
      case 'Q': Q = 1; break;
      case 'r':
        num_restarts = atoi(optarg);
        if (num_restarts < 1) {
          fprintf(stderr, "Algorithm must run at least once\n");
          return 0;
        }
        break;
      case 's':
        seed = atoi(optarg); break;
      case 't':
        timeout_secs = atoi(optarg);
        if (timeout_secs < 1) {
          fprintf(stderr, "Algorithm must run for at least one second\n");
          return 0;
        }

        if (timeout_method == DEFAULT) {
          timeout_method = TIME;
        } else {
          timeout_method = BOTH;
        }
        break;
      case 'T':
        timeout_flips = atoi(optarg);
        if (timeout_flips < 1) {
          fprintf(stderr, "Algorithm must run for at least one flip\n");
          return 0;
        }

        if (timeout_method == DEFAULT) {
          timeout_method = FLIPS;
        } else {
          timeout_method = BOTH;
        }
        break;
      case 'v': v = 1; break;
      case 'w':
        w = 1;
        init_clause_weight = atof(optarg);
        if (init_clause_weight < 0.0) {
          fprintf(stderr, "Initial weight must be a positive number\n");
          return 0;
        }
        break;
      default: print_usage(argv[0]); return 0;
    }
  }

  // Check that no configuration settings are conflicting
  if (d == 1) {
    if (ao + Ao + co + Co != 0) {
      fprintf(stderr, "Cannot supply -d and -aAcC flags at the same time\n");
      return 0;
    }

    // Original DDFW configuration
    init_clause_weight = 8.0;
    mult_a = 1.0;
    mult_A = 1.0;
    add_c = -2.0;
    add_C = -1.0;
  }

  if (q + Q + v > 1) {
    fprintf(stderr, "More than one verbosity argument provided\n");
    return 0;
  }

  if (filename == NULL) {
    fprintf(stderr, "No input CNF filename provided\n");
    return 0;
  }

  // All command-line arguments have been validated, proceed with setup
  // Print banner
  log_str("c ------------------------------------------------------------\n");
  log_str("c          DDFW - Divide and Distribute Fixed Weights\n");
  log_str("c                  Implemented by Cayden Codel\n");
  log_str("c ------------------------------------------------------------\n");

  // Print version number and current configuration
  log_str("c  Version %d\nc\n", VERSION);
  log_str("c  Run configuration:\n");
  log_str("c  -f %s\n", filename);
  log_str("c  -w %lf\n", init_clause_weight);
  log_str("c  -a %lf\n", mult_a);
  log_str("c  -A %lf\n", mult_A);
  log_str("c  -c %lf\n", add_c);
  log_str("c  -C %lf\n", add_C);
  log_str("c  -r %d\n", num_restarts);
  switch (timeout_method) {
    case DEFAULT:
    case TIME:
      log_str("c  -t %d\n", timeout_secs);
      log_str("c  -T infty\n");
      break;
    case FLIPS:
      log_str("c  -t infty\n");
      log_str("c  -T %d\n", timeout_flips);
      break;
    case BOTH:
      log_str("c  -t %d\n", timeout_secs);
      log_str("c  -T %d\n", timeout_flips);
      break;
    default:
      fprintf(stderr, "Internal error\n");
      return 0;
  }
  log_str("c  -s %d\n", seed);

  // Set verbosity
  if (q == 1) {
    log_str("c  -q\nc\n");
    set_verbosity(SILENT);
  } else if (Q == 1) {
    log_str("c  -Q\nc\n");
    set_verbosity(SILENT);
    suppress_solution = 1;
  } else if (v == 1) {
    log_str("c  -v\nc\n");
    set_verbosity(VERBOSE);
  } else {
    log_str("c\n");
    set_verbosity(NORMAL);
  }
  printf("c\n");

  log_str("c All done parsing CLI args, opening file %s\n", filename);
  parse_cnf_file(filename);

  // Run the algorithm as many times as specified at the command-line
  srand(seed);
  for (algorithm_run = 1; algorithm_run <= num_restarts; algorithm_run++) {
    run_ddfw_algorithm();
    reset_data_structures();
  }

  return 0;
}
