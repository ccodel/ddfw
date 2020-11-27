/** @file logger.c
 *  @brief Implementations of logging functions.
 *
 *  This file provides the implementations of functions for logging various
 *  aspects of the DDFW algorithm. Some functions are called every loop, while
 *  others are only called when the algorithm is run with the -v flag.
 *
 *  The functions implemented below are of two types:
 *
 *   - Printing functions: always print to stdout, regardless of verbosity.
 *   - Logging functions: print to stdout, according to verbosity.
 *
 *  Printing functions are prefixed by "print_", while logging functions
 *  are prefixed by "log_".
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>

#include "logger.h"
#include "clause.h"

#ifndef ABS
#define ABS(x)   (((x) < 0) ? -(x) : (x))
#endif

#ifndef MIN
#define MIN(x, y)  (((x) < (y)) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x, y)  (((x) > (y)) ? (x) : (y))
#endif

/** The verbosity level maintained by the logger. Set by set_verbosity(). 
 *  The default level is NORMAL.
 */
static verbose_t vlevel = NORMAL;

/** @brief Emits an error message and exits. */
static void unrecognized_verbosity_level() {
  printf("c ERROR: UNRECOGNIZED VERBOSITY LEVEL, EXITING...\n");
  exit(-1);
}

/** @brief Print help information when the -h flag is used.
 *
 *  When the -h flag is specified at the command line, a help message is
 *  printed by this function. Included in the help message are the various
 *  flags the algorithm allows.
 *
 *  @param runtime_path The string name of the running executable, e.g. ./ddfw
 */
void print_help(char *runtime_path) {
  printf("\n%s: Divide and Distribute Fixed Weights\n\n", runtime_path);
  printf("  -a <double>         Provide");
  printf(" a multiplicative constant for weight distribution.\n");
  printf("  -c <double>         Provide");
  printf(" an additive constant for weight distribution.\n");
  printf("  -A <double>         Provide");
  printf(" a multiplicative constant for WD under W_init.\n");
  printf("  -C <double>         Provide");
  printf(" an additive constant for WD under W_init.\n");
  printf("  -d                  Use DDFW original paper settings.\n");
  printf("  -f <filename>       Provide a .cnf file.\n");
  printf("  -h                  Display this help message.\n");
  printf("  -m <method>         Selection method");
  printf(" ([U]niform, [W]eighted, [B]est)\n");
  printf("  -q                  Quiet mode, no printing.\n");
  printf("  -r <runs>           Provide an optional number of runs\n");
  printf("  -s <seed>           Provide an optional randomization seed.\n");
  printf("  -t <timeout>        Provide");
  printf(" an optional number of seconds until timeout.\n");
  printf("  -T <flips>          Provide");
  printf(" an optional number of flips until timeout.\n");
  printf("  -v                  Turn on verbose printing.\n");
  printf("  -w <double>         Initial weight for all clauses.\n");
  printf("\n");

}

/** @brief Print usage information when the provided arguments aren't
 *         well-formed.
 *
 *  This implementation of the DDFW allows for several parameters to be
 *  specified or tuned by the user at the command line. If these arguments
 *  aren't well-formed, usage information is printed instead of execution
 *  continuing like normal. This function only prints the usage information.
 *
 *  @param runtime_path The string name of the running executable, e.g. ./ddfw
 */
void print_usage(char *runtime_path) {
  printf("Usage: %s [-hqv] -f <filename> [-s <seed>] [-t <timeout>]\n", 
      runtime_path);
  printf("                [-w <double>] [-aAcC <double>]\n");
}

/** @brief Returns the logger's verbosity level.
 *
 *  @return The logger's current verbosity level.
 */
verbose_t get_verbosity() {
  return vlevel;
}

/** @brief Sets the logger's verbosity level.
 *
 *  When the verbosity level is set, regardless of the level it is set to,
 *  a comment line is emitted by the logger noting the new level of the logger.
 *
 *  @param verbosity_level The verbosity level the logger should assume.
 */
void set_verbosity(verbose_t verbosity_level) {
  vlevel = verbosity_level;
  switch (verbosity_level) {
    case SILENT:
      printf("c Verbosity level set to silent.\n");
      break;
    case NORMAL:
      printf("c Verbosity level set to normal.\n");
      break;
    case VERBOSE:
      printf("c Verbosity level set to verbose.\n");
      break;
    default:
      unrecognized_verbosity_level();
  }
}

/** @brief Logs clause information according to verbosity level.
 *
 *  @param c The index to the clause to log.
 */
void log_clause(int c_idx) {
  switch (vlevel) {
    case SILENT:
      return;
    case NORMAL:
    case VERBOSE:
      printf("c Clause %d, weight: %.2f, sz: %d, sat lits: %d\n",
          c_idx, clause_weights[c_idx], clause_sizes[c_idx], 
          clause_num_true_lits[c_idx]);
      break;
    default:
      unrecognized_verbosity_level();
  }
}

/** @brief Logs literal information according to verbosity level.
 *
 *  @param l A pointer to the literal to log.
 */
void log_literal(int l_idx) {
  switch (vlevel) {
    case SILENT:
      return;
    case NORMAL:
      printf("c Literal %d (var %d), occurrences: %d\n",
          l_idx, VAR_IDX(l_idx), literal_occ[l_idx]);
      break;
    case VERBOSE:
      printf("c Literal %d (var %d), occurrences: %dc    In clauses  \n",
          l_idx, VAR_IDX(l_idx), literal_occ[l_idx]);
      const int occ = literal_occ[l_idx];
      int *l_to_clauses = literal_clauses[l_idx];
      for (int i = 0; i < occ; i++) {
        printf("%d ", l_to_clauses[i]);
      }
      printf("\n");
      break;
    default:
      unrecognized_verbosity_level();
  }
}

/** @brief Logs the provided string according to the specified format string.
 *
 *  If the verbosity level is SILENT, the function has no effect.
 *
 *  @param format A format string, such as those passed to printf().
 *  @param ... A variable list of tokens for the format string.
 */
void log_str(const char *format, ...) {
  switch (vlevel) {
    case SILENT:
      return;
    case NORMAL:
    case VERBOSE:;
      va_list ap;
      va_start(ap, format);
      vprintf(format, ap);
      break;
    default:
      unrecognized_verbosity_level();
  }
}

/** @brief Logs the provided string according to the specified format string.
 *         Prints to standard error.
 *
 *  Verbosity has no effect on error logging.
 *
 *  @param format A format string, such as those passed to printf().
 *  @param ... A variable list of tokens for the format string.
 */
void log_err(const char *format, ...) {
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr, format, ap);
}

/** @brief Logs the weights associated with each clause.
 *
 */
void log_weights() {
  if (vlevel == SILENT)
    return;

  for (int i = 0; i < num_clauses; i++) {
    if (clause_num_true_lits[i] <= 1) {
      log_clause(i);
    }
  }
}

/** @brief Logs the reducing cost literals.
 *
 *  If verbosity level is set to NORMAL, just those literals in the
 *  reducing cost lits array are logged.
 *
 *  If verbosity level is set to VERBOSE, then more information is
 *  output, such as how the calculation was made.
 */
void log_reducing_cost_lits() {
  if (vlevel == SILENT)
    return;

  // printf("c Found %d cost reducing literals\n", num_cost_reducing_vars);

  /*
  for (int i = 0; i < num_cost_reducing_vars; i++) {
    int l_idx = LIT_IDX(cost_reducing_vars[i]);
    int not_l_idx = NEGATED_IDX(l_idx);
    int assigned = ASSIGNMENT(l_idx);

    printf("c %d (var %d) is cost reducing, current truth value %d\n",
        l_idx, VAR_IDX(l_idx), assigned);
    if (vlevel == VERBOSE) {
      printf("c   Critical clauses for this literal:\n");

      const int occ = literal_occ[l_idx];
      int *l_to_clauses = literal_clauses[l_idx];
      for (int c = 0; c < occ; c++) {
        const int c_idx = *l_to_clauses;
        if (!assigned && clause_num_true_lits[c_idx] == 0) {
          printf("c     %d has 0 sat lit, weight: (+) %.4f\n", 
              c_idx, clause_weights[c_idx]);
        } else if (assigned && clause_num_true_lits[c_idx] == 1) {
          printf("c     %d has 1 sat lits, weight: (-) %.4f\n", 
              c_idx, clause_weights[c_idx]);
        }

        l_to_clauses++;
      }

      const int not_occ = literal_occ[not_l_idx];
      l_to_clauses = literal_clauses[not_l_idx];
      for (int c = 0; c < not_occ; c++) {
        const int c_idx = *l_to_clauses;
        if (assigned && clause_num_true_lits[c_idx] == 0) {
          printf("c      %d has 0 sat lits, weight: (+) %.4f\n", 
              c_idx, clause_weights[c_idx]);
        } else if (!assigned && clause_num_true_lits[c_idx] == 1) {
          printf("c     %d has 1 sat lit, weight: (-) %.4f\n",
              c_idx, clause_weights[c_idx]);
        }

        l_to_clauses++;
      }
    }
  }
  */
}

/** @brief Logs common statistics collected throughout the algorithm.
 *  
 *  Prints number of flips, time taken, best found so far.
 *
 *  Prints 
 */
void log_statistics(int run, struct timeval *start, struct timeval *stop) {
  // run | Solved? | flips | best | best_step | time
  printf("c Stats: %d | %d | %ld | %d | %d | %.3f\n",
      run, (num_unsat_clauses == 0), num_flips, lowest_unsat_clauses,
      lowest_unsat_step,
      ((double) (stop->tv_usec - start->tv_usec) / 1000000) + 
      ((double) (stop->tv_sec - start->tv_sec)));
}

/** @brief Logs the current assignment in the global formula variable.
 *
 *  The logged output prints the assignment in a grid 50 variables wide.
 *  Guiding numbers appear at the top and sides to help grid readability.
 */
void log_assignment() {
  switch (vlevel) {
    case SILENT:
      return;
    case VERBOSE:
      // Print out the top banner
      printf("c Assign:   1   5    10   5    20   5    30   5    40\n");
      printf("c           --------------------------------------------------");

      char *a = assignment + 1;
      for (int i = 0; i < num_vars; i++) {
        if (i % 50 == 0) {
          printf("\nc %06d    ", i);
        }

        char bit = *a;
        printf("%d", bit);
        a++;
      }
      printf("\n");
    case NORMAL:
      printf("c There are %d unsatisfied clauses after %ld flips\n", 
          num_unsat_clauses, num_flips);
      break;
    default:
      unrecognized_verbosity_level();
  }
}

/** @brief Prints the assignment in the satisfiable case.
 *
 *  Outputs assignment with 15 variables per line.
 *  TODO better formatting later
 */
void output_assignment() {
  if (num_unsat_clauses != 0) {
    printf("s UNSATISFIABLE\n");
  } else {
    printf("s SATISFIABLE\nv ");
    char *a = assignment + 1;
    for (int i = 1; i <= num_vars; i++) {
      char bit = *a;
      a++;

      if (bit == 0) {
        printf("-%d ", i);
      } else {
        printf("%d ", i);
      }

      // Newline character for formatting
      if (i % 10 == 0) {
        printf("\nv ");
      }
    }

    printf("0\n");
  }
}
