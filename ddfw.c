/** @file ddfw.c - Divide and Distribute Fixed Weights
 *  @brief An implementation of the DDFW algorithm, as presented in 
 *         the paper: "Neighbourhood Clause Weight Redistribution in 
 *                     Local Search for SAT" [Ishtaiwi, Thornton, Sattar, Pham]
 *
 *  Note that this implementation of DDFW assumes that:
 *   - No clause has more than MAX_CLAUSE_LENGTH literals in it
 *   - No clause contains a tautology (e.g. (x v !x))
 *  
 *  ///////////////////////////////////////////////////////////////////////////
 *  // HIGH LEVEL OVERVIEW OF DDFW
 *  ///////////////////////////////////////////////////////////////////////////
 *  
 *  TODO High level discussion here.
 *
 *  ///////////////////////////////////////////////////////////////////////////
 *  // IMPLEMENTATION
 *  ///////////////////////////////////////////////////////////////////////////
 *  
 *  TODO Implementation discussion here.
 *
 *  ///////////////////////////////////////////////////////////////////////////
 *  // QUIRKY IMPLEMENTATION DETAILS
 *  ///////////////////////////////////////////////////////////////////////////
 *
 *  According to DIMACS format, variables are specified by nonzero integers.
 *  A positive integer corresponds to a positive literal l, whereas a negative
 *  integer corrsponds to a negated literal, !l. Note that the problem line
 *  defines the total number of *variables*, not the total number of literals.
 *  Therefore, there are twice as many literals as there are variables.
 *
 *  In this implementation, each literal has a unique literal_t struct that
 *  is allocated contiguously in the literals array. Literals are stored
 *  adjacently to their negated form. So, if the literal is (l = 5), then
 *  the index into the literals array to access the literal_t struct for
 *  l is (2 * 5), and the index for !l is (2 * 5) + 1. Because DIMACS variables
 *  are 1-indexed, the literals array has "two too many" to preserve the
 *  1-indexing of DIMACS format.
 *
 *  The clauses are 0-indexed, however, as the clause_t array, clauses, is
 *  built internally.
 *
 *  The motivation behind using a [l1, !l1, l2, !l2, ...] literals storage
 *  structure, as opposed to all positive literals on one side of a pointer
 *  and all negated literals on the other, is to take advantage of caching
 *  when accessing l and !l together. For example, when flipping a literal
 *  in the assignment, clauses for both l and !l must be updated, and so
 *  both literal_t structs are accessed together. Whether or not such caching
 *  effects can be felt during runtime is something the author has not
 *  investigated, and stands as a potential optimization.
 *
 *  ///////////////////////////////////////////////////////////////////////////
 *  // CREDITS AND ACKNOWLEDGEMENTS
 *  ///////////////////////////////////////////////////////////////////////////
 *
 *  @author   Cayden Codel (ccodel@andrew.cmu.edu)
 *  @advisor  Marijn Heule (mheule@andrew.cmu.edu)
 *
 *  Original DDFW paper authored by
 *  @author   Abdelraouf Ishtaiwi (a.ishtaiwi@giffith.edu.au)
 *  @author   John Thornton       (j.thornton@griffith.edu.au)
 *  @author   Abdul Sattar        (a.sattar@griffith.edu.au)
 *  @author Duc Nghia Pham        (d.n.pham@griffith.edu.au)
 *  
 *  ///////////////////////////////////////////////////////////////////////////
 *  // BUGS AND TODOS
 *  ///////////////////////////////////////////////////////////////////////////
 *
 *  @bug No known bugs.
 *  @todos Many todos.
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <getopt.h>
#include <time.h>

#include "clause.h"
#include "cnf_parser.h"
#include "logger.h"

///////////////////////////////////////////////////////////////////////////////
// CONFIGURATION MACROS
///////////////////////////////////////////////////////////////////////////////

/** @brief The maximum number of command line arguments allowed.
 *
 *  While this number alone will not catch ill-formed command line arguments,
 *  if does provide a quick way to short-circuit in main().
 *
 *  UPDATE THIS NUMBER IF NEW COMMAND LINE ARGUMENTS ARE IMPLEMENTED.
 */
#define MAX_CLI_ARGS  10

/** Default number of seconds to loop before timing out. 
 *  Can be toggled with the -t <timeout_secs> flag at the command line.
 **/
#define DEFAULT_TIMEOUT_SECS              (100)

/** Default number of times a loop body is run before the number of elapsed
 *  seconds is checked. Not currently configurable.
 *
 *  TODO maybe define this based on some computed runtime quantity? Make
 *       smaller based on number of clauses/literals?
 */
#define DEFAULT_LOOPS_PER_TIMEOUT_CHECK   (100000)

/** Default randomization seed. */
#define DEFAULT_SEED                      (0xdeadd00d)

/** The transfer weight equation is
 *
 *  clause_to weight += a (clause_from weight) + c
 *
 *  where a is a multiplicative scalar and c a constant.
 *
 *  As the initial weight in the DDFW paper was 8, then we set
 *
 *  a = 0.75, c = 0.25
 */
#define A   (0.75)
#define C   (0.25)

/** @brief The default probability that the maximum weight neighboring clause
 *         is disregarded for a random clause.
 *
 *  TODO configurable?
 */
#define DEFAULT_RANDOM_CLAUSE_PROB          (0.01)

// TODO versioning

///////////////////////////////////////////////////////////////////////////////
// HELPER MACROS
///////////////////////////////////////////////////////////////////////////////

/** Calculates the absolute value of the number passed in. */
#ifndef ABS
#define ABS(x)     (((x) < 0) ? -(x) : (x))
#endif

/** Finds the minimum of two numbers. */
#ifndef MIN
#define MIN(x, y)  (((x) < (y)) ? (x) : (y))
#endif

/** Finds the maximum of two numbers. */
#ifndef MAX
#define MAX(x, y)  (((x) > (y)) ? (x) : (y))
#endif

///////////////////////////////////////////////////////////////////////////////
// STATIC GLOBAL VARIABLE DECLARATIONS
///////////////////////////////////////////////////////////////////////////////

static int timeout_secs = DEFAULT_TIMEOUT_SECS;  // Seconds until timeout

static int random_clause_prob = (int) (1.0 / DEFAULT_RANDOM_CLAUSE_PROB);

////////////////////////////////////////////////////////////////////////////////
// DDFW algorithm implementation
////////////////////////////////////////////////////////////////////////////////

/** @brief Finds those literals that cause a positive change in weighted
 *         cost if flipped and stores them in reducing_cost_lits.
 *         Stores the number of such literals in num_reducing_cost_lits.
 *
 *  According to the DDFW algorithm, the first thing done per each loop body
 *  is to "find and return a list L of literals causing the greatest reduction
 *  in weighted cost delta(w) when flipped."
 *
 *  TODO CURRENTLY INEFFICIENT
 *  To do so, the literals are looped over. The change in weighted cost for
 *  literal l is defined as the cost lost minus the cost gained from making
 *  clauses unsatisfied.
 *
 *  Let W(l) be the total weight of clauses with l or !l, sat or unsat.
 *  Then if delta(w) for l is W, delta(w) for !l is "break." TODO
 *
 */
static void find_cost_reducing_literals() {
  // Loop through those variables in the cost_compute_vars to compute cost
  // Most efficient to do from the back of the cost compute vars array
  int *cc_vars = cost_compute_vars + num_cost_compute_vars - 1;
  const int num_cc_vars = num_cost_compute_vars;
  for (int i = 0; i < num_cc_vars; i++) {
    const int v_idx = *cc_vars;
    const int l_idx = LIT_IDX(v_idx);
    const int assigned = ASSIGNMENT(l_idx);
    double satisfied_weight = 0.0;
    double unsatisfied_weight = 0.0;

    // log_str("Checking var %d, truth value %d;", v_idx, assigned);

    // Determine the index of the true literal
    int true_idx, false_idx;
    if (assigned) {
      true_idx = l_idx;
      false_idx = NEGATED_IDX(l_idx);
    } else {
      true_idx = NEGATED_IDX(l_idx);
      false_idx = l_idx;
    }

    const int true_occ = literal_occ[true_idx];
    const int false_occ = literal_occ[false_idx];

    // Loop over satisfied clauses containing the true literal
    int *l_to_clauses = literal_clauses[true_idx];
    for (int c = 0; c < true_occ; c++) {
      const int c_idx = *l_to_clauses;
      if (clause_num_true_lits[c_idx] == 1) {
        unsatisfied_weight += clause_weights[c_idx];
      }

      l_to_clauses++;
    }

    // Loop over unsatisfied clauses containing the false literal
    l_to_clauses = literal_clauses[false_idx];
    for (int c = 0; c < false_occ; c++) {
      const int c_idx = *l_to_clauses;
      if (clause_num_true_lits[c_idx] == 0) {
        satisfied_weight += clause_weights[c_idx];
      }

      l_to_clauses++;
    }

    // Determine if flipping the truth value of the true literal
    //   would result in more satisfied weight than unsatisfied weight
    const double diff = satisfied_weight - unsatisfied_weight;
    // log_str(" diff %.4f\n", diff);
    if (diff > 0.0) {
      add_cost_reducing_var(v_idx); // Filter out 0 cost variables?
    } else if (diff <= 0.0) {
      remove_cost_reducing_var(v_idx);
    }

    // Remove var from compute list
    remove_cost_compute_var(v_idx);
    cc_vars--;
  }
}

/** @brief Transfers weight from one clause to another in the distribute step.
 *
 *  Since the transferral was specified one way for the discrete case, I
 *  will be playing around with the continuous case here. For now, the
 *  weight is distributed by
 *   
 *   - halving the weight of the satisfied clause
 *   - adding that weight to the unsatisfied clause
 *
 *  @param from_idx The clause index that weight is taken from.
 *  @param to_idx   The clause index that weight is given to.
 */
static inline void transfer_weight(int from_idx, int to_idx) {
  double orig = clause_weights[from_idx];
  clause_weights[from_idx] *= A;
  clause_weights[from_idx] += C;
  double diff = orig - clause_weights[from_idx];
  clause_weights[to_idx] += diff;

  const int from_size = clause_sizes[from_idx];
  const int to_size = clause_sizes[to_idx];

  // In the satisfied "from" clause, add vars to cost compute
  //   only if the clause is "critical"
  if (clause_num_true_lits[from_idx] <= 1) {
    int *from_lits = clause_literals[from_idx];
    for (int l = 0; l < from_size; l++) {
      const int l_idx = *from_lits;
      add_cost_compute_var(VAR_IDX(l_idx));
      from_lits++;
    }
  }

  // In the unsatisfied "to" clause, add vars to cost compute
  int *to_lits = clause_literals[to_idx];
  for (int l = 0; l < to_size; l++) {
    const int l_idx = *to_lits;
    add_cost_compute_var(VAR_IDX(l_idx));
    to_lits++;
  }
}

/** @brief Distribute weights from satisfied to unsatisfied clauses.
 *  
 *  In the case where no literals may be flipped to decrease the weight
 *  of the unsat clauses, the weight will be re-distributed from satisfied
 *  to unsatisfied clauses according to the following rule:
 *
 *  for each unsatisfied clause c,
 *    select a satisfied same sign neighboring clause cn with max weight wn
 *
 */
static void distribute_weights() {
  log_str("c Distributing weights\n");

  // Loop over all clauses, picking out those that are false
  const int uc = num_unsat_clauses;
  int *false_clauses = false_clause_members;
  for (int c = 0; c < uc; c++) {
    const int c_idx = *false_clauses;
    const int size = clause_sizes[c_idx];

    int max_neighbor_idx = -1;
    double max_neighbor_weight = -1.0;

    // Loop over the literals in the clause to search for neighboring sat clause
    int *cl_lits = clause_literals[c_idx];
    for (int l = 0; l < size; l++) {
      const int l_idx = *cl_lits;
      const int occ = literal_occ[l_idx];

      // For each literal, search neighbors for a satisfied clause
      int *l_clauses = literal_clauses[l_idx];
      for (int cn = 0; cn < occ; cn++) {
        const int cn_idx = *l_clauses;
        if (clause_num_true_lits[cn_idx] > 0 
            && clause_weights[cn_idx] > max_neighbor_weight) {
          max_neighbor_idx = cn_idx;
          max_neighbor_weight = clause_weights[cn_idx];
        }

        l_clauses++;
      }

      if (max_neighbor_idx != -1) {
        if (rand() % random_clause_prob != 0) {
          transfer_weight(max_neighbor_idx, c_idx);
        } else {
          // Loop over clauses randomly until satisfied, good weight clause
          unsigned int r_idx;
          do {
            r_idx = ((unsigned int) rand()) % occ;
          } while (clause_num_true_lits[r_idx] == 0); // TODO weight case

          transfer_weight(r_idx, c_idx);
        }
      }

      cl_lits++;
    }

    false_clauses++;
  }
}

/** @brief Runs the DDFW algorithm.
 *
 *  Main loop.
 */
static void run_algorithm() {
  generate_random_assignment();

  // Record the time to ensure no timeout
  int timeout_loop_counter = DEFAULT_LOOPS_PER_TIMEOUT_CHECK;
  int seconds_at_start = time(NULL);

  int lit_to_flip;
  while (num_unsat_clauses > 0) {
    // log_weights();
    // TODO extreme reset
    /*
    for (int i = 1; i <= num_vars; i++) {
      add_cost_compute_var(i);
    }
    */

    find_cost_reducing_literals();
    log_reducing_cost_lits();
    
    // See if any literals will reduce the weight
    // TODO no check against delta(W) = 0
    if (num_cost_reducing_vars > 0) {
      unsigned int rand_lit = ((unsigned int) rand()) % num_cost_reducing_vars;
      lit_to_flip = cost_reducing_vars[rand_lit];
    }// else if (num_zero_cost_lits > 0) {
      // TODO
      //printf(stderr, "No weights\n");
    //} 
    else {
      distribute_weights();
      continue;
    }

    flip_variable(lit_to_flip);

    // Determine if enough loops have passed to update time variable
    timeout_loop_counter--;
    if (timeout_loop_counter == 0) {
      printf("c There are now %d unsat clauses after %d flips\n",
          num_unsat_clauses, num_flips);
      timeout_loop_counter = DEFAULT_LOOPS_PER_TIMEOUT_CHECK;

      // Check for timeout
      if (time(NULL) - seconds_at_start >= timeout_secs)
        break;
    }
  }

  // Print solution
  if (num_unsat_clauses == 0) {
    output_assignment();
    printf("\n\nc Satisfied in %ld seconds\n", time(NULL) - seconds_at_start);
  } else {
    printf("s UNSATISFIABLE\n");
    printf("c Could not solve due to timeout.\n");
  }
}

////////////////////////////////////////////////////////////////////////////////
// Main execution
////////////////////////////////////////////////////////////////////////////////

/** @brief Main function. Processes input arguments and kicks off solver.
 *
 *  @param argc The number of arguments to the command line.
 *  @param argv The array of string arguments given on the command line.
 */
int main(int argc, char *argv[]) {
  if (argc == 1 || argc > MAX_CLI_ARGS) {
    print_usage(argv[0]);
    exit(1);
  }

  int seed = DEFAULT_SEED;
  char *filename = NULL;
  extern char *optarg;
  char opt;
  while ((opt = getopt(argc, argv, "hvqf:s:t:")) != -1) {
    switch (opt) { 
      case 'f':
        filename = optarg;
        break;
      case 'h':
        print_help(argv[0]);
        return 0;
      case 'q':
        set_verbosity(SILENT);
        break;
      case 's':
        seed = atoi(optarg);
        if (seed != 0) {
          log_str("c Using randomization seed %d", seed);
          srand(seed);
        }
        break;
      case 't':
        timeout_secs = atoi(optarg);
        break;
      case 'v':
        set_verbosity(VERBOSE);
        break;
      default:
        print_usage(argv[0]);
    }
  }

  printf("c ------------------------------------------------------------\n");
  printf("c          DDFW - Divide and Distribute Fixed Weights\n");
  printf("c                  Implemented by Cayden Codel\n");
  printf("c                          Version 0.1\n");
  printf("c ------------------------------------------------------------\n");

  if (filename == NULL) {
    fprintf(stderr, "c No filename provided, exiting.\n");
    exit(1);
  }

  log_str("c Running with timeout of %d seconds\n", timeout_secs);

  if (seed == DEFAULT_SEED) {
    log_str("c Using default randomization seed %d\n", seed);
    srand(seed);
  }

  log_str("c All done parsing CLI args, opening file %s\n", filename);
  parse_cnf_file(filename);

  /*
  // TODO remove when verify copy is good TODO calloc to array compare
  reducing_cost_copy = calloc(num_literals, sizeof(int));
  if (reducing_cost_copy == NULL) {
    fprintf(stderr, "c Ran out of memory on copy, exiting\n");
    exit(-1);
  }
  */

  run_algorithm();

  return 0;
}
