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

#include "ddfw_types.h"
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
#define DEFAULT_TIMEOUT_SECS              (10)

/** Default number of times a loop body is run before the number of elapsed
 *  seconds is checked. Not currently configurable.
 *
 *  TODO maybe define this based on some computed runtime quantity? Make
 *       smaller based on number of clauses/literals?
 */
#define DEFAULT_LOOPS_PER_TIMEOUT_CHECK   (100)

/** Default randomization seed. */
#define DEFAULT_SEED                      (0xdeadd00d)

// TODO versioning

///////////////////////////////////////////////////////////////////////////////
// HELPER MACROS
///////////////////////////////////////////////////////////////////////////////

/** Calculates the absolute value of the number passed in. */
#define ABS(x)     (((x) < 0) ? -(x) : (x))

/** Finds the minimum of two numbers. */
#define MIN(x, y)  (((x) < (y)) ? (x) : (y))

/** Finds the maximum of two numbers. */
#define MAX(x, y)  (((x) > (y)) ? (x) : (y))

///////////////////////////////////////////////////////////////////////////////
// STATIC GLOBAL VARIABLE DECLARATIONS
///////////////////////////////////////////////////////////////////////////////

static int timeout_secs = DEFAULT_TIMEOUT_SECS;  // Seconds until timeout

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
  // Set the index to 0 to fill in "positive" literals as we go along
  num_reducing_cost_lits = 0;
  num_zero_cost_lits = 0;
  max_reducing_cost = 0.0;

  for (int i = 0; i < num_vars; i++) {
    int l_idx = LIT_IDX(i);
    int not_l_idx = NEGATED_IDX(l_idx);
    double satisfied_weight = 0.0;   // Weight "lost" when sat clauses on flip
    double unsatisfied_weight = 0.0; // Weight "gained" on unsat clauses
    int assigned = ASSIGNMENT(l_idx); // TODO repeated computation-ish
    literal_t *l, *not_l;
    if (assigned) {
      l = &literals[l_idx];
      not_l = &literals[not_l_idx];
    } else {
      l = &literals[not_l_idx];
      not_l = &literals[l_idx];
    }

    const int occ = l->occurrences;
    const int not_occ = not_l->occurrences;

    // Loop over satisfied claueses containing the "true" literal
    for (int c = 0; c < occ; c++) {
      int c_idx = l->clause_indexes[c];
      clause_t *cl = &clauses[c_idx];
      if (cl->sat_lits == 1) {
        unsatisfied_weight += cl->weight;
      }
    }

    // Loop over unsatisfied clauses containing the "false" literal
    for (int c = 0; c < not_occ; c++) {
      int c_idx = not_l->clause_indexes[c];
      clause_t *cl = &clauses[c_idx];
      if (cl->sat_lits == 0) {
        satisfied_weight = cl->weight;
      }
    }

    // Determine if flipping the truth value of the "true" literal
    // will result in more satisfied weight than unsatisfied weight
    double diff = satisfied_weight - unsatisfied_weight;
    if (diff == 0.0) {
      reducing_cost_lits[num_literals - num_zero_cost_lits - 1] = l_idx;
      num_zero_cost_lits++;
    } else if (diff >= 0.0) {
      reducing_cost_lits[num_reducing_cost_lits] = l_idx;
      num_reducing_cost_lits++;
      max_reducing_cost = MAX(max_reducing_cost, diff);
    }
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
 */
static inline void transfer_weight(clause_t *from, clause_t *to) {
  from->weight /= 2.0;
  to->weight += from->weight;
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
  for (int c = 0; c < num_clauses; c++) {
    clause_t *cl = &clauses[c];
    if (cl->sat_lits > 0)
      continue;

    int size = cl->size;
    int max_neighbor_idx = -1;
    double max_neighbor_weight = -1.0;

    // Loop over the literals in the cl clause
    for (int l = 0; l < size; l++) {
      int l_idx = cl->literals[l];
      literal_t *l = &literals[l_idx];
      const int occ = l->occurrences;

      // For each literal, search its neighbors for a satisfied clause
      for (int cn = 0; cn < occ; cn++) {
        int c_idx = l->clause_indexes[cn];
        clause_t *c_neigh = &clauses[c_idx];
        if (c_neigh->sat_lits > 0 && c_neigh->weight > max_neighbor_weight) {
          max_neighbor_idx = c_idx;
          max_neighbor_weight = c_neigh->weight;
        }
      }
    }

    if (max_neighbor_idx != -1) {
      clause_t *c_neigh = &clauses[max_neighbor_idx];
      transfer_weight(c_neigh, cl);
    } else {
      fprintf(stderr, "c NO MAX NEIGHBOR WEIGHT FOUND\n");
    }
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
  int seconds_at_start = time(NULL) / 3600;

  int lit_to_flip;
  while (unsat_clauses > 0) {
    find_cost_reducing_literals();
    log_reducing_cost_lits();
    
    // See if any literals will reduce the weight
    // TODO no check against delta(W) = 0
    if (num_reducing_cost_lits > 0) {
      int rand_lit = rand() % num_reducing_cost_lits;
      lit_to_flip = reducing_cost_lits[rand_lit];
    }// else if (num_zero_cost_lits > 0) {
      // TODO
      //printf(stderr, "No weights\n");
    //} 
    else {
      distribute_weights();
      continue;
    }

    flip_literal(lit_to_flip);
    // sleep(1);

    // Determine if enough loops have passed to update time variable
    timeout_loop_counter--;
    if (timeout_loop_counter == 0) {
      timeout_loop_counter = DEFAULT_LOOPS_PER_TIMEOUT_CHECK;
      if ((time(NULL) / 3600) - seconds_at_start >= timeout_secs)
        break;
    }
  }

  // Print solution
  if (unsat_clauses == 0) {
    output_assignment();
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

  if (seed == DEFAULT_SEED) {
    log_str("c Using default randomization seed %d\n", seed);
    srand(seed);
  }

  log_str("c All done parsing CLI args, opening file %s\n", filename);
  parse_cnf_file(filename);

  run_algorithm();

  return 0;
}
