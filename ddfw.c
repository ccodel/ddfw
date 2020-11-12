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
 *  Potential improvement according to
 *
 *   https://www.keithschwarz.com/darts-dice-coins/
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
 *  @author   Duc Nghia Pham      (d.n.pham@griffith.edu.au)
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
#include <float.h>

#include "clause.h"
#include "cnf_parser.h"
#include "logger.h"
#include "xmalloc.h"

///////////////////////////////////////////////////////////////////////////////
// CONFIGURATION MACROS
///////////////////////////////////////////////////////////////////////////////

/** Default number of seconds to loop before timing out. 
 *  Can be toggled with the -t <timeout_secs> flag at the command line.
 **/
#define DEFAULT_TIMEOUT_SECS              (100)

/** Default number of flips before timing out. */
#define DEFAULT_TIMEOUT_FLIPS             1000000

/** Default number of times a loop body is run before the number of elapsed
 *  seconds is checked. Not currently configurable.
 *
 *  TODO maybe define this based on some computed runtime quantity? Make
 *       smaller based on number of clauses/literals?
 */
#define DEFAULT_LOOPS_PER_TIMEOUT_CHECK   (1000000)

/** Default randomization seed. */
#define DEFAULT_SEED                      (0xdeadd00d)

/** When transferring weights, weight is transferred from a clause cf to a
 *  clause ct. If we let the weight of cf to be wf and the weight of ct to
 *  be wt, then the weights of the two are adjusted according to:
 *
 *    cf = a * cf + c
 *    ct = cf - (a * cf + c)
 *
 *  To more closely reflect the original paper, a is set to 1 and c is -2.
 */
#define A   (1.0)
#define C   (-2.0)

/** @brief Numerator to random walk if no cost reducing literals */
#define DEFAULT_RANDOM_WALK_NUM              (15)

/** @brief Denominator to random walk if no cost reducing literals */
#define DEFAULT_RANDOM_WALK_DEN              (100)

/** @brief The default probability that the maximum weight neighboring clause
 *         is disregarded for a random clause.
 *
 *  TODO configurable?
 */
#define DEFAULT_RANDOM_CLAUSE_PROB          (0.01)

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


/** @brief Defines whether the algorithm will timeout due to time or
 *         number flips.
 */
typedef enum timeout_method {
  TIME, FLIPS
} timeout_method_t;

/** @brief Defines the selection method for the flipped variable.
 *
 *  UNIFORM is a uniform distribution across the candidate variables.
 *  WEIGHTED is a weighted probability distribution according to weights.
 *  BEST is the best weighted variable.
 */
typedef enum variable_selection_method {
  UNIFORM, WEIGHTED, BEST
} selection_method_t;

///////////////////////////////////////////////////////////////////////////////
// STATIC GLOBAL VARIABLE DECLARATIONS
///////////////////////////////////////////////////////////////////////////////

/** @brief What run number the program is on. */
static int runs;

/** @brief Timeout method. By default, times out by time. */
static timeout_method_t timeout_method = TIME;

/** @brief Variable selection method. By default, weighted probabilities. */
static selection_method_t selection_method = WEIGHTED;

/** @brief Number of seconds before instance timeout.
 *
 *  Can be toggled with the -t <secs> flag when running the executable.
 *  The amount of time elapsed is checked every DEFAULT_LOOPS_PER_TIMEOUT_CHECK 
 *  loops of the main loop body.
 */
static int timeout_secs = DEFAULT_TIMEOUT_SECS;

/** @brief Number of flips before instance timeout.
 *
 *  Can be toggled with the -T <flips> flag when running the executable.
 */
static int timeout_flips = -1;

/** @brief The probability that, instead of transferring weight from the
 *         maximum weighted neighbor, weight is transferred from a randomly
 *         chosen satisfied neighboring clause.
 */
static int random_clause_prob = (int) (1.0 / DEFAULT_RANDOM_CLAUSE_PROB);

/** @brief The multiplicative constant in weight transferral. */
static double a = A;

/** @brief The additive constant in weight transferral. */
static double c = C;

/** @brief The multiplicative constant when underneath w_init. */
static double a_prime = A;

/** @brief The additive constant when underneath w_init. */
static double c_prime = C;

#ifdef DEBUG
static int *cost_reducing_idxs_copy = NULL;
#endif

////////////////////////////////////////////////////////////////////////////////
// DDFW algorithm implementation
////////////////////////////////////////////////////////////////////////////////

#ifdef DEBUG
/** @brief Verifies the computed cost reducing vars by looping over all
 *         variables, instead of taking "cost_compute_vars" at its word.
 */
static void verify_cost_reducing_vars() {
  // First step, copy over the cost reducing indexes array to compare later
  memcpy(cost_reducing_idxs_copy, cost_reducing_idxs, 
      (num_vars + 1) * sizeof(int));
  memset(cost_reducing_idxs, 0xff, (num_vars + 1) * sizeof(int));

  // Next, clear number of cost reducing vars and loop through all vars
  // TODO this is just copied from below, pull out into helper function?
  const int prev_num_cost_reducing = num_cost_reducing_vars;
  num_cost_reducing_vars = 0;
  for (int v_idx = 1; v_idx <= num_vars; v_idx++) {
    const int l_idx = LIT_IDX(v_idx);
    const int assigned = ASSIGNMENT(l_idx);
    double satisfied_weight = 0.0;
    double unsatisfied_weight = 0.0;

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
    if (diff > 0.0) {
      add_cost_reducing_var(v_idx, diff);
    } else if (diff <= 0.0) {
      remove_cost_reducing_var(v_idx);
    }
  }

  // Now compare the two indexes array
  if (prev_num_cost_reducing != num_cost_reducing_vars) {
    printf("Number of cost reducing variables not the same, prev %d, now %d\n",
        prev_num_cost_reducing, num_cost_reducing_vars);
    exit(-1);
  }

  for (int i = 1; i <= num_vars; i++) {
    if ((cost_reducing_idxs[i] == -1 && cost_reducing_idxs_copy[i] != -1) ||
        (cost_reducing_idxs[i] != -1 && cost_reducing_idxs_copy[i] == -1)) {
      printf("Arrays differed at var %d\n", i);
      exit(-1);
    }
  }
}
#endif

/** @brief Finds those literals that cause a positive change in weighted
 *         cost if flipped and stores them in reducing_cost_lits.
 *         Stores the number of such literals in num_reducing_cost_lits.
 *
 *  According to the DDFW algorithm, the first thing done per each loop body
 *  is to "find and return a list L of literals causing the greatest reduction
 *  in weighted cost delta(w) when flipped."
 *
 *  While the list L could be computed by looping over every literal, every
 *  loop, that would result in much redundant computation, as only those
 *  literals involved in a clause with a flipped literal would have their
 *  cost calculations change.
 *
 *  As a result, flipping literals causes literal indexes to be stored in 
 *  cost_compute_vars, which acts as a sort of stack. This function pops
 *  every variable index in the stack and re-computes its delta(w). If
 *  delta(w) is strictly positive, the variable is added to the list
 *  of cost reducing variables.
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
      add_cost_reducing_var(v_idx, diff); // Filter out 0 cost variables?
    } else if (diff <= 0.0) {
      remove_cost_reducing_var(v_idx);
    }

    // Remove var from compute list
    remove_cost_compute_var(v_idx);
    cc_vars--;
  }

#ifdef DEBUG
  verify_cost_reducing_vars();
#endif
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
  // Check if weight is under init_clause_weight
  double orig = clause_weights[from_idx];

  // Use prime version of variables if under init weight
  if (orig < init_clause_weight) {
    clause_weights[from_idx] *= a_prime;
    clause_weights[from_idx] += c_prime;
  } else {
    clause_weights[from_idx] *= a;
    clause_weights[from_idx] += c;
  }

  // Whatever weight was taken away is given to unsat clause
  double diff = orig - clause_weights[from_idx];
  clause_weights[to_idx] += diff;

  // Assume the "from" clause is satisfied - then only the satisfied
  //   literal inside must have its cost change re-computed
  // assert(clause_num_true_lits[from_idx] > 0);
  if (clause_num_true_lits[from_idx] == 1) {
    const int mask_lit = clause_lit_masks[from_idx];
    add_cost_compute_var(VAR_IDX(mask_lit));
  }

  // In the unsatisfied "to" clause, add vars to cost compute
  const int to_size = clause_sizes[to_idx];
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
  // Loop over all clauses, picking out those that are false
  const int uc = num_unsat_clauses;
  int *false_clauses = false_clause_members;
  for (int c = 0; c < uc; c++) {
    const int c_idx = *false_clauses;
    const int size = clause_sizes[c_idx];

    // Scan for any neighboring clause that is satisfied and has max weight
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

        // Move to the next clause index for this literal
        l_clauses++;
      }

      // Move to the next literal in this clause
      cl_lits++;
    }

    // If a maximum weight neighbor has been found for this clause
    if (max_neighbor_idx != -1) {
      // Transfer weight with almost 1 probability
      if (rand() % random_clause_prob != 0) {
        transfer_weight(max_neighbor_idx, c_idx);
      } else {
        // Select random satisfied clause instead by choosing random indexes
        // NOTE: This was not mentioned in the original paper, but instead
        // found in the original code by the authors. The form it takes in
        // the UBCSAT code is to take any random clause, rather than a
        // neighboring same-sign clause. Potential TODO
        unsigned int r_idx;
        do {
          r_idx = ((unsigned int) rand()) % num_clauses;
        } while (clause_num_true_lits[r_idx] == 0); // TODO weight case

        transfer_weight(r_idx, c_idx);
      }
    }

    // Move to the next false clause in the array
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
  struct timeval start_time, stop_time;

  gettimeofday(&start_time, NULL);

  int var_to_flip;
  while (num_unsat_clauses > 0) {
    find_cost_reducing_literals();
    
    // See if any literals will reduce the weight
    if (num_cost_reducing_vars > 0) {
      unsigned int rand_var = 0;
      switch (selection_method) {
        case UNIFORM:
          rand_var = ((unsigned int) rand()) % num_cost_reducing_vars;
          break;
        case WEIGHTED:;
          double total_so_far = 0.0;
          double rand_double = ABS(((double) rand()) 
              / ((double) RAND_MAX) * total_cost_reducing_weight);
          double *ws = cost_reducing_weights;

          // Loop through the weights until found total greater than rand_d
          for (int i = 0; i < num_cost_reducing_vars; i++) {
            total_so_far += *ws;
            if (rand_double <= total_so_far) {
              rand_var = i;
              break;
            }

            ws++;
          }
          break;
        case BEST:;
          double best_weight = DBL_MAX;
          int best_idx = 0;
          ws = cost_reducing_weights;
          for (int i = 0; i < num_cost_reducing_vars; i++) {
            const double w = *ws;
            if (w < best_weight) {
              best_weight = w;
              best_idx = i;
            }
          }

          rand_var = best_idx;
          break;
        default:
          fprintf(stderr, "Unrecognized selection method\n");
          exit(-1);
      }

      var_to_flip = cost_reducing_vars[rand_var];
    } else if (((unsigned int) rand()) % DEFAULT_RANDOM_WALK_DEN <
        DEFAULT_RANDOM_WALK_NUM) {
      var_to_flip = ((unsigned int) rand()) % num_vars;
    } else {
      distribute_weights();
      continue;
    }

    flip_variable(var_to_flip);

    // Determine if enough flips/loops have passed to update time variable
    if (timeout_method == FLIPS && num_flips >= timeout_flips) {
      break;
    } else {
      timeout_loop_counter--;
      if (timeout_loop_counter == 0) {
        timeout_loop_counter = DEFAULT_LOOPS_PER_TIMEOUT_CHECK;
        gettimeofday(&stop_time, NULL);

        log_str("c %d unsatisfied clauses after %d flips\n", 
            num_unsat_clauses, num_flips);

        // Check for timeout
        if (stop_time.tv_sec - start_time.tv_sec >= timeout_secs)
          break;
      }
    }
  }

  gettimeofday(&stop_time, NULL);

  // Print solution
  output_assignment();
  log_statistics(runs, &start_time, &stop_time);
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
  if (argc == 1) {
    print_usage(argv[0]);
    exit(1);
  }
 
  int seed = DEFAULT_SEED;
  char *filename = NULL;
  extern char *optarg;
  char opt;
  while ((opt = getopt(argc, argv, "dhvqa:A:c:C:f:m:r:s:t:T:w:")) != -1) {
    switch (opt) { 
      case 'a':
        a = atof(optarg);
        log_str("c Multiplicative scalar is %lf\n", a);
        break;
      case 'A':
        a_prime = atof(optarg);
        log_str("c Alternative multiplicative scalar is %lf\n", a_prime);
        break;
      case 'c':
        c = atof(optarg);
        log_str("c Additive constant is %lf\n", c);
        break;
      case 'C':
        c_prime = atof(optarg);
        log_str("c Alternative additive constant is %lf\n", c_prime);
        break;
      case 'd':
        // Original DDFW paper settings
        init_clause_weight = 8.0;
        a = 1.0;
        a_prime = 1.0;
        c = -2.0;
        c_prime = -1.0;
        break;
      case 'f':
        filename = optarg;
        break;
      case 'h':
        print_help(argv[0]);
        return 0;
      case 'm':
        switch (optarg[0]) {
          case 'U':
            selection_method = UNIFORM;
            log_str("c Selection method is uniform distribution\n");
            break;
          case 'W':
            selection_method = WEIGHTED;
            log_str("c Selection method is weighted distribution\n");
            break;
          case 'B':
            selection_method = BEST;
            log_str("c Selection method is best\n");
            break;
          default:
            selection_method = WEIGHTED;
            log_str("c Unrecognized selection method, default is weighted\n");
            break;
        }
        break;
      case 'q':
        set_verbosity(SILENT);
        break;
      case 'r':
        num_restarts = atoi(optarg);
        log_str("c Will run the algorithm %d times\n", num_restarts);
        break;
      case 's':
        seed = atoi(optarg);
        if (seed != 0) {
          log_str("c Using randomization seed %d\n", seed);
          srand(seed);
        }
        break;
      case 't':
        timeout_secs = atoi(optarg);
        log_str("c Timeout set to %d\n", timeout_secs);
        break;
      case 'T':
        timeout_flips = atoi(optarg);
        timeout_method = FLIPS;
        log_str("c Timeout set to %d flips\n", timeout_flips);
        break;
      case 'v':
        set_verbosity(VERBOSE);
        break;
      case 'w':
        init_clause_weight = atof(optarg);
        log_str("c Default clause weight set to %lf\n", init_clause_weight);
        break;
      default:
        print_usage(argv[0]);
    }
  }

  log_str("c ------------------------------------------------------------\n");
  log_str("c          DDFW - Divide and Distribute Fixed Weights\n");
  log_str("c                  Implemented by Cayden Codel\n");
  log_str("c                          Version 0.1\n");
  log_str("c ------------------------------------------------------------\n");

  if (filename == NULL) {
    fprintf(stderr, "c No filename provided, exiting.\n");
    exit(1);
  }

  if (timeout_method == TIME) {
    log_str("c Running with timeout of %d seconds\n", timeout_secs);
  } else {
    log_str("c Running with timeout of %d flips\n", timeout_flips);
  }

  switch (selection_method) {
    case BEST:
      log_str("c Running with selection method BEST\n");
      break;
    case UNIFORM:
      log_str("c Running with selection method UNIFORM\n");
      break;
    case WEIGHTED:
    default:
      log_str("c Running with selection method WEIGHTED\n");
      break;
  }

  if (seed == DEFAULT_SEED) {
    log_str("c Using default randomization seed %d\n", seed);
    srand(seed);
  }

  log_str("c All done parsing CLI args, opening file %s\n", filename);
  parse_cnf_file(filename);

#ifdef DEBUG
  cost_reducing_idxs_copy = xmalloc((num_vars + 1) * sizeof(int));
#endif

  for (runs = 1; runs <= num_restarts; runs++) {
    run_algorithm();
    reset_data_structures();
  }

  return 0;
}
