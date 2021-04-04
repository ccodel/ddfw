/** @file ddfw.c - Divide and Distribute Fixed Weights
 *  @brief An implementation of the DDFW algorithm, as presented in 
 *         the paper: "Neighbourhood Clause Weight Redistribution in 
 *         Local Search for SAT" [Ishtaiwi, Thornton, Sattar, Pham].
 *
 *  Note that this implementation of DDFW assumes that:
 *   - No clause has more than MAX_CLAUSE_LENGTH literals in it
 *   - No holes exist in the variables in the CNF (undefined behavior)
 *
 *  Note that this implementation may perform best after:
 *   - Preprocessing of single literal clauses
 *  
 *  ///////////////////////////////////////////////////////////////////////////
 *  // HIGH LEVEL OVERVIEW OF DDFW
 *  ///////////////////////////////////////////////////////////////////////////
 *  
 *  DDFW is a stochastic local search (SLS) SAT solver that attempts to solve
 *  satisfiable CNF boolean formulas by iteratively flipping boolean variables
 *  until either a solution is found or timeout is reached. When the CNF
 *  formula is read in, an initial starting weight is given to each clause.
 *  For each variable flip, a variable that can reduce the total amount of
 *  weight held by unsatisfiable clauses is chosen. If no such variable exists,
 *  then weight is distributed from satisfied clauses to unsatisfied clauses.
 *  
 *  Unlike the original DDFW implementation, which used integers for clause
 *  weights only, this implementation of DDFW uses floating-point nubmers
 *  for clause weights. Floating-point numbers allow for greater expressivity
 *  in the weight re-distribution policy used when no weight-reducing variables
 *  are left to flip.
 *
 *  The body of the algorithm is implemented in this file, but helper files
 *  contain important data structures, parsing functions, and logging utilities.
 *  Below is a listing of each helper file and its supporting purpose:
 *
 *    - clause:  Data structures for literals, clauses, assignment, etc.
 *    - cnf_parser: Parses the provided CNF file and pulls out clause info.
 *    - logger:  Logs strings and data structures according to verbosity level.
 *    - xmalloc: Guaranteed memory allocation. Wrappers around exit().
 *  
 *  See below for a high-level description of various implementation features.
 *  See each file listed above for a more in-depth description and for code.
 *
 *  ///////////////////////////////////////////////////////////////////////////
 *  // IMPLEMENTATION
 *  ///////////////////////////////////////////////////////////////////////////
 *  
 *  DDFW solves boolean formulas in conjunctive normal form (CNF), meaning that
 *  boolean variables are grouped in disjunctive clauses (ORs) joined by ANDs.
 *  The format of the CNF formula provides many opportunities for optimization
 *  and "tricks of the trade" to minimize the amount of work done at each
 *  variable flip, keeping the implementation efficient.
 *
 *  In spirit, literals are represented by the struct
 *
 *    typedef struct boolean_literal {
 *      int var_symbol;
 *      int occurrences;
 *      clause_t **clauses;
 *    } literal_t;
 *
 *  where each literal stores its variable symbol from the CNF file, the 
 *  number of clauses it occurs in, and an array of clauses that contain 
 *  it. Similarly, clauses are represented, in spirit, by the struct
 *
 *    typedef struct disjunctive_clause {
 *      int id;
 *      int size;
 *      int num_true_lits;
 *      int sat_lit_mask;
 *      double weight;
 *      literal_t **literals;
 *    } clause_t;
 *
 *  where each clause stores which line from the CNF file it lies on, the
 *  number of literals it contains, how many of those literals are true under 
 *  the current assignment, an XOR of its satisfying literals to easily
 *  identify the symbol of the last remaining satisfying literal, its weight,
 *  and an array of literals contained in it.
 *
 *  However, implementing the above data structures caused two orders of
 *  magnitude slowdown due to pointer arithmetic inside and across structs.
 *  Therefore, the above data structures are "split up" into their component
 *  fields as individual arrays indexed by literal or clause ID. While doing
 *  so makes the codebase harder to manage, the speedup is necessary to have
 *  a practical SAT solver useful for anything more than trivial formulas.
 *  The breakup of the data structures is handled by clause.h/.c.
 *
 *  The main body of DDFW is implemented in run_ddfw_algorithm(). After noting
 *  the system time, the loop will continually flip variables and distribute
 *  weight between clauses until a solution is found or until one of the two
 *  timeout conditions is met (CPU clock time or number of flips). To minimize
 *  the number of system calls, the time is not checked on every loop, but
 *  once every LOOPS_PER_TIMEOUT_CHECK loops.
 *
 *  In choosing which variable to flip, one of three selection methods are
 *  used: UNIFORM, WEIGHTED, and BEST. The UNIFORM method will select any
 *  cost-reducing varaible with equal probability. The WEIGHTED method will
 *  select variables with a weighted probability according to the amount of
 *  weight reduced by flipping that variable. The BEST method will select
 *  the variable with the highest weight reduction.
 *
 *  If no cost-reducing variables are present, then for each unsatisfied
 *  clause, weight is distributed from one of its adjacent satisfied clauses
 *  with the hightest weight. Occasionally, a random satisfied clause is
 *  chosen instead. Weight is distributed according to a linear rule:
 *
 *    Let c be the current unsatisfied clause, and let n be its highest weight
 *    neighbor. Also, INIT_WEIGHT is the initial weight given to each clause,
 *    "a" and "A" are mutliplicative constants, and "c" and "C" are additive
 *    constants. Then weight is re-assigned via
 *
 *    double diff;
 *    if (n->weight < INIT_WEIGHT) diff = n->weight - (a * n->weight + c);
 *    else                         diff = n->weight - (A * n->weight + C);
 *    c->weight += diff;
 *    n->weight -= diff;
 *
 *  By adjusting a, A, c, and C at the command line, different linear rules
 *  can be applied to the weights of the clauses. Setting a = A and c = C
 *  removes the effect of the if-statement.
 *
 *  Additional details regarding command-line arguments, data structures,
 *  and other implementation details can be found in the appropriate files
 *  or above each function.
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
 *  In this implementation, each literal has a literal_t struct "in spirit" that
 *  is allocated across the split-up arrays in clause.h/.c. Literals are stored
 *  adjacently to their negated form. So, if the literal is (l = 5), then
 *  the index into the literals array to access the literal_t struct for
 *  l is (2 * 5), and the index for !l is (2 * 5) + 1. Because DIMACS variables
 *  are 1-indexed, each array containing literal information has "two too many" 
 *  elements to preserve the 1-indexing of DIMACS format.
 *
 *  The clauses are 0-indexed, however, across their split-up arrays.
 *
 *  The motivation behind using a [l1, !l1, l2, !l2, ...] literals storage
 *  structure, as opposed to all positive literals on one side of a pointer
 *  and all negated literals on the other, is to take advantage of caching
 *  when accessing l and !l together. For example, when flipping a literal
 *  in the assignment, clauses for both l and !l must be updated, and so
 *  both literal_t structs are accessed together "in spirit." Whether or not 
 *  such caching effects can be felt during runtime is something the author 
 *  has not investigated, and stands as a potential optimization.
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
 *  Potential improvement according to
 *   https://www.keithschwarz.com/darts-dice-coins/
 *
 *  The WEIGHTED and BEST probability distributions are found linearly, which
 *  can cause slowdown. A tree or segtree data structure may be best here.
 *
 *  @bugs No known bugs.
 *
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
#include <signal.h>

#include "ddfw.h"
#include "clause.h"
#include "cnf_parser.h"
#include "logger.h"
#include "xmalloc.h"
#include "neighborhood.h"
#include "weight_transfer.h"
#include "verifier.h"

///////////////////////////////////////////////////////////////////////////////
// CONFIGURATION MACROS
///////////////////////////////////////////////////////////////////////////////

/** @brief Default number of seconds to loop before timing out. 
 *  
 *  Can be toggled with the -t <timeout_secs> flag at the command line.
 **/
#define DEFAULT_TIMEOUT_SECS        100

/** @brief The number of times the main DDFW loop body is run before the number 
 *         of elapsed seconds is checked. Not currently configurable.
 *
 *  TODO maybe define this based on some computed runtime quantity? Make
 *       smaller based on number of clauses/literals?
 */
#ifdef DEBUG
#define LOOPS_PER_TIMEOUT_CHECK     1000
#else
#define LOOPS_PER_TIMEOUT_CHECK     1000
#endif

/** @brief Default multiplicative and additive constants.
 *
 *  When transferring weights, weight is transferred according to a linear
 *  rule, highlighted at the top of the file. DEFAULT_A is the default 
 *  multiplicative constant for both "a" and "A," and DEFAULT_C is the
 *  default additive constant for both "c" and "C."
 *
 *  The values of 1.0 and -2.0 were chosen to match the original DDFW paper.
 */
#define DEFAULT_A   (1.0)
#define DEFAULT_C   (-2.0)

/** @brief The probability with which to select a random variable to flip 
 *         when no cost-reducing variables remain.
 *
 *  Because the pseudo-random number generator only returns integers, an
 *  integer representation for the random walk probability was chosen. The
 *  random number is taken mod DEN, and compared to be less than NUM.
 *
 *  In the original DDFW paper, the random walk probability was 15%.
 */
#define RANDOM_WALK_NUM       15
#define RANDOM_WALK_DEN       100

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

/** @brief What run number the algorithm is on. */
int algorithm_run;

/** @brief Timeout method. Times out after DEFAULT_TIMEOUT_SECS seconds.
 *
 *  Can be toggled with the -t and -T command-line options.
 */
timeout_method_t timeout_method = DEFAULT;

/** @brief Variable selection method. Selects variable by weighted probability.
 *
 *  Can be toggled with the -m command-line option.
 */
selection_method_t selection_method = WEIGHTED;

/** @brief Number of seconds before instance timeout.
 *
 *  Can be toggled with the -t <secs> command-line option. The elapsed time
 *  is checked every LOOPS_PER_TIMEOUT_CHECK loops of the main loop body.
 */
int timeout_secs = DEFAULT_TIMEOUT_SECS;

/** @brief Number of flips before instance timeout.
 *
 *  Can be toggled with the -T <flips> flag when running the executable.
 *  When the timeout_method is DEFAULT, the number of flips is not checked.
 */
int timeout_flips = -1;

/** @brief The multiplicative constant in weight transferral. 
 *
 *  When the weight of the neighboring clause is at least init_clause_weight,
 *  then "a" is used. Can be specified with the -a <a> command-line option.
 */
double mult_a = DEFAULT_A;

/** @brief The additive constant in weight transferral.
 *
 *  When the weight of the neighboring clause is at least init_clause_weight,
 *  then "c" is used. Can be specified with the -c <c> command-line option.
 */
double add_c = DEFAULT_C;

/** @brief The multiplicative constant for when beneath init_clause_weight.
 *  
 *  When the weight of the neighboring clause is less than init_clause_weight,
 *  then "A" is used. Can be specified with the -A <A> command-line option.
 */
double mult_A = DEFAULT_A;

/** @brief The additive constant for when beneath init_clause_weight.
 *
 *  When the weight of the neighboring clause is less than init_clause_weight,
 *  then "C" is used. Can be specified with the -C <C> command-line option.
 */
double add_C = DEFAULT_C;

/** @brief Suppresses the printing of a solution if found. 
 * 
 *  TODO move this into the logger
 */
int suppress_solution = 0;

/** @brief The number of flips between logging weight statistics.
 *
 *  By default, weight statistics are not logged while DDFW is running. If the
 *  -l <flips> option is specified at the command-line, then every <flips>
 *  flips, the weight statistics are logged with a call to the function
 *  log_weight_statistics().
 */
int weight_statistics_log_rate = 0;

// Additional options for weight transfer rules
transfer_group_t transfer_grp = SINGULAR;
rule_group_t rule_grp = RULE_SINGULAR;
available_weight_option_t available_weight_opt = RAW;

///////////////////////////////////////////////////////////////////////////////
// DDFW algorithm implementation
///////////////////////////////////////////////////////////////////////////////

/** @brief Finds the literals that cause a decrease in the total amount of
 *         weight held by unsatisfied clauses if flipped and stores their
 *         indexes into the array cost_reducing_vars.
 *
 *  According to the DDFW algorithm, the first thing done per each loop body
 *  is to "find and return a list L of literals causing the greatest reduction
 *  in weighted cost delta(w) when flipped." While the list L could be computed 
 *  by looping over every literal, every loop, that would result in much 
 *  redundant computation, as only those literals involved in a clause with a 
 *  flipped literal would have their cost calculations change.
 *
 *  As a result, flipping literals causes variable indexes to be stored in 
 *  cost_compute_vars, which acts as a sort of stack. This function pops off
 *  every variable index in the stack and re-computes its delta(w). If
 *  delta(w) is strictly positive, the variable index is added to the list
 *  of cost reducing variables.
 */
static void find_cost_reducing_literals() {
  // Loop through those variables in the cost_compute_vars
  // Most efficient to do from the back of the cost compute vars array
  int *cc_vars = cost_compute_vars + num_cost_compute_vars - 1;
  const int num_cc_vars = num_cost_compute_vars;
  for (int i = 0; i < num_cc_vars; i++) {
    const int v_idx = *cc_vars;
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

    // Check clauses containing the true literal for ones that can become unsat
    unsatisfied_weight = literal_crit_sat_weights[true_idx];

    // Check clauses containing the false literal for those that can become sat
    satisfied_weight = literal_unsat_weights[false_idx];

    // Determine if flipping the truth value of the true literal
    //   would result in more satisfied weight than unsatisfied weight
    const double diff = satisfied_weight - unsatisfied_weight;
    if (diff > 0.0) {
      add_cost_reducing_var(v_idx, diff);
    } else if (diff <= 0.0) {
      remove_cost_reducing_var(v_idx);
    }

    // Remove var from compute list
    remove_cost_compute_var(v_idx);
    cc_vars--;
  }

  verify_cost_reducing_vars();
}


/** @brief Runs the DDFW algorithm.
 *
 *  Main loop.
 */
void run_ddfw_algorithm() {

  generate_random_assignment();

  // Record the time to ensure no timeout
  int timeout_loop_counter = LOOPS_PER_TIMEOUT_CHECK;
  struct timeval start_time, stop_time;
  gettimeofday(&start_time, NULL);

  // Spit out statistics before any flips, if enabled
  if (weight_statistics_log_rate > 0) {
    log_weight_statistics(algorithm_run, 0, 0);
  }

  int var_to_flip;
  while (num_unsat_clauses > 0) {
    // As we loop, verify various components and data structures
    verify_clauses_and_assignment();
    verify_crit_sat_unsat_weights();

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
    } else if (((unsigned int) rand()) % RANDOM_WALK_DEN < RANDOM_WALK_NUM) {
      var_to_flip = ((unsigned int) rand()) % num_vars;
    } else {
      distribute_weights();
      goto check_timeout;
    }

    flip_variable(var_to_flip);

    // Log in-run statistics if the rate is reached
    if (weight_statistics_log_rate > 0 && 
        num_flips % weight_statistics_log_rate == 0) {
      double average = get_transfer_average();
      log_weight_statistics(algorithm_run, num_weight_transfers, average);
    }

    // Determine if enough flips/loops have passed to update time variable
    if ((timeout_method == FLIPS || timeout_method == BOTH) 
        && num_flips >= timeout_flips) {
      break;
    }

check_timeout:
    if (timeout_method != FLIPS) {
      timeout_loop_counter--;
      if (timeout_loop_counter == 0) {
        timeout_loop_counter = LOOPS_PER_TIMEOUT_CHECK;
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
  if (!suppress_solution) {
    output_assignment();
  }

  log_statistics(algorithm_run, &start_time, &stop_time);
}
