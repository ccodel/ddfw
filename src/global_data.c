/** @file global_data.h
 *  @brief Global data structures shared across all files.
 *
 *  A compilation of all global data shared between the DDFW files. File-
 *  specific data structures and supporting variables should be declared
 *  within the .c file.
 *
 *  The data structures presented here should be data structures that are
 *  common across all DDFW implementations, i.e. structures that optimize
 *  should be relegated to other files.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#include <string.h>

#include "global_data.h"
#include "initializer.h"
#include "xmalloc.h"

#include <stdio.h>

/** @brief The number of variables in the CNF formula.
 *
 *  Each CNF formula (in DIMACS format) is comprised of clauses of disjunctions
 *  joined by conjuctions. The Boolean literals of each clause are either
 *  positive or negative. The number of distinct variable symbols (with
 *  positive and negative literals being considered as the same variable)
 *  is stored here when parsing the specific CNF file (see cnf_parser.c).
 */
int num_vars;

/** @brief The number of literals in the CNF formula.
 *  
 *  Generally, this is twice the number of variables (2 * num_vars).
 */
int num_literals;

/** @brief The number of clauses in the CNF formula.
 *
 *  Each CNF formula (in DIMACS format) is comprised of disjunctive clauses.
 *  This variable stores how many clauses there are in the formula when
 *  being parsed (see cnf_parser.c).
 */
int num_clauses;

/** @brief The initial weight to give all clauses.
 *
 *  The DDFW algorithm assigns a weight to each clause and distributes weights
 *  during the course of the algorithm. This value, specified at the command
 *  line, is the value of initial weight assigned to all clauses at the start.
 */
weight init_clause_weight;

/** @brief The total amount of weight held by satisfied clauses. */
weight total_sat_clause_weight;

/** @brief The total amount of weight held by unsatisfied clauses. */
weight total_unsat_clause_weight;

/** @brief The number of Boolean variable flips performed by the current run.
 *
 *  The DDFW algorithm incrementally flips Boolean variables until timeout or
 *  until a satisfying assignment is found. This value stores the current
 *  number of flips the algorithm has reached.
 */
long num_flips = 0;

/** @brief The number of Boolean variable flips since getting a better
 *         assignment (a lower number of unsatisfied clauses).
 *
 *  The DDFW+ algorithm uses this values to do re-weighting.
 */
long num_flips_since_improvement = 0;

/** @brief The true/false values to the Boolean variables.
 *
 *  Indexed by variable number (VAR_IDX).
 *
 *  TODO make into a more compact bit array?
 */
char *assignment = NULL;

/** @brief The number of unsatisfied clauses under the current assignment. */
int num_unsat_clauses;

/** @brief 0-indexed array that stores claues IDs which have no true lit. */
int *unsat_clauses;

/** @brief The best true/false values to the Boolean variables found so far.
 *
 *  Whenever an assignment is found which lowers the number of unsatisfied
 *  clauses below the best found so far, that assignment is copied here.
 *  Some DDFW algorithms may restart to the best assignment previously found.
 */
char *best_assignment = NULL;

/** @brief The number of unsatisfied clauses in the best assignment found. */
int best_num_unsat_clauses;

/** @brief The flip number on which the best assignment was found. */
long best_flip_num;

/** @brief Stores the sizes of each clause.
 *
 *  The array is 0-indexed by clause ID.
 */
int *clause_sizes = NULL;

/** @brief Stores the weights held by each clause.
 *
 *  The array is 0-indexed by clause ID. Initialized with all the same weight.
 */
weight *clause_weights = NULL;

/** @brief Stores the number of true literals in each clause.
 *
 *  The array is 0-indexed by clause ID.
 */
int *clause_num_true_lits = NULL;

/** @brief Stores the mask of true literal IDs for each clause.
 *
 *  When a clause only has a single literal that is true, it is often important
 *  to know which literal that is. XORing all the literal IDs together
 *  provides for a way to recover the ID of the last remaining satisfying
 *  literal.
 *
 *  The array is 0-indexed by clause ID.
 */
int *clause_lit_masks = NULL;

/** @brief Stores the literals in each clause.
*
*  The array is 0-indexed by clause ID. The values stored in the array
*  are literal IDs (see LIT_IDX).
*/
int **clause_literals = NULL;

/** @brief Stores how many times each literal occurs.
*
*  The array is 1-indexed by literal ID (see LIT_IDX).
*/
int *literal_occ = NULL;

/** @brief Stores the clauses that each literal is a member in.
*
*  The array is 1-indexed by literal ID (see LIT_IDX). The values in the
*  array are clause IDs.
*/
int **literal_clauses = NULL;

/*************************** FUNCTIONS ****************************************/

/** @brief Allocates memory needed by global data structures.
 *
 *  This function allocates memory needed based on the size of the CNF formula
 *  being considered, and so should be called after the header line of the
 *  DIMACS CNF file has been parsed by cnf_parser.c.
 */
void allocate_global_memory(void) {
  const int var_arr_size = num_vars + 1;
  const int lit_arr_size = num_literals + 2;

  assignment = xmalloc(var_arr_size * sizeof(char));
  best_assignment = xmalloc(var_arr_size * sizeof(char));
  unsat_clauses = xmalloc(num_clauses * sizeof(int));

  clause_sizes = xmalloc(num_clauses * sizeof(int));
  clause_weights = xmalloc(num_clauses * sizeof(weight));
  clause_num_true_lits = xmalloc(num_clauses * sizeof(int));
  clause_lit_masks = xcalloc(num_clauses, sizeof(int));
  clause_literals = xcalloc(num_clauses, sizeof(int *));

  literal_occ = xcalloc(lit_arr_size, sizeof(int));
  literal_clauses = xcalloc(lit_arr_size, sizeof(int *));
}

void initialize_global_literals_to_clauses(void) {
  // Allocate a helper count buffer
  int *clause_counts = xcalloc((num_literals + 2), sizeof(int));

  // Loop through the clauses to add clause IDs to each literal in the clause
  for (int c = 0; c < num_clauses; c++) {
    const int clause_size = clause_sizes[c];
    int *lits = clause_literals[c];

    // For each literal in the clause, add the clause ID to the array
    for (int l = 0; l < clause_size; l++) {
      const int l_idx = lits[l];
      if (literal_clauses[l_idx] == NULL) {
        literal_clauses[l_idx] = xmalloc(literal_occ[l_idx] * sizeof(int));
      }

      literal_clauses[l_idx][clause_counts[l_idx]] = c;
      clause_counts[l_idx]++;
    }
  }

  xfree(clause_counts);
}

static void reset_clause_weights_local(void) {
  weight *end = clause_weights + num_clauses;
  for (weight *w = clause_weights; w < end; w++) {
    *w = init_clause_weight;
  }
}

/** @brief Resets the weights held by the clauses to the initial weight.
 *
 *  Some DDFW implementations reset the weights held by the clauses (e.g. in
 *  order to escape local minimums). Clause weights are set to the initial
 *  weight, as specified at the command line.
 */
void reset_clause_weights(void) {
  reset_clause_weights_local();
  initialize_structures_after_reweighting();
}

/** @brief Sets the weights held by all unsatisfied clauses to a provided value.
 * 
 *  Once done, the initializer handles recomputing structures with new weights.
 *
 *  @param w  The weight value that all unsatisfied clauses should take.
 */
void set_unsat_weights(weight w) {
  total_unsat_clause_weight = w * num_unsat_clauses;
  weight *end = clause_weights + num_clauses;
  int *true_lits_ptr = clause_num_true_lits;
  for (weight *wptr = clause_weights; wptr < end; wptr++) {
    if (*true_lits_ptr == 0) {
      *wptr = w;
    }

    true_lits_ptr++;
  }

  initialize_structures_after_reweighting();
}

/** @brief Increases the weights held by all unsatisfied clauses by a
 *         provided value. If the value is negative, the weights are instead
 *         decreased.
 *
 *  If decreasing a clause's weight would take it below 0, weight is
 *  decreased to 0 instead.
 *
 *  Once done, the initializer handles recomputing structures with new weights.
 *
 *  @param w  The weight value that all unsatisfied clauses should increase by.
 */
void increase_unsat_weights(weight w) {
  if (w == 0)
    return;

  weight *end = clause_weights + num_clauses;
  int *true_lits_ptr = clause_num_true_lits;
  for (weight *wptr = clause_weights; wptr < end; wptr++) {
    if (*true_lits_ptr == 0) {
      const int new_weight = *wptr + w;
      if (new_weight < 0) {
        total_unsat_clause_weight -= *wptr;
        *wptr = 0;
      } else {
        total_unsat_clause_weight += w;
        *wptr = new_weight;
      }
    }

    true_lits_ptr++;
  }

  initialize_structures_after_reweighting();
}

/** @brief Sets the weights held by all satisfied clauses to a provided value.
 *
 *  Once done, the initializer handles recomputing structures with new weights.
 *
 *  @param w  The weight value that all satisfied clauses should take.
 */
void set_sat_weights(weight w) {
  total_sat_clause_weight = w * num_unsat_clauses;
  weight *end = clause_weights + num_clauses;
  int *true_lits_ptr = clause_num_true_lits;
  for (weight *wptr = clause_weights; wptr < end; wptr++) {
    if (*true_lits_ptr > 0) {
      *wptr = w;
    }

    true_lits_ptr++;
  }

  initialize_structures_after_reweighting();
}

/** @brief Increases the weights held by all satisfied clauses by a
 *         provided value. If the value is negative, the weights are instead
 *         decreased.
 *
 *  Once done, the initializer handles recomputing structures with new weights.
 *
 *  @param w  The weight value that all satisfied clauses should increase by.
 */
void increase_sat_weights(weight w) {
  if (w == 0)
    return;

  weight *end = clause_weights + num_clauses;
  int *true_lits_ptr = clause_num_true_lits;
  for (weight *wptr = clause_weights; wptr < end; wptr++) {
    if (*true_lits_ptr > 0) {
      const int new_weight = *wptr + w;
      if (new_weight < 0) {
        total_sat_clause_weight -= *wptr;
        *wptr = 0;
      } else {
        total_sat_clause_weight += w;
        *wptr = new_weight;
      }
    }

    true_lits_ptr++;
  }

  initialize_structures_after_reweighting();
}

/** @brief Resets all data structures managed by global_data.c.
 *
 *  Resets all data structures managed by global_data.c to default values, e.g.
 *  values fresh for a new algorithm run.
 */
void initialize_global_data_for_alg_run(void) {
  // Set single variables first
  total_sat_clause_weight = 0;
  total_unsat_clause_weight = 0;
  num_flips = 0;
  num_flips_since_improvement = 0;
  num_unsat_clauses = 0;
  best_num_unsat_clauses = 0;
  best_flip_num = 0;

  // Clear arrays
  memset(assignment, 0, (num_vars + 1) * sizeof(char));
  memset(best_assignment, 0, (num_vars + 1) * sizeof(char));
  memset(clause_num_true_lits, 0, num_clauses * sizeof(int));
  memset(clause_lit_masks, 0, num_clauses * sizeof(int));
  reset_clause_weights_local();
}
