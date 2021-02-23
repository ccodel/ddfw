/** @file clause.c
 *  @brief Implements functions to create and interact with CNf DDFW clauses.
 *
 *  TODO
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#include <stdlib.h>
#include <string.h>

#include "clause.h"
#include "logger.h"
#include "xmalloc.h"
#include "neighborhood.h"

#include <stdio.h>

#ifndef ABS
#define ABS(x)   (((x) < 0) ? -(x) : (x))
#endif

#ifndef MIN
#define MIN(x, y)  (((x) < (y)) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x, y)  (((x) > (y)) ? (x) : (y))
#endif

// CNF information
int num_vars = 0;
int num_literals = 0;
int num_clauses = 0;
double init_clause_weight = DEFAULT_CLAUSE_WEIGHT;
double unsat_clause_weight = 0.0;

// Statistics
int num_restarts = 1;
long num_flips = 0;
int lowest_unsat_clauses = 0;
int lowest_unsat_step = 0;

// Formula information
char *assignment = NULL;

// Clause information - 0 indexed
int *clause_sizes = NULL;
double *clause_weights = NULL;
int *clause_num_true_lits = NULL;
int *clause_lit_masks = NULL;
int **clause_literals = NULL;

// Literal information - 2-indexed (use LIT_IDX)
int *literal_occ = NULL;
int **literal_clauses = NULL;

// Bookkeeping structures
// Membership struct for false clauses
int *false_clause_members = NULL;
int *false_clause_indexes = NULL;
int num_unsat_clauses = 0;

// Membership struct for cost reducing literals
/** @brief The array containing indexes of literals that, if flipped,
 *         reduce the cost of the unsatisfied literals.
 *
 *  The exact way that the list is calculated and a literal is chosen from
 *  within is uncertain in the paper. I have determined that those literals
 *  that reduce cost a strictly positive amount should be prioritized over
 *  those that result in no weight change. Therefore, those literals that
 *  reduce cost are placed "in the front" of the array, while those literals
 *  that result in no weight change on a flip are placed "in the back" of
 *  the array. The variables num_reducing_cost_lits and num_zero_cost_lits
 *  are updated accordingly.
 *
 *  TODO update for _lits and _membership
 */
int *cost_reducing_vars = NULL;
int *cost_reducing_idxs = NULL;
double *cost_reducing_weights = NULL;
double total_cost_reducing_weight = 0.0;
int num_cost_reducing_vars = 0;

int *cost_compute_vars = NULL;
int *cost_compute_idxs = NULL;
int num_cost_compute_vars = 0;

// TODO make this general / into a macro
static inline void add_false_clause(int c_idx) {
  if (false_clause_indexes[c_idx] == -1) {
    false_clause_indexes[c_idx] = num_unsat_clauses;
    false_clause_members[num_unsat_clauses] = c_idx;
    num_unsat_clauses++;
    unsat_clause_weight += clause_weights[c_idx];
  }
}

static inline void remove_false_clause(int c_idx) {
  const int idx = false_clause_indexes[c_idx];
  if (idx != -1) {
    num_unsat_clauses--;

    // If idx is not at the end, move the end of members to this index
    if (idx != num_unsat_clauses) {
      const int end_clause = false_clause_members[num_unsat_clauses];
      false_clause_members[idx] = end_clause;
      false_clause_indexes[end_clause] = idx;
    }
    
    false_clause_indexes[c_idx] = -1;
    unsat_clause_weight -= clause_weights[c_idx];
  }
}

inline void add_cost_reducing_var(const int v_idx, const double w) {
  const int idx = cost_reducing_idxs[v_idx];
  if (idx == -1) {
    cost_reducing_idxs[v_idx] = num_cost_reducing_vars;
    cost_reducing_vars[num_cost_reducing_vars] = v_idx;
    cost_reducing_weights[num_cost_reducing_vars] = w;
    total_cost_reducing_weight += w;
    num_cost_reducing_vars++;
  } else {
    // Update weight, change nothing else
    total_cost_reducing_weight -= cost_reducing_weights[idx];
    cost_reducing_weights[idx] = w;
    total_cost_reducing_weight += w;
  }
}

inline void remove_cost_reducing_var(const int v_idx) {
  const int idx = cost_reducing_idxs[v_idx];
  if (idx != -1) {
    num_cost_reducing_vars--;
    total_cost_reducing_weight -= cost_reducing_weights[idx];

    // If idx is not at the end, swap the end with this index
    if (idx != num_cost_reducing_vars) {
      const int end_lit = cost_reducing_vars[num_cost_reducing_vars];
      const double w = cost_reducing_weights[num_cost_reducing_vars];
      cost_reducing_vars[idx] = end_lit;
      cost_reducing_weights[idx] = w;
      cost_reducing_idxs[end_lit] = idx;
    }

    cost_reducing_idxs[v_idx] = -1;
  }
}

inline void add_cost_compute_var(const int v_idx) {
  if (cost_compute_idxs[v_idx] == -1) {
    cost_compute_idxs[v_idx] = num_cost_compute_vars;
    cost_compute_vars[num_cost_compute_vars] = v_idx;
    num_cost_compute_vars++;
  }
}

inline void remove_cost_compute_var(const int v_idx) {
  const int idx = cost_compute_idxs[v_idx];
  if (idx != -1) {
    num_cost_compute_vars--;

    // If idx is not at the end, swap the end with this index
    if (idx != num_cost_compute_vars) {
      const int end_var = cost_compute_vars[num_cost_compute_vars];
      cost_compute_vars[idx] = end_var;
      cost_compute_idxs[end_var] = idx;
    }

    cost_compute_idxs[v_idx] = -1;
  }
}


/** @brief Initializes the global formula to get ready for the specified
 *         number of clauses and variables/literals. Allocates the necessary
 *         memory.
 *
 *  @param num_cs The number of clauses in the CNF input file.
 *  @param num_vs The number of variables in the CNF input file.
 */
void initialize_formula(int num_cs, int num_vs) {
  // Set global formula variables
  // Other global variables were given their default values above
  num_vars = num_vs;
  num_literals = 2 * num_vars;
  num_clauses = num_cs;

  const int lit_arr = num_literals + 2;

  log_str("c Allocating memory for %d variables and %d clauses\n",
      num_vars, num_clauses);

  // Allocate all global pointer structures
  assignment = xmalloc((num_vars + 1) * sizeof(char));

  // TODO which need to be calloced?
  clause_sizes = xmalloc(num_clauses * sizeof(int));
  clause_weights = xmalloc(num_clauses * sizeof(double));
  clause_num_true_lits =  xmalloc(num_clauses * sizeof(int));
  clause_lit_masks = xmalloc(num_clauses * sizeof(int));
  clause_literals = xcalloc(num_clauses, sizeof(int *));

  literal_occ = xcalloc(lit_arr, sizeof(int));
  literal_clauses = xcalloc(lit_arr, sizeof(int *));

  false_clause_members = xmalloc(num_clauses * sizeof(int));
  false_clause_indexes = xmalloc(num_clauses * sizeof(int));
  memset(false_clause_indexes, 0xff, num_clauses * sizeof(int));
  
  cost_reducing_vars = xmalloc(num_vars * sizeof(int));
  cost_reducing_idxs = xmalloc((num_vars + 1) * sizeof(int));
  cost_reducing_weights = xmalloc(num_vars * sizeof(double));
  memset(cost_reducing_idxs, 0xff, (num_vars + 1) * sizeof(int));

  cost_compute_vars = xmalloc(num_vars * sizeof(int));
  cost_compute_idxs = xmalloc((num_vars + 1) * sizeof(int));
  memset(cost_compute_idxs, 0xff, (num_vars + 1) * sizeof(int));

  // Allocate memory for neighborhood structures
  allocate_neighborhoods();
  
  log_str("c Allocated memory successfully\n");
}

/** @brief Initializes the clause at the specified index.
 *
 *  TODO
 *
 *  @param clause_idx  The index into the formula->clauses array.
 *  @param size        The number of literals in the clause.
 *  @param lit_idxs    The indexes into formula->literals for the clause's lits.
 */
void initialize_clause(int clause_idx, int size, int *lit_idxs) {
  // Note that clause_num_true_lits and clause_lit_masks will be set
  // when the assignment is randomized, and so do not need to be init. here
  clause_sizes[clause_idx] = size;
  clause_weights[clause_idx] = init_clause_weight;
  clause_literals[clause_idx] = xmalloc(size * sizeof(int));
  memcpy(clause_literals[clause_idx], lit_idxs, size * sizeof(int));
}


/** @brief Processes the clauses read in from the CNF input file and before
 *         the DDFW algorithm is run.
 *
 *  Some post-processing that is done is
 *    - Calculates neighboring clauses for each clause (used in wght dist)
 *    - TODO fill in
 *
 *  @return void, but calls exit(-1) on memory allocation failure.
 */
void process_clauses(void) {
  // Alias variables
  const int nc = num_clauses;
  int *sizes = clause_sizes;
  int **literals = clause_literals;

  int *clause_counts = xcalloc((num_literals + 2), sizeof(int));

  // Loop through the clauses to add clause index to each literal in the clause
  for (int c = 0; c < nc; c++) {
    const int size = *sizes;
    int *lits = *literals;

    // For each literal in the clause, add the clause to the literal
    for (int l = 0; l < size; l++) {
      const int l_idx = *lits;
      if (literal_clauses[l_idx] == NULL) {
        literal_clauses[l_idx] = xmalloc(literal_occ[l_idx] * sizeof(int));
      }

      literal_clauses[l_idx][clause_counts[l_idx]] = c;
      clause_counts[l_idx] = clause_counts[l_idx] + 1;
      lits++;
    }

    sizes++;
    literals++;
  }

  free(clause_counts);

  // Add all literals to cost compute structure
  for (int i = 1; i <= num_vars; i++) {
    add_cost_compute_var(i);
  }

  // After processing everything above, initialize neighborhood structures
  compute_neighborhoods();
}


/** @brief Resets the various data structures in use by clause.c
 *         so the algorithm can run again on the same CNF file.
 */
void reset_data_structures(void) {
  num_flips = 0;
  lowest_unsat_clauses = num_clauses;
  lowest_unsat_step = 0;
  
  num_unsat_clauses = 0;
  memset(false_clause_indexes, 0xff, num_clauses * sizeof(int));

  num_cost_reducing_vars = 0;
  memset(cost_reducing_idxs, 0xff, (num_vars + 1) * sizeof(int));

  num_cost_compute_vars = 0;
  memset(cost_compute_idxs, 0xff, (num_vars + 1) * sizeof(int));

  // Redistribute weights to the clauses
  double *weights = clause_weights;
  for (int c = 0; c < num_clauses; c++) {
    *weights = init_clause_weight;
    weights++;
  }
}

/** @brief Generates a random variable assignment for the global formula.
 *
 *  A random assignment is chosen for the num_vars variables. The random
 *  assignment is generated by storing random numbers into the assignment
 *  bitvector in the global formula. A bit of 1 indicates that the variable
 *  is set to true and its negation is false.
 *
 *  The clauses and literals are updated to reflect the random assignment,
 *  which includes
 *    
 *    - Sat
 *    - Mask
 *    - TODO
 *
 *  Note that the function does not initialize the srand() randomization
 *  library. Such randomization should be handled at command-line parsing.
 *
 *
 *   TODO clean up docs
 *  1. A random assignment is chosen for the num_vars variables.
 *     The random assignment is generated by looping over each
 *     "slot" in the assignment variable and assigning it a random
 *     bit. A 1 in slot i indicates that x_i is assigned to true.
 *
 *  2. The clauses are looped over to determine satisfiability.
 *     After a random assignment is generated, many clauses may evaluate
 *     to true. Looping over the clauses determines how many literals in
 *     each clause are true, and overall how many clauses are true.
 *     Note the usage of an XOR satisfiability mask in each clause to help
 *     determine which literal is the last one keeping a clause true in
 *     the case that clause->sat_lits is 1.
 */
void generate_random_assignment(void) {
  log_str("c Randomizing assignment\n");

  unsat_clause_weight = 0.0;

  // Give random bits to assignment
  // TODO cast to int* for BITS_IN_BYTE fewer calls to rand()
  //const int slots = (num_vars + BITS_IN_BYTE) / BITS_IN_BYTE;
  for (int i = 1; i <= num_vars; i++) {
    assignment[i] = (char) (rand() & 0x1);
  }

  // TODO think about whether clause -> literal is faster than literal -> clause
  // Iterate through the clauses to see which are satisfied
  const int nc = num_clauses;
  int *mask = clause_lit_masks;
  int *sizes = clause_sizes;
  int *true_lits = clause_num_true_lits;
  int **literals = clause_literals;
  int new_mask, sat_lits;

  for (int i = 0; i < nc; i++) {
    sat_lits = 0;
    new_mask = 0;
    int *ls = *literals;
    const int size = *sizes;

    // Loop through the literals to see how many are satisfied
    for (int l = 0; l < size; l++) {
      int lit_idx = *ls;
      int negated = IS_NEGATED(lit_idx);
      int assigned = ASSIGNMENT(lit_idx);

      // If the literal evaluates to true, update clause information
      if ((assigned && !negated) || (!assigned && negated)) {
        sat_lits++;
        new_mask ^= lit_idx;
      }

      ls++;
    }

    // If not satisfied, update membership in false clauses
    if (sat_lits == 0) {
      add_false_clause(i);
    }

    // Store lit and mask information back to the arrays
    *true_lits = sat_lits;
    *mask = new_mask;

    // Move pointers to next clause
    sizes++;
    mask++;
    true_lits++;
    literals++;
  }

  lowest_unsat_clauses = num_unsat_clauses;

  // Once all sat/unsat clauses have been computed, update neighborhood structs
  initialize_neighborhoods();

  if (get_verbosity() == VERBOSE) {
    log_assignment();
  }
}

/** @brief Takes an index of a literal to flip and flips it in the assignment.
 *
 *  The function will change the truth value of the variable corresponding to
 *  the provided lit_idx.
 *
 *  Any internal changes as a result of flipping the variable (weights, marking
 *  literals, etc.) also occur here.
 *
 *  @param lit_idx The index into the literals array of the literal to flip.
 */
void flip_variable(const int var_idx) {
  const int pos_lit_idx = LIT_IDX(var_idx);
  const int not_lit_idx = NEGATED_IDX(pos_lit_idx);
  const int assigned = ASSIGNMENT(pos_lit_idx);

  add_cost_compute_var(var_idx);

  // Determine which literal has the truth value, to flip correct clauses
  int l_idx, not_l_idx;
  if (assigned) {
    l_idx = pos_lit_idx;
    not_l_idx = not_lit_idx;
  } else {
    l_idx = not_lit_idx;
    not_l_idx = pos_lit_idx;
  }

  const int l_occ = literal_occ[l_idx];
  const int not_occ = literal_occ[not_l_idx];

  // For each clause containing l, set l to false
  int *l_to_clauses = literal_clauses[l_idx];
  for (int c = 0; c < l_occ; c++) {
    const int c_idx = *l_to_clauses;

    // Remove literal from sat XOR mask
    clause_num_true_lits[c_idx]--;
    const int true_lits = clause_num_true_lits[c_idx];
    clause_lit_masks[c_idx] ^= pos_lit_idx;

    // In the case where the clause now has 0 true literals,
    //   all literals become critical, and may be added to cost compute
    if (true_lits == 0) {
      // Also, the clause is false, and may be added to that list
      add_false_clause(c_idx);

      const int size = clause_sizes[c_idx];
      int *c_to_lits = clause_literals[c_idx];
      for (int cl_lit = 0; cl_lit < size; cl_lit++) {
        const int cl_var_idx = VAR_IDX(*c_to_lits);
        add_cost_compute_var(cl_var_idx); 
        c_to_lits++;
      }

      // The clause affects the neighborhood on flipping sign
      update_neighborhood_on_flip(c_idx);
    } else if (true_lits == 1) {
      // If instead the clause has 1 true literal, only the last true
      // literal is critical and can make the clause false if flipped
      // Get the last true literal from the literal mask
      const int mask_lit = clause_lit_masks[c_idx];
      add_cost_compute_var(VAR_IDX(mask_lit));
    }

    l_to_clauses++;
  }

  // For each clause containing !l, set !l to true
  int *not_l_to_clauses = literal_clauses[not_l_idx];
  for (int c = 0; c < not_occ; c++) {
    const int c_idx = *not_l_to_clauses;

    // Add literal to sat XOR mask
    clause_num_true_lits[c_idx]++;
    const int true_lits = clause_num_true_lits[c_idx];
    const int mask_before = clause_lit_masks[c_idx]; // See else if case
    clause_lit_masks[c_idx] ^= not_l_idx;

    // If the clause was just made true, then all literals may now reduce
    //   cost less, and so must be re-computed
    if (true_lits == 1) {
      // Also, the clause is no longer false, and may be removed
      remove_false_clause(c_idx);

      const int size = clause_sizes[c_idx];
      int *c_to_lits = clause_literals[c_idx];
      for (int cl_lit = 0; cl_lit < size; cl_lit++) {
        const int cl_var_idx = VAR_IDX(*c_to_lits);
        add_cost_compute_var(cl_var_idx);
        c_to_lits++;
      }

      // The clause affects the neighborhood on flipping sign
      update_neighborhood_on_flip(c_idx);
    } else if (true_lits == 2) {
      // If instead the literal is made "non-critical," only the other true
      //   literal must be re-computed
      add_cost_compute_var(VAR_IDX(mask_before));
    }

    not_l_to_clauses++;
  }

  // Affect assignment bitvector
  assignment[var_idx] = !assignment[var_idx];

  if (get_verbosity() == VERBOSE) {
    log_assignment();
  }

  if (num_unsat_clauses < lowest_unsat_clauses) {
    log_str("c Record %d after %d flips\n", 
        num_unsat_clauses, num_flips);
    lowest_unsat_clauses = num_unsat_clauses;
    lowest_unsat_step = num_flips;
  }

  num_flips++;
}
