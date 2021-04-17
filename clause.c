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
#include "verifier.h"

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
char *best_assignment = NULL;
long num_flips_since_improvement = 0;

// Clause information - 0 indexed
int *clause_sizes = NULL;
double *clause_weights = NULL;
int *clause_num_true_lits = NULL;
int *clause_lit_masks = NULL;
int **clause_literals = NULL;

// Literal information - 2-indexed (use LIT_IDX)
int *literal_occ = NULL;
int **literal_clauses = NULL;

// Breakdown of literal clause weight info
double *literal_unsat_weights = NULL;
double *literal_crit_sat_weights = NULL;

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

// TODO experimental
long ddfw_plus_counter = 0;
int ddfw_plus_boolean = 0;
int ddfw_reweighting_enabled = 0;

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
  best_assignment = xmalloc((num_vars + 1) * sizeof(char));

  // TODO which need to be calloced?
  clause_sizes = xmalloc(num_clauses * sizeof(int));
  clause_weights = xmalloc(num_clauses * sizeof(double));
  clause_num_true_lits =  xmalloc(num_clauses * sizeof(int));
  clause_lit_masks = xmalloc(num_clauses * sizeof(int));
  clause_literals = xcalloc(num_clauses, sizeof(int *));

  literal_occ = xcalloc(lit_arr, sizeof(int));
  literal_clauses = xcalloc(lit_arr, sizeof(int *));

  literal_unsat_weights = xcalloc(lit_arr, sizeof(double));
  literal_crit_sat_weights = xcalloc(lit_arr, sizeof(double));

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
  num_flips_since_improvement = 0;
  lowest_unsat_clauses = num_clauses;
  lowest_unsat_step = 0;
  
  num_unsat_clauses = 0;
  memset(false_clause_indexes, 0xff, num_clauses * sizeof(int));

  total_cost_reducing_weight = 0.0;
  num_cost_reducing_vars = 0;
  memset(cost_reducing_idxs, 0xff, (num_vars + 1) * sizeof(int));
  memset(cost_reducing_weights, 0, num_vars * sizeof(double));

  num_cost_compute_vars = 0;
  memset(cost_compute_idxs, 0xff, (num_vars + 1) * sizeof(int));

  // Add all literals to cost compute structure
  for (int i = 1; i <= num_vars; i++) {
    add_cost_compute_var(i);
  }

  // Redistribute weights to the clauses
  double *weights = clause_weights;
  for (int c = 0; c < num_clauses; c++) {
    *weights = init_clause_weight;
    weights++;
  }

  memset(literal_crit_sat_weights, 0, (num_literals + 2) * sizeof(double));
  memset(literal_unsat_weights, 0, (num_literals + 2) * sizeof(double));
}

/** @brief Resets the weights held by all clauses. Sets equal to W_init. */
void reset_clause_weights(void) {
  double *weights = clause_weights;
  for (int c = 0; c < num_clauses; c++) {
    *weights = init_clause_weight;
    weights++;
  }
}

/** @brief Resets the structure that records which variables, when flipped,
 *         reduce the amount of weight held by unsatisfied clauses.
 */
void reset_cost_reducing_struct(void) {
  num_cost_reducing_vars = 0;
  total_cost_reducing_weight = 0.0;
  memset(cost_reducing_idxs, 0xff, (num_vars + 1) * sizeof(int));
  memset(cost_reducing_weights, 0, num_vars * sizeof(double));
  // This isn't strictly necessary, hence why it's commented out
  // memset(cost_reducing_vars, 0, num_vars * sizeof(int));
}

/** @brief Resets the structure that records which variables need their
 *         cost reduction values computed.
 *
 *  Adds all variables to the structure to be computed.
 */
void reset_cost_compute_struct(void) {
  num_cost_compute_vars = 0;
  memset(cost_compute_idxs, 0xff, (num_vars + 1) * sizeof(int));
  // This isn't strictly necessary, hence why it's commented out
  // memset(cost_compute_vars, 0, num_vars * sizeof(int));
  for (int i = 1; i <= num_vars; i++) {
    add_cost_compute_var(i);
  }
}

// TODO experimental
static void compute_critical_weights_struct(void) {
  // Compute unsat weight and critical sat weight for all literals
  // TODO alias pointers later
  for (int c = 0; c < num_clauses; c++) {
    const int clause_size = clause_sizes[c];
    int *literals = clause_literals[c];
    const double w = clause_weights[c];
    if (clause_num_true_lits[c] == 0) {
      for (int l = 0; l < clause_size; l++) {
        literal_unsat_weights[literals[l]] += w;
      }
    } else if (clause_num_true_lits[c] == 1) {
      const int mask_lit = clause_lit_masks[c];
      literal_crit_sat_weights[mask_lit] += w;
    }
  }
}

/** @brief Resets the structure that keeps track of the unsat and critical
 *         weights associated with each literal.
 */
void reset_critical_weights_struct(void) {
  memset(literal_crit_sat_weights, 0, (num_literals + 2) * sizeof(double));
  memset(literal_unsat_weights, 0, (num_literals + 2) * sizeof(double));
}

/** @brief Resets the structure that keeps track of the false clauses. */
static void reset_false_clause_struct(void) {
  num_unsat_clauses = 0;
  memset(false_clause_indexes, 0xff, num_clauses * sizeof(int));
  // memset(false_clause_members, 0, num_clauses * sizeof(int));
}

// TODO experimental
void reset_to_ddfw_plus_weightings(void) {
  ddfw_plus_counter++;
  if (ddfw_plus_counter >= num_clauses) {
    // printf("Resetting to ddfw plus at %ld flips\n", num_flips);
    ddfw_plus_counter = 0;
    if (!ddfw_plus_boolean) {
      ddfw_plus_boolean = 1;
      const double w_to_add = init_clause_weight / 2.0;
      for (int c = 0; c < num_clauses; c++) {
        clause_weights[c] += w_to_add;
      }
    } else {
      ddfw_plus_boolean = 0;
      const double sw = init_clause_weight * 1.5;
      for (int c = 0; c < num_clauses; c++) {
        if (clause_num_true_lits[c] == 0) {
          clause_weights[c] = init_clause_weight;
        } else {
          clause_weights[c] = sw;
        }
      }
    }

    reset_cost_reducing_struct();
    reset_cost_compute_struct();
    reset_critical_weights_struct();
    compute_critical_weights_struct();
    initialize_neighborhoods();
  }
}

/** @brief Resets all structures that the clause file is in charge of
 *         to prepare for a new run of the DDFW algorithm.
 *
 *  The clause file is in charge of:
 *    - clause weights
 *    - cost reducing variables
 *    - cost compute variables
 *    - (best) assignment
 *    - Various run statistics
 */
void reset_clause_for_alg_run(void) {
  reset_clause_weights();
  reset_cost_reducing_struct();
  reset_cost_compute_struct();
  reset_critical_weights_struct();
  reset_false_clause_struct();
  
  // Reset various other statistics
  unsat_clause_weight = 0.0;
  num_flips = 0;
  lowest_unsat_clauses = num_clauses;
  lowest_unsat_step = 0;

  // TODO experimental
  ddfw_plus_counter = 0;
  ddfw_plus_boolean = 0;
}

/** @brief Resets the structures the clause file is in charge of when the
 *         assignment changes by more than one flip at a time.
 *
 *  Called in generate_new_assignment() and when restoring to best assign.
 */
static void reset_clause_structures_for_new_assignment(void) {
  reset_false_clause_struct();

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
      int assigned = ASSIGNMENT(lit_idx);
      if (assigned) {
        sat_lits++;
        new_mask ^= lit_idx;
      }

      ls++;
    }

    if (sat_lits == 0) {
      add_false_clause(i);
    }

    // Store lit and mask information back to the arrays
    *true_lits = sat_lits;
    *mask = new_mask;

    sizes++;
    mask++;
    true_lits++;
    literals++;
  }

  lowest_unsat_clauses = num_unsat_clauses;

  // Once all sat/unsat clauses have been computed, update neighborhood structs
  initialize_neighborhoods();
  compute_critical_weights_struct();
}

/** @brief Restores the structures that the clause file is in charge of
 *         back to the best assignment found so far.
 */
void restore_to_best_assignment(void) {
  reset_clause_weights();
  reset_cost_reducing_struct();
  reset_cost_compute_struct();
  memcpy(assignment, best_assignment, (num_vars + 1) * sizeof(char));
  num_flips_since_improvement = 0;
  reset_clause_structures_for_new_assignment();
  verify_clauses_and_assignment();
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

  // Give random bits to assignment
  for (int i = 1; i <= num_vars; i++) {
    assignment[i] = (char) (rand() & 0x1);
  }

  reset_clause_structures_for_new_assignment(); 
  verify_clauses_and_assignment();
}

/** @brief Takes an index of a variable to flip and flips it in the assignment.
 *
 *  The function will change the truth value of the variable and update the
 *  appropriate structures.
 *
 *  @param lit_idx The index into the literals array of the literal to flip.
 */
void flip_variable(const int var_idx) {
  const int pos_lit_idx = LIT_IDX(var_idx);
  const int not_lit_idx = NEGATED_IDX(pos_lit_idx);
  const int assigned = ASSIGNMENT(pos_lit_idx);

  // Affect assignment bitvector
  assignment[var_idx] = !assignment[var_idx];
  add_cost_compute_var(var_idx);

  // Determine which literal has the truth value, to flip correct clauses
  int l_idx = (assigned) ? pos_lit_idx : not_lit_idx;
  int not_l_idx = (assigned) ? not_lit_idx : pos_lit_idx;

  const int l_occ = literal_occ[l_idx];
  const int not_occ = literal_occ[not_l_idx];

  // Since l is being set to false, then by definition, it cannot be a sat lit
  literal_crit_sat_weights[l_idx] = 0;

  // For each clause containing l, set l to false
  int *l_to_clauses = literal_clauses[l_idx];
  for (int c = 0; c < l_occ; c++) {
    const int c_idx = *l_to_clauses;

    // Remove literal from sat XOR mask
    clause_num_true_lits[c_idx]--;
    const int true_lits = clause_num_true_lits[c_idx];
    clause_lit_masks[c_idx] ^= l_idx;

    // In the case where the clause now has 0 true literals,
    //   all literals become critical, and may be added to cost compute
    if (true_lits == 0) {
      // Also, the clause is false, and may be added to that list
      add_false_clause(c_idx);

      // TODO document later
      const int size = clause_sizes[c_idx];
      int *c_to_lits = clause_literals[c_idx];
      const double w = clause_weights[c_idx];
      for (int cl_lit = 0; cl_lit < size; cl_lit++) {
        const int lit = *c_to_lits;
        const int cl_var_idx = VAR_IDX(lit);
        add_cost_compute_var(cl_var_idx); 

        // Take weight to unsat
        literal_unsat_weights[lit] += w;

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

      // Add weight from all literals to crit sat weights
      const double w = clause_weights[c_idx];
      literal_crit_sat_weights[mask_lit] += w;
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
      const double w = clause_weights[c_idx];
      literal_crit_sat_weights[not_l_idx] += w;
      for (int cl_lit = 0; cl_lit < size; cl_lit++) {
        const int lit = *c_to_lits;
        const int cl_var_idx = VAR_IDX(lit);
        add_cost_compute_var(cl_var_idx);

        // Move weight from unsat to crit sat
        literal_unsat_weights[lit] -= w;

        c_to_lits++;
      }

      // The clause affects the neighborhood on flipping sign
      update_neighborhood_on_flip(c_idx);
    } else if (true_lits == 2) {
      // If instead the literal is made "non-critical," only the other true
      //   literal must be re-computed
      add_cost_compute_var(VAR_IDX(mask_before));

      // Now that clause is no longer critical, remove crit sat weight
      const double w = clause_weights[c_idx];
      literal_crit_sat_weights[mask_before] -= w;
    }

    not_l_to_clauses++;
  }
 
  if (get_verbosity() == VERBOSE) {
    log_assignment();
  }

  if (num_unsat_clauses < lowest_unsat_clauses) {
    lowest_unsat_clauses = num_unsat_clauses;
    lowest_unsat_step = num_flips;
    memcpy(best_assignment, assignment, (num_vars + 1) * sizeof(char));
    num_flips_since_improvement = 0;
  } else {
    ddfw_plus_counter++;
    num_flips_since_improvement++;
    if (ddfw_reweighting_enabled) {
      reset_to_ddfw_plus_weightings();
    }
  }

  num_flips++;
}
