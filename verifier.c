/** @file verifier.c
 *  @brief A collection of verification methods that, unless specifically
 *         compiled in, are not included in the compilation of DDFW.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#include "xmalloc.h"
#include "verifier.h"
#include "clause.h"
#include "neighborhood.h"
#include "weight_transfer.h"

/** If debug is enabled, implement the verifier functions */
#ifdef DEBUG

#ifndef ABS
#define ABS(x)     (((x) < 0) ? -(x) : (x))
#endif

#ifndef MIN
#define MIN(x, y)  (((x) < (y)) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x, y)  (((x) > (y)) ? (x) : (y))
#endif

#define ERR_IF(cond, msg)  if (cond) { fprintf(stderr, msg); exit(-1); }

// Declare some variables used by the verifier

/** Copy of membership array data structure used for cost reducing calcs. */
static int *cost_reducing_idxs_copy = NULL;
static int *cost_reducing_vars_copy = NULL;
static double *cost_reducing_weights_copy = NULL;

void initialize_verifier(void) {
  if (cost_reducing_idxs_copy == NULL) {
    cost_reducing_idxs_copy = xmalloc((num_vars + 1) * sizeof(int));
    cost_reducing_vars_copy = xmalloc(num_vars * sizeof(int));
    cost_reducing_weights_copy = xmalloc(num_vars * sizeof(double));
  }
}

/** @brief Verifies the consistency of maximum and weights for neighborhoods. */
void verify_neighborhoods(void) {
  for (int c = 0; c < num_clauses; c++) {
    int num_sat = 0;
    double w = 0.0;
    double max_w = 0.0;

    // Traverse the clauses in the neighborhood
    int *cs = neigh_clauses[c];
    const int size = neigh_sizes[c];
    for (int i = 0; i < size; i++) {
      const int c_idx = *cs;
      if (clause_num_true_lits[c_idx] > 0) {
        num_sat++;
        w += clause_weights[c_idx];
        if (clause_weights[c_idx] > max_w) {
          max_w = clause_weights[c_idx];
        }
      }

      cs++;
    }

    // Once done, check the tracked stats
    ERR_IF(neigh_num_sat[c] != num_sat, "Neigh num sat\n")
    ERR_IF(ABS(neigh_weights[c] - w) > 0.1, "Neigh weights\n")
    if (neigh_max_idxs[c] != -1) {
      ERR_IF(ABS(neigh_max_weights[c] - max_w) > 0.1, "Neigh max weights\n")
      ERR_IF(clause_num_true_lits[neigh_max_idxs[c]] == 0, "Neigh max idx\n")
      ERR_IF(ABS(clause_weights[neigh_max_idxs[c]] - max_w) > 0.1, 
          "Verification of neigh max weight\n")
    }
  }
}

/** @brief Verifies the weights stored for critically satisfied clauses. */
void verify_crit_sat_unsat_weights(void) {
  for (int l = 2; l < num_literals + 2; l++) {
    double sat_w = 0.0;
    double unsat_w = 0.0;
    const int is_true = ASSIGNMENT(l);
   
    // Loop through clauses, add weights of those with one sat lit
    int *cs = literal_clauses[l];
    const int occ = literal_occ[l];
    for (int c = 0; c < occ; c++) {
      const int c_idx = *cs;
      if (is_true && clause_num_true_lits[c_idx] == 1) { 
        ERR_IF(clause_lit_masks[c_idx] != l, "Masks\n")
        sat_w += clause_weights[c_idx];
      } else if (clause_num_true_lits[c_idx] == 0) {
        unsat_w += clause_weights[c_idx];
      }

      cs++;
    }

    // Now compare the computed value with the stored value
    if (is_true) {
      ERR_IF(ABS(literal_crit_sat_weights[l] - sat_w) > 0.1, "Crit sat w\n")
    }
    
    ERR_IF(ABS(literal_unsat_weights[l] - unsat_w) > 0.1, "Unsat w\n")
  }
}

/** @brief Verifies counts of satisfying literals and clauses.
 *
 *  Loops through all clauses to check that the number of satisfying
 *  literals for that clause is correct. Also keeps track of the number
 *  of unsatisfied clauses in the formula and compares against the
 *  stored value "num_unsat_clauses".
 */
void verify_clauses_and_assignment(void) {
  // Alias pointers used for iteration
  int *num_true_lits = clause_num_true_lits;
  int *sizes = clause_sizes;
  int **cl_to_lits = clause_literals;
  int unsat_clause_counter = 0;

  for (int i = 0; i < num_clauses; i++) {
    const int size = *sizes;
    const int true_lits = *num_true_lits;
    int *lits = *cl_to_lits;
    int true_lit = -1;
    int true_mask = 0;

    // Loop through the literals and check against assignment
    int true_lit_counter = 0;
    for (int j = 0; j < size; j++) {
      int assigned = ASSIGNMENT(*lits);
      if (assigned) {
        true_lit = *lits;
        true_lit_counter++;
        true_mask ^= *lits;
      }

      lits++;
    }

    ERR_IF(true_lit_counter != true_lits, "True literal count\n")

    if (true_lit_counter == 0) {
      unsat_clause_counter++;
    } else {
      ERR_IF(clause_lit_masks[i] != true_mask, "Mask lit\n")
      if (true_lit_counter == 1) {
        ERR_IF(clause_lit_masks[i] != true_lit, "Single mask lit\n")
      }
    }

    num_true_lits++;
    sizes++;
    cl_to_lits++;
  }

  ERR_IF(unsat_clause_counter != num_unsat_clauses, "Unsat clause number\n")
}


/** @brief Verifies the reported cost reducing vars by looping over all
 *         variables, rather than taking the "cost_compute_vars" array
 *         at its word.
 */
void verify_cost_reducing_vars(void) {
  // First step, copy over the cost reducing indexes array to compare later
  memcpy(cost_reducing_idxs_copy, cost_reducing_idxs,
      (num_vars + 1) * sizeof(int));
  memset(cost_reducing_idxs, 0xff, (num_vars + 1) * sizeof(int));

  // Store previous indexes, as they get messed up below
  memcpy(cost_reducing_vars_copy, cost_reducing_vars,
      num_cost_reducing_vars * sizeof(int));

  // Next, copy over the weights
  memcpy(cost_reducing_weights_copy, cost_reducing_weights,
      num_cost_reducing_vars * sizeof(double));
  memset(cost_reducing_weights, 0, num_cost_reducing_vars * sizeof(int));

  // Alias arrays
  double *crw = cost_reducing_weights;
  double *crwc = cost_reducing_weights_copy;

  // Next, clear number of cost reducing vars and loop through all vars
  // TODO this is just copied from below, pull out into helper function?
  const int prev_num_cost_reducing = num_cost_reducing_vars;
  const double prev_total_cost = total_cost_reducing_weight;
  num_cost_reducing_vars = 0;

  // Due to floating point error, we keep an error bound for "close calls"
  int error_bound = 0;

  total_cost_reducing_weight = 0.0;
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

    if (ABS(diff) < 0.1) {
      error_bound++;
    }

    if (diff > 0.0) {
      add_cost_reducing_var(v_idx, diff);
    } else if (diff <= 0.0) {
      remove_cost_reducing_var(v_idx);
    }
  }

  ERR_IF((prev_num_cost_reducing > num_cost_reducing_vars + error_bound) ||
         (prev_num_cost_reducing < num_cost_reducing_vars - error_bound), 
         "num crvs\n")
  ERR_IF(ABS(total_cost_reducing_weight - prev_total_cost) > 0.01, "crw diff\n")

  // For each variable, check the internal consistency of both arrays
  // Then cross-check arrays for membership and stored weight
  for (int i = 1; i <= num_vars; i++) {
    const int crix = cost_reducing_idxs[i];
    const int cricx = cost_reducing_idxs_copy[i];

    // If both are present, check that their arrays are interally consistent
    if (crix != -1 && cricx != -1) {
      ERR_IF(cost_reducing_vars_copy[cricx] != i, "Crvs not right\n")
      ERR_IF(cost_reducing_vars[crix] != i, "Crvs copy not right\n")
      ERR_IF(ABS(crw[crix] - crwc[cricx]) > 0.01, "Crws differed\n")
    }
  }

  // Copy the copies back
  memcpy(cost_reducing_idxs, cost_reducing_idxs_copy,
      (num_vars + 1) * sizeof(int));
  memcpy(cost_reducing_weights, cost_reducing_weights_copy,
      num_cost_reducing_vars * sizeof(double));
  memcpy(cost_reducing_vars, cost_reducing_vars_copy,
      num_cost_reducing_vars * sizeof(int));
  total_cost_reducing_weight = prev_total_cost;
  num_cost_reducing_vars = prev_num_cost_reducing;
}

#else

void initialize_verifier(void) { (void) 0; }
void verify_neighborhoods(void) { (void) 0; }
void verify_crit_sat_unsat_weights(void) { (void) 0; }
void verify_clauses_and_assignment(void) { (void) 0; }
void verify_cost_reducing_vars(void) { (void) 0; }

#endif /* DEBUG */
