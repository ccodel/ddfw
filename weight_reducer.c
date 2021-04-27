/** @file weight_reducer.c
 *  @brief Computes the Boolean variables that, when flipped, reduce the clause
 *         weight held by the unsatisfied clauses.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#include <string.h>
#include <stdio.h>

#include "weight_reducer.h"
#include "global_data.h"
#include "xmalloc.h"

/** @brief Stores the variables which, when flipped, reduce the weight held
 *         by the unsatisfied clauses.
 *
 *  The array is 0-indexed, and its values are literal indexes (see LIT_IDX).
 */
int *unsat_weight_reducing_vars;

/** @brief Stores the indexes of the literal indexes in the _vars array.
 *
 *  The array is 1-indexed by literal index (see LIT_IDX).
 */
static int *unsat_weight_reducing_idxs;

/** @brief The number of variables which, when flipped, reduce the weight held
 *         by the unsatisfied clauses.
 */
int num_unsat_weight_reducing_vars;

/** Indexed like in _vars, not in _idxs */
weight *unsat_weight_reducing_weights;
weight total_unsat_weight_reducing_weight;

static int *unsat_weight_compute_vars;
static int *unsat_weight_compute_idxs;
static int num_unsat_weight_compute_vars;

weight *literal_unsat_weights;
weight *literal_critical_sat_weights;

void allocate_weight_reducer_memory(void) {
  unsat_weight_reducing_vars = xmalloc(num_vars * sizeof(int));
  unsat_weight_reducing_weights = xmalloc(num_vars * sizeof(weight));
  unsat_weight_reducing_idxs = xmalloc((num_vars + 1) * sizeof(int));
  unsat_weight_compute_vars = xmalloc(num_vars * sizeof(int));
  unsat_weight_compute_idxs = xmalloc((num_vars + 1) * sizeof(int));
  literal_unsat_weights = xmalloc((num_literals + 2) * sizeof(weight));
  literal_critical_sat_weights = xmalloc((num_literals + 2) * sizeof(weight));
}

static void add_weight_reducing_var(const int v_idx, const weight w) {
  const int idx = unsat_weight_reducing_idxs[v_idx];
  // If not present, add to structure. If present, update weight
  if (idx == -1) {
    // Add to end - acts like a queue structure
    unsat_weight_reducing_idxs[v_idx] = num_unsat_weight_reducing_vars;
    unsat_weight_reducing_vars[num_unsat_weight_reducing_vars] = v_idx;
    unsat_weight_reducing_weights[num_unsat_weight_reducing_vars] = w;
    total_unsat_weight_reducing_weight += w;
    num_unsat_weight_reducing_vars++;
  } else {
    total_unsat_weight_reducing_weight -= unsat_weight_reducing_weights[idx];
    unsat_weight_reducing_weights[idx] = w;
    total_unsat_weight_reducing_weight += w;
  }
}

static void remove_weight_reducing_var(const int v_idx) {
  const int idx = unsat_weight_reducing_idxs[v_idx];
  // If the variable is weight reducing, remove it and swap with end variable
  if (idx != -1) {
    num_unsat_weight_reducing_vars--;
    total_unsat_weight_reducing_weight -= unsat_weight_reducing_weights[idx];

    // If idx is not at the end, swap with ending
    if (idx != num_unsat_weight_reducing_vars) {
      int end_var = unsat_weight_reducing_vars[num_unsat_weight_reducing_vars];
      weight w = unsat_weight_reducing_weights[num_unsat_weight_reducing_vars];
      unsat_weight_reducing_vars[idx] = end_var;
      unsat_weight_reducing_weights[idx] = w;
      unsat_weight_reducing_idxs[end_var] = idx;
    }

    unsat_weight_reducing_idxs[v_idx] = -1;
  }
}

void add_weight_compute_var(const int v_idx) {
  // If the variable is not queued, add to end
  if (unsat_weight_compute_idxs[v_idx] == -1) {
    unsat_weight_compute_idxs[v_idx] = num_unsat_weight_compute_vars;
    unsat_weight_compute_vars[num_unsat_weight_compute_vars] = v_idx;
    num_unsat_weight_compute_vars++;
  }
}

static void clear_weight_reducing_structures(void) {
  memset(unsat_weight_reducing_idxs, 0xff, (num_vars + 1) * sizeof(int));
  memset(unsat_weight_compute_idxs, 0xff, (num_vars + 1) * sizeof(int));
  num_unsat_weight_reducing_vars = 0;
  total_unsat_weight_reducing_weight = 0.0;
  num_unsat_weight_compute_vars = 0;
  for (int v = 1; v <= num_vars; v++) {
    add_weight_compute_var(v);
  }
}

static void compute_literal_critical_weights(void) {
  memset(literal_unsat_weights, 0, (num_literals + 2) * sizeof(weight));
  memset(literal_critical_sat_weights, 0, (num_literals + 2) * sizeof(weight));

  for (int i = 0; i < num_clauses; i++) {
    const int clause_size = clause_sizes[i];
    const weight w = clause_weights[i];
    const int num_true_lits = clause_num_true_lits[i];
    int *literals = clause_literals[i];

    // If the clause has no literals which evaluate to true, then the clause
    // contributes to the unsat weight, and flipping any literal inside will
    // make the clause true, and thus the weight "becomes satisfied"
    if (num_true_lits == 0) {
      for (int l = 0; l < clause_size; l++) {
        literal_unsat_weights[literals[l]] += w;
      }
    } else if (num_true_lits == 1) {
      // Flipping the true variable to false would decrease the satisfied weight
      literal_critical_sat_weights[clause_lit_masks[i]] += w;
    }
  }
}

/** @brief Computes the variables which, when flipped, reduces the weight held
 *         by the unsatisfied clauses.
 *
 *  The results are stored in the unsat_weight_reducing_vars array.
 *
 *  Typically, when a single clause's weight is changed or when a variable is
 *  flipped, only a small number of variables must have their weight reducing
 *  status re-computed. This function takes those variables from the 
 *  unsat_weight_compute_vars array and processes them each in turn. The array
 *  is then cleared at the end.
 */
void compute_weight_reducing_variables(void) {
  for (int i = 0; i < num_unsat_weight_compute_vars; i++) {
    const int v_idx = unsat_weight_compute_vars[i];
    const int l_idx = LIT_IDX(v_idx);
    const int assigned = ASSIGNMENT(l_idx);
    const int true_idx = (assigned) ? l_idx : NEGATED_IDX(l_idx);
    const int false_idx = NEGATED_IDX(true_idx);

    // Compare the unsat weight gained if the true literal were flipped
    //  versus the unsat weight lost if the false literal were flipped
    const weight diff = literal_unsat_weights[false_idx] - 
                        literal_critical_sat_weights[true_idx];
    if (diff > 0) {
      add_weight_reducing_var(v_idx, diff);
    } else {
      remove_weight_reducing_var(v_idx);
    }

    // Remove the variable from the weight compute array
    unsat_weight_compute_idxs[v_idx] = -1;
  }

  num_unsat_weight_compute_vars = 0;

  // This does nothing if DEBUG is not on
  // TODO give this responsibility to the verifier.c
  verify_weight_reducer();
}

void compute_weight_reducing_after_assignment(void) {
  clear_weight_reducing_structures();
  compute_literal_critical_weights();
  compute_weight_reducing_variables();
}

void compute_weight_reducing_after_reweighting(void) {
  clear_weight_reducing_structures();
  compute_literal_critical_weights();
  compute_weight_reducing_variables();
}

/** @brief Verifies the internal consistency of the weight reducer structs.
 *
 *  First calculates the critical weight values, then verifies that the current
 *  state of weight reducing variables agrees with hard-computed values.
 */
#ifdef DEBUG
void verify_weight_reducer(void) {
  // First, we check that unsat/critical weights are correct, by literal
  for (int l = 2; l <= num_literals + 2; l++) {
    weight sat_w = 0;
    weight unsat_w = 0;
    const int is_true = ASSIGNMENT(l);

    // Loop through the clauses for this literal, add appropriate weights
    const int occ = literal_occ[l];
    int *clauses = literal_clauses[l];
    for (int c = 0; c < occ; c++) {
      const int c_idx = clauses[c];
      if (is_true && clause_num_true_lits[c_idx] == 1) {
        sat_w += clause_weights[c_idx];
      } else if (clause_num_true_lits[c_idx] == 0) {
        unsat_w += clause_weights[c_idx];
      }
    }

    if (is_true) {
      ERR_IF(ABS(literal_critical_sat_weights[l] - sat_w) > 0.1, "Crit sat\n")
    }

    ERR_IF(ABS(literal_unsat_weights[l] - unsat_w) > 0.1, "Unsat w\n")
  }

  // Next, compute weight reducing all by hand
  // NOTE: We can use the critical weight structs since we just checked those
  weight total_wrw = 0;
  int num_weight_reducing = 0;

  for (int v = 0; v <= num_vars; v++) {
    const int l_idx = LIT_IDX(v);
    const int assigned = ASSIGNMENT(l_idx);
    const int true_idx = (assigned) ? l_idx : NEGATED_IDX(l_idx);
    const int false_idx = NEGATED_IDX(true_idx);
    const int true_occ = literal_occ[true_idx];
    const int false_occ = literal_occ[false_idx];

    const weight diff = literal_unsat_weights[false_idx] - 
                        literal_critical_sat_weights[true_idx];
    const weight curr_diff = unsat_weight_reducing_weights[v];
    if (diff > 0) {
      // Check to see if in the list
      ERR_IF(unsat_weight_reducing_idxs[v] == -1, "Weight reducing idxs\n")
      ERR_IF(ABS(diff - curr_diff) > 1, "Weight reducing diffs\n")
      num_weight_reducing++;
      total_wrw += diff;
    } else {
      ERR_IF(unsat_weight_reducing_idxs[v] != -1, "!Weight reducing idxs\n")
      ERR_IF(ABS(diff - curr_diff) > 1, "!Weight reducing diff\n")
    }
  }

  ERR_IF(num_weight_reducing != num_weight_reducing_vars, "Num wr vars\n")
  ERR_IF(ABS(total_weight_reducing_weight - total_wrw) > 1, "Total wrw\n")
}

#else

void verify_weight_reducer(void) { (void) 0; }

#endif /* DEBUG */
