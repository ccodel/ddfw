/** @file weight_reducer.h
 *  @brief Computes the Boolean variables that, when flipped, reduce the clause
 *         weight held by the unsatisfied clauses.
 *
 *  See weight_reducer.c for implementation details.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#ifndef _WEIGHT_REDUCER_H_
#define _WEIGHT_REDUCER_H_

#include "global_data.h"

extern int *unsat_weight_reducing_vars;
extern int num_unsat_weight_reducing_vars;

extern weight *unsat_weight_reducing_weights;
extern weight total_unsat_weight_reducing_weight;

extern weight *literal_unsat_weights;
extern weight *literal_critical_sat_weights;

/****************************** FUNCTIONS *************************************/

void allocate_weight_reducer_memory(void);

// To incrementally add variables when weight/assignment changed
void add_weight_compute_var(const int v_idx);

// To be called in normal operation
void compute_weight_reducing_variables(void);

// To be called after a larger change in algorithmic state
void compute_weight_reducing_after_assignment(void);
void compute_weight_reducing_after_reweighting(void);

void verify_weight_reducer(void);

#endif /* _WEIGHT_REDUCER_H_ */
