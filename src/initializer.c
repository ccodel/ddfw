/** @file initializer.c
 *  @brief Handles initializing and resetting various structures.
 *
 *  The normal operation of a stochastic local search algorithm proceeds in
 *  steps, e.g. flipping the sign of a singular Boolean variable. However,
 *  if the algorithm decides to change the assignment or the clause weights
 *  more than one "atom" at a time, then structures must be recomputed.
 *  The initializer file handles the recomputing needed across all files.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#include "initializer.h"
#include "weight_reducer.h"
#include "assignment.h"
#include "neighborhood.h"
#include "verifier.h"

#include <stdio.h>

/** @brief Compute structures after the CNF file is completely parsed. */
void initialize_structures_after_parsing_CNF(void) {
  initialize_global_literals_to_clauses();
  initialize_neighborhoods();
}

/** @brief Compute structures after the assignment changes by more than one
 *         flip.
 */
void initialize_structures_after_assignment(void) {
  compute_unsat_after_assignment();
  compute_weight_reducing_after_assignment();
  initialize_neighborhoods_after_assignment();

  // Does not run anything if DEBUG is not enabled
  verify_after_assignment();
}

/** @brief Compute structures after the clause weights change by more than one
 *         weight transfer at a time.
 */
void initialize_structures_after_reweighting(void) {
  compute_weight_reducing_after_reweighting();
  initialize_neighborhoods_after_reweighting();

  // Does not run anything if DEBUG is not enabled
  verify_after_reweighting();
}
