/** @file neighborhood.c
 *  @brief Implements functions to manage clause neighborhoods and track
 *         statistics about the neighborhoods.
 *
 *  DDFW transfers weight from neighboring clauses to unsatisfied clauses when
 *  no more cost-reducing variables exist. Variants on the original DDFW
 *  algorithm will take weight from greater subsets of neighboring clauses
 *  in various ways, and so standard tracking of neighborhood relations and
 *  aggregate data speeds up computation (e.g. by not having to recompute
 *  the sum of weights held in a single neighborhood).
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include "neighborhood.h"
#include "global_data.h"
#include "verifier.h"
#include "xmalloc.h"
#include "logger.h"

/** @brief The size of each neighborhood.
 *
 *  The array is 0-indexed by clause ID. The values are the number of clauses
 *  other than the considered clause that contain same-sign literals that are
 *  also present in the considered clause.
 */
int *neigh_sizes = NULL;

/** @brief The number of neighbors that are satisfied.
 *  
 *  The array is 0-indexed by clause ID.
 */
int *neigh_num_sat = NULL;

/** @brief The clause IDs that are neighbors.
 *
 *  The array is 0-indexed by clause ID.
 */
int **neigh_clauses = NULL;

/** @brief The sum total weight held by the neighborhood.
 *
 *  The array is 0-indexed by clause ID.
 */
weight *neigh_weights = NULL;

/** @brief The maximum weight held by any single clause in the neighborhood.
 *  
 *  The array is 0-indexed by clause ID.
 */
weight *neigh_max_weights = NULL;

/** @brief @brief The clause ID of the neighbor with maximum weight.
 *
 *  The array is 0-indexed by clause ID.
 */
int *neigh_max_idxs = NULL;

/** @brief Allocates the memory needed for neighborhood tracking structures.
 *
 *  Should be called after the CNF file has been parsed.
 */
void allocate_neighborhood_memory(void) {
  neigh_sizes = xmalloc(num_clauses * sizeof(int));
  neigh_num_sat = xmalloc(num_clauses * sizeof(int));
  neigh_clauses = xmalloc(num_clauses * sizeof(int *));
  neigh_weights = xmalloc(num_clauses * sizeof(weight));
  neigh_max_weights = xmalloc(num_clauses * sizeof(weight));
  neigh_max_idxs = xmalloc(num_clauses * sizeof(int));
}

/** @brief Initializes the neighborhood structures.
 *
 *  Should be called after the entirety of the CNF file has been parsed.
 *
 *  To compute the neighborhoods, an outer loop examines each clause, then an 
 *  inner loop goes through each literal and then those clauses associated 
 *  with each literal. If two literals are neighbors of the same clause, 
 *  then it is not added twice.
 */
void initialize_neighborhoods(void) {
  // Helper buffer for clause indices for neighbors
  int *neigh_buf = xmalloc(num_clauses * sizeof(int));
  int buf_idx = 0;

  // Alias some pointers
  int **literals = clause_literals;
  int *sizes = clause_sizes;

  // Loop through each clause and examine all literals in that clause
  for (int c = 0; c < num_clauses; c++) {
    const int clause_size = *sizes;
    int *lits = *literals;
    buf_idx = 0;

    // For each literal, check its neighboring clauses for neighbors
    for (int l = 0; l < clause_size; l++) {
      const int l_idx = *lits;
      const int occ = literal_occ[l_idx];
      int *lc = literal_clauses[l_idx];

      // Loop through all neighboring clauses for this lit and add if new
      for (int nc = 0; nc < occ; nc++) {
        const int neigh_idx = *lc;
        if (neigh_idx == c) {
          lc++;
          continue;
        }

        // Add to the neighbors found if new
        int already_found = 0;
        for (int i = 0; i < buf_idx; i++) {
          if (neigh_buf[i] == neigh_idx) {
            already_found = 1;
            break;
          }
        }

        if (!already_found) {
          neigh_buf[buf_idx] = neigh_idx;
          buf_idx++;
        }

        lc++;
      }
      
      lits++;
    }

    // Once done looping through literals in this clause, copy over neighbors
    neigh_sizes[c] = buf_idx;
    neigh_clauses[c] = xmalloc(buf_idx * sizeof(int));
    memcpy(neigh_clauses[c], neigh_buf, buf_idx * sizeof(int));

    literals++;
    sizes++;
  }

  xfree(neigh_buf);
}

void initialize_neighborhoods_after_assignment(void) {
  // Recompute which clauses in the neighborhoods are satisfied
  for (int c = 0; c < num_clauses; c++) {
    int num_sat_neigh = 0;
    int *nc = neigh_clauses[c];
    const int neigh_size = neigh_sizes[c];
    for (int n = 0; n < neigh_size; n++) {
      const int n_idx = nc[n];
      if (clause_num_true_lits[n_idx] > 0) {
        num_sat_neigh++;
      }
    }

    neigh_num_sat[c] = num_sat_neigh;
  }

  // Now compute all other information
  initialize_neighborhoods_after_reweighting();
}

void initialize_neighborhoods_after_reweighting(void) {
  // Recompute the weights and max weights of the neighborhoods
  for (int c = 0; c < num_clauses; c++) {
    int *nc = neigh_clauses[c];
    const int neigh_size = neigh_sizes[c];
    weight sum_neigh_weights = 0;
    weight max_neigh_weight = 0;
    int idx_max_weight = -1;

    for (int n = 0; n < neigh_size; n++) {
      const int n_idx = nc[n];
      if (clause_num_true_lits[n_idx] > 0) {
        const weight w = clause_weights[n_idx];
        sum_neigh_weights += w;

        if (w > max_neigh_weight) {
          max_neigh_weight = w;
          idx_max_weight = n_idx;
        }
      }
    }

    neigh_weights[c] = sum_neigh_weights;
    neigh_max_weights[c] = max_neigh_weight;
    neigh_max_idxs[c] = idx_max_weight;
  }
}

/** @brief Updates the neighborhood structures when a clause has a variable
 *         that's flipped.
 *
 *  Should only be called if the clause flips between 0 and 1 satisfied lits.
 */
void update_neighborhoods_on_flip(const int c_idx) {
  int *neighs = neigh_clauses[c_idx];
  const int neigh_size = neigh_sizes[c_idx];
  const weight w = clause_weights[c_idx];

  if (clause_num_true_lits[c_idx] == 0) {
    // Assume that the clause flipped from sat to unsat - remove weight, etc.
    for (int nc = 0; nc < neigh_size; nc++) {
      const int nc_idx = *neighs;
      neigh_num_sat[nc_idx]--;
      neigh_weights[nc_idx] -= w;
      if (neigh_max_idxs[nc_idx] == c_idx) {
        neigh_max_idxs[nc_idx] = -1;
      }

      neighs++;
    }
  } else {
    // Assume that the clause flipped from unsat to sat - add weight, etc.
    for (int nc = 0; nc < neigh_size; nc++) {
      const int nc_idx = *neighs;
      neigh_num_sat[nc_idx]++;
      neigh_weights[nc_idx] += w;
      if (neigh_max_weights[nc_idx] < w || 
          (neigh_max_weights[nc_idx] == w && neigh_max_idxs[nc_idx] == -1)) {
        neigh_max_weights[nc_idx] = w;
        neigh_max_idxs[nc_idx] = c_idx;
      }

      neighs++;
    }
  }

  verify_neighborhoods();
}


/** @brief Updates the neighborhood structures when a clause has weight
 *         transferred from it.
 *
 *  Should only be called on satisfied clauses when weight is taken away.
 *  Should only be called after the weight has been transferred.
 *
 *  @param c_idx  The index number of the clause
 *  @param diff   Positive weight, Weight taken away from c_idx.
 */
void update_neighborhoods_on_weight_transfer(const int c_idx, weight diff) {
  int *neighs = neigh_clauses[c_idx];
  const int neigh_size = neigh_sizes[c_idx];
  for (int nc = 0; nc < neigh_size; nc++) {
    const int nc_idx = *neighs;
    neigh_weights[nc_idx] -= diff;
    if (neigh_max_idxs[nc_idx] == c_idx) {
      neigh_max_idxs[nc_idx] = -1;
    }

    neighs++;
  }

  verify_neighborhoods();
}
