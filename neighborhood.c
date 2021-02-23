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

#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include "clause.h"
#include "xmalloc.h"
#include "neighborhood.h"


int *neigh_sizes;
int *neigh_num_sat;
int **neigh_clauses;
double *neigh_weights;
double *neigh_max_weights;
int *neigh_max_idxs;


/** @brief Allocates the memory needed for neighborhood tracking structures.
 *
 *  Should be called after the CNF file has been parsed for the number of
 *  clauses and variables.
 */
void allocate_neighborhoods(void) {
  neigh_sizes = xmalloc(num_clauses * sizeof(int));
  neigh_num_sat = xmalloc(num_clauses * sizeof(int));
  neigh_clauses = xmalloc(num_clauses * sizeof(int *));
  neigh_weights = xmalloc(num_clauses * sizeof(double));
  neigh_max_weights = xmalloc(num_clauses * sizeof(double));
  neigh_max_idxs = xmalloc(num_clauses * sizeof(int));
}


/** @brief Computes the neighborhood structure between clauses.
 *
 *  Should be called after the entirety of the CNF file has been parsed, and
 *  the clause and literal information has been initialized by other modules.
 *
 *  TODO Currently, to compute the neighborhoods, an outer loop examines each
 *  clause, then an inner loop goes through each literal and then those clauses
 *  associated with each literal. If two literals are neighbors of the same
 *  clause, then it is not added twice.
 *
 *  TODO because the operation is O(n^2) (or abouts). However, since this is
 *  called once each run, it's excusable.
 */
void compute_neighborhoods(void) {
  // Helper buffer for clause indices for neighbors
  int *neigh_buf = xmalloc(num_clauses * sizeof(int));
  int buf_idx = 0;

  // Alias some pointers for speeup in the loops
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

  free(neigh_buf);
}


/** @brief Initializes the neighborhood statistics, given an initial assignment.
 *  
 *  Assumes that generate_initial_assignment(), or something similar, has been
 *  called, and the numbers of unsatisfied literals and such have been placed
 *  into their appropriate global variables.
 *
 *  TODO another expensive O(n^2) computation, see above
 */
void initialize_neighborhoods(void) {
  // Zero out some fields before computation
  memset(neigh_num_sat, 0, num_clauses * sizeof(int));
  memset(neigh_weights, 0, num_clauses * sizeof(double));
  memset(neigh_max_weights, 0, num_clauses * sizeof(double));
  memset(neigh_max_idxs, UINT8_MAX, num_clauses * sizeof(int));

  // Alias pointers to be incremented in the loop for speedup
  int **ncs = neigh_clauses;
  int *ns = neigh_sizes;
  int *nns = neigh_num_sat;
  double *nw = neigh_weights;
  double *nmw = neigh_max_weights;
  int *nmi = neigh_max_idxs;

  // Loop through each clause and raw compute
  for (int c = 0; c < num_clauses; c++) {
    int *clauses = *ncs;
    int neigh_size = *ns;
    if (neigh_size == 0) goto next_clause; // Skip stats if no neighborhood

    // Computed stats are: num sat, sum of weights, max weight
    int num_sat_neigh = 0;
    double sum_neigh_weights = 0.0;
    double max_neigh_weight = 0.0;
    int idx_max_weight = -1;

    for (int nc = 0; nc < neigh_size; nc++) {
      const int nc_idx = *clauses;
      if (clause_num_true_lits[nc_idx] > 0) {
        num_sat_neigh++;
        const double w = clause_weights[nc_idx];
        sum_neigh_weights += w;
        if (w > max_neigh_weight) {
          max_neigh_weight = w;
          idx_max_weight = nc_idx;
        }
      }

      clauses++;
    }

    // Take the information learned from the neighboring clauses and update
    *nns = num_sat_neigh;
    *nw = sum_neigh_weights;
    *nmw = max_neigh_weight;
    *nmi = idx_max_weight;

next_clause: // Increment all the pointers to move onto the next clause
    ncs++; ns++; nns++; nw++; nmw++; nmi++;
  }
}


/** @brief Updates the neighborhood structures when a clause has a variable
 *         that's flipped.
 *
 *  Should only be called if the clause flips between 0 and 1 satisfied lits.
 */
void update_neighborhood_on_flip(const int c_idx) {
  int *neighs = neigh_clauses[c_idx];
  const int neigh_size = neigh_sizes[c_idx];
  const double w = clause_weights[c_idx];

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
      if (neigh_max_weights[nc_idx] < w) {
        neigh_max_weights[nc_idx] = w;
        neigh_max_idxs[nc_idx] = c_idx;
      }

      neighs++;
    }
  }
}


/** @brief Updates the neighborhood structures when a clause has weight
 *         transferred from it.
 *
 *  Should only be called on satisfied clauses when weight is taken away.
 *  Should only be called after the weight has been transferred.
 *
 *  @param c_idx  The index number of the clause
 *  @param diff   Positive double, Weight taken away from c_idx.
 */
void update_neighborhood_on_weight_transfer(const int c_idx, double diff) {
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
}
