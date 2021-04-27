/** @file neighborhood.h
 *  @brief Implements functions to manage clause neighborhoods and track
 *         statistics about the neighborhoods.
 *
 *  See neighborhood.c for implementation details.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#ifndef _NEIGHBORHOOD_H_
#define _NEIGHBORHOOD_H_

#include "global_data.h"

// All of the below are 0-indexed by clause ID
extern int *neigh_sizes;
extern int *neigh_num_sat;      
extern int **neigh_clauses;       // Indexes (unsorted) of neigh clauses
extern weight *neigh_weights;     // Sum of weights of sat neigh clauses
extern weight *neigh_max_weights; // Weight of max sat neigh clause
extern int *neigh_max_idxs;       // Index number of max, -1 if changed

void allocate_neighborhood_memory(void); // Call after header parsed
void initialize_neighborhoods(void);     // Call after entire CNF is parsed

void initialize_neighborhoods_after_assignment(void);
void initialize_neighborhoods_after_reweighting(void);

void update_neighborhoods_on_flip(const int c_idx);
void update_neighborhoods_on_weight_transfer(const int c_idx, const weight d);

void verify_neighborhoods();

#endif /* _NEIGHBORHOOD_H_ */
