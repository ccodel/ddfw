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

// All of the below are 0-indexed by clause number
extern int *neigh_sizes;
extern int *neigh_num_sat;      
extern int **neigh_clauses;       // Indexes (unsorted) of neigh clauses
extern double *neigh_weights;     // Sum of weights of sat neigh clauses
extern double *neigh_max_weights; // Weight of max sat neigh clause
extern int *neigh_max_idxs;       // Index number of max, -1 if changed

void allocate_neighborhoods(void);     // Call after header parsed in
void compute_neighborhoods(void);      // Call after clauses parsed in
void initialize_neighborhoods(void);   // Call after initial var assignment

void update_neighborhood_on_flip(const int c_idx);
void update_neighborhood_on_weight_transfer(const int c_idx, const double diff);

#endif /* _NEIGHBORHOOD_H_ */
