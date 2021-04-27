/** @file assignment.h
 *  @brief Implements functions to initialize and change the Boolean variable
 *         truth assignment while running DDFW.
 *
 *  See assignment.c for implementation details.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#ifndef _ASSIGNMENT_H_
#define _ASSIGNMENT_H_

void allocate_assignment_memory(void);

void compute_unsat_after_assignment(void);

void generate_random_assignment(void);
void flip_variable(const int v_idx);
void restore_to_best_assignment(void);

void verify_clauses_and_assignment(void);

/** Implementation of DDFW+ weight toggle */
/*
extern long ddfw_plus_counter;
extern int ddfw_plus_boolean;
extern int ddfw_reweighting_enabled;
void reset_to_ddfw_plus_weightings(void);
*/

#endif /* _CLAUSE_H_ */
