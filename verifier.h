/** @file verifier.h
 *  @brief A collection of verification methods that, unless specifically
 *         compiled in, are not included in the compilation of DDFW.
 *
 *  See verifier.c for implementation details.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#ifndef _VERIFIER_H_
#define _VERIFIER_H_

// Call during initialization stage
void initialize_verifier(void);

// For neighborhood.c
void verify_neighborhoods(void);

// For ddfw.c
void verify_crit_sat_unsat_weights(void);
void verify_clauses_and_assignment(void);
void verify_cost_reducing_vars(void);

#endif /* _VERIFIER_H_ */
