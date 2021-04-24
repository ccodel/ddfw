/** @file initializer.h
 *  @brief Handles initializing and resetting various structures.
 *
 *  See initializer.c for implementation details.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#ifndef _INITIALIZER_H_
#define _INITIALIZER_H_

void initialize_structures_after_parsing_CNF(void);
void initialize_structures_after_assignment(void);
void initialize_structures_after_reweighting(void);

#endif /* _INITIALIZER_H_ */
