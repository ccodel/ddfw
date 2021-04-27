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

#define ERR_IF(cond, msg)  if (cond) { fprintf(stderr, msg); exit(-1); }

void verify_after_flip(void);
void verify_after_weight_transfer(void);
void verify_after_assignment(void);
void verify_after_reweighting(void);

#endif /* _VERIFIER_H_ */
