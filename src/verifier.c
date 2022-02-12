/** @file verifier.c
 *  @brief A collection of verification methods that, unless specifically
 *         compiled in, are not included in the compilation of DDFW.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#include "xmalloc.h"
#include "verifier.h"
#include "assignment.h"
#include "neighborhood.h"
#include "weight_transfer.h"

/** If debug is enabled, implement the verifier functions */
#ifdef DEBUG

void verify_after_flip(void) {
  verify_clauses_and_assignment();
  verify_neighborhoods();
}

void verify_after_weight_transfer(void) {
  verify_neighborhoods();
}

void verify_after_assignment(void) {
  verify_clauses_and_assignment();
  verify_neighborhoods();
}

void verify_after_reweighting(void) {
  verify_neighborhoods();
}

#else

void verify_after_flip(void) { (void) 0; }
void verify_after_weight_transfer(void) { (void) 0; }
void verify_after_assignment(void) { (void) 0; }
void verify_after_reweighting(void) { (void) 0; }

#endif /* DEBUG */
