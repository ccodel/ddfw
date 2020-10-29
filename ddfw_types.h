/** @file ddfw_types.h
 *  @brief Common types used by the DDFW algorithm.
 *
 *  NOTE: These are "spiritual" types, in the sense that
 *  accessing everything as a struct is more inefficient.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *  
 *  @bug No known bugs.
 */

#ifndef _DDFW_TYPES_H_
#define _DDFW_TYPES_H_

typedef struct sat_clause {
  int size;
  int sat_lits;
  int sat_mask;
  double weight;
  int *literals;
} clause_t;

typedef struct sat_literal {
  int occurrences;
  int marking;
  int *clause_indexes;
} literal_t;

#endif /* _DDFW_TYPES_H_ */
