/** @file cnf_parser.h
 *  @brief Implements DIMACS CNF file parsing.
 *
 *  See cnf_parser.c for implementation details.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#ifndef _CNF_PARSER_H_
#define _CNF_PARSER_H_

#include "clause.h"

/** DIMACS character that starts a comment line. */
#define DIMACS_COMMENT_LINE ('c')

/** DIMACS character that starts a problem line. */
#define DIMACS_PROBLEM_LINE ('p')

void parse_cnf_file(char *filename);

#endif /* _CNF_PARSER_H_ */
