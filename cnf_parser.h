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

/** @brief The maximum number of variables any one clause is assumed to have.
 *
 *  A cap on each clause length allows a finite buffer to be allocated at
 *  compile time to simplify the parsing process. See cnf_parser.c for how
 *  the buffer is used.
 *
 *  TODO remove this limit later
 */
#define MAX_CLAUSE_LENGTH 10000

void parse_cnf_file(char *filename);

#endif /* _CNF_PARSER_H_ */
