/** @file cnf_parser.c
 *  @brief Implements DIMACS CNF file parsing.
 *
 *  When DDFW is run, a DIMACS CNF file is required. cnf_parser.c will
 *  parse the literals and clauses in the file and pass them off to
 *  clause.c for allocation and initialization.
 *
 *  No file format beyond DIMACS CNF is currently supported by this
 *  implementation of DDFW.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#include <stdlib.h>
#include <stdio.h>

#include "logger.h"
#include "cnf_parser.h"

/** A buffer for literal indices. Used to store the literals for each
 *  clause when parsing the CNF file.
 *
 *  TODO consider hijacking the reducing_cost_lits extern variable
 */
static int clause_idx_buf[MAX_CLAUSE_LENGTH];

/** @brief Parses a CNF file with the provided file name.
 *
 *  @param filename The name of the CNF file to open.
 *  @return void, but will call exit(-1) if an error is encountered.
 */
void parse_cnf_file(char *filename) {
  FILE *f = fopen(filename, "r");
  if (f == NULL) {
    fprintf(stderr, "c Error: Unable to open the file: %s\n", filename);
    exit(-1);
  } else {
    log_str("c Successfully opened file %s\n", filename);
  }

  // Keep scanning the file until the problem line is found
  int problem_found = 0;
  while (!problem_found) {
    char c = fgetc(f);
    switch (c) {
      case DIMACS_COMMENT_LINE:
        log_str("c Found comment line, skipping\n");
        fscanf(f, "%*[^\n]");
        continue;
      case DIMACS_PROBLEM_LINE:
        log_str("c Found problem line, ");
        int num_clauses, num_vars;
        fscanf(f, "%*s %d %d\n", &num_vars, &num_clauses);
        log_str("file claims to have %d variables and %d clauses\n",
            num_vars, num_clauses);
        initialize_formula(num_clauses, num_vars);
        problem_found = 1;
        break;
    }
  }

  // Now scan in all the clauses
  log_str("c Scanning clauses\n");
  int scanned_clauses = 0;
  while (scanned_clauses < num_clauses) {
    int var = 0, clause_size = 0;
    do {
      fscanf(f, "%d", &var);
      if (var != 0) {
        // New literal for the clause - count occurrence, place into buffer
        int lit_idx = LIT_IDX(var);
        literals[lit_idx].occurrences++;
        clause_idx_buf[clause_size++] = lit_idx;
      } else {
        // Read in a 0, so we are done readig this clause
        // Initialize a new clause
        initialize_clause(scanned_clauses++, clause_size, clause_idx_buf);
      }
    } while (var != 0);
  }

  log_str("c Done with I/O, post-processing CNF file\n");
  process_clauses();
  log_str("c Done post-processing, onto the algorithm!\n");
  fclose(f);
}
