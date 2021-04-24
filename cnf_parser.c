/** @file cnf_parser.c
 *  @brief Implements DIMACS CNF file parsing.
 *
 *  When DDFW is run, a DIMACS CNF file is required. cnf_parser.c will
 *  parse the literals and clauses in the file and pass them off to
 *  clause.c for allocation and initialization.
 *
 *  No other CNF file format is currently supported by this implementation.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "global_data.h"
#include "assignment.h"
#include "neighborhood.h"
#include "weight_reducer.h"
#include "initializer.h"
#include "xmalloc.h"
#include "cnf_parser.h"
#include "logger.h"

/** DIMACS character that starts a comment line. */
#define DIMACS_COMMENT_LINE ('c')

/** DIMACS character that starts a problem line. */
#define DIMACS_PROBLEM_LINE ('p')

/** @brief A buffer for literal IDs. Used to store the literals for each
 *         clause when parsing the CNF file.
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
        fscanf(f, "%*s %d %d\n", &num_vars, &num_clauses);
        num_literals = num_vars * 2;
        log_str("c File claims to have %d variables and %d clauses\n",
            num_vars, num_clauses);
        problem_found = 1;
        break;
    }
  }

  // Allocate memory after parsing header
  allocate_global_memory();
  allocate_assignment_memory();
  allocate_neighborhood_memory();
  allocate_weight_reducer_memory();

  // Now scan in all the clauses
  int scanned_clauses = 0;
  while (scanned_clauses < num_clauses) {
    int var = 0, clause_size = 0;
    do {
      fscanf(f, "%d", &var);
      if (var != 0) {
        // New literal for the clause - count occurrence, place into buffer
        int lit_idx = LIT_IDX(var);
        literal_occ[lit_idx]++;
        clause_idx_buf[clause_size++] = lit_idx;
      } else {
        // We have read in a 0, so we are done reading this clause
        // Set its information (size, literals, etc.)
        clause_sizes[scanned_clauses] = clause_size;
        int *lits = xmalloc(clause_size * sizeof(int));
        memcpy(lits, clause_idx_buf, clause_size * sizeof(int));
        clause_literals[scanned_clauses] = lits;
        scanned_clauses++;
      }
    } while (var != 0);
  }

  fclose(f);
  log_str("c Done parsing, now post-processing CNF file\n");

  initialize_structures_after_parsing_CNF();

  log_str("c Done post-processing, onto the algorithm!\n");
}
