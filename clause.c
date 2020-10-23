/** @file clause.c
 *  @brief Implements functions to create and interact with CNf DDFW clauses.
 *
 *  TODO
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#include <stdlib.h>
#include <string.h>

#include "ddfw_types.h"
#include "clause.h"
#include "logger.h"

#ifndef ABS
#define ABS(x)   (((x) < 0) ? -(x) : (x))
#endif

#ifndef MIN
#define MIN(x, y)  (((x) < (y)) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x, y)  (((x) > (y)) ? (x) : (y))
#endif

int num_literals = 0;
int num_clauses = 0;
int num_vars = 0;
int num_reducing_cost_lits = 0;
int num_zero_cost_lits = 0;
double max_reducing_cost = 0.0;
int flips = 0;
int unsat_clauses = 0;

char *assignment = NULL;
literal_t *literals = NULL;
clause_t *clauses = NULL;

/** @brief The array containing indexes of literals that, if flipped,
 *         reduce the cost of the unsatisfied literals.
 *
 *  The exact way that the list is calculated and a literal is chosen from
 *  within is uncertain in the paper. I have determined that those literals
 *  that reduce cost a strictly positive amount should be prioritized over
 *  those that result in no weight change. Therefore, those literals that
 *  reduce cost are placed "in the front" of the array, while those literals
 *  that result in no weight change on a flip are placed "in the back" of
 *  the array. The variables num_reducing_cost_lits and num_zero_cost_lits
 *  are updated accordingly.
 */
int *reducing_cost_lits = NULL;

/** @brief Initializes the global formula to get ready for the specified
 *         number of clauses and variables/literals. Allocates the necessary
 *         memory.
 *
 *  @param num_cs The number of clauses in the CNF input file.
 *  @param num_vs The number of variables in the CNF input file.
 */
void initialize_formula(int num_cs, int num_vs) {
  // Set global formula variables
  num_vars = num_vs;
  num_literals = 2 * num_vars;
  num_clauses = num_cs;
  num_reducing_cost_lits = 0;
  flips = 0;
  unsat_clauses = 0;

  log_str("c Allocating memory for %d variables and %d clauses\n",
      num_vars, num_clauses);

  // The assignment WAS a bitvector, so round up to a multiple of BITS_IN_BYTE
  // The addition of BITS_IN_BYTE ensures that integer division does not
  // "cut off" var truth values. Note that assignment is 1-indexed.
  //assignment = malloc(
  //    ((num_vars + BITS_IN_BYTE) / BITS_IN_BYTE) * sizeof(char));
  assignment = malloc((num_vars + 1) * sizeof(char));
  if (assignment == NULL) {
    log_str("c Ran out of memory when allocating assignment, exiting\n");
    exit(-1);
  }

  // TODO repeat code
  // Need calloc to zero out pointers, see process_clauses()
  // Need to add two to allow for 1-indexing of the literals array
  literals = calloc((num_literals + 2), sizeof(literal_t));
  if (literals == NULL) {
    log_str("c Ran out of memory when allocating literals, exiting\n");
    exit(-1);
  }

  clauses = malloc(num_clauses * sizeof(clause_t));
  if (clauses == NULL) {
    log_str("c Ran out of memory when allocating clauses, exiting\n");
    exit(-1);
  }

  // Zero out because array comparison, TODO to malloc later
  reducing_cost_lits = calloc(num_literals, sizeof(int));
  if (reducing_cost_lits == NULL) {
    log_str("c Ran out of memory when allocating extras, exiting\n");
    exit(-1);
  }

  log_str("c Allocated memory successfully\n");
}

/** @brief Initializes the clause at the specified index.
 *
 *  A clause is defined as a clause_t struct in the clauses array.
 *  TODO
 *
 *  @param clause_idx  The index into the formula->clauses array.
 *  @param size        The number of literals in the clause.
 *  @param lit_idxs    The indexes into formula->literals for the clause's lits.
 */
void initialize_clause(int clause_idx, int size, int *lit_idxs) {
  clause_t *c = &clauses[clause_idx];
  c->size = size;
  c->weight = DEFAULT_CLAUSE_WEIGHT;
  c->literals = malloc(size * sizeof(int));
  if (c->literals == NULL) {
    log_str("c Ran out of memory when allocating clause literals, exiting\n");
    exit(-1);
  }

  memcpy(c->literals, lit_idxs, size * sizeof(int));
}

/** @brief Processes the clauses read in from the CNF input file and before
 *         the DDFW algorithm is run.
 *
 *  Some post-processing that is done is
 *    - Set up clause pointers in each literal TODO efficiency?
 *  
 *  @return void, but calls exit(-1) on memory allocation failure.
 */
void process_clauses() {
  // Alias variables
  const int nc = num_clauses; // TODO remove?
  const int nl = num_literals; // TODO remove?

  // Loop through each literal in each clause and add the clause to the literal
  // Store literal in clause information in l->clause_indexes.
  // Hijack the l->marking int for clause index, reset to 0 at end
  for (int c = 0; c < nc; c++) {
    clause_t *cl = &clauses[c];
    const int size = cl->size;
    int *lits = cl->literals;
    for (int l = 0; l < size; l++) {
      int l_idx = lits[l];
      literal_t *lit = &literals[l_idx];
      
      // If the literal hasn't been touched yet, allocate space for clause idxs
      if (lit->clause_indexes == NULL) {
        lit->clause_indexes = malloc(lit->occurrences * sizeof(int));
        if (lit->clause_indexes == NULL) {
          log_str("c Error: Unable to allocate memory for clause indexes\n");
          exit(-1);
        }
        lit->marking = 0;
      }

      // Add the clause index to the literal
      lit->clause_indexes[lit->marking] = c;
      lit->marking++;
    }
  }

  // Set all lit->marking variables to 0
  for (int l = 0; l <= nl; l++) {
    literal_t *lit = &literals[l];
    lit->marking = 0;
  }
}

/** @brief Generates a random variable assignment for the global formula.
 *
 *  A random assignment is chosen for the num_vars variables. The random
 *  assignment is generated by storing random numbers into the assignment
 *  bitvector in the global formula. A bit of 1 indicates that the variable
 *  is set to true and its negation is false.
 *
 *  The clauses and literals are updated to reflect the random assignment,
 *  which includes
 *    
 *    - Sat
 *    - Mask
 *    - TODO
 *
 *  Note that the function does not initialize the srand() randomization
 *  library. Such randomization should be handled at command-line parsing.
 *
 *
 *   TODO clean up docs
 *  1. A random assignment is chosen for the num_vars variables.
 *     The random assignment is generated by looping over each
 *     "slot" in the assignment variable and assigning it a random
 *     bit. A 1 in slot i indicates that x_i is assigned to true.
 *
 *  2. The clauses are looped over to determine satisfiability.
 *     After a random assignment is generated, many clauses may evaluate
 *     to true. Looping over the clauses determines how many literals in
 *     each clause are true, and overall how many clauses are true.
 *     Note the usage of an XOR satisfiability mask in each clause to help
 *     determine which literal is the last one keeping a clause true in
 *     the case that clause->sat_lits is 1.
 */
void generate_random_assignment() {
  log_str("c Randomizing assignment\n");

  // Give random bits to assignment
  // TODO cast to int* for BITS_IN_BYTE fewer calls to rand()
  //const int slots = (num_vars + BITS_IN_BYTE) / BITS_IN_BYTE;
  for (int i = 1; i <= num_vars; i++) {
    assignment[i] = (char) (rand() & 0x1);
  }

  // TODO think about whether clause -> literal is faster than literal -> clause
  // Iterate through the clauses to see which are satisfied
  const int nc = num_clauses;
  for (int i = 0; i < nc; i++) {
    clause_t *c = &clauses[i];
    c->sat_lits = 0;
    c->sat_mask = 0;

    // Loop through the literals to see how many are satisfied
    int *lits = c->literals;
    const int size = c->size;
    for (int l = 0; l < size; l++) {
      int lit_idx = lits[l];
      int negated = IS_NEGATED(lit_idx);
      int assigned = ASSIGNMENT(lit_idx);

      // If the literal evaluates to true, update clause information
      if ((assigned && !negated) || (!assigned && negated)) {
        c->sat_lits++;
        c->sat_mask ^= lit_idx;
      }
    }

    if (c->sat_lits == 0) {
      unsat_clauses++;
    }
  }

  log_assignment();
}

/** @brief Takes an index of a literal to flip and flips it in the assignment.
 *
 *  The function will change the truth value of the variable corresponding to
 *  the provided lit_idx.
 *
 *  Any internal changes as a result of flipping the variable (weights, marking
 *  literals, etc.) also occur here.
 *
 *  @param lit_idx The index into the literals array of the literal to flip.
 */
void flip_literal(int lit_idx) {
  log_str("c Flipping %d\n", lit_idx);

  const int pos_lit_idx = POS_LIT_IDX(lit_idx);
  const int not_lit_idx = NEGATED_IDX(pos_lit_idx);
  const int var_idx = VAR_IDX(lit_idx);
  int assigned = ASSIGNMENT(lit_idx);

  // TODO Clear markings in the positive literal
  literal_t *l, *not_l;
  if (assigned) {
    l = &literals[pos_lit_idx];
    CLEAR_MARKING(l);
    not_l = &literals[not_lit_idx];
  } else {
    l = &literals[not_lit_idx];
    not_l = &literals[pos_lit_idx];
    CLEAR_MARKING(not_l);
  }

  int occ = l->occurrences;
  int not_occ = not_l->occurrences;

  // For each clause containing l, set l to false
  for (int c = 0; c < occ; c++) {
    int c_idx = l->clause_indexes[c];
    clause_t *cl = &clauses[c_idx];

    // Remove ltieral from sat XOR mask
    cl->sat_lits--;
    cl->sat_mask ^= lit_idx;

    // Clear markings of literals involved in the clause
    const int size = cl->size;
    for (int cl_lit = 0; cl_lit < size; cl_lit++) {
      int l_idx = POS_LIT_IDX(cl->literals[cl_lit]);
      literal_t *clause_literal = &literals[l_idx];
      CLEAR_MARKING(clause_literal);
    }

    if (cl->sat_lits == 0) {
      unsat_clauses++;
    }
  }

  // For each clause containing !l, set !l to true
  for (int c = 0; c < not_occ; c++) {
    int c_idx = not_l->clause_indexes[c];
    clause_t *cl = &clauses[c_idx];

    // Add literal to sat XOR mask
    cl->sat_lits++;
    cl->sat_mask ^= lit_idx;

    // Clear markings of literals involved in the clause
    const int size = cl->size;
    for (int cl_lit = 0; cl_lit < size; cl_lit++) {
      int l_idx = POS_LIT_IDX(cl->literals[cl_lit]);
      literal_t *clause_literal = &literals[l_idx];
      CLEAR_MARKING(clause_literal);
    }

    if (cl->sat_lits == 1) {
      unsat_clauses--;
    }
  }

  // Affect assignment bitvector
  assignment[var_idx] = !assignment[var_idx];
  // int shift = var_idx & BYTE_MASK;
  // int bit = 1 << shift;
  // assignment[var_idx / BITS_IN_BYTE] ^= bit;

  log_str("c There are %d remaining unsatisfied clauses\n", unsat_clauses);
  if (get_verbosity() == VERBOSE) {
    log_assignment();
  }

  flips++;
}
