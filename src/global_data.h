/** @file global_data.h
 *  @brief Global data structures shared across all files.
 *
 *  A compilation of all global data shared between the DDFW files. File-
 *  specific data structures and supporting variables should be declared
 *  within the .c file.
 *
 *  The data structures presented here should be data structures that are
 *  common across all DDFW implementations, i.e. structures that optimize
 *  should be relegated to other files.
 *
 *  See global_data.c for implementation details and documentation for each
 *  global variable and data structure.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#ifndef _GLOBAL_DATA_H_
#define _GLOBAL_DATA_H_

/** @brief Type of clause weights.
 *
 *  In the original DDFW paper, the weights assigned to each clause were
 *  integers. Here, they are taken as floating point values. However, to
 *  change back to integers, change the type below and recompile.
 *
 *  TODO potential optimization with integer computations.
 */
typedef double weight;

/******************************** GLOBAL MACROS *******************************/

/** @brief Takes a variable symbol from the DIMACS CNF input file and outputs
 *         an index into the literals array.
 *
 *  The literals array stores literal_t structs for l and !l. This macro
 *  indexes into the literals array by doubling positive variable symbols
 *  and adding one to the positive double of negated variable symbols.
 *
 *  @param  A variable symbol from the DIMACS CNF input file.
 *  @return An index into the literals array for that variable symbol.
 */
#define LIT_IDX(x)      (((x) < 0) ? (2 * -(x) + 1) : (2 * (x)))

/** @brief Takes an index into the literals array and determines if the
 *         literal index corresponds to a positive or a negated literal.
 *
 *  Because LIT_IDX doubles positive variable symbols and adds one to doubled
 *  negative variable symbols, a check against the least-significant bit of
 *  any literal index is sufficient to determine if the index corresponds to
 *  a literal or its negated twin.
 *
 *  @param  An index into the literals array.
 *  @return 1 if the literal index is a negated literal, 0 otherwise.
 */
#define IS_NEGATED(x)   ((x) & 0x1)

/** @brief Takes an index into the literals array and returns the index
 *         corresponding to the negated form of the literal.
 *
 *  To switch between negated forms of literal indices, even indices should
 *  have 1 added to themm, while odd indices should have 1 subtracted. An
 *  XOR against 0x1 does the trick.
 *
 *  @param  An index into the literals array.
 *  @return The index corresponding to the negated form of the literal
 *          index passed in.
 */
#define NEGATED_IDX(x)  ((x) ^ 0x1)

/** @brief Takes an index into the literals array and returns the literal
 *         index corresponding to the positive form of the literal.
 *
 *  For example, if 5 := !l were passed in, then 4 would be returned, but
 *  if 4 := l were passed in, then 4 would also be returned.
 *
 *  @param  An index into the literals array.
 *  @return The index corresponding to the positive form of the literal.
 */
#define POS_LIT_IDX(x)  ((x) & ~(0x1))

/** @brief Takes an index into the literals array and returns the positive
 *         variable index for that literal, e.g. to index into *assignment.
 *
 *  Because the literal indexes are just doubled variable symbols, just
 *  integer divide by 2 to get the positive variable value desired.
 *
 *  @param  An index into the literals array.
 *  @return The variable index for that literal, always positive.
 */
#define VAR_IDX(x)      ((x) / 2)

/** @brief Takes an index into the literals array and returns the variable
 *         symbol for that index, as in the one from the DIMACS CNF input file.
 *
 *  Just undo the mapping from LIT_IDX. Integer division shaves off any +1s.
 *
 *  @param  An index into the literals array.
 *  @return The variable symbol for that literal, as it appears in the DIMACS
 *          CNF input file.
 */
#define VAR_SYM(x)      (((x) & 0x1) ? -((x) / 2) : (x) / 2)

/** @brief Takes a lit_idx and retrieves the assignment for that literal.
 *
 *  Must divide by 2 to get rid of the doubling of l and !l in literals array.
 *
 *  @param x A literal index.
 *  @return 1 if the variable is set to true, and 0 otherwise.
 */
#define ASSIGNMENT(x)      (((x) & 0x1) ^ assignment[(x) / 2])

/****************************** GLOBAL DATA STRUCTURES ************************/
/** Various CNF information */
extern int num_vars;
extern int num_literals;
extern int num_clauses;

/** Global values about clause weights */
extern weight init_clause_weight;
extern weight total_sat_clause_weight;
extern weight total_unsat_clause_weight;

/** Information about the formula DDFW operates on */
extern long num_flips;
extern long num_flips_since_improvement;

extern char *assignment;
extern int *unsat_clauses;
extern int num_unsat_clauses;

extern char *best_assignment;
extern int best_num_unsat_clauses;
extern long best_flip_num;

/** Clause information, 0-indexed by clause ID (assigned at CNF parse) */
extern int *clause_sizes;
extern weight *clause_weights;
extern int *clause_num_true_lits;
extern int *clause_lit_masks;
extern int **clause_literals; // Stores literal indexes involved in the clause

/** Literal information, 1-indexed by way of LIT_IDX */
extern int *literal_occ;
extern int **literal_clauses; // Stores clause IDs the literal is a part of

/*************************** FUNCTIONS ****************************************/
/** Allocate the memory needed for global data structures */
void allocate_global_memory(void);

/** Allocates and computes literal -> clause structures above after CNF parse */
void initialize_global_literals_to_clauses(void);

/** Various functions to reset global data structures to default values */
void reset_clause_weights(void);

void set_unsat_weights(weight w);
void increase_unsat_weights(weight w);
void set_sat_weights(weight w);
void increase_sat_weights(weight w);

void initialize_global_data_for_alg_run(void);

#endif /* _GLOBAL_DATA_H_ */
