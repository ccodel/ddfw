/** @file clause.h
 *  @brief Implements functions to create and interact with CNF DDFW clauses.
 *
 *  See clause.c for implementation details.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#ifndef _CLAUSE_H_
#define _CLAUSE_H_

/** A cap on each clause length allows a finite buffer to be allocated at
 *  compile time to simplify the parsing process. See cnf_parser.c for how
 *  the buffer is used.
 *
 *  TODO remove this limit later
 */
#define MAX_CLAUSE_LENGTH 10000

/** @brief Default starting weight for each clause.
 *  TODO make it so this can be toggled with a command line argument.
 */
#define DEFAULT_CLAUSE_WEIGHT (100.0)

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

/** @brief Defines the number of bites in a byte.
 *
 *  The standard 8 bits in a byte are used to build the assignment bitvector,
 *  with a truth value for each variable corresponding to a single bit in
 *  the assignment.
 *
 *  Note that if the bitvector is instead changed to a char-vector, this
 *  macro is (probably) no longer needed.
 */
// #define BITS_IN_BYTE    8

/** @brief Bitmask to get a number in the range [0, BITS_IN_BYTE).
 *
 *  The byte mask is used to find the position in the assignment bitvector
 *  to flip. If the bitvector changes to a char-vector, this macro is
 *  (probably) no longer needed.
 */
// #define BYTE_MASK      0x7

/** @brief Takes a lit_idx and retrieves the assignment for that literal.
 *  
 *  Must divide by 2 to get rid of the doubling of l and !l in literals array.
 *
 *  @param x A literal index.
 *  @return 1 if the variable is set to true, and 0 otherwise.
 */
#define ASSIGNMENT(x)      (assignment[(x) / 2])
// #define ASSIGNMENT(x)  (((assignment[(x) / (2 * BITS_IN_BYTE)]) \
//    >> (((x) / 2) & BYTE_MASK)) & 0x1)

/** Global variables for the single formula DDFW is solving. */
// TODO does packaging into a struct so there is one extern too slow?

/** 
 *
 *  In order to facilitate fast DDFW solving, a CNF formula includes not just
 *  the literals and clauses involved, but also miscellaneous other structures
 *  that help amortize repeated helper function calls in the loop body.
 */
// CNF information
extern int num_vars;
extern int num_literals;
extern int num_clauses;

// Statistics
extern int num_restarts;
extern int num_flips;
extern double unsat_weight;

// Formula information - 1-indexed (VAR_IDX indexed)
extern char *assignment;

// Clause information - 0 indexed
extern int *clause_sizes;
extern double *clause_weights;
extern int *clause_num_true_lits;
extern int *clause_lit_masks;
extern int **clause_literals; // Index by clause num, then by literal num

// Literal information - use LIT_IDX
extern int *literal_occ;
extern int **literal_clauses; // Index by literal num, then by clause num

// Bookkeeping structures
// Membership struct for false clauses
extern int *false_clause_members;
extern int *false_clause_indexes;
extern int num_unsat_clauses;

// Membership struct for cost reducing literals
extern int *cost_reducing_lits; // Make into vars?
extern int *cost_reducing_idxs;
extern int num_cost_reducing_lits;
void add_cost_reducing_lit(int l_idx); // Helper

// Membership struct for variables to compute cost reduction
extern int *cost_compute_vars;
extern int *cost_compute_idxs;
extern int num_cost_compute_vars;
void add_cost_compute_var(int v_idx);
void remove_cost_compute_var(int v_idx); // Helper

// Functions
void initialize_formula(int num_cs, int num_vs);
void initialize_clause(int clause_idx, int size, int *lit_idxs);
void process_clauses();

void generate_random_assignment();
void flip_literal(int lit_idx);

#endif /* _CLAUSE_H_ */
