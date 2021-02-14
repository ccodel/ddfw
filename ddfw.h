/** @file ddfw.h
 *  @brief An implementation of the DDFW algorithm, as presented in
 *         the paper: "Neighbourhood Clause Weight Redistribution in
 *         Local Search for SAT" [Ishtaiwi, Thornton, Sattar, Pham].
 *
 *  See ddfw.c for implementation details.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug See ddfw.c for bug listing.
 */

#ifndef _DDFW_H_
#define _DDFW_H_

/** @brief The version number for the current implementation of DDFW. */
#define VERSION 3

/** @brief Default randomization seed.
 *
 *  Can be toggled with the -s <seed> flag at the command line.
 */
#define DEFAULT_RAND_SEED  0xdeadd00d

/** @brief Defines how the algorithm times out.
 *
 *  The algorithm can time out based on CPU time or number of flips, or both.
 *
 *  DEFAULT is the default behavior, which times out after DEFAULT_TIMEOUT_SECS.
 *  TIME makes the algorithm time out after a number of seconds.
 *  FLIPS makes the algorithm time out after a number of flips.
 *  BOTH makes the algorithm time out after flips or time, whichever is first.
 */
typedef enum timeout_method {
  DEFAULT, TIME, FLIPS, BOTH
} timeout_method_t;

/** @brief Defines the selection method for the flipped variable.
 *
 *  UNIFORM is a uniform distribution across the candidate variables.
 *  WEIGHTED is a weighted probability distribution according to weights.
 *  BEST is the best cost-reducing variable.
 */
typedef enum variable_selection_method {
  UNIFORM, WEIGHTED, BEST
} selection_method_t;

/** @brief Defines how many clauses receive a weight transfer for each
 *         unsatisfied clause, when no more cost-reducing variables exist.
 *
 *  SINGULAR indicates that weight is transferred from a single clause.
 *  ABOVE_INIT indicates that only those clauses with weight above the initial
 *    weight have weight transferred from them.
 *  ALL indicates that weight is transferred from all neighboring clauses.
 */
typedef enum weight_transfer_group {
  SINGULAR, ABOVE_INIT, ALL
} transfer_group_t;

/** @brief Defines how much weight the weight transfer rule is applied to.
 *
 *  SINGLUAR indicates that the rule is applied on a per-clause basis.
 *  SUM indicates that the rule is applied on the sum of candidate weights.
 *    The weight is then transferred proportional to the fraction of the
 *    total weight each candidate clause has.
 *  AVG indicates that the rule is applied on the average of the candidate
 *    weights, and distributed evenly (not proportionally).
 */
typedef enum weight_rule_group {
  RULE_SINGULAR, RULE_SUM, RULE_AVG
} rule_group_t;

/** @brief Defines how much weight is available for transferral.
 *
 *  RAW indicates that all weight in the clause is available.
 *  MINUS_INIT indicates that only weight above the initial weight is available.
 */
typedef enum available_weight_option {
  RAW, MINUS_INIT
} available_weight_option_t;

extern int algorithm_run;
extern timeout_method_t timeout_method;
extern selection_method_t selection_method;
extern int timeout_secs;
extern int timeout_flips;
extern double mult_a;
extern double mult_A;
extern double add_c;
extern double add_C;
extern int suppress_solution;
extern int weight_statistics_log_rate;

extern transfer_group_t transfer_grp;
extern rule_group_t rule_grp;
extern available_weight_option_t available_weight_opt;

void run_ddfw_algorithm();

#endif /* _DDFW_H_ */
