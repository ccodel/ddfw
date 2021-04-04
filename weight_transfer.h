/** @file weight_transfer.h
 *  @brief Handles the transfer of weight from satisfied to unsatisfied
 *         clauses when no more cost-reducing variables exist.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#ifndef _WEIGHT_TRANSFER_H_
#define _WEIGHT_TRANSFER_H_

/** @brief Defines how many clauses receive a weight transfer for each
 *         unsatisfied clause, when no more cost-reducing variables exist.
 *
 *  SINGULAR indicates that weight is transferred from a single clause.
 *  ABOVE_INIT indicates that only those clauses with weight above the initial
 *    weight have weight transferred from them.
 *  ALL indicates that weight is transferred from all neighboring clauses.
 */
typedef enum weight_transfer_group {
  SINGULAR, ABOVE_INIT, ALL, ACTIVE
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


// Public options - toggled at command line
extern transfer_group_t transfer_grp;
extern rule_group_t rule_grp;
extern available_weight_option_t available_weight_opt;

extern int num_weight_transfers;

// Functions
void distribute_weights(void);
double get_transfer_average(void);

#endif /* _WEIGHT_TRANSFER_H_ */
