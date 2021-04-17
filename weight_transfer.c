/** @file weight_transfer.c
 *  @brief Handles the transfer of weight from satisfied to unsatisfied
 *         clauses when no more cost-reducing variables exist.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 */

#include <stdlib.h>
#include <stdio.h>

#include "weight_transfer.h"
#include "neighborhood.h"
#include "clause.h"
#include "ddfw.h"

/** @brief The probability that the maximum weight neighboring clause is
 *         disregarded for a random clause when distributing weights.
 *
 *  Because the pseudo-random number generator only returns integers, an
 *  integer representation for the probability was chosen. The random number
 *  is taken mod DEN, and compared to be less than NUM.
 *
 *  While not mentioned in the original DDFW paper, the implementation provided
 *  by the authors included a probability to disregard the maximum weight
 *  neighbor. The value from the original implementation was 1%.
 */
#define RAND_NEIGH_NUMERATOR    1
#define RAND_NEIGH_DENOMINATOR  100

transfer_group_t transfer_grp;
rule_group_t rule_grp;
available_weight_option_t available_weight_opt;

int num_weight_transfers = 0;

#define WEIGHT_TRANSFER_MEMORY 1000
static double recent_transfer_weights[WEIGHT_TRANSFER_MEMORY];
static int recent_transfer_idx = 0;

/** @brief Returns the index of a random satisfied clause.
 *  
 *  The satisfied clause has at least the initial weight.
 */
static int get_random_satisfied_clause() {
  int loop_counter = 0;
  int loop_stop = 100 * num_clauses;
  int idx;
  do {
    idx = ((unsigned int) rand()) % num_clauses;
    loop_counter++;
  } while ((clause_num_true_lits[idx] == 0 ||
        clause_weights[idx] < init_clause_weight) && loop_counter < loop_stop);

  return idx;
}

/** @brief Applies the linear rule for weight transferral. Returns the amount
 *         of weight to be transfered.
 *  
 *  Applies the linear rule of
 *    
 *    mult_a * w + add_c
 *
 *  Depending on if w < init_clause_weight or w >= init_clause_weight.
 *  Also if the RAW or MINUS_INIT option is specified.
 *
 *  @param Raw weight potentiall available for a transfer.
 */
static double calculate_transfer_weight(double w) {
  switch (available_weight_opt) {
    case RAW:
      if (w <= init_clause_weight) {
        return w - ((mult_A * w) + add_C);
      } else {
        return w - ((mult_a * w) + add_c);
      }
    case MINUS_INIT:
      if (w <= init_clause_weight) {
        return (w - ((mult_A * w) + add_C)) - init_clause_weight;
      } else {
        return (w - ((mult_a * w) + add_c)) - init_clause_weight;
      }
    default:
      fprintf(stderr, "Unrecognized available weight option\n"); exit(-1);
  }
}


/** @brief Transfers weight from a neighboring satisfied clause to an
 *         unsatisfied clause when no more cost-reducing variables exist.
 *
 *  @param from_idx  The clause index that weight is taken from.
 *  @param to_idx    The clause index that weight is given to.
 *  @param w         The amount of weight to transfer. Note that if not enough
 *                   weight exists for a transfer, then the transfer does not
 *                   take place.
 */
static void transfer_weight(int from_idx, int to_idx, double w) {
  // Abort if there isn't enough weight for a transfer
  if (clause_weights[from_idx] < w) return;

  clause_weights[from_idx] -= w;
  clause_weights[to_idx] += w;
  unsat_clause_weight += w;

  // Store the change in weight in the transfer memory
  num_weight_transfers++;
  recent_transfer_weights[recent_transfer_idx++] = w;
  if (recent_transfer_idx == WEIGHT_TRANSFER_MEMORY) {
    recent_transfer_idx = 0;
  }

  // If the satisfied "from" literal is critical, then a change in weight could
  // update its cost-reducing status
  if (clause_num_true_lits[from_idx] == 1) {
    const int mask_lit = clause_lit_masks[from_idx];
    add_cost_compute_var(VAR_IDX(mask_lit));

    // Since critical, update the critical weight for the remaining var
    literal_crit_sat_weights[mask_lit] -= w;
  }

  // All literals in the "to" clause can be cost-reducing, add to cost compute
  const int to_size = clause_sizes[to_idx];
  int *to_lits = clause_literals[to_idx];
  for (int l = 0; l < to_size; l++) {
    const int l_idx = *to_lits;
    add_cost_compute_var(VAR_IDX(l_idx));

    // Also add weight to the unsat weights since more weight on the clause
    literal_unsat_weights[l_idx] += w;

    to_lits++;
  }

  // Since weight has been taken from a sat clause, update its neighborhood
  update_neighborhood_on_weight_transfer(from_idx, w);
}


/** @brief Distributes weight from a singular satisfied to a specified 
 *         unsatisfied clause.
 *
 *  @param c_idx  The index number of the unsatisfied clause to give weight to.
 */
static void distribute_singular(const int c_idx) {
  const int neigh_size = neigh_sizes[c_idx];
  int *neighborhood = neigh_clauses[c_idx];

  // Goal is to find a neighbor with maximum weight
  int max_idx = -1;
  double max_weight = -1; // Positive weights required, so OK as min
  int is_in_neighborhood = 0;

  int second_max_idx = -1;
  double second_max_weight = -1; // Used as an optimization technique

  // Check max store first
  if (neigh_max_idxs[c_idx] != -1) {
    // There is a maximum neighbor - use it and clear it
    max_idx = neigh_max_idxs[c_idx];
    max_weight = neigh_max_weights[c_idx];
    neigh_max_idxs[c_idx] = -1;
  } else {
    // Scan neighbor weights for max weight
    for (int c = 0; c < neigh_size; c++) {
      const int cn_idx = *neighborhood;
      if (clause_num_true_lits[cn_idx] > 0) {
        double w = clause_weights[cn_idx];
        if (w > max_weight) {
          second_max_idx = max_idx;
          second_max_weight = max_weight;
          max_idx = cn_idx;
          max_weight = w;
        } else if (w > second_max_weight) {
          second_max_idx = cn_idx;
          second_max_weight = w;
        }
      }

      neighborhood++;
    }
  }

  /* While not mentioned in the original DDFW paper, it was found in the
   * original authors' code: with small probability, toss the max weight
   * neighbor found and choose a random satisfied neighbor instead.
   */
  if (max_idx != -1) {
    if (max_weight < init_clause_weight || ((unsigned int) rand()) 
        % RAND_NEIGH_DENOMINATOR < RAND_NEIGH_NUMERATOR) {
      max_idx = -1;
    } else {
      is_in_neighborhood = 1;
    }
  }

  if (max_idx == -1) {
    max_idx = get_random_satisfied_clause();
  }

  // Transfer the weight, then update the max if in neighborhood
  max_weight = clause_weights[max_idx];
  double w = calculate_transfer_weight(max_weight);

  if (is_in_neighborhood) {
    if (second_max_idx != -1 && second_max_weight > clause_weights[max_idx]) {
      neigh_max_weights[c_idx] = second_max_weight;
      neigh_max_idxs[c_idx] = second_max_idx;
    } else {
      neigh_max_weights[c_idx] = clause_weights[max_idx];
      neigh_max_idxs[c_idx] = max_idx;
    }
  }

  transfer_weight(max_idx, c_idx, w);
}


/** @brief Distributes weight from all clauses with weight above the initial
 *         weight.
 *
 *  @param c_idx  The index number of the clause to transfer weight to.
 */
static void distribute_above_init(const int c_idx) {
  // If the randomness says so, take weight from random neighbor instead
  if (((unsigned int) rand()) % RAND_NEIGH_DENOMINATOR < RAND_NEIGH_NUMERATOR) {
    int idx = get_random_satisfied_clause();
    double cw = clause_weights[idx];
    double w = calculate_transfer_weight(cw);
    transfer_weight(idx, c_idx, w);
    return;
  }

  const int neigh_size = neigh_sizes[c_idx];
  int *neighborhood = neigh_clauses[c_idx];

  if (neigh_size == 0) return; // Strange short-circuit, just in case
  if (neigh_num_sat[c_idx] == 0) return;

  // If the transfer rule is SUM or AVG, then a total sum is needed
  double sum = neigh_weights[c_idx];
  double avg_weight = sum / ((double) neigh_size);
  double avg = calculate_transfer_weight(avg_weight); // TODO no?

  // Loop over all neighbors and transfer if weight above init
  for (int c = 0; c < neigh_size; c++) {
    const int cn_idx = *neighborhood;
    const double w = clause_weights[cn_idx];
    if (clause_num_true_lits[cn_idx] > 0 && w >= init_clause_weight) {
      // Apply the rule, one way or another
      switch (rule_grp) {
        case RULE_SINGULAR:
          transfer_weight(cn_idx, c_idx, calculate_transfer_weight(w)); break;
        case RULE_SUM:
          transfer_weight(cn_idx, c_idx, (w / sum) * avg);              break;
        case RULE_AVG:
          transfer_weight(cn_idx, c_idx, avg);                          break;
        default:
          fprintf(stderr, "Unrecognized rule group\n"); exit(-1);
      }
    }
    neighborhood++;
  }
}


/** @brief Distributes weight from all satisfied clauses it can.
 *
 *  @param c_idx  The index number of the clause to transfer weight to.
 */
static void distribute_all(const int c_idx) {
  
  const int neigh_size = neigh_sizes[c_idx];
  int *neighborhood = neigh_clauses[c_idx];

  if (neigh_size == 0) return; // Short-circuit to prevent div by 0
  if (neigh_num_sat[c_idx] == 0) return;

  // If the transfer rule is SUM or AVG, then a total sum is needed
  double sum = neigh_weights[c_idx];
  double avg_weight = sum / ((double) neigh_size);
  double avg = calculate_transfer_weight(avg_weight); // TODO no?

  // Loop over all neighbors and transfer weight if satisfied
  for (int c = 0; c < neigh_size; c++) {
    const int cn_idx = *neighborhood;
    const double w = clause_weights[cn_idx];
    if (clause_num_true_lits[cn_idx] > 0) {
      // If the randomness says so, take weight from random neighbor instead
      if (((unsigned int) rand()) % RAND_NEIGH_DENOMINATOR < RAND_NEIGH_NUMERATOR) {
        int idx = get_random_satisfied_clause();
        double cw = clause_weights[idx];
        double w = calculate_transfer_weight(cw);
        transfer_weight(idx, c_idx, w);
        continue;
      }

      switch (rule_grp) {
        case RULE_SINGULAR:
          transfer_weight(cn_idx, c_idx, calculate_transfer_weight(w)); break;
        case RULE_SUM:
          transfer_weight(cn_idx, c_idx, (w / sum) * avg);              break;
        case RULE_AVG:
          transfer_weight(cn_idx, c_idx, avg);                          break;
        default:
          fprintf(stderr, "Unrecognized rule group\n"); exit(-1);
      }
    }

    neighborhood++;
  }
}


/** @brief Applies the rule to the sum of satisfied weight.
 *
 *  @param c_idx  The index number of the clause to transfer weight to.
 */
static void distribute_active(const int c_idx) {
  const int neigh_size = neigh_sizes[c_idx];
  int *neighborhood = neigh_clauses[c_idx];

  if (neigh_size == 0) return; // Short-circuit to prevent div by 0
  if (neigh_num_sat[c_idx] == 0) return;

  // If the transfer rule is SUM or AVG, then a total sum is needed
  double sum = neigh_weights[c_idx];
  double weight_diff = calculate_transfer_weight(sum);
  double avg = weight_diff / ((double) neigh_num_sat[c_idx]);

  // Loop over all neighbors and transfer weight if satisfied
  for (int c = 0; c < neigh_size; c++) {
    const int cn_idx = *neighborhood;
    const double w = clause_weights[cn_idx];
    if (clause_num_true_lits[cn_idx] > 0) {
      switch (rule_grp) {
        case RULE_SUM:
          transfer_weight(cn_idx, c_idx, (w / sum) * avg);              break;
        case RULE_AVG:
          transfer_weight(cn_idx, c_idx, avg);                          break;
        default:
          fprintf(stderr, "Unrecognized rule group\n"); exit(-1);
      }
    }

    neighborhood++;
  }
}



/** @brief Distribute weights from satisfied to unsatisfied clauses.
 *
 *  In the case where no literals may be flipped to decrease the weight
 *  held by the unsatisfied clauses, the weight will be re-distributed from
 *  satisfied to unsatisfied clauses according to the settings specified
 *  at the command line.
 */
void distribute_weights(void) {
  // Loop over all false clauses
  int *false_clauses = false_clause_members;
  for (int c = 0; c < num_unsat_clauses; c++) {
    const int c_idx = *false_clauses;
    switch (transfer_grp) {
      case SINGULAR:
        distribute_singular(c_idx);      break;
      case ABOVE_INIT:
        distribute_above_init(c_idx);    break;
      case ALL:
        distribute_all(c_idx);           break;
      case ACTIVE:
        distribute_active(c_idx);        break;
      default:
        fprintf(stderr, "Unrecognized transfer group\n"); exit(-1);
    }

    false_clauses++;
  }
}


/** @brief Get the average of the most recent weight transfer amounts.
 *
 *  Calculates the average of the last ~1000 weight transfer amounts.
 */
double get_transfer_average(void) {
  double weight_sum = 0.0;
  for (int i = 0; i < WEIGHT_TRANSFER_MEMORY; i++) {
    weight_sum += recent_transfer_weights[i];
  }

  return weight_sum / ((double) WEIGHT_TRANSFER_MEMORY);
}
