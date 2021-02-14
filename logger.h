/** @file logger.h
 *  @brief Implementations of logging functions.
 *
 *  See logger.c for implementation details.
 *
 *  @author Cayden Codel (ccodel@andrew.cmu.edu)
 *
 *  @bug No known bugs.
 *  @TODO Think about this file's name - printer or logger?
 */

#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <sys/time.h>

/** @brief Defines the verbosity level for the algorithm.
 *
 *  A single verbosity level is maintained by the logger. Set the verbosity
 *  level with a call to set_verbosity().
 *
 *  The three levels are:
 *
 *    SILENT:  Logging functions make no output.
 *    NORMAL:  Logging functions print some amount of information.
 *    VERBOSE: Logging functions print plenty of information.
 */
typedef enum verbosity_level {
  SILENT, NORMAL, VERBOSE
} verbose_t;

/** -- Printing functions: print to stdout, regardless of verbosity. -- */

/** Print help information when the -h flag is used. */
void print_help(char *runtime_path);

/** Print usage information when the provided arguments aren't well-formed. */
void print_usage(char *runtime_path);

/** -- Logging functions: print to stdout, according to verbosity level. -- */

verbose_t get_verbosity(void);
void set_verbosity(verbose_t verbosity_level);

void log_clause(int c_idx);
void log_literal(int l_idx);
void log_str(const char *format, ...);
void log_err(const char *format, ...);

/** Loop information */
void log_weights(void);
void log_reducing_cost_lits(void);

void log_weight_statistics(int run, int transfers, double weight_transfer_avg);
void log_statistics(int run, struct timeval *start, struct timeval *stop);

void log_assignment(void);
void output_assignment(void);

#endif /* _LOGGER_H_ */
