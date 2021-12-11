#ifndef __POPULATION_H__
#define __POPULATION_H__
#include <include/chromosome_base.h>
#include <include/common.h>
#include <include/job.h>
#include <include/job_base.h>
#include <include/machine.h>
#include <vector>
#include "include/lot.h"
#include "include/machine_base.h"

typedef struct round_t {
    int round_no;
    int AMOUNT_OF_JOBS;
    int AMOUNT_OF_MACHINES;
    job_t *jobs;  // sample;
    // machine_t ** machines;
    std::map<unsigned int, machine_t *> machines;
    tool_t **tools;
    wire_t **wires;
    process_time_t **process_times;
    int *size_of_process_times;

} round_t;

typedef struct chromosomes_t {
    double velocity;
    chromosome_base_t current_chromosome;
    chromosome_base_t history_best_chromosome;
} chromosomes_t;

struct population_t {
    unsigned int no;
    struct {
        int AMOUNT_OF_CHROMOSOMES;
        int AMOUNT_OF_R_CHROMOSOMES;
        double EVOLUTION_RATE;
        double SELECTION_RATE;
        int GENERATIONS;
    } parameters;

    struct {
        list_operations_t *list_ops;
        job_base_operations_t *job_ops;
        machine_base_operations_t *machine_ops;
    } operations;

    std::vector<std::vector<lot_group_t> > groups;

    unsigned int current_round_no;
    round_t round;

    // chromosome_base_t * chromosomes;
    chromosomes_t *chromosomes;
    chromosome_base_t global_best_chromosome;
};

void initializePopulation(population_t *pop);

void clearARound(population_t *pop, int round);

void algorithm(population_t *pop);

#endif
