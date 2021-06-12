#ifndef __POPULATION_H__
#define __POPULATION_H__
#include <vector>
#include "include/lot.h"
#include <include/common.h>
#include <include/job.h>
#include <include/job_base.h>
#include <include/machine.h>
#include <include/chromosome_base.h>

typedef struct round_t{
    int round_no;
    int AMOUNT_OF_JOBS;
    int AMOUNT_OF_MACHINES;
    job_t * jobs; // sample;
    // machine_t ** machines;
    std::map<unsigned int, machine_t *> machines;
    tool_t ** tools;
    wire_t ** wires;
    process_time_t  ** process_times;
    int * size_of_process_times;
     
}round_t;


struct population_t{
    unsigned int no;
    struct {
        int AMOUNT_OF_MACHINES;
        int AMOUNT_OF_CHROMOSOMES;
        int AMOUNT_OF_R_CHROMOSOMES;
        double EVOLUTION_RATE;
        double SELECTION_RATE;
        int GENERATIONS;
    }parameters;
    
    std::vector<std::vector<lot_group_t> > groups;

    unsigned int current_round_no;
    round_t round;

    chromosome_base_t * chromosomes; 

    // struct {
    //     job_t ** jobs;
    //     machine_t **machines;
    //     tool_t ** tools;
    //     wire_t ** wires;
    // }device_objects;

    // struct {
    //     job_t ** address_of_cujobs;
    //     machine_t ** address_of_cumachines;
    //     tool_t ** address_of_tools;
    //     wire_t ** address_of_wires;
    //     process_time_t **address_of_process_time_entry;
    // }host_objects;

};

void initializePopulation(population_t *pop);

void clearARound(population_t *pop, int round);

void algorithm(population_t *pop);

#endif
