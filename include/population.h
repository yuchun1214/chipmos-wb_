#ifndef __POPULATION_H__
#define __POPULATION_H__
#include <vector>
#include "include/lot.h"
#include "include/machine_base.h"
#include <include/common.h>
#include <include/job.h>
#include <include/job_base.h>
#include <include/machine.h>
#include <include/chromosome_base.h>

typedef struct task_t{
    int task_id;
    int AMOUNT_OF_JOBS;
    int AMOUNT_OF_MACHINES;
    // sample
    job_t * jobs; 
    machine_t * machines;
    tool_t * tools;
    wire_t * wires;
    process_time_t  ** process_times;
    int * size_of_process_times;
}task_t;

typedef struct evolution_factors_t{
    unsigned int *c_selected1;
    unsigned int *c_selected2;
    unsigned int *cut_points;
    unsigned int *range;

    unsigned int *m_selected;
    unsigned int *gene_idx;
    double *new_genes;
}evolution_factors;


struct population_t{
    unsigned int no;
    struct {
        int AMOUNT_OF_CHROMOSOMES;
        int AMOUNT_OF_R_CHROMOSOMES;
        double EVOLUTION_RATE;
        double SELECTION_RATE;
        int GENERATIONS;
        int SWAP_CHROMOSOMES;
    }parameters;

    task_t task;

    struct {
        job_t ** jobs;
        machine_t ** machines;
        tool_t ** tools;
        wire_t ** wires;
        process_time_t ** process_times;
        int * size_of_process_times;

        list_operations_t * list_ops;
        job_base_operations_t * job_ops;
        machine_base_operations_t * machine_ops;
    } device_objects;

    struct {
        job_t ** address_of_cujobs;
        machine_t ** address_of_cumachines;
        tool_t ** address_of_tools;
        wire_t ** address_of_wires;
        process_time_t ** address_of_process_times_entry;
    }host_objects;

    struct {
        chromosome_base_t * chromosomes;
        double ** genes;
        double ** address_of_cugenes;
        
        chromosome_base_t * host_chromosomes;
    }chromosomes;



    struct{
        evolution_factors_t device;
        evolution_factors_t host;
    }evolution_factors;

};


void initializePopulation(population_t * pop);

void cpyResult(population_t * pop, char *fileanme);

void *geneticAlgorithm(void *pop);


void swapPopulation(population_t pops[], const int AMOUNT_OF_POPULATIONS);

#endif
