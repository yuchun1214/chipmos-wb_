#include <cuda.h>
#include <cuda_runtime_api.h>
#include <include/population.h>

void clearARound(population_t *pop, int _round) {}

void initializeARound(population_t *pop, int _round)
{
    // round_t round = pop->samples.rounds[_round];
    // int AMOUNT_OF_R_CHROMOSOMES = pop->parameters.AMOUNT_OF_R_CHROMOSOMES;
    // int AMOUNT_OF_JOBS = round.AMOUNT_OF_JOBS;

    // // prepare jobs
    // job_t * job_entry;
    // job_t ** jobs;
    // job_t ** address_of_job_entry;
    // cudaCheck(cudaMallocHost((void**)&address_of_job_entry,
    // sizeof(job_t*)*AMOUNT_OF_R_CHROMOSOMES), "cudaMallocHost for addresses of
    // entries of jobs"); cudaCheck(cudaMallocHost((void**)&job_entry,
    // sizeof(job_t)*AMOUNT_OF_JOBS), "cudaMallocHost for an entry of jobs");
    // for(int i = 0; i < AMOUNT_OF_R_CHROMOSOMES; ++i){
    //     job_t * temp;
    //     cudaCheck(cudaMalloc((void**)&temp, sizeof(job_t)*AMOUNT_OF_JOBS),
    //     "cudaMalloc for an entry of jobs"); cudaCheck(cudaMemcpy(temp,
    //     job_entry, sizeof(job_t)*AMOUNT_OF_JOBS, cudaMemcpyHostToDevice),
    //     "cudaMemcpy for an entry of jobs");
    // }
    // cudaCheck(cudaMalloc((void**)&jobs, sizeof(job_t
    // *)*AMOUNT_OF_R_CHROMOSOMES), "cudaMalloc for a round of jobs");
    // cudaCheck(cudaMemcpy(jobs, address_of_job_entry, sizeof(job_t
    // *)*AMOUNT_OF_R_CHROMOSOMES, cudaMemcpyHostToDevice), "cudaMemcpy for
    // address of job entry from host to device");

    // prepare tool
}

void initializePopulation(population_t *pop) {}

void algorithm(population_t *pop) {}
