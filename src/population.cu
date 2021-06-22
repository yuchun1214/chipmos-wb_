#include <include/population.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

void initializePopulation(population_t * pop){
    
    int AMOUNT_OF_R_CHROMOSOMES = pop->parameters.AMOUNT_OF_R_CHROMOSOMES;
    int AMOUNT_OF_JOBS = pop->task.AMOUNT_OF_JOBS;
    int AMOUNT_OF_MACHINES = pop->task.AMOUNT_OF_MACHINES;

    job_t ** jobs;
    job_t ** address_of_jobs;
    job_t * job_entry;
    
    // malloc for jobs;
    cudaCheck(cudaMallocHost((void**)&address_of_jobs, sizeof(job_t*)*AMOUNT_OF_R_CHROMOSOMES), "cudaMallocHost for address_of_jobs");
    // malloc for entry
    for(int i = 0; i < AMOUNT_OF_R_CHROMOSOMES; ++i){
        cudaCheck(cudaMalloc((void**)&job_entry, sizeof(job_t)*AMOUNT_OF_JOBS), "cudaMalloc for entry of jobs");
        cudaCheck(cudaMemcpy(job_entry, pop->task.jobs, AMOUNT_OF_JOBS, cudaMemcpyHostToDevice), "cudaMemcpy entry of jobs from host to device");
        address_of_jobs[i] = job_entry;
    }
    cudaCheck(cudaMalloc((void**)&jobs, sizeof(job_t *)*AMOUNT_OF_R_CHROMOSOMES), "cudaMalloc for jobs");
    cudaCheck(cudaMemcpy(jobs, address_of_jobs, sizeof(job_t*)*AMOUNT_OF_R_CHROMOSOMES, cudaMemcpyHostToDevice), "cudaMemcpy jobs from host to device");
    
   
}

