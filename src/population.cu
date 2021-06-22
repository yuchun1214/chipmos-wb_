#include <include/population.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

#define cudamalloc_objects(type, dvc_objects, address_of_objects, sample, amount_of_objects, amount_of_chromosomes, obj_name)\
    cudaCheck(cudaMallocHost((void**)&address_of_objects, sizeof(type*)*amount_of_chromosomes), "cudaMallocHost for address_of_address" obj_name);\
    for(int i = 0; i < amount_of_chromosomes; ++i){\
        type * entry;\
        cudaCheck(cudaMalloc((void**)&entry, sizeof(type)*amount_of_objects), "cudaMalloc for entry of " obj_name);\
        cudaCheck(cudaMemcpy(entry, sample, sizeof(type)*amount_of_objects, cudaMemcpyHostToDevice), "cudaMemcpy entry of " obj_name "from host to device");\
        address_of_objects[i] = entry;\
    }\
    cudaCheck(cudaMalloc((void**)&dvc_objects, sizeof(type*)*amount_of_chromosomes), "cudaMalloc for " obj_name);\
    cudaCheck(cudaMemcpy(dvc_objects, address_of_objects, sizeof(type*)*amount_of_chromosomes, cudaMemcpyHostToDevice), "cudaMemcpy" obj_name "from host to device");



// void __global__ initJobs(job_t ** jobs, process_time_t ** process_times, 
    
void initializePopulation(population_t * pop){
    
    int AMOUNT_OF_R_CHROMOSOMES = pop->parameters.AMOUNT_OF_R_CHROMOSOMES;
    int AMOUNT_OF_JOBS = pop->task.AMOUNT_OF_JOBS;
    int AMOUNT_OF_MACHINES = pop->task.AMOUNT_OF_MACHINES;

    // malloc for jobs;
    job_t ** jobs;
    job_t ** address_of_jobs;
    // job_t * job_entry;
    cudamalloc_objects(job_t, jobs, address_of_jobs, pop->task.jobs, AMOUNT_OF_JOBS, AMOUNT_OF_R_CHROMOSOMES, "jobs"); 
    // cudaCheck(cudaMallocHost((void**)&address_of_jobs, sizeof(job_t*)*AMOUNT_OF_R_CHROMOSOMES), "cudaMallocHost for address_of_jobs");
    // // malloc for entry
    // for(int i = 0; i < AMOUNT_OF_R_CHROMOSOMES; ++i){
    //     cudaCheck(cudaMalloc((void**)&job_entry, sizeof(job_t)*AMOUNT_OF_JOBS), "cudaMalloc for entry of jobs");
    //     cudaCheck(cudaMemcpy(job_entry, pop->task.jobs, AMOUNT_OF_JOBS, cudaMemcpyHostToDevice), "cudaMemcpy entry of jobs from host to device");
    //     address_of_jobs[i] = job_entry;
    // }
    // cudaCheck(cudaMalloc((void**)&jobs, sizeof(job_t *)*AMOUNT_OF_R_CHROMOSOMES), "cudaMalloc for jobs");
    // cudaCheck(cudaMemcpy(jobs, address_of_jobs, sizeof(job_t*)*AMOUNT_OF_R_CHROMOSOMES, cudaMemcpyHostToDevice), "cudaMemcpy jobs from host to device");

    pop->device_objects.jobs = jobs;
    pop->host_objects.address_of_cujobs = address_of_jobs; 
    
    // malloc for machines
    machine_t **machines;
    machine_t **address_of_machines;
    // machine_t *machine_entry;
    cudamalloc_objects(machine_t,  machines, address_of_machines, pop->task.machines, AMOUNT_OF_MACHINES, AMOUNT_OF_R_CHROMOSOMES, "machines");
    // cudaCheck(cudaMallocHost((void**)&address_of_machines, sizeof(machine_t *)*AMOUNT_OF_R_CHROMOSOMES), "cudaMallocHost for address_of_machines"); 
    // for(int i = 0; i < AMOUNT_OF_R_CHROMOSOMES; ++i){
    //     cudaCheck(cudaMalloc((void**)&machine_entry, sizeof(machine_t)*AMOUNT_OF_MACHINES), "cudamalloc for entry of machines");
    //     cudaCheck(cudaMemcpy(machine_entry, pop->task.machines, sizeof(machine_t)*AMOUNT_OF_MACHINES, cudaMemcpyHostToDevice), "cudaMemcpy entry of machines from host to device");
    //     address_of_machines[i] = machine_entry;
    // }
    // cudaCheck(cudaMalloc((void**)&machines, sizeof(machine_t*)*AMOUNT_OF_R_CHROMOSOMES), "cudaMalloc for machines");
    // cudaCheck(cudaMemcpy(machines, address_of_machines, sizeof(machine_t*)*AMOUNT_OF_R_CHROMOSOMES, cudaMemcpyHostToDevice), "cudaMemcpy for machines from hst to dvc");
    
    pop->device_objects.machines = machines;
    pop->host_objects.address_of_cumachines = address_of_machines;

    // malloc for tools and wires;
    tool_t ** tools;
    tool_t ** address_of_tools;
    cudamalloc_objects(tool_t, tools, address_of_tools, pop->task.tools, AMOUNT_OF_MACHINES, AMOUNT_OF_R_CHROMOSOMES, "tools");

    wire_t ** wires;
    wire_t ** address_of_wires;
    cudamalloc_objects(wire_t, wires, address_of_wires, pop->task.wires, AMOUNT_OF_MACHINES, AMOUNT_OF_R_CHROMOSOMES, "wires");
}


void geneticAlgorithm(population_t *pop){

}
