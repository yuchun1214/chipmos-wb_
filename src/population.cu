#include <cuda.h>
#include <cuda_runtime_api.h>
#include <ctime>
#include <curand.h>
#include <curand_kernel.h>

#include "include/population.h"
#include "include/chromosome_base.h"
#include "include/entity.h"
#include "include/infra.h"
#include "include/job.h"
#include "include/job_base.h"
#include "include/linked_list.h"
#include "include/machine_base.h"

#define cudamalloc_objects(type, dvc_objects, address_of_objects, sample,      \
                           amount_of_objects, amount_of_chromosomes, obj_name) \
    cudaCheck(cudaMallocHost((void **) &address_of_objects,                    \
                             sizeof(type *) * amount_of_chromosomes),          \
              "cudaMallocHost for address_of_address" obj_name);               \
    for (int i = 0; i < amount_of_chromosomes; ++i) {                          \
        type *entry;                                                           \
        cudaCheck(                                                             \
            cudaMalloc((void **) &entry, sizeof(type) * amount_of_objects),    \
            "cudaMalloc for entry of " obj_name);                              \
        cudaCheck(cudaMemcpy(entry, sample, sizeof(type) * amount_of_objects,  \
                             cudaMemcpyHostToDevice),                          \
                  "cudaMemcpy entry of " obj_name "from host to device");      \
        address_of_objects[i] = entry;                                         \
    }                                                                          \
    cudaCheck(cudaMalloc((void **) &dvc_objects,                               \
                         sizeof(type *) * amount_of_chromosomes),              \
              "cudaMalloc for " obj_name);                                     \
    cudaCheck(cudaMemcpy(dvc_objects, address_of_objects,                      \
                         sizeof(type *) * amount_of_chromosomes,               \
                         cudaMemcpyHostToDevice),                              \
              "cudaMemcpy" obj_name "from host to device");


#define device_malloc_evolution_factor(factor, type, amount_of_chromosomes) \
    cudaCheck(cudaMalloc((void **) &(factor),                               \
                         sizeof(type) * (amount_of_chromosomes)),           \
              "cudaMalloc for " #factor);

#define host_malloc_evolution_factor(factor, type, amount_of_chromosomes) \
    cudaCheck(cudaMallocHost((void **) &(factor),                         \
                             sizeof(type) * (amount_of_chromosomes)),     \
              "cudaMalloc for " #factor);

#define cpy_factor_h2d(dest, src, type, amount_of_chromosomes)              \
    cudaCheck(cudaMemcpy(dest, src, sizeof(type) * (amount_of_chromosomes), \
                         cudaMemcpyHostToDevice),                           \
              "cudaMemcpy for " #src " from host to device");

__global__ void setup_rand_state(curandState *state, unsigned long long *seed){
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    curand_init(seed[idx], idx, 0, &state[idx]);
}

__device__ unsigned int generate_random_uint(curandState * state, int lower, int upper, int difference){
    float f_rnd;
    int i_rnd;
    do{
        f_rnd = curand_uniform(state);
        i_rnd = f_rnd * (upper - lower) + lower;
    }while(difference > 0 && difference == i_rnd);
    
    return i_rnd;
}


__global__ void generateCrossoverFactors(curandState *states, unsigned int *c_selected1, unsigned int *c_selected2, unsigned int *cut_points, unsigned int *range, int GENE_SIZE, int FACTOR_SIZE, int AMOUNT_OF_CHROMOSOMES){
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if(idx < FACTOR_SIZE){
        unsigned int r1, r2, tmp;

        r1 = generate_random_uint(&states[idx], 0, AMOUNT_OF_CHROMOSOMES, -1);
        r2 = generate_random_uint(&states[idx], 0, AMOUNT_OF_CHROMOSOMES, r1); 
        c_selected1[idx] = r1;
        c_selected2[idx] = r2;

        tmp = cut_points[idx] = generate_random_uint(&states[idx], 0, GENE_SIZE, -1);
        range[idx] = generate_random_uint(&states[idx], 0, GENE_SIZE - tmp, -1);

    }
}

__global__ void generateMutationFactors(curandState *states, unsigned int *m_selected, unsigned int *gene_idx, double *ngenes, int GENE_SIZE, int FACTOR_SIZE, int AMOUNT_OF_CHROMOSOMES){
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    if(idx < FACTOR_SIZE){
        m_selected[idx] = generate_random_uint(&states[idx], 0, AMOUNT_OF_CHROMOSOMES, -1);
        gene_idx[idx] = generate_random_uint(&states[idx], 0, GENE_SIZE, -1);
        ngenes[idx] = curand_uniform_double(&states[idx]);
    }
}


__global__ void initOps(list_operations_t *list_ops,
                        job_base_operations_t *job_ops,
                        machine_base_operations_t *machine_ops)
{
    *list_ops = LINKED_LIST_OPS;
    *job_ops = JOB_BASE_OPS;
    *machine_ops = MACHINE_BASE_OPS;
    machine_ops->setup_times[0] = setup_time_CWN;
    machine_ops->setup_times[1] = setup_time_CK;
    machine_ops->setup_times[2] = setup_time_EU;
    machine_ops->setup_times[3] = setup_time_MC_SC;
    machine_ops->setup_times[4] = setup_time_CSC;
    machine_ops->setup_times[5] = setup_time_USC;
    machine_ops->sizeof_setup_time_function_array = 6;
    machine_ops->reset = machine_reset;
}


__global__ void initializeJobs(job_t **jobs,
                               process_time_t **process_times,
                               job_base_operations_t *job_ops,
                               chromosome_base_t *chromosomes,
                               int AMOUNT_OF_JOBS,
                               int AMOUNT_OF_R_CHROMOSOMES)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x < AMOUNT_OF_R_CHROMOSOMES && y < AMOUNT_OF_JOBS) {
        job_initialize(&jobs[x][y]);
        job_ops->set_process_time(&jobs[x][y].base, process_times[y],
                                  jobs[x][y].base.size_of_process_time);
        job_ops->set_ms_gene_addr(&jobs[x][y].base,
                                  chromosomes[x].ms_genes + y);
        job_ops->set_os_gene_addr(&jobs[x][y].base,
                                  chromosomes[x].os_genes + y);
    }
}

__global__ void binding(job_t **jobs,
                        chromosome_base_t *chromosomes,
                        job_base_operations_t *ops,
                        const int AMOUNT_OF_JOBS,
                        const int AMOUNT_OF_CHROMOSOMES)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x < AMOUNT_OF_CHROMOSOMES && y < AMOUNT_OF_JOBS) {
        ops->set_ms_gene_addr(&jobs[x][y].base, chromosomes[x].ms_genes + y);
        ops->set_os_gene_addr(&jobs[x][y].base, chromosomes[x].os_genes + y);
    }
}

__global__ void initializeMachines(machine_t **machines,
                                   tool_t **tools,
                                   wire_t **wires,
                                   const int AMOUNT_OF_MACHINES,
                                   const int AMOUNT_OF_CHROMOSOMES)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    if (x < AMOUNT_OF_CHROMOSOMES && y < AMOUNT_OF_MACHINES) {
        initMachine(&machines[x][y]);
        machines[x][y].tool = &tools[x][y];
        machines[x][y].wire = &wires[x][y];
    }
}

__global__ void initializeChromosomes(chromosome_base_t *chromosomes,
                                      double **genes,
                                      const int AMOUNT_OF_JOBS,
                                      const int AMOUNT_OF_CHROMOSOMES)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if (x < AMOUNT_OF_CHROMOSOMES) {
        chromosomes[x].chromosome_no = x;
        chromosomes[x].fitnessValue = 0;
        chromosomes[x].gene_size = AMOUNT_OF_JOBS << 1;
        chromosome_base_init(chromosomes + x, genes[x]);
    }
}

__global__ void initializeChromosomes(chromosome_base_t *chromosomes,
                                      const int AMOUNT_OF_CHROMOSOMES)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if (x < AMOUNT_OF_CHROMOSOMES) {
        chromosomes[x].chromosome_no = x;
        chromosomes[x].fitnessValue = 0;
        chromosome_base_init(&chromosomes[x], chromosomes[x].genes);
    }
}

__global__ void resetMachines(machine_t **machines,
                              machine_base_operations_t *ops,
                              int AMOUNT_OF_MACHINES,
                              int AMOUNT_OF_CHROMOSOMES)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x < AMOUNT_OF_CHROMOSOMES && y < AMOUNT_OF_MACHINES) {
        ops->reset(&machines[x][y].base);
        machines[x][y].makespan = 0;
    }
}

__global__ void machineSelection(job_t **jobs,
                                 job_base_operations_t *jbops,
                                 int AMOUNT_OF_JOBS,
                                 int AMOUNT_OF_R_CHROMOSOMES)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x < AMOUNT_OF_R_CHROMOSOMES && y < AMOUNT_OF_JOBS) {
        int machine_idx = jbops->machine_selection(&jobs[x][y].base);
        jobs[x][y].base.machine_no =
            jobs[x][y].base.process_time[machine_idx].machine_no;
        jobs[x][y].base.ptime =
            jobs[x][y].base.process_time[machine_idx].process_time;
    }
}

__global__ void machineSelection2(job_t **jobs,
                                  machine_t **machines,
                                  machine_base_operations_t *ops,
                                  int AMOUNT_OF_JOBS,
                                  int AMOUNT_OF_MACHINES,
                                  int AMOUNT_OF_R_CHROMOSOMES)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x < AMOUNT_OF_R_CHROMOSOMES && y < AMOUNT_OF_MACHINES) {
        for (int i = 0; i < AMOUNT_OF_JOBS; ++i) {
            if (jobs[x][i].base.machine_no == machines[x][y].base.machine_no) {
                ops->add_job(&machines[x][y].base, &jobs[x][i].list);
            }
        }
    }
}

__global__ void getMachineInformation(machine_t *machines,
                                      unsigned int *machine_numbers,
                                      unsigned int *sizeof_jobs,
                                      int AMOUNT_OF_MACHINES)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if (x < AMOUNT_OF_MACHINES) {
        machine_numbers[x] = machines[x].base.machine_no;
        sizeof_jobs[x] = machines[x].base.size_of_jobs;
    }
}

__global__ void getMachineJobs(machine_t *machines,
                               double **data,
                               int AMOUNT_OF_MACHINES)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if (x < AMOUNT_OF_MACHINES) {
        list_ele_t *list = machines[x].base.root;
        job_t *job;
        int i = 0;
        while (list) {
            job = (job_t *) list->ptr_derived_object;
            data[x][i] = job->base.start_time;
            list = list->next;
            ++i;
        }
    }
}

__global__ void sortJob(machine_t **machines,
                        list_operations_t *list_ops,
                        machine_base_operations_t *mbops,
                        int AMOUNT_OF_MACHINES,
                        int AMOUNT_OF_R_CHROMOSOMES)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x < AMOUNT_OF_R_CHROMOSOMES && y < AMOUNT_OF_MACHINES) {
        mbops->sort_job(&machines[x][y].base, list_ops);
    }
}

__global__ void scheduling(machine_t **machines,
                           job_base_operations_t *job_ops,
                           machine_base_operations_t *machine_ops,
                           int AMOUNT_OF_MACHINES,
                           int AMOUNT_OF_R_CHROMOSOMES)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x < AMOUNT_OF_R_CHROMOSOMES && y < AMOUNT_OF_MACHINES) {
        scheduling(&machines[x][y], job_ops, machine_ops);
    }
}

__global__ void computeFitnessValue(machine_t **machines,
                                    chromosome_base_t *chromosomes,
                                    int AMOUNT_OF_MACHINES,
                                    int AMOUNT_OF_R_CHROMOSOMES)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if (x < AMOUNT_OF_R_CHROMOSOMES) {
        double total_completion_time = 0;
        for (int i = 0; i < AMOUNT_OF_MACHINES; ++i) {
            // if(machines[x][i].total_completion_time > worst){
            total_completion_time += machines[x][i].total_completion_time;
            // }
        }
        chromosomes[x].fitnessValue = total_completion_time;
    }
}

__global__ void sortChromosomes(chromosome_base_t *chromosomes,
                                int AMOUNT_OF_CHROMOSOMES)
{
    chromosome_base_t temp;
    for (int i = 0; i < AMOUNT_OF_CHROMOSOMES - 1; ++i) {
        for (int j = 0; j < AMOUNT_OF_CHROMOSOMES - 1; ++j) {
            if (chromosomes[j].fitnessValue > chromosomes[j + 1].fitnessValue) {
                temp = chromosomes[j];
                chromosomes[j] = chromosomes[j + 1];
                chromosomes[j + 1] = temp;
            }
        }
    }
}

__global__ void crossover(chromosome_base_t *chromosomes,
                          unsigned int *selected1,
                          unsigned int *selected2,
                          unsigned int *cut_points,
                          unsigned int *range,
                          unsigned int offset,
                          const int AMOUNT_OF_JOBS,
                          const int AMOUNT_OF_FACTORS,
                          const int AMOUNT_OF_CHROMOSOMES)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    double *gene1, *gene2, *r_gene1, *r_gene2;
    if (x < AMOUNT_OF_FACTORS) {
        gene1 = chromosomes[selected1[x]].genes;
        gene2 = chromosomes[selected2[x]].genes;
        r_gene1 = chromosomes[2 * x + offset + AMOUNT_OF_CHROMOSOMES].genes;
        r_gene2 = chromosomes[2 * x + offset + AMOUNT_OF_CHROMOSOMES + 1].genes;
        memcpy(r_gene1, gene1, sizeof(double) * AMOUNT_OF_JOBS * 2);
        memcpy(r_gene2, gene2, sizeof(double) * AMOUNT_OF_JOBS * 2);
        memcpy(r_gene1 + cut_points[x], gene2 + cut_points[x],
               sizeof(double) * range[x]);
        memcpy(r_gene2 + cut_points[x], gene1 + cut_points[x],
               sizeof(double) * range[x]);
    }
}

__global__ void mutation(chromosome_base_t *chromosomes,
                         unsigned int *selected,
                         unsigned int *gene_idx,
                         double *ngenes,
                         unsigned int offset,
                         const int AMOUNT_OF_JOBS,
                         const int AMOUNT_OF_MUTATIONS,
                         const int AMOUNT_OF_CHROMOSOMES)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    double *gene, *r_gene;
    if (x < AMOUNT_OF_MUTATIONS) {
        gene = chromosomes[selected[x]].genes;
        r_gene = chromosomes[x + offset + AMOUNT_OF_CHROMOSOMES].genes;
        memcpy(r_gene, gene, sizeof(double) * AMOUNT_OF_JOBS * 2);
        r_gene[gene_idx[x]] = ngenes[x];
    }
}

void generateCrossoverFactors(struct evolution_factors_t *factors,
                              int factor_size,
                              int gene_size,
                              const int AMOUNT_OF_R_CHROMOSOMES)
{
    int tmp;
    for (int i = 0; i < factor_size; ++i) {
        tmp = factors->c_selected1[i] =
            random_range(0, AMOUNT_OF_R_CHROMOSOMES, -1);
        factors->c_selected2[i] = random_range(0, AMOUNT_OF_R_CHROMOSOMES, tmp);
        tmp = factors->cut_points[i] = random_range(0, gene_size, -1);
        factors->range[i] = random_range(0, gene_size - tmp, -1);
    }
}

void generateMutationFactors(struct evolution_factors_t *factors,
                             int factor_size,
                             int gene_size,
                             const int AMOUNT_OF_R_CHROMOSOMES)
{
    for (int i = 0; i < factor_size; ++i) {
        factors->m_selected[i] = random_range(0, AMOUNT_OF_R_CHROMOSOMES, -1);
        factors->gene_idx[i] = random_range(0, gene_size, -1);
    }
    random(factors->new_genes, factor_size);
}

void cpyEvolutionFactors(struct evolution_factors_t *dest,
                         struct evolution_factors_t *src,
                         int AMOUNT_OF_CHROMOSOMES)
{
    cpy_factor_h2d(dest->c_selected1, src->c_selected1, unsigned int,
                   AMOUNT_OF_CHROMOSOMES);
    cpy_factor_h2d(dest->c_selected2, src->c_selected2, unsigned int,
                   AMOUNT_OF_CHROMOSOMES);
    cpy_factor_h2d(dest->m_selected, src->m_selected, unsigned int,
                   AMOUNT_OF_CHROMOSOMES);
    cpy_factor_h2d(dest->cut_points, src->cut_points, unsigned int,
                   AMOUNT_OF_CHROMOSOMES);
    cpy_factor_h2d(dest->range, src->range, unsigned int,
                   AMOUNT_OF_CHROMOSOMES);
    cpy_factor_h2d(dest->new_genes, src->new_genes, double,
                   AMOUNT_OF_CHROMOSOMES);
    cpy_factor_h2d(dest->gene_idx, src->gene_idx, unsigned int,
                   AMOUNT_OF_CHROMOSOMES);
}

void cpyResult(population_t *pop, char *filename)
{
    int AMOUNT_OF_JOBS = pop->task.AMOUNT_OF_JOBS;

    job_t *jobs;
    cudaCheck(cudaMallocHost((void **) &jobs, sizeof(job_t) * AMOUNT_OF_JOBS),
              "cudaMallocHost for jobs");
    cudaCheck(
        cudaMemcpy(jobs, pop->host_objects.address_of_cujobs[0],
                   sizeof(job_t) * AMOUNT_OF_JOBS, cudaMemcpyDeviceToHost),
        "cudaMemcpy jobs from d 2 h");

    FILE *file = fopen(filename, "w");
    for (int i = 0; i < AMOUNT_OF_JOBS; ++i) {
        fprintf(file, "%s, %s, %s, %s, %.3f, %.3f\n",
                jobs[i].base.job_info.data.text, jobs[i].part_no.data.text,
                jobs[i].part_id.data.text,
                convertUIntToEntityName(jobs[i].base.machine_no).c_str(),
                jobs[i].base.start_time, jobs[i].base.end_time);
    }
    cudaCheck(cudaFreeHost(jobs), "cudaFreeHost for jobs in cpyResult");
    return;
}

void initializePopulation(population_t *pop)
{
    int AMOUNT_OF_R_CHROMOSOMES = pop->parameters.AMOUNT_OF_R_CHROMOSOMES;
    int AMOUNT_OF_JOBS = pop->task.AMOUNT_OF_JOBS;
    int AMOUNT_OF_MACHINES = pop->task.AMOUNT_OF_MACHINES;

    /*******************************malloc for jobs;***************/
    job_t **jobs;
    job_t **address_of_jobs;
    // job_t * job_entry;
    cudamalloc_objects(job_t, jobs, address_of_jobs, pop->task.jobs,
                       AMOUNT_OF_JOBS, AMOUNT_OF_R_CHROMOSOMES, "jobs");
    // cudaCheck(cudaMallocHost((void**)&address_of_jobs,
    // sizeof(job_t*)*AMOUNT_OF_R_CHROMOSOMES), "cudaMallocHost for
    // address_of_jobs");
    // // malloc for entry
    // for(int i = 0; i < AMOUNT_OF_R_CHROMOSOMES; ++i){
    //     cudaCheck(cudaMalloc((void**)&job_entry,
    //     sizeof(job_t)*AMOUNT_OF_JOBS), "cudaMalloc for entry of jobs");
    //     cudaCheck(cudaMemcpy(job_entry, pop->task.jobs, AMOUNT_OF_JOBS,
    //     cudaMemcpyHostToDevice), "cudaMemcpy entry of jobs from host to
    //     device"); address_of_jobs[i] = job_entry;
    // }
    // cudaCheck(cudaMalloc((void**)&jobs, sizeof(job_t
    // *)*AMOUNT_OF_R_CHROMOSOMES), "cudaMalloc for jobs");
    // cudaCheck(cudaMemcpy(jobs, address_of_jobs,
    // sizeof(job_t*)*AMOUNT_OF_R_CHROMOSOMES, cudaMemcpyHostToDevice),
    // "cudaMemcpy jobs from host to device");
    pop->device_objects.jobs = jobs;
    pop->host_objects.address_of_cujobs = address_of_jobs;

    /**************************malloc for process_times**********/
    process_time_t **process_times;
    process_time_t **address_of_process_times;
    cudaCheck(cudaMallocHost((void **) &address_of_process_times,
                             sizeof(process_time_t *) * AMOUNT_OF_JOBS),
              "cudaMallocHost for address of process times");
    for (int i = 0; i < AMOUNT_OF_JOBS; ++i) {
        process_time_t *entry;
        cudaCheck(cudaMalloc((void **) &entry,
                             sizeof(process_time_t) *
                                 pop->task.size_of_process_times[i]),
                  "cudaMalloc for an entry of process times");
        cudaCheck(cudaMemcpy(entry, pop->task.process_times[i],
                             sizeof(process_time_t) *
                                 pop->task.size_of_process_times[i],
                             cudaMemcpyHostToDevice),
                  "cudaMemcpy for an entry of process time");
        address_of_process_times[i] = entry;
    }
    cudaCheck(cudaMalloc((void **) &process_times,
                         sizeof(process_time_t *) * AMOUNT_OF_JOBS),
              "cudaMalloc for process times");
    cudaCheck(cudaMemcpy(process_times, address_of_process_times,
                         sizeof(process_time_t *) * AMOUNT_OF_JOBS,
                         cudaMemcpyHostToDevice),
              "cudaMemcpy for process times");
    pop->device_objects.process_times = process_times;
    pop->host_objects.address_of_process_times_entry = address_of_process_times;

    /***************************malloc for machines***************/
    machine_t **machines;
    machine_t **address_of_machines;
    cudamalloc_objects(machine_t, machines, address_of_machines,
                       pop->task.machines, AMOUNT_OF_MACHINES,
                       AMOUNT_OF_R_CHROMOSOMES, "machines");
    // cudaCheck(cudaMallocHost((void**)&address_of_machines, sizeof(machine_t
    // *)*AMOUNT_OF_R_CHROMOSOMES), "cudaMallocHost for address_of_machines");
    // for(int i = 0; i < AMOUNT_OF_R_CHROMOSOMES; ++i){
    //     cudaCheck(cudaMalloc((void**)&machine_entry,
    //     sizeof(machine_t)*AMOUNT_OF_MACHINES), "cudamalloc for entry of
    //     machines"); cudaCheck(cudaMemcpy(machine_entry, pop->task.machines,
    //     sizeof(machine_t)*AMOUNT_OF_MACHINES, cudaMemcpyHostToDevice),
    //     "cudaMemcpy entry of machines from host to device");
    //     address_of_machines[i] = machine_entry;
    // }
    // cudaCheck(cudaMalloc((void**)&machines,
    // sizeof(machine_t*)*AMOUNT_OF_R_CHROMOSOMES), "cudaMalloc for machines");
    // cudaCheck(cudaMemcpy(machines, address_of_machines,
    // sizeof(machine_t*)*AMOUNT_OF_R_CHROMOSOMES, cudaMemcpyHostToDevice),
    // "cudaMemcpy for machines from hst to dvc");
    pop->device_objects.machines = machines;
    pop->host_objects.address_of_cumachines = address_of_machines;

    /*********************malloc for tools and wires***************/
    tool_t **tools;
    tool_t **address_of_tools;
    cudamalloc_objects(tool_t, tools, address_of_tools, pop->task.tools,
                       AMOUNT_OF_MACHINES, AMOUNT_OF_R_CHROMOSOMES, "tools");
    pop->device_objects.tools = tools;
    pop->host_objects.address_of_tools = address_of_tools;

    wire_t **wires;
    wire_t **address_of_wires;
    cudamalloc_objects(wire_t, wires, address_of_wires, pop->task.wires,
                       AMOUNT_OF_MACHINES, AMOUNT_OF_R_CHROMOSOMES, "wires");
    pop->device_objects.wires = wires;
    pop->host_objects.address_of_wires = address_of_wires;

    // setup chromosomes
    chromosome_base_t *chromosomes;
    double **genes;
    double **address_of_cu_genes;
    double *tmp_genes;
    cudaCheck(cudaMallocHost((void **) &tmp_genes,
                             sizeof(double) * (AMOUNT_OF_JOBS << 1)),
              "cudaMallocHost for genes sample");
    cudaCheck(cudaMalloc((void **) &chromosomes,
                         sizeof(chromosome_base_t) * AMOUNT_OF_R_CHROMOSOMES),
              "cudaMalloc for chromosomes");
    cudaCheck(cudaMallocHost((void **) &address_of_cu_genes,
                             sizeof(double *) * AMOUNT_OF_R_CHROMOSOMES),
              "cudaMallocHost for address of genes");
    for (int i = 0; i < AMOUNT_OF_R_CHROMOSOMES; ++i) {
        double *tmp;
        random(tmp_genes, AMOUNT_OF_JOBS << 1);
        cudaCheck(
            cudaMalloc((void **) &tmp, sizeof(double) * (AMOUNT_OF_JOBS << 1)),
            "cudaMalloc for entry of genes");
        cudaCheck(
            cudaMemcpy(tmp, tmp_genes, sizeof(double) * (AMOUNT_OF_JOBS << 1),
                       cudaMemcpyHostToDevice),
            "cudaMemcpy for genes");
        address_of_cu_genes[i] = tmp;
    }
    cudaCheck(cudaMalloc((void **) &genes,
                         sizeof(double *) * AMOUNT_OF_R_CHROMOSOMES),
              "cudaMalloc for genes");
    cudaCheck(cudaMemcpy(genes, address_of_cu_genes,
                         sizeof(double *) * AMOUNT_OF_R_CHROMOSOMES,
                         cudaMemcpyHostToDevice),
              "cudaMemcpy genes from host to device");

    // setup host chromosomes
    chromosome_base_t *host_chromosomes;
    cudaCheck(cudaMallocHost(
                  (void **) &host_chromosomes,
                  sizeof(chromosome_base_t) * pop->parameters.SWAP_CHROMOSOMES),
              "cudaMallocHost for host chromosomes");


    pop->chromosomes.chromosomes = chromosomes;
    pop->chromosomes.address_of_cugenes = address_of_cu_genes;
    pop->chromosomes.genes = genes;
    pop->chromosomes.host_chromosomes = host_chromosomes;

    //==================prepare evolution factors===========================//
    // device
    device_malloc_evolution_factor(pop->evolution_factors.device.c_selected1,
                                   unsigned int, AMOUNT_OF_R_CHROMOSOMES);
    device_malloc_evolution_factor(pop->evolution_factors.device.c_selected2,
                                   unsigned int, AMOUNT_OF_R_CHROMOSOMES);
    device_malloc_evolution_factor(pop->evolution_factors.device.cut_points,
                                   unsigned int, AMOUNT_OF_R_CHROMOSOMES);
    device_malloc_evolution_factor(pop->evolution_factors.device.range,
                                   unsigned int, AMOUNT_OF_R_CHROMOSOMES);
    device_malloc_evolution_factor(pop->evolution_factors.device.m_selected,
                                   unsigned int, AMOUNT_OF_R_CHROMOSOMES);
    device_malloc_evolution_factor(pop->evolution_factors.device.new_genes,
                                   double, AMOUNT_OF_R_CHROMOSOMES);
    device_malloc_evolution_factor(pop->evolution_factors.device.gene_idx,
                                   unsigned int, AMOUNT_OF_R_CHROMOSOMES);
    // host
    host_malloc_evolution_factor(pop->evolution_factors.host.c_selected1,
                                 unsigned int, AMOUNT_OF_R_CHROMOSOMES);
    host_malloc_evolution_factor(pop->evolution_factors.host.c_selected2,
                                 unsigned int, AMOUNT_OF_R_CHROMOSOMES);
    host_malloc_evolution_factor(pop->evolution_factors.host.cut_points,
                                 unsigned int, AMOUNT_OF_R_CHROMOSOMES);
    host_malloc_evolution_factor(pop->evolution_factors.host.range,
                                 unsigned int, AMOUNT_OF_R_CHROMOSOMES);
    host_malloc_evolution_factor(pop->evolution_factors.host.m_selected,
                                 unsigned int, AMOUNT_OF_R_CHROMOSOMES);
    host_malloc_evolution_factor(pop->evolution_factors.host.new_genes, double,
                                 AMOUNT_OF_R_CHROMOSOMES);
    host_malloc_evolution_factor(pop->evolution_factors.host.gene_idx,
                                 unsigned int, AMOUNT_OF_R_CHROMOSOMES);



    // setup ops
    list_operations_t *list_ops;
    cudaCheck(cudaMalloc((void **) &list_ops, sizeof(list_operations_t)),
              "cudaMalloc for list operations");
    job_base_operations_t *job_ops;
    cudaCheck(cudaMalloc((void **) &job_ops, sizeof(job_base_operations_t)),
              "cudaMalloc for job operations");
    machine_base_operations_t *machine_ops;
    cudaCheck(
        cudaMalloc((void **) &machine_ops, sizeof(machine_base_operations_t) +
                                               sizeof(setup_time_func_t) * 7),
        "cudaMalloc for machine operations");
    pop->device_objects.list_ops = list_ops;
    pop->device_objects.job_ops = job_ops;
    pop->device_objects.machine_ops = machine_ops;

    // initialization
    dim3 machine_chromosome_thread(10, 100);
    dim3 machine_chromosome_block(
        AMOUNT_OF_R_CHROMOSOMES / machine_chromosome_thread.x,
        AMOUNT_OF_MACHINES / machine_chromosome_thread.y + 1);

    dim3 job_chromosome_thread(10, 100);
    dim3 job_chromosome_block(AMOUNT_OF_R_CHROMOSOMES / job_chromosome_thread.x,
                              AMOUNT_OF_JOBS / job_chromosome_thread.y + 1);

    initOps<<<1, 1>>>(list_ops, job_ops, machine_ops);
    initializeChromosomes<<<AMOUNT_OF_R_CHROMOSOMES, 1>>>(
        chromosomes, genes, AMOUNT_OF_JOBS, AMOUNT_OF_R_CHROMOSOMES);
    initializeJobs<<<job_chromosome_block, job_chromosome_thread>>>(
        jobs, process_times, job_ops, chromosomes, AMOUNT_OF_JOBS,
        AMOUNT_OF_R_CHROMOSOMES);
    initializeMachines<<<machine_chromosome_block, machine_chromosome_thread>>>(
        machines, tools, wires, AMOUNT_OF_MACHINES, AMOUNT_OF_R_CHROMOSOMES);

    cudaDeviceSynchronize();
}

void testMachineInformation(population_t *pop, int n)
{
    int AMOUNT_OF_MACHINES = pop->task.AMOUNT_OF_MACHINES;

    unsigned int *d_machine_nums;
    unsigned int *h_machine_nums;
    unsigned int *d_sizeof_jobs;
    unsigned int *h_sizeof_jobs;

    size_t array_size = sizeof(unsigned int) * AMOUNT_OF_MACHINES;

    cudaCheck(cudaMalloc((void **) &d_machine_nums, array_size),
              "cudaMalloc for d_machine_nums");
    cudaCheck(cudaMallocHost((void **) &h_machine_nums, array_size),
              "cudaMalloc for h_machine_nums");
    cudaCheck(cudaMalloc((void **) &d_sizeof_jobs, array_size),
              "cudaMalloc for d_sizeof_jobs");
    cudaCheck(cudaMallocHost((void **) &h_sizeof_jobs, array_size),
              "cudaMalloc for h_sizeof_jobs");

    getMachineInformation<<<AMOUNT_OF_MACHINES, 1>>>(
        pop->host_objects.address_of_cumachines[n], d_machine_nums,
        d_sizeof_jobs, AMOUNT_OF_MACHINES);

    cudaCheck(cudaMemcpy(h_machine_nums, d_machine_nums, array_size,
                         cudaMemcpyDeviceToHost),
              "cudaMemcpy from device to host for machine_num");
    cudaCheck(cudaMemcpy(h_sizeof_jobs, d_sizeof_jobs, array_size,
                         cudaMemcpyDeviceToHost),
              "cudaMemcpy from device to host for sizeof_jobs");

    double **genes;
    double **address_of_genes;
    cudaCheck(cudaMallocHost((void **) &address_of_genes,
                             sizeof(double *) * AMOUNT_OF_MACHINES),
              "cudaMallocHost for entry of genes in testMachineInformation");
    cudaCheck(
        cudaMalloc((void **) &genes, sizeof(double *) * AMOUNT_OF_MACHINES),
        "cudaMallocHost for entries of genes in testMachineInformation");
    for (int i = 0; i < AMOUNT_OF_MACHINES; ++i) {
        double *entry;
        cudaCheck(
            cudaMalloc((void **) &entry, sizeof(double) * h_sizeof_jobs[i]),
            "cudaMalloc for entry of genes");
        address_of_genes[i] = entry;
    }
    cudaCheck(cudaMemcpy(genes, address_of_genes,
                         sizeof(double *) * AMOUNT_OF_MACHINES,
                         cudaMemcpyHostToDevice),
              "cudaMemcpy for genes from host to device");

    // for(int i = 0; i < AMOUNT_OF_MACHINES; ++i){
    //     printf("%s : %u\n",
    //     convertUIntToEntityName(h_machine_nums[i]).c_str(),
    //     h_sizeof_jobs[i]);
    // }
    getMachineJobs<<<AMOUNT_OF_MACHINES, 1>>>(
        pop->host_objects.address_of_cumachines[n], genes, AMOUNT_OF_MACHINES);
    for (int i = 0; i < AMOUNT_OF_MACHINES; ++i) {
        printf("%s[%d] : ", convertUIntToEntityName(h_machine_nums[i]).c_str(),
               h_sizeof_jobs[i]);
        double *entry;
        if (h_sizeof_jobs[i] == 0) {
            printf("\n");
            continue;
        }
        cudaCheck(
            cudaMallocHost((void **) &entry, sizeof(double) * h_sizeof_jobs[i]),
            "cudaMallocHost for entry of genes");
        cudaCheck(cudaMemcpy(entry, address_of_genes[i],
                             sizeof(double) * h_sizeof_jobs[i],
                             cudaMemcpyDeviceToHost),
                  "cudaMemcpy entry from device to host");
        for (int j = 0; j < h_sizeof_jobs[i]; ++j) {
            printf("%.2f ", entry[j]);
        }
        cudaCheck(cudaFreeHost(entry), "cudaFreeHost for entry");
        printf("\n");
    }

    cudaFree(d_machine_nums);
    cudaFree(d_sizeof_jobs);
    cudaFreeHost(h_machine_nums);
    cudaFreeHost(h_sizeof_jobs);
}


double geneticAlgorithm(void *_pop)
{
    population_t *pop = (population_t *) _pop;

    int AMOUNT_OF_R_CHROMOSOMES = pop->parameters.AMOUNT_OF_R_CHROMOSOMES;
    int AMOUNT_OF_JOBS = pop->task.AMOUNT_OF_JOBS;
    int AMOUNT_OF_MACHINES = pop->task.AMOUNT_OF_MACHINES;
    int AMOUNT_OF_CHROMOSOMES = pop->parameters.AMOUNT_OF_CHROMOSOMES;
    int CROSSOVER_AMOUNT;
    int MUTATION_AMOUNT;

    CROSSOVER_AMOUNT = pop->parameters.AMOUNT_OF_CHROMOSOMES *
                           pop->parameters.EVOLUTION_RATE;

    MUTATION_AMOUNT =
            pop->parameters.AMOUNT_OF_CHROMOSOMES - CROSSOVER_AMOUNT;

    int n = max(CROSSOVER_AMOUNT, MUTATION_AMOUNT);


    // setup random states and seed
    curandState_t *states;
    cudaCheck(cudaMalloc((void**)&states, sizeof(curandState_t)*n), "cudaMalloc for curandState");
    unsigned long long *h_seeds, *d_seeds;
    cudaMalloc((void**)&d_seeds, sizeof(unsigned long long)*n);
    cudaMallocHost((void**)&h_seeds, sizeof(unsigned long long)*n);
    for(int i = 0; i < n; ++i){
        h_seeds[i] = time(NULL);
    }
    cudaMemcpy(d_seeds, h_seeds, sizeof(unsigned long long)*n, cudaMemcpyHostToDevice);
    setup_rand_state<<<1, n>>>(states, d_seeds); 


    chromosome_base_t *chrs;
    cudaCheck(cudaMallocHost((void **) &chrs, sizeof(chromosome_base_t) *
                                                  AMOUNT_OF_R_CHROMOSOMES),
              "cudaMalloc for chromosomes for testing");


    dim3 machine_chromosome_thread(10, 100);
    dim3 machine_chromosome_block(
        AMOUNT_OF_R_CHROMOSOMES / machine_chromosome_thread.x,
        AMOUNT_OF_MACHINES / machine_chromosome_thread.y + 1);

    dim3 job_chromosome_thread(10, 100);
    dim3 job_chromosome_block(AMOUNT_OF_R_CHROMOSOMES / job_chromosome_thread.x,
                              AMOUNT_OF_JOBS / job_chromosome_thread.y + 1);

    initializeChromosomes<<<AMOUNT_OF_R_CHROMOSOMES, 1>>>(
        pop->chromosomes.chromosomes, AMOUNT_OF_R_CHROMOSOMES);

    clock_t q1 = clock();
    clock_t q2 = q1 + pop->parameters.GENERATIONS * CLOCKS_PER_SEC;
    for (int i = 0; q1 < q2;  ++i, q1 = clock()) {
        // printf("generation %d\n",i);
        binding<<<job_chromosome_block, job_chromosome_thread>>>(
            pop->device_objects.jobs, pop->chromosomes.chromosomes,
            pop->device_objects.job_ops, AMOUNT_OF_JOBS,
            AMOUNT_OF_R_CHROMOSOMES);

        resetMachines<<<machine_chromosome_block, machine_chromosome_thread>>>(
            pop->device_objects.machines, pop->device_objects.machine_ops,
            AMOUNT_OF_MACHINES, AMOUNT_OF_R_CHROMOSOMES);

        machineSelection<<<job_chromosome_block, job_chromosome_thread>>>(
            pop->device_objects.jobs, pop->device_objects.job_ops,
            AMOUNT_OF_JOBS, AMOUNT_OF_R_CHROMOSOMES);

        machineSelection2<<<machine_chromosome_block,
                            machine_chromosome_thread>>>(
            pop->device_objects.jobs, pop->device_objects.machines,
            pop->device_objects.machine_ops, AMOUNT_OF_JOBS, AMOUNT_OF_MACHINES,
            AMOUNT_OF_R_CHROMOSOMES);

        sortJob<<<machine_chromosome_block, machine_chromosome_thread>>>(
            pop->device_objects.machines, pop->device_objects.list_ops,
            pop->device_objects.machine_ops, AMOUNT_OF_MACHINES,
            AMOUNT_OF_R_CHROMOSOMES);

        scheduling<<<machine_chromosome_block, machine_chromosome_thread>>>(
            pop->device_objects.machines, pop->device_objects.job_ops,
            pop->device_objects.machine_ops, AMOUNT_OF_MACHINES,
            AMOUNT_OF_R_CHROMOSOMES);

        computeFitnessValue<<<AMOUNT_OF_R_CHROMOSOMES, 1>>>(
            pop->device_objects.machines, pop->chromosomes.chromosomes,
            AMOUNT_OF_MACHINES, AMOUNT_OF_R_CHROMOSOMES);

        sortChromosomes<<<1, 1>>>(pop->chromosomes.chromosomes,
                                  AMOUNT_OF_R_CHROMOSOMES);

//        cudaCheck(cudaMemcpy(chrs, pop->chromosomes.chromosomes, sizeof(chromosome_base_t)*AMOUNT_OF_R_CHROMOSOMES, cudaMemcpyDeviceToHost), "cudaMemcpy for chromosomes for testing"); 
//        printf("fitness_value = %f\n", chrs[0].fitnessValue);
//

        generateCrossoverFactors<<<1, CROSSOVER_AMOUNT>>>(
                states, 
                pop->evolution_factors.device.c_selected1,
                pop->evolution_factors.device.c_selected2,
                pop->evolution_factors.device.cut_points,
                pop->evolution_factors.device.range,
                AMOUNT_OF_JOBS << 1,
                CROSSOVER_AMOUNT >> 1,
                AMOUNT_OF_CHROMOSOMES);

        generateMutationFactors<<<1, MUTATION_AMOUNT>>>(states,
                                pop->evolution_factors.device.m_selected,
                                pop->evolution_factors.device.gene_idx,
                                pop->evolution_factors.device.new_genes,
                                AMOUNT_OF_JOBS << 1, // gene size
                                MUTATION_AMOUNT, // factor size
                                AMOUNT_OF_CHROMOSOMES);



        // generateCrossoverFactors(&pop->evolution_factors.host,
        //                          CROSSOVER_AMOUNT >> 1, AMOUNT_OF_JOBS << 1,
        //                          AMOUNT_OF_CHROMOSOMES);
        // generateMutationFactors(&pop->evolution_factors.host, MUTATION_AMOUNT,
        //                         AMOUNT_OF_JOBS << 1, AMOUNT_OF_CHROMOSOMES);
        // cpyEvolutionFactors(&pop->evolution_factors.device,
        //                     &pop->evolution_factors.host,
        //                     AMOUNT_OF_CHROMOSOMES);
        // cudaStreamSynchronize(0);
        // cudaDeviceSynchronize();

        crossover<<<1, CROSSOVER_AMOUNT>>>(
            pop->chromosomes.chromosomes,
            pop->evolution_factors.device.c_selected1,
            pop->evolution_factors.device.c_selected2,
            pop->evolution_factors.device.cut_points,
            pop->evolution_factors.device.range, 0, AMOUNT_OF_JOBS,
            CROSSOVER_AMOUNT >> 1, AMOUNT_OF_CHROMOSOMES);

        mutation<<<1, MUTATION_AMOUNT>>>(
            pop->chromosomes.chromosomes,
            pop->evolution_factors.device.m_selected,
            pop->evolution_factors.device.gene_idx,
            pop->evolution_factors.device.new_genes, CROSSOVER_AMOUNT,
            AMOUNT_OF_JOBS, MUTATION_AMOUNT, AMOUNT_OF_CHROMOSOMES);

        cudaStreamSynchronize(0);
    }

    cudaCheck(cudaMemcpy(chrs, pop->chromosomes.chromosomes, sizeof(chromosome_base_t)*AMOUNT_OF_R_CHROMOSOMES, cudaMemcpyDeviceToHost), "cudaMemcpy for chromosomes for testing"); 
    // printf("fitness_value = %f\n", chrs[0].fitnessValue);
    double fitnessValue = chrs[0].fitnessValue;
    cudaCheck(cudaFreeHost(chrs), "cudaFree chrs");
    // pthread_exit(NULL);
    return fitnessValue;
}


void swapPopulation(population_t pops[], const int AMOUNT_OF_POPULATIONS)
{
    for (int i = 0; i < AMOUNT_OF_POPULATIONS - 1; ++i) {
        cudaCheck(cudaMemcpy(pops[i + 1].chromosomes.host_chromosomes,
                             pops[i].chromosomes.chromosomes,
                             sizeof(chromosome_base_t) *
                                 pops[i].parameters.SWAP_CHROMOSOMES,
                             cudaMemcpyDeviceToHost),
                  "cudaMemcpy d2h for chromosomes");
    }

    cudaCheck(
        cudaMemcpy(
            pops[0].chromosomes.host_chromosomes,
            pops[AMOUNT_OF_POPULATIONS - 1].chromosomes.chromosomes,
            sizeof(chromosome_base_t) * pops[0].parameters.SWAP_CHROMOSOMES,
            cudaMemcpyDeviceToHost),
        "cudaMemcpy d2h for first population");

    // copy to device
    for (int i = 0; i < AMOUNT_OF_POPULATIONS; ++i) {
        cudaCheck(cudaMemcpy(pops[i].chromosomes.chromosomes,
                             pops[i].chromosomes.host_chromosomes,
                             sizeof(chromosome_base_t) *
                                 pops[i].parameters.SWAP_CHROMOSOMES,
                             cudaMemcpyHostToDevice),
                  "cudaMemcpy h2d for chromosomes");
    }
}
