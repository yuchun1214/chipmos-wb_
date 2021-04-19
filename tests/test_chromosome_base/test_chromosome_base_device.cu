#include <bits/floatn.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <driver_types.h>
#include <include/linked_list.h>
#include <include/machine_base.h>
#include <include/chromosome_base.h>
#include <tests/include/test_chromosome_base.h>
#include <tests/include/test_machine_base.h>
#include <tests/include/def.h>

size_t memory_usage = 0;

cudaError_t test_cudaMalloc(void ** ptr, size_t size){
	memory_usage += size;
	return cudaMalloc(ptr, size);
}

#define cudaMalloc test_cudaMalloc 
		

extern int JOB_AMOUNT;
extern int MACHINE_AMOUNT;
extern int CHROMOSOME_AMOUNT;
extern int GENERATIONS;

class TestChromosomeBaseDevice : public testing::Test{
public:
	Machine ** machines;
	Machine ** address_machines_arr;
	job_t ** jobs;
	job_t ** address_jobs_arr;
	Chromosome * chromosomes;

	process_time_t **processTimes;
	process_time_t **address_process_time_arr;

	double *genes;
	double *host_genes;
	unsigned int *device_can_run_machine_size;
	unsigned int *host_can_run_machine_size;

	size_t gene_size;

	list_operations_t  *ops;
	
	int R_JOB_AMOUNT;
	int R_MACHINE_AMOUNT;
	int R_CHROMOSOME_AMOUNT;

	void random_shuffle(double *genes, size_t size);

	void SetUp() override;
	void TearDown() override;
};

void TestChromosomeBaseDevice::random_shuffle(double *genes, size_t size){
	for(unsigned int i = 0; i < size; ++i){
		genes[i] = (double)rand() / (double)RAND_MAX;	
	}
}

void TestChromosomeBaseDevice::SetUp(){
	R_JOB_AMOUNT = JOB_AMOUNT * (CHROMOSOME_AMOUNT<<1);
	R_MACHINE_AMOUNT = MACHINE_AMOUNT * (CHROMOSOME_AMOUNT<<1);
	R_CHROMOSOME_AMOUNT = CHROMOSOME_AMOUNT << 1;
	
	// allocating jobs
	cudaCheck(cudaMalloc((void**)&jobs, sizeof(job_t*) * R_CHROMOSOME_AMOUNT), "allocating jobs...");
	cudaCheck(cudaMallocHost((void**)&address_jobs_arr, sizeof(job_t*) * R_CHROMOSOME_AMOUNT), "allocating host address_job_arr");
	


	job_t * tmp;
	for(int i = 0; i < R_CHROMOSOME_AMOUNT; ++i){
		cudaCheck(cudaMalloc((void**)&tmp, sizeof(job_t) * JOB_AMOUNT), "allocating jobs for a chromosome");
		address_jobs_arr[i] = tmp;
	}
	cudaCheck(cudaMemcpy(jobs, address_jobs_arr, sizeof(job_t*) * R_CHROMOSOME_AMOUNT, cudaMemcpyHostToDevice), "copy jobs from host to device");
	
	
	// allocating machines
	cudaCheck( cudaMalloc((void**)&machines, sizeof(Machine*)*R_CHROMOSOME_AMOUNT), "alloating machines...");
	cudaCheck( cudaMallocHost((void**)&address_machines_arr, sizeof(Machine*)*R_CHROMOSOME_AMOUNT), "allocating host address_machines_arr");
	Machine *machines_tmp;
	for(int i = 0; i < R_CHROMOSOME_AMOUNT; ++i){
		cudaCheck( cudaMalloc((void**)&machines_tmp, sizeof(Machine)*MACHINE_AMOUNT), "allocating machines for a chromosome");
		address_machines_arr[i] = machines_tmp;
	}
	cudaCheck( cudaMemcpy(machines, address_machines_arr, sizeof(Machine*)*R_CHROMOSOME_AMOUNT, cudaMemcpyHostToDevice), "copy machines from host to device");

	// allocating chromosomes
	cudaCheck( cudaMalloc((void**)&chromosomes, sizeof(Chromosome)*R_CHROMOSOME_AMOUNT), "allocating chromosomes");

	
	// prepare host_can_run_machine_size
	cudaCheck( cudaMallocHost((void**)&host_can_run_machine_size, sizeof(unsigned int)*JOB_AMOUNT), "allocating host_can_run_machine_size on host");
	cudaCheck( cudaMalloc((void**)&device_can_run_machine_size, sizeof(unsigned int)*JOB_AMOUNT), "allocating device_can_run_machines_size on device");
	for(int i = 0; i < JOB_AMOUNT; ++i){
		host_can_run_machine_size[i] = rand() % 200 + 400;
	}
	cudaCheck(cudaMemcpy(device_can_run_machine_size, host_can_run_machine_size, sizeof(unsigned int)*JOB_AMOUNT, cudaMemcpyHostToDevice), "copy can run tool");

	// prepare process_time
	cudaCheck( cudaMallocHost((void**)&address_process_time_arr, sizeof(process_time_t *)*JOB_AMOUNT), "allocating process time on host");
	cudaCheck( cudaMalloc((void**)&processTimes, sizeof(process_time_t *)*JOB_AMOUNT), "allocating process time on device");
	process_time_t *process_time_tmp_host;
	process_time_t *process_time_tmp;
	for(int i = 0; i < JOB_AMOUNT; ++i){
		cudaCheck(cudaMalloc((void**)&process_time_tmp, sizeof(process_time_t) * host_can_run_machine_size[i]), "allocating process time on device");
		cudaCheck(cudaMallocHost((void**)&process_time_tmp_host, sizeof(process_time_t) * host_can_run_machine_size[i]), "allocating process time on host");

		for(unsigned int j = 0; j < host_can_run_machine_size[i]; ++j){
			process_time_tmp_host[j].machine_no = rand() % MACHINE_AMOUNT;
			process_time_tmp_host[j].process_time = rand() % 1000;
		}

		cudaCheck(cudaMemcpy(process_time_tmp, process_time_tmp_host, sizeof(process_time_t) * host_can_run_machine_size[i], cudaMemcpyHostToDevice), "copy process time from host to deivce");
		cudaCheck(cudaFreeHost(process_time_tmp_host), "cuda free process_time_tmp_host");
		address_process_time_arr[i] = process_time_tmp;	
	}
	cudaCheck( cudaMemcpy(processTimes, address_process_time_arr, sizeof(process_time_t *)*JOB_AMOUNT, cudaMemcpyHostToDevice), "copy can run tool from host to device");

	// alloc genes
	cudaCheck(cudaMalloc((void**)&genes, sizeof(double)*(JOB_AMOUNT<<1)*(CHROMOSOME_AMOUNT<<1)),"cuda alloc genes");
	cudaCheck(cudaMallocHost((void**)&host_genes, sizeof(double)*(JOB_AMOUNT<<1)*(CHROMOSOME_AMOUNT<<1)),"cuda malloc hostld");

	// alloc ops
	cudaCheck(cudaMalloc((void**)&ops, sizeof(list_operations_t)), "alloc ops");
}

void TestChromosomeBaseDevice::TearDown(){
	// free ops
	cudaCheck(cudaFree(ops), "Free ops...");

	// free jobs
	cudaCheck(cudaFree(jobs), "Free jobs");
	for(int i = 0; i < R_CHROMOSOME_AMOUNT; ++i){
		cudaCheck(cudaFree(address_jobs_arr[i]), "Free an array of jobs");
	}
	cudaCheck( cudaFreeHost(address_jobs_arr), "Free address_job_arr");

	// free machines
	cudaCheck(cudaFree(machines), "Free machines");
	for(int i = 0; i < R_CHROMOSOME_AMOUNT; ++i){
		cudaCheck(cudaFree(address_machines_arr[i]), "Free an array of machines");
	}
	cudaCheck(cudaFreeHost(address_machines_arr), "Free addres_machines_arr");

	// free can_run_machine
	cudaCheck(cudaFree(device_can_run_machine_size), "Free device_can_run_machine_size");
	cudaCheck(cudaFreeHost(host_can_run_machine_size), "Free host_can_run_mahcine_size");

	// free chromosomes
	cudaCheck(cudaFree(chromosomes), "Free chromosomes");

	// free process time
	cudaCheck(cudaFree(processTimes), "Free processTimes");
	for(int i = 0; i < JOB_AMOUNT; ++i){
		cudaCheck(cudaFree(address_process_time_arr[i]), "Free an array of process_time");
	}
	cudaCheck(cudaFreeHost(address_process_time_arr), "Free address_process_time_arr");
}

__global__ void machineSetup(Machine **machines, int MACHINE_AMOUNT, int CHROMOSOME_AMOUNT){
	int x = threadIdx.x + blockDim.x * blockIdx.x;
	int y = threadIdx.y + blockDim.y * blockIdx.y; 
	if(x < CHROMOSOME_AMOUNT && y < MACHINE_AMOUNT){
		machines[x][y].base.init = initMachineBase;
		initMachine(&machines[x][y]);
	}
}

__global__ void chromosomeSetup(Chromosome *chromosomes, double * genes, int JOB_AMOUNT, int CHROMOSOME_AMOUNT){
	int idx = threadIdx.x + blockDim.x * blockIdx.x;
	if(idx < CHROMOSOME_AMOUNT){
		chromosomes[idx].val = idx;
		chromosomes[idx].base.gene_size = JOB_AMOUNT<<1;
		chromosomes[idx].base.chromosome_no = idx;
		initChromosomeBase(&chromosomes[idx].base, genes + idx*(JOB_AMOUNT<<1));
	}
}

__global__ void jobSetup(job_t ** jobs, unsigned int *can_run_tool_size, process_time_t ** process_times, int JOB_AMOUNT, int CHROMOSOME_AMOUNT){
	int x = threadIdx.x + blockDim.x * blockIdx.x;
	int y = threadIdx.y + blockDim.y * blockIdx.y;
	if(x < CHROMOSOME_AMOUNT && y < JOB_AMOUNT){
		initJob(&jobs[x][y]);
		jobs[x][y].base.job_no = y;
		jobs[x][y].base.setProcessTime(&jobs[x][y].base, process_times[y], can_run_tool_size[y]);
	}
}

__global__ void jobBindGenes(job_t **jobs, Chromosome * chromosomes, int JOB_AMOUNT, int R_CHROMOSOME_AMOUNT){
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y; 
	if(x < R_CHROMOSOME_AMOUNT && y < JOB_AMOUNT){
		jobs[x][y].base.setMsGenePointer(&jobs[x][y].base, chromosomes[x].base.ms_genes + y);
		jobs[x][y].base.setOsSeqGenePointer(&jobs[x][y].base, chromosomes[x].base.os_genes + y);
	}
}

__global__ void machineSelection(job_t **jobs, int JOB_AMOUNT, int R_CHROMOSOME_AMOUNT){
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int machine_idx;
	if(x < R_CHROMOSOME_AMOUNT && y < JOB_AMOUNT){
		machine_idx = jobs[x][y].base.machineSelection(&jobs[x][y].base);
		jobs[x][y].base.machine_no = jobs[x][y].base.process_time[machine_idx].machine_no;
	}
}

__global__ void machineSelection2(job_t **jobs, Machine **machines, int JOB_AMOUNT, int MACHINE_AMOUNT, int R_CHROMOSOME_AMOUNT){
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	if( x < R_CHROMOSOME_AMOUNT && y < MACHINE_AMOUNT){
		for(int i = 0; i < JOB_AMOUNT; ++i){
			if(jobs[x][i].base.machine_no == y){
				machines[x][y].base.addJob(&machines[x][y].base, &jobs[x][i]);
			}
		}
	}
}

__global__ void sortJob(Machine **machines, list_operations_t *ops, int MACHINE_AMOUNT, int R_CHROMOSOME_AMOUNT){
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	if( x < R_CHROMOSOME_AMOUNT && y < MACHINE_AMOUNT){
		machines[x][y].base.sortJob(&machines[x][y].base, ops);	
	}
}

// __global__ void machineSelectJob(Machine *machines, job_t *jobs,

__global__ void operationSetup(list_operations_t *ops){
	ops->init = initList;
	ops->setNext = __listEleSetNext;
	ops->setPrev = __listEleSetPrev;
}

TEST_F(TestChromosomeBaseDevice, test_chromosome_base_device){
	// setup grid dimension
	dim3 machine_chromosome_thread(32, 32);
	dim3 machine_chromosome_block(R_CHROMOSOME_AMOUNT >> 5, MACHINE_AMOUNT >> 5);
	
	dim3 job_chromosome_thread(32, 32);
	dim3 job_chromosome_block(R_CHROMOSOME_AMOUNT >> 5, JOB_AMOUNT >> 5); // (R_CHROMOSOME_AMOUNT / 32, )

	// setup kernel
	jobSetup<<<job_chromosome_block, job_chromosome_thread>>>(jobs, device_can_run_machine_size, processTimes , JOB_AMOUNT, R_CHROMOSOME_AMOUNT);
	machineSetup<<<machine_chromosome_block, machine_chromosome_thread>>>(machines, MACHINE_AMOUNT, R_CHROMOSOME_AMOUNT);
	chromosomeSetup<<<100, 100>>>(chromosomes, genes, JOB_AMOUNT, R_CHROMOSOME_AMOUNT);

	jobBindGenes<<<job_chromosome_block, job_chromosome_thread>>>(jobs, chromosomes, JOB_AMOUNT, R_CHROMOSOME_AMOUNT);
	operationSetup<<<1, 1>>>(ops);
	cudaDeviceSynchronize();

	PRINTF("Device Memory Usage = %lu\n", memory_usage);
	cudaEvent_t startEvent, stopEvent;
	cudaCheck(cudaEventCreate(&startEvent), "create start event");
	cudaCheck(cudaEventCreate(&stopEvent), "create stop event");
	// start computing...
	// machine selection
	// PRINTF("Start Computing...\n");
	cudaCheck(cudaEventRecord(startEvent, 0), "cuda event record start event");
//	cudaProfilerStart();
	machineSelection<<<job_chromosome_block, job_chromosome_thread>>>(jobs, JOB_AMOUNT, R_CHROMOSOME_AMOUNT);
	// cudaDeviceSynchronize();
	// PRINTF("Finish machine selection part 1\n");
	// PRINTF("Start machine selection part2\n");
	machineSelection2<<<machine_chromosome_block, machine_chromosome_thread>>>(jobs, machines, JOB_AMOUNT, MACHINE_AMOUNT, R_CHROMOSOME_AMOUNT);
	// cudaDeviceSynchronize();
	// PRINTF("Finish machine selection part2\n");
	sortJob<<<machine_chromosome_block, machine_chromosome_thread>>>(machines, ops, MACHINE_AMOUNT, R_CHROMOSOME_AMOUNT);
//	cudaProfilerStop();
	// cudaDeviceSynchronize();
	cudaCheck(cudaEventRecord(stopEvent, 0), "cuda event record stop event");
	cudaCheck(cudaEventSynchronize(stopEvent), "cuda event sync stop event");

	float ms;
	cudaCheck(cudaEventElapsedTime(&ms, startEvent, stopEvent), "get elapsed time");

	PRINTF("Elapsed Time : %.3f\n", ms / 1000.0);

	// PRINTF("Finish sorting jobs\n");
}
