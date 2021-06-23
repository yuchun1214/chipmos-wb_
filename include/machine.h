#ifndef __MACHINE_H__
#define __MACHINE_H__

#include <include/job.h>
#include <include/machine_base.h>
#include <include/infra.h>
#include "include/def.h"
#include "include/job_base.h"

typedef struct __info_t machine_info_t;
typedef struct __info_t tool_info_t;
typedef struct __info_t wire_info_t;

typedef struct ancillary_resource_t{
    unsigned int no;
    struct __info_t name;
    unsigned int machine_no;
    double time;
}ares_t;

typedef ares_t tool_t;
typedef ares_t wire_t;


typedef struct __machine_t{
    machine_base_t base;
    ares_t * tool;
    ares_t * wire;
    double total_completion_time;
    double makespan;
}machine_t;

bool ares_ptr_comp(ares_t *, ares_t *);
bool ares_comp(ares_t , ares_t);

__qualifier__ void machine_reset(machine_base_t * base);

__qualifier__ double setup_time_CWN(job_base_t *, job_base_t*);
__qualifier__ double setup_time_CK(job_base_t*, job_base_t *);
__qualifier__ double setup_time_EU(job_base_t*, job_base_t *);
__qualifier__ double setup_time_MC_SC(job_base_t*, job_base_t *);
__qualifier__ double setup_time_CSC(job_base_t*, job_base_t *);
__qualifier__ double setup_time_USC(job_base_t*, job_base_t *);

__qualifier__ void initMachine(machine_t *machine);

__qualifier__ void scheduling(machine_t *machine, job_base_operations_t * jbops, machine_base_operations_t *ops);

#endif
