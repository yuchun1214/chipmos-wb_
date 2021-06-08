#ifndef __JOB_H__
#define __JOB_H__

#include <include/job_base.h>
#include <include/linked_list.h>
#include <include/infra.h>

typedef struct __info_t job_info_t;

typedef struct job_t{
    job_info_t part_no;
    job_info_t pin_package;
    char urgent_code;
    job_base_t base;
    list_ele_t list;
}job_t;


#endif
