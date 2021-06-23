#include <include/job.h>
#include "include/def.h"

__qualifier__ double jobGetValue(void* _self){
    list_ele_t *self = (list_ele_t *) _self;
    job_t *j = (job_t *) self->ptr_derived_object;
    return *(j->base.os_seq_gene);
}

__qualifier__ void job_initialize(job_t *self){
    _list_init(&self->list);
    self->list.ptr_derived_object = self;
    self->list.get_value = jobGetValue;
    
    job_base_init(&self->base);
    self->base.ptr_derived_object = self;
}
