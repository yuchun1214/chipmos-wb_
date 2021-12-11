#include <include/job.h>

double jobGetValue(void *_self)
{
    list_ele_t *self = (list_ele_t *) _self;
    job_t *j = (job_t *) self->ptr_derived_object;
    return *(j->base.os_seq_gene);
}
