#include <include/machine.h>
#include <string.h>

bool ares_ptr_comp(ares_t *a1, ares_t *a2)
{
    return a1->time < a2->time;
}

bool ares_comp(ares_t a1, ares_t a2)
{
    return a1.time < a2.time;
}

void machine_reset(machine_base_t *base)
{
    machine_t *m = (machine_t *) base->ptr_derived_object;
    machine_base_reset(base);
    m->makespan = 0;
    m->total_completion_time = 0;
    m->tool->time = 0;
    m->wire->time = 0;
}

double setup_time_CWN(job_base_t *_prev, job_base_t *_next)
{
    if (_prev) {
        job_t *prev = (job_t *) _prev->ptr_derived_object;
        job_t *next = (job_t *) _next->ptr_derived_object;
        if (isSameInfo(prev->part_no, next->part_no))
            return 0.0;
        else
            return 30.0;
    }
    return 0;
}

double setup_time_CK(job_base_t *_prev, job_base_t *_next)
{
    if (_prev) {
        job_t *prev = (job_t *) _prev->ptr_derived_object;
        job_t *next = (job_t *) _next->ptr_derived_object;
        if (prev->part_no.data.text[5] == next->part_no.data.text[5]) {
            return 0.0;
        } else {
            return 0;
        }
    }
    return 0;
}

double setup_time_EU(job_base_t *_prev, job_base_t *_next)
{
    if (_next) {
        job_t *next = (job_t *) _next->ptr_derived_object;
        if (next->urgent_code)
            return 0;  // FIXME
        else
            return 0;
    }
    return -1;
}

double setup_time_MC_SC(job_base_t *_prev, job_base_t *_next)
{
    if (_prev && _next) {
        job_t *prev = (job_t *) _prev->ptr_derived_object;
        job_t *next = (job_t *) _next->ptr_derived_object;
        if (prev->pin_package.data.number[0] ==
            next->pin_package.data.number[0]) {
            return 84;
        } else {
            return 90;
        }
    }
    return 0;
}

double setup_time_CSC(job_base_t *_prev, job_base_t *_next)
{
    if (_next) {
        job_t *next = (job_t *) _next->ptr_derived_object;
        if (next->part_no.data.text[5] != 'A')
            return 0;
    }
    return 0;
}

double setup_time_USC(job_base_t *_prev, job_base_t *_next)
{
    if (_next) {
        job_t *next = (job_t *) _next->ptr_derived_object;
        if (next->part_no.data.text[5] != 'A' && next->urgent_code == 'Y')
            return 72;
        else
            return 0;
    }
    return -1;
}

double calculateSetupTime(job_t *prev,
                          job_t *next,
                          machine_base_operations_t *ops)
{
    // for all setup time function
    double time = 0;
    for (unsigned int i = 0; i < ops->sizeof_setup_time_function_array; ++i) {
        if (prev)
            time += ops->setup_times[i](&prev->base, &next->base);
        else
            time += ops->setup_times[i](NULL, &next->base);
    }
    return time;
}

void scheduling(machine_t *machine, machine_base_operations_t *ops)
{
    list_ele_t *it;
    job_t *job;
    job_t *prev_job = NULL;
    it = machine->base.root;
    double arrival_time;
    double setup_time;
    bool hasICSI = false;
    double start_time = machine->base.avaliable_time;
    double total_completion_time = 0;
    while (it) {
        job = (job_t *) it->ptr_derived_object;
        arrival_time = job->base.arriv_t;
        setup_time = calculateSetupTime(prev_job, job, ops);
        if (!hasICSI) {
            if (strncmp(job->customer.data.text, "ICSI", 4) == 0) {
                setup_time += 54;
                hasICSI = true;
            }
        }

        start_time = (start_time + setup_time) > arrival_time
                         ? start_time + setup_time
                         : arrival_time;
        set_start_time(&job->base, start_time);
        start_time = get_end_time(&job->base);
        total_completion_time += start_time;
        prev_job = job;
        it = it->next;
    }
    machine->makespan = start_time;
    machine->tool->time = start_time;
    machine->wire->time = start_time;
    machine->total_completion_time = total_completion_time;
}
