#include <include/chromosome.h>
#include <cstdlib>
#include "include/infra.h"
#include "include/job_base.h"
#include "include/linked_list.h"

using namespace std;

void copyChromosome(chromosome_base_t c1, chromosome_base_t c2)
{
    memcpy(c1.genes, c2.genes, sizeof(double) * c1.gene_size);
}


void crossover(chromosome_base_t p1,
               chromosome_base_t p2,
               chromosome_base_t c1,
               chromosome_base_t c2)
{
    memcpy(c1.genes, p1.genes, sizeof(double) * p1.gene_size);
    memcpy(c2.genes, p2.genes, sizeof(double) * p2.gene_size);
    int cutpoint1 = random_range(0, p1.gene_size, -1);
    int cutpoint2 = random_range(0, p2.gene_size, cutpoint1);

    if (cutpoint1 > cutpoint2)
        std::swap(cutpoint1, cutpoint2);

    int size = cutpoint2 - cutpoint1;

    memcpy(c1.genes + cutpoint1, p2.genes + cutpoint1, sizeof(double) * size);
    memcpy(c2.genes + cutpoint1, p1.genes + cutpoint1, sizeof(double) * size);
}


void mutation(chromosome_base_t p, chromosome_base_t c)
{
    memcpy(c.genes, p.genes, sizeof(double) * p.gene_size);
    int pos = random_range(0, p.gene_size, -1);
    double rnd = (double) rand() / (double) RAND_MAX;
    c.genes[pos] = rnd;
}

double decoding(chromosome_base_t chromosome,
                job_t *jobs,
                std::map<unsigned int, machine_t *> machines,
                machine_base_operations_t *machine_ops,
                list_operations_t *list_ops,
                job_base_operations_t *job_ops,
                int AMOUNT_OF_JOBS)
{
    unsigned int machine_idx, machine_no;
    for (map<unsigned int, machine_t *>::iterator it = machines.begin();
         it != machines.end(); it++) {
        machine_ops->reset(&(it->second->base));
    }
    // machine selection
    for (int j = 0; j < AMOUNT_OF_JOBS; ++j) {
        job_ops->set_ms_gene_addr(&jobs[j].base, chromosome.ms_genes + j);
        job_ops->set_os_gene_addr(&jobs[j].base, chromosome.os_genes + j);
        machine_idx = job_ops->machine_selection(&jobs[j].base);
        machine_no = jobs[j].base.machine_no;
        machine_ops->add_job(&machines[machine_no]->base, &jobs[j].list);
    }

    // sorting;
    for (map<unsigned int, machine_t *>::iterator it = machines.begin();
         it != machines.end(); it++) {
        machine_ops->sort_job(&it->second->base, list_ops);
    }

    // scheduling
    double value = 0;
    for (map<unsigned int, machine_t *>::iterator it = machines.begin();
         it != machines.end(); it++) {
        scheduling(it->second, machine_ops);
        // if(it->second->total_completion_time > value)
        value += it->second->total_completion_time;
    }

    return value;
}

int chromosomeCmpr(const void *_c1, const void *_c2)
{
    chromosome_base_t *c1 = (chromosome_base_t *) _c1;
    chromosome_base_t *c2 = (chromosome_base_t *) _c2;
    if (c1->fitnessValue > c2->fitnessValue)
        return true;
    return false;
}
