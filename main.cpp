#include <ctime>
#include <exception>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include "include/chromosome_base.h"
#include "include/job_base.h"
#include "include/linked_list.h"
#include "include/machine_base.h"

#include <include/arrival.h>
#include <include/chromosome.h>
#include <include/condition_card.h>
#include <include/csv.h>
#include <include/da.h>
#include <include/entity.h>
#include <include/infra.h>
#include <include/lot.h>
#include <include/population.h>
#include <include/route.h>


using namespace std;

round_t createARound(vector<lot_group_t> group,
                     machines_t &machines,
                     ancillary_resources_t &tools,
                     ancillary_resources_t &wires);


void initializeARound(round_t *r);

void initializePopulation(population_t *pop,
                          machines_t &machines,
                          ancillary_resources_t &tools,
                          ancillary_resources_t &wires,
                          int round);

void initializeOperations(population_t *pop);


double geneticAlgorithm(population_t *pop);

void freeJobs(round_t *round);
void freeResources(round_t *round);
void freeChromosomes(population_t *pop);

int main(int argc, const char *argv[])
{
    csv_t arguments(argv[1], "r", true, true);
    map<string, string> elements = arguments.getElements(0);

    csv_t lot_csv(elements["nlots"], "r", true, true);
    lot_csv.setHeaders(
        map<string, string>({{"route", "route"},
                             {"lot_number", "lot_number"},
                             {"pin_package", "pin_package"},
                             {"bd_id", "recipe"},
                             {"prod_id", "prod_id"},
                             {"part_id", "part_id"},
                             {"part_no", "part_no"},
                             {"urgent_code", "urgent_code"},
                             {"qty", "qty"},
                             {"dest_oper", "dest_oper"},
                             {"oper", "dest_oper"},
                             {"hold", "hold"},
                             {"mvin", "mvin"},
                             {"queue_time", "queue_time"},
                             {"fcst_time", "fcst_time"},
                             {"amount_of_tools", "amount_of_tools"},
                             {"amount_of_wires", "amount_of_wires"},
                             {"CAN_RUN_MODELS", "CAN_RUN_MODELS"},
                             {"PROCESS_TIME", "PROCESS_TIME"},
                             {"uphs", "uphs"},
                             {"customer", "customer"}}));
    lot_csv.trim(" ");
    vector<lot_t> allots;
    for (int i = 0, nrows = lot_csv.nrows(); i < nrows; ++i) {
        allots.push_back(lot_t(lot_csv.getElements(i)));
    }

    lots_t lots;
    lots.addLots(allots);


    ancillary_resources_t tools(lots.amountOfTools());
    ancillary_resources_t wires(lots.amountOfWires());

    csv_t machine_csv(elements["nmachines"], "r", true, true);
    machine_csv.trim(" ");
    machine_csv.setHeaders(map<string, string>({{"entity", "ENTITY"},
                                                {"model", "MODEL"},
                                                {"recover_time", "OUTPLAN"}}));

    csv_t location_csv(elements["locations"], "r", true, true);
    location_csv.trim(" ");
    location_csv.setHeaders(
        map<string, string>({{"entity", "Entity"}, {"location", "Location"}}));


    char *text = strdup(elements["std_time"].c_str());
    entities_t entities(text);
    entities.addMachines(machine_csv, location_csv);
    machines_t machines;
    machines.addMachines(entities.getAllEntity());

    // vector<entity_t *> test_ents = entities.getAllEntity();
    // iter(test_ents, i){
    //     printf("%s : %f\n", test_ents[i]->entity_name.c_str(),
    //     test_ents[i]->recover_time);
    // }

    vector<vector<lot_group_t> > round_groups = lots.rounds(entities);

    double tme = 60;

    population_t pop =
        population_t{.parameters = {.AMOUNT_OF_CHROMOSOMES = 100,
                                    .AMOUNT_OF_R_CHROMOSOMES = 200,
                                    .EVOLUTION_RATE = 0.8,
                                    .SELECTION_RATE = 0.2,
                                    .GENERATIONS = (tme / round_groups.size())},
                     .groups = round_groups,
                     .current_round_no = 0};

    // printf("GENERATIONS = %d\n", pop.parameters.GENERATIONS);
    double total_fitness_val = 0;
    srand(time(NULL));
    initializeOperations(&pop);
    iter(pop.groups, i)
    {
        initializePopulation(&pop, machines, tools, wires, i);
        // printf("amount of jobs = %d\n", pop.round.AMOUNT_OF_JOBS);
        total_fitness_val += geneticAlgorithm(&pop);
        freeJobs(&pop.round);
        freeResources(&pop.round);
        freeChromosomes(&pop);
    }
    printf("%f\n", total_fitness_val);
    return 0;
}

void initializeOperations(population_t *pop)
{
    machine_base_operations_t *machine_ops;
    machine_ops = (machine_base_operations_t *) malloc(
        sizeof(machine_base_operations_t) + sizeof(setup_time_t) * 7);
    machine_ops->add_job = _machine_base_add_job;
    machine_ops->sort_job = _machine_base_sort_job;
    machine_ops->setup_times[0] = setup_time_CWN;
    machine_ops->setup_times[1] = setup_time_CK;
    machine_ops->setup_times[2] = setup_time_EU;
    machine_ops->setup_times[3] = setup_time_MC_SC;
    machine_ops->setup_times[4] = setup_time_CSC;
    machine_ops->setup_times[5] = setup_time_USC;
    machine_ops->sizeof_setup_time_function_array = 6;
    machine_ops->sort_job = _machine_base_sort_job;
    machine_ops->add_job = _machine_base_add_job;
    machine_ops->reset = machine_reset;

    list_operations_t *list_ops;
    list_ops = (list_operations_t *) malloc(sizeof(list_operations_t));
    *list_ops = LINKED_LIST_OPS;

    job_base_operations_t *job_ops;
    job_ops = (job_base_operations_t *) malloc(sizeof(job_base_operations_t));
    *job_ops = JOB_BASE_OPS;

    pop->operations.machine_ops = machine_ops;
    pop->operations.list_ops = list_ops;
    pop->operations.job_ops = job_ops;
}

void initializePopulation(population_t *pop,
                          machines_t &machines,
                          ancillary_resources_t &tools,
                          ancillary_resources_t &wires,
                          int round)
{
    pop->round = createARound(pop->groups[round], machines, tools, wires);
    pop->chromosomes = (chromosomes_t *) malloc(
        sizeof(chromosomes_t) * pop->parameters.AMOUNT_OF_CHROMOSOMES);

    for (int i = 0; i < pop->parameters.AMOUNT_OF_CHROMOSOMES; ++i) {
        pop->chromosomes[i].velocity = 0.0;
        pop->chromosomes[i].current_chromosome.chromosome_no = i;
        pop->chromosomes[i].history_best_chromosome.chromosome_no = i;
        pop->chromosomes[i].history_best_chromosome.fitnessValue = 0.0;

        pop->chromosomes[i].current_chromosome.gene_size =
            pop->chromosomes[i].history_best_chromosome.gene_size =
                pop->round.AMOUNT_OF_JOBS << 1;

        pop->chromosomes[i].current_chromosome.genes = (double *) malloc(
            sizeof(double) * (pop->round.AMOUNT_OF_JOBS << 1));

        pop->chromosomes[i].history_best_chromosome.genes = (double *) malloc(
            sizeof(double) * (pop->round.AMOUNT_OF_JOBS << 1));

        pop->chromosomes[i].current_chromosome.ms_genes =
            pop->chromosomes[i].current_chromosome.genes;

        pop->chromosomes[i].current_chromosome.os_genes =
            pop->chromosomes[i].current_chromosome.genes +
            pop->round.AMOUNT_OF_JOBS;

        random(pop->chromosomes[i].current_chromosome.genes,
               pop->chromosomes[i].current_chromosome.gene_size);
    }

    pop->global_best_chromosome.genes =
        (double *) malloc(sizeof(double) * (pop->round.AMOUNT_OF_JOBS << 1));
    pop->global_best_chromosome.fitnessValue = 0.0;
}


round_t prepareResources(vector<lot_group_t> group,
                         machines_t &machines,
                         ancillary_resources_t &tools,
                         ancillary_resources_t &wires)
{
    int AMOUNT_OF_MACHINES = 0;
    vector<tool_t *> alltools;
    vector<wire_t *> allwires;
    map<string, machine_t *> allmachines = machines.getMachines();
    vector<machine_t *> selected_machines;
    map<unsigned int, machine_t *> machines_map;


    vector<tool_t *> ts;
    vector<wire_t *> ws;
    iter(group, i)
    {
        ts = tools.aRound(group[i].tool_name, group[i].machine_amount);
        ws = wires.aRound(group[i].wire_name, group[i].machine_amount);
        iter(group[i].entities, j)
        {
            machine_t *m = allmachines[group[i].entities[j]->entity_name];
            m->base.ptr_derived_object = m;
            m->tool = ts[j];
            m->wire = ws[j];
            strcpy(ts[j]->name.data.text, group[i].tool_name.c_str());
            strcpy(ws[j]->name.data.text, group[i].wire_name.c_str());
            ts[j]->machine_no =
                convertEntityNameToUInt(group[i].entities[j]->entity_name);
            ws[j]->machine_no =
                convertEntityNameToUInt(group[i].entities[j]->entity_name);

            // set the recover time max(machine, tool, wire);
            double mx = max(m->base.avaliable_time, m->tool->time);
            mx = max(mx, m->wire->time);
            // TODO: calculate setup time
            m->base.avaliable_time = mx;

            machines_map[m->base.machine_no] = m;
            alltools.push_back(ts[j]);
            allwires.push_back(ws[j]);
            AMOUNT_OF_MACHINES += 1;
        }
    }

    // prepare tool
    tool_t **tools_ar;
    wire_t **wires_ar;
    tools_ar = (tool_t **) malloc(sizeof(tool_t *) * AMOUNT_OF_MACHINES);
    wires_ar = (wire_t **) malloc(sizeof(wire_t *) * AMOUNT_OF_MACHINES);

    iter(alltools, i)
    {
        tools_ar[i] = alltools[i];
        wires_ar[i] = allwires[i];
    }


    return round_t{.AMOUNT_OF_MACHINES = AMOUNT_OF_MACHINES,
                   .machines = machines_map,
                   .tools = tools_ar,
                   .wires = wires_ar};
}

round_t prepareJobs(vector<lot_group_t> group)
{
    int AMOUNT_OF_JOBS = 0;
    vector<lot_t *> lots;
    iter(group, i) { lots += group[i].lots; }

    AMOUNT_OF_JOBS = lots.size();
    job_t *jobs = (job_t *) malloc(sizeof(job_t) * AMOUNT_OF_JOBS);
    process_time_t **pts =
        (process_time_t **) malloc(sizeof(process_time_t *) * AMOUNT_OF_JOBS);
    int *size_of_pt = (int *) malloc(sizeof(int) * AMOUNT_OF_JOBS);

    // prepare jobs
    // set process time
    iter(lots, i)
    {
        jobs[i] = lots[i]->job();
        vector<string> can_run_ents = lots[i]->getCanRunEntities();
        map<string, double> ent_process_time = lots[i]->getEntitiyProcessTime();
        pts[i] = (process_time_t *) malloc(sizeof(process_time_t) *
                                           can_run_ents.size());
        size_of_pt[i] = can_run_ents.size();
        iter(can_run_ents, j)
        {
            pts[i][j].machine_no = convertEntityNameToUInt(can_run_ents[j]);
            pts[i][j].process_time = ent_process_time[can_run_ents[j]];
        }
        set_process_time(&jobs[i].base, pts[i], size_of_pt[i]);
        jobs[i].base.ptr_derived_object = &jobs[i];
        jobs[i].list.ptr_derived_object = &jobs[i];
    }

    return round_t{.AMOUNT_OF_JOBS = AMOUNT_OF_JOBS,
                   .jobs = jobs,
                   .process_times = pts,
                   .size_of_process_times = size_of_pt};
}

void freeJobs(round_t *round)
{
    free(round->jobs);
    for (int i = 0; i < round->AMOUNT_OF_JOBS; ++i) {
        free(round->process_times[i]);
    }
    free(round->process_times);
    free(round->size_of_process_times);

    round->jobs = NULL;
    round->process_times = NULL;
    round->size_of_process_times = NULL;
}

void freeResources(round_t *round)
{
    free(round->tools);
    free(round->wires);
    round->tools = NULL;
    round->wires = NULL;
}

void freeChromosomes(population_t *pop)
{
    for (int i = 0; i < pop->parameters.AMOUNT_OF_CHROMOSOMES; ++i) {
        free(pop->chromosomes[i].current_chromosome.genes);
        free(pop->chromosomes[i].history_best_chromosome.genes);
    }
    free(pop->chromosomes);
    pop->chromosomes = NULL;
}

round_t createARound(vector<lot_group_t> group,
                     machines_t &machines,
                     ancillary_resources_t &tools,
                     ancillary_resources_t &wires)
{
    round_t round_res = prepareResources(group, machines, tools, wires);
    round_t round_jobs = prepareJobs(group);
    round_res.AMOUNT_OF_JOBS = round_jobs.AMOUNT_OF_JOBS;
    round_res.jobs = round_jobs.jobs;
    round_res.process_times = round_jobs.process_times;
    round_res.size_of_process_times = round_jobs.size_of_process_times;
    return round_res;
}

void compareAndAssignChromosomes(chromosome_base_t *c1, chromosome_base_t *c2)
{
    if (c1->fitnessValue < c2->fitnessValue || c2->fitnessValue == 0.0) {
        memcpy(c2->genes, c1->genes, sizeof(double) * c1->gene_size);
        c2->fitnessValue = c1->fitnessValue;
    }
}

double geneticAlgorithm(population_t *pop)
{
    int AMOUNT_OF_JOBS = pop->round.AMOUNT_OF_JOBS;
    job_t *jobs = pop->round.jobs;
    chromosomes_t *chromosomes = pop->chromosomes;

    map<unsigned int, machine_t *> machines = pop->round.machines;

    // ops
    machine_base_operations_t *machine_ops = pop->operations.machine_ops;
    list_operations_t *list_ops = pop->operations.list_ops;
    job_base_operations_t *job_ops = pop->operations.job_ops;

    // initialize machine_op
    int k;
    clock_t end = clock() + pop->parameters.GENERATIONS * CLOCKS_PER_SEC;
    clock_t c = clock();
    for (k = 0; c < end; ++k, c = clock()) {
        for (int i = 0; i < pop->parameters.AMOUNT_OF_CHROMOSOMES; ++i) {
            for (int j = 0, size = pop->round.AMOUNT_OF_JOBS << 1; j < size;
                 ++j) {
                int integer_part =
                    pop->chromosomes[i].current_chromosome.genes[j] +=
                    pop->chromosomes[i].velocity;
                pop->chromosomes[i].current_chromosome.genes[j] =
                    abs(pop->chromosomes[i].current_chromosome.genes[j] -
                        integer_part);
            }
        }

        for (int i = 0; i < pop->parameters.AMOUNT_OF_CHROMOSOMES;
             ++i) {  // for all chromosomes
            chromosomes[i].current_chromosome.fitnessValue =
                decoding(chromosomes[i].current_chromosome, jobs, machines,
                         machine_ops, list_ops, job_ops, AMOUNT_OF_JOBS);
        }

        // find the global best
        double best_fitness_value = 1000000000000;
        int idx;
        for (int i = 0; i < pop->parameters.AMOUNT_OF_CHROMOSOMES; ++i) {
            if (chromosomes[i].current_chromosome.fitnessValue <
                best_fitness_value) {
                best_fitness_value =
                    chromosomes[i].current_chromosome.fitnessValue;
                idx = i;
            }
        }

        // chromosomes[idx] is the best in this generation and try to set the
        // real global best
        compareAndAssignChromosomes(&chromosomes[idx].current_chromosome,
                                    &pop->global_best_chromosome);
        memcpy(chromosomes[0].history_best_chromosome.genes,
               chromosomes[idx].current_chromosome.genes,
               sizeof(double) * AMOUNT_OF_JOBS);

        // calculate velocity
        for (int i = 0; i < pop->parameters.AMOUNT_OF_CHROMOSOMES; ++i) {
            for (int j = 0, size = AMOUNT_OF_JOBS << 1; j < size; ++j) {
                chromosomes[i].velocity =
                    0.5 * chromosomes[i].velocity +
                    2 * randomDouble() *
                        (chromosomes[0].history_best_chromosome.genes[j] -
                         chromosomes[i].current_chromosome.genes[j]) +
                    2 * randomDouble() *
                        (pop->global_best_chromosome.genes[j] -
                         chromosomes[i].current_chromosome.genes[j]);
            }
        }


        // printf("%d,%f\n",k, pop->global_best_chromosome.fitnessValue);
    }

    return pop->global_best_chromosome.fitnessValue;
}
