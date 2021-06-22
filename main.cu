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
#include <include/condition_card.h>
#include <include/csv.h>
#include <include/da.h>
#include <include/entity.h>
#include <include/infra.h>
#include <include/lot.h>
#include <include/route.h>
#include <include/population.h>
#include <include/chromosome.h>

#include <cuda.h>


using namespace std;

// round_t createARound(vector<lot_group_t> group, machines_t &machines, ancillary_resources_t & tools, ancillary_resources_t & wires);
// 
// 
// void initializeARound(round_t * r);
// 
// void initializePopulation(population_t *pop, machines_t & machines, ancillary_resources_t & tools, ancillary_resources_t & wires);
// 
// 
// void geneticAlgorithm(population_t * pop);


task_t createTaskFromLotGroups(vector<lot_group_t> groups, ancillary_resources_t & tools, ancillary_resources_t & wires, machines_t & machines);

int main(int argc, const char *argv[])
{
    srand(time(NULL));
    csv_t lot_csv("lots.csv", "r", true, true);
    lot_csv.setHeaders(map<string, string>({
                    {"route", "route"},
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
                    {"customer", "customer"}
                }));
    lot_csv.trim(" ");
    vector<lot_t> allots;
    for(int i = 0, nrows = lot_csv.nrows(); i < nrows; ++i){
        allots.push_back(lot_t(lot_csv.getElements(i)));
    } 

    lots_t lots;
    lots.addLots(allots);

        
    // lots.addLots(createLots("WipOutPlanTime_.csv", "product_find_process_id.csv",
    //                "process_find_lot_size_and_entity.csv", "fcst.csv",
    //                "routelist.csv", "newqueue_time.csv",
    //                "BOM List_20210521.csv", "Process find heatblock.csv",
    //                "EMS Heatblock data.csv", "GW Inventory.csv"));


    ancillary_resources_t tools(lots.amountOfTools());
    ancillary_resources_t wires(lots.amountOfWires());

    csv_t machine_csv("machines.csv", "r", true, true);
    machine_csv.trim(" ");
    machine_csv.setHeaders(map<string, string>({{"entity", "ENTITY"},
                                                {"model", "MODEL"},
                                                {"recover_time", "OUTPLAN"}}));

    csv_t location_csv("locations.csv", "r", true, true);
    location_csv.trim(" ");
    location_csv.setHeaders(
        map<string, string>({{"entity", "Entity"}, {"location", "Location"}}));


    char *text = strdup("2020/12/19 10:50");
    entities_t entities(text);
    entities.addMachines(machine_csv, location_csv);
    machines_t machines;
    machines.addMachines(entities.getAllEntity());

    vector<vector<lot_group_t> > round_groups = lots.rounds(entities);

    task_t t = createTaskFromLotGroups(round_groups[0], tools, wires, machines);
    
    population_t pop{
        .no = 0,
        .parameters = {
            .AMOUNT_OF_CHROMOSOMES = 100,
            .AMOUNT_OF_R_CHROMOSOMES = 200,
            .EVOLUTION_RATE = 0.2,
            .SELECTION_RATE = 0.3,
            .GENERATIONS = 10,
        },
        .task = t
    };
    
    initializePopulation(&pop);
    return 0;
}


task_t createTaskFromLotGroups(vector<lot_group_t> groups, ancillary_resources_t & tools, ancillary_resources_t & wires, machines_t & machines){
    // setup jobs 
    int AMOUNT_OF_JOBS = 0;
    int AMOUNT_OF_MACHINES = 0;
    int k = 0;
    iter(groups, i){
        AMOUNT_OF_JOBS += groups[i].lots.size();
        AMOUNT_OF_MACHINES += groups[i].machine_amount;
    }
    
    job_t * jobs;
    cudaCheck(cudaMallocHost(&jobs, sizeof(job_t)*AMOUNT_OF_JOBS), "cudasMallocHost for jobs");
    
    // setup jobs data
    iter(groups, i){
        iter(groups[i].lots, j){
            jobs[k] = groups[i].lots[j]->job(); 
            ++k;
        }
    }

    // setup process time
    int *size_of_process_times;
    cudaCheck(cudaMallocHost(&size_of_process_times, sizeof(int)*AMOUNT_OF_JOBS), "cudaMallocHost for size_of_process_times");
    process_time_t **process_times;
    cudaCheck(cudaMallocHost(&process_times, sizeof(process_time_t*)*AMOUNT_OF_JOBS), "cudaMallocHost for process times");
    process_time_t *process_times_entry;
    k = 0;
    iter(groups, i){
        iter(groups[i].lots, j){
            map<string, double> mpts = groups[i].lots[j]->getEntitiyProcessTime();
            
            cudaCheck(cudaMallocHost(&process_times_entry, sizeof(process_time_t)*mpts.size()), "cudaMallocHost for process time entry");
            int l = 0;
            for(map<string, double>::iterator it = mpts.begin(); it != mpts.end(); it++){
                process_times_entry[l] = process_time_t{
                    .machine_no = convertEntityNameToUInt(it->first),
                    .process_time = it->second
                };  
                ++l;
            }
            process_times[k] = process_times_entry;
            size_of_process_times[k] = mpts.size();
            ++k;
        }
    }

    // setup machines
    map<string, machine_t *> allmachines = machines.getMachines();
    machine_t * ms;
    cudaCheck(cudaMallocHost(&ms, sizeof(machine_t)*AMOUNT_OF_MACHINES), "cudaMallocHost for machines");
    k = 0;
    iter(groups, i){
        iter(groups[i].entities, j){
            string ent_name = groups[i].entities[j]->entity_name;
            ms[k] = *(allmachines[ent_name]);
            ++k;
        }
    }

    // setup tools and machines
    tool_t * ts;
    wire_t * ws;
    cudaCheck(cudaMallocHost(&ts, sizeof(tool_t)*AMOUNT_OF_MACHINES), "cudaMallocHost for tools");
    cudaCheck(cudaMallocHost(&ws, sizeof(wire_t)*AMOUNT_OF_MACHINES), "cudaMallocHost for wires");
    k = 0;
    iter(groups, i){
        vector<tool_t *> v_ts = tools.aRound(groups[i].tool_name, groups[i].machine_amount);
        vector<wire_t *> v_ws = wires.aRound(groups[i].wire_name, groups[i].machine_amount);
        for(int j = 0; j < groups[i].machine_amount; ++j){
            ts[k] = *(v_ts[j]);
            ws[k] = *(v_ws[j]);
            ++k;
        }

    }
    return task_t {
        .AMOUNT_OF_JOBS = AMOUNT_OF_JOBS,
        .AMOUNT_OF_MACHINES = AMOUNT_OF_MACHINES,
        .jobs = jobs,
        .machines = ms,
        .tools = ts,
        .wires = ws,
        .process_times = process_times,
        .size_of_process_times = size_of_process_times
    };
}
