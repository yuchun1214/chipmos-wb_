#include <gtest/gtest.h>

#define private public
#define protected public

#include "include/machines.h"

#undef private
#undef protected

#include <map>
#include <string>
#include <vector>

using namespace std;

struct test_case_t {
    map<string, string> input;
    bool out;
};

class test_machines_t_canJobRunOnTheMachine : public testing::Test
{
public:
    map<string, string> case1;
    map<string, string> case2;
    map<string, string> case3;
    job_t *job;
    machines_t *machines;
    machine_t *machine;

    void createTestingObjects(job_t **job,
                              machine_t **machine,
                              machines_t **machines,
                              map<string, string> cs);

    void SetUp() override;
    void TearDown() override;
};

void test_machines_t_canJobRunOnTheMachine::createTestingObjects(
    job_t **job,
    machine_t **machine,
    machines_t **machines,
    map<string, string> cs)
{
    job_t *j = new job_t();
    j->base.job_info = stringToInfo(cs["lot_number"]);

    machine_t *m = new machine_t();
    m->model_name = stringToInfo(cs["machine_model"]);
    m->location = stringToInfo(cs["machine_location"]);

    machines_t *ms = new machines_t();
    char *text = strdup(cs["can_run_models"].c_str());
    vector<string> can_run_models = split(text, ',');
    free(text);

    text = strdup(cs["process_time"].c_str());
    vector<string> process_time_str = split(text, ',');
    free(text);

    map<string, double> process_times;
    foreach (process_time_str, i) {
        process_times[can_run_models[i]] = stod(process_time_str[i]);
    }

    text = strdup(cs["can_run_location"].c_str());
    vector<string> can_run_location = split(text, ',');
    ms->_job_can_run_locations[cs["lot_number"]] = can_run_location;
    ms->_job_process_times[cs["lot_number"]] = process_times;

    if (cs["is_automotive"].compare("Y") == 0) {
        ms->setAutomotiveLotNumber(set<string>{cs["lot_number"]});
    }

    *job = j;
    *machine = m;
    *machines = ms;
}

void test_machines_t_canJobRunOnTheMachine::TearDown()
{
    delete job;
    delete machine;
    delete machines;
}

void test_machines_t_canJobRunOnTheMachine::SetUp()
{
    case1 = map<string, string>(
        {{"lot_number", "L1"},
         {"cust", "MTI"},
         {"machine_location", "TA"},
         {"machine_model", "UTC1000"},
         {"machine_name", "BB001"} {"can_run_location", "TA,TB"},
         {"can_run_models", "UTC1000,UTC2000,UTC3000"},
         {"process_time", "1,2,3"},
         {"is_automotive", "Y"}});

    // case2 = map<string, string>({{"lot_number", "L1"},
    //                              {"machine_location", "TA"},
    //                              {"machine_model", "UTC5000"},
    //                              {"can_run_location", "TA,TB"},
    //                              {"can_run_models",
    //                              "UTC1000,UTC2000,UTC3000"},
    //                              {"process_time", "1,2,3"}});
    // case3 = map<string, string>({{"lot_number", "L1"},
    //                              {"machine_location", "TC"},
    //                              {"machine_model", "UTC5000"},
    //                              {"can_run_location", "TA,TB"},
    //                              {"can_run_models",
    //                              "UTC1000,UTC2000,UTC3000"},
    //                              {"process_time", "1,2,3"}});
}

TEST_F(test_machines_t_canJobRunOnTheMachine, test_correctness1)
{
    createTestingObjects(&job, &machine, &machines, case1);
    EXPECT_TRUE(machines->_canJobRunOnTheMachine(job, machine));
}

// TEST_F(test_machines_t_canJobRunOnTheMachine, test_correctness2)
// {
//     createTestingObjects(&job, &machine, &machines, case2);
//     EXPECT_FALSE(machines->_canJobRunOnTheMachine(job, machine));
// }
//
// TEST_F(test_machines_t_canJobRunOnTheMachine, test_correctness3)
// {
//     createTestingObjects(&job, &machine, &machines, case3);
//     EXPECT_FALSE(machines->_canJobRunOnTheMachine(job, machine));
// }
