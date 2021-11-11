#include <gtest/gtest.h>
#include <ctime>
#include <string>

#define private public
#define protected public

#include "include/entities.h"
#include "include/entity.h"
#include "include/lot.h"
#include "include/lots.h"
#include "include/machine.h"
#include "include/machines.h"

#undef private
#undef protected

using namespace std;

struct test_case_t {
    int _number_of_resources;
    map<string, int> _case;
    map<string, int> _ans;
};

class test_machines_t_distributeAResource : public testing::Test
{
public:
    vector<test_case_t> cases;
    machines_t *machines;

    void SetUp() override;
    void TearDown() override;
};

void test_machines_t_distributeAResource::SetUp()
{
    machines = nullptr;
    machines = new machines_t();
    if (machines == nullptr)
        exit(EXIT_FAILURE);
    cases.push_back(test_case_t{100,
                                {
                                    {"g1"s, 40},
                                    {"g2"s, 60},
                                },
                                {{"g1"s, 40}, {"g2"s, 60}}});

    // cases.push_back((struct test_case_t){
    //     ._number_of_resources = 50,
    //     ._case = map<string, int>({
    //         {"g1", 40},
    //         {"g2", 60},
    //     }),
    //     ._ans = map<string, int>({{"g1", 20}, {"g2", 30}})});

    // cases.push_back((struct test_case_t){
    //     ._number_of_resources = 48,
    //     ._case = map<string, int>({
    //         {"g1", 40},
    //         {"g2", 60},
    //     }),
    //     ._ans = map<string, int>({{"g1", 19}, {"g2", 29}})});

    // cases.push_back((struct test_case_t){
    //     ._number_of_resources = 29,
    //     ._case = map<string, int>({
    //         {"g1", 43},
    //         {"g2", 31},
    //     }),
    //     ._ans = map<string, int>({{"g1", 17}, {"g2", 12}})});

    // cases.push_back(
    //     (struct test_case_t){._number_of_resources = 3,
    //                          ._case = map<string, int>({
    //                              {"g1", 43},
    //                              {"g2", 31},
    //                          }),
    //                          ._ans = map<string, int>({{"g1", 2}, {"g2",
    //                          1}})});

    // cases.push_back(
    //     (struct test_case_t){._number_of_resources = 2,
    //                          ._case = map<string, int>({
    //                              {"g1", 43},
    //                              {"g2", 31},
    //                          }),
    //                          ._ans = map<string, int>({{"g1", 1}, {"g2",
    //                          1}})});

    // cases.push_back(
    //     (struct test_case_t){._number_of_resources = 1,
    //                          ._case = map<string, int>({
    //                              {"g1", 43},
    //                              {"g2", 31},
    //                          }),
    //                          ._ans = map<string, int>({{"g1", 0}, {"g2",
    //                          1}})});
    // cases.push_back(
    //     (struct test_case_t){._number_of_resources = 6,
    //                          ._case = map<string, int>({
    //                              {"g1", 1},
    //                              {"g2", 12},
    //                          }),
    //                          ._ans = map<string, int>({{"g1", 1}, {"g2",
    //                          5}})});



    // cases.push_back((struct test_case_t){
    //     ._number_of_resources = 5,
    //     ._case = map<string, int>({{"g1", 43}, {"g2", 31}, {"g3", 25}}),
    //     ._ans = map<string, int>({{"g1", 3}, {"g2", 1}, {"g3", 1}})});

    // cases.push_back((struct test_case_t){
    //     ._number_of_resources = 3,
    //     ._case = map<string, int>({{"g1", 43}, {"g2", 1}, {"g3", 5}}),
    //     ._ans = map<string, int>({{"g1", 1}, {"g2", 1}, {"g3", 1}})});
}

void test_machines_t_distributeAResource::TearDown()
{
    delete machines;
}

TEST_F(test_machines_t_distributeAResource, test_correctness)
{
    for (unsigned int i = 0; i < cases.size(); ++i) {
        map<string, int> ans = machines->_distributeAResource(
            cases[i]._number_of_resources, cases[i]._case);
        // for (map<string, int>::iterator it = cases[i]._ans.begin();
        //      it != cases[i]._ans.end(); ++it) {
        //     EXPECT_NO_THROW(ans.at(it->first));
        //     EXPECT_EQ(ans[it->first], cases[i]._ans[it->first]);
        // }
    }
}
