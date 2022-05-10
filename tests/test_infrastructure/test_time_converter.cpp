#include <gtest/gtest.h>
#include <iostream>
#include <string>

#define private public
#include "include/time_converter.h"
#undef private

using namespace std;

struct string_to_time_test_case {
    string input;
    time_t ans;
    string base_time;
    time_t difference;
};

ostream &operator<<(ostream &out, string_to_time_test_case cs)
{
    return out << "(" << cs.input << "," << cs.ans << "," << cs.base_time << ","
               << cs.difference << ")";
}

class test_string_to_time_t
    : public testing::TestWithParam<string_to_time_test_case>
{
};

TEST_P(test_string_to_time_t, test_correctness)
{
    auto cs = GetParam();
    EXPECT_EQ(cs.ans, timeConverter()(cs.input))
        << "Failed with case : " << cs << endl;
}

TEST_P(test_string_to_time_t, test_ctor_setup_base_time_by_str)
{
    auto cs = GetParam();
    timeConverter tc(cs.input);
    EXPECT_EQ(cs.ans, tc._base_time);
    EXPECT_EQ(tc.getBaseTime(), tc._base_time);
}

TEST_P(test_string_to_time_t, test_ctor_setup_base_time_by_value)
{
    auto cs = GetParam();
    timeConverter tc(cs.ans);
    EXPECT_EQ(cs.ans, tc._base_time);
    EXPECT_EQ(tc.getBaseTime(), tc._base_time);
}

TEST_P(test_string_to_time_t, test_time_difference)
{
    auto cs = GetParam();
    EXPECT_EQ(cs.difference, timeConverter(cs.base_time)(cs.input));
    EXPECT_EQ((-1) * cs.difference, timeConverter(cs.input)(cs.base_time));
    EXPECT_EQ((-1) * cs.difference, timeConverter(cs.ans)(cs.base_time));
}

TEST_P(test_string_to_time_t, test_set_base_time)
{
    auto cs = GetParam();
    timeConverter tc;
    tc.setBaseTime(cs.input);
    EXPECT_EQ(cs.ans, tc._base_time);
}

TEST_P(test_string_to_time_t, test_time_difference_operator)
{
    auto cs = GetParam();
    timeConverter tc_base(cs.base_time);
    timeConverter tc_input(cs.input);
    EXPECT_EQ(cs.difference, tc_base - tc_input);
}

TEST_P(test_string_to_time_t, test_time_difference_operator2)
{
    auto cs = GetParam();
    EXPECT_EQ(cs.difference,
              timeConverter(cs.base_time) - timeConverter(cs.input));
}


INSTANTIATE_TEST_SUITE_P(
    string_to_time,
    test_string_to_time_t,
    testing::Values(
        string_to_time_test_case{"21-12-31 20:45", 1640954700, "22-03-03 03:16",
                                 -5293860},
        string_to_time_test_case{"21-12-31 20:45", 1640954700, "22-03-03 03:16",
                                 -5293860},
        string_to_time_test_case{"22-03-03 03:16", 1646248560, "21-12-31 20:45",
                                 5293860},
        string_to_time_test_case{"21-12-31 20:45", 1640954700,
                                 "2022/03/03 03:16", -5293860},
        string_to_time_test_case{"22-03-03 03:16", 1646248560,
                                 "2021/12/31 20:45", 5293860},
        string_to_time_test_case{"22-3-03 03:16", 1646248560,
                                 "2021/11/23 20:45", 8577060},
        string_to_time_test_case{"22-3-03 03:16", 1646248560, "21-11-23 20:45",
                                 8577060},
        string_to_time_test_case{"22-03-3 03:16", 1646248560, "", 1646248560},
        string_to_time_test_case{"22-3-3 03:16", 1646248560, "", 1646248560},
        string_to_time_test_case{"20-02-29 03:16", 1582917360, "", 1582917360},
        string_to_time_test_case{"22-02-08 19:25", 1644319500, "", 1644319500},
        string_to_time_test_case{"22-02-08 19:25:01", 1644319501, "",
                                 1644319501},
        string_to_time_test_case{"22-2-8 19:25:01", 1644319501, "", 1644319501},
        string_to_time_test_case{"22-02-8 19:25:01", 1644319501, "",
                                 1644319501},
        string_to_time_test_case{"22-2-08 19:25:01", 1644319501, "",
                                 1644319501},
        string_to_time_test_case{"2022/2/8 05:40:01", 1644270001, "",
                                 1644270001},
        string_to_time_test_case{"2022/02/08 05:40:01", 1644270001, "",
                                 1644270001},
        string_to_time_test_case{"2022/2/08 05:40:01", 1644270001, "",
                                 1644270001},
        string_to_time_test_case{"2022/02/8 05:40:01", 1644270001, "",
                                 1644270001},
        string_to_time_test_case{"2021/11/23 20:45", 1637671500, "",
                                 1637671500},
        string_to_time_test_case{"2021/12/31 20:45", 1640954700, "",
                                 1640954700},
        string_to_time_test_case{"2021/5/1 20:45", 1619873100, "", 1619873100},
        string_to_time_test_case{"2021/05/1 20:45", 1619873100, "", 1619873100},
        string_to_time_test_case{"2021/5/01 20:45", 1619873100, "", 1619873100},
        string_to_time_test_case{"2021/05/01 20:45", 1619873100, "",
                                 1619873100},

        string_to_time_test_case{"22-03-03 03:03", 1646247780, "", 1646247780},
        string_to_time_test_case{"22-03-03 3:03", 1646247780, "", 1646247780},
        string_to_time_test_case{"22-03-03 03:3", 1646247780, "", 1646247780},
        string_to_time_test_case{"22-03-03 3:3", 1646247780, "", 1646247780},
        string_to_time_test_case{"22-03-03 3:3:3", 1646247783, "", 1646247783},

        string_to_time_test_case{"2022/03/03 03:03", 1646247780, "",
                                 1646247780},
        string_to_time_test_case{"2022/03/03 3:03", 1646247780, "", 1646247780},
        string_to_time_test_case{"2022/03/03 03:3", 1646247780, "", 1646247780},
        string_to_time_test_case{"2022/03/03 3:3", 1646247780, "", 1646247780},
        string_to_time_test_case{"2022/03/03 3:3:3", 1646247783, "",
                                 1646247783},
        string_to_time_test_case{"2022/3/14 1:59:26", 1647194366, "",
                                 1647194366},
        string_to_time_test_case{"2022/3/14 15:9:26", 1647241766, "",
                                 1647241766},

        string_to_time_test_case{"2022/05/01", 1651334400, "", 1651334400},
        string_to_time_test_case{"2022/5/01", 1651334400, "", 1651334400},
        string_to_time_test_case{"2022/05/1", 1651334400, "", 1651334400},
        string_to_time_test_case{"2022/5/1", 1651334400, "", 1651334400},

        string_to_time_test_case{"22-05-01", 1651334400, "", 1651334400},
        string_to_time_test_case{"22-5-01", 1651334400, "", 1651334400},
        string_to_time_test_case{"22-05-1", 1651334400, "", 1651334400},
        string_to_time_test_case{"22-5-1", 1651334400, "", 1651334400},
        string_to_time_test_case{"22-05-10 02:06:57", 1652119617, "",
                                 1652119617},
        string_to_time_test_case{"22-05-10 02:06", 1652119560, "", 1652119560},
        string_to_time_test_case{"2022/05/10 02:06:57", 1652119617, "",
                                 1652119617},
        string_to_time_test_case{"2022/05/10 02:06", 1652119560, "",
                                 1652119560},



        string_to_time_test_case{"22-005-01", 0, "", 0},
        string_to_time_test_case{"", 0, "", 0},
        string_to_time_test_case{"2021\\05\\01 20:45", 0, "", 0},
        string_to_time_test_case{"2021-05-1", 0, "", 0}));
