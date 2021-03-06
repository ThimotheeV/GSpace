#define CATCH_CONFIG_MAIN unit_test::CoalescenceTableTest
#include "catch.hpp"
#include <iostream>

#include "coalescence_table.hpp"
#include "settings.hpp"

namespace unit_test
{
    struct CoalescenceTableTest
    {
        // unit_test::CoalescenceTableTest() {}
        // ~unit_test::CoalescenceTableTest() {}

        static std::vector<coa_event_t> &getCoalescenceTable(coa_table_c &c)
        {
            return c.Coalescence_table;
        }

        static std::vector<std::tuple<std::vector<int>, double, int>> &getAncestryTree(tree_gen_c &tg)
        {
            return tg.Ancestry_tree;
        }

        static coa_table_c sort_table_I(coa_table_c &c)
        {
            return tree_gen_c::sort_table_I(c);
        }

        static coa_table_c sort_table_R(coa_table_c &c)
        {
            return tree_gen_c::sort_table_R(c);
        }
    };
} // namespace unit_test

//Caution : exemple are the same as the article of Kellerher but with all the node num  -1

/**************************************************************/
/*                 Coalescence table                          */
/**************************************************************/

TEST_CASE_METHOD(unit_test::CoalescenceTableTest, "coalescence_table_test")
{
    SECTION("get_coa_event")
    {
        coa_event_t event(2, 10, 4, {2, 3}, 0.071, -1);

        REQUIRE(coa_table_c::get_left_brkpt_l(event) == 2);
        REQUIRE(coa_table_c::get_right_brkpt_r(event) == 10);
        REQUIRE(coa_table_c::get_num_node_u(event) == 4);
        REQUIRE((coa_table_c::get_childs_node_c(event))[0] == 2);
        REQUIRE((coa_table_c::get_childs_node_c(event))[1] == 3);
        REQUIRE(coa_table_c::get_time_t(event) == 0.071);
    }

    SECTION("set_coalescance_table")
    {
        coa_table_c table_event;

        table_event.set_coa_event(2, 10, 4, {2, 3}, 0.071, -1);
        table_event.set_coa_event(0, 2, 5, {1, 3}, 0.090, -1);
        table_event.set_coa_event(2, 10, 5, {1, 4}, 0.090, 3);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[2]) == 2);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[2]) == 10);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[2]) == 5);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[0] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[1] == 4);
        REQUIRE(coa_table_c::get_time_t(table_event.get_coa_table()[2]) == 0.090);
        REQUIRE(coa_table_c::get_gen(table_event.get_coa_table()[2]) == 3);
    }

    SECTION("size_coalescance_table")
    {
        coa_table_c table_event;

        table_event.set_coa_event(2, 10, 4, {2, 3}, 0.071, -1);
        table_event.set_coa_event(0, 2, 5, {1, 3}, 0.090, -1);
        table_event.set_coa_event(2, 10, 5, {1, 4}, 0.090, -1);

        REQUIRE(table_event.size() == 3);
    }

    SECTION("sort_table_I_normal_order_generate_tree_algorithm")
    {
        coa_event_t event1(2, 10, 4, {2, 3}, 0.071, -1);
        coa_event_t event2(0, 2, 5, {1, 3}, 0.090, -1);
        coa_event_t event3(2, 10, 5, {1, 4}, 0.090, -1);
        coa_event_t event4(0, 7, 6, {0, 5}, 0.170, -1);
        coa_event_t event5(7, 10, 7, {0, 5}, 0.202, -1);
        coa_event_t event6(0, 2, 8, {2, 6}, 0.253, -1);

        coa_table_c table_event;

        getCoalescenceTable(table_event).push_back(event1);
        getCoalescenceTable(table_event).push_back(event2);
        getCoalescenceTable(table_event).push_back(event3);
        getCoalescenceTable(table_event).push_back(event4);
        getCoalescenceTable(table_event).push_back(event5);
        getCoalescenceTable(table_event).push_back(event6);

        table_event = sort_table_I(table_event);
        REQUIRE(coa_table_c::get_left_brkpt_l(table_event[0]) == 0);
        REQUIRE(coa_table_c::get_time_t(table_event[0]) == 0.090);
        REQUIRE(coa_table_c::get_left_brkpt_l(table_event[2]) == 0);
        REQUIRE(coa_table_c::get_time_t(table_event[2]) == 0.253);
        REQUIRE(coa_table_c::get_left_brkpt_l(table_event[5]) == 7);
        REQUIRE(coa_table_c::get_time_t(table_event[5]) == 0.202);
    }

    SECTION("sort_table_R_reverse_order_generate_tree_algorithm")
    {
        coa_event_t event1(2, 10, 4, {2, 3}, 0.071, -1);
        coa_event_t event2(0, 2, 5, {1, 3}, 0.090, -1);
        coa_event_t event3(2, 10, 5, {1, 4}, 0.090, -1);
        coa_event_t event4(0, 7, 6, {0, 5}, 0.170, -1);
        coa_event_t event5(7, 10, 7, {0, 5}, 0.202, -1);
        coa_event_t event6(0, 2, 8, {2, 6}, 0.253, -1);

        coa_table_c table_event;

        getCoalescenceTable(table_event).push_back(event1);
        getCoalescenceTable(table_event).push_back(event2);
        getCoalescenceTable(table_event).push_back(event3);
        getCoalescenceTable(table_event).push_back(event4);
        getCoalescenceTable(table_event).push_back(event5);
        getCoalescenceTable(table_event).push_back(event6);

        table_event = sort_table_R(table_event);

        REQUIRE(coa_table_c::get_right_brkpt_r(table_event[0]) == 2);
        REQUIRE(coa_table_c::get_time_t(table_event[0]) == 0.253);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event[3]) == 10);
        REQUIRE(coa_table_c::get_time_t(table_event[3]) == 0.202);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event[5]) == 10);
        REQUIRE(coa_table_c::get_time_t(table_event[5]) == 0.071);
    }

    SECTION("multi_coa_coalescance_table")
    {
        coa_table_c table_event;

        table_event.set_coa_event(0, 1, 4, {3, 2}, 0, 1);
        table_event.set_coa_event(0, 1, 5, {4, 1}, 0, 1);

        table_event.group_multi_coa(0);

        table_event.set_coa_event(0, 1, 6, {5, 0}, 0, 1);

        table_event.group_multi_coa(0);

        REQUIRE(table_event.size() == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[0]) == 6);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[1] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[2] == 3);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[3] == 2);
        REQUIRE(coa_table_c::get_time_t(table_event.get_coa_table()[0]) == 0);
        REQUIRE(coa_table_c::get_gen(table_event.get_coa_table()[0]) == 1);
    }

    SECTION("group_multi_coa_coalescance_table")
    {
        coa_table_c table_event;

        table_event.set_coa_event(2, 10, 4, {2, 3}, 0.071, 1);
        table_event.set_coa_event(0, 2, 5, {1, 3}, 0.090, 2);
        table_event.set_coa_event(2, 8, 6, {1, 4}, 0.090, 2);
        table_event.set_coa_event(2, 8, 7, {0, 6}, 0.090, 2);

        table_event.group_multi_coa(1);

        REQUIRE(table_event.size() == 3);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[1]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[1]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[1]) == 5);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[0] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[1] == 3);
        REQUIRE(coa_table_c::get_time_t(table_event.get_coa_table()[1]) == 0.090);
        REQUIRE(coa_table_c::get_gen(table_event.get_coa_table()[1]) == 2);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[2]) == 2);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[2]) == 8);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[2]) == 7);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[1] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[2] == 4);
        REQUIRE(coa_table_c::get_time_t(table_event.get_coa_table()[2]) == 0.090);
        REQUIRE(coa_table_c::get_gen(table_event.get_coa_table()[2]) == 2);
    }

    SECTION("simple_group_multi_coa_coalescance_table")
    {
        coa_table_c table_event;

        table_event.set_coa_event(0, 1, 10, {4, 7}, 0, 1);
        table_event.set_coa_event(0, 1, 11, {0, 2}, 0, 1);
        table_event.set_coa_event(0, 1, 12, {1, 3}, 0, 1);
        table_event.set_coa_event(0, 1, 13, {12, 9}, 0, 1);

        table_event.group_multi_coa(0);

        REQUIRE(table_event.size() == 3);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[0]) == 10);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[0] == 4);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[1] == 7);
        REQUIRE(coa_table_c::get_time_t(table_event.get_coa_table()[0]) == 0);
        REQUIRE(coa_table_c::get_gen(table_event.get_coa_table()[0]) == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[1]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[1]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[1]) == 11);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[1] == 2);
        REQUIRE(coa_table_c::get_time_t(table_event.get_coa_table()[1]) == 0);
        REQUIRE(coa_table_c::get_gen(table_event.get_coa_table()[1]) == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[2]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[2]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[2]) == 13);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[0] == 9);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[1] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[2] == 3);
        REQUIRE(coa_table_c::get_time_t(table_event.get_coa_table()[2]) == 0);
        REQUIRE(coa_table_c::get_gen(table_event.get_coa_table()[2]) == 1);
    }

    SECTION("intercal_other_coa_group_multi_coa_coalescance_table")
    {
        coa_table_c table_event;

        table_event.set_coa_event(0, 1, 100, {0, 2}, 0, 1);
        table_event.set_coa_event(0, 1, 101, {100, 3}, 0, 1);
        table_event.set_coa_event(0, 1, 102, {101, 5}, 0, 1);

        table_event.group_multi_coa(0);

        table_event.set_coa_event(0, 1, 103, {102, 7}, 0, 1);

        table_event.group_multi_coa(0);

        table_event.set_coa_event(0, 1, 104, {103, 83}, 0, 1);

        table_event.group_multi_coa(0);

        table_event.set_coa_event(0, 1, 105, {4, 94}, 0, 1);
        table_event.set_coa_event(0, 1, 106, {104, 99}, 0, 1);

        table_event.group_multi_coa(0);

        REQUIRE(table_event.size() == 2);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[0]) == 105);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[0] == 4);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[1] == 94);
        REQUIRE(coa_table_c::get_time_t(table_event.get_coa_table()[0]) == 0);
        REQUIRE(coa_table_c::get_gen(table_event.get_coa_table()[0]) == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[1]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[1]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[1]) == 106);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[0] == 99);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[1] == 83);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[2] == 7);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[3] == 5);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[4] == 3);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[5] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[6] == 2);
        REQUIRE(coa_table_c::get_time_t(table_event.get_coa_table()[1]) == 0);
        REQUIRE(coa_table_c::get_gen(table_event.get_coa_table()[1]) == 1);
    }

    SECTION("group_multi_coa_partial_clean_coalescance_table")
    {
        coa_table_c table_event;

        table_event.set_coa_event(0, 3, 4, {2, 3}, 0.090, 1);
        table_event.set_coa_event(0, 3, 5, {1, 4}, 0.090, 1);

        table_event.group_multi_coa(0);

        REQUIRE(table_event.size() == 1);

        table_event.set_coa_event(2, 3, 6, {0, 5}, 0.090, 1);

        table_event.group_multi_coa(0);

        REQUIRE(table_event.size() == 2);
        //Three coalescences in part 0 to 2
        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[0]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[0]) == 5);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[0] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[1] == 2);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[2] == 3);
        REQUIRE(coa_table_c::get_time_t(table_event.get_coa_table()[0]) == 0.090);
        REQUIRE(coa_table_c::get_gen(table_event.get_coa_table()[0]) == 1);
        //Four coalescences in part 2 to 3
        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[1]) == 2);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[1]) == 3);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[1]) == 6);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[1] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[2] == 2);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[3] == 3);
        REQUIRE(coa_table_c::get_time_t(table_event.get_coa_table()[1]) == 0.090);
        REQUIRE(coa_table_c::get_gen(table_event.get_coa_table()[1]) == 1);
    }

    SECTION("group_multi_coa_complexe_coalescance_table")
    {
        coa_table_c table_event;

        table_event.set_coa_event(0, 3, 5, {1, 4}, 0.077, 0);
        table_event.set_coa_event(2, 3, 6, {0, 5}, 0.077, 0);
        
        table_event.set_coa_event(5, 9, 4, {2, 1}, 0.090, 1);
        table_event.set_coa_event(6, 8, 5, {4, 3}, 0.090, 1);
        table_event.set_coa_event(4, 5, 6, {1, 0}, 0.090, 1);
        table_event.set_coa_event(5, 6, 7, {4, 0}, 0.090, 1);
        table_event.set_coa_event(6, 8, 8, {5, 0}, 0.090, 1);
        table_event.set_coa_event(8, 9, 9, {4, 0}, 0.090, 1);
        table_event.set_coa_event(9, 10, 10, {2, 0}, 0.090, 1);

        table_event.group_multi_coa(2);

        REQUIRE(table_event.size() == 7);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[2]) == 4);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[2]) == 5);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[2]) == 6);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[0] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[1] == 0);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[3]) == 5);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[3]) == 6);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[3]) == 7);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[3]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[3]))[1] == 2);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[3]))[2] == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[4]) == 6);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[4]) == 8);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[4]) == 8);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[4]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[4]))[1] == 3);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[4]))[2] == 2);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[4]))[3] == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[5]) == 8);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[5]) == 9);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[5]) == 9);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[5]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[5]))[1] == 2);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[5]))[2] == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[6]) == 9);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[6]) == 10);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[6]) == 10);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[6]))[0] == 2);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[6]))[1] == 0);
    }

    /************EDGES CASES******************/

    SECTION("NO_group_multi_coalescance_table")
    {
        coa_table_c table_event;

        table_event.set_coa_event(0, 3, 12, {1, 4}, 0.077, 0);
        table_event.set_coa_event(2, 3, 13, {0, 5}, 0.077, 0);
        table_event.set_coa_event(5, 9, 14, {2, 1}, 0.090, 1);
        table_event.set_coa_event(6, 8, 15, {4, 3}, 0.090, 1);
        table_event.set_coa_event(4, 5, 16, {1, 0}, 0.090, 1);
        table_event.set_coa_event(5, 6, 17, {4, 0}, 0.090, 1);
        table_event.set_coa_event(6, 8, 18, {5, 0}, 0.090, 1);
        table_event.set_coa_event(8, 9, 19, {4, 0}, 0.090, 1);
        table_event.set_coa_event(9, 10, 20, {2, 0}, 0.090, 1);

        table_event.group_multi_coa(2);

        REQUIRE(table_event.size() == 9);
    }

    SECTION("group_multi_coa_strange_case_coalescance_table")
    {
        coa_table_c table_event;

        table_event.set_coa_event(0, 10, 15, {11, 12}, 0, 25);
        table_event.set_coa_event(3, 6, 16, {15, 13}, 0, 25);
        table_event.set_coa_event(6, 9, 16, {15, 14}, 0, 25);

        table_event.group_multi_coa(0);

        REQUIRE(table_event.size() == 4);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[0]) == 3);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[0]) == 15);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[0] == 11);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[1] == 12);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[1]) == 3);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[1]) == 6);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[1]) == 16);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[0] == 13);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[1] == 11);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[2] == 12);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[2]) == 6);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[2]) == 9);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[2]) == 16);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[0] == 14);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[1] == 11);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[2] == 12);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[3]) == 9);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[3]) == 10);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[3]) == 15);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[3]))[0] == 11);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[3]))[1] == 12);
    }

    SECTION("group_multi_coa_strange_case_II_coalescance_table")
    {
        coa_table_c table_event;

        table_event.set_coa_event(2, 4, 57, {3, 17}, 0, 19);
        table_event.set_coa_event(2, 3, 58, {57, 9}, 0, 19);

        table_event.group_multi_coa(0);

        REQUIRE(table_event.size() == 2);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[0]) == 2);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[0]) == 3);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[0]) == 58);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[0] == 9);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[1] == 3);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[2] == 17);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[1]) == 3);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[1]) == 4);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[1]) == 57);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[0] == 3);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[1] == 17);
    }

    //Special case du to reshape_segs in multi coa process.
    //Not a problem for GSpace simu but need be handle in algorithms, quick patch in add_new_chrom
    SECTION("group_multi_coa_strange_case_III_coalescance_table")
    {
        coa_table_c table_event;

        table_event.set_coa_event(4, 5, 1008, {120, 554}, 0, 1);
        table_event.set_coa_event(5, 6, 1008, {554, 120, 200}, 0, 1);
        table_event.set_coa_event(4, 6, 1009, {1008, 791}, 0, 1);

        table_event.group_multi_coa(0);

        REQUIRE(table_event.size() == 2);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[0]) == 4);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[0]) == 5);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[0]) == 1009);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[0] == 791);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[1] == 120);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[2] == 554);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[1]) == 5);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[1]) == 6);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[1]) == 1009);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[0] == 791);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[1] == 554);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[2] == 120);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[3] == 200);
    }

    // -.-
    SECTION("group_multi_coa_strange_case_IV_coalescance_table")
    {
        coa_table_c table_event;

        table_event.set_coa_event(6, 8, 1298, {1191, 967}, 0, 3);
        table_event.set_coa_event(8, 9, 1298, {782, 967}, 0, 3);
        table_event.set_coa_event(7, 9, 1299, {1298, 778}, 0, 3);

        table_event.group_multi_coa(0);

        REQUIRE(table_event.size() == 3);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[0]) == 6);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[0]) == 7);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[0]) == 1298);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[0] == 1191);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[1] == 967);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[1]) == 7);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[1]) == 8);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[1]) == 1299);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[0] == 778);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[1] == 1191);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[2] == 967);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[2]) == 8);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[2]) == 9);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[2]) == 1299);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[0] == 778);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[1] == 782);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[2] == 967);
    }

    // HAAAAAAAAAAAAAAAAAAAAA !!!!!!!!!!!!!!!!!!!!!
    SECTION("group_multi_coa_strange_case_V_coalescance_table")
    {
        coa_table_c table_event;

        table_event.set_coa_event(0, 2, 2254, {1565, 177}, 0, 6);
        table_event.set_coa_event(2, 4, 2254, {334, 177}, 0, 6);
        table_event.set_coa_event(1, 3, 2256, {2254, 936}, 0, 6);

        table_event.group_multi_coa(0);

        REQUIRE(table_event.size() == 4);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[0]) == 2254);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[0] == 1565);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[0]))[1] == 177);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[1]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[1]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[1]) == 2256);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[0] == 936);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[1] == 1565);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[1]))[2] == 177);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[2]) == 2);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[2]) == 3);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[2]) == 2256);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[0] == 936);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[1] == 334);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[2]))[2] == 177);

        REQUIRE(coa_table_c::get_left_brkpt_l(table_event.get_coa_table()[3]) == 3);
        REQUIRE(coa_table_c::get_right_brkpt_r(table_event.get_coa_table()[3]) == 4);
        REQUIRE(coa_table_c::get_num_node_u(table_event.get_coa_table()[3]) == 2254);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[3]))[0] == 334);
        REQUIRE((coa_table_c::get_childs_node_c(table_event.get_coa_table()[3]))[1] == 177);
    }

    /**************************************************************/
    /*                 Tree Generator                             */
    /**************************************************************/

    SECTION("tree_generator")
    {
        coa_table_c table_event;

        table_event.set_coa_event(2, 10, 4, {2, 3}, 0.071, -1);
        table_event.set_coa_event(0, 2, 5, {1, 3}, 0.090, -1);
        table_event.set_coa_event(2, 10, 5, {1, 4}, 0.090, -1);
        table_event.set_coa_event(0, 7, 6, {0, 5}, 0.170, -1);
        table_event.set_coa_event(7, 10, 7, {0, 5}, 0.202, -1);
        table_event.set_coa_event(0, 2, 8, {2, 6}, 0.253, -1);

        tree_gen_c generator(table_event, 9, 4);

        REQUIRE(coa_table_c::get_num_node_u(generator.I_insert_order_recomb[0]) == 5);
        REQUIRE(coa_table_c::get_num_node_u(generator.R_remove_order_recomb[0]) == 8);
        REQUIRE(generator.Max_u == 9);
        REQUIRE(generator.Leaf_number == 4);
    }

    SECTION("find_MRCA")
    {
        std::vector<std::tuple<std::vector<int>, double, int>> ancestry_tree(9);
        ancestry_tree[0] = std::make_tuple(std::vector<int>(), 0, 0);
        ancestry_tree[1] = std::make_tuple(std::vector<int>(), 0, 0);
        ancestry_tree[2] = std::make_tuple(std::vector<int>(), 0, 0);
        ancestry_tree[3] = std::make_tuple(std::vector<int>(), 0, 0);
        ancestry_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1);
        ancestry_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1);
        ancestry_tree[6] = std::make_tuple(std::vector<int>(), 0, 0);
        ancestry_tree[7] = std::make_tuple(std::vector<int>{4, 5}, 2, -1);
        ancestry_tree[8] = std::make_tuple(std::vector<int>(), 0, 0);

        auto num_node_MRCA = find_MRCA(ancestry_tree);

        REQUIRE(num_node_MRCA == 7);
    }

    SECTION("tree_generator_next_tree_one_tree_case")
    {
        coa_table_c table_event;
        //one tree case
        table_event.set_coa_event(2, 10, 4, {2, 3}, 0.071, -1);
        table_event.set_coa_event(0, 2, 5, {1, 3}, 0.090, -1);
        table_event.set_coa_event(2, 10, 5, {1, 4}, 0.090, -1);
        table_event.set_coa_event(0, 7, 6, {0, 5}, 0.170, -1);
        table_event.set_coa_event(7, 10, 7, {0, 5}, 0.202, -1);
        table_event.set_coa_event(0, 2, 8, {2, 6}, 0.253, -1);

        tree_gen_c generator(table_event, 9, 4);

        auto result = generator.set_next_tree();

        auto tree = getAncestryTree(generator);

        REQUIRE(std::get<0>(result) == 0);
        REQUIRE(std::get<1>(result) == 2);

        REQUIRE((std::get<0>(tree[0])).empty());
        REQUIRE(std::get<1>(tree[0]) == 0);

        REQUIRE((std::get<0>(tree[5]))[0] == 1);
        REQUIRE((std::get<0>(tree[5]))[1] == 3);
        REQUIRE(std::get<1>(tree[5]) == 0.090);

        REQUIRE((std::get<0>(tree[8]))[0] == 2);
        REQUIRE((std::get<0>(tree[8]))[1] == 6);
        REQUIRE(std::get<1>(tree[8]) == 0.253);
    }

    SECTION("tree_generator_next_two_three_tree_case")
    {
        coa_table_c table_event;

        table_event.set_coa_event(2, 10, 4, {2, 3}, 0.071, -1);
        table_event.set_coa_event(0, 2, 5, {1, 3}, 0.090, -1);
        table_event.set_coa_event(2, 10, 5, {1, 4}, 0.090, -1);
        table_event.set_coa_event(0, 7, 6, {0, 5}, 0.170, -1);
        table_event.set_coa_event(7, 10, 7, {0, 5}, 0.202, -1);
        table_event.set_coa_event(0, 2, 8, {2, 6}, 0.253, -1);

        tree_gen_c generator(table_event, 9, 4);

        generator.set_next_tree();
        //View second tree
        auto result2 = generator.set_next_tree();

        auto tree2 = getAncestryTree(generator);

        REQUIRE(generator.MRCA_node_nbr == 6);

        REQUIRE(std::get<0>(result2) == 2);
        REQUIRE(std::get<1>(result2) == 7);

        REQUIRE((std::get<0>(tree2[3])).empty());
        REQUIRE(std::get<1>(tree2[3]) == 0);

        REQUIRE((std::get<0>(tree2[5]))[0] == 1);
        REQUIRE((std::get<0>(tree2[5]))[1] == 4);
        REQUIRE(std::get<1>(tree2[5]) == 0.090);

        REQUIRE((std::get<0>(tree2[8])).empty());
        REQUIRE(std::get<1>(tree2[8]) == 0);
    }

    SECTION("tree_generator_next_tree_three_tree_case")
    {
        coa_table_c table_event;

        table_event.set_coa_event(2, 10, 4, {2, 3}, 0.071, -1);
        table_event.set_coa_event(0, 2, 5, {1, 3}, 0.090, -1);
        table_event.set_coa_event(2, 10, 5, {1, 4}, 0.090, -1);
        table_event.set_coa_event(0, 7, 6, {0, 5}, 0.170, -1);
        table_event.set_coa_event(7, 10, 7, {0, 5}, 0.202, -1);
        table_event.set_coa_event(0, 2, 8, {2, 6}, 0.253, -1);

        tree_gen_c generator(table_event, 9, 4);

        generator.set_next_tree();
        generator.set_next_tree();
        //View second tree
        auto result2 = generator.set_next_tree();

        auto tree2 = getAncestryTree(generator);

        REQUIRE(generator.MRCA_node_nbr == 7);

        REQUIRE(std::get<0>(result2) == 7);
        REQUIRE(std::get<1>(result2) == 10);

        REQUIRE((std::get<0>(tree2[3])).empty());
        REQUIRE(std::get<1>(tree2[3]) == 0);

        REQUIRE((std::get<0>(tree2[6])).empty());
        REQUIRE(std::get<1>(tree2[6]) == 0);

        REQUIRE((std::get<0>(tree2[7]))[0] == 0);
        REQUIRE((std::get<0>(tree2[7]))[1] == 5);
        REQUIRE(std::get<1>(tree2[7]) == 0.202);
    }

    //Kellerher exemple with sample 4 at begin and 5 at end
    SECTION("tree_generator_multi_coa_case")
    {
        coa_table_c table_event;

        table_event.set_coa_event(2, 10, 6, {3, 4, 5}, 0.071, -1);
        table_event.set_coa_event(0, 2, 7, {1, 3, 5}, 0.090, -1);
        table_event.set_coa_event(2, 10, 7, {1, 2, 6}, 0.090, -1);
        table_event.set_coa_event(0, 7, 8, {0, 2, 7}, 0.170, -1);
        table_event.set_coa_event(7, 10, 9, {0, 7}, 0.202, -1);
        table_event.set_coa_event(0, 2, 10, {4, 8}, 0.253, -1);

        tree_gen_c generator(table_event, 11, 6);

        generator.set_next_tree();
        generator.set_next_tree();
        //View second tree
        auto result2 = generator.set_next_tree();

        auto tree2 = getAncestryTree(generator);

        REQUIRE(generator.MRCA_node_nbr == 9);

        REQUIRE(std::get<0>(result2) == 7);
        REQUIRE(std::get<1>(result2) == 10);

        REQUIRE((std::get<0>(tree2[3])).empty());
        REQUIRE(std::get<1>(tree2[3]) == 0);

        REQUIRE((std::get<0>(tree2[6]))[0] == 3);
        REQUIRE((std::get<0>(tree2[6]))[1] == 4);
        REQUIRE((std::get<0>(tree2[6]))[2] == 5);
        REQUIRE(std::get<1>(tree2[6]) == 0.071);

        REQUIRE((std::get<0>(tree2[7]))[0] == 1);
        REQUIRE((std::get<0>(tree2[7]))[1] == 2);
        REQUIRE((std::get<0>(tree2[7]))[2] == 6);
        REQUIRE(std::get<1>(tree2[7]) == 0.090);

        REQUIRE((std::get<0>(tree2[8])).empty());
        REQUIRE(std::get<1>(tree2[8]) == 0);

        REQUIRE((std::get<0>(tree2[9]))[0] == 0);
        REQUIRE((std::get<0>(tree2[9]))[1] == 7);
        REQUIRE(std::get<1>(tree2[9]) == 0.202);
    }
}
