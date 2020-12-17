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

    static std::vector<std::tuple<std::vector<int>, double, int, int>> &getAncestryTree(tree_gen_c &tg)
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
        table_event.set_coa_event(0, 1, 103, {102, 7}, 0, 1);
        table_event.set_coa_event(0, 1, 104, {103, 83}, 0, 1);
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

    SECTION("nbr_leaf_for_node")
    {
        coa_table_c table_event;

        tree_gen_c generator(table_event, 9, 4);

        getAncestryTree(generator)[5] = std::make_tuple(std::vector<int>(1, 3), 0, 0, 2);
        getAncestryTree(generator)[6] = std::make_tuple(std::vector<int>(0, 5), 0, 0, 3);
        getAncestryTree(generator)[8] = std::make_tuple(std::vector<int>(2, 6), 0, 0, 4);

        REQUIRE(generator.nbr_leaf_for_node({1, 3}) == 2);
        REQUIRE(generator.nbr_leaf_for_node({0, 5}) == 3);
        REQUIRE(generator.nbr_leaf_for_node({2, 6}) == 4);
    }

    SECTION("find_MRCA")
    {
        std::vector<std::tuple<std::vector<int>, double, int, int>> ancestry_tree(9);
        ancestry_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        ancestry_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        ancestry_tree[2] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        ancestry_tree[3] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        ancestry_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1, 2);
        ancestry_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1, 2);
        ancestry_tree[6] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        ancestry_tree[7] = std::make_tuple(std::vector<int>{4, 5}, 2, -1, 4);
        ancestry_tree[8] = std::make_tuple(std::vector<int>(), 0, 0, 1);

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
        REQUIRE(std::get<3>(tree[0]) == 1);

        REQUIRE((std::get<0>(tree[5]))[0] == 1);
        REQUIRE((std::get<0>(tree[5]))[1] == 3);
        REQUIRE(std::get<1>(tree[5]) == 0.090);
        REQUIRE(std::get<3>(tree[5]) == 2);

        REQUIRE((std::get<0>(tree[8]))[0] == 2);
        REQUIRE((std::get<0>(tree[8]))[1] == 6);
        REQUIRE(std::get<1>(tree[8]) == 0.253);
        REQUIRE(std::get<3>(tree[8]) == 4);
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
        REQUIRE(std::get<3>(tree2[3]) == 1);

        REQUIRE((std::get<0>(tree2[5]))[0] == 1);
        REQUIRE((std::get<0>(tree2[5]))[1] == 4);
        REQUIRE(std::get<1>(tree2[5]) == 0.090);
        REQUIRE(std::get<3>(tree2[5]) == 3);

        REQUIRE((std::get<0>(tree2[8])).empty());
        REQUIRE(std::get<1>(tree2[8]) == 0);
        REQUIRE(std::get<3>(tree2[8]) == 0);
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
        REQUIRE(std::get<3>(tree2[3]) == 1);

        REQUIRE((std::get<0>(tree2[6])).empty());
        REQUIRE(std::get<1>(tree2[6]) == 0);
        REQUIRE(std::get<3>(tree2[6]) == 0);

        REQUIRE((std::get<0>(tree2[7]))[0] == 0);
        REQUIRE((std::get<0>(tree2[7]))[1] == 5);
        REQUIRE(std::get<1>(tree2[7]) == 0.202);
        REQUIRE(std::get<3>(tree2[7]) == 4);
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
        REQUIRE(std::get<3>(tree2[3]) == 1);

        REQUIRE((std::get<0>(tree2[6]))[0] == 3);
        REQUIRE((std::get<0>(tree2[6]))[1] == 4);
        REQUIRE((std::get<0>(tree2[6]))[2] == 5);
        REQUIRE(std::get<1>(tree2[6]) == 0.071);
        REQUIRE(std::get<3>(tree2[6]) == 3);

        REQUIRE((std::get<0>(tree2[7]))[0] == 1);
        REQUIRE((std::get<0>(tree2[7]))[1] == 2);
        REQUIRE((std::get<0>(tree2[7]))[2] == 6);
        REQUIRE(std::get<1>(tree2[7]) == 0.090);
        REQUIRE(std::get<3>(tree2[7]) == 5);

        REQUIRE((std::get<0>(tree2[8])).empty());
        REQUIRE(std::get<1>(tree2[8]) == 0);
        REQUIRE(std::get<3>(tree2[8]) == 0);

        REQUIRE((std::get<0>(tree2[9]))[0] == 0);
        REQUIRE((std::get<0>(tree2[9]))[1] == 7);
        REQUIRE(std::get<1>(tree2[9]) == 0.202);
        REQUIRE(std::get<3>(tree2[9]) == 6);
    }
}
