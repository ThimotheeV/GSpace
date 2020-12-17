
#define CATCH_CONFIG_MAIN FenwickTreeTest
#include "catch.hpp"
#include <iostream>

#include "fenwick_tree.hpp"

namespace unit_test
{
struct FenwickTreeTest
{
    // FenwickTreeTest() {}
    // ~FenwickTreeTest() {}

    static void increase_tree_size(fenwick_tree_c &f, std::size_t index)
    {
        f.increase_tree_size(index);
    }

    static std::vector<long int> &getFenwickTree(fenwick_tree_c &f)
    {
        return f.Fenwick_tree;
    }

    static std::size_t getMidldle_index_value(fenwick_tree_c &f)
    {
        return f.Middle_index_value;
    }
};
} // namespace unit_test

TEST_CASE_METHOD(unit_test::FenwickTreeTest, "fenwick_tree_test")
{
    SECTION("basic_constructor")
    {
        fenwick_tree_c tree_one;

        REQUIRE(!(getFenwickTree(tree_one)[0]));
        REQUIRE(getMidldle_index_value(tree_one) == 1);
    }

    SECTION("advance_constructor")
    {
        fenwick_tree_c tree_twelve;
        fenwick_tree_c tree_forty_four;
        tree_twelve.create_fenwick_tree(12);
        tree_forty_four.create_fenwick_tree(44);

        REQUIRE(!(getFenwickTree(tree_twelve)[0]));
        REQUIRE(getMidldle_index_value(tree_twelve) == 8);
        REQUIRE(!(getFenwickTree(tree_forty_four)[0]));
        REQUIRE(getMidldle_index_value(tree_forty_four) == 32);
    }

    SECTION("increase_size")
    {
        fenwick_tree_c tree;
        tree.create_fenwick_tree(7);

        increase_tree_size(tree, 12);

        REQUIRE(!(getFenwickTree(tree)[0]));
        REQUIRE(getMidldle_index_value(tree) == 8);
        REQUIRE(getFenwickTree(tree).size() == 17);
    }

    SECTION("increment")
    {
        fenwick_tree_c tree;
        tree.create_fenwick_tree(7);
        tree.up_at_index_with_value(0, 5);

        REQUIRE(!(getFenwickTree(tree)[0]));
        REQUIRE(getFenwickTree(tree)[1] == 5);
        REQUIRE(getFenwickTree(tree)[2] == 5);
        REQUIRE(getFenwickTree(tree)[4] == 5);
        REQUIRE(getFenwickTree(tree)[3] == 0);
    }

    SECTION("increment_increase_size")
    {
        fenwick_tree_c tree;
        tree.create_fenwick_tree(7);
        tree.up_at_index_with_value(4, 1);
        tree.up_at_index_with_value(8, 5);

        REQUIRE(!(getFenwickTree(tree)[0]));
        REQUIRE(getFenwickTree(tree)[5] == 1);
        REQUIRE(getFenwickTree(tree)[8] == 1);
        REQUIRE(getFenwickTree(tree)[9] == 5);
        REQUIRE(getFenwickTree(tree)[10] == 5);
        REQUIRE(getFenwickTree(tree)[11] == 0);
    }

    SECTION("get_frequency")
    {
        fenwick_tree_c tree;
        tree.create_fenwick_tree(7);
        tree.up_at_index_with_value(0, 5);

        REQUIRE(!(getFenwickTree(tree)[0]));
        REQUIRE(tree.get_frequency_index(0) == 5);
        REQUIRE(tree.get_frequency_index(1) == 0);
    }

    SECTION("set_value")
    {
        fenwick_tree_c tree;
        tree.create_fenwick_tree(7);
        tree.set_value_at_index(0, 5);

        REQUIRE(!(getFenwickTree(tree)[0]));
        REQUIRE(getFenwickTree(tree)[1] == 5);
        REQUIRE(getFenwickTree(tree)[2] == 5);
        REQUIRE(getFenwickTree(tree)[4] == 5);
        REQUIRE(getFenwickTree(tree)[3] == 0);

        tree.set_value_at_index(0, 3);

        REQUIRE(!(getFenwickTree(tree)[0]));
        REQUIRE(getFenwickTree(tree)[4] == 3);
    }

    SECTION("set_value_increase_size")
    {
        fenwick_tree_c tree;
        tree.create_fenwick_tree(7);
        tree.set_value_at_index(8, 5);

        REQUIRE(!(getFenwickTree(tree)[0]));
        REQUIRE(getFenwickTree(tree)[8] == 0);
        REQUIRE(getFenwickTree(tree)[9] == 5);
        REQUIRE(getFenwickTree(tree)[10] == 5);
        REQUIRE(getFenwickTree(tree)[11] == 0);

        tree.set_value_at_index(11, 3);

        REQUIRE(!(getFenwickTree(tree)[0]));
        REQUIRE(getFenwickTree(tree)[12] == 8);
        REQUIRE(getFenwickTree(tree).size() == 17);
    }

    SECTION("get_cumulative_frequency")
    {
        fenwick_tree_c tree;
        tree.create_fenwick_tree(7);
        tree.set_value_at_index(0, 5);
        tree.set_value_at_index(3, 2);

        REQUIRE(!(getFenwickTree(tree)[0]));
        REQUIRE(tree.get_cumul_freq_at_index(2) == 5);
        REQUIRE(tree.get_cumul_freq_at_index(4) == 7);
    }

    SECTION("get_tot_cumul_freq")
    {
        fenwick_tree_c tree;
        tree.create_fenwick_tree(7);
        tree.set_value_at_index(3, 2);
        tree.set_value_at_index(6, 3);

        REQUIRE(!(getFenwickTree(tree)[0]));
        REQUIRE(tree.get_tot_cumul_freq() == 5);
    }

    SECTION("seg_index")
    {
        fenwick_tree_c tree;
        tree.create_fenwick_tree(7);

        tree.set_value_at_index(0, 2);
        tree.set_value_at_index(1, 3);
        tree.set_value_at_index(2, 4);
        tree.set_value_at_index(4, 5);

        REQUIRE(!(getFenwickTree(tree)[0]));
        REQUIRE(tree.seg_index(10) == 4);
    }

    SECTION("FWT massive amount of value")
    {
        
    }
}