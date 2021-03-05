#define CATCH_CONFIG_MAIN Debug
#include "catch.hpp"
#include <iostream>

#include "debug.hpp"

TEST_CASE("check_tree_test")
{
    std::vector<std::tuple<std::vector<int>, double, int>> gene_tree(7);
    gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0);
    gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0);
    gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0);
    gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0);
    gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1);
    gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1);
    gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, 2, -1);

    check_tree(gene_tree, 6, 4, 0);

    REQUIRE(1);
}
