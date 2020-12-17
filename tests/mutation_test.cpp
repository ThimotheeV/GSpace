#define CATCH_CONFIG_MAIN MutTest
#include "catch.hpp"
#include <iostream>

#include "mutation.hpp"
#include "settings.hpp"
#include "summary_stat.hpp"

TEST_CASE("nucl_distrib_test")
{
    SECTION("cumul_distrib")
    {
        auto result = cumul_distrib({0, 1, 1, 1});

        REQUIRE(result[0] == Approx(0).margin(1));
        REQUIRE(result[1] == Approx(333333).margin(1));
        REQUIRE(result[2] == Approx(666666).margin(1));
        REQUIRE(result[3] == Approx(1000000).margin(1));
    }
}
TEST_CASE("mutation_test")
{
    //2 indiv_c, 3 nuceotides => indiv_c 0 : 2 muts ; indiv_c 1 : 1 muts ; ancestre 2 (no mut)
    //indiv_c 0 mut 2
    SECTION("apply_mut_to_site_mut")
    {
        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        auto &simu_param = singleton_c<simu_param_c>::instance();
        simu_param.Continuous_time_approxim = true;
        rand_gen.put_seed(3);

        mut_genenerator_c mut_generator;

        std::map<int, int> sample_mutated_state;
        sample_mutated_state.emplace(2, 1);
        auto mut_jcm = sample_mutated_state.find(2);

        mut_generator.Mut_model.apply_mut_to_site(mut_model_enum::jcm, mut_jcm);

        REQUIRE(mut_jcm->second == 3);
        REQUIRE(sizeof(mut_jcm->second) == 4);
    }

    //indiv_c 0
    SECTION("mut_events")
    {
        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        auto &simu_param = singleton_c<simu_param_c>::instance();
        simu_param.Continuous_time_approxim = true;
        rand_gen.put_seed(3);

        auto &muta_param = singleton_c<muta_param_c>::instance();
        muta_param.Mod_mut_name = mut_model_enum::jcm;
        muta_param.Scaled_mut_rate_theta = 1;

        mut_genenerator_c mut_generator;

        std::vector<int> ancestry_seq(3, 2);

        mut_generator.All_indivs_mut_state = std::vector<std::map<int, int>>(3);

        int node_index = 0;

        mut_generator.Begin_end_sequence = {0, 3};

        double long_branch = 1;

        mut_generator.mut_events(true, node_index, long_branch, muta_param.Scaled_mut_rate_theta, ancestry_seq, rand_gen);

        REQUIRE(mut_generator.All_indivs_mut_state[0].find(0) == mut_generator.All_indivs_mut_state[0].end());

        REQUIRE(mut_generator.All_indivs_mut_state[0].find(1) != mut_generator.All_indivs_mut_state[0].end());
        REQUIRE(mut_generator.All_indivs_mut_state[0].find(1)->second == 3);

        REQUIRE(mut_generator.All_indivs_mut_state[0].find(2) != mut_generator.All_indivs_mut_state[0].end());
        REQUIRE(mut_generator.All_indivs_mut_state[0].find(2)->second == 1);
    }

    SECTION("browse_tree")
    {
        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        auto &simu_param = singleton_c<simu_param_c>::instance();
        simu_param.Continuous_time_approxim = true;
        rand_gen.put_seed(3);

        auto &muta_param = singleton_c<muta_param_c>::instance();
        muta_param.Mod_mut_name = mut_model_enum::jcm;
        muta_param.Scaled_mut_rate_theta = 1;

        mut_genenerator_c mut_generator;

        auto gene_tree = std::vector<std::tuple<std::vector<int>, double, int, int>>(3);
        gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[2] = std::make_tuple(std::vector<int>{0, 1}, 1, -1, 2);

        std::vector<int> ancestry_seq(3, 2);

        mut_generator.All_indivs_mut_state = std::vector<std::map<int, int>>(3);
        mut_generator.Begin_end_sequence = {0, 3};

        std::vector<int> node_file;
        node_file.reserve(3);
        node_file.push_back(2);

        mut_generator.browse_tree(true, gene_tree, node_file, ancestry_seq);

        REQUIRE(mut_generator.All_indivs_mut_state[0].find(0) == mut_generator.All_indivs_mut_state[0].end());
        REQUIRE(mut_generator.All_indivs_mut_state[0].find(1) != mut_generator.All_indivs_mut_state[0].end());
        REQUIRE(mut_generator.All_indivs_mut_state[0].find(1)->second == 3);
        REQUIRE(mut_generator.All_indivs_mut_state[0].find(2) != mut_generator.All_indivs_mut_state[0].end());
        REQUIRE(mut_generator.All_indivs_mut_state[0].find(2)->second == 1);

        REQUIRE(mut_generator.All_indivs_mut_state[1].find(0) != mut_generator.All_indivs_mut_state[1].end());
        REQUIRE(mut_generator.All_indivs_mut_state[1].find(0)->second == 1);
        REQUIRE(mut_generator.All_indivs_mut_state[1].find(1) == mut_generator.All_indivs_mut_state[1].end());
        REQUIRE(mut_generator.All_indivs_mut_state[1].find(2) != mut_generator.All_indivs_mut_state[1].end());
        REQUIRE(mut_generator.All_indivs_mut_state[1].find(2)->second == 1);
    }

    SECTION("mut_generation")
    {
        auto &muta_param = singleton_c<muta_param_c>::instance();
        muta_param.Mod_mut_name = mut_model_enum::jcm;
        muta_param.Scaled_mut_rate_theta = 1;

        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.n_total_sample_size = 4;

        std::vector<std::tuple<std::vector<int>, double, int, int>> gene_tree(7);
        gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1, 2);
        gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1, 2);
        gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, 2, -1, 4);

        int nbr_of_recomb_site_m = 3;
        std::vector<int> ancestry_seq(nbr_of_recomb_site_m, 2);

        std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};
        int MRCA_index = 6;

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        auto &simu_param = singleton_c<simu_param_c>::instance();
        simu_param.Continuous_time_approxim = true;
        rand_gen.put_seed(3);

        mut_genenerator_c mut_generator;

        auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, MRCA_index);

        REQUIRE(mut_generator.All_indivs_mut_state[0].at(0) == 1);
        REQUIRE(mut_generator.All_indivs_mut_state[0].at(1) == 3);
        REQUIRE(mut_generator.All_indivs_mut_state[0].at(2) == 2);
        REQUIRE(mut_generator.All_indivs_mut_state[0].size() == 3);

        REQUIRE(mut_generator.All_indivs_mut_state[1].at(1) == 3);
        REQUIRE(mut_generator.All_indivs_mut_state[1].at(2) == 1);
        REQUIRE(mut_generator.All_indivs_mut_state[1].size() == 2);

        REQUIRE(mut_generator.All_indivs_mut_state[2].at(0) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[2].at(1) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[2].at(2) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[2].size() == 3);

        REQUIRE(mut_generator.All_indivs_mut_state[3].at(0) == 3);
        REQUIRE(mut_generator.All_indivs_mut_state[3].at(1) == 2);
        REQUIRE(mut_generator.All_indivs_mut_state[3].at(2) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[3].size() == 3);
    }

    SECTION("mut_generation_gen_by_gen")
    {
        auto &simu_param = singleton_c<simu_param_c>::instance();
        simu_param.Continuous_time_approxim = false;

        auto &muta_param = singleton_c<muta_param_c>::instance();
        muta_param.Mod_mut_name = mut_model_enum::jcm;
        muta_param.Unscaled_mut_rate_mu = 1;

        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.n_total_sample_size = 4;

        std::vector<std::tuple<std::vector<int>, double, int, int>> gene_tree(7);
        gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, -10, 1, 2);
        gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, -10, 1, 2);
        gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, -10, 2, 4);

        int nbr_of_recomb_site_m = 3;
        std::vector<int> ancestry_seq(nbr_of_recomb_site_m, 2);

        std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};
        int MRCA_index = 6;

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        mut_genenerator_c mut_generator;

        auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, MRCA_index);

        REQUIRE(mut_generator.All_indivs_mut_state[0].at(0) == 1);
        REQUIRE(mut_generator.All_indivs_mut_state[0].at(1) == 2);
        REQUIRE(mut_generator.All_indivs_mut_state[0].at(2) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[0].size() == 3);

        REQUIRE(mut_generator.All_indivs_mut_state[1].at(0) == 1);
        REQUIRE(mut_generator.All_indivs_mut_state[1].at(1) == 2);
        REQUIRE(mut_generator.All_indivs_mut_state[1].at(2) == 1);
        REQUIRE(mut_generator.All_indivs_mut_state[1].size() == 3);

        REQUIRE(mut_generator.All_indivs_mut_state[2].at(0) == 3);
        REQUIRE(mut_generator.All_indivs_mut_state[2].at(1) == 3);
        REQUIRE(mut_generator.All_indivs_mut_state[2].at(2) == 1);
        REQUIRE(mut_generator.All_indivs_mut_state[2].size() == 3);

        REQUIRE(mut_generator.All_indivs_mut_state[3].at(0) == 3);
        REQUIRE(mut_generator.All_indivs_mut_state[3].at(1) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[3].at(2) == 2);
        REQUIRE(mut_generator.All_indivs_mut_state[3].size() == 3);
    }

    SECTION("multi_descendant_init")
    {
        auto &simu_param = singleton_c<simu_param_c>::instance();
        simu_param.Continuous_time_approxim = true;

        auto &muta_param = singleton_c<muta_param_c>::instance();
        muta_param.Mod_mut_name = mut_model_enum::jcm;
        muta_param.Scaled_mut_rate_theta = 1;

        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.n_total_sample_size = 4;

        std::vector<std::tuple<std::vector<int>, double, int, int>> gene_tree(7);
        gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[4] = std::make_tuple(std::vector<int>{0, 1, 2}, 1, -1, 3);
        gene_tree[5] = std::make_tuple(std::vector<int>{3, 4}, 2, -1, 2);

        int nbr_of_recomb_site_m = 3;
        std::vector<int> ancestry_seq(nbr_of_recomb_site_m, 2);

        std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};
        int MRCA_index = 5;

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        mut_genenerator_c mut_generator;

        auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, MRCA_index);

        REQUIRE(mut_generator.All_indivs_mut_state[0].at(1) == 2);
        REQUIRE(mut_generator.All_indivs_mut_state[0].size() == 1);

        REQUIRE(mut_generator.All_indivs_mut_state[2].at(0) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[2].at(1) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[2].at(2) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[2].size() == 3);

        REQUIRE(mut_generator.All_indivs_mut_state[3].at(0) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[3].at(1) == 1);
        REQUIRE(mut_generator.All_indivs_mut_state[3].at(2) == 3);
        REQUIRE(mut_generator.All_indivs_mut_state[3].size() == 3);
    }
    /******************/

    SECTION("mut_with_mod_iam")
    {
        std::vector<std::tuple<std::vector<int>, double, int, int>> gene_tree(7);
        gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1, 2);
        gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1, 2);
        gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, 2, -1, 4);

        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sequence_length = 3;
        std::vector<int> ancestry_seq(samp_param.Sequence_length, 1);

        std::array<int, 2> begin_end_seq_current_tree = {0, samp_param.Sequence_length};

        auto &muta_param = singleton_c<muta_param_c>::instance();
        muta_param.Mod_mut_name = mut_model_enum::iam;
        muta_param.Scaled_mut_rate_theta = 0.5;

        samp_param.n_total_sample_size = 4;

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        auto &simu_param = singleton_c<simu_param_c>::instance();
        simu_param.Continuous_time_approxim = true;
        rand_gen.put_seed(3);

        int MRCA_index = 6;

        mut_genenerator_c mut_generator;

        auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, MRCA_index);

        REQUIRE(mut_generator.All_indivs_mut_state[0].at(0) == 3);
        REQUIRE(mut_generator.All_indivs_mut_state[0].at(1) == 2);
        REQUIRE(mut_generator.All_indivs_mut_state[0].size() == 2);

        REQUIRE(mut_generator.All_indivs_mut_state[1].at(0) == 5);
        REQUIRE(mut_generator.All_indivs_mut_state[1].at(1) == 2);
        REQUIRE(mut_generator.All_indivs_mut_state[1].at(2) == 2);
        REQUIRE(mut_generator.All_indivs_mut_state[1].size() == 3);

        REQUIRE(mut_generator.All_indivs_mut_state[2].at(0) == 7);
        REQUIRE(mut_generator.All_indivs_mut_state[2].at(1) == 3);
        REQUIRE(mut_generator.All_indivs_mut_state[2].size() == 2);

        REQUIRE(mut_generator.All_indivs_mut_state[3].at(1) == 3);
        REQUIRE(mut_generator.All_indivs_mut_state[3].size() == 1);
    }

    SECTION("mut_mod_c_kam")
    {
        std::vector<std::tuple<std::vector<int>, double, int, int>> gene_tree(7);
        gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1, 2);
        gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1, 2);
        gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, 2, -1, 4);

        int nbr_of_recomb_site_m = 3;
        std::vector<int> ancestry_seq(3, 5);

        auto &muta_param = singleton_c<muta_param_c>::instance();
        muta_param.Mod_mut_name = mut_model_enum::kam;
        muta_param.Scaled_mut_rate_theta = 1;
        muta_param.K_min = 1;
        muta_param.K_max = 10;

        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.n_total_sample_size = 4;

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        auto &simu_param = singleton_c<simu_param_c>::instance();
        simu_param.Continuous_time_approxim = true;
        rand_gen.put_seed(3);

        int MRCA_index = 6;

        mut_genenerator_c mut_generator;

        std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};

        auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, MRCA_index);

        REQUIRE(mut_generator.All_indivs_mut_state[0].at(0) == 7);
        REQUIRE(mut_generator.All_indivs_mut_state[0].at(1) == 6);
        REQUIRE(mut_generator.All_indivs_mut_state[0].at(2) == 8);
        REQUIRE(mut_generator.All_indivs_mut_state[0].size() == 3);

        REQUIRE(mut_generator.All_indivs_mut_state[1].at(0) == 9);
        REQUIRE(mut_generator.All_indivs_mut_state[1].at(1) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[1].at(2) == 8);
        REQUIRE(mut_generator.All_indivs_mut_state[1].size() == 3);

        REQUIRE(mut_generator.All_indivs_mut_state[2].at(0) == 8);
        REQUIRE(mut_generator.All_indivs_mut_state[2].at(1) == 10);
        REQUIRE(mut_generator.All_indivs_mut_state[2].at(2) == 8);
        REQUIRE(mut_generator.All_indivs_mut_state[2].size() == 3);

        REQUIRE(mut_generator.All_indivs_mut_state[3].at(0) == 2);
        REQUIRE(mut_generator.All_indivs_mut_state[3].at(1) == 6);
        REQUIRE(mut_generator.All_indivs_mut_state[3].at(2) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[3].size() == 3);
    }

    SECTION("mut_mod_c_smm")
    {
        std::vector<std::tuple<std::vector<int>, double, int, int>> gene_tree(7);
        gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1, 2);
        gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1, 2);
        gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, 2, -1, 4);

        int nbr_of_recomb_site_m = 3;
        std::vector<int> ancestry_seq(3, 5);

        std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};

        auto &muta_param = singleton_c<muta_param_c>::instance();
        muta_param.Mod_mut_name = mut_model_enum::smm;
        muta_param.Scaled_mut_rate_theta = 1;
        muta_param.K_min = 4;
        muta_param.K_max = 6;

        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.n_total_sample_size = 4;

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        auto &simu_param = singleton_c<simu_param_c>::instance();
        simu_param.Continuous_time_approxim = true;
        rand_gen.put_seed(3);

        int MRCA_index = 6;

        mut_genenerator_c mut_generator;

        auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, MRCA_index);

        REQUIRE(mut_generator.All_indivs_mut_state[0].at(0) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[0].at(1) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[0].at(2) == 6);
        REQUIRE(mut_generator.All_indivs_mut_state[0].size() == 3);

        REQUIRE(mut_generator.All_indivs_mut_state[1].at(1) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[1].at(2) == 5);
        REQUIRE(mut_generator.All_indivs_mut_state[1].size() == 2);

        REQUIRE(mut_generator.All_indivs_mut_state[2].at(0) == 5);
        REQUIRE(mut_generator.All_indivs_mut_state[2].at(1) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[2].at(2) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[2].size() == 3);

        REQUIRE(mut_generator.All_indivs_mut_state[3].at(0) == 6);
        REQUIRE(mut_generator.All_indivs_mut_state[3].at(1) == 5);
        REQUIRE(mut_generator.All_indivs_mut_state[3].at(2) == 5);
        REQUIRE(mut_generator.All_indivs_mut_state[3].size() == 3);
    }
    /******************/
    SECTION("mut_mod_c_gsm")
    {
        std::vector<std::tuple<std::vector<int>, double, int, int>> gene_tree(7);
        gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0, 1);
        gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1, 2);
        gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1, 2);
        gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, 2, -1, 4);

        int nbr_of_recomb_site_m = 3;
        std::vector<int> ancestry_seq(3, 5);

        std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};

        auto &muta_param = singleton_c<muta_param_c>::instance();
        muta_param.Mod_mut_name = mut_model_enum::gsm;
        muta_param.Scaled_mut_rate_theta = 1;
        muta_param.K_min = 1;
        muta_param.K_max = 10;
        muta_param.P_gsm = 0.5;

        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.n_total_sample_size = 4;

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        auto &simu_param = singleton_c<simu_param_c>::instance();
        simu_param.Continuous_time_approxim = true;
        rand_gen.put_seed(3);

        int MRCA_index = 6;

        mut_genenerator_c mut_generator;

        auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, MRCA_index);

        REQUIRE(mut_generator.All_indivs_mut_state[0].at(0) == 3);
        REQUIRE(mut_generator.All_indivs_mut_state[0].at(1) == 8);
        REQUIRE(mut_generator.All_indivs_mut_state[0].size() == 2);

        REQUIRE(mut_generator.All_indivs_mut_state[1].at(0) == 3);
        REQUIRE(mut_generator.All_indivs_mut_state[1].at(1) == 8);
        REQUIRE(mut_generator.All_indivs_mut_state[1].size() == 2);

        REQUIRE(mut_generator.All_indivs_mut_state[2].at(0) == 8);
        REQUIRE(mut_generator.All_indivs_mut_state[2].at(1) == 7);
        REQUIRE(mut_generator.All_indivs_mut_state[2].at(2) == 9);
        REQUIRE(mut_generator.All_indivs_mut_state[2].size() == 3);

        REQUIRE(mut_generator.All_indivs_mut_state[3].at(0) == 8);
        REQUIRE(mut_generator.All_indivs_mut_state[3].at(1) == 9);
        REQUIRE(mut_generator.All_indivs_mut_state[3].at(2) == 4);
        REQUIRE(mut_generator.All_indivs_mut_state[3].size() == 3);
    }

    //     SECTION("mut_mod_c_k80")
    //     {
    //         std::vector<std::tuple<std::vector<int>, double, int, int>> gene_tree(7);
    //         gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1, 2);
    //         gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1, 2);
    //         gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, 2, -1, 4);

    //         int nbr_of_recomb_site_m = 3;
    //         std::vector<int> ancestry_seq(nbr_of_recomb_site_m, 0);

    //         std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};

    //         auto &muta_param = singleton_c<muta_param_c>::instance();
    //         muta_param.Mod_mut_name = mut_model_enum::k80;
    //         muta_param.Scaled_mut_rate_theta = 0.5;
    //         muta_param.Ratio_transi_transver = {0.5, 0.7};

    //         auto &samp_param = singleton_c<samp_param_c>::instance();
    //         samp_param.n_total_sample_size = 4;

    //         auto &rand_gen = singleton_c<rand_gen_c>::instance();
    // auto &simu_param = singleton_c<simu_param_c>::instance();
    // simu_param.Continuous_time_approxim = true;
    //         rand_gen.put_seed(3);

    //         int MRCA_index = 6;

    //         mut_genenerator_c mut_generator;

    //         auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, MRCA_index);

    //         REQUIRE(mut_generator.All_indivs_mut_state[0].at(0) == 2);
    //         REQUIRE(mut_generator.All_indivs_mut_state[0].at(1) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[0].size() == 2);

    //         REQUIRE(mut_generator.All_indivs_mut_state[1].at(0) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[1].at(1) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[1].at(2) == 2);
    //         REQUIRE(mut_generator.All_indivs_mut_state[1].size() == 3);

    //         REQUIRE(mut_generator.All_indivs_mut_state[2].at(0) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[2].at(2) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[2].size() == 2);

    //         REQUIRE(mut_generator.All_indivs_mut_state[3].at(2) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[3].size() == 1);
    //     }

    //     SECTION("mut_mod_c_f81")
    //     {
    //         std::vector<std::tuple<std::vector<int>, double, int, int>> gene_tree(7);
    //         gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1, 2);
    //         gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1, 2);
    //         gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, 2, -1, 4);

    //         int nbr_of_recomb_site_m = 3;
    //         std::vector<int> ancestry_seq(nbr_of_recomb_site_m, 0);

    //         std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};

    //         auto &muta_param = singleton_c<muta_param_c>::instance();
    //         muta_param.Mod_mut_name = mut_model_enum::f81;
    //         muta_param.Scaled_mut_rate_theta = 0.5;
    //         muta_param.Equi_base_freq.A = 0.1;
    //         muta_param.Equi_base_freq.T = 0.4;
    //         muta_param.Equi_base_freq.C = 0.2;
    //         muta_param.Equi_base_freq.G = 0.3;

    //         auto &samp_param = singleton_c<samp_param_c>::instance();
    //         samp_param.n_total_sample_size = 4;

    //         auto &rand_gen = singleton_c<rand_gen_c>::instance();
    // auto &simu_param = singleton_c<simu_param_c>::instance();
    // simu_param.Continuous_time_approxim = true;
    //         rand_gen.put_seed(3);

    //         int MRCA_index = 6;

    //         mut_genenerator_c mut_generator;

    //         auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, MRCA_index);

    //         REQUIRE(mut_generator.All_indivs_mut_state[0].at(0) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[0].at(1) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[0].size() == 2);

    //         REQUIRE(mut_generator.All_indivs_mut_state[1].at(0) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[1].at(1) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[1].at(2) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[1].size() == 3);

    //         REQUIRE(mut_generator.All_indivs_mut_state[2].at(0) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[2].at(2) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[2].size() == 2);

    //         REQUIRE(mut_generator.All_indivs_mut_state[3].at(2) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[3].size() == 1);
    //     }

    //     SECTION("mut_mod_c_hky")
    //     {
    //         std::vector<std::tuple<std::vector<int>, double, int, int>> gene_tree(7);
    //         gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1, 2);
    //         gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1, 2);
    //         gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, 2, -1, 4);

    //         int nbr_of_recomb_site_m = 3;
    //         std::vector<int> ancestry_seq(nbr_of_recomb_site_m, 0);

    //         std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};

    //         auto &muta_param = singleton_c<muta_param_c>::instance();
    //         muta_param.Mod_mut_name = mut_model_enum::hky;
    //         muta_param.Scaled_mut_rate_theta = 0.5;
    //         muta_param.Ratio_transi_transver = {0.5, 0.7};
    //         muta_param.Equi_base_freq.A = 0.6;
    //         muta_param.Equi_base_freq.T = 0.1;
    //         muta_param.Equi_base_freq.C = 0.2;
    //         muta_param.Equi_base_freq.G = 0.1;

    //         auto &samp_param = singleton_c<samp_param_c>::instance();
    //         samp_param.n_total_sample_size = 4;

    //         auto &rand_gen = singleton_c<rand_gen_c>::instance();
    // auto &simu_param = singleton_c<simu_param_c>::instance();
    // simu_param.Continuous_time_approxim = true;
    //         rand_gen.put_seed(3);

    //         int MRCA_index = 6;

    //         mut_genenerator_c mut_generator;

    //         auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, MRCA_index);

    //         REQUIRE(mut_generator.All_indivs_mut_state[0].at(0) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[0].at(1) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[0].size() == 2);

    //         REQUIRE(mut_generator.All_indivs_mut_state[1].at(0) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[1].at(1) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[1].at(2) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[1].size() == 3);

    //         REQUIRE(mut_generator.All_indivs_mut_state[2].at(0) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[2].at(1) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[2].size() == 2);

    //         REQUIRE(mut_generator.All_indivs_mut_state[3].at(0) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[3].size() == 1);
    //     }

    //     SECTION("mut_mod_c_tn93")
    //     {
    //         std::vector<std::tuple<std::vector<int>, double, int, int>> gene_tree(7);
    //         gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1, 2);
    //         gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1, 2);
    //         gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, 2, -1, 4);

    //         int nbr_of_recomb_site_m = 3;
    //         std::vector<int> ancestry_seq(nbr_of_recomb_site_m, 0);

    //         std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};

    //         auto &muta_param = singleton_c<muta_param_c>::instance();
    //         muta_param.Mod_mut_name = mut_model_enum::tn93;
    //         muta_param.Scaled_mut_rate_theta = 0.5;
    //         muta_param.Ratio_transi_transver = {1.5, 0.5};
    //         muta_param.Equi_base_freq.A = 0.6;
    //         muta_param.Equi_base_freq.T = 0.1;
    //         muta_param.Equi_base_freq.C = 0.2;
    //         muta_param.Equi_base_freq.G = 0.1;

    //         auto &samp_param = singleton_c<samp_param_c>::instance();
    //         samp_param.n_total_sample_size = 4;

    //         auto &rand_gen = singleton_c<rand_gen_c>::instance();
    // auto &simu_param = singleton_c<simu_param_c>::instance();
    // simu_param.Continuous_time_approxim = true;
    //         rand_gen.put_seed(3);

    //         int MRCA_index = 6;

    //         mut_genenerator_c mut_generator;

    //         auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, MRCA_index);

    //         REQUIRE(mut_generator.All_indivs_mut_state[0].at(0) == 2);
    //         REQUIRE(mut_generator.All_indivs_mut_state[0].at(1) == 2);
    //         REQUIRE(mut_generator.All_indivs_mut_state[0].size() == 2);

    //         REQUIRE(mut_generator.All_indivs_mut_state[1].at(0) == 2);
    //         REQUIRE(mut_generator.All_indivs_mut_state[1].at(1) == 2);
    //         REQUIRE(mut_generator.All_indivs_mut_state[1].at(2) == 3);
    //         REQUIRE(mut_generator.All_indivs_mut_state[1].size() == 3);

    //         REQUIRE(mut_generator.All_indivs_mut_state[2].at(0) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[2].at(1) == 2);
    //         REQUIRE(mut_generator.All_indivs_mut_state[2].size() == 2);

    //         REQUIRE(mut_generator.All_indivs_mut_state[3].at(0) == 2);
    //         REQUIRE(mut_generator.All_indivs_mut_state[3].size() == 1);
    //     }

    //     /******************/

    //     SECTION("mut_with_mod_ism")
    //     {
    //         std::vector<std::tuple<std::vector<int>, double, int, int>> gene_tree(7);
    //         gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //         gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1, 2);
    //         gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1, 2);
    //         gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, 2, -1, 4);

    //         int nbr_of_recomb_site_m = 3;
    //         std::vector<int> ancestry_seq(nbr_of_recomb_site_m, 0);

    //         std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};

    //         auto &muta_param = singleton_c<muta_param_c>::instance();
    //         muta_param.Mod_mut_name = mut_model_enum::ism;
    //         muta_param.Scaled_mut_rate_theta = 0.5;

    //         auto &samp_param = singleton_c<samp_param_c>::instance();
    //         samp_param.n_total_sample_size = 4;
    //         samp_param.Sequence_length = nbr_of_recomb_site_m;

    //         auto &rand_gen = singleton_c<rand_gen_c>::instance();
    // auto &simu_param = singleton_c<simu_param_c>::instance();
    // simu_param.Continuous_time_approxim = true;
    //         rand_gen.put_seed(3);

    //         int MRCA_index = 6;

    //         mut_genenerator_c mut_generator;

    //         auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, MRCA_index);

    //         REQUIRE(mut_generator.All_indivs_mut_state[0].at(0) == 3);
    //         REQUIRE(mut_generator.All_indivs_mut_state[0].at(1) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[0].size() == 2);

    //         REQUIRE(mut_generator.All_indivs_mut_state[1].at(0) == 12);
    //         REQUIRE(mut_generator.All_indivs_mut_state[1].at(1) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[1].at(2) == 1);
    //         REQUIRE(mut_generator.All_indivs_mut_state[1].size() == 3);

    //         REQUIRE(mut_generator.All_indivs_mut_state[2].at(0) == 48);
    //         REQUIRE(mut_generator.All_indivs_mut_state[2].at(1) == 2);
    //         REQUIRE(mut_generator.All_indivs_mut_state[2].size() == 2);

    //         REQUIRE(mut_generator.All_indivs_mut_state[3].at(1) == 2);
    //         REQUIRE(mut_generator.All_indivs_mut_state[3].size() == 1);
    //     }
    // }

    // TEST_CASE("apply_mut_test")
    // {
    //     SECTION("haploid")
    //     {
    //         auto &simu_param = singleton_c<simu_param_c>::instance();
    //         simu_param.Continuous_time_approxim = false;

    //         //call in mut_process
    //         muta_param_c &muta_param = singleton_c<muta_param_c>::instance();
    //         muta_param.Mod_mut_name = mut_model_enum::iam;
    //         muta_param.Ancestry_seq = {1, 1};
    //         muta_param.Unscaled_mut_rate_mu = 1;

    //         auto &samp_param = singleton_c<samp_param_c>::instance();
    //         samp_param.Sequence_length = 2;
    //         samp_param.n_total_sample_size = 4;
    //         samp_param.Ploidy = 1;

    //         coa_table_c coa_table;
    //         coa_table.set_coa_event(0, 1, 4, {0, 1}, 0.090, 1);
    //         coa_table.set_coa_event(0, 1, 5, {2, 3}, 0.090, 1);
    //         coa_table.set_coa_event(1, 2, 4, {1, 3}, 0.090, 1);
    //         coa_table.set_coa_event(1, 2, 5, {0, 2}, 0.090, 1);
    //         coa_table.set_coa_event(0, 2, 6, {4, 5}, 0.090, 2);

    //         auto result = apply_mut(coa_table, 7, muta_param, samp_param);
    //         sample_mutated_state_type expect = {{{0, 3}, {1, 11}},
    //                                             {{0, 4}, {1, 9}},
    //                                             {{0, 5}, {1, 12}},
    //                                             {{0, 6}, {1, 10}}};
    //         REQUIRE(expect == result);
    //     }

    //     SECTION("diploid")
    //     {
    //         auto &simu_param = singleton_c<simu_param_c>::instance();
    //         simu_param.Continuous_time_approxim = false;

    //         //call in mut_process
    //         muta_param_c &muta_param = singleton_c<muta_param_c>::instance();
    //         muta_param.Mod_mut_name = mut_model_enum::iam;
    //         muta_param.Ancestry_seq = {1, 1};
    //         muta_param.Unscaled_mut_rate_mu = 1;

    //         auto &samp_param = singleton_c<samp_param_c>::instance();
    //         samp_param.Sequence_length = 2;
    //         samp_param.n_total_sample_size = 2;
    //         samp_param.Ploidy = 2;

    //         coa_table_c coa_table;
    //         coa_table.set_coa_event(0, 1, 4, {0, 1}, 0.090, 1);
    //         coa_table.set_coa_event(0, 1, 5, {2, 3}, 0.090, 1);
    //         coa_table.set_coa_event(1, 2, 4, {1, 3}, 0.090, 1);
    //         coa_table.set_coa_event(1, 2, 5, {0, 2}, 0.090, 1);
    //         coa_table.set_coa_event(0, 2, 6, {4, 5}, 0.090, 2);

    //         auto result = apply_mut(coa_table, 7, muta_param, samp_param);
    //         sample_mutated_state_type expect = {{{0, 3}, {1, 11}},
    //                                             {{0, 4}, {1, 9}},
    //                                             {{0, 5}, {1, 12}},
    //                                             {{0, 6}, {1, 10}}};
    //         REQUIRE(expect == result);
    //     }
}