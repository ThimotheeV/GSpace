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

        mut_jcm->second = mut_generator.Mut_model.apply_mut_to_site(mut_model_enum::jcm, mut_jcm->second, mut_jcm->first);

        REQUIRE(mut_jcm->second == 3);
        REQUIRE(sizeof(mut_jcm->second) == 4);
    }

    //situation : first branch with mutation in the tree
    SECTION("gen/gen add_mut")
    {
        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        auto &tskit = singleton_c<tskit_struct_c>::instance();
        tskit.Tskit_output = false;

        mut_genenerator_c mut_gen;

        int node_index = 0;
        int nbr_mut = 1000;
        int nbr_site = 10;
        std::vector<int> ancester_state = std::vector<int>(nbr_site, 1);

        mut_gen.Mod_mut_name = mut_model_enum::jcm;

        //Draw a site
        mut_gen.Site_distri = std::uniform_int_distribution<long int>(0, nbr_site - 1);
        //For each site, the markov chain mutation
        mut_gen.Mut_chain = std::map<int, std::list<int>>();
        //Recap all mut emplacement by branch
        mut_gen.Mut_register = std::list<std::list<int>>{std::list<int>{}};

        mut_gen.add_mut(ancester_state, rand_gen, node_index, nbr_mut);

        //Because nbr_of_mut >>> nbr_site
        REQUIRE(mut_gen.Mut_chain.size() == static_cast<std::size_t>(nbr_site));
        REQUIRE(mut_gen.Mut_register.size() == static_cast<std::size_t>(1));

        REQUIRE(mut_gen.Mut_register.back().size() == static_cast<std::size_t>(nbr_site));
    }

    SECTION("remove_mut")
    {
        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        auto &tskit = singleton_c<tskit_struct_c>::instance();
        tskit.Tskit_output = false;

        mut_genenerator_c mut_gen;

        int node_index = 0;
        int nbr_mut = 1000;
        int nbr_site = 10;
        std::vector<int> ancester_state = std::vector<int>(nbr_site, 1);

        mut_gen.Mod_mut_name = mut_model_enum::jcm;

        //Draw a site
        mut_gen.Site_distri = std::uniform_int_distribution<long int>(0, nbr_site - 1);
        //For each site, the markov chain mutation
        mut_gen.Mut_chain = std::map<int, std::list<int>>();
        //Recap all mut emplacement by branch
        mut_gen.Mut_register = std::list<std::list<int>>{std::list<int>{}};

        mut_gen.add_mut(ancester_state, rand_gen, node_index, nbr_mut);
        mut_gen.remove_mut();

        REQUIRE(mut_gen.Mut_chain.size() == 0);
        REQUIRE(mut_gen.Mut_register.size() == 0);
    }

    SECTION("mut_generation")
    {
        auto &muta_param = singleton_c<muta_param_c>::instance();
        muta_param.Mod_mut_name = mut_model_enum::jcm;
        muta_param.Scaled_mut_rate_theta = 1;

        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.n_total_sample_size = 4;

        std::vector<std::tuple<std::vector<int>, double, int>> gene_tree(7);
        gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0);
        gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0);
        gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0);
        gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0);
        gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1);
        gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1);
        gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, 10, -1);

        int nbr_of_recomb_site_m = 4;
        std::vector<int> ancestry_seq{1, 2, 3, 4};

        std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};
        int MRCA_index = 6;

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        auto &simu_param = singleton_c<simu_param_c>::instance();
        simu_param.Continuous_time_approxim = true;
        rand_gen.put_seed(3);

        std::vector<std::array<int, 5>> result = std::vector<std::array<int, 5>>{static_cast<std::size_t>(nbr_of_recomb_site_m), std::array<int, 5>{0, 0, 0, 0, 0}};

        //jcm => 1/4 of each nucl if mut rate and time large
        int repet = 10000;
        for (int i = 0; i < repet; ++i)
        {
            mut_genenerator_c mut_generator;

            auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, muta_param, MRCA_index, samp_param.n_total_sample_size);

            for (std::size_t indiv = 0; indiv < All_indivs_mut_state.size(); ++indiv)
            {
                for (auto const &mut : All_indivs_mut_state[indiv])
                {
                    ++result[mut.first][mut.second];
                }
            }
        }

        int equi = repet / 4 * samp_param.n_total_sample_size;
        for (std::size_t site = 0; site < result.size(); ++site)
        {
            for (std::size_t nucl = 1; nucl <= 4; ++nucl)
            {
                REQUIRE(result[site][nucl] == Approx(equi).margin(equi * 0.05));
            }
        }
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

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        std::vector<std::tuple<std::vector<int>, double, int>> gene_tree(7);
        gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0);
        gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0);
        gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0);
        gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0);
        gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, -10, 3);
        gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, -10, 3);
        gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, -10, 8);

        int nbr_of_recomb_site_m = 4;
        std::vector<int> ancestry_seq{1, 2, 3, 4};

        std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};
        int MRCA_index = 6;

        std::vector<std::array<int, 5>> result = std::vector<std::array<int, 5>>{static_cast<std::size_t>(nbr_of_recomb_site_m), std::array<int, 5>{0, 0, 0, 0, 0}};

        //jcm => 1/4 of each nucl if mut rate and time large
        int repet = 10000;
        for (int i = 0; i < repet; ++i)
        {
            mut_genenerator_c mut_generator;

            auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, muta_param, MRCA_index, samp_param.n_total_sample_size);

            for (std::size_t indiv = 0; indiv < All_indivs_mut_state.size(); ++indiv)
            {
                for (auto const &mut : All_indivs_mut_state[indiv])
                {
                    ++result[mut.first][mut.second];
                }
            }
        }

        int equi = repet / 4 * samp_param.n_total_sample_size;
        for (std::size_t site = 0; site < result.size(); ++site)
        {
            for (std::size_t nucl = 1; nucl <= 4; ++nucl)
            {
                REQUIRE(result[site][nucl] == Approx(equi).margin(equi * 0.05));
            }
        }
    }

    SECTION("multi_descendant_init")
    {
        auto &simu_param = singleton_c<simu_param_c>::instance();
        simu_param.Continuous_time_approxim = false;

        auto &muta_param = singleton_c<muta_param_c>::instance();
        muta_param.Mod_mut_name = mut_model_enum::jcm;
        muta_param.Unscaled_mut_rate_mu = 1;

        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.n_total_sample_size = 4;

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        std::vector<std::tuple<std::vector<int>, double, int>> gene_tree(6);
        gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0);
        gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0);
        gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0);
        gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0);
        gene_tree[4] = std::make_tuple(std::vector<int>{0, 1, 2}, 0, 3);
        gene_tree[5] = std::make_tuple(std::vector<int>{3, 4}, 0, 8);

        int nbr_of_recomb_site_m = 4;
        std::vector<int> ancestry_seq{1, 2, 3, 4};

        std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};
        int MRCA_index = 5;

        std::vector<std::array<int, 5>> result = std::vector<std::array<int, 5>>{static_cast<std::size_t>(nbr_of_recomb_site_m), std::array<int, 5>{0, 0, 0, 0, 0}};

        //jcm => 1/4 of each nucl if mut rate and time large
        int repet = 10000;
        for (int i = 0; i < repet; ++i)
        {
            mut_genenerator_c mut_generator;

            auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, muta_param, MRCA_index, samp_param.n_total_sample_size);

            for (std::size_t indiv = 0; indiv < All_indivs_mut_state.size(); ++indiv)
            {
                for (auto const &mut : All_indivs_mut_state[indiv])
                {
                    ++result[mut.first][mut.second];
                }
            }
        }

        int equi = repet / 4 * samp_param.n_total_sample_size;
        for (std::size_t site = 0; site < result.size(); ++site)
        {
            for (std::size_t nucl = 1; nucl <= 4; ++nucl)
            {
                REQUIRE(result[site][nucl] == Approx(equi).margin(equi * 0.05));
            }
        }
    }
    // // /******************/

    // SECTION("mut_with_mod_iam")
    // {
    //     std::vector<std::tuple<std::vector<int>, double, int>> gene_tree(7);
    //     gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0);
    //     gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0);
    //     gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0);
    //     gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0);
    //     gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1);
    //     gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1);
    //     gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, 2, -1);

    //     auto &samp_param = singleton_c<samp_param_c>::instance();
    //     samp_param.Sequence_length = 5;
    //     std::vector<int> ancestry_seq(samp_param.Sequence_length, 1);

    //     std::array<int, 2> begin_end_seq_current_tree = {0, samp_param.Sequence_length};

    //     auto &muta_param = singleton_c<muta_param_c>::instance();
    //     muta_param.Mod_mut_name = mut_model_enum::iam;
    //     muta_param.Scaled_mut_rate_theta = 0.5;

    //     samp_param.n_total_sample_size = 4;

    //     auto &rand_gen = singleton_c<rand_gen_c>::instance();
    //     auto &simu_param = singleton_c<simu_param_c>::instance();
    //     simu_param.Continuous_time_approxim = true;
    //     rand_gen.put_seed(3);

    //     int MRCA_index = 6;

    //     int repet = 1;
    //     //site
    //     std::vector<std::map<int, int>> result = std::vector<std::map<int, int>>(samp_param.Sequence_length);
    //     for (auto & ap : result)
    //     {
    //         ap.emplace(1, repet*samp_param.n_total_sample_size );
    //     }

    //     for (int i = 0; i < repet; ++i)
    //     {
    //         mut_genenerator_c mut_generator;

    //         auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, muta_param, MRCA_index, samp_param.n_total_sample_size);

    //         for (std::size_t indiv = 0; indiv < All_indivs_mut_state.size(); ++indiv)
    //         {
    //             for (auto const &mut : All_indivs_mut_state[indiv])
    //             {
    //                 auto temp = result.at(mut.first).emplace(mut.second, 0);
    //                 ++(temp.first->second);
    //                 --result.at(mut.first).at(1);
    //             }
    //         }
    //     }

    //     for (int site = 0; site < samp_param.Sequence_length; ++site)
    //     {
    //         int sum = 0;
    //         std::cout << "##############SITE##############" << std::endl;
    //         for (auto map_state : result.at(site))
    //         {
    //             std::cout << "state : " << map_state.first << ",nbr : " << map_state.second << std::endl;
    //             sum += map_state.second;
    //         }
    //         std::cout << "SUM : " << sum << std::endl;
    //     }

    //     REQUIRE(1 == 1);
    // }

    // SECTION("mut_mod_c_kam")
    // {
    //     std::vector<std::tuple<std::vector<int>, double, int, int>> gene_tree(7);
    //     gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //     gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //     gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //     gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //     gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1, 2);
    //     gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1, 2);
    //     gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, 2, -1, 4);

    //     int nbr_of_recomb_site_m = 3;
    //     std::vector<int> ancestry_seq(3, 5);

    //     auto &muta_param = singleton_c<muta_param_c>::instance();
    //     muta_param.Mod_mut_name = mut_model_enum::kam;
    //     muta_param.Scaled_mut_rate_theta = 1;
    //     muta_param.K_min = 1;
    //     muta_param.K_max = 10;

    //     auto &samp_param = singleton_c<samp_param_c>::instance();
    //     samp_param.n_total_sample_size = 4;

    //     auto &rand_gen = singleton_c<rand_gen_c>::instance();
    //     auto &simu_param = singleton_c<simu_param_c>::instance();
    //     simu_param.Continuous_time_approxim = true;
    //     rand_gen.put_seed(3);

    //     int MRCA_index = 6;

    //     mut_genenerator_c mut_generator;

    //     std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};

    //     auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, MRCA_index);

    //     REQUIRE(mut_generator.All_indivs_mut_state[0].at(0) == 7);
    //     REQUIRE(mut_generator.All_indivs_mut_state[0].at(1) == 6);
    //     REQUIRE(mut_generator.All_indivs_mut_state[0].at(2) == 8);
    //     REQUIRE(mut_generator.All_indivs_mut_state[0].size() == 3);

    //     REQUIRE(mut_generator.All_indivs_mut_state[1].at(0) == 9);
    //     REQUIRE(mut_generator.All_indivs_mut_state[1].at(1) == 4);
    //     REQUIRE(mut_generator.All_indivs_mut_state[1].at(2) == 8);
    //     REQUIRE(mut_generator.All_indivs_mut_state[1].size() == 3);

    //     REQUIRE(mut_generator.All_indivs_mut_state[2].at(0) == 8);
    //     REQUIRE(mut_generator.All_indivs_mut_state[2].at(1) == 10);
    //     REQUIRE(mut_generator.All_indivs_mut_state[2].at(2) == 8);
    //     REQUIRE(mut_generator.All_indivs_mut_state[2].size() == 3);

    //     REQUIRE(mut_generator.All_indivs_mut_state[3].at(0) == 2);
    //     REQUIRE(mut_generator.All_indivs_mut_state[3].at(1) == 6);
    //     REQUIRE(mut_generator.All_indivs_mut_state[3].at(2) == 4);
    //     REQUIRE(mut_generator.All_indivs_mut_state[3].size() == 3);
    // }

    // SECTION("mut_mod_c_smm")
    // {
    //     std::vector<std::tuple<std::vector<int>, double, int, int>> gene_tree(7);
    //     gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //     gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //     gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //     gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //     gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1, 2);
    //     gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1, 2);
    //     gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, 2, -1, 4);

    //     int nbr_of_recomb_site_m = 3;
    //     std::vector<int> ancestry_seq(3, 5);

    //     std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};

    //     auto &muta_param = singleton_c<muta_param_c>::instance();
    //     muta_param.Mod_mut_name = mut_model_enum::smm;
    //     muta_param.Scaled_mut_rate_theta = 1;
    //     muta_param.K_min = 4;
    //     muta_param.K_max = 6;

    //     auto &samp_param = singleton_c<samp_param_c>::instance();
    //     samp_param.n_total_sample_size = 4;

    //     auto &rand_gen = singleton_c<rand_gen_c>::instance();
    //     auto &simu_param = singleton_c<simu_param_c>::instance();
    //     simu_param.Continuous_time_approxim = true;
    //     rand_gen.put_seed(3);

    //     int MRCA_index = 6;

    //     mut_genenerator_c mut_generator;

    //     auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, MRCA_index);

    //     REQUIRE(mut_generator.All_indivs_mut_state[0].at(0) == 4);
    //     REQUIRE(mut_generator.All_indivs_mut_state[0].at(1) == 4);
    //     REQUIRE(mut_generator.All_indivs_mut_state[0].at(2) == 6);
    //     REQUIRE(mut_generator.All_indivs_mut_state[0].size() == 3);

    //     REQUIRE(mut_generator.All_indivs_mut_state[1].at(1) == 4);
    //     REQUIRE(mut_generator.All_indivs_mut_state[1].at(2) == 5);
    //     REQUIRE(mut_generator.All_indivs_mut_state[1].size() == 2);

    //     REQUIRE(mut_generator.All_indivs_mut_state[2].at(0) == 5);
    //     REQUIRE(mut_generator.All_indivs_mut_state[2].at(1) == 4);
    //     REQUIRE(mut_generator.All_indivs_mut_state[2].at(2) == 4);
    //     REQUIRE(mut_generator.All_indivs_mut_state[2].size() == 3);

    //     REQUIRE(mut_generator.All_indivs_mut_state[3].at(0) == 6);
    //     REQUIRE(mut_generator.All_indivs_mut_state[3].at(1) == 5);
    //     REQUIRE(mut_generator.All_indivs_mut_state[3].at(2) == 5);
    //     REQUIRE(mut_generator.All_indivs_mut_state[3].size() == 3);
    // }
    // /******************/
    // SECTION("mut_mod_c_gsm")
    // {
    //     std::vector<std::tuple<std::vector<int>, double, int, int>> gene_tree(7);
    //     gene_tree[0] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //     gene_tree[1] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //     gene_tree[2] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //     gene_tree[3] = std::make_tuple(std::vector<int>(), 0, 0, 1);
    //     gene_tree[4] = std::make_tuple(std::vector<int>{0, 1}, 1, -1, 2);
    //     gene_tree[5] = std::make_tuple(std::vector<int>{2, 3}, 1, -1, 2);
    //     gene_tree[6] = std::make_tuple(std::vector<int>{4, 5}, 2, -1, 4);

    //     int nbr_of_recomb_site_m = 3;
    //     std::vector<int> ancestry_seq(3, 5);

    //     std::array<int, 2> begin_end_seq_current_tree = {0, nbr_of_recomb_site_m};

    //     auto &muta_param = singleton_c<muta_param_c>::instance();
    //     muta_param.Mod_mut_name = mut_model_enum::gsm;
    //     muta_param.Scaled_mut_rate_theta = 1;
    //     muta_param.K_min = 1;
    //     muta_param.K_max = 10;
    //     muta_param.P_gsm = 0.5;

    //     auto &samp_param = singleton_c<samp_param_c>::instance();
    //     samp_param.n_total_sample_size = 4;

    //     auto &rand_gen = singleton_c<rand_gen_c>::instance();
    //     auto &simu_param = singleton_c<simu_param_c>::instance();
    //     simu_param.Continuous_time_approxim = true;
    //     rand_gen.put_seed(3);

    //     int MRCA_index = 6;

    //     mut_genenerator_c mut_generator;

    //     auto All_indivs_mut_state = mut_generator.mut_generation(gene_tree, ancestry_seq, begin_end_seq_current_tree, MRCA_index);

    //     REQUIRE(mut_generator.All_indivs_mut_state[0].at(0) == 3);
    //     REQUIRE(mut_generator.All_indivs_mut_state[0].at(1) == 8);
    //     REQUIRE(mut_generator.All_indivs_mut_state[0].size() == 2);

    //     REQUIRE(mut_generator.All_indivs_mut_state[1].at(0) == 3);
    //     REQUIRE(mut_generator.All_indivs_mut_state[1].at(1) == 8);
    //     REQUIRE(mut_generator.All_indivs_mut_state[1].size() == 2);

    //     REQUIRE(mut_generator.All_indivs_mut_state[2].at(0) == 8);
    //     REQUIRE(mut_generator.All_indivs_mut_state[2].at(1) == 7);
    //     REQUIRE(mut_generator.All_indivs_mut_state[2].at(2) == 9);
    //     REQUIRE(mut_generator.All_indivs_mut_state[2].size() == 3);

    //     REQUIRE(mut_generator.All_indivs_mut_state[3].at(0) == 8);
    //     REQUIRE(mut_generator.All_indivs_mut_state[3].at(1) == 9);
    //     REQUIRE(mut_generator.All_indivs_mut_state[3].at(2) == 4);
    //     REQUIRE(mut_generator.All_indivs_mut_state[3].size() == 3);
    // }

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