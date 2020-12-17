
#define CATCH_CONFIG_MAIN HudsonHTest
#include "catch.hpp"
#include <iostream>

#include "kelleher_algo.hpp"

TEST_CASE()
{
    REQUIRE(1);
}

TEST_CASE("algo_h_monolocus_test")
{
    SECTION("haploid_init")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                    {0, 0},
                                    {0, 0},
                                    {1, 1}};
        samp_param.Ploidy = 1;

        samp_param.n_total_sample_size = 4;
        samp_param.Chr_nbr = 1;
        samp_param.Sequence_length = 10;

        struct_arg_c main;
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Population_size_N = samp_param.n_total_sample_size;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        indiv_vec.ini(lat, samp_param);
        auto chr_index = 0;
        main.ini(indiv_vec, samp_param, chr_index);

        REQUIRE(indiv_vec.size() == static_cast<std::size_t>(samp_param.n_total_sample_size));

        REQUIRE(indiv_vec[0]->Ident == 0);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 1) == nullptr);
        REQUIRE((indiv_vec[0]->get_chr(chr_index, 0))->U_ident_node == 0);
        REQUIRE((indiv_vec[0]->get_chr(chr_index, 0))->L_left_brkpt == 0);
        REQUIRE((indiv_vec[0]->get_chr(chr_index, 0))->R_right_brkpt == 10);
        REQUIRE((indiv_vec[0]->get_chr(chr_index, 0))->get_indiv() == indiv_vec[0]);

        REQUIRE(indiv_vec[samp_param.n_total_sample_size - 1]->Ident == samp_param.n_total_sample_size - 1);
        REQUIRE(indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 1) == nullptr);
        REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 0))->U_ident_node == samp_param.n_total_sample_size - 1);
        REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 0))->L_left_brkpt == 0);
        REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 0))->R_right_brkpt == 10);
        REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 0))->get_indiv() == indiv_vec[samp_param.n_total_sample_size - 1]);

        REQUIRE(main.L_cumul_nbr_brkpt_per_seg.get_frequency_index(0) == samp_param.Sequence_length - 1);
        REQUIRE(main.L_cumul_nbr_brkpt_per_seg.get_frequency_index(samp_param.n_total_sample_size - 1) == samp_param.Sequence_length - 1);

        REQUIRE(lat.Lattice[0].Indivs_in_pop.size() == static_cast<std::size_t>(samp_param.n_total_sample_size - 1));

        auto iter = lat.Lattice[0].Indivs_in_pop.begin();
        REQUIRE(std::get<1>(*iter) == indiv_vec[0]);
        REQUIRE(std::get<1>(*++iter) == indiv_vec[1]);
        REQUIRE(std::get<1>(*++iter) == indiv_vec[2]);

        REQUIRE(lat.Lattice[3].Indivs_in_pop.size() == 1);
        REQUIRE(std::get<1>(*(lat.Lattice[3].Indivs_in_pop.begin())) == indiv_vec[3]);
        REQUIRE(*(lat.Nodes_with_lineage.begin()) == &(lat.Lattice[0]));

        REQUIRE(main.S_intersection_count.find(0)->second == samp_param.n_total_sample_size);
        REQUIRE(main.S_intersection_count.find(samp_param.Sequence_length)->second == -1);

        REQUIRE(main.W_next_node_nbr == samp_param.n_total_sample_size);
        REQUIRE(main.T_time == 0);
    }

    SECTION("diploid_init")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                    {0, 0},
                                    {0, 0},
                                    {1, 1}};

        samp_param.n_total_sample_size = 4;
        samp_param.Chr_nbr = 1;
        samp_param.Sequence_length = 10;

        struct_arg_c main;
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Population_size_N = samp_param.n_total_sample_size;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        indiv_vec.ini(lat, samp_param);
        samp_param.Ploidy = 2;
        auto chr_index = 0;
        main.ini(indiv_vec, samp_param, chr_index);

        REQUIRE(indiv_vec.size() == static_cast<std::size_t>(samp_param.n_total_sample_size));

        REQUIRE(indiv_vec[0]->Ident == 0);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 1) != nullptr);

        REQUIRE((indiv_vec[0]->get_chr(chr_index, 0))->U_ident_node == 0);
        REQUIRE((indiv_vec[0]->get_chr(chr_index, 0))->L_left_brkpt == 0);
        REQUIRE((indiv_vec[0]->get_chr(chr_index, 0))->R_right_brkpt == 10);
        REQUIRE((indiv_vec[0]->get_chr(chr_index, 0))->get_indiv() == indiv_vec[0]);

        REQUIRE((indiv_vec[0]->get_chr(chr_index, 1))->U_ident_node == 1);
        REQUIRE((indiv_vec[0]->get_chr(chr_index, 1))->L_left_brkpt == 0);
        REQUIRE((indiv_vec[0]->get_chr(chr_index, 1))->R_right_brkpt == 10);
        REQUIRE((indiv_vec[0]->get_chr(chr_index, 1))->get_indiv() == indiv_vec[0]);

        REQUIRE(indiv_vec[samp_param.n_total_sample_size - 1]->Ident == samp_param.n_total_sample_size - 1);
        REQUIRE(indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 1) != nullptr);

        REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 0))->U_ident_node == 2 * samp_param.n_total_sample_size - 2);
        REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 0))->L_left_brkpt == 0);
        REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 0))->R_right_brkpt == 10);
        REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 0))->get_indiv() == indiv_vec[samp_param.n_total_sample_size - 1]);

        REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 1))->U_ident_node == 2 * samp_param.n_total_sample_size - 1);
        REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 1))->L_left_brkpt == 0);
        REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 1))->R_right_brkpt == 10);
        REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 1))->get_indiv() == indiv_vec[samp_param.n_total_sample_size - 1]);

        REQUIRE(main.L_cumul_nbr_brkpt_per_seg.get_frequency_index(0) == samp_param.Sequence_length - 1);
        REQUIRE(main.L_cumul_nbr_brkpt_per_seg.get_frequency_index(2 * samp_param.n_total_sample_size - 1) == samp_param.Sequence_length - 1);

        REQUIRE(lat.Lattice[0].Indivs_in_pop.size() == static_cast<std::size_t>(samp_param.n_total_sample_size - 1));

        auto iter = lat.Lattice[0].Indivs_in_pop.begin();
        REQUIRE(std::get<1>(*iter) == indiv_vec[0]);
        REQUIRE(std::get<1>(*++iter) == indiv_vec[1]);
        REQUIRE(std::get<1>(*++iter) == indiv_vec[2]);

        REQUIRE(lat.Lattice[3].Indivs_in_pop.size() == 1);
        REQUIRE(std::get<1>(*(lat.Lattice[3].Indivs_in_pop.begin())) == indiv_vec[3]);
        REQUIRE(*(lat.Nodes_with_lineage.begin()) == &(lat.Lattice[0]));

        REQUIRE(main.S_intersection_count.find(0)->second == 2 * samp_param.n_total_sample_size);
        REQUIRE(main.S_intersection_count.find(samp_param.Sequence_length)->second == -1);

        REQUIRE(main.W_next_node_nbr == 2 * samp_param.n_total_sample_size);
        REQUIRE(main.T_time == 0);
    }

    SECTION("choose_event")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                    {0, 0},
                                    {0, 0},
                                    {0, 0}};
        samp_param.Ploidy = 1;

        samp_param.n_total_sample_size = 4;
        samp_param.Chr_nbr = 1;
        samp_param.Sequence_length = 10;

        struct_arg_c main;
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {0, 0};
        demo_param.Population_size_N = samp_param.n_total_sample_size;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        indiv_vec.ini(lat, samp_param);
        auto chr_index = 0;
        main.ini(indiv_vec, samp_param, chr_index);
        //TODO : regarder les prob ;p
        REQUIRE(!(main.choose_event(0, indiv_vec, rand_gen)));
    }
}

TEST_CASE("algo_h_multilocus_test")
{
    SECTION("haploid_init")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                    {0, 0},
                                    {0, 0},
                                    {1, 1}};
        samp_param.Ploidy = 1;

        samp_param.n_total_sample_size = 4;
        samp_param.Chr_nbr = 5;
        samp_param.Sequence_length = 10;

        indiv_stock_c indiv_vec(samp_param.Chr_nbr);

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Population_size_N = samp_param.n_total_sample_size;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        indiv_vec.ini(lat, samp_param);

        for (int chr_index = 0; chr_index < samp_param.Chr_nbr; ++chr_index)
        {
            struct_arg_c main;
            main.ini(indiv_vec, samp_param, chr_index);

            REQUIRE(indiv_vec.size() == static_cast<std::size_t>(samp_param.n_total_sample_size));

            REQUIRE(indiv_vec[0]->Ident == 0);
            REQUIRE(indiv_vec[0]->get_chr(chr_index, 1) == nullptr);
            REQUIRE((indiv_vec[0]->get_chr(chr_index, 0))->U_ident_node == 0);
            REQUIRE((indiv_vec[0]->get_chr(chr_index, 0))->L_left_brkpt == 0);
            REQUIRE((indiv_vec[0]->get_chr(chr_index, 0))->R_right_brkpt == 10);
            REQUIRE((indiv_vec[0]->get_chr(chr_index, 0))->get_indiv() == indiv_vec[0]);

            REQUIRE(indiv_vec[samp_param.n_total_sample_size - 1]->Ident == samp_param.n_total_sample_size - 1);
            REQUIRE(indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 1) == nullptr);
            REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 0))->U_ident_node == samp_param.n_total_sample_size - 1);
            REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 0))->L_left_brkpt == 0);
            REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 0))->R_right_brkpt == 10);
            REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 0))->get_indiv() == indiv_vec[samp_param.n_total_sample_size - 1]);

            REQUIRE(main.L_cumul_nbr_brkpt_per_seg.get_frequency_index(0) == samp_param.Sequence_length - 1);
            REQUIRE(main.L_cumul_nbr_brkpt_per_seg.get_frequency_index(samp_param.n_total_sample_size - 1) == samp_param.Sequence_length - 1);

            REQUIRE(main.S_intersection_count.find(0)->second == samp_param.n_total_sample_size);
            REQUIRE(main.S_intersection_count.find(samp_param.Sequence_length)->second == -1);

            REQUIRE(main.W_next_node_nbr == samp_param.n_total_sample_size);
            REQUIRE(main.T_time == 0);
        }

        REQUIRE(lat.Lattice[0].Indivs_in_pop.size() == static_cast<std::size_t>(samp_param.n_total_sample_size - 1));

        auto iter = lat.Lattice[0].Indivs_in_pop.begin();
        REQUIRE(std::get<1>(*iter) == indiv_vec[0]);
        REQUIRE(std::get<1>(*++iter) == indiv_vec[1]);
        REQUIRE(std::get<1>(*++iter) == indiv_vec[2]);

        REQUIRE(lat.Lattice[3].Indivs_in_pop.size() == 1);
        REQUIRE(std::get<1>(*(lat.Lattice[3].Indivs_in_pop.begin())) == indiv_vec[3]);
        REQUIRE(*(lat.Nodes_with_lineage.begin()) == &(lat.Lattice[0]));
    }

    SECTION("diploid_init")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                    {0, 0},
                                    {0, 0},
                                    {1, 1}};

        samp_param.n_total_sample_size = 4;
        samp_param.Chr_nbr = 1;
        samp_param.Sequence_length = 10;

        indiv_stock_c indiv_vec(samp_param.Chr_nbr);

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Population_size_N = samp_param.n_total_sample_size;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        indiv_vec.ini(lat, samp_param);
        samp_param.Ploidy = 2;

        for (int chr_index = 0; chr_index < samp_param.Chr_nbr; ++chr_index)
        {
            struct_arg_c main;
            main.ini(indiv_vec, samp_param, chr_index);

            REQUIRE(indiv_vec.size() == static_cast<std::size_t>(samp_param.n_total_sample_size));

            REQUIRE(indiv_vec[0]->Ident == 0);
            REQUIRE(indiv_vec[0]->get_chr(chr_index, 1) != nullptr);

            REQUIRE((indiv_vec[0]->get_chr(chr_index, 0))->U_ident_node == 0);
            REQUIRE((indiv_vec[0]->get_chr(chr_index, 0))->L_left_brkpt == 0);
            REQUIRE((indiv_vec[0]->get_chr(chr_index, 0))->R_right_brkpt == 10);
            REQUIRE((indiv_vec[0]->get_chr(chr_index, 0))->get_indiv() == indiv_vec[0]);

            REQUIRE((indiv_vec[0]->get_chr(chr_index, 1))->U_ident_node == 1);
            REQUIRE((indiv_vec[0]->get_chr(chr_index, 1))->L_left_brkpt == 0);
            REQUIRE((indiv_vec[0]->get_chr(chr_index, 1))->R_right_brkpt == 10);
            REQUIRE((indiv_vec[0]->get_chr(chr_index, 1))->get_indiv() == indiv_vec[0]);

            REQUIRE(indiv_vec[samp_param.n_total_sample_size - 1]->Ident == samp_param.n_total_sample_size - 1);
            REQUIRE(indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 1) != nullptr);

            REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 0))->U_ident_node == 2 * samp_param.n_total_sample_size - 2);
            REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 0))->L_left_brkpt == 0);
            REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 0))->R_right_brkpt == 10);
            REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 0))->get_indiv() == indiv_vec[samp_param.n_total_sample_size - 1]);

            REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 1))->U_ident_node == 2 * samp_param.n_total_sample_size - 1);
            REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 1))->L_left_brkpt == 0);
            REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 1))->R_right_brkpt == 10);
            REQUIRE((indiv_vec[samp_param.n_total_sample_size - 1]->get_chr(chr_index, 1))->get_indiv() == indiv_vec[samp_param.n_total_sample_size - 1]);

            REQUIRE(main.L_cumul_nbr_brkpt_per_seg.get_frequency_index(0) == samp_param.Sequence_length - 1);
            REQUIRE(main.L_cumul_nbr_brkpt_per_seg.get_frequency_index(2 * samp_param.n_total_sample_size - 1) == samp_param.Sequence_length - 1);

            REQUIRE(main.S_intersection_count.find(0)->second == 2 * samp_param.n_total_sample_size);
            REQUIRE(main.S_intersection_count.find(samp_param.Sequence_length)->second == -1);

            REQUIRE(main.W_next_node_nbr == 2 * samp_param.n_total_sample_size);
            REQUIRE(main.T_time == 0);
        }

        REQUIRE(lat.Lattice[0].Indivs_in_pop.size() == static_cast<std::size_t>(samp_param.n_total_sample_size - 1));

        auto iter = lat.Lattice[0].Indivs_in_pop.begin();
        REQUIRE(std::get<1>(*iter) == indiv_vec[0]);
        REQUIRE(std::get<1>(*++iter) == indiv_vec[1]);
        REQUIRE(std::get<1>(*++iter) == indiv_vec[2]);

        REQUIRE(lat.Lattice[3].Indivs_in_pop.size() == 1);
        REQUIRE(std::get<1>(*(lat.Lattice[3].Indivs_in_pop.begin())) == indiv_vec[3]);
        REQUIRE(*(lat.Nodes_with_lineage.begin()) == &(lat.Lattice[0]));
    }
}
