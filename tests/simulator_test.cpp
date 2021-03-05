#define CATCH_CONFIG_MAIN SimulatorTest
#include "catch.hpp"
#include <iostream>

#include "simulator.hpp"

TEST_CASE()
{
    REQUIRE(1);
}

TEST_CASE("simulator_test")
{
    SECTION("apply_first_recomb")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};
        samp_param.Ploidy = 1;
        samp_param.n_total_sample_size = 10;
        samp_param.Sequence_length = 100;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

        std::vector<struct_arg_c> struct_arg_vec(1);
        auto &struct_arg = struct_arg_vec.at(0);
        ;
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = samp_param.n_total_sample_size;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        indiv_vec.ini(lat, samp_param);
        struct_arg.ini(indiv_vec, samp_param, chr_index);

        std::uniform_int_distribution<long int> uniform_distribution(1, struct_arg.L_cumul_nbr_brkpt_per_seg.get_tot_cumul_freq());
        long int h_choosen_brkpt = uniform_distribution(rand_gen.Seed_gen);
        int choosen_seg_index = struct_arg.L_cumul_nbr_brkpt_per_seg.seg_index(h_choosen_brkpt);
        auto y_choosen_seg = struct_arg.Pool_seg.get_obj_from_index(choosen_seg_index);

        auto k_choose_brkpt = y_choosen_seg->R_right_brkpt - struct_arg.L_cumul_nbr_brkpt_per_seg.get_cumul_freq_at_index(choosen_seg_index) + h_choosen_brkpt - 1;

        //Random draw for brkpt number : 559
        auto new_seg = apply_first_recomb(struct_arg.Pool_seg, struct_arg.L_cumul_nbr_brkpt_per_seg, y_choosen_seg, k_choose_brkpt);
        auto indiv = new indiv_c(samp_param.Chr_nbr, lat.node({0, 0}));
        indiv->update_indiv(chr_index, new_seg.at(1), nullptr);
        indiv_vec.add_indiv(indiv);

        REQUIRE(indiv_vec.size() == 11);

        REQUIRE(indiv_vec[5]->get_chr(chr_index, 0)->L_left_brkpt == 0);
        REQUIRE(indiv_vec[5]->get_chr(chr_index, 0)->R_right_brkpt == 59);
        REQUIRE(indiv_vec[5]->get_chr(chr_index, 0)->U_ident_node == 5);

        REQUIRE(indiv_vec[10]->get_chr(chr_index, 0)->L_left_brkpt == 59);
        REQUIRE(indiv_vec[10]->get_chr(chr_index, 0)->R_right_brkpt == 100);
        REQUIRE(indiv_vec[10]->get_chr(chr_index, 0)->U_ident_node == 5);

        REQUIRE(struct_arg.S_intersection_count.size() == 2);
    }

    SECTION("algo_coa")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};

        samp_param.n_total_sample_size = 10;
        samp_param.Ploidy = 1;
        samp_param.Sequence_length = 100;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

        std::vector<struct_arg_c> struct_arg_vec(1);
        auto &struct_arg = struct_arg_vec.at(0);
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        coa_table_c coa_table;

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = samp_param.n_total_sample_size;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        indiv_vec.ini(lat, samp_param);
        struct_arg.ini(indiv_vec, samp_param, chr_index);

        //Random draw for coalescence => 5 and 1
        int num_indiv_1, num_indiv_2;
        std::uniform_int_distribution<long int> uniform_distribution(0, indiv_vec.size() - 1);
        num_indiv_1 = uniform_distribution(rand_gen.Seed_gen);
        num_indiv_2 = uniform_distribution(rand_gen.Seed_gen);

        //Choose two differents seg
        while (num_indiv_1 == num_indiv_2)
        {
            num_indiv_2 = uniform_distribution(rand_gen.Seed_gen);
        }

        //Choose the num_indiv_1-iem node in the map indiv_vec
        auto seg1 = indiv_vec[num_indiv_1]->get_chr(chr_index, 0);
        //Choose the num_indiv_2-iem node in the map indiv_vec
        auto seg2 = indiv_vec[num_indiv_2]->get_chr(chr_index, 0);

        auto new_seg = apply_coa_two_lineage(coa_table, struct_arg.Pool_seg, struct_arg.L_cumul_nbr_brkpt_per_seg, struct_arg.S_intersection_count, struct_arg.W_next_node_nbr, struct_arg.T_time, -1, seg1, seg2);

        //No update for Indiv_vec in coa
        REQUIRE(indiv_vec.size() == 10);

        REQUIRE(new_seg->L_left_brkpt == 0);
        REQUIRE(new_seg->R_right_brkpt == 100);
        REQUIRE(new_seg->U_ident_node == 10);

        auto coa_tab_event = coa_table.get_coa_table().at(0);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_tab_event) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_tab_event) == 100);
        REQUIRE(coa_table_c::get_num_node_u(coa_tab_event) == 10);
        REQUIRE((coa_table_c::get_childs_node_c(coa_tab_event))[0] == 5);
        REQUIRE((coa_table_c::get_childs_node_c(coa_tab_event))[1] == 1);
        REQUIRE(coa_table_c::get_time_t(coa_tab_event) == 0);
        REQUIRE(coa_table_c::get_gen(coa_tab_event) == -1);

        REQUIRE(struct_arg.S_intersection_count.size() == 2);
    }

    SECTION("ditrib_and_coa_lign_btw_1_chr")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};

        samp_param.n_total_sample_size = 10;
        samp_param.Ploidy = 1;
        samp_param.Sequence_length = 100;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

        std::vector<struct_arg_c> struct_arg_vec(1);
        auto &struct_arg = struct_arg_vec.at(0);
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(1);

        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = samp_param.n_total_sample_size;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        indiv_vec.ini(lat, samp_param);
        struct_arg.ini(indiv_vec, samp_param, chr_index);

        auto indiv_1 = indiv_vec[1];
        auto indiv_2 = indiv_vec[2];

        ditrib_and_coa_lign_btw_1_chr(coa_table_vec, indiv_vec, struct_arg_vec, -1, indiv_1, indiv_2);

        REQUIRE(indiv_vec.size() == 9);

        REQUIRE(indiv_vec[1]->get_chr(chr_index, 0)->L_left_brkpt == 0);
        REQUIRE(indiv_vec[1]->get_chr(chr_index, 0)->R_right_brkpt == 100);
        REQUIRE(indiv_vec[1]->get_chr(chr_index, 0)->U_ident_node == 10);
        REQUIRE(indiv_vec[1]->get_chr(chr_index, 0)->get_indiv() == indiv_vec[1]);

        REQUIRE(indiv_vec[1]->get_chr(chr_index, 1) == nullptr);

        auto coa_tab_event = coa_table.get_coa_table().at(0);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_tab_event) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_tab_event) == 100);
        REQUIRE(coa_table_c::get_num_node_u(coa_tab_event) == 10);
        REQUIRE((coa_table_c::get_childs_node_c(coa_tab_event))[0] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(coa_tab_event))[1] == 2);
        REQUIRE(coa_table_c::get_time_t(coa_tab_event) == 0);
        REQUIRE(coa_table_c::get_gen(coa_tab_event) == -1);

        REQUIRE(struct_arg.S_intersection_count.size() == 2);
        REQUIRE((struct_arg.S_intersection_count.begin())->second == 9);
    }

    SECTION("ditrib_and_coa_lign_btw_1_chr_simple_coa")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 0},
                                       {0, 0}};

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(1);

        int gen = 1;
        samp_param.Ploidy = 1;
        samp_param.n_total_sample_size = 3;
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 3;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        std::vector<struct_arg_c> algo_gen_vec(1);
        auto &algo_gen = algo_gen_vec.at(0);
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        indiv_vec.ini(lat, samp_param);
        algo_gen.ini(indiv_vec, samp_param, chr_index);

        auto new_cell = ditrib_and_coa_lign_btw_1_chr(coa_table_vec, indiv_vec, algo_gen_vec, gen, indiv_vec[0], indiv_vec[1]);

        REQUIRE(new_cell == indiv_vec[0]);

        REQUIRE(coa_table.size() == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 3);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 1);
        REQUIRE(coa_table_c::get_gen(coa_table[0]) == 1);

        REQUIRE(algo_gen.S_intersection_count.size() == 2);
        REQUIRE(algo_gen.S_intersection_count.at(0) == 2);
        REQUIRE(algo_gen.S_intersection_count.at(1) == -1);
    }

    SECTION("ditrib_and_coa_lign_btw_1_chr_mrca")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        //Two independante coa one with mrca and the other with no real coa    auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 0}};

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(1);

        int gen = 1;
        samp_param.Ploidy = 1;
        samp_param.n_total_sample_size = 2;
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        std::vector<struct_arg_c> algo_gen_vec(1);
        auto &algo_gen = algo_gen_vec.at(0);
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        indiv_vec.ini(lat, samp_param);
        algo_gen.ini(indiv_vec, samp_param, chr_index);

        auto new_cell = ditrib_and_coa_lign_btw_1_chr(coa_table_vec, indiv_vec, algo_gen_vec, gen, indiv_vec[0], indiv_vec[1]);

        REQUIRE(new_cell == nullptr);
        REQUIRE(indiv_vec.size() == 0);

        REQUIRE(coa_table.size() == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 1);
        REQUIRE(coa_table_c::get_gen(coa_table[0]) == 1);

        REQUIRE(algo_gen.S_intersection_count.size() == 2);
        REQUIRE(algo_gen.S_intersection_count.at(0) == 0);
        REQUIRE(algo_gen.S_intersection_count.at(1) == -1);
    }

    SECTION("ditrib_and_coa_lign_btw_1_chr_after_partial_mrca")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        //Two independante coa one with mrca and the other with no real coa
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 0},
                                       {0, 0}};

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(1);

        int gen = 1;
        samp_param.Ploidy = 1;
        samp_param.n_total_sample_size = 3;
        samp_param.Sequence_length = 2;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        std::vector<struct_arg_c> algo_gen_vec(1);
        auto &algo_gen = algo_gen_vec.at(0);
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        indiv_vec.ini(lat, samp_param);
        algo_gen.ini(indiv_vec, samp_param, chr_index);

        //Bidouillage ON

        indiv_vec[0]->get_chr(chr_index, 0)->R_right_brkpt = 1;
        indiv_vec[1]->get_chr(chr_index, 0)->R_right_brkpt = 1;

        algo_gen.S_intersection_count.at(0) = 2;
        algo_gen.S_intersection_count.emplace(1, 3);

        //Bidouillage OFF

        auto new_cell = ditrib_and_coa_lign_btw_1_chr(coa_table_vec, indiv_vec, algo_gen_vec, gen, indiv_vec[0], indiv_vec[1]);

        new_cell = ditrib_and_coa_lign_btw_1_chr(coa_table_vec, indiv_vec, algo_gen_vec, gen, new_cell, indiv_vec[2]);

        REQUIRE(new_cell == indiv_vec[2]);

        REQUIRE(coa_table.size() == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 3);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 1);
        REQUIRE(coa_table_c::get_gen(coa_table[0]) == 1);

        REQUIRE(algo_gen.S_intersection_count.size() == 3);
        REQUIRE(algo_gen.S_intersection_count.at(0) == 0);
        REQUIRE(algo_gen.S_intersection_count.at(1) == 3);
        REQUIRE(algo_gen.S_intersection_count.at(2) == -1);
    }

    SECTION("ditrib_and_coa_lign_btw_1_chr_multi_coa")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 0},
                                       {0, 0},
                                       {0, 0}};

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(1);

        int gen = 1;
        samp_param.Ploidy = 1;
        samp_param.n_total_sample_size = 4;
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        std::vector<struct_arg_c> algo_gen_vec(1);
        auto &algo_gen = algo_gen_vec.at(0);
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        indiv_vec.ini(lat, samp_param);
        algo_gen.ini(indiv_vec, samp_param, chr_index);

        auto new_cell = ditrib_and_coa_lign_btw_1_chr(coa_table_vec, indiv_vec, algo_gen_vec, gen, indiv_vec[0], indiv_vec[1]);

        new_cell = ditrib_and_coa_lign_btw_1_chr(coa_table_vec, indiv_vec, algo_gen_vec, gen, new_cell, indiv_vec[1]);

        REQUIRE(new_cell == indiv_vec[0]);
        REQUIRE(coa_table.size() == 2);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 4);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 1);
        REQUIRE(coa_table_c::get_gen(coa_table[0]) == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[1]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[1]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[1]) == 5);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[0] == 4);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[1] == 3);
        REQUIRE(coa_table_c::get_gen(coa_table[1]) == 1);

        REQUIRE(algo_gen.S_intersection_count.size() == 2);
        REQUIRE(algo_gen.S_intersection_count.at(0) == 2);
        REQUIRE(algo_gen.S_intersection_count.at(1) == -1);
    }

    SECTION("ditrib_and_coa_lign_btw_2_chr")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0}, {0, 0}};

        samp_param.n_total_sample_size = 2;
        samp_param.Ploidy = 2;
        samp_param.Sequence_length = 100;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

        std::vector<struct_arg_c> struct_arg_vec(1);
        auto &struct_arg = struct_arg_vec.at(0);
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = samp_param.n_total_sample_size;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        indiv_vec.ini(lat, samp_param);
        struct_arg.ini(indiv_vec, samp_param, chr_index);

        auto num_indiv_1 = indiv_vec[0];
        auto num_indiv_2 = indiv_vec[1];

        REQUIRE(struct_arg.S_intersection_count.at(0) == 4);

        ditrib_and_coa_lign_btw_2_chr(coa_table_vec, indiv_vec, struct_arg_vec, 1, num_indiv_1, num_indiv_2);

        //Total mrca
        REQUIRE(indiv_vec.size() == 1);

        REQUIRE(indiv_vec[0]->get_chr(chr_index, 0)->L_left_brkpt == 0);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 0)->R_right_brkpt == 100);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 0)->U_ident_node == 4);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 0)->get_indiv() == indiv_vec[0]);

        REQUIRE(indiv_vec[0]->get_chr(chr_index, 1)->L_left_brkpt == 0);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 1)->R_right_brkpt == 100);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 1)->U_ident_node == 5);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 1)->get_indiv() == indiv_vec[0]);

        REQUIRE(coa_table.get_coa_table().size() == 2);

        REQUIRE(struct_arg.S_intersection_count.size() == 2);
        REQUIRE(struct_arg.S_intersection_count.at(0) == 2);
        REQUIRE(struct_arg.S_intersection_count.at(100) == -1);
    }

    SECTION("ditrib_and_coa_lign_btw_2_chr_simple_coa")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 0}};

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(1);

        int gen = 1;
        samp_param.Ploidy = 2;
        samp_param.n_total_sample_size = 2;
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 3;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        std::vector<struct_arg_c> algo_gen_vec(1);
        auto &algo_gen = algo_gen_vec.at(0);
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        indiv_vec.ini(lat, samp_param);
        algo_gen.ini(indiv_vec, samp_param, chr_index);

        auto new_cell = ditrib_and_coa_lign_btw_2_chr(coa_table_vec, indiv_vec, algo_gen_vec, gen, indiv_vec[0], indiv_vec[1]);

        REQUIRE(new_cell == indiv_vec[0]);

        REQUIRE(coa_table.size() == 2);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 4);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 2);
        REQUIRE(coa_table_c::get_gen(coa_table[0]) == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[1]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[1]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[1]) == 5);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[0] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[1] == 3);
        REQUIRE(coa_table_c::get_gen(coa_table[1]) == 1);

        REQUIRE(algo_gen.S_intersection_count.size() == 2);
        REQUIRE(algo_gen.S_intersection_count.at(0) == 2);
        REQUIRE(algo_gen.S_intersection_count.at(1) == -1);
    }

    SECTION("ditrib_and_coa_lign_btw_2_chr_partial_mrca")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        //Two independante coa one with mrca and the other with no real coa    auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 0}};

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(1);

        int gen = 1;
        samp_param.Ploidy = 1;
        samp_param.n_total_sample_size = 2;
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        std::vector<struct_arg_c> algo_gen_vec(1);
        auto &algo_gen = algo_gen_vec.at(0);
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        indiv_vec.ini(lat, samp_param);
        algo_gen.ini(indiv_vec, samp_param, chr_index);

        //Bidouillage ON
        auto seg0 = algo_gen.Pool_seg.get_new_obj();
        *seg0 = seg(1, 2, 0, nullptr, nullptr);
        indiv_vec[0]->update_homolog_seg(chr_index, seg0);

        auto seg1 = algo_gen.Pool_seg.get_new_obj();
        *seg1 = seg(2, 3, 1, nullptr, nullptr);
        indiv_vec[1]->update_homolog_seg(chr_index, seg1);

        algo_gen.S_intersection_count.at(1) = 2;
        algo_gen.S_intersection_count.emplace(2, 2);
        algo_gen.S_intersection_count.emplace(3, -1);
        //Bidouillage OFF

        auto new_cell = ditrib_and_coa_lign_btw_2_chr(coa_table_vec, indiv_vec, algo_gen_vec, gen, indiv_vec[0], indiv_vec[1]);

        REQUIRE(new_cell == indiv_vec[0]);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 0) == seg0);

        REQUIRE(coa_table.size() == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 1);
        REQUIRE(coa_table_c::get_gen(coa_table[0]) == 1);

        REQUIRE(algo_gen.S_intersection_count.size() == 3);
        REQUIRE(algo_gen.S_intersection_count.at(0) == 0);
        REQUIRE(algo_gen.S_intersection_count.at(1) == 2);
        REQUIRE(algo_gen.S_intersection_count.at(3) == -1);
    }

    SECTION("ditrib_and_coa_lign_btw_2_chr_after_partial_mrca")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        //Two independante coa one with mrca and the other with no real coa    auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 0},
                                       {0, 0}};

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(1);

        int gen = 1;
        samp_param.Ploidy = 1;
        samp_param.n_total_sample_size = 3;
        samp_param.Sequence_length = 2;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        std::vector<struct_arg_c> algo_gen_vec(1);
        auto &algo_gen = algo_gen_vec.at(0);
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        indiv_vec.ini(lat, samp_param);
        algo_gen.ini(indiv_vec, samp_param, chr_index);

        //Bidouillage ON

        indiv_vec[0]->get_chr(chr_index, 0)->R_right_brkpt = 1;
        indiv_vec[1]->get_chr(chr_index, 0)->R_right_brkpt = 1;

        algo_gen.S_intersection_count.at(0) = 2;
        algo_gen.S_intersection_count.emplace(1, 3);

        //Bidouillage OFF

        auto new_cell = ditrib_and_coa_lign_btw_2_chr(coa_table_vec, indiv_vec, algo_gen_vec, gen, indiv_vec[0], indiv_vec[1]);

        new_cell = ditrib_and_coa_lign_btw_2_chr(coa_table_vec, indiv_vec, algo_gen_vec, gen, new_cell, indiv_vec[2]);

        REQUIRE(new_cell == indiv_vec[2]);

        REQUIRE(coa_table.size() == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 3);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 1);
        REQUIRE(coa_table_c::get_gen(coa_table[0]) == 1);

        REQUIRE(algo_gen.S_intersection_count.size() == 3);
        REQUIRE(algo_gen.S_intersection_count.at(0) == 0);
        REQUIRE(algo_gen.S_intersection_count.at(1) == 3);
        REQUIRE(algo_gen.S_intersection_count.at(2) == -1);
    }

    SECTION("ditrib_and_coa_lign_btw_2_chr_multi_coa")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 0},
                                       {0, 0},
                                       {0, 0}};

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(1);

        int gen = 1;
        samp_param.Ploidy = 1;
        samp_param.n_total_sample_size = 4;
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        std::vector<struct_arg_c> algo_gen_vec(1);
        auto &algo_gen = algo_gen_vec.at(0);
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        indiv_vec.ini(lat, samp_param);
        algo_gen.ini(indiv_vec, samp_param, chr_index);

        auto new_cell = ditrib_and_coa_lign_btw_2_chr(coa_table_vec, indiv_vec, algo_gen_vec, gen, indiv_vec[0], indiv_vec[1]);

        new_cell = ditrib_and_coa_lign_btw_2_chr(coa_table_vec, indiv_vec, algo_gen_vec, gen, new_cell, indiv_vec[1]);

        REQUIRE(new_cell == indiv_vec[0]);
        REQUIRE(coa_table.size() == 2);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 4);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 1);
        REQUIRE(coa_table_c::get_gen(coa_table[0]) == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[1]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[1]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[1]) == 5);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[0] == 4);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[1] == 3);
        REQUIRE(coa_table_c::get_gen(coa_table[1]) == 1);

        REQUIRE(algo_gen.S_intersection_count.size() == 2);
        REQUIRE(algo_gen.S_intersection_count.at(0) == 2);
        REQUIRE(algo_gen.S_intersection_count.at(1) == -1);
    }

    SECTION("simulator_hudson_approximate_time")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0}, {0, 0}};
        samp_param.n_total_sample_size = 2;
        samp_param.Chr_nbr = 1;
        samp_param.Sequence_length = 2;

        auto &recomb_param = singleton_c<recomb_param_c>::instance();
        recomb_param.Scaled_recomb_rate_rho = 1;

        std::vector<coa_table_c> coa_table_vec(samp_param.Chr_nbr);

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(5);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;

        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        auto next_max_node_nbr = ARG_simulator_continuous_time(lat, coa_table_vec, samp_param, recomb_param, rand_gen);
        auto coa_table = coa_table_vec[0];

        REQUIRE(next_max_node_nbr[0] == 4);
        REQUIRE(coa_table.size() == 2);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 0);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[1]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[1]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[1]) == 3);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[1] == 1);
    }

    SECTION("crossing_over")
    {
        int chr_index = 0;
        int chr_nbr = 1;
        std::vector<struct_arg_c> algo_gen_vec(1);
        auto &algo_gen = algo_gen_vec.at(0);
        indiv_stock_c indiv_vec(chr_nbr);

        seg *seg_chr_1 = algo_gen.Pool_seg.get_new_obj();
        seg *seg_chr_2 = algo_gen.Pool_seg.get_new_obj();

        *seg_chr_1 = seg(0, 5, 0, nullptr, nullptr);
        *seg_chr_2 = seg(5, 10, 0, nullptr, nullptr);

        algo_gen.L_cumul_nbr_brkpt_per_seg.set_value_at_index(algo_gen.Pool_seg.get_index_from_obj(seg_chr_1), 5);
        algo_gen.L_cumul_nbr_brkpt_per_seg.set_value_at_index(algo_gen.Pool_seg.get_index_from_obj(seg_chr_2), 5);

        algo_gen.S_intersection_count.emplace(0, 2);
        algo_gen.S_intersection_count.emplace(10, -1);

        algo_gen.T_time = 0;

        seg *y_choosen_seg = seg_chr_1;
        int k_choose_brkpt = 4;
        seg *homolog_seg = seg_chr_2;

        seg *first_y_choosen_seg = y_choosen_seg;
        indiv_c indiv(chr_nbr, nullptr);
        indiv.update_indiv(chr_index, seg_chr_1, seg_chr_2);

        crossing_over(algo_gen, y_choosen_seg, k_choose_brkpt, first_y_choosen_seg, homolog_seg, &indiv, chr_index);

        REQUIRE(indiv.get_chr(chr_index, 0)->L_left_brkpt == 0);
        REQUIRE(indiv.get_chr(chr_index, 0)->R_right_brkpt == 4);
        REQUIRE(indiv.get_chr(chr_index, 0)->U_ident_node == 0);
        REQUIRE(indiv.get_chr(chr_index, 0)->Previous_seg == nullptr);

        REQUIRE(indiv.get_chr(chr_index, 0)->Next_seg->L_left_brkpt == 5);
        REQUIRE(indiv.get_chr(chr_index, 0)->Next_seg->R_right_brkpt == 10);
        REQUIRE(indiv.get_chr(chr_index, 0)->Next_seg->U_ident_node == 0);
        REQUIRE(indiv.get_chr(chr_index, 0)->Next_seg->Next_seg == nullptr);

        REQUIRE(indiv.get_chr(chr_index, 1)->L_left_brkpt == 4);
        REQUIRE(indiv.get_chr(chr_index, 1)->R_right_brkpt == 5);
        REQUIRE(indiv.get_chr(chr_index, 1)->U_ident_node == 0);
        REQUIRE(indiv.get_chr(chr_index, 1)->Previous_seg == nullptr);
        REQUIRE(indiv.get_chr(chr_index, 1)->Next_seg == nullptr);

        REQUIRE(algo_gen.S_intersection_count.size() == 2);
        REQUIRE(algo_gen.S_intersection_count.at(0) == 2);
        REQUIRE(algo_gen.S_intersection_count.at(10) == -1);
    }

    SECTION("crossing_over_inter_seg")
    {
        int chr_index = 0;
        int chr_nbr = 1;
        std::vector<struct_arg_c> algo_gen_vec(1);
        auto &algo_gen = algo_gen_vec.at(0);
        indiv_stock_c indiv_vec(chr_nbr);

        seg *seg_chr_1 = algo_gen.Pool_seg.get_new_obj();
        seg *seg_chr_2 = algo_gen.Pool_seg.get_new_obj();
        seg *seg2_chr_2 = algo_gen.Pool_seg.get_new_obj();

        *seg_chr_1 = seg(0, 5, 0, nullptr, nullptr);
        *seg_chr_2 = seg(0, 3, 0, nullptr, seg2_chr_2);
        *seg2_chr_2 = seg(5, 10, 0, seg_chr_2, nullptr);

        algo_gen.L_cumul_nbr_brkpt_per_seg.set_value_at_index(algo_gen.Pool_seg.get_index_from_obj(seg_chr_1), 5);
        algo_gen.L_cumul_nbr_brkpt_per_seg.set_value_at_index(algo_gen.Pool_seg.get_index_from_obj(seg_chr_2), 3);
        algo_gen.L_cumul_nbr_brkpt_per_seg.set_value_at_index(algo_gen.Pool_seg.get_index_from_obj(seg2_chr_2), 7);

        //FALSE S_intersection_count but not use so osef
        algo_gen.S_intersection_count.emplace(0, 2);
        algo_gen.S_intersection_count.emplace(10, -1);

        algo_gen.T_time = 0;

        seg *y_choosen_seg = seg_chr_1;
        int k_choose_brkpt = 4;
        seg *homolog_seg = seg_chr_2;

        seg *first_y_choosen_seg = y_choosen_seg;
        indiv_c indiv(chr_nbr, nullptr);
        indiv.update_indiv(chr_index, seg_chr_1, seg_chr_2);

        crossing_over(algo_gen, y_choosen_seg, k_choose_brkpt, first_y_choosen_seg, homolog_seg, &indiv, chr_index);

        REQUIRE(indiv.get_chr(chr_index, 0)->L_left_brkpt == 0);
        REQUIRE(indiv.get_chr(chr_index, 0)->R_right_brkpt == 4);
        REQUIRE(indiv.get_chr(chr_index, 0)->U_ident_node == 0);
        REQUIRE(indiv.get_chr(chr_index, 0)->Previous_seg == nullptr);

        REQUIRE(indiv.get_chr(chr_index, 0)->Next_seg->L_left_brkpt == 5);
        REQUIRE(indiv.get_chr(chr_index, 0)->Next_seg->R_right_brkpt == 10);
        REQUIRE(indiv.get_chr(chr_index, 0)->Next_seg->U_ident_node == 0);
        REQUIRE(indiv.get_chr(chr_index, 0)->Next_seg->Next_seg == nullptr);

        REQUIRE(indiv.get_chr(chr_index, 1)->L_left_brkpt == 0);
        REQUIRE(indiv.get_chr(chr_index, 1)->R_right_brkpt == 3);
        REQUIRE(indiv.get_chr(chr_index, 1)->U_ident_node == 0);
        REQUIRE(indiv.get_chr(chr_index, 1)->Previous_seg == nullptr);

        REQUIRE(indiv.get_chr(chr_index, 1)->Next_seg->L_left_brkpt == 4);
        REQUIRE(indiv.get_chr(chr_index, 1)->Next_seg->R_right_brkpt == 5);
        REQUIRE(indiv.get_chr(chr_index, 1)->Next_seg->U_ident_node == 0);
        REQUIRE(indiv.get_chr(chr_index, 1)->Next_seg->Next_seg == nullptr);

        REQUIRE(algo_gen.S_intersection_count.size() == 2);
        REQUIRE(algo_gen.S_intersection_count.at(0) == 2);
        REQUIRE(algo_gen.S_intersection_count.at(10) == -1);
    }

    SECTION("haploid_first_recomb_gen_by_gen")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0}};

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        samp_param.Ploidy = 1;
        rand_gen.put_seed(3);

        samp_param.n_total_sample_size = 1;
        samp_param.Sequence_length = 2;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;
        double unscaled_recomb_rate = 1;

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        std::vector<struct_arg_c> algo_gen_vec(1);
        auto &algo_gen = algo_gen_vec.at(0);
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);

        indiv_vec.ini(lat, samp_param);
        algo_gen.ini(indiv_vec, samp_param, chr_index);

        recomb_gen_by_gen(algo_gen, unscaled_recomb_rate, chr_index, rand_gen);

        REQUIRE(indiv_vec.size() == 1);
        REQUIRE(algo_gen.L_cumul_nbr_brkpt_per_seg.get_tot_cumul_freq() == 0);

        REQUIRE(indiv_vec[0]->get_chr(chr_index, 0) == algo_gen.Pool_seg.get_obj_from_index(0));
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 0)->get_indiv() == indiv_vec[0]);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 0)->L_left_brkpt == 0);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 0)->R_right_brkpt == 1);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 0)->U_ident_node == 0);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 0)->Previous_seg == nullptr);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 0)->Next_seg == nullptr);

        REQUIRE(indiv_vec[0]->get_chr(chr_index, 1) == algo_gen.Pool_seg.get_obj_from_index(1));
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 1)->get_indiv() == indiv_vec[0]);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 1)->L_left_brkpt == 1);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 1)->R_right_brkpt == 2);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 1)->U_ident_node == 0);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 1)->Previous_seg == nullptr);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 1)->Next_seg == nullptr);

        REQUIRE(algo_gen.S_intersection_count.size() == 2);
        REQUIRE(algo_gen.S_intersection_count.at(0) == 1);
        REQUIRE(algo_gen.S_intersection_count.at(2) == -1);
    }

    SECTION("haploid_second_recomb_gen_by_gen")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0}};

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        samp_param.Ploidy = 1;
        rand_gen.put_seed(3);

        samp_param.n_total_sample_size = 1;
        samp_param.Sequence_length = 3;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;
        double unscaled_recomb_rate = 1;

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        std::vector<struct_arg_c> algo_gen_vec(1);
        auto &algo_gen = algo_gen_vec.at(0);
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);

        indiv_vec.ini(lat, samp_param);
        algo_gen.ini(indiv_vec, samp_param, chr_index);

        recomb_gen_by_gen(algo_gen, unscaled_recomb_rate, chr_index, rand_gen);

        REQUIRE(indiv_vec.size() == 1);
        REQUIRE(algo_gen.L_cumul_nbr_brkpt_per_seg.get_tot_cumul_freq() == 2);

        auto chr1_1 = indiv_vec[0]->get_chr(chr_index, 1);
        auto chr1_2 = indiv_vec[0]->get_chr(chr_index, 1)->Next_seg;
        auto chr2 = indiv_vec[0]->get_chr(chr_index, 0);

        REQUIRE(chr1_1 == algo_gen.Pool_seg.get_obj_from_index(0));
        REQUIRE(chr1_1->get_indiv() == indiv_vec[0]);
        REQUIRE(chr1_1->L_left_brkpt == 0);
        REQUIRE(chr1_1->R_right_brkpt == 1);
        REQUIRE(chr1_1->U_ident_node == 0);
        REQUIRE(chr1_1->Previous_seg == nullptr);
        REQUIRE(chr1_1->Next_seg == chr1_2);

        REQUIRE(chr1_2 == algo_gen.Pool_seg.get_obj_from_index(2));
        REQUIRE(chr1_2->get_indiv() == nullptr);
        REQUIRE(chr1_2->L_left_brkpt == 2);
        REQUIRE(chr1_2->R_right_brkpt == 3);
        REQUIRE(chr1_2->U_ident_node == 0);
        REQUIRE(chr1_2->Previous_seg == chr1_1);
        REQUIRE(chr1_2->Next_seg == nullptr);

        REQUIRE(chr2 == algo_gen.Pool_seg.get_obj_from_index(1));
        REQUIRE(chr2->get_indiv() == indiv_vec[0]);
        REQUIRE(chr2->L_left_brkpt == 1);
        REQUIRE(chr2->R_right_brkpt == 2);
        REQUIRE(chr2->U_ident_node == 0);
        REQUIRE(chr2->Previous_seg == nullptr);
        REQUIRE(chr2->Next_seg == nullptr);

        REQUIRE(algo_gen.S_intersection_count.size() == 2);
        REQUIRE(algo_gen.S_intersection_count.at(0) == 1);
        REQUIRE(algo_gen.S_intersection_count.at(3) == -1);
    }

    SECTION("haploid_coa_gen_by_gen_one_pop_haploid")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 0},
                                       {0, 0},
                                       {0, 0}};

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        samp_param.Ploidy = 1;
        rand_gen.put_seed(2);

        int gen = 1;
        samp_param.n_total_sample_size = 4;
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        std::vector<struct_arg_c> algo_gen_vec(1);
        auto &algo_gen = algo_gen_vec.at(0);
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        indiv_vec.ini(lat, samp_param);
        algo_gen.ini(indiv_vec, samp_param, chr_index);

        coa_gen_by_gen(indiv_vec, algo_gen_vec, lat, coa_table_vec, gen, samp_param, rand_gen);

        REQUIRE(coa_table.size() == 1);
        //Maternal chr
        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 6);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 3);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[2] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[3] == 1);
        REQUIRE(coa_table_c::get_time_t(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_gen(coa_table[0]) == 1);
    }

    SECTION("haploid_coa_gen_by_gen_two_pop_haploid")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 0},
                                       {1, 1},
                                       {1, 1}};

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        samp_param.Ploidy = 1;
        rand_gen.put_seed(5);

        int gen = 1;
        samp_param.n_total_sample_size = 4;
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {2, 2};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        std::vector<struct_arg_c> algo_gen_vec(1);
        auto &algo_gen = algo_gen_vec.at(0);
        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        indiv_vec.ini(lat, samp_param);
        algo_gen.ini(indiv_vec, samp_param, chr_index);

        coa_gen_by_gen(indiv_vec, algo_gen_vec, lat, coa_table_vec, gen, samp_param, rand_gen);

        REQUIRE(lat.Lattice[0].Indivs_in_pop.size() == 2);
        REQUIRE(lat.Lattice[1].Indivs_in_pop.size() == 0);

        REQUIRE(coa_table.size() == 2);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 4);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 1);
        REQUIRE(coa_table_c::get_gen(coa_table[0]) == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[1]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[1]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[1]) == 5);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[0] == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[1] == 3);
        REQUIRE(coa_table_c::get_gen(coa_table[1]) == 1);
    }

    SECTION("haploid_simulator_gen_by_gen_simple_coa")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0}, {0, 0}};
        samp_param.Ploidy = 1;

        samp_param.n_total_sample_size = 2;

        auto &recomb_param = singleton_c<recomb_param_c>::instance();
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 1;
        recomb_param.Unscaled_recomb_rate = 0;

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 1;
        demo_param.Population_size_N = demo_param.Pop_size_per_node * demo_param.Lattice_size.at(0) * demo_param.Lattice_size.at(1);
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        auto next_max_node_nbr = ARG_simulator_gen_by_gen(lat, rmap, coa_table_vec, rand_gen);

        REQUIRE(next_max_node_nbr.at(0) == 3);
        REQUIRE(coa_table.size() == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 1);
        REQUIRE(coa_table_c::get_gen(coa_table[0]) == 1);
    }

    SECTION("haploid_simulator_gen_by_gen_multi_coa")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0}, {0, 0}, {0, 0}};
        samp_param.Ploidy = 1;

        samp_param.n_total_sample_size = 3;

        auto &recomb_param = singleton_c<recomb_param_c>::instance();
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 1;
        recomb_param.Unscaled_recomb_rate = 0;

        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(1);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 1;
        demo_param.Population_size_N = demo_param.Pop_size_per_node * demo_param.Lattice_size.at(0) * demo_param.Lattice_size.at(1);
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        auto next_max_node_nbr = ARG_simulator_gen_by_gen(lat, rmap, coa_table_vec, rand_gen);

        REQUIRE(next_max_node_nbr.at(0) == 5);
        REQUIRE(coa_table.size() == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 4);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[2] == 1);
        REQUIRE(coa_table_c::get_gen(coa_table[0]) == 1);
    }

    SECTION("haploid_simulator_gen_by_gen_simple_recomb")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0}, {0, 0}};
        samp_param.Ploidy = 1;

        samp_param.n_total_sample_size = 2;

        auto &recomb_param = singleton_c<recomb_param_c>::instance();
        samp_param.Sequence_length = 2;
        samp_param.Chr_nbr = 1;
        recomb_param.Unscaled_recomb_rate = 1;

        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 2;
        demo_param.Population_size_N = demo_param.Pop_size_per_node * demo_param.Lattice_size.at(0) * demo_param.Lattice_size.at(1);
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        auto next_max_node_nbr = ARG_simulator_gen_by_gen(lat, rmap, coa_table_vec, rand_gen);

        REQUIRE(next_max_node_nbr.at(0) == 4);
        REQUIRE(coa_table.size() == 2);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 1);
        REQUIRE(coa_table_c::get_gen(coa_table[0]) == 2);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[1]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[1]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[1]) == 3);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[1] == 1);
        REQUIRE(coa_table_c::get_gen(coa_table[1]) == 3);
    }

    SECTION("haploid_simulator_gen_by_gen_multi_recomb")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0}, {0, 0}};
        samp_param.Ploidy = 1;

        samp_param.n_total_sample_size = 2;

        auto &recomb_param = singleton_c<recomb_param_c>::instance();
        samp_param.Sequence_length = 3;
        samp_param.Chr_nbr = 1;
        recomb_param.Unscaled_recomb_rate = 1;

        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 10;
        demo_param.Population_size_N = demo_param.Pop_size_per_node * demo_param.Lattice_size.at(0) * demo_param.Lattice_size.at(1);
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);
        extend_lattice_c rmap(lat);

        auto next_max_node_nbr = ARG_simulator_gen_by_gen(lat, rmap, coa_table_vec, rand_gen);

        REQUIRE(next_max_node_nbr.at(0) == 5);
        REQUIRE(coa_table.size() == 3);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 0);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[1]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[1]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[1]) == 3);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[1] == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[2]) == 2);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[2]) == 3);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[2]) == 4);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[2]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[2]))[1] == 1);
    }

    SECTION("diploid_simulator_gen_by_gen_simple_recomb")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0}};
        samp_param.Ploidy = 2;

        samp_param.n_total_sample_size = 1;

        auto &recomb_param = singleton_c<recomb_param_c>::instance();
        samp_param.Sequence_length = 2;
        samp_param.Chr_nbr = 1;
        recomb_param.Unscaled_recomb_rate = 1;

        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 1;
        demo_param.Population_size_N = demo_param.Pop_size_per_node * demo_param.Lattice_size.at(0) * demo_param.Lattice_size.at(1);
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        extend_lattice_c rmap(lat);

        auto next_max_node_nbr = ARG_simulator_gen_by_gen(lat, rmap, coa_table_vec, rand_gen);

        REQUIRE(next_max_node_nbr.at(0) == 4);
        REQUIRE(coa_table.size() == 2);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 1);
        REQUIRE(coa_table_c::get_gen(coa_table[0]) == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[1]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[1]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[1]) == 3);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[1] == 1);
        REQUIRE(coa_table_c::get_gen(coa_table[1]) == 1);
    }

    SECTION("diploid_simulator_gen_by_gen_multi_recomb")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0}, {0, 0}};
        samp_param.Ploidy = 2;

        samp_param.n_total_sample_size = 2;

        auto &recomb_param = singleton_c<recomb_param_c>::instance();
        samp_param.Sequence_length = 2;
        samp_param.Chr_nbr = 1;
        recomb_param.Unscaled_recomb_rate = 1;

        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        rand_gen_c rand_gen;
        rand_gen.put_seed(3);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 10;
        demo_param.Population_size_N = demo_param.Pop_size_per_node * demo_param.Lattice_size.at(0) * demo_param.Lattice_size.at(1);
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);
        extend_lattice_c rmap(lat);

        auto next_max_node_nbr = ARG_simulator_gen_by_gen(lat, rmap, coa_table_vec, rand_gen);

        REQUIRE(next_max_node_nbr.at(0) == 10);
        REQUIRE(coa_table.size() == 6);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 4);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 3);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[1]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[1]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[1]) == 5);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[0] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[1] == 3);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[2]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[2]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[2]) == 6);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[2]))[0] == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[2]))[1] == 5);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[3]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[3]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[3]) == 7);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[3]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[3]))[1] == 4);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[4]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[4]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[4]) == 8);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[4]))[0] == 7);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[4]))[1] == 2);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[5]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[5]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[5]) == 9);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[5]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[5]))[1] == 6);
    }

    SECTION("diploid_simulator_gen_by_gen_2_subpop_simple_coa")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 1}};
        samp_param.Ploidy = 2;

        samp_param.n_total_sample_size = 2;

        auto &recomb_param = singleton_c<recomb_param_c>::instance();
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 1;
        recomb_param.Unscaled_recomb_rate = 0;

        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        auto &rand_gen = singleton_c<rand_gen_c>::instance();

        rand_gen.put_seed(3);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 2};
        demo_param.Pop_size_per_node = 1;
        demo_param.Population_size_N = demo_param.Pop_size_per_node * demo_param.Lattice_size.at(0) * demo_param.Lattice_size.at(1);
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.5;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);
        extend_lattice_c rmap(lat);

        auto next_max_node_nbr = ARG_simulator_gen_by_gen(lat, rmap, coa_table_vec, rand_gen);

        REQUIRE(next_max_node_nbr.at(0) == 7);
        REQUIRE(coa_table.size() == 3);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 4);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[1]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[1]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[1]) == 5);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[0] == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[1] == 3);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[2]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[2]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[2]) == 6);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[2]))[0] == 4);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[2]))[1] == 5);
    }

    SECTION("diploid_simulator_gen_by_gen_2_subpop")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 0},
                                       {0, 1},
                                       {0, 1}};
        samp_param.Ploidy = 2;

        samp_param.n_total_sample_size = 4;

        auto &recomb_param = singleton_c<recomb_param_c>::instance();
        samp_param.Sequence_length = 2;
        samp_param.Chr_nbr = 1;
        recomb_param.Unscaled_recomb_rate = 0.5;

        std::vector<coa_table_c> coa_table_vec(1);
        auto &coa_table = coa_table_vec.at(0);

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 2};
        demo_param.Pop_size_per_node = 2;
        demo_param.Population_size_N = demo_param.Pop_size_per_node * demo_param.Lattice_size.at(0) * demo_param.Lattice_size.at(1);
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.5;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);
        extend_lattice_c rmap(lat);

        auto next_max_node_nbr = ARG_simulator_gen_by_gen(lat, rmap, coa_table_vec, rand_gen);

        REQUIRE(next_max_node_nbr.at(0) == 20);
        REQUIRE(coa_table.size() == 11);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[0]) == 9);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0])).size() == 3);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[0] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[1] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[0]))[2] == 2);
        REQUIRE(coa_table_c::get_gen(coa_table[0]) == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[1]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[1]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[1]) == 12);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[0] == 3);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[1]))[1] == 7);
        REQUIRE(coa_table_c::get_gen(coa_table[1]) == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[4]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[4]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[4]) == 13);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[4])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[4]))[0] == 4);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[4]))[1] == 6);
        REQUIRE(coa_table_c::get_gen(coa_table[4]) == 3);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[7]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[7]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[7]) == 16);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[7])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[7]))[0] == 9);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[7]))[1] == 12);
        REQUIRE(coa_table_c::get_gen(coa_table[7]) == 7);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[8]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[8]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[8]) == 17);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[8])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[8]))[0] == 16);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[8]))[1] == 5);
        REQUIRE(coa_table_c::get_gen(coa_table[8]) == 10);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[9]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[9]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[9]) == 18);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[9])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[9]))[0] == 17);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[9]))[1] == 13);
        REQUIRE(coa_table_c::get_gen(coa_table[9]) == 16);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[2]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[2]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[2]) == 11);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[2])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[2]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[2]))[1] == 3);
        REQUIRE(coa_table_c::get_gen(coa_table[2]) == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[3]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[3]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[3]) == 12);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[3])).size() == 3);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[3]))[0] == 7);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[3]))[1] == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[3]))[2] == 1);
        REQUIRE(coa_table_c::get_gen(coa_table[3]) == 1);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[5]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[5]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[5]) == 14);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[5])).size() == 3);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[5]))[0] == 12);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[5]))[1] == 4);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[5]))[2] == 6);
        REQUIRE(coa_table_c::get_gen(coa_table[5]) == 3);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[6]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[6]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[6]) == 15);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[6])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[6]))[0] == 11);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[6]))[1] == 5);
        REQUIRE(coa_table_c::get_gen(coa_table[6]) == 4);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table[10]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table[10]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table[10]) == 19);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[10])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[10]))[0] == 14);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table[10]))[1] == 15);
        REQUIRE(coa_table_c::get_gen(coa_table[10]) == 22);
    }

    SECTION("diploid_simulator_gen_by_gen_2_subpop_2_chr")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 0},
                                       {0, 1},
                                       {0, 1}};
        samp_param.Ploidy = 1;

        samp_param.n_total_sample_size = 4;

        auto &recomb_param = singleton_c<recomb_param_c>::instance();
        samp_param.Sequence_length = 2;
        samp_param.Chr_nbr = 2;
        recomb_param.Unscaled_recomb_rate = 0.5;

        std::vector<coa_table_c> coa_table_vec(samp_param.Chr_nbr);

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 2};
        demo_param.Pop_size_per_node = 10;
        demo_param.Population_size_N = demo_param.Pop_size_per_node * demo_param.Lattice_size.at(0) * demo_param.Lattice_size.at(1);
        demo_param.Disp_dist_max = {1, 1};
        demo_param.Proba_migr = 0.5;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);
        extend_lattice_c rmap(lat);

        auto next_max_node_nbr = ARG_simulator_gen_by_gen(lat, rmap, coa_table_vec, rand_gen);

        auto &coa_table1 = coa_table_vec.at(0);
        REQUIRE(next_max_node_nbr.at(0) == 10);
        REQUIRE(coa_table1.size() == 5);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table1[1]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table1[1]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table1[1]) == 6);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table1[1])).size() == 3);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table1[1]))[0] == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table1[1]))[1] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table1[1]))[2] == 1);
        REQUIRE(coa_table_c::get_gen(coa_table1[1]) == 4);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table1[4]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table1[4]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table1[4]) == 9);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table1[4])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table1[4]))[0] == 6);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table1[4]))[1] == 3);
        REQUIRE(coa_table_c::get_gen(coa_table1[4]) == 36);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table1[0]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table1[0]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table1[0]) == 4);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table1[0])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table1[0]))[0] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table1[0]))[1] == 0);
        REQUIRE(coa_table_c::get_gen(coa_table1[0]) == 2);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table1[2]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table1[2]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table1[2]) == 7);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table1[2])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table1[2]))[0] == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table1[2]))[1] == 4);
        REQUIRE(coa_table_c::get_gen(coa_table1[2]) == 8);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table1[3]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table1[3]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table1[3]) == 8);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table1[3])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table1[3]))[0] == 7);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table1[3]))[1] == 3);
        REQUIRE(coa_table_c::get_gen(coa_table1[3]) == 30);

        auto &coa_table2 = coa_table_vec.at(1);
        REQUIRE(next_max_node_nbr.at(1) == 10);
        REQUIRE(coa_table2.size() == 6);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table2[0]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table2[0]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table2[0]) == 4);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[0])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[0]))[0] == 3);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[0]))[1] == 0);
        REQUIRE(coa_table_c::get_gen(coa_table2[0]) == 4);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table2[4]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table2[4]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table2[4]) == 8);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[4])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[4]))[0] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[4]))[1] == 4);
        REQUIRE(coa_table_c::get_gen(coa_table2[4]) == 17);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table2[5]) == 0);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table2[5]) == 1);
        REQUIRE(coa_table_c::get_num_node_u(coa_table2[5]) == 9);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[5])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[5]))[0] == 8);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[5]))[1] == 2);
        REQUIRE(coa_table_c::get_gen(coa_table2[5]) == 18);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table2[1]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table2[1]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table2[1]) == 5);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[1])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[1]))[0] == 0);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[1]))[1] == 3);
        REQUIRE(coa_table_c::get_gen(coa_table2[1]) == 7);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table2[2]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table2[2]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table2[2]) == 6);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[2])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[2]))[0] == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[2]))[1] == 5);
        REQUIRE(coa_table_c::get_gen(coa_table2[2]) == 9);

        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table2[3]) == 1);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table2[3]) == 2);
        REQUIRE(coa_table_c::get_num_node_u(coa_table2[3]) == 7);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[3])).size() == 2);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[3]))[0] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table2[3]))[1] == 6);
        REQUIRE(coa_table_c::get_gen(coa_table2[3]) == 14);
    }
}

TEST_CASE("segregate_indiv_simulator_test")
{
    SECTION("segregate_indiv_mono_locus")
    {
        int chr_nbr = 1;
        int chr_index = 0;
        indiv_stock_c indiv_vec(chr_nbr);
        obj_mem_manag_c<seg> pool_seg;

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        seg *seg0 = pool_seg.get_new_obj();
        seg *seg1 = pool_seg.get_new_obj();
        seg *seg2 = pool_seg.get_new_obj();

        indiv_vec.new_indiv(lat.node({0, 0}));
        indiv_vec.back()->update_indiv(chr_index, seg0, seg1);
        lat.add_indiv(indiv_vec.back(), {0, 0});
        seg0->set_indiv(indiv_vec[0]);
        seg1->set_indiv(indiv_vec[0]);

        indiv_vec.new_indiv(lat.node({0, 0}));
        indiv_vec.back()->update_indiv(chr_index, seg2, nullptr);
        lat.add_indiv(indiv_vec.back(), {0, 0});
        seg2->set_indiv(indiv_vec[1]);

        segregate_chr_btw_indiv(lat, indiv_vec, chr_nbr, rand_gen);

        REQUIRE(indiv_vec.size() == 3);

        REQUIRE(indiv_vec[0]->get_chr(chr_index, 0) == seg0);
        REQUIRE(seg0->get_indiv() == indiv_vec[0]);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 1) == nullptr);

        REQUIRE(indiv_vec[1]->get_chr(chr_index, 0) == seg2);
        REQUIRE(seg2->get_indiv() == indiv_vec[1]);
        REQUIRE(indiv_vec[1]->get_chr(chr_index, 1) == nullptr);

        REQUIRE(indiv_vec[2]->get_chr(chr_index, 0) == seg1);
        REQUIRE(seg1->get_indiv() == indiv_vec[2]);
        REQUIRE(indiv_vec[2]->get_chr(chr_index, 1) == nullptr);

        REQUIRE(lat.Lattice[0].Indivs_in_pop.size() == 3);
    }

    SECTION("segregate_indiv_multi_locus")
    {
        int chr_nbr = 3;
        indiv_stock_c indiv_vec(chr_nbr);
        obj_mem_manag_c<seg> pool_seg;

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(55);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        std::vector<seg *> seg_result(chr_nbr * 2, nullptr);
        indiv_vec.new_indiv(lat.node({0, 0}));
        int seg_index = 0;
        for (int chr_index = 0; chr_index < chr_nbr; ++chr_index)
        {
            seg_result[seg_index] = pool_seg.get_new_obj();
            seg_result[seg_index + 1] = pool_seg.get_new_obj();
            indiv_vec.back()->update_indiv(chr_index, seg_result[seg_index], seg_result[seg_index + 1]);
            seg_result[seg_index]->set_indiv(indiv_vec[0]);
            seg_result[seg_index + 1]->set_indiv(indiv_vec[0]);
            seg_index += 2;
        }
        lat.add_indiv(indiv_vec.back(), {0, 0});

        std::vector<seg *> seg_result1(chr_nbr, nullptr);
        indiv_vec.new_indiv(lat.node({0, 0}));
        for (int chr_index = 0; chr_index < chr_nbr; ++chr_index)
        {
            seg_result1[chr_index] = pool_seg.get_new_obj();
            indiv_vec.back()->update_indiv(chr_index, seg_result1[chr_index], nullptr);
            seg_result1[chr_index]->set_indiv(indiv_vec[1]);
        }
        lat.add_indiv(indiv_vec.back(), {0, 0});
        segregate_chr_btw_indiv(lat, indiv_vec, chr_nbr, rand_gen);

        REQUIRE(indiv_vec.size() == 4);

        REQUIRE(indiv_vec[0]->get_chr(0, 0) == seg_result[0]);
        REQUIRE(seg_result[0]->get_indiv() == indiv_vec[0]);
        REQUIRE(indiv_vec[0]->get_chr(0, 1) == nullptr);

        REQUIRE(indiv_vec[0]->get_chr(1, 0) == seg_result[3]);
        REQUIRE(seg_result[3]->get_indiv() == indiv_vec[0]);
        REQUIRE(indiv_vec[0]->get_chr(1, 1) == nullptr);

        REQUIRE(indiv_vec[0]->get_chr(2, 0) == seg_result[5]);
        REQUIRE(seg_result[5]->get_indiv() == indiv_vec[0]);
        REQUIRE(indiv_vec[0]->get_chr(2, 1) == nullptr);

        REQUIRE(indiv_vec[2]->get_chr(0, 0) == seg_result[1]);
        REQUIRE(seg_result[1]->get_indiv() == indiv_vec[2]);
        REQUIRE(indiv_vec[2]->get_chr(0, 1) == nullptr);

        REQUIRE(indiv_vec[2]->get_chr(1, 0) == seg_result[2]);
        REQUIRE(seg_result[2]->get_indiv() == indiv_vec[2]);
        REQUIRE(indiv_vec[2]->get_chr(1, 1) == nullptr);

        REQUIRE(indiv_vec[2]->get_chr(2, 0) == seg_result[4]);
        REQUIRE(seg_result[4]->get_indiv() == indiv_vec[2]);
        REQUIRE(indiv_vec[2]->get_chr(2, 1) == nullptr);

        REQUIRE(indiv_vec[1]->get_chr(0, 0) == seg_result1[0]);
        REQUIRE(seg_result1[0]->get_indiv() == indiv_vec[1]);
        REQUIRE(indiv_vec[1]->get_chr(0, 1) == nullptr);

        REQUIRE(indiv_vec[1]->get_chr(1, 0) == nullptr);
        REQUIRE(indiv_vec[1]->get_chr(1, 1) == nullptr);

        REQUIRE(indiv_vec[1]->get_chr(2, 0) == seg_result1[2]);
        REQUIRE(seg_result1[2]->get_indiv() == indiv_vec[1]);
        REQUIRE(indiv_vec[1]->get_chr(2, 1) == nullptr);

        REQUIRE(indiv_vec[3]->get_chr(0, 0) == nullptr);
        REQUIRE(indiv_vec[3]->get_chr(0, 1) == nullptr);

        REQUIRE(indiv_vec[3]->get_chr(1, 0) == seg_result1[1]);
        REQUIRE(seg_result1[1]->get_indiv() == indiv_vec[3]);
        REQUIRE(indiv_vec[3]->get_chr(1, 1) == nullptr);

        REQUIRE(indiv_vec[3]->get_chr(2, 0) == nullptr);
        REQUIRE(indiv_vec[3]->get_chr(2, 1) == nullptr);

        REQUIRE(lat.Lattice[0].Indivs_in_pop.size() == 4);
    }

    SECTION("segregate_indiv_multi_locus_with_empty_first_locus")
    {
        int chr_nbr = 2;
        indiv_stock_c indiv_vec(chr_nbr);
        obj_mem_manag_c<seg> pool_seg;

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(55);

        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Edge_effects = edge_effect_enum::absorbing;
        demo_param.Lattice_size = {1, 1};
        demo_param.Pop_size_per_node = 1;
        demo_param.Disp_dist_max = {0, 0};
        demo_param.Proba_migr = 0.2;
        double g_geo_param = 0.2;
        demo_param.Disp_func[0] = fwd_disp_distrib_c::geometric(g_geo_param);
        demo_param.Disp_func[1] = fwd_disp_distrib_c::geometric(g_geo_param);

        fwd_disp_distrib_c fwd_distrib(demo_param);
        lattice_c lat(rand_gen, fwd_distrib, demo_param);

        indiv_vec.new_indiv(lat.node({0, 0}));
        indiv_vec.back()->update_indiv(0, nullptr, nullptr);
        seg *seg1 = pool_seg.get_new_obj();
        seg *seg2 = pool_seg.get_new_obj();
        indiv_vec.back()->update_indiv(1, seg1, seg2);
        seg1->set_indiv(indiv_vec[0]);
        seg2->set_indiv(indiv_vec[0]);
        lat.add_indiv(indiv_vec.back(), {0, 0});

        indiv_vec.new_indiv(lat.node({0, 0}));
        seg *seg3 = pool_seg.get_new_obj();
        indiv_vec.back()->update_indiv(1, seg3, nullptr);
        seg3->set_indiv(indiv_vec[1]);
        lat.add_indiv(indiv_vec.back(), {0, 0});

        segregate_chr_btw_indiv(lat, indiv_vec, chr_nbr, rand_gen);

        REQUIRE(indiv_vec.size() == 3);

        REQUIRE(indiv_vec[0]->get_chr(0, 0) == nullptr);
        REQUIRE(indiv_vec[0]->get_chr(0, 1) == nullptr);

        REQUIRE(indiv_vec[0]->get_chr(1, 0) == seg1);
        REQUIRE(seg1->get_indiv() == indiv_vec[0]);
        REQUIRE(indiv_vec[0]->get_chr(1, 1) == nullptr);

        REQUIRE(indiv_vec[2]->get_chr(0, 0) == nullptr);
        REQUIRE(indiv_vec[2]->get_chr(0, 1) == nullptr);

        REQUIRE(indiv_vec[2]->get_chr(1, 0) == seg2);
        REQUIRE(seg2->get_indiv() == indiv_vec[2]);
        REQUIRE(indiv_vec[2]->get_chr(1, 1) == nullptr);

        REQUIRE(indiv_vec[1]->get_chr(0, 0) == nullptr);
        REQUIRE(indiv_vec[1]->get_chr(0, 1) == nullptr);

        REQUIRE(indiv_vec[1]->get_chr(1, 0) == seg3);
        REQUIRE(seg3->get_indiv() == indiv_vec[1]);
        REQUIRE(indiv_vec[1]->get_chr(1, 1) == nullptr);

        REQUIRE(lat.Lattice[0].Indivs_in_pop.size() == 3);
    }
}