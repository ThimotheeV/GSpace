#define CATCH_CONFIG_MAIN AlgoRecombinationTest
#include "catch.hpp"
#include <iostream>

#include "kelleher_algo.hpp"

TEST_CASE()
{
    REQUIRE(1);
}

TEST_CASE("algo_r_test")
{
    SECTION("choose_brkpt_one_seg")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0}};
        samp_param.Ploidy = 1;
        samp_param.n_total_sample_size = 1;
        samp_param.Sequence_length = 10;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

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
        main.ini(indiv_vec, samp_param, chr_index);

        alg_recomb_c recomb;

        std::uniform_int_distribution<long int> uniform_distribution(1, main.L_cumul_nbr_brkpt_per_seg.get_tot_cumul_freq());
        long int h_choosen_brkpt = uniform_distribution(rand_gen.Seed_gen);

        int choosen_seg_index = main.L_cumul_nbr_brkpt_per_seg.seg_index(h_choosen_brkpt);
        auto y_choosen_seg = main.Pool_seg.get_obj_from_index(choosen_seg_index);

        auto k_choose_brkpt = y_choosen_seg->R_right_brkpt - main.L_cumul_nbr_brkpt_per_seg.get_cumul_freq_at_index(choosen_seg_index) + h_choosen_brkpt - 1;

        bool result = recomb.ini(y_choosen_seg, k_choose_brkpt);

        REQUIRE(result == 1);
        REQUIRE(main.Pool_seg.get_index_from_obj(recomb.Y_choosen_seg) == 0);
        REQUIRE(!(recomb.X_prev_seg));
        REQUIRE(!(recomb.Z_new_seg));
    }

    SECTION("choose_brkpt_two_segs")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                    {0, 0}};
        samp_param.Ploidy = 1;
        samp_param.n_total_sample_size = 2;
        samp_param.Sequence_length = 10;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

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
        main.ini(indiv_vec, samp_param, chr_index);
        alg_recomb_c recomb;

        std::uniform_int_distribution<long int> uniform_distribution(1, main.L_cumul_nbr_brkpt_per_seg.get_tot_cumul_freq());
        long int h_choosen_brkpt = uniform_distribution(rand_gen.Seed_gen);

        int choosen_seg_index = main.L_cumul_nbr_brkpt_per_seg.seg_index(h_choosen_brkpt);
        auto y_choosen_seg = main.Pool_seg.get_obj_from_index(choosen_seg_index);

        auto k_choose_brkpt = y_choosen_seg->R_right_brkpt - main.L_cumul_nbr_brkpt_per_seg.get_cumul_freq_at_index(choosen_seg_index) + h_choosen_brkpt - 1;

        bool result = recomb.ini(y_choosen_seg, k_choose_brkpt);

        REQUIRE(result == 1);
        REQUIRE(!(recomb.X_prev_seg));
        REQUIRE(!(recomb.Z_new_seg));
    }

    SECTION("choose_brkpt_two_sub_segs")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                    {0, 0}};
        samp_param.Ploidy = 1;
        samp_param.n_total_sample_size = 2;
        samp_param.Sequence_length = 10;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

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
        main.ini(indiv_vec, samp_param, chr_index);

        seg seg0;
        seg seg1;

        seg0 = seg(0, 1, 0, nullptr, &seg1);
        seg1 = seg(6, 7, 1, &seg0, nullptr);

        (indiv_vec[0])->update_indiv(chr_index, &seg0, nullptr);
        main.L_cumul_nbr_brkpt_per_seg.up_at_index_with_value(0, -8);

        alg_recomb_c recomb;

        std::uniform_int_distribution<long int> uniform_distribution(1, main.L_cumul_nbr_brkpt_per_seg.get_tot_cumul_freq());
        long int h_choosen_brkpt = uniform_distribution(rand_gen.Seed_gen);

        int choosen_seg_index = main.L_cumul_nbr_brkpt_per_seg.seg_index(h_choosen_brkpt);
        auto y_choosen_seg = main.Pool_seg.get_obj_from_index(choosen_seg_index);

        auto k_choose_brkpt = y_choosen_seg->R_right_brkpt - main.L_cumul_nbr_brkpt_per_seg.get_cumul_freq_at_index(choosen_seg_index) + h_choosen_brkpt - 1;

        bool result = recomb.ini(y_choosen_seg, k_choose_brkpt);

        REQUIRE(result == 1);
        REQUIRE(main.Pool_seg.get_index_from_obj(recomb.Y_choosen_seg) == 1);
        REQUIRE(!(recomb.X_prev_seg));
        REQUIRE(!(recomb.Z_new_seg));
    }

    SECTION("break_between_seg_two_sub_seg")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                    {0, 0}};
        samp_param.Ploidy = 1;
        samp_param.n_total_sample_size = 2;
        samp_param.Sequence_length = 10;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

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
        main.ini(indiv_vec, samp_param, chr_index);

        indiv_vec[0]->get_chr(chr_index, 0)->Next_seg = indiv_vec[1]->get_chr(chr_index, 0);
        main.L_cumul_nbr_brkpt_per_seg.up_at_index_with_value(0, -8);
        indiv_vec[1]->get_chr(chr_index, 0)->Previous_seg = indiv_vec[0]->get_chr(chr_index, 0);

        alg_recomb_c recomb;

        std::uniform_int_distribution<long int> uniform_distribution(1, main.L_cumul_nbr_brkpt_per_seg.get_tot_cumul_freq());
        long int h_choosen_brkpt = uniform_distribution(rand_gen.Seed_gen);

        int choosen_seg_index = main.L_cumul_nbr_brkpt_per_seg.seg_index(h_choosen_brkpt);
        auto y_choosen_seg = main.Pool_seg.get_obj_from_index(choosen_seg_index);

        auto k_choose_brkpt = y_choosen_seg->R_right_brkpt - main.L_cumul_nbr_brkpt_per_seg.get_cumul_freq_at_index(choosen_seg_index) + h_choosen_brkpt - 1;

        recomb.ini(y_choosen_seg, k_choose_brkpt);
        recomb.break_between_segment();

        REQUIRE(main.Pool_seg.get_index_from_obj(recomb.Y_choosen_seg) == 1);
        REQUIRE(recomb.Z_new_seg == recomb.Y_choosen_seg);
        REQUIRE(!(indiv_vec[0]->get_chr(chr_index, 0)->Next_seg));
        REQUIRE(!(indiv_vec[1]->get_chr(chr_index, 0)->Previous_seg));
    }

    SECTION("break_within_seg_one_seg")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0}};
        samp_param.Ploidy = 1;
        samp_param.n_total_sample_size = 1;
        samp_param.Sequence_length = 10;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

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
        main.ini(indiv_vec, samp_param, chr_index);
        alg_recomb_c recomb;

        std::uniform_int_distribution<long int> uniform_distribution(1, main.L_cumul_nbr_brkpt_per_seg.get_tot_cumul_freq());
        long int h_choosen_brkpt = uniform_distribution(rand_gen.Seed_gen);

        int choosen_seg_index = main.L_cumul_nbr_brkpt_per_seg.seg_index(h_choosen_brkpt);
        auto y_choosen_seg = main.Pool_seg.get_obj_from_index(choosen_seg_index);

        auto k_choose_brkpt = y_choosen_seg->R_right_brkpt - main.L_cumul_nbr_brkpt_per_seg.get_cumul_freq_at_index(choosen_seg_index) + h_choosen_brkpt - 1;

        recomb.ini(y_choosen_seg, k_choose_brkpt);
        recomb.break_within_segment(main.Pool_seg, main.L_cumul_nbr_brkpt_per_seg);

        REQUIRE(main.Pool_seg.get_index_from_obj(recomb.Y_choosen_seg) == 0);
        REQUIRE(indiv_vec[0]->get_chr(chr_index, 0)->L_left_brkpt == 0);
        REQUIRE(recomb.Z_new_seg->R_right_brkpt == 10);
        REQUIRE(!(indiv_vec[0]->get_chr(chr_index, 0)->Next_seg));
        REQUIRE(!(recomb.Z_new_seg->Previous_seg));
        REQUIRE(recomb.X_prev_seg == recomb.Y_choosen_seg);
    }

    SECTION("break_within_seg_two_sub_seg")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                    {0, 0}};
        samp_param.Ploidy = 1;
        samp_param.n_total_sample_size = 2;
        samp_param.Sequence_length = 10;
        samp_param.Chr_nbr = 1;
        int chr_index = 0;

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
        main.ini(indiv_vec, samp_param, chr_index);

        seg seg0;
        seg seg1;

        seg0 = seg(0, 5, 0, nullptr, &seg1);
        seg1 = seg(5, 10, 1, &seg0, nullptr);

        (indiv_vec[0])->update_indiv(chr_index, &seg0, nullptr);
        main.L_cumul_nbr_brkpt_per_seg.up_at_index_with_value(0, -10);

        alg_recomb_c recomb;

        std::uniform_int_distribution<long int> uniform_distribution(1, main.L_cumul_nbr_brkpt_per_seg.get_tot_cumul_freq());
        long int h_choosen_brkpt = uniform_distribution(rand_gen.Seed_gen);

        int choosen_seg_index = main.L_cumul_nbr_brkpt_per_seg.seg_index(h_choosen_brkpt);
        auto y_choosen_seg = main.Pool_seg.get_obj_from_index(choosen_seg_index);

        auto k_choose_brkpt = y_choosen_seg->R_right_brkpt - main.L_cumul_nbr_brkpt_per_seg.get_cumul_freq_at_index(choosen_seg_index) + h_choosen_brkpt - 1;

        bool result = recomb.ini(y_choosen_seg, k_choose_brkpt);

        if (result)
        {
            recomb.break_within_segment(main.Pool_seg, main.L_cumul_nbr_brkpt_per_seg);

            if (recomb.Y_choosen_seg == indiv_vec[0]->get_chr(chr_index, 0))
            {
                REQUIRE(recomb.Z_new_seg->R_right_brkpt == 5);
                REQUIRE(!(recomb.Z_new_seg->Previous_seg));
                REQUIRE(recomb.Z_new_seg->Next_seg == indiv_vec[1]->get_chr(chr_index, 0));
            }
            else
            {
                REQUIRE(recomb.Z_new_seg->R_right_brkpt == 10);
                REQUIRE(!(recomb.Z_new_seg->Previous_seg));
                REQUIRE(!(recomb.Z_new_seg->Next_seg));
            }
        }
    }
}
