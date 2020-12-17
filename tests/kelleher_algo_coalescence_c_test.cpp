#define CATCH_CONFIG_MAIN AlgoCoalescenceCTest

#include <iostream>

#include "kelleher_algo.hpp"
#include "catch.hpp"

// TEST_CASE()
// {
//     REQUIRE(1);
// }

TEST_CASE("coa_algo_c_test")
{
    SECTION("simple_coa_algo")
    {
        struct_arg_c coa;

        seg *seg1 = coa.Pool_seg.get_new_obj();
        seg *seg2 = coa.Pool_seg.get_new_obj();
        seg *seg3 = coa.Pool_seg.get_new_obj();

        *seg1 = seg(0, 1, 0, nullptr, nullptr, nullptr);
        *seg2 = seg(2, 3, 0, nullptr, nullptr, nullptr);
        *seg3 = seg(5, 6, 0, nullptr, nullptr, nullptr);

        seg1->Next_seg = seg2;
        seg2->Previous_seg = seg1;

        simple_coa(coa, seg1, seg3);

        REQUIRE(seg2->Next_seg == seg3);
        REQUIRE(seg3->Previous_seg == seg2);
        REQUIRE(coa.L_cumul_nbr_brkpt_per_seg.get_frequency_index(coa.Pool_seg.get_index_from_obj(seg3)) == 3);

        REQUIRE_THROWS(simple_coa(coa, seg3, seg1));

        auto seg4 = new seg(5, 6, 0, seg1, nullptr, nullptr);

        REQUIRE_THROWS(simple_coa(coa, seg1, seg4));
    }

    SECTION("choose_ancestors")
    {
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                    {0, 0}};
        samp_param.Ploidy = 1;

        samp_param.n_total_sample_size = 2;
        samp_param.Chr_nbr = 1;
        samp_param.Sequence_length = 10;

        struct_arg_c main;
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

        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        indiv_vec.ini(lat, samp_param);
        
        auto chr_index = 0;
        main.ini(indiv_vec, samp_param, chr_index);

        alg_coa_c coa;

        coa.choose_ancestors(indiv_vec[0]->get_chr(0, 0), indiv_vec[1]->get_chr(0, 0));

        REQUIRE(coa.X_prev_seg);
        REQUIRE(coa.Y_choosen_seg);
        REQUIRE(coa.X_prev_seg != coa.Y_choosen_seg);
    }

    SECTION("choose_case_echange_place_and_no_overlap_between_seg")
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

        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        indiv_vec.ini(lat, samp_param);
        
        auto chr_index = 0;
        main.ini(indiv_vec, samp_param, chr_index);

        alg_coa_c coa;

        coa.X_prev_seg = main.Pool_seg.get_new_obj();
        *coa.X_prev_seg = seg(6, 10, 1);
        coa.Y_choosen_seg = main.Pool_seg.get_new_obj();
        *coa.Y_choosen_seg = seg(0, 5, 2);

        REQUIRE(!(coa.choose_case(main.Pool_seg)));
        REQUIRE(!(coa.X_prev_seg));
        REQUIRE(coa.Y_choosen_seg->R_right_brkpt == 10);
    }

    SECTION("choose_case_partial_overleap_seg")
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

        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        indiv_vec.ini(lat, samp_param);
        
        auto chr_index = 0;
        main.ini(indiv_vec, samp_param, chr_index);

        alg_coa_c coa;

        coa.X_prev_seg = main.Pool_seg.get_new_obj();
        *coa.X_prev_seg = seg(0, 8, 1);
        coa.Y_choosen_seg = main.Pool_seg.get_new_obj();
        *coa.Y_choosen_seg = seg(4, 10, 2);

        REQUIRE(!(coa.choose_case(main.Pool_seg)));
        REQUIRE(coa.Alpha_new_seg->R_right_brkpt == 4);
        REQUIRE(coa.Alpha_new_seg->U_ident_node == 1);
        REQUIRE(coa.X_prev_seg->L_left_brkpt == 4);
    }

    SECTION("choose_case_total_overleap_seg")
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

        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        indiv_vec.ini(lat, samp_param);
        
        auto chr_index = 0;
        main.ini(indiv_vec, samp_param, chr_index);

        alg_coa_c coa;

        coa.X_prev_seg = main.Pool_seg.get_new_obj();
        *coa.X_prev_seg = seg(0, 10, 1);
        coa.Y_choosen_seg = main.Pool_seg.get_new_obj();
        *coa.Y_choosen_seg = seg(0, 10, 2);

        REQUIRE(coa.choose_case(main.Pool_seg));
    }

    SECTION("coa")
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

        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        indiv_vec.ini(lat, samp_param);
        
        auto chr_index = 0;
        main.ini(indiv_vec, samp_param, chr_index);

        alg_coa_c coa;

        coa.X_prev_seg = main.Pool_seg.get_new_obj();
        *coa.X_prev_seg = seg(5, 10, 1);
        coa.Y_choosen_seg = main.Pool_seg.get_new_obj();
        *coa.Y_choosen_seg = seg(5, 9, 2);

        REQUIRE(coa.coa(main.S_intersection_count, main.W_next_node_nbr));
        REQUIRE(main.W_next_node_nbr == 5);
        REQUIRE(coa.L_left_brkpt == 5);
        REQUIRE(coa.Rp_min_right_brkpt == 9);
        REQUIRE(main.S_intersection_count.count(5) == 1);
        REQUIRE(main.S_intersection_count.count(9) == 1);
    }

    SECTION("seg_mrca")
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

        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        indiv_vec.ini(lat, samp_param);
        
        auto chr_index = 0;
        main.ini(indiv_vec, samp_param, chr_index);

        alg_coa_c coa;

        coa.X_prev_seg = main.Pool_seg.get_new_obj();
        *coa.X_prev_seg = seg(5, 10, 1);
        coa.Y_choosen_seg = main.Pool_seg.get_new_obj();
        *coa.Y_choosen_seg = seg(5, 9, 2);

        coa.coa(main.S_intersection_count, main.W_next_node_nbr);
        coa.seg_mrca(main.S_intersection_count);

        REQUIRE(main.S_intersection_count.find(5)->second == 0);
        REQUIRE(coa.R_right_brkpt == 9);
    }

    SECTION("decrement_overlaps")
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

        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        indiv_vec.ini(lat, samp_param);
        
        auto chr_index = 0;
        main.ini(indiv_vec, samp_param, chr_index);

        alg_coa_c coa;

        coa.X_prev_seg = main.Pool_seg.get_new_obj();
        *coa.X_prev_seg = seg(5, 10, 1);
        coa.Y_choosen_seg = main.Pool_seg.get_new_obj();
        *coa.Y_choosen_seg = seg(5, 9, 2);

        coa.coa(main.S_intersection_count, main.W_next_node_nbr);
        coa.decrement_overlaps(main.Pool_seg, main.S_intersection_count);

        REQUIRE(main.S_intersection_count.find(5)->second == 3);
        REQUIRE(main.S_intersection_count.find(9)->second == 4);
        REQUIRE(coa.Alpha_new_seg->L_left_brkpt == 5);
        REQUIRE(coa.Alpha_new_seg->R_right_brkpt == 9);
        REQUIRE(coa.Alpha_new_seg->U_ident_node == 4);
    }

    SECTION("update_segs")
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

        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        indiv_vec.ini(lat, samp_param);
        
        auto chr_index = 0;
        main.ini(indiv_vec, samp_param, chr_index);

        alg_coa_c coa;
        coa_table_c coa_table;

        coa.X_prev_seg = main.Pool_seg.get_new_obj();
        *coa.X_prev_seg = seg(5, 10, 1);
        coa.Y_choosen_seg = main.Pool_seg.get_new_obj();
        *coa.Y_choosen_seg = seg(5, 9, 2);

        coa.coa(main.S_intersection_count, main.W_next_node_nbr);
        coa.decrement_overlaps(main.Pool_seg, main.S_intersection_count);
        coa.update_segs(main.Pool_seg, main.L_cumul_nbr_brkpt_per_seg, coa_table, main.T_time, -1);

        REQUIRE(coa.X_prev_seg->L_left_brkpt == 9);
        REQUIRE(!(coa.Y_choosen_seg));
        REQUIRE(coa_table.get_coa_table().size() == 1);
        REQUIRE(coa_table_c::get_left_brkpt_l(coa_table.get_coa_table()[0]) == 5);
        REQUIRE(coa_table_c::get_right_brkpt_r(coa_table.get_coa_table()[0]) == 9);
        REQUIRE(coa_table_c::get_num_node_u(coa_table.get_coa_table()[0]) == 4);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table.get_coa_table()[0]))[0] == 1);
        REQUIRE((coa_table_c::get_childs_node_c(coa_table.get_coa_table()[0]))[1] == 2);
        REQUIRE(coa_table_c::get_time_t(coa_table.get_coa_table()[0]) == 0);
        REQUIRE(coa_table_c::get_gen(coa_table.get_coa_table()[0]) == -1);
    }

    SECTION("update_brkpt_no_alpha")
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

        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        indiv_vec.ini(lat, samp_param);
        
        auto chr_index = 0;
        main.ini(indiv_vec, samp_param, chr_index);

        alg_coa_c coa;

        coa.update_brkpts(main.Pool_seg, main.L_cumul_nbr_brkpt_per_seg);

        REQUIRE(coa.New_first_seg == nullptr);
    }

    SECTION("update_brkpt_no_seg_z")
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

        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        indiv_vec.ini(lat, samp_param);
        
        auto chr_index = 0;
        main.ini(indiv_vec, samp_param, chr_index);

        alg_coa_c coa;
        coa_table_c coa_table;

        coa.X_prev_seg = main.Pool_seg.get_new_obj();
        *coa.X_prev_seg = seg(5, 10, 1);
        coa.Y_choosen_seg = main.Pool_seg.get_new_obj();
        *coa.Y_choosen_seg = seg(5, 9, 2);

        coa.coa(main.S_intersection_count, main.W_next_node_nbr);
        coa.decrement_overlaps(main.Pool_seg, main.S_intersection_count);
        coa.update_segs(main.Pool_seg, main.L_cumul_nbr_brkpt_per_seg, coa_table, main.T_time, -1);
        coa.update_brkpts(main.Pool_seg, main.L_cumul_nbr_brkpt_per_seg);

        REQUIRE(indiv_vec.size() == 4);
        REQUIRE(coa.New_first_seg->L_left_brkpt == 5);                       //index(6) = 4 in init + 1coa.Seg_* + 1 in decrement_overlaps : last one in update_brkpts
        REQUIRE(main.L_cumul_nbr_brkpt_per_seg.get_frequency_index(6) == 3); //number of brkpt in the seg
        REQUIRE(!(coa.Z_new_seg->Previous_seg));
    }

    SECTION("update_brkpt_with_seg_z")
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

        indiv_stock_c indiv_vec(samp_param.Chr_nbr);
        indiv_vec.ini(lat, samp_param);
        
        auto chr_index = 0;
        main.ini(indiv_vec, samp_param, chr_index);

        alg_coa_c coa;
        coa_table_c coa_table;

        coa.X_prev_seg = main.Pool_seg.get_new_obj();
        *coa.X_prev_seg = seg(5, 10, 1);
        coa.Y_choosen_seg = main.Pool_seg.get_new_obj();
        *coa.Y_choosen_seg = seg(5, 9, 2);
        coa.Z_new_seg = main.Pool_seg.get_new_obj();
        *coa.Z_new_seg = seg(0, 5, 3);

        auto save_Z_new_seg = coa.Z_new_seg;

        coa.coa(main.S_intersection_count, main.W_next_node_nbr);
        coa.decrement_overlaps(main.Pool_seg, main.S_intersection_count);
        coa.update_segs(main.Pool_seg, main.L_cumul_nbr_brkpt_per_seg, coa_table, main.T_time, -1);
        coa.update_brkpts(main.Pool_seg, main.L_cumul_nbr_brkpt_per_seg);

        REQUIRE(coa.Alpha_new_seg->Previous_seg == save_Z_new_seg);
        REQUIRE(main.L_cumul_nbr_brkpt_per_seg.get_frequency_index(7) == 4); //index(7) = 4 in init + 2coa.Seg_* + 1 in decrement_overlaps
        REQUIRE(save_Z_new_seg->Next_seg == coa.Alpha_new_seg);              //4 = number of brkpt between right(Z_new_seg) and right(Y_choosen_seg)
        REQUIRE(!(coa.Bool_defrag));
    }

    SECTION("prune_avl_tree")
    {

        alg_coa_c coa;

        std::map<int, int> S_intersection_count;

        S_intersection_count.insert(std::pair<int, int>(0, 5));
        S_intersection_count.insert(std::pair<int, int>(1, 10));
        S_intersection_count.insert(std::pair<int, int>(2, 10));
        S_intersection_count.insert(std::pair<int, int>(3, 6));
        S_intersection_count.insert(std::pair<int, int>(4, 6));
        S_intersection_count.insert(std::pair<int, int>(5, 10));

        coa.prune_tree(S_intersection_count);

        REQUIRE(S_intersection_count.find(0)->second == 5);
        REQUIRE(S_intersection_count.find(1)->second == 10);
        REQUIRE(S_intersection_count.count(2) == 0);
        REQUIRE(S_intersection_count.find(3)->second == 6);
        REQUIRE(S_intersection_count.count(4) == 0);
        REQUIRE(S_intersection_count.find(5)->second == 10);
    }

    SECTION("reshape_segs")
    {

        obj_mem_manag_c<seg> Pool_seg;
        fenwick_tree_c cumul_nbr_brkpt_per_seg;
        cumul_nbr_brkpt_per_seg.create_fenwick_tree(6);

        seg *Seg_one = Pool_seg.get_new_obj();
        seg *Seg_two = Pool_seg.get_new_obj();
        seg *Seg_three = Pool_seg.get_new_obj();
        seg *Seg_four = Pool_seg.get_new_obj();
        seg *Seg_five = Pool_seg.get_new_obj();
        seg *Seg_six = Pool_seg.get_new_obj();

        *Seg_one = seg(0, 10, 1, nullptr, Seg_two);
        *Seg_two = seg(10, 20, 2, Seg_one, Seg_three);
        *Seg_three = seg(21, 29, 3, Seg_two, Seg_four);
        *Seg_four = seg(30, 40, 3, Seg_three, Seg_five);
        *Seg_five = seg(40, 50, 3, Seg_four, Seg_six);
        *Seg_six = seg(50, 60, 3, Seg_five, nullptr);

        std::size_t index = Pool_seg.get_index_from_obj(Seg_one);
        cumul_nbr_brkpt_per_seg.set_value_at_index(index, 9);
        index = Pool_seg.get_index_from_obj(Seg_two);
        cumul_nbr_brkpt_per_seg.set_value_at_index(index, 9);
        index = Pool_seg.get_index_from_obj(Seg_three);
        cumul_nbr_brkpt_per_seg.set_value_at_index(index, 7);
        index = Pool_seg.get_index_from_obj(Seg_four);
        cumul_nbr_brkpt_per_seg.set_value_at_index(index, 9);
        index = Pool_seg.get_index_from_obj(Seg_five);
        cumul_nbr_brkpt_per_seg.set_value_at_index(index, 9);
        index = Pool_seg.get_index_from_obj(Seg_six);
        cumul_nbr_brkpt_per_seg.set_value_at_index(index, 9);

        REQUIRE(cumul_nbr_brkpt_per_seg.get_tot_cumul_freq() == 52);

        alg_coa_c coa;
        coa.Z_new_seg = Seg_six;
        coa.reshape_segs(Pool_seg, cumul_nbr_brkpt_per_seg);

        // REQUIRE_THROWS(Pool_seg.get_index_from_obj(Seg_five));
        // REQUIRE_THROWS(Pool_seg.get_index_from_obj(Seg_six));
        REQUIRE(Seg_four->L_left_brkpt == 30);
        REQUIRE(Seg_four->R_right_brkpt == 60);
        REQUIRE(Seg_four->U_ident_node == 3);
        REQUIRE(Seg_four->Previous_seg == Seg_three);
        REQUIRE(!(Seg_four->Next_seg));
    }
}
