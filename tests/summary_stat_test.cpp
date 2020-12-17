#define CATCH_CONFIG_MAIN ProbaIndentest
#include "catch.hpp"
#include <iostream>

#include "summary_stat.hpp"

TEST_CASE("summary_stat")
{
    SECTION("combination_sum_stat")
    {
        auto result = combination(2, 4);
        REQUIRE(result == 6);

        result = combination(3, 3);
        REQUIRE(result == 1);

        REQUIRE_THROWS(combination(8, 7));
    }

    SECTION("moy_cumul_prob_id_1_2_loc_c")
    {
        std::array<double, 2> accu_proba_moy = {0, 0};

        std::vector<double> vec = {0.8, 0.83, 0.83, 0.73};
        accu_proba_moy = calc_cumul_mean_var(accu_proba_moy, vec, 0);

        vec = {0.33, 0.5, 0.4, 0.6};
        accu_proba_moy = calc_cumul_mean_var(accu_proba_moy, vec, 1);

        REQUIRE(std::get<0>(accu_proba_moy) == Approx(0.6275).margin(0.0001));
        REQUIRE(std::get<1>(accu_proba_moy) == Approx(0.0349).margin(0.0001));
    }

    SECTION("moy_cumul_one_value")
    {
        std::array<double, 2> accu_proba_moy = {0, 0};

        std::vector<double> vec = {0.8, 0.83, 0.83, 0.73, 0.33, 0.5, 0.4, 0.6};

        for (std::size_t i = 0; i < vec.size(); ++i)
        {
            accu_proba_moy = calc_cumul_mean_var(accu_proba_moy, {vec[i]}, i);
        }

        REQUIRE(std::get<0>(accu_proba_moy) == Approx(0.6275).margin(0.0001));
        REQUIRE(std::get<1>(accu_proba_moy) == Approx(0.0349).margin(0.0001));
    }

    SECTION("2_haplo_samples_coa_tree_metrics_c")
    {
        int population_size_N = 10;
        int n_total_sample_size = 2;
        coa_tree_metrics_c coa_tree_metric;
        coa_tree_metric.calcul_coa_tree_metrics(false, 1, population_size_N, n_total_sample_size);

        REQUIRE(coa_tree_metric.Theo_MRCA_mean == 20);
    }

    SECTION("10_haplo_samples_coa_tree_metrics_c")
    {
        int population_size_N = 10;
        int n_total_sample_size = 10;
        coa_tree_metrics_c coa_tree_metric;
        coa_tree_metric.calcul_coa_tree_metrics(false, 1, population_size_N, n_total_sample_size);

        REQUIRE(coa_tree_metric.Theo_MRCA_mean == 36);
    }

    SECTION("calc_prob_id_1_loc_simulator")
    {
        std::vector<std::array<int, 2>> coord_vec{{0, 0}, {0, 0}, {0, 0}, {0, 1}, {0, 1}, {0, 1}};

        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = coord_vec;
        samp_param.Sample_size_per_node = {3};

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        samp_param.Sequence_length = 6;
        samp_param.Chr_nbr = 1;
        auto &muta_param = singleton_c<muta_param_c>::instance();
        muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));
        samp_param.Dist_class_nbr = 2;

        Prob_id_1_loc_Qr_res_c q_result;

        std::vector<std::pair<int, int>> indiv1 = {{0, 2}, {2, 2}, {3, 0}, {4, 0}};
        std::vector<std::pair<int, int>> indiv2 = {{0, 2}, {1, 2}, {2, 2}, {4, 2}, {5, 0}};
        std::vector<std::pair<int, int>> indiv3 = {{0, 2}, {1, 0}, {2, 0}};
        std::vector<std::pair<int, int>> indiv4 = {{2, 2}, {5, 2}};
        std::vector<std::pair<int, int>> indiv5 = {{1, 2}, {2, 2}, {3, 2}};
        std::vector<std::pair<int, int>> indiv6 = {{1, 0}, {3, 2}, {5, 2}};

        std::vector<std::vector<std::vector<std::pair<int, int>>>> sample_mutated_state_vec = {{indiv1, indiv2, indiv3, indiv4, indiv5, indiv6}};

        data_plane_vec_c data_plane_vec(sample_mutated_state_vec, samp_param, muta_param.Ancestry_seq);
        q_result.calc_Qr(data_plane_vec, 1);

        REQUIRE(q_result.Qr[0] == Approx(0.4167).epsilon(0.001));
    }

    SECTION("calc_prob_id_1_loc_simulator_at_dist_2")
    {
        std::vector<std::array<int, 2>> coord_vec{{0, 0}, {0, 0}, {0, 1}, {0, 1}, {0, 2}, {0, 2}};

        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = coord_vec;
        samp_param.Sample_size_per_node = {2};

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        samp_param.Sequence_length = 6;
        samp_param.Chr_nbr = 1;
        auto &muta_param = singleton_c<muta_param_c>::instance();
        muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));
        samp_param.Dist_class_nbr = 3;

        Prob_id_1_loc_Qr_res_c q_result;

        std::vector<std::pair<int, int>> indiv1 = {{0, 2}, {2, 2}, {3, 0}, {4, 0}};
        std::vector<std::pair<int, int>> indiv2 = {{0, 2}, {1, 2}, {2, 2}, {4, 2}, {5, 0}};
        std::vector<std::pair<int, int>> indiv3 = {{0, 2}, {1, 0}, {2, 0}};
        std::vector<std::pair<int, int>> indiv4 = {{2, 2}, {5, 2}};
        std::vector<std::pair<int, int>> indiv5 = {{1, 2}, {2, 2}, {3, 2}};
        std::vector<std::pair<int, int>> indiv6 = {{1, 0}, {3, 2}, {5, 2}};

        std::vector<std::vector<std::vector<std::pair<int, int>>>> sample_mutated_state_vec = {{indiv1, indiv2, indiv3, indiv4, indiv5, indiv6}};

        data_plane_vec_c data_plane_vec(sample_mutated_state_vec, samp_param, muta_param.Ancestry_seq);
        q_result.calc_Qr(data_plane_vec, 1);

        REQUIRE(q_result.Qr[0] == Approx(0.3889).epsilon(0.0001));
        REQUIRE(q_result.Qr[1] == Approx(0.3750).epsilon(0.0001));
        REQUIRE(q_result.Qr[2] == Approx(0.1667).epsilon(0.001));
    }

    SECTION("cumul_prob_id_1_2_loc_Q_result_c")
    {
        std::vector<std::array<int, 2>> coord_vec{{0, 0}, {0, 0}, {0, 1}, {0, 1}};

        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = coord_vec;
        samp_param.Sample_size_per_node = {2};

        auto &rand_gen = singleton_c<rand_gen_c>::instance();
        rand_gen.put_seed(3);

        samp_param.Sequence_length = 3;
        samp_param.Chr_nbr = 1;
        auto &muta_param = singleton_c<muta_param_c>::instance();
        muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));
        samp_param.Dist_class_nbr = 2;

        auto &info_collector = singleton_c<info_collector_c>::instance();
        info_collector.Iterative_stats = true;

        Prob_id_1_loc_Qr_res_c q_result;

        std::vector<std::pair<int, int>> indiv1 = {{0, 2}, {2, 2}};
        std::vector<std::pair<int, int>> indiv2 = {{0, 2}, {1, 2}, {2, 2}};
        std::vector<std::pair<int, int>> indiv3 = {{0, 2}, {1, 0}};
        std::vector<std::pair<int, int>> indiv4 = {{0, 2}, {2, 0}};

        std::vector<std::vector<std::vector<std::pair<int, int>>>> sample_mutated_state_vec = {{indiv1, indiv2, indiv3, indiv4}};

        data_plane_vec_c data_plane_vec(sample_mutated_state_vec, samp_param, muta_param.Ancestry_seq);
        q_result.calc_Qr(data_plane_vec, 0);

        REQUIRE(q_result.Qr_cumul_m_v[0].at(0) == Approx(0.5).epsilon(0.1));

        std::vector<std::pair<int, int>> indiv5 = {{1, 0}};
        std::vector<std::pair<int, int>> indiv6 = {{0, 0}, {1, 0}, {2, 0}};
        std::vector<std::pair<int, int>> indiv7 = {{0, 0}};
        std::vector<std::pair<int, int>> indiv8 = {{0, 2}, {2, 0}};

        std::vector<std::vector<std::vector<std::pair<int, int>>>> sample_mutated_state_vec2 = {{indiv5, indiv6, indiv7, indiv8}};

        data_plane_vec.update_data_plane_vec(sample_mutated_state_vec2, muta_param.Ancestry_seq, samp_param);
        q_result.calc_Qr(data_plane_vec, 1);

        REQUIRE(q_result.Qr_cumul_m_v[0].at(0) == Approx(0.4167).epsilon(0.0001));
    }
}
