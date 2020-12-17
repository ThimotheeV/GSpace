#define CATCH_CONFIG_MAIN FunctionnalTestMutation
#include "catch.hpp"

#include "simulator.hpp"
#include <iostream>

/*********IAM**********/

TEST_CASE("Hudson__IAM__Wright_Fisher_model__1_pop__2_theta__1_site__4000000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);
    // 7500000
    simu_param.Repetition_nbr = 4000000;
    simu_param.Continuous_time_approxim = true;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 50000;
    demo_param.Population_size_N = 50000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.00001;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.333333300000;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Hudson__IAM__Wright_Fisher_model__1_pop__2_theta__10_indepsite__400000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);
    // 7500000
    simu_param.Repetition_nbr = 400000;
    simu_param.Continuous_time_approxim = true;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 50000;
    demo_param.Population_size_N = 50000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Chr_nbr = 10;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.00001;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.333333300000;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__IAM__Wright_Fisher_model__1_pop__5000_diplo_indivs_per_pop__0.0001_mut_rate__1_site__400000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 400000;
    simu_param.Continuous_time_approxim = false;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 5000;
    demo_param.Population_size_N = 5000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.0001;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.3332999989;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__IAM__Wright_Fisher_model__1_pop__5000_diplo_indivs_per_pop__0.0001_mut_rate__10_indepsite__30000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 30000;
    simu_param.Continuous_time_approxim = false;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 5000;
    demo_param.Population_size_N = 5000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Chr_nbr = 10;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.0001;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.3332999989;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

/*********KAM**********/

TEST_CASE("Hudson__KAM_1_5__Wright_Fisher_model__1_pop__2_theta__1_site__1500000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 1500000;
    simu_param.Continuous_time_approxim = true;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 50000;
    demo_param.Population_size_N = 50000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::kam;
    muta_param.K_min = 1;
    muta_param.K_max = 5;
    muta_param.Unscaled_mut_rate_mu = 0.00001;

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    std::uniform_int_distribution distrib(muta_param.K_min, muta_param.K_max);
    muta_param.Ancestry_seq.resize(samp_param.Chr_nbr);

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
        {
            muta_param.Ancestry_seq[chr].resize(samp_param.Sequence_length);
            for (int i = 0; i < samp_param.Sequence_length; ++i)
            {
                muta_param.Ancestry_seq[chr][i] = distrib(rand_gen.Seed_gen);
            }
        }
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.428571397959;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Hudson__KAM_1_5__Wright_Fisher_model__1_pop__2_theta__10_indepsite__600000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 600000;
    simu_param.Continuous_time_approxim = true;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 50000;
    demo_param.Population_size_N = 50000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Chr_nbr = 10;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::kam;
    muta_param.K_min = 1;
    muta_param.K_max = 5;
    muta_param.Unscaled_mut_rate_mu = 0.00001;

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    std::uniform_int_distribution distrib(muta_param.K_min, muta_param.K_max);
    muta_param.Ancestry_seq.resize(samp_param.Chr_nbr);

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
        {
            muta_param.Ancestry_seq[chr].resize(samp_param.Sequence_length);
            for (int i = 0; i < samp_param.Sequence_length; ++i)
            {
                muta_param.Ancestry_seq[chr][i] = distrib(rand_gen.Seed_gen);
            }
        }
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.428571397959;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__KAM_1_5__Wright_Fisher_model__1_pop__5000_diplo_indivs_per_pop__0.0001_mut_rate__1_site__500000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 500000;
    simu_param.Continuous_time_approxim = false;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 5000;
    demo_param.Population_size_N = 5000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::kam;
    muta_param.K_min = 1;
    muta_param.K_max = 5;
    muta_param.Unscaled_mut_rate_mu = 0.0001;

    output_stat_c output;
    output.Stat_out = false;

    std::uniform_int_distribution distrib(3, 3);
    muta_param.Ancestry_seq.resize(samp_param.Chr_nbr);

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
        {
            muta_param.Ancestry_seq[chr].resize(samp_param.Sequence_length);
            for (int i = 0; i < samp_param.Sequence_length; ++i)
            {
                muta_param.Ancestry_seq[chr][i] = distrib(rand_gen.Seed_gen);
            }
        }
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.4285408153;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__KAM_1_5__Wright_Fisher_model__1_pop__5000_diplo_indivs_per_pop__0.0001_mut_rate__10_indepsite__75000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 75000;
    simu_param.Continuous_time_approxim = false;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 5000;
    demo_param.Population_size_N = 5000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Chr_nbr = 10;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::kam;
    muta_param.K_min = 1;
    muta_param.K_max = 5;
    muta_param.Unscaled_mut_rate_mu = 0.0001;

    output_stat_c output;
    output.Stat_out = false;

    std::uniform_int_distribution distrib(muta_param.K_min, muta_param.K_max);
    muta_param.Ancestry_seq.resize(samp_param.Chr_nbr);

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
        {
            muta_param.Ancestry_seq[chr].resize(samp_param.Sequence_length);
            for (int i = 0; i < samp_param.Sequence_length; ++i)
            {
                muta_param.Ancestry_seq[chr][i] = distrib(rand_gen.Seed_gen);
            }
        }
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.4285408153;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

/********SMM***********/

TEST_CASE("Hudson__SMM_1_10__Wright_Fisher_model__1_pop__2_theta__1_site__1500000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 1500000;
    simu_param.Continuous_time_approxim = true;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 50000;
    demo_param.Population_size_N = 50000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::smm;
    muta_param.K_min = 1;
    muta_param.K_max = 10;
    muta_param.Unscaled_mut_rate_mu = 0.00001;

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    std::uniform_int_distribution distrib(muta_param.K_min, muta_param.K_max);
    muta_param.Ancestry_seq.resize(samp_param.Chr_nbr);

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
        {
            muta_param.Ancestry_seq[chr].resize(samp_param.Sequence_length);
            for (int i = 0; i < samp_param.Sequence_length; ++i)
            {
                muta_param.Ancestry_seq[chr][i] = distrib(rand_gen.Seed_gen);
            }
        }
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.487188035355;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Hudson__SMM_1_10000__Wright_Fisher_model__1_pop__2_theta__1_site__3000000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 3000000;
    simu_param.Continuous_time_approxim = true;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 50000;
    demo_param.Population_size_N = 50000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::smm;
    muta_param.K_min = 1;
    muta_param.K_max = 10000;
    muta_param.Unscaled_mut_rate_mu = 0.00001;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 5000));

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.447213564274;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__SMM_1_10__Wright_Fisher_model__1_pop__5000_diplo_indivs_per_pop__0.0001_mut_rate__1_site__1000000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 1000000;
    simu_param.Continuous_time_approxim = false;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 5000;
    demo_param.Population_size_N = 5000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::smm;
    muta_param.K_min = 1;
    muta_param.K_max = 10;
    muta_param.Unscaled_mut_rate_mu = 0.0001;

    output_stat_c output;
    output.Stat_out = false;

    std::uniform_int_distribution distrib(muta_param.K_min, muta_param.K_max);
    muta_param.Ancestry_seq.resize(samp_param.Chr_nbr);

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
        {
            muta_param.Ancestry_seq[chr].resize(samp_param.Sequence_length);
            for (int i = 0; i < samp_param.Sequence_length; ++i)
            {
                muta_param.Ancestry_seq[chr][i] = distrib(rand_gen.Seed_gen);
            }
        }
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.487188035355;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

/********GSM***********/

TEST_CASE("Hudson__GSM_1_10__Wright_Fisher_model__1_pop__2_theta__1_site__1000000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 1000000;
    simu_param.Continuous_time_approxim = true;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 50000;
    demo_param.Population_size_N = 50000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::gsm;
    muta_param.K_min = 1;
    muta_param.K_max = 10;
    muta_param.P_gsm = 0.22;
    muta_param.Unscaled_mut_rate_mu = 0.00001;

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    std::uniform_int_distribution distrib(muta_param.K_min, muta_param.K_max);
    muta_param.Ancestry_seq.resize(samp_param.Chr_nbr);

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
        {
            muta_param.Ancestry_seq[chr].resize(samp_param.Sequence_length);
            for (int i = 0; i < samp_param.Sequence_length; ++i)
            {
                muta_param.Ancestry_seq[chr][i] = distrib(rand_gen.Seed_gen);
            }
        }
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.497653221634;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

// TODO : vérifier le pGSM => problème de définition
TEST_CASE("Hudson__GSM_1_10000__Wright_Fisher_model__1_pop__2_theta__1_site__1000000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 1000000;
    simu_param.Continuous_time_approxim = true;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 5000;
    demo_param.Population_size_N = 5000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::gsm;
    muta_param.K_min = 1;
    muta_param.K_max = 10000;
    muta_param.P_gsm = 1 - 0.74;
    muta_param.Unscaled_mut_rate_mu = 0.0001;

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 5000));

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.35607658248;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__GSM_1_10__Wright_Fisher_model__1_pop__5000_diplo_indivs_per_pop__0.0001_mut_rate__1_site__750000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 750000;
    simu_param.Continuous_time_approxim = false;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 5000;
    demo_param.Population_size_N = 5000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::gsm;
    muta_param.K_min = 1;
    muta_param.K_max = 10;
    muta_param.P_gsm = 0.22;
    muta_param.Unscaled_mut_rate_mu = 0.0001;

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    std::uniform_int_distribution distrib(muta_param.K_min, muta_param.K_max);
    muta_param.Ancestry_seq.resize(samp_param.Chr_nbr);

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
        {
            muta_param.Ancestry_seq[chr].resize(samp_param.Sequence_length);
            for (int i = 0; i < samp_param.Sequence_length; ++i)
            {
                muta_param.Ancestry_seq[chr][i] = distrib(rand_gen.Seed_gen);
            }
        }
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.497653221634;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

/********JCM***********/

TEST_CASE("Hudson__JCM__Wright_Fisher_model__1_pop__2_theta__1_site__3000000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 3000000;
    simu_param.Continuous_time_approxim = true;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 50000;
    demo_param.Population_size_N = 50000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::jcm;
    muta_param.Unscaled_mut_rate_mu = 0.00001;

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 2));

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    //From IBDSim with fixe MRCA_seq of size 1
    auto Q0 = 0.4545074723;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Hudson__JCM__Wright_Fisher_model__1_pop__2_theta__100_indepsite__65000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 65000;
    simu_param.Continuous_time_approxim = true;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 50000;
    demo_param.Population_size_N = 50000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Chr_nbr = 100;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::jcm;
    muta_param.Unscaled_mut_rate_mu = 0.00001;

    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    output_stat_c output;
    output.Stat_out = false;

    std::uniform_int_distribution distrib(1, 4);
    muta_param.Ancestry_seq.resize(samp_param.Chr_nbr);

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
        {
            muta_param.Ancestry_seq[chr].resize(samp_param.Sequence_length);
            for (int i = 0; i < samp_param.Sequence_length; ++i)
            {
                muta_param.Ancestry_seq[chr][i] = distrib(rand_gen.Seed_gen);
            }
        }
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    //From IBDSim
    auto Q0 = 0.4545074723;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__JCM__Wright_Fisher_model__1_pop__5000_diplo_indivs_per_pop__0.0001_mut_rate__1_site__300000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 300000;
    simu_param.Continuous_time_approxim = false;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 5000;
    demo_param.Population_size_N = 5000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::kam;
    muta_param.K_min = 1;
    muta_param.K_max = 4;
    muta_param.Unscaled_mut_rate_mu = 0.0001;

    output_stat_c output;
    output.Stat_out = false;

    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 2));

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    //From IBDSim with fixe MRCA_seq of size 1
    auto Q0 = 0.4545074723;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__JCM__Wright_Fisher_model__1_pop__5000_diplo_indivs_per_pop__0.0001_mut_rate__10_indepsite__45000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 45000;
    simu_param.Continuous_time_approxim = false;

    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;

    demo_param.Lattice_size = {0, 0};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    samp_param.Sample_size_per_node = std::vector<int>(1, 100);
    samp_param.n_total_sample_size = 100;

    demo_param.Pop_size_per_node = 5000;
    demo_param.Population_size_N = 5000;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x <= demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y <= demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 1;
    samp_param.Chr_nbr = 10;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::jcm;
    muta_param.Unscaled_mut_rate_mu = 0.0001;

    output_stat_c output;
    output.Stat_out = false;

    std::uniform_int_distribution distrib(1, 4);
    muta_param.Ancestry_seq.resize(samp_param.Chr_nbr);

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
        {
            muta_param.Ancestry_seq[chr].resize(samp_param.Sequence_length);
            for (int i = 0; i < samp_param.Sequence_length; ++i)
            {
                muta_param.Ancestry_seq[chr][i] = distrib(rand_gen.Seed_gen);
            }
        }
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    //From IBDSim with fixe MRCA_seq of size 1
    auto Q0 = 0.4545074723;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

/********K80***********/

/**********F81*********/

/*********HKY**********/

/*********TN9**********/