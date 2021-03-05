#define CATCH_CONFIG_MAIN FunctionnalTestMutation
#include "catch.hpp"

#include "simulator.hpp"
#include <iostream>


TEST_CASE("Gen_per_gen__IAM__Wright_Fisher_model__1_pop__20_haplo_indivs_per_pop__0_05_mut_rate__1_site__2000000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 2000000;
    simu_param.Continuous_time_approxim = false;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {1, 1};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    demo_param.Pop_size_per_node = 20;
    demo_param.Population_size_N = 20;

    samp_param.Sample_size_per_node = std::vector<int>(1, 10);
    samp_param.n_total_sample_size = 10;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x < demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y < demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 1;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 1;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.05;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.316389132340052;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__IAM__Wright_Fisher_model__1_pop__10_diplo_indivs_per_pop__0_05_mut_rate__1_site__2000000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 2000000;
    simu_param.Continuous_time_approxim = false;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {1, 1};
    demo_param.Disp_dist_max = {0, 0};
    demo_param.Proba_migr = 0;

    demo_param.Pop_size_per_node = 10;
    demo_param.Population_size_N = 10;

    samp_param.Sample_size_per_node = std::vector<int>(1, 5);
    samp_param.n_total_sample_size = 5;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x < demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y < demo_param.Lattice_size.at(1); ++y)
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
    muta_param.Unscaled_mut_rate_mu = 0.05;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.316389132340052;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__IAM__Island_model__2_pop__20_haplo_indivs_per_pop__0_0051_prob_migration__0_05_mut_rate__1_site__30000000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 30000000;
    simu_param.Continuous_time_approxim = false;

    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {1, 2};
    demo_param.Disp_dist_max = {0, 1};
    demo_param.Edge_effects = edge_effect_enum::circular;
    demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
    demo_param.Disp_func[1] = fwd_disp_distrib_c::uniform();
    demo_param.Proba_migr = 0.0051;

    demo_param.Pop_size_per_node = 20;
    demo_param.Population_size_N = 40;

    samp_param.Sample_size_per_node = std::vector<int>(2, 10);
    samp_param.n_total_sample_size = 20;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x < demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y < demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 2;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 1;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.05;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(1);
    auto Q0 = 0.296897551216174;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));

    auto m_Q1_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(0);
    auto v_Q1_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(1);
    auto Q1 = 0.0285126846953135;

    REQUIRE(m_Q1_seq == Approx(Q1).margin(0.0001));
    REQUIRE(Q1 < m_Q1_seq + 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q1 > m_Q1_seq - 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__IAM__Island_model__2_pop__10_diplo_indivs_per_pop__0_0051_prob_migration__0_05_mut_rate__1_site__20000000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 20000000;
    simu_param.Continuous_time_approxim = false;

    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {1, 2};
    demo_param.Disp_dist_max = {0, 1};
    demo_param.Edge_effects = edge_effect_enum::circular;
    demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
    demo_param.Disp_func[1] = fwd_disp_distrib_c::uniform();
    demo_param.Proba_migr = 0.0051;

    demo_param.Pop_size_per_node = 10;
    demo_param.Population_size_N = 20;

    samp_param.Sample_size_per_node = std::vector<int>(2, 5);
    samp_param.n_total_sample_size = 10;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x < demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y < demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 2;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.05;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Qw_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwi_cumul_m_v.at(0);
    auto v_Qw_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwi_cumul_m_v.at(1);
    auto Qw = 0.300513;

    auto m_Q0_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(1);
    auto Q0 = 0.297725;

    auto m_Q1_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(0);
    auto v_Q1_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(1);
    auto Q1 = 0.0285922;

    // std::cout << "m_Qw_seq : " << m_Qw_seq << std::endl;
    // std::cout << "v_Qw_seq : " << v_Qw_seq << std::endl;
    // std::cout << "*************" << std::endl;
    // std::cout << "m_Q0_seq : " << m_Q0_seq << std::endl;
    // std::cout << "v_Q0_seq : " << v_Q0_seq << std::endl;
    // std::cout << "*************" << std::endl;
    // std::cout << "m_Q1_seq : " << m_Q1_seq << std::endl;
    // std::cout << "v_Q1_seq : " << v_Q1_seq << std::endl;

    REQUIRE(m_Qw_seq == Approx(Qw).margin(0.0001));
    REQUIRE(Qw < m_Qw_seq + 1.96 * (std::sqrt(v_Qw_seq / simu_param.Repetition_nbr)));
    REQUIRE(Qw > m_Qw_seq - 1.96 * (std::sqrt(v_Qw_seq / simu_param.Repetition_nbr)));

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));

    REQUIRE(m_Q1_seq == Approx(Q1).margin(0.0001));
    REQUIRE(Q1 < m_Q1_seq + 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q1 > m_Q1_seq - 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__IAM__Island_model__2_pop__200_haplo_indivs_per_pop__0_0051_prob_migration__0_005_mut_rate__1_site__3000000_rep__1/10000_precision")
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
    simu_param.Continuous_time_approxim = false;

    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {1, 2};
    demo_param.Disp_dist_max = {0, 1};
    demo_param.Edge_effects = edge_effect_enum::circular;
    demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
    demo_param.Disp_func[1] = fwd_disp_distrib_c::uniform();
    demo_param.Proba_migr = 0.0051;

    demo_param.Pop_size_per_node = 200;
    demo_param.Population_size_N = 400;

    samp_param.Sample_size_per_node = std::vector<int>(2, 10);
    samp_param.n_total_sample_size = 20;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x < demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y < demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 2;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 1;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.005;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(1);
    auto Q0 = 0.247424328890905;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));

    auto m_Q1_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(0);
    auto v_Q1_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(1);
    auto Q1 = 0.126043705590188;

    REQUIRE(m_Q1_seq == Approx(Q1).margin(0.0001));
    REQUIRE(Q1 < m_Q1_seq + 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q1 > m_Q1_seq - 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__IAM__Island_model__2_pop__100_diplo_indivs_per_pop__0_0051_prob_migration__0_005_mut_rate__1_site__2000000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 2500000;
    simu_param.Continuous_time_approxim = false;

    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {1, 2};
    demo_param.Disp_dist_max = {0, 1};
    demo_param.Edge_effects = edge_effect_enum::circular;
    demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
    demo_param.Disp_func[1] = fwd_disp_distrib_c::uniform();
    demo_param.Proba_migr = 0.0051;

    demo_param.Pop_size_per_node = 100;
    demo_param.Population_size_N = 200;

    samp_param.Sample_size_per_node = std::vector<int>(2, 10);
    samp_param.n_total_sample_size = 20;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x < demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y < demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 2;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.005;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Qw_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwi_cumul_m_v.at(0);
    auto v_Qw_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwi_cumul_m_v.at(1);
    auto Qw = 0.248995;

    auto m_Q0_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(1);
    auto Q0 = 0.247736;

    auto m_Q1_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(0);
    auto v_Q1_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(1);
    auto Q1 = 0.126202;

    // std::cout << "m_Qw_seq : " << m_Qw_seq << std::endl;
    // std::cout << "v_Qw_seq : " << v_Qw_seq << std::endl;
    // std::cout << "*************" << std::endl;
    // std::cout << "m_Q0_seq : " << m_Q0_seq << std::endl;
    // std::cout << "v_Q0_seq : " << v_Q0_seq << std::endl;
    // std::cout << "*************" << std::endl;
    // std::cout << "m_Q1_seq : " << m_Q1_seq << std::endl;
    // std::cout << "v_Q1_seq : " << v_Q1_seq << std::endl;

    REQUIRE(m_Qw_seq == Approx(Qw).margin(0.0001));
    REQUIRE(Qw < m_Qw_seq + 1.96 * (std::sqrt(v_Qw_seq / simu_param.Repetition_nbr)));
    REQUIRE(Qw > m_Qw_seq - 1.96 * (std::sqrt(v_Qw_seq / simu_param.Repetition_nbr)));

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));

    REQUIRE(m_Q1_seq == Approx(Q1).margin(0.0001));
    REQUIRE(Q1 < m_Q1_seq + 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q1 > m_Q1_seq - 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__IAM__Island_model__4_pop__20_haplo_indivs_per_pop__0_0051_prob_migration__0_05_mut_rate__1_site__4000000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 4000000;
    simu_param.Continuous_time_approxim = false;

    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {1, 4};
    demo_param.Disp_dist_max = {0, 3};
    demo_param.Edge_effects = edge_effect_enum::circular;
    demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
    demo_param.Disp_func[1] = fwd_disp_distrib_c::uniform();
    demo_param.Proba_migr = 0.0051;

    demo_param.Pop_size_per_node = 20;
    demo_param.Population_size_N = 80;

    samp_param.Sample_size_per_node = std::vector<int>(4, 5);
    samp_param.n_total_sample_size = 20;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x < demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y < demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 4;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 1;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.05;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(1);
    auto Q0 = 0.295747861645301;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));

    auto m_Q1_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(0);
    auto v_Q1_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(1);
    auto Q1 = 0.0100648247276543;

    REQUIRE(m_Q1_seq == Approx(Q1).margin(0.0001));
    REQUIRE(Q1 < m_Q1_seq + 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q1 > m_Q1_seq - 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__IAM__Island_model__4_pop__10_diplo_indivs_per_pop__0_0051_prob_migration__0_05_mut_rate__1_site__20000000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 20000000;
    simu_param.Continuous_time_approxim = false;

    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {1, 4};
    demo_param.Disp_dist_max = {0, 3};
    demo_param.Edge_effects = edge_effect_enum::circular;
    demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
    demo_param.Disp_func[1] = fwd_disp_distrib_c::uniform();
    demo_param.Proba_migr = 0.0051;

    demo_param.Pop_size_per_node = 10;
    demo_param.Population_size_N = 40;

    samp_param.Sample_size_per_node = std::vector<int>(4, 2);
    samp_param.n_total_sample_size = 8;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x < demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y < demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 4;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 2;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.05;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Qw_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwi_cumul_m_v.at(0);
    auto v_Qw_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwi_cumul_m_v.at(1);
    auto Qw = 0.299574;

    auto m_Q0_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(1);
    auto Q0 = 0.296621;

    auto m_Q1_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(0);
    auto v_Q1_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(1);
    auto Q1 = 0.0100945;

    // std::cout << "m_Qw_seq : " << m_Qw_seq << std::endl;
    // std::cout << "v_Qw_seq : " << v_Qw_seq << std::endl;
    // std::cout << "*************" << std::endl;
    // std::cout << "m_Q0_seq : " << m_Q0_seq << std::endl;
    // std::cout << "v_Q0_seq : " << v_Q0_seq << std::endl;
    // std::cout << "*************" << std::endl;
    // std::cout << "m_Q1_seq : " << m_Q1_seq << std::endl;
    // std::cout << "v_Q1_seq : " << v_Q1_seq << std::endl;

    REQUIRE(m_Qw_seq == Approx(Qw).margin(0.0001));
    REQUIRE(Qw < m_Qw_seq + 1.96 * (std::sqrt(v_Qw_seq / simu_param.Repetition_nbr)));
    REQUIRE(Qw > m_Qw_seq - 1.96 * (std::sqrt(v_Qw_seq / simu_param.Repetition_nbr)));

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));

    REQUIRE(m_Q1_seq == Approx(Q1).margin(0.0001));
    REQUIRE(Q1 < m_Q1_seq + 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q1 > m_Q1_seq - 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
}

// WARNING : Too much variance for testing in moderage time (less than 1 day)
// *************
// Q0_goal : 0.927009974733469
// m_Q0_seq : 0.927164
// v_Q0_seq : 0.0047823
// *************
// Q1_goal : 0.591377342954598
// m_Q1_seq : 0.589927
// v_Q1_seq : 0.0565423
// *************
// TEST_CASE("Gen_per_gen__IAM__Island_model__20_pop__30_haplo_indivs_per_pop__0_0051_prob_migration__0_0001_mut_rate__1_site__10000000_rep__1/10000_precision")
// {
//     auto &rand_gen = singleton_c<rand_gen_c>::instance();
//     auto &simu_param = singleton_c<simu_param_c>::instance();
//     auto &demo_param = singleton_c<demo_param_c>::instance();
//     auto &samp_param = singleton_c<samp_param_c>::instance();
//     auto &recomb_param = singleton_c<recomb_param_c>::instance();
//     auto &muta_param = singleton_c<muta_param_c>::instance();
//     auto &info_collect = singleton_c<info_collector_c>::instance();

//     rand_gen.put_seed(72);

//     simu_param.Repetition_nbr = 10000000;
//     simu_param.Continuous_time_approxim = false;

//     info_collect.Stats = true;
//     info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

//     demo_param.Lattice_size = {1, 20};
//     demo_param.Disp_dist_max = {0, 19};
//     demo_param.Edge_effects = edge_effect_enum::circular;
//     demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
//     demo_param.Disp_func[1] = fwd_disp_distrib_c::uniform();
//     demo_param.Proba_migr = 0.00333;

//     demo_param.Pop_size_per_node = 30;
//     demo_param.Population_size_N = 600;

//     samp_param.Sample_size_per_node = std::vector<int>(20, 2);
//     samp_param.n_total_sample_size = 40;

//     samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
//     for (int x = 0; x < demo_param.Lattice_size.at(0); ++x)
//     {
//         for (int y = 0; y < demo_param.Lattice_size.at(1); ++y)
//         {
//             auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
//             samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
//         }
//     }

//     samp_param.Sequence_length = 1;
//     samp_param.Dist_class_nbr = 20;
//     recomb_param.Unscaled_recomb_rate = 0;

//     samp_param.Ploidy = 1;
//     muta_param.Mod_mut_name = mut_model_enum::iam;
//     muta_param.Unscaled_mut_rate_mu = 0.0001;
//     muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

//     output_stat_c output;
//     output.Stat_out = false;

//     for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
//     {
//         sample_simulator(output, rep);
//     }

//     auto m_Q0_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(0);
//     auto v_Q0_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(1);
//     auto Q0 = 0.927009974733469;

//     auto m_Q1_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(0);
//     auto v_Q1_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(1);
//     auto Q1 = 0.591377342954598;

//     std::cout << "m_Q0_seq : " << m_Q0_seq << std::endl;
//     std::cout << "v_Q0_seq : " << v_Q0_seq << std::endl;
//     std::cout << "*************" << std::endl;
//     std::cout << "m_Q1_seq : " << m_Q1_seq << std::endl;
//     std::cout << "v_Q1_seq : " << v_Q1_seq << std::endl;

//     REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
//     REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
//     REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));

//     REQUIRE(m_Q1_seq == Approx(Q1).margin(0.0001));
//     REQUIRE(Q1 < m_Q1_seq + 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
//     REQUIRE(Q1 > m_Q1_seq - 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
// }

// WARNING : Too much variance for testing in moderage time (less than 1 day)
// *************
// Qw_goal : 0.93135
// m_Qw_seq : 0.931488
// v_Qw_seq : 0.00271667
// *************
// Q0_goal : 0.929098
// m_Q0_seq : 0.929255
// v_Q0_seq : 0.00263081
// *************
// Q1_goal : 0.592709
// m_Q1_seq : 0.591267
// v_Q1_seq : 0.0559707
// *************
// TEST_CASE("Gen_per_gen__IAM__Island_model__20_pop__15_diplo_indivs_per_pop__0_0051_prob_migration__0_0001_mut_rate__1_site__10000000_rep__1/10000_precision")
// {
//     auto &rand_gen = singleton_c<rand_gen_c>::instance();
//     auto &simu_param = singleton_c<simu_param_c>::instance();
//     auto &demo_param = singleton_c<demo_param_c>::instance();
//     auto &samp_param = singleton_c<samp_param_c>::instance();
//     auto &recomb_param = singleton_c<recomb_param_c>::instance();
//     auto &muta_param = singleton_c<muta_param_c>::instance();
//     auto &info_collect = singleton_c<info_collector_c>::instance();

//     rand_gen.put_seed(72);

//     simu_param.Repetition_nbr = 10000000;
//     simu_param.Continuous_time_approxim = false;

//     info_collect.Stats = true;
//     info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

//     demo_param.Lattice_size = {1, 20};
//     demo_param.Disp_dist_max = {0, 19};
//     demo_param.Edge_effects = edge_effect_enum::circular;
//     demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
//     demo_param.Disp_func[1] = fwd_disp_distrib_c::uniform();
//     demo_param.Proba_migr = 0.00333;

//     demo_param.Pop_size_per_node = 15;
//     demo_param.Population_size_N = 300;

//     samp_param.Sample_size_per_node = std::vector<int>(20, 5);
//     samp_param.n_total_sample_size = 100;

//     samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
//     for (int x = 0; x < demo_param.Lattice_size.at(0); ++x)
//     {
//         for (int y = 0; y < demo_param.Lattice_size.at(1); ++y)
//         {
//             auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
//             samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
//         }
//     }

//     samp_param.Sequence_length = 1;
//     samp_param.Dist_class_nbr = 20;
//     recomb_param.Unscaled_recomb_rate = 0;

//     samp_param.Ploidy = 2;
//     muta_param.Mod_mut_name = mut_model_enum::iam;
//     muta_param.Unscaled_mut_rate_mu = 0.0001;
//     muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

//     output_stat_c output;
//     output.Stat_out = false;

//     for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
//     {
//         sample_simulator(output, rep);
//     }

//     auto m_Qw_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwi_cumul_m_v.at(0);
//     auto v_Qw_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwi_cumul_m_v.at(1);
//     auto Qw = 0.93135;

//     auto m_Q0_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(0);
//     auto v_Q0_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(1);
//     auto Q0 = 0.929098;

//     auto m_Q1_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(0);
//     auto v_Q1_seq = output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(1);
//     auto Q1 = 0.592709;

//     std::cout << "m_Qw_seq : " << m_Qw_seq << std::endl;
//     std::cout << "v_Qw_seq : " << v_Qw_seq << std::endl;
//     std::cout << "*************" << std::endl;
//     std::cout << "m_Q0_seq : " << m_Q0_seq << std::endl;
//     std::cout << "v_Q0_seq : " << v_Q0_seq << std::endl;
//     std::cout << "*************" << std::endl;
//     std::cout << "m_Q1_seq : " << m_Q1_seq << std::endl;
//     std::cout << "v_Q1_seq : " << v_Q1_seq << std::endl;

//     REQUIRE(m_Qw_seq == Approx(Qw).margin(0.0001));
//     REQUIRE(Qw < m_Qw_seq + 1.96 * (std::sqrt(v_Qw_seq / simu_param.Repetition_nbr)));
//     REQUIRE(Qw > m_Qw_seq - 1.96 * (std::sqrt(v_Qw_seq / simu_param.Repetition_nbr)));

//     REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
//     REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
//     REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));

//     REQUIRE(m_Q1_seq == Approx(Q1).margin(0.0001));
//     REQUIRE(Q1 < m_Q1_seq + 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
//     REQUIRE(Q1 > m_Q1_seq - 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
// }

/****************************/
/*            IBD           */
/****************************/

TEST_CASE("Gen_per_gen__IAM__IBD_stepping_stone__10_pop__2_haplo_indivs_per_pop__0_6666_prob_migration__0_0005_mut_rate__1_site__10000000_rep__1/10000_precision")
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();

    rand_gen.put_seed(3);

    simu_param.Repetition_nbr = 10000000;
    simu_param.Continuous_time_approxim = false;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {1, 10};
    demo_param.Disp_dist_max = {0, 1};
    demo_param.Edge_effects = edge_effect_enum::circular;
    demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
    demo_param.Disp_func[1] = fwd_disp_distrib_c::uniform();
    demo_param.Proba_migr = 0.6666;

    demo_param.Pop_size_per_node = 2;
    demo_param.Population_size_N = 20;

    samp_param.Sample_size_per_node = std::vector<int>(10, 2);
    samp_param.n_total_sample_size = 20;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x < demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y < demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 10;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 1;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.0005;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.98052;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));

    auto m_Q1_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[1].at(0);
    auto v_Q1_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[1].at(1);
    auto Q1 = 0.978344;

    REQUIRE(m_Q1_seq == Approx(Q1).margin(0.0001));
    REQUIRE(Q1 < m_Q1_seq + 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q1 > m_Q1_seq - 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));

    auto m_Q2_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[2].at(0);
    auto v_Q2_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[2].at(1);
    auto Q2 = 0.974677;

    REQUIRE(m_Q2_seq == Approx(Q2).margin(0.0001));
    REQUIRE(Q2 < m_Q2_seq + 1.96 * (std::sqrt(v_Q2_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q2 > m_Q2_seq - 1.96 * (std::sqrt(v_Q2_seq / simu_param.Repetition_nbr)));

    auto m_Q3_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[3].at(0);
    auto v_Q3_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[3].at(1);
    auto Q3 = 0.970656;

    REQUIRE(m_Q3_seq == Approx(Q3).margin(0.0001));
    REQUIRE(Q3 < m_Q3_seq + 1.96 * (std::sqrt(v_Q3_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q3 > m_Q3_seq - 1.96 * (std::sqrt(v_Q3_seq / simu_param.Repetition_nbr)));

    auto m_Q4_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[4].at(0);
    auto v_Q4_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[4].at(1);
    auto Q4 = 0.968581;

    REQUIRE(m_Q4_seq == Approx(Q4).margin(0.0001));
    REQUIRE(Q4 < m_Q4_seq + 1.96 * (std::sqrt(v_Q4_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q4 > m_Q4_seq - 1.96 * (std::sqrt(v_Q4_seq / simu_param.Repetition_nbr)));

    auto m_Q5_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[5].at(0);
    auto v_Q5_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[5].at(1);
    auto Q5 = 0.96782;

    REQUIRE(m_Q5_seq == Approx(Q5).margin(0.0001));
    REQUIRE(Q5 < m_Q5_seq + 1.96 * (std::sqrt(v_Q5_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q5 > m_Q5_seq - 1.96 * (std::sqrt(v_Q5_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__IAM__IBD_stepping_stone__100_pop__2_haplo_indivs_per_pop__0_6666_prob_migration__0_0005_mut_rate__1_site__1000000_rep__1/10000_precision")
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
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {10,10};
    demo_param.Disp_dist_max = {1, 1};
    demo_param.Edge_effects = edge_effect_enum::circular;
    demo_param.Disp_func[0] = fwd_disp_distrib_c::uniform();
    demo_param.Disp_func[1] = fwd_disp_distrib_c::uniform();
    demo_param.Proba_migr = 0.6666;

    demo_param.Pop_size_per_node = 2;
    demo_param.Population_size_N = 200;

    samp_param.Sample_size_per_node = std::vector<int>(100, 2);
    samp_param.n_total_sample_size = 200;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x < demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y < demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 20;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 1;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.0005;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.836510;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));

    auto m_Q1_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[1].at(0);
    auto v_Q1_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[1].at(1);
    auto Q1 = 0.831826;

    REQUIRE(m_Q1_seq == Approx(Q1).margin(0.0001));
    REQUIRE(Q1 < m_Q1_seq + 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q1 > m_Q1_seq - 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));

    auto m_Q2_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[2].at(0);
    auto v_Q2_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[2].at(1);
    auto Q2 = 0.824988;

    // Q not just compute by axes
    // TODO : Option pour calculer par axe
    // REQUIRE(m_Q2_seq == Approx(Q2).margin(0.0001));
    // REQUIRE(Q2 < m_Q2_seq + 1.96 * (std::sqrt(v_Q2_seq / simu_param.Repetition_nbr)));
    // REQUIRE(Q2 > m_Q2_seq - 1.96 * (std::sqrt(v_Q2_seq / simu_param.Repetition_nbr)));

    // auto m_Q3_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[3].at(0);
    // auto v_Q3_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[3].at(1);
    // auto Q3 = 0.817718;

    // REQUIRE(m_Q3_seq == Approx(Q3).margin(0.0001));
    // REQUIRE(Q3 < m_Q3_seq + 1.96 * (std::sqrt(v_Q3_seq / simu_param.Repetition_nbr)));
    // REQUIRE(Q3 > m_Q3_seq - 1.96 * (std::sqrt(v_Q3_seq / simu_param.Repetition_nbr)));

    // auto m_Q4_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[4].at(0);
    // auto v_Q4_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[4].at(1);
    // auto Q4 = 0.814846;

    // REQUIRE(m_Q4_seq == Approx(Q4).margin(0.0001));
    // REQUIRE(Q4 < m_Q4_seq + 1.96 * (std::sqrt(v_Q4_seq / simu_param.Repetition_nbr)));
    // REQUIRE(Q4 > m_Q4_seq - 1.96 * (std::sqrt(v_Q4_seq / simu_param.Repetition_nbr)));

    // auto m_Q5_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[5].at(0);
    // auto v_Q5_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[5].at(1);
    // auto Q5 = 0.813805;

    // REQUIRE(m_Q5_seq == Approx(Q5).margin(0.0001));
    // REQUIRE(Q5 < m_Q5_seq + 1.96 * (std::sqrt(v_Q5_seq / simu_param.Repetition_nbr)));
    // REQUIRE(Q5 > m_Q5_seq - 1.96 * (std::sqrt(v_Q5_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__IAM__IBD_pareto_n_3_79809435_dist_max_49__20_pop__1_haplo_indivs_per_pop__0_6_prob_migration__0_0005_mut_rate__1_site__3000000_rep__1/10000_precision")
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
    simu_param.Continuous_time_approxim = false;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {1, 10};
    demo_param.Disp_dist_max = {0, 49};
    demo_param.Edge_effects = edge_effect_enum::circular;
    demo_param.Disp_func[0] = fwd_disp_distrib_c::pareto(3.79809435);
    demo_param.Disp_func[1] = fwd_disp_distrib_c::pareto(3.79809435);
    demo_param.Proba_migr = 0.6;

    demo_param.Pop_size_per_node = 2;
    demo_param.Population_size_N = 20;

    samp_param.Sample_size_per_node = std::vector<int>(10, 2);
    samp_param.n_total_sample_size = 20;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x < demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y < demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 10;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 1;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.0005;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.980492;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));

    auto m_Q1_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[1].at(0);
    auto v_Q1_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[1].at(1);
    auto Q1 = 0.978869;

    REQUIRE(m_Q1_seq == Approx(Q1).margin(0.0001));
    REQUIRE(Q1 < m_Q1_seq + 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q1 > m_Q1_seq - 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));

    auto m_Q2_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[2].at(0);
    auto v_Q2_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[2].at(1);
    auto Q2 = 0.975505;

    REQUIRE(m_Q2_seq == Approx(Q2).margin(0.0001));
    REQUIRE(Q2 < m_Q2_seq + 1.96 * (std::sqrt(v_Q2_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q2 > m_Q2_seq - 1.96 * (std::sqrt(v_Q2_seq / simu_param.Repetition_nbr)));

    auto m_Q3_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[3].at(0);
    auto v_Q3_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[3].at(1);
    auto Q3 = 0.972546;

    REQUIRE(m_Q3_seq == Approx(Q3).margin(0.0001));
    REQUIRE(Q3 < m_Q3_seq + 1.96 * (std::sqrt(v_Q3_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q3 > m_Q3_seq - 1.96 * (std::sqrt(v_Q3_seq / simu_param.Repetition_nbr)));

    auto m_Q4_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[4].at(0);
    auto v_Q4_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[4].at(1);
    auto Q4 = 0.9709;

    REQUIRE(m_Q4_seq == Approx(Q4).margin(0.0001));
    REQUIRE(Q4 < m_Q4_seq + 1.96 * (std::sqrt(v_Q4_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q4 > m_Q4_seq - 1.96 * (std::sqrt(v_Q4_seq / simu_param.Repetition_nbr)));

    auto m_Q5_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[5].at(0);
    auto v_Q5_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[5].at(1);
    auto Q5 = 0.970346;

    REQUIRE(m_Q5_seq == Approx(Q5).margin(0.0001));
    REQUIRE(Q5 < m_Q5_seq + 1.96 * (std::sqrt(v_Q5_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q5 > m_Q5_seq - 1.96 * (std::sqrt(v_Q5_seq / simu_param.Repetition_nbr)));
}

TEST_CASE("Gen_per_gen__IAM__IBD_pareto_n_2_74376568017753_dist_max_48__10_pop__2_haplo_indivs_per_pop__0_7_prob_migration__0_0005_mut_rate__1_site__3000000_rep__1/10000_precision")
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
    simu_param.Continuous_time_approxim = false;
    info_collect.Stats = true;
    info_collect.Prob_id_1_loc_Qr = true;
    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;

    demo_param.Lattice_size = {1, 10};
    demo_param.Disp_dist_max = {0, 48};
    demo_param.Edge_effects = edge_effect_enum::circular;
    demo_param.Disp_func[0] = fwd_disp_distrib_c::pareto(2.74376568017753);
    demo_param.Disp_func[1] = fwd_disp_distrib_c::pareto(2.74376568017753);
    demo_param.Proba_migr = 0.7;

    demo_param.Pop_size_per_node = 2;
    demo_param.Population_size_N = 20;

    samp_param.Sample_size_per_node = std::vector<int>(10, 2);
    samp_param.n_total_sample_size = 20;

    samp_param.Sample_coord_vec = std::vector<std::array<int, 2>>();
    for (int x = 0; x < demo_param.Lattice_size.at(0); ++x)
    {
        for (int y = 0; y < demo_param.Lattice_size.at(1); ++y)
        {
            auto coord01 = std::vector<std::array<int, 2>>(samp_param.Sample_size_per_node[0], std::array<int, 2>({x, y}));
            samp_param.Sample_coord_vec.insert(samp_param.Sample_coord_vec.end(), coord01.begin(), coord01.end());
        }
    }

    samp_param.Sequence_length = 1;
    samp_param.Dist_class_nbr = 10;
    recomb_param.Unscaled_recomb_rate = 0;

    samp_param.Ploidy = 1;
    muta_param.Mod_mut_name = mut_model_enum::iam;
    muta_param.Unscaled_mut_rate_mu = 0.0005;
    muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));

    output_stat_c output;
    output.Stat_out = false;

    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        sample_simulator(output, rep);
    }

    auto m_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0);
    auto v_Q0_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1);
    auto Q0 = 0.980435;

    REQUIRE(m_Q0_seq == Approx(Q0).margin(0.0001));
    REQUIRE(Q0 < m_Q0_seq + 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q0 > m_Q0_seq - 1.96 * (std::sqrt(v_Q0_seq / simu_param.Repetition_nbr)));

    auto m_Q1_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[1].at(0);
    auto v_Q1_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[1].at(1);
    auto Q1 = 0.979572;

    REQUIRE(m_Q1_seq == Approx(Q1).margin(0.0001));
    REQUIRE(Q1 < m_Q1_seq + 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q1 > m_Q1_seq - 1.96 * (std::sqrt(v_Q1_seq / simu_param.Repetition_nbr)));

    auto m_Q2_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[2].at(0);
    auto v_Q2_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[2].at(1);
    auto Q2 = 0.978036;

    REQUIRE(m_Q2_seq == Approx(Q2).margin(0.0001));
    REQUIRE(Q2 < m_Q2_seq + 1.96 * (std::sqrt(v_Q2_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q2 > m_Q2_seq - 1.96 * (std::sqrt(v_Q2_seq / simu_param.Repetition_nbr)));

    auto m_Q3_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[3].at(0);
    auto v_Q3_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[3].at(1);
    auto Q3 = 0.976451;

    REQUIRE(m_Q3_seq == Approx(Q3).margin(0.0001));
    REQUIRE(Q3 < m_Q3_seq + 1.96 * (std::sqrt(v_Q3_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q3 > m_Q3_seq - 1.96 * (std::sqrt(v_Q3_seq / simu_param.Repetition_nbr)));

    auto m_Q4_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[4].at(0);
    auto v_Q4_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[4].at(1);
    auto Q4 = 0.975642;

    REQUIRE(m_Q4_seq == Approx(Q4).margin(0.0001));
    REQUIRE(Q4 < m_Q4_seq + 1.96 * (std::sqrt(v_Q4_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q4 > m_Q4_seq - 1.96 * (std::sqrt(v_Q4_seq / simu_param.Repetition_nbr)));

    auto m_Q5_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[5].at(0);
    auto v_Q5_seq = output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[5].at(1);
    auto Q5 = 0.975374;

    REQUIRE(m_Q5_seq == Approx(Q5).margin(0.0001));
    REQUIRE(Q5 < m_Q5_seq + 1.96 * (std::sqrt(v_Q5_seq / simu_param.Repetition_nbr)));
    REQUIRE(Q5 > m_Q5_seq - 1.96 * (std::sqrt(v_Q5_seq / simu_param.Repetition_nbr)));
}