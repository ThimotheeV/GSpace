#define CATCH_CONFIG_MAIN SettingsTest
#include "catch.hpp"
#include <iostream>

#include "settings.hpp"
#include "mutation.hpp"
#include "output.hpp"

TEST_CASE("settings_test")
{
    SECTION("singleton_instance")
    {
        {
            auto &simu_param = singleton_c<simu_param_c>::instance();
            simu_param.Repetition_nbr = 1;
        }
        auto const &sim_p = singleton_c<simu_param_c>::instance();
        REQUIRE(sim_p.Repetition_nbr == 1);
    }

    SECTION("file_reading")
    {
        const auto res_str = read_file("testfile.ini");
        const auto vecres = slice_unix_windows_file_by_line(res_str);

        REQUIRE(vecres.size() == 3);
        REQUIRE(vecres[0] == "One sTrinG=machin");
        REQUIRE(vecres[1] == "one inT =   2  ");
        REQUIRE(vecres[2] == "one vecT_flot    = 4.5 , 5.0 ,   6");

        const auto pairvec = slice_by_char(vecres[2], '=');
        REQUIRE(pairvec[0] == "one vecT_flot    ");
        const auto clean_pair_vec = remove_spaces_tab_underscores(pairvec[1]);
        REQUIRE(double_vector_parse_by_comma(clean_pair_vec) == std::vector<double>{4.5, 5.0, 6.0});
    }

    SECTION("parse_file ")
    {
        const auto file_str = read_file("IbdSettings.txt");
        parser_str(file_str);
        check_param();
        apply_param();
        auto const &simu_param = singleton_c<simu_param_c>::instance();
        auto const &demo_param = singleton_c<demo_param_c>::instance();
        auto const &samp_param = singleton_c<samp_param_c>::instance();
        auto const &recomb_param = singleton_c<recomb_param_c>::instance();
        auto const &muta_param = singleton_c<muta_param_c>::instance();

        REQUIRE(simu_param.Repetition_nbr == 555);
        REQUIRE(simu_param.Continuous_time_approxim == false);
        REQUIRE(simu_param.Setting_filename == "GSpaceSettings.txt");
        REQUIRE(simu_param.Output_dir == ".");

        std::array<int, 2> temp = {2, 2};
        REQUIRE(demo_param.Lattice_size == temp);
        temp = {1, 1};
        REQUIRE(demo_param.Disp_dist_max == temp);
        REQUIRE(demo_param.Edge_effects == edge_effect_enum::absorbing);
        REQUIRE(demo_param.Dispersal_distrib == dispersal_distrib_enum::geometric);
        REQUIRE(demo_param.Proba_migr == 0.2);
        REQUIRE(demo_param.G_geo_param == 0.2);
        REQUIRE(demo_param.Population_size_N == 90);
        REQUIRE(demo_param.Pop_size_per_node == 10);

        REQUIRE(samp_param.n_total_sample_size == 12);
        REQUIRE(samp_param.Sample_size_per_node[0] == 3);

        REQUIRE(recomb_param.Unscaled_recomb_rate == 0.01);
        REQUIRE(samp_param.Chr_nbr == 2);
        REQUIRE(samp_param.Sequence_length == 2);
        REQUIRE(samp_param.Ploidy == 2);

        REQUIRE(muta_param.Mod_mut_name == mut_model_enum::kam);
        REQUIRE(muta_param.Unscaled_mut_rate_mu == 0.00001);
        //Random nucl
        std::vector<std::vector<int>> temp3 = {{5, 6}, {7, 4}};
        REQUIRE(muta_param.Ancestry_seq == temp3);
        REQUIRE(muta_param.K_min == 4);
        REQUIRE(muta_param.K_max == 10);
        REQUIRE(muta_param.P_gsm == 0.38);
        auto temp4 = std::make_tuple(0.7, 1.5);
        REQUIRE(std::get<0>(muta_param.Ratio_transi_transver) == std::get<0>(temp4));
        REQUIRE(std::get<1>(muta_param.Ratio_transi_transver) == std::get<1>(temp4));

        REQUIRE(muta_param.Equi_base_freq.A == 0.1);
        REQUIRE(muta_param.Equi_base_freq.G == 0.3);
        REQUIRE(muta_param.Equi_base_freq.C == 0.4);
        REQUIRE(muta_param.Equi_base_freq.T == 0.2);
    }

    SECTION("sampling_organisation")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        demo_param.Nbr_node_sampled_x = 1;
        demo_param.Nbr_node_sampled_y = 1;
        demo_param.Min_sample_coord_x = 0;
        demo_param.Min_sample_coord_y = 0;

        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_size_per_node = {100};

        compute_sample_coord(samp_param, demo_param);

        REQUIRE(samp_param.n_total_sample_size == 100);

        samp_param.Sample_coord_vec.clear();

        demo_param.Nbr_node_sampled_x = 2;
        demo_param.Nbr_node_sampled_y = 2;
        demo_param.Min_sample_coord_x = 1;
        demo_param.Min_sample_coord_y = 1;
        samp_param.Sample_size_per_node = {3};

        compute_sample_coord(samp_param, demo_param);

        std::vector<std::array<int, 2>> result = {{1, 1}, {1, 1}, {1, 1}, {1, 2}, {1, 2}, {1, 2}, {2, 1}, {2, 1}, {2, 1}, {2, 2}, {2, 2}, {2, 2}};

        REQUIRE(samp_param.Sample_coord_vec.size() == result.size());
        REQUIRE(samp_param.Sample_coord_vec == result);
    }

    SECTION("sampling_organisation_with_sample_coord")
    {
        auto &demo_param = singleton_c<demo_param_c>::instance();
        auto &samp_param = singleton_c<samp_param_c>::instance();
        auto &muta_param = singleton_c<muta_param_c>::instance();

        muta_param.Ancestry_seq.clear();

        auto const &file_str = read_file("IbdSettings2.txt");
        parser_str(file_str);
        check_param();
        compute_sample_coord(samp_param, demo_param);
        std::vector<std::array<int, 2>> result = {{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}};

        REQUIRE(samp_param.Sample_coord_vec.size() == result.size());
        REQUIRE(samp_param.Sample_coord_vec == result);
    }

    SECTION("diagnostic_tables")
    {
        auto const &file_str = read_file("IbdSettings.txt");
        parser_str(file_str);

        auto const &info_collect = singleton_c<info_collector_c>::instance();

        REQUIRE(info_collect.Prob_id_1_loc_Qr == true);
        REQUIRE(info_collect.Prob_id_1_loc_Qwi_wd_bd == true);
        REQUIRE(info_collect.Prob_id_1_2_loc == true);
    }
}

TEST_CASE("cmd_line_settings_test")
{
    SECTION("cmd_reading")
    {
        int argc = 1;
        char gs[] = "./GSpace";
        char *argv[4] = {gs};
        std::string file_str = read_write_cmdline(argc, argv);

        REQUIRE(file_str.empty());
    }

    SECTION("cmd_reading")
    {
        int argc = 4;
        char gs[] = "./GSpace";
        char in[] = "Input=Settings.txt";
        char out[] = "Output=Trash";
        char opn[] = "OutputPointNbr=10";
        char *argv[4] = {gs, in, out, opn};
        std::string file_str = read_write_cmdline(argc, argv);

        REQUIRE(file_str == "Input=Settings.txt\nOutput=Trash\nOutputPointNbr=10");
    }

    SECTION("cmd_parsing")
    {
        int argc = 3;
        char gs[] = "./GSpace";
        char in[] = "Setting_filename=Settings.txt";
        char out[] = "Output_dir=Trash";
        char *argv[3] = {gs, in, out};
        std::string file_str = read_write_cmdline(argc, argv);
        parser_str(file_str);

        auto const &simu_param = singleton_c<simu_param_c>::instance();
        REQUIRE(simu_param.Setting_filename == "Settings.txt");
        REQUIRE(simu_param.Output_dir == "Trash");
    }
}

TEST_CASE("check_ancestry_settings_test")
{
    SECTION("iam special case")
    {
        const auto file_str = read_file("IbdSettings.txt");
        parser_str(file_str);

        auto &samp_param = singleton_c<samp_param_c>::instance();
        auto &muta_param = singleton_c<muta_param_c>::instance();

        samp_param.Chr_nbr = 1;
        samp_param.Sequence_length = 10;
        muta_param.Mod_mut_name = mut_model_enum::iam;
        muta_param.K_min = 0;
        muta_param.K_max = 0;

        muta_param.Ancestry_seq = {};

        check_param();

        std::vector<std::vector<int>> temp = {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1}};
        REQUIRE(muta_param.Ancestry_seq == temp);
    }

    SECTION("allelic 0 case")
    {
        const auto file_str = read_file("IbdSettings.txt");
        parser_str(file_str);

        auto &samp_param = singleton_c<samp_param_c>::instance();
        auto &muta_param = singleton_c<muta_param_c>::instance();

        samp_param.Chr_nbr = 2;
        samp_param.Sequence_length = 2;
        muta_param.Mod_mut_name = mut_model_enum::kam;
        muta_param.K_min = 52;
        muta_param.K_max = 52;

        muta_param.Ancestry_seq = {{0}};

        check_param();

        std::vector<std::vector<int>> temp = {{52, 52}, {52, 52}};
        REQUIRE(muta_param.Ancestry_seq == temp);
    }

    SECTION("allelic 0 and values case")
    {
        const auto file_str = read_file("IbdSettings.txt");
        parser_str(file_str);

        auto &samp_param = singleton_c<samp_param_c>::instance();
        auto &muta_param = singleton_c<muta_param_c>::instance();

        samp_param.Chr_nbr = 2;
        samp_param.Sequence_length = 2;
        muta_param.Mod_mut_name = mut_model_enum::kam;
        muta_param.K_min = 51;
        muta_param.K_max = 52;

        muta_param.Ancestry_seq = {{0}, {51, 51}};

        check_param();

        std::vector<std::vector<int>> temp = {{52, 52}, {51, 51}};
        REQUIRE(muta_param.Ancestry_seq == temp);
    }

    SECTION("allelic with 0 in values case")
    {
        const auto file_str = read_file("IbdSettings.txt");
        parser_str(file_str);

        auto &samp_param = singleton_c<samp_param_c>::instance();
        auto &muta_param = singleton_c<muta_param_c>::instance();

        samp_param.Chr_nbr = 2;
        samp_param.Sequence_length = 2;
        muta_param.Mod_mut_name = mut_model_enum::kam;
        muta_param.K_min = 51;
        muta_param.K_max = 52;

        muta_param.Ancestry_seq = {{0, 51}, {51, 51}};

        check_param();

        std::vector<std::vector<int>> temp = {{52, 51}, {51, 51}};
        REQUIRE(muta_param.Ancestry_seq == temp);
    }

    SECTION("allelic values case")
    {
        const auto file_str = read_file("IbdSettings.txt");
        parser_str(file_str);

        auto &samp_param = singleton_c<samp_param_c>::instance();
        auto &muta_param = singleton_c<muta_param_c>::instance();

        samp_param.Chr_nbr = 2;
        samp_param.Sequence_length = 2;
        muta_param.Mod_mut_name = mut_model_enum::kam;
        muta_param.K_min = 1;
        muta_param.K_max = 50;

        muta_param.Ancestry_seq = {{22, 23}, {11, 12}};

        check_param();

        std::vector<std::vector<int>> temp = {{22, 23}, {11, 12}};
        REQUIRE(muta_param.Ancestry_seq == temp);
    }

    SECTION("nucleotidic 0 case")
    {
        const auto file_str = read_file("IbdSettings.txt");
        parser_str(file_str);

        auto &samp_param = singleton_c<samp_param_c>::instance();
        auto &muta_param = singleton_c<muta_param_c>::instance();

        samp_param.Chr_nbr = 2;
        samp_param.Sequence_length = 2;
        muta_param.Mod_mut_name = mut_model_enum::jcm;

        muta_param.Ancestry_seq = {{0}};

        check_param();

        REQUIRE(muta_param.Ancestry_seq.size() == 2);
        REQUIRE(muta_param.Ancestry_seq.at(0).size() == 2);
        REQUIRE(muta_param.Ancestry_seq.at(1).size() == 2);
    }

    SECTION("nucleotidic 0 and values case")
    {
        const auto file_str = read_file("IbdSettings.txt");
        parser_str(file_str);

        auto &samp_param = singleton_c<samp_param_c>::instance();
        auto &muta_param = singleton_c<muta_param_c>::instance();

        samp_param.Chr_nbr = 2;
        samp_param.Sequence_length = 2;
        muta_param.Mod_mut_name = mut_model_enum::jcm;

        muta_param.Ancestry_seq = {{0}, {2, 3}};

        check_param();

        REQUIRE(muta_param.Ancestry_seq.size() == 2);
        REQUIRE(muta_param.Ancestry_seq.at(0).size() == 2);
        REQUIRE(muta_param.Ancestry_seq.at(1).size() == 2);
    }

    SECTION("nucleotidic with 0 in values case")
    {
        const auto file_str = read_file("IbdSettings.txt");
        parser_str(file_str);

        auto &samp_param = singleton_c<samp_param_c>::instance();
        auto &muta_param = singleton_c<muta_param_c>::instance();

        samp_param.Chr_nbr = 2;
        samp_param.Sequence_length = 2;
        muta_param.Mod_mut_name = mut_model_enum::jcm;

        muta_param.Ancestry_seq = {{0, 1}, {2, 3}};

        check_param();

        REQUIRE(muta_param.Ancestry_seq.size() == 2);
        REQUIRE(muta_param.Ancestry_seq.at(0).size() == 2);
        REQUIRE(muta_param.Ancestry_seq.at(0).at(0) != 0);
        REQUIRE(muta_param.Ancestry_seq.at(0).at(1) == 1);
        REQUIRE(muta_param.Ancestry_seq.at(1).size() == 2);
    }

    SECTION("nucleotidic values case")
    {
        const auto file_str = read_file("IbdSettings.txt");
        parser_str(file_str);

        auto &samp_param = singleton_c<samp_param_c>::instance();
        auto &muta_param = singleton_c<muta_param_c>::instance();

        samp_param.Chr_nbr = 2;
        samp_param.Sequence_length = 2;
        muta_param.Mod_mut_name = mut_model_enum::jcm;

        muta_param.Ancestry_seq = {{1, 4}, {1, 2}};

        check_param();

        REQUIRE(muta_param.Ancestry_seq.size() == 2);
        REQUIRE(muta_param.Ancestry_seq.at(0).size() == 2);
        REQUIRE(muta_param.Ancestry_seq.at(1).size() == 2);
    }
}
