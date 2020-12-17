#define CATCH_CONFIG_MAIN ProbaIndentest
#include "catch.hpp"
#include <iostream>

#include "summary_stat.hpp"

TEST_CASE("summary_stat_prob_id_c_test")
{
    // SECTION("constructor")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     REQUIRE(proba.Ancestry_seq.size() == number_of_site);
    // }

    // SECTION("recombnstitution_full_seq_indiv")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {1, 0};
    //     indiv1[1] = {3, 2};

    //     auto full_seq_indiv1 = reconstitution_full_int_seq_indiv(ancestry_seq, indiv1);

    //     REQUIRE(full_seq_indiv1[0] == 1);
    //     REQUIRE(full_seq_indiv1[1] == 0);
    //     REQUIRE(full_seq_indiv1[2] == 1);
    //     REQUIRE(full_seq_indiv1[3] == 2);
    //     REQUIRE(full_seq_indiv1[4] == 1);
    // }

    // SECTION("phiA_comparison_2_indiv_1_site_all_different_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {1, 0};
    //     indiv1[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv2(3);
    //     indiv2[0] = {0, 2};
    //     indiv2[1] = {2, 0};
    //     indiv2[2] = {4, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum);

    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[0] == 1);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[1] == 1);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[2] == 1);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[3] == 1);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[4] == 1);
    //     REQUIRE(proba.Nbr_dif_for_phiA_seq == 1);
    // }

    // SECTION("phiA_comparison_2_indiv_1_site_all_same_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {1, 0};
    //     indiv1[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {3, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum);

    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[0] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[1] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[2] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[3] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[4] == 0);
    //     REQUIRE(proba.Nbr_dif_for_phiA_seq == 0);
    // }

    // SECTION("phiA_comparison_2_indiv_1_site_rand_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(4);
    //     indiv1[0] = {1, 1};
    //     indiv1[1] = {2, 2};
    //     indiv1[2] = {3, 2};
    //     indiv1[3] = {4, 1};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 1};
    //     indiv2[1] = {3, 1};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum);

    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[0] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[1] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[2] == 1);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[3] == 1);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[4] == 0);
    //     REQUIRE(proba.Nbr_dif_for_phiA_seq == 1);
    // }

    // SECTION("phiA_comparison_3_indiv_1_site_rand_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(4);
    //     indiv1[0] = {1, 1};
    //     indiv1[1] = {2, 2};
    //     indiv1[2] = {3, 2};
    //     indiv1[3] = {4, 1};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 1};
    //     indiv2[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {3, 1};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum);

    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[0] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[1] == 2);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[2] == 3);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[3] == 2);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiA[4] == 0);
    //     REQUIRE(proba.Total_comparison_phiAB == 3);
    //     REQUIRE(proba.Nbr_dif_for_phiA_seq == 3);
    // }

    // SECTION("phiA_moy")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(4);
    //     indiv1[0] = {1, 1};
    //     indiv1[1] = {2, 2};
    //     indiv1[2] = {3, 2};
    //     indiv1[3] = {4, 1};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 1};
    //     indiv2[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {3, 1};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum);
    //     std::array<double, 2> proba_moy = calc_mean_var(proba.PhiA);

    //     REQUIRE(std::get<0>(proba_moy) == Approx(1 - 0.46).margin(0.01));
    //     REQUIRE(std::get<1>(proba_moy) == Approx(0.16).margin(0.01));
    // }

    // /************************************/

    // SECTION("comparison_one_indiv_one_indiv_at_one_site")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     Prob_id_1_loc_Qr_res_c proba;

    //     std::vector<std::array<int, 2>> indiv1(3);
    //     indiv1[0] = {0, 2};
    //     indiv1[1] = {3, 2};
    //     indiv1[2] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv2(4);
    //     indiv2[0] = {0, 0};
    //     indiv2[1] = {1, 2};
    //     indiv2[2] = {2, 2};
    //     indiv2[3] = {4, 0};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv2};

    //     std::array<int, 2> value_comparison = proba.comparison_one_indiv_one_group_at_one_site(indiv1, resum, ancestry_seq);

    //     REQUIRE(value_comparison[0] == 1);
    //     REQUIRE(value_comparison[1] == 5);
    // }

    // SECTION("comparison_one_indiv_one_group_at_one_site")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     Prob_id_1_loc_Qr_res_c proba;

    //     std::vector<std::array<int, 2>> indiv1(3);
    //     indiv1[0] = {0, 2};
    //     indiv1[1] = {3, 2};
    //     indiv1[2] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv2(4);
    //     indiv2[0] = {0, 0};
    //     indiv2[1] = {1, 2};
    //     indiv2[2] = {2, 2};
    //     indiv2[3] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv3(4);
    //     indiv3[0] = {0, 0};
    //     indiv3[1] = {1, 2};
    //     indiv3[2] = {2, 2};
    //     indiv3[3] = {4, 2};

    //     std::vector<std::array<int, 2>> indiv4(3);
    //     indiv4[0] = {0, 0};
    //     indiv4[1] = {3, 2};
    //     indiv4[2] = {4, 0};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv2, indiv3, indiv4};

    //     std::array<int, 2> value_comparison = proba.comparison_one_indiv_one_group_at_one_site(indiv1, resum, ancestry_seq);

    //     REQUIRE(value_comparison[0] == 5);
    //     REQUIRE(value_comparison[1] == 15);
    // }

    // SECTION("comparison_two_group_at_one_site(one_indiv_vs_one_group)")
    // {
    //     auto &muta_param = singleton_c<muta_param_c>::instance();
    //     int number_of_site = 5;
    //     muta_param.Ancestry_seq = std::vector<int>(number_of_site, 1);

    //     Prob_id_1_loc_Qr_res_c proba;

    //     std::vector<std::array<int, 2>> indiv1(3);
    //     indiv1[0] = {0, 2};
    //     indiv1[1] = {1, 2};
    //     indiv1[2] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv2(4);
    //     indiv2[0] = {0, 0};
    //     indiv2[1] = {1, 0};
    //     indiv2[2] = {3, 2};
    //     indiv2[3] = {4, 2};

    //     std::vector<std::array<int, 2>> indiv3(3);
    //     indiv3[0] = {0, 2};
    //     indiv3[1] = {2, 2};
    //     indiv3[2] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv4(2);
    //     indiv4[0] = {0, 0};
    //     indiv4[1] = {4, 2};

    //     std::vector<std::array<int, 2>> indiv5(2);
    //     indiv5[0] = {2, 2};
    //     indiv5[1] = {3, 0};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1};
    //     std::vector<std::vector<std::array<int, 2>>> resum1 = {indiv2, indiv3, indiv4, indiv5};

    //     std::array<int, 2> value_comparison = proba.comparison_two_group_at_one_site({resum, resum1}, 0, 1);

    //     REQUIRE(value_comparison[0] == 4);
    //     REQUIRE(value_comparison[1] == 20);
    // }

    // SECTION("comparison_one_group_at_one_site")
    // {
    //     auto &muta_param = singleton_c<muta_param_c>::instance();
    //     int number_of_site = 5;
    //     muta_param.Ancestry_seq = std::vector<int>(number_of_site, 1);

    //     Prob_id_1_loc_Qr_res_c proba;

    //     std::vector<std::array<int, 2>> indiv1(3);
    //     indiv1[0] = {0, 2};
    //     indiv1[1] = {1, 2};
    //     indiv1[2] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv2(4);
    //     indiv2[0] = {0, 0};
    //     indiv2[1] = {1, 0};
    //     indiv2[2] = {3, 2};
    //     indiv2[3] = {4, 2};

    //     std::vector<std::array<int, 2>> indiv3(3);
    //     indiv3[0] = {0, 2};
    //     indiv3[1] = {2, 2};
    //     indiv3[2] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv4(2);
    //     indiv4[0] = {0, 0};
    //     indiv4[1] = {4, 2};

    //     std::vector<std::array<int, 2>> indiv5(2);
    //     indiv5[0] = {2, 2};
    //     indiv5[1] = {3, 0};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3, indiv4, indiv5};

    //     std::array<int, 2> value_comparison = proba.comparison_two_group_at_one_site({resum}, 0, 0);

    //     REQUIRE(value_comparison[0] == 12);
    //     REQUIRE(value_comparison[1] == 50);
    // }

    // SECTION("comparison_two_group_at_one_site")
    // {
    //     auto &muta_param = singleton_c<muta_param_c>::instance();
    //     int number_of_site = 5;
    //     muta_param.Ancestry_seq = std::vector<int>(number_of_site, 1);

    //     Prob_id_1_loc_Qr_res_c proba;

    //     std::vector<std::array<int, 2>> indiv1(3);
    //     indiv1[0] = {0, 2};
    //     indiv1[1] = {1, 2};
    //     indiv1[2] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv2(4);
    //     indiv2[0] = {0, 0};
    //     indiv2[1] = {1, 0};
    //     indiv2[2] = {3, 2};
    //     indiv2[3] = {4, 2};

    //     std::vector<std::array<int, 2>> indiv3(3);
    //     indiv3[0] = {0, 2};
    //     indiv3[1] = {2, 2};
    //     indiv3[2] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv4(2);
    //     indiv4[0] = {0, 0};
    //     indiv4[1] = {4, 2};

    //     std::vector<std::array<int, 2>> indiv5(2);
    //     indiv5[0] = {2, 2};
    //     indiv5[1] = {3, 0};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2};
    //     std::vector<std::vector<std::array<int, 2>>> resum1 = {indiv3, indiv4, indiv5};

    //     std::array<int, 2> value_comparison = proba.comparison_two_group_at_one_site({resum, resum1}, 0, 1);

    //     REQUIRE(value_comparison[0] == 7);
    //     REQUIRE(value_comparison[1] == 30);
    // }

    // SECTION("comparison_node_sample_2_by_2_at_1_and_2_dist")
    // {
    //     auto &muta_param = singleton_c<muta_param_c>::instance();
    //     int number_of_site = 6;
    //     muta_param.Ancestry_seq = std::vector<int>(number_of_site, 1);

    //     Prob_id_1_loc_Qr_res_c proba;

    //     std::vector<std::array<int, 2>> indiv1(4);
    //     indiv1[0] = {0, 2};
    //     indiv1[1] = {2, 2};
    //     indiv1[2] = {3, 0};
    //     indiv1[3] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv2(5);
    //     indiv2[0] = {0, 2};
    //     indiv2[1] = {1, 2};
    //     indiv2[2] = {2, 2};
    //     indiv2[3] = {4, 2};
    //     indiv2[4] = {5, 0};

    //     std::vector<std::array<int, 2>> indiv3(3);
    //     indiv3[0] = {0, 2};
    //     indiv3[1] = {1, 0};
    //     indiv3[2] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv4(2);
    //     indiv4[0] = {2, 2};
    //     indiv4[1] = {5, 2};

    //     std::vector<std::array<int, 2>> indiv5(3);
    //     indiv5[0] = {1, 2};
    //     indiv5[1] = {2, 2};
    //     indiv5[2] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv6(3);
    //     indiv6[0] = {1, 0};
    //     indiv6[1] = {3, 2};
    //     indiv6[2] = {5, 2};

    //     std::vector<std::vector<std::vector<std::array<int, 2>>>> resum = {{indiv1, indiv2, indiv3}, {indiv4, indiv5, indiv6}};

    //     std::array<int, 2> value_comparison_Q0 = proba.comparison_node_sample_2_by_2(resum, {{0, 0}, {1, 1}});

    //     REQUIRE(value_comparison_Q0[0] == 15);
    //     REQUIRE(value_comparison_Q0[1] == 36);

    //     std::array<int, 2> value_comparison_Q1 = proba.comparison_node_sample_2_by_2(resum, {{0, 1}});

    //     REQUIRE(value_comparison_Q1[0] == 14);
    //     REQUIRE(value_comparison_Q1[1] == 54);
    // }
    // /************************************/

    // SECTION("phiAB_comparison_2_indiv_2_site_dist_1_all_different_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);
    //     ;

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {1, 0};
    //     indiv1[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv2(3);
    //     indiv2[0] = {0, 2};
    //     indiv2[1] = {2, 0};
    //     indiv2[2] = {4, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 1);

    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[1][0] == 1);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[1][1] == 1);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[1][2] == 1);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[1][3] == 1);
    // }

    // SECTION("phiAB_comparison_2_indiv_2_site_dist_1_all_same_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {1, 0};
    //     indiv1[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {3, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 1);

    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[1][0] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[1][1] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[1][2] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[1][3] == 0);
    // }

    // SECTION("phiAB_comparison_2_indiv_2_site_dist_1_all_rand_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(4);
    //     indiv1[0] = {1, 1};
    //     indiv1[1] = {2, 2};
    //     indiv1[2] = {3, 2};
    //     indiv1[3] = {4, 1};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 1};
    //     indiv2[1] = {3, 1};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 1);

    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[1][0] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[1][1] == 1);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[1][2] == 1);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[1][3] == 1);
    // }

    // SECTION("phiAB_comparison_3_indiv_2_site_dist_1_all_different_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(4);
    //     indiv1[0] = {1, 1};
    //     indiv1[1] = {2, 2};
    //     indiv1[2] = {3, 2};
    //     indiv1[3] = {4, 1};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 1};
    //     indiv2[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {3, 1};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 1);

    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[1][0] == 2);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[1][1] == 3);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[1][2] == 3);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[1][3] == 2);
    //     REQUIRE(proba.Total_comparison_phiAB == 3);
    // }

    // SECTION("phiAB_dist_1_moy")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(4);
    //     indiv1[0] = {1, 1};
    //     indiv1[1] = {2, 2};
    //     indiv1[2] = {3, 2};
    //     indiv1[3] = {4, 1};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 1};
    //     indiv2[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {3, 1};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 1);
    //     std::array<double, 2> proba_moy = calc_mean_var(proba.PhiAB[1]);

    //     //tolerance in percent
    //     REQUIRE(std::get<0>(proba_moy) == Approx(1 - 0.83).margin(0.01));
    //     REQUIRE(std::get<1>(proba_moy) == Approx(0.0278).margin(0.0001));
    // }

    // /************************************/

    // SECTION("phiAB_comparison_2_indiv_2_site_dist_3_all_different_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);
    //     ;

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {1, 0};
    //     indiv1[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv2(3);
    //     indiv2[0] = {0, 2};
    //     indiv2[1] = {2, 0};
    //     indiv2[2] = {4, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 3);

    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[3][0] == 1);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[3][1] == 1);
    // }

    // SECTION("phiAB_comparison_2_indiv_2_site_dist_3_all_same_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {1, 0};
    //     indiv1[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {3, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 3);

    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[3][0] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[3][1] == 0);
    // }

    // SECTION("phiAB_comparison_2_indiv_2_site_dist_3_all_rand_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(4);
    //     indiv1[0] = {1, 1};
    //     indiv1[1] = {2, 2};
    //     indiv1[2] = {3, 2};
    //     indiv1[3] = {4, 1};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 1};
    //     indiv2[1] = {3, 1};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 3);

    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[3][0] == 1);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[3][1] == 0);
    // }

    // SECTION("phiAB_comparison_3_indiv_2_site_dist_3_all_different_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(4);
    //     indiv1[0] = {1, 1};
    //     indiv1[1] = {2, 2};
    //     indiv1[2] = {3, 2};
    //     indiv1[3] = {4, 1};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 1};
    //     indiv2[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {3, 1};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 3);

    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[3][0] == 2);
    //     REQUIRE(proba.Nbr_dif_per_site_for_phiAB[3][1] == 2);
    // }

    // SECTION("phiAB_dist_3_moy")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(4);
    //     indiv1[0] = {1, 1};
    //     indiv1[1] = {2, 2};
    //     indiv1[2] = {3, 2};
    //     indiv1[3] = {4, 1};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 1};
    //     indiv2[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {3, 1};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 3);
    //     std::array<double, 2> proba_moy = calc_mean_var(proba.PhiAB[3]);

    //     REQUIRE(std::get<0>(proba_moy) == Approx(1 - 0.66).margin(0.01));
    //     REQUIRE(std::get<1>(proba_moy) == Approx(0).margin(0.0001));
    // }

    // /************************************/

    // SECTION("gammaAB_comparison_3_indiv_2_site_dist_1_all_different_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(3);
    //     indiv1[0] = {0, 0};
    //     indiv1[1] = {2, 2};
    //     indiv1[2] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv3(3);
    //     indiv3[0] = {0, 2};
    //     indiv3[1] = {2, 0};
    //     indiv3[2] = {4, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 1);

    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][0] == 6);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][1] == 6);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][2] == 6);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][3] == 6);
    // }

    // SECTION("gammaAB_comparison_3_indiv_2_site_dist_1_all_same_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(3);
    //     indiv1[0] = {0, 0};
    //     indiv1[1] = {1, 1};
    //     indiv1[2] = {4, 2};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {0, 0};
    //     indiv2[1] = {4, 2};

    //     std::vector<std::array<int, 2>> indiv3(3);
    //     indiv3[0] = {0, 0};
    //     indiv3[1] = {1, 1};
    //     indiv3[2] = {4, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 1);

    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][0] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][1] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][2] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][3] == 0);
    // }

    // SECTION("gammaAB_comparison_3_indiv_2_site_dist_1_all_rand_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(4);
    //     indiv1[0] = {1, 1};
    //     indiv1[1] = {2, 2};
    //     indiv1[2] = {3, 2};
    //     indiv1[3] = {4, 1};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 1};
    //     indiv2[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {3, 1};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 1);

    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][0] == 4);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][1] == 6);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][2] == 6);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][3] == 4);
    // }

    // SECTION("gammaAB_comparison_3_indiv_2_site_dist_3_all_different_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(3);
    //     indiv1[0] = {0, 0};
    //     indiv1[1] = {2, 2};
    //     indiv1[2] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv3(3);
    //     indiv3[0] = {0, 2};
    //     indiv3[1] = {2, 0};
    //     indiv3[2] = {4, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 3);

    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[3][0] == 6);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[3][1] == 6);
    // }

    // SECTION("gammaAB_comparison_3_indiv_2_site_dist_3_all_same_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(3);
    //     indiv1[0] = {0, 0};
    //     indiv1[1] = {1, 1};
    //     indiv1[2] = {4, 2};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {0, 0};
    //     indiv2[1] = {4, 2};

    //     std::vector<std::array<int, 2>> indiv3(3);
    //     indiv3[0] = {0, 0};
    //     indiv3[1] = {1, 1};
    //     indiv3[2] = {4, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 3);

    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][0] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][1] == 0);
    // }

    // SECTION("gammaAB_comparison_3_indiv_2_site_dist_3_all_rand_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(4);
    //     indiv1[0] = {1, 1};
    //     indiv1[1] = {2, 2};
    //     indiv1[2] = {3, 2};
    //     indiv1[3] = {4, 1};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 1};
    //     indiv2[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {3, 1};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 3);

    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[3][0] == 4);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[3][1] == 4);
    // }

    // SECTION("gammaAB_comparison_4_indiv_2_site_dist_1_all_rand_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(4);
    //     indiv1[0] = {1, 1};
    //     indiv1[1] = {2, 2};
    //     indiv1[2] = {3, 2};
    //     indiv1[3] = {4, 1};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {2, 1};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {3, 0};

    //     std::vector<std::array<int, 2>> indiv4(4);
    //     indiv4[0] = {0, 0};
    //     indiv4[1] = {1, 2};
    //     indiv4[2] = {2, 2};
    //     indiv4[3] = {4, 0};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3, indiv4};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 1);

    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][0] == 22);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][1] == 24);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][2] == 22);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][3] == 22);
    // }

    // SECTION("gammaAB_comparison_4_indiv_2_site_dist_1_all_rand_mut_bis")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {2, 2};
    //     indiv1[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv4(3);
    //     indiv4[0] = {1, 0};
    //     indiv4[1] = {2, 2};
    //     indiv4[2] = {3, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3, indiv4};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 1);

    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][0] == 12);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][1] == 22);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][2] == 24);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[1][3] == 20);
    //     REQUIRE(proba.Nbr_dif_for_phiA_seq == 6);
    // }

    // SECTION("gammaAB_comparison_4_indiv_2_site_dist_3_all_rand_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {2, 2};
    //     indiv1[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv4(3);
    //     indiv4[0] = {1, 0};
    //     indiv4[1] = {2, 2};
    //     indiv4[2] = {3, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3, indiv4};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 3);

    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[3][0] == 16);
    //     REQUIRE(proba.Nbr_dif_per_site_for_gammaAB[3][1] == 18);
    //     REQUIRE(proba.Total_comparison_gammaAB == 24);
    // }

    // SECTION("gammaAB_dist_1_moy")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {2, 2};
    //     indiv1[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv4(3);
    //     indiv4[0] = {1, 0};
    //     indiv4[1] = {2, 2};
    //     indiv4[2] = {3, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3, indiv4};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 1);
    //     std::array<double, 2> proba_moy = calc_mean_var(proba.GammaAB[1]);

    //     REQUIRE(std::get<0>(proba_moy) == Approx(1 - 0.8125).margin(0.0001));
    //     REQUIRE(std::get<1>(proba_moy) == Approx(0.036).margin(0.001));
    // }

    // SECTION("gammaAB_dist_3_moy")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {2, 2};
    //     indiv1[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv4(3);
    //     indiv4[0] = {1, 0};
    //     indiv4[1] = {2, 2};
    //     indiv4[2] = {3, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3, indiv4};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 3);
    //     std::array<double, 2> proba_moy = calc_mean_var(proba.GammaAB[3]);

    //     REQUIRE(std::get<0>(proba_moy) == Approx(1 - 0.708).margin(0.001));
    //     REQUIRE(std::get<1>(proba_moy) == Approx(0.0017).margin(0.0001));
    // }

    // /*********************************/

    // SECTION("deltaAB_comparison_4_indiv_2_site_dist_1_all_different_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(4);
    //     indiv1[0] = {0, 2};
    //     indiv1[1] = {1, 0};
    //     indiv1[2] = {3, 2};
    //     indiv1[3] = {4, 2};

    //     std::vector<std::array<int, 2>> indiv2(4);
    //     indiv2[0] = {0, 0};
    //     indiv2[1] = {1, 2};
    //     indiv2[2] = {2, 2};
    //     indiv2[3] = {3, 0};

    //     std::vector<std::array<int, 2>> indiv3(1);
    //     indiv3[0] = {4, 2};

    //     std::vector<std::array<int, 2>> indiv4(2);
    //     indiv4[0] = {0, 0};
    //     indiv4[1] = {2, 0};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3, indiv4};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 1);

    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[1][0] == 6);
    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[1][1] == 6);
    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[1][2] == 6);
    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[1][3] == 6);
    // }

    // SECTION("deltaAB_comparison_4_indiv_2_site_dist_1_all_same_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(4);
    //     std::vector<std::array<int, 2>> indiv2(4);
    //     std::vector<std::array<int, 2>> indiv3(1);
    //     std::vector<std::array<int, 2>> indiv4(2);

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3, indiv4};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 1);

    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[1][0] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[1][1] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[1][2] == 0);
    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[1][3] == 0);
    // }

    // SECTION("deltaAB_comparison_4_indiv_2_site_dist_1_all_rand_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {2, 2};
    //     indiv1[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv4(3);
    //     indiv4[0] = {1, 0};
    //     indiv4[1] = {2, 2};
    //     indiv4[2] = {3, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3, indiv4};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 1);

    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[1][0] == 3);
    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[1][1] == 5);
    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[1][2] == 5);
    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[1][3] == 5);
    // }

    // SECTION("deltaAB_comparison_4_indiv_2_site_dist_3_all_different_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(4);
    //     indiv1[0] = {0, 2};
    //     indiv1[1] = {1, 0};
    //     indiv1[2] = {3, 2};
    //     indiv1[3] = {4, 2};

    //     std::vector<std::array<int, 2>> indiv2(4);
    //     indiv2[0] = {0, 0};
    //     indiv2[1] = {1, 2};
    //     indiv2[2] = {2, 2};
    //     indiv2[3] = {3, 0};

    //     std::vector<std::array<int, 2>> indiv3(1);
    //     indiv3[0] = {4, 2};

    //     std::vector<std::array<int, 2>> indiv4(2);
    //     indiv4[0] = {0, 0};
    //     indiv4[1] = {2, 0};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3, indiv4};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 3);

    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[3][0] == 6);
    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[3][1] == 6);
    // }

    // SECTION("deltaAB_comparison_4_indiv_2_site_dist_3_all_rand_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {2, 2};
    //     indiv1[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv4(3);
    //     indiv4[0] = {1, 0};
    //     indiv4[1] = {2, 2};
    //     indiv4[2] = {3, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3, indiv4};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 3);

    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[3][0] == 4);
    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[3][1] == 4);
    // }

    // SECTION("deltaAB_comparison_5_indiv_2_site_dist_1_all_rand_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {2, 2};
    //     indiv1[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv4(3);
    //     indiv4[0] = {1, 0};
    //     indiv4[1] = {2, 2};
    //     indiv4[2] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv5(2);
    //     indiv5[0] = {0, 0};
    //     indiv5[1] = {2, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3, indiv4, indiv5};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 1);

    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[1][0] == 24);
    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[1][1] == 25);
    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[1][2] == 25);
    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[1][3] == 22);
    // }

    // SECTION("deltaAB_comparison_5_indiv_2_site_dist_3_all_rand_mut")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {2, 2};
    //     indiv1[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv4(3);
    //     indiv4[0] = {1, 0};
    //     indiv4[1] = {2, 2};
    //     indiv4[2] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv5(2);
    //     indiv5[0] = {0, 0};
    //     indiv5[1] = {2, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3, indiv4, indiv5};
    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 3);

    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[3][0] == 22);
    //     REQUIRE(proba.Nbr_dif_per_site_for_deltaAB[3][1] == 22);
    //     REQUIRE(proba.Total_comparison_deltaAB == 30);
    //     REQUIRE(proba.Nbr_dif_for_phiA_seq == 10);
    // }

    // SECTION("deltaAB_dist_1_moy")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {2, 2};
    //     indiv1[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv4(3);
    //     indiv4[0] = {1, 0};
    //     indiv4[1] = {2, 2};
    //     indiv4[2] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv5(2);
    //     indiv5[0] = {0, 0};
    //     indiv5[1] = {2, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3, indiv4, indiv5};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 1);
    //     std::array<double, 2> proba_moy = calc_mean_var(proba.DeltaAB[1]);

    //     REQUIRE(std::get<0>(proba_moy) == Approx(1 - 0.8).margin(0.1));
    //     REQUIRE(std::get<1>(proba_moy) == Approx(0.0016).margin(0.0001));
    // }

    // SECTION("deltaAB_dist_3_moy")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);

    //     prob_id_c proba(ancestry_seq);

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {2, 2};
    //     indiv1[1] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 0};
    //     indiv3[1] = {2, 0};

    //     std::vector<std::array<int, 2>> indiv4(3);
    //     indiv4[0] = {1, 0};
    //     indiv4[1] = {2, 2};
    //     indiv4[2] = {3, 2};

    //     std::vector<std::array<int, 2>> indiv5(2);
    //     indiv5[0] = {0, 0};
    //     indiv5[1] = {2, 2};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3, indiv4, indiv5};

    //     proba.comparison_indivs_in_one_group_at_two_sites(resum, 3);
    //     std::array<double, 2> proba_moy = calc_mean_var(proba.DeltaAB[3]);

    //     REQUIRE(std::get<0>(proba_moy) == Approx(1 - 0.73).margin(0.01));
    //     REQUIRE(std::get<1>(proba_moy) == Approx(0).margin(0.0001));
    // }
}
