#define CATCH_CONFIG_MAIN ProbaFrequencyAlleletest
#include "catch.hpp"
#include <iostream>

#include "summary_stat.hpp"

TEST_CASE("summary_stat_freq_allele_c_test")
{
    // SECTION("constructor")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);
    //     int sample_size = 4;

    //     freq_allele_c frequency_proba(ancestry_seq, sample_size);

    //     REQUIRE(frequency_proba.Ancestry_seq.size() == number_of_site);
    //     REQUIRE(frequency_proba.Sample_size == 4);
    // }

    // SECTION("metering_mut_for_frequency_simple_case")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);
    //     int sample_size = 3;

    //     freq_allele_c frequency_proba(ancestry_seq, sample_size);

    //     std::vector<std::array<int, 2>> indiv1(2);
    //     indiv1[0] = {0, 0};
    //     indiv1[1] = {3, 3};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv3(1);
    //     indiv3[0] = {1, 3};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3};

    //     frequency_proba.metering_mut_for_frequency(resum);

    //     REQUIRE(frequency_proba.Mut_count_by_site.size() == 5);
    //     REQUIRE(frequency_proba.Mut_count_by_site[0].size() == 1);
    //     REQUIRE(frequency_proba.Mut_count_by_site[1].size() == 2);
    //     REQUIRE(frequency_proba.Mut_count_by_site[2].size() == 0);
    //     REQUIRE(frequency_proba.Mut_count_by_site[3].size() == 1);
    //     REQUIRE(frequency_proba.Mut_count_by_site[4].size() == 1);

    //     REQUIRE(std::get<0>(frequency_proba.Mut_count_by_site[0][0]) == 0);
    //     REQUIRE(std::get<1>(frequency_proba.Mut_count_by_site[0][0]) == 1);

    //     REQUIRE(std::get<0>(frequency_proba.Mut_count_by_site[1][0]) == 0);
    //     REQUIRE(std::get<1>(frequency_proba.Mut_count_by_site[1][0]) == 1);

    //     REQUIRE(std::get<0>(frequency_proba.Mut_count_by_site[1][1]) == 3);
    //     REQUIRE(std::get<1>(frequency_proba.Mut_count_by_site[1][1]) == 1);

    //     REQUIRE(std::get<0>(frequency_proba.Mut_count_by_site[3][0]) == 3);
    //     REQUIRE(std::get<1>(frequency_proba.Mut_count_by_site[3][0]) == 1);

    //     REQUIRE(std::get<0>(frequency_proba.Mut_count_by_site[4][0]) == 0);
    //     REQUIRE(std::get<1>(frequency_proba.Mut_count_by_site[4][0]) == 1);
    // }

    // SECTION("metering_mut_for_frequency_complexe_case")
    // {
    //     int number_of_site = 5;
    //     std::vector<int> ancestry_seq(number_of_site, 1);
    //     int sample_size = 5;

    //     freq_allele_c frequency_proba(ancestry_seq, sample_size);

    //     std::vector<std::array<int, 2>> indiv1(3);
    //     indiv1[0] = {0, 0};
    //     indiv1[1] = {2, 2};
    //     indiv1[2] = {3, 3};

    //     std::vector<std::array<int, 2>> indiv2(2);
    //     indiv2[0] = {1, 0};
    //     indiv2[1] = {4, 0};

    //     std::vector<std::array<int, 2>> indiv3(2);
    //     indiv3[0] = {1, 3};
    //     indiv3[1] = {2, 2};

    //     std::vector<std::array<int, 2>> indiv4(3);
    //     indiv4[0] = {1, 2};
    //     indiv4[1] = {2, 0};
    //     indiv4[2] = {3, 3};

    //     std::vector<std::array<int, 2>> indiv5(3);
    //     indiv5[0] = {0, 2};
    //     indiv5[1] = {2, 3};
    //     indiv5[2] = {4, 0};

    //     std::vector<std::vector<std::array<int, 2>>> resum = {indiv1, indiv2, indiv3, indiv4, indiv5};

    //     frequency_proba.metering_mut_for_frequency(resum);

    //     REQUIRE(frequency_proba.Mut_count_by_site.size() == 5);
    //     REQUIRE(frequency_proba.Mut_count_by_site[0].size() == 2);
    //     REQUIRE(frequency_proba.Mut_count_by_site[1].size() == 3);
    //     REQUIRE(frequency_proba.Mut_count_by_site[2].size() == 3);
    //     REQUIRE(frequency_proba.Mut_count_by_site[3].size() == 1);
    //     REQUIRE(frequency_proba.Mut_count_by_site[4].size() == 1);

    //     REQUIRE(std::get<0>(frequency_proba.Mut_count_by_site[0][0]) == 0);
    //     REQUIRE(std::get<1>(frequency_proba.Mut_count_by_site[0][0]) == 1);

    //     REQUIRE(std::get<0>(frequency_proba.Mut_count_by_site[0][1]) == 2);
    //     REQUIRE(std::get<1>(frequency_proba.Mut_count_by_site[0][1]) == 1);

    //     REQUIRE(std::get<0>(frequency_proba.Mut_count_by_site[1][0]) == 0);
    //     REQUIRE(std::get<1>(frequency_proba.Mut_count_by_site[1][0]) == 1);

    //     REQUIRE(std::get<0>(frequency_proba.Mut_count_by_site[1][1]) == 2);
    //     REQUIRE(std::get<1>(frequency_proba.Mut_count_by_site[1][1]) == 1);

    //     REQUIRE(std::get<0>(frequency_proba.Mut_count_by_site[1][2]) == 3);
    //     REQUIRE(std::get<1>(frequency_proba.Mut_count_by_site[1][2]) == 1);

    //     REQUIRE(std::get<0>(frequency_proba.Mut_count_by_site[2][0]) == 0);
    //     REQUIRE(std::get<1>(frequency_proba.Mut_count_by_site[2][0]) == 1);

    //     REQUIRE(std::get<0>(frequency_proba.Mut_count_by_site[2][1]) == 2);
    //     REQUIRE(std::get<1>(frequency_proba.Mut_count_by_site[2][1]) == 2);

    //     REQUIRE(std::get<0>(frequency_proba.Mut_count_by_site[2][2]) == 3);
    //     REQUIRE(std::get<1>(frequency_proba.Mut_count_by_site[2][2]) == 1);

    //     REQUIRE(std::get<0>(frequency_proba.Mut_count_by_site[3][0]) == 3);
    //     REQUIRE(std::get<1>(frequency_proba.Mut_count_by_site[3][0]) == 2);

    //     REQUIRE(std::get<0>(frequency_proba.Mut_count_by_site[4][0]) == 0);
    //     REQUIRE(std::get<1>(frequency_proba.Mut_count_by_site[4][0]) == 2);
    // }
}
