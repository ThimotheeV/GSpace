#define CATCH_CONFIG_MAIN LibTest
#include "catch.hpp"
#include "data_plane_vec_lib.hpp"

TEST_CASE("data_plane_vec_lib_interface_test")
{
    SECTION("3 haplo indiv, 3 deme, 1 locus")
    {
        //Simple case : 3 haplo indivs, 3 deme, 1 locus
        std::vector<std::vector<std::vector<std::pair<int, int>>>> sample_mutated_state_vec = {{{{0, 2}},
                                                                                                      {{0, 3}},
                                                                                                      {{0, 4}}}};
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 1},
                                       {0, 2}};
        samp_param.Ploidy = 1;
        samp_param.Sample_size_per_node = {1};
        samp_param.n_total_sample_size = 3;
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 1;
        samp_param.Dist_class_nbr = 5;

        std::vector<std::vector<int>> ancestry_seq = {{1}};
        data_plane_vec_c data_vec(sample_mutated_state_vec, samp_param, ancestry_seq);

        REQUIRE(data_vec.get_Ploidy() == 1);
        REQUIRE(data_vec.size() == 3);
        REQUIRE(data_vec.nbr_of_deme() == 3);
        REQUIRE(data_vec.base_nbr_locus_per_indiv() == 1); //number of locus in a specifique indiv
        REQUIRE(data_vec.nbr_of_locus_tot() == 3);
        REQUIRE(data_vec.nbr_of_indiv() == 3);

        REQUIRE(data_vec.nbr_of_indiv_per_deme(0) == 1); //deme size
        REQUIRE(data_vec.nbr_of_indiv_per_deme(1) == 1);
        REQUIRE(data_vec.nbr_of_indiv_per_deme(2) == 1);

        REQUIRE(data_vec.cumul_nbr_of_indiv_per_deme() == std::vector<int>{0, 1, 2});

        REQUIRE(data_vec.get_indiv(0) == 0);
        REQUIRE(data_vec.get_indiv(1) == 1);
        REQUIRE(data_vec.get_indiv(2) == 2);

        REQUIRE(data_vec.get_feature(0).Deme == 0);
        REQUIRE(data_vec.get_feature(1).Deme == 1);
        REQUIRE(data_vec.get_feature(2).Deme == 2);

        REQUIRE(data_vec.nomiss_nbr_of_gene_per_loc(0) == 3);

        REQUIRE(data_vec.nomiss_nbr_of_indiv_per_loc(0) == 3);

        REQUIRE(data_vec.nomiss_nbr_of_gene_per_loc_per_deme(0) == std::vector<int>{1, 1, 1});

        REQUIRE(data_vec.nomiss_nbr_of_indiv_per_loc_per_deme(0) == std::vector<int>{1, 1, 1});

        REQUIRE(data_vec.nomiss_nbr_of_deme_per_loc(0) == 3);

        REQUIRE(data_vec.nbr_allele_per_loc(0) == 3);
        REQUIRE(data_vec.allele_state_per_loc(0) == std::vector<std::array<int, 2>>{{2, 1}, {3, 1}, {4, 1}});

        REQUIRE(data_vec.get_plane_vec() == std::vector<int>{2, 3, 4});

        REQUIRE(data_vec.same_loc_in_indiv(0, 0));
        REQUIRE(!data_vec.same_loc_in_indiv(0, 1));

        REQUIRE(data_vec.same_deme(0, 0));
        REQUIRE(!data_vec.same_deme(0, 1));

        REQUIRE(data_vec.nomiss_data_indiv(0)[0] == data_vec.nomiss_data_indiv(1)[0]);
        REQUIRE(data_vec.nomiss_data_indiv(1)[0] == data_vec.nomiss_data_indiv(2)[0]);

        REQUIRE(data_vec.dist_btw_deme(0, 0) == 0);
        REQUIRE(data_vec.dist_btw_deme(0, 1) == 1);
        REQUIRE(data_vec.dist_btw_deme(0, 2) == 2);
        REQUIRE(data_vec.dist_btw_deme(1, 0) == 1);
        REQUIRE(data_vec.dist_btw_deme(1, 2) == 1);

        REQUIRE(data_vec.dist_class_btw_deme(0, 0) == 0);
        REQUIRE(data_vec.dist_class_btw_deme(0, 1) == 2);
        REQUIRE(data_vec.dist_class_btw_deme(0, 2) == 4);
        REQUIRE(data_vec.dist_class_btw_deme(1, 0) == 2);
        REQUIRE(data_vec.dist_class_btw_deme(1, 2) == 2);
    }

    SECTION("3 diploid indiv, 3 deme, 1 locus")
    {
        //Simple case : 3 diploid indiv, 3 deme, 1 locus
        std::vector<std::vector<std::vector<std::pair<int, int>>>> sample_mutated_state_vec = {{{{0, 2}},
                                                                                                      {{0, 2}},
                                                                                                      {{0, 3}},
                                                                                                      {{0, 4}},
                                                                                                      {{0, 4}},
                                                                                                      {{0, 5}}}};
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 1},
                                       {0, 2}};
        samp_param.Ploidy = 2;
        samp_param.Sample_size_per_node = {1};
        samp_param.n_total_sample_size = 3;
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 1;
        samp_param.Dist_class_nbr = 3;

        std::vector<std::vector<int>> ancestry_seq = {{1}};
        data_plane_vec_c data_vec(sample_mutated_state_vec, samp_param, ancestry_seq);

        REQUIRE(data_vec.get_Ploidy() == 2);
        REQUIRE(data_vec.size() == 6);
        REQUIRE(data_vec.nbr_of_deme() == 3);
        REQUIRE(data_vec.base_nbr_locus_per_indiv() == 1); //number of locus in a specifique indiv
        REQUIRE(data_vec.nbr_of_locus_tot() == 6);
        REQUIRE(data_vec.nbr_of_indiv() == 3);

        REQUIRE(data_vec.nbr_of_indiv_per_deme(0) == 1); //deme size
        REQUIRE(data_vec.nbr_of_indiv_per_deme(1) == 1);
        REQUIRE(data_vec.nbr_of_indiv_per_deme(2) == 1);

        REQUIRE(data_vec.cumul_nbr_of_indiv_per_deme() == std::vector<int>{0, 1, 2});

        REQUIRE(data_vec.get_indiv(0) == 0);
        REQUIRE(data_vec.get_indiv(1) == 0);
        REQUIRE(data_vec.get_indiv(2) == 1);
        REQUIRE(data_vec.get_indiv(3) == 1);
        REQUIRE(data_vec.get_indiv(4) == 2);
        REQUIRE(data_vec.get_indiv(5) == 2);

        REQUIRE(data_vec.get_feature(0).Deme == 0);
        REQUIRE(data_vec.get_feature(1).Deme == 1);
        REQUIRE(data_vec.get_feature(2).Deme == 2);

        REQUIRE(data_vec.nomiss_nbr_of_gene_per_loc(0) == 6);

        REQUIRE(data_vec.nomiss_nbr_of_indiv_per_loc(0) == 3);

        REQUIRE(data_vec.nomiss_nbr_of_gene_per_loc_per_deme(0) == std::vector<int>{2, 2, 2});

        REQUIRE(data_vec.nomiss_nbr_of_indiv_per_loc_per_deme(0) == std::vector<int>{1, 1, 1});

        REQUIRE(data_vec.nomiss_nbr_of_deme_per_loc(0) == 3);

        REQUIRE(data_vec.nbr_allele_per_loc(0) == 4);
        REQUIRE(data_vec.allele_state_per_loc(0) == std::vector<std::array<int, 2>>{{2, 2}, {3, 1}, {4, 2}, {5, 1}});

        REQUIRE(data_vec.get_plane_vec() == std::vector<int>{2, 2, 3, 4, 4, 5});

        REQUIRE(data_vec.same_loc_in_indiv(0, 1));
        REQUIRE(!data_vec.same_loc_in_indiv(1, 2));
        REQUIRE(data_vec.same_loc_in_indiv(2, 3));
        REQUIRE(!data_vec.same_loc_in_indiv(2, 4));

        REQUIRE(data_vec.same_deme(0, 1));
        REQUIRE(!data_vec.same_deme(1, 2));
        REQUIRE(data_vec.same_deme(2, 3));
        REQUIRE(!data_vec.same_deme(2, 4));

        REQUIRE(data_vec.nomiss_data_indiv(0)[0] == data_vec.nomiss_data_indiv(1)[0]);
        REQUIRE(data_vec.nomiss_data_indiv(1)[0] == data_vec.nomiss_data_indiv(2)[0]);

        REQUIRE(data_vec.dist_btw_deme(0, 0) == 0);
        REQUIRE(data_vec.dist_btw_deme(0, 1) == 0);
        REQUIRE(data_vec.dist_btw_deme(0, 2) == 1);
        REQUIRE(data_vec.dist_btw_deme(1, 4) == 2);
        REQUIRE(data_vec.dist_btw_deme(1, 2) == 1);

        REQUIRE(data_vec.dist_class_btw_deme(0, 0) == 0);
        REQUIRE(data_vec.dist_class_btw_deme(0, 1) == 0);
        REQUIRE(data_vec.dist_class_btw_deme(0, 2) == 1);
        REQUIRE(data_vec.dist_class_btw_deme(1, 4) == 2);
        REQUIRE(data_vec.dist_class_btw_deme(1, 2) == 1);
    }

    SECTION("3 diploid indiv, 3 deme, 2 independant locus")
    {
        //Case : 2 independant locus, 3 diploid indiv, 3 deme, 1 locus
        std::vector<std::vector<std::vector<std::pair<int, int>>>> sample_mutated_state_vec = {{{{0, 1}},
                                                                                                      {{0, 1}},
                                                                                                      {{0, 6}},
                                                                                                      {{0, 7}},
                                                                                                      {{0, 7}},
                                                                                                      {{0, 4}}},
                                                                                                     {{{0, 2}},
                                                                                                      {{0, 2}},
                                                                                                      {{0, 3}},
                                                                                                      {{0, 4}},
                                                                                                      {{0, 4}},
                                                                                                      {{0, 5}}}};
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 1},
                                       {0, 2}};
        samp_param.Ploidy = 2;
        samp_param.Sample_size_per_node = {1};
        samp_param.n_total_sample_size = 3;
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 2;
        samp_param.Dist_class_nbr = 3;

        std::vector<std::vector<int>> ancestry_seq = {{1}};
        data_plane_vec_c data_vec(sample_mutated_state_vec, samp_param, ancestry_seq);

        REQUIRE(data_vec.get_Ploidy() == 2);
        REQUIRE(data_vec.size() == 12);
        REQUIRE(data_vec.nbr_of_deme() == 3);
        REQUIRE(data_vec.base_nbr_locus_per_indiv() == 2); //number of locus in a specifique indiv
        REQUIRE(data_vec.nbr_of_locus_tot() == 12);
        REQUIRE(data_vec.nbr_of_indiv() == 3);

        REQUIRE(data_vec.nbr_of_indiv_per_deme(0) == 1); //deme size
        REQUIRE(data_vec.nbr_of_indiv_per_deme(1) == 1);
        REQUIRE(data_vec.nbr_of_indiv_per_deme(2) == 1);

        REQUIRE(data_vec.cumul_nbr_of_indiv_per_deme() == std::vector<int>{0, 1, 2});

        REQUIRE(data_vec.get_indiv(0) == 0);
        REQUIRE(data_vec.get_indiv(1) == 0);
        REQUIRE(data_vec.get_indiv(2) == 1);
        REQUIRE(data_vec.get_indiv(3) == 1);
        REQUIRE(data_vec.get_indiv(4) == 2);
        REQUIRE(data_vec.get_indiv(5) == 2);

        REQUIRE(data_vec.get_indiv(6) == 0);
        REQUIRE(data_vec.get_indiv(7) == 0);
        REQUIRE(data_vec.get_indiv(8) == 1);
        REQUIRE(data_vec.get_indiv(9) == 1);
        REQUIRE(data_vec.get_indiv(10) == 2);
        REQUIRE(data_vec.get_indiv(11) == 2);

        REQUIRE(data_vec.get_feature(0).Deme == 0);
        REQUIRE(data_vec.get_feature(1).Deme == 1);
        REQUIRE(data_vec.get_feature(2).Deme == 2);

        REQUIRE(data_vec.nomiss_nbr_of_gene_per_loc(0) == 6);
        REQUIRE(data_vec.nomiss_nbr_of_indiv_per_loc(0) == 3);
        REQUIRE(data_vec.nomiss_nbr_of_gene_per_loc_per_deme(0) == std::vector<int>{2, 2, 2});
        REQUIRE(data_vec.nomiss_nbr_of_indiv_per_loc_per_deme(0) == std::vector<int>{1, 1, 1});
        REQUIRE(data_vec.nomiss_nbr_of_deme_per_loc(0) == 3);
        REQUIRE(data_vec.nbr_allele_per_loc(0) == 4);
        REQUIRE(data_vec.allele_state_per_loc(0) == std::vector<std::array<int, 2>>{{1, 2}, {4, 1}, {6, 1}, {7, 2}});

        REQUIRE(data_vec.nomiss_nbr_of_gene_per_loc(1) == 6);
        REQUIRE(data_vec.nomiss_nbr_of_indiv_per_loc(1) == 3);
        REQUIRE(data_vec.nomiss_nbr_of_gene_per_loc_per_deme(1) == std::vector<int>{2, 2, 2});
        REQUIRE(data_vec.nomiss_nbr_of_indiv_per_loc_per_deme(1) == std::vector<int>{1, 1, 1});
        REQUIRE(data_vec.nomiss_nbr_of_deme_per_loc(1) == 3);
        REQUIRE(data_vec.nbr_allele_per_loc(1) == 4);
        REQUIRE(data_vec.allele_state_per_loc(1) == std::vector<std::array<int, 2>>{{2, 2}, {3, 1}, {4, 2}, {5, 1}});

        REQUIRE(data_vec.get_plane_vec() == std::vector<int>{1, 1, 6, 7, 7, 4, 2, 2, 3, 4, 4, 5});

        REQUIRE(data_vec.same_loc_in_indiv(0, 1));
        REQUIRE(!data_vec.same_loc_in_indiv(0, 6));
        REQUIRE(!data_vec.same_loc_in_indiv(1, 2));
        REQUIRE(!data_vec.same_loc_in_indiv(1, 7));
        REQUIRE(data_vec.same_loc_in_indiv(2, 3));
        REQUIRE(!data_vec.same_loc_in_indiv(2, 8));
        REQUIRE(!data_vec.same_loc_in_indiv(2, 4));
        REQUIRE(!data_vec.same_loc_in_indiv(2, 9));

        REQUIRE(data_vec.same_deme(0, 1));
        REQUIRE(data_vec.same_deme(0, 6));
        REQUIRE(!data_vec.same_deme(1, 2));
        REQUIRE(!data_vec.same_deme(1, 8));
        REQUIRE(data_vec.same_deme(2, 3));
        REQUIRE(data_vec.same_deme(2, 8));
        REQUIRE(!data_vec.same_deme(2, 4));
        REQUIRE(!data_vec.same_deme(2, 10));

        REQUIRE(data_vec.nomiss_data_indiv(0)[0] == data_vec.nomiss_data_indiv(1)[0]);
        REQUIRE(data_vec.nomiss_data_indiv(1)[0] == data_vec.nomiss_data_indiv(2)[0]);

        REQUIRE(data_vec.dist_btw_deme(0, 0) == 0);
        REQUIRE(data_vec.dist_btw_deme(0, 1) == 0);
        REQUIRE(data_vec.dist_btw_deme(0, 2) == 1);
        REQUIRE(data_vec.dist_btw_deme(1, 4) == 2);
        REQUIRE(data_vec.dist_btw_deme(1, 2) == 1);

        REQUIRE(data_vec.dist_class_btw_deme(0, 0) == 0);
        REQUIRE(data_vec.dist_class_btw_deme(0, 1) == 0);
        REQUIRE(data_vec.dist_class_btw_deme(0, 2) == 1);
        REQUIRE(data_vec.dist_class_btw_deme(1, 4) == 2);
        REQUIRE(data_vec.dist_class_btw_deme(1, 2) == 1);
    }

    SECTION("4 diploid indiv, 2 deme, 2 locus")
    {
        std::vector<std::vector<std::vector<std::pair<int, int>>>> sample_mutated_state_vec = {{{{0, 2}},
                                                                                                      {{0, 2}, {1, 1}},
                                                                                                      {{1, 2}},
                                                                                                      {{0, 2}, {1, 2}},
                                                                                                      {{0, 0}},
                                                                                                      {{0, 0}, {1, 2}},
                                                                                                      {{1, 0}},
                                                                                                      {{0, 3}, {1, 2}}}};
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {1, 0}};
        samp_param.Ploidy = 2;
        samp_param.Sample_size_per_node = {2};
        samp_param.n_total_sample_size = 4;
        samp_param.Sequence_length = 2;
        samp_param.Chr_nbr = 1;
        samp_param.Dist_class_nbr = 10;

        std::vector<std::vector<int>> ancestry_seq = {{1, 1}};
        data_plane_vec_c data_vec(sample_mutated_state_vec, samp_param, ancestry_seq);

        REQUIRE(data_vec.get_Ploidy() == 2);
        REQUIRE(data_vec.size() == 16);
        REQUIRE(data_vec.nbr_of_deme() == 2);
        REQUIRE(data_vec.base_nbr_locus_per_indiv() == 2); //number of locus in a specifique indiv
        REQUIRE(data_vec.nbr_of_locus_tot() == 16);
        REQUIRE(data_vec.nbr_of_indiv() == 4);

        REQUIRE(data_vec.nbr_of_indiv_per_deme(0) == 2); //deme size
        REQUIRE(data_vec.nbr_of_indiv_per_deme(1) == 2);

        REQUIRE(data_vec.cumul_nbr_of_indiv_per_deme() == std::vector<int>{0, 2});

        REQUIRE(data_vec.get_indiv(0) == 0);
        REQUIRE(data_vec.get_indiv(1) == 0);
        REQUIRE(data_vec.get_indiv(2) == 1);
        REQUIRE(data_vec.get_indiv(3) == 1);
        REQUIRE(data_vec.get_indiv(4) == 2);
        REQUIRE(data_vec.get_indiv(5) == 2);
        REQUIRE(data_vec.get_indiv(6) == 3);
        REQUIRE(data_vec.get_indiv(7) == 3);
        REQUIRE(data_vec.get_indiv(8) == 0);
        REQUIRE(data_vec.get_indiv(9) == 0);
        REQUIRE(data_vec.get_indiv(10) == 1);
        REQUIRE(data_vec.get_indiv(11) == 1);
        REQUIRE(data_vec.get_indiv(12) == 2);
        REQUIRE(data_vec.get_indiv(13) == 2);
        REQUIRE(data_vec.get_indiv(14) == 3);
        REQUIRE(data_vec.get_indiv(15) == 3);

        REQUIRE(data_vec.get_feature(0).Deme == 0);
        REQUIRE(data_vec.get_feature(1).Deme == 0);
        REQUIRE(data_vec.get_feature(2).Deme == 1);
        REQUIRE(data_vec.get_feature(3).Deme == 1);

        REQUIRE(data_vec.nomiss_nbr_of_gene_per_loc(0) == 8);

        REQUIRE(data_vec.nomiss_nbr_of_indiv_per_loc(0) == 4);

        REQUIRE(data_vec.nomiss_nbr_of_gene_per_loc_per_deme(0) == std::vector<int>{4, 4});

        REQUIRE(data_vec.nomiss_nbr_of_indiv_per_loc_per_deme(0) == std::vector<int>{2, 2});

        REQUIRE(data_vec.nomiss_nbr_of_deme_per_loc(0) == 2);

        REQUIRE(data_vec.nbr_allele_per_loc(0) == 4);
        REQUIRE(data_vec.nbr_allele_per_loc(1) == 3);

        REQUIRE(data_vec.allele_state_per_loc(0) == std::vector<std::array<int, 2>>{{0, 2}, {1, 2}, {2, 3}, {3, 1}});
        REQUIRE(data_vec.allele_state_per_loc(1) == std::vector<std::array<int, 2>>{{0, 1}, {1, 3}, {2, 4}});

        REQUIRE(data_vec.get_plane_vec() == std::vector<int>{2, 2, 1, 2, 0, 0, 1, 3, 1, 1, 2, 2, 1, 2, 0, 2});

        REQUIRE(data_vec.same_loc_in_indiv(0, 1));
        REQUIRE(data_vec.same_loc_in_indiv(8, 9));
        REQUIRE(!data_vec.same_loc_in_indiv(9, 10));
        //Not in the same locus
        REQUIRE(!data_vec.same_loc_in_indiv(2, 10));

        REQUIRE(data_vec.same_deme(0, 1));
        REQUIRE(data_vec.same_deme(1, 10));
        REQUIRE(data_vec.same_deme(2, 11));
        REQUIRE(!data_vec.same_deme(2, 12));

        REQUIRE(data_vec.nomiss_data_indiv(0)[0] == data_vec.nomiss_data_indiv(1)[0]);
        REQUIRE(data_vec.nomiss_data_indiv(1)[0] == data_vec.nomiss_data_indiv(2)[0]);

        REQUIRE(data_vec.dist_btw_deme(0, 0) == 0);
        REQUIRE(data_vec.dist_btw_deme(0, 1) == 0);
        REQUIRE(data_vec.dist_btw_deme(0, 2) == 0);
        REQUIRE(data_vec.dist_btw_deme(1, 4) == 1);
        REQUIRE(data_vec.dist_btw_deme(2, 11) == 0);
    }

    SECTION("update 3 diploid indiv, 3 deme, 1 locus")
    {
        //Simple case : 3 diploid indiv, 3 deme, 1 locus
        std::vector<std::vector<std::vector<std::pair<int, int>>>> sample_mutated_state_vec = {{{{0, 2}},
                                                                                                      {{0, 2}},
                                                                                                      {{0, 3}},
                                                                                                      {{0, 4}},
                                                                                                      {{0, 4}},
                                                                                                      {{0, 5}}}};
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 1},
                                       {0, 2}};
        samp_param.Ploidy = 2;
        samp_param.Sample_size_per_node = {1};
        samp_param.n_total_sample_size = 3;
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 1;
        samp_param.Dist_class_nbr = 3;

        std::vector<std::vector<int>> ancestry_seq = {{1}};
        data_plane_vec_c data_vec(sample_mutated_state_vec, samp_param, ancestry_seq);

        std::vector<std::vector<std::vector<std::pair<int, int>>>> sample_mutated_state_vec2 = {{{{0, 0}},
                                                                                                       {{0, 0}},
                                                                                                       {{0, 1}},
                                                                                                       {{0, 1}},
                                                                                                       {{0, 2}},
                                                                                                       {{0, 2}}}};

        data_vec.update_data_plane_vec(sample_mutated_state_vec2, ancestry_seq, samp_param);

        REQUIRE(data_vec.get_Ploidy() == 2);
        REQUIRE(data_vec.size() == 6);
        REQUIRE(data_vec.nbr_of_deme() == 3);
        REQUIRE(data_vec.base_nbr_locus_per_indiv() == 1); //number of locus in a specifique indiv
        REQUIRE(data_vec.nbr_of_locus_tot() == 6);
        REQUIRE(data_vec.nbr_of_indiv() == 3);

        REQUIRE(data_vec.nbr_of_indiv_per_deme(0) == 1); //deme size
        REQUIRE(data_vec.nbr_of_indiv_per_deme(1) == 1);
        REQUIRE(data_vec.nbr_of_indiv_per_deme(2) == 1);

        REQUIRE(data_vec.cumul_nbr_of_indiv_per_deme() == std::vector<int>{0, 1, 2});

        REQUIRE(data_vec.get_indiv(0) == 0);
        REQUIRE(data_vec.get_indiv(1) == 0);
        REQUIRE(data_vec.get_indiv(2) == 1);
        REQUIRE(data_vec.get_indiv(3) == 1);
        REQUIRE(data_vec.get_indiv(4) == 2);
        REQUIRE(data_vec.get_indiv(5) == 2);

        REQUIRE(data_vec.get_feature(0).Deme == 0);
        REQUIRE(data_vec.get_feature(1).Deme == 1);
        REQUIRE(data_vec.get_feature(2).Deme == 2);

        REQUIRE(data_vec.nomiss_nbr_of_gene_per_loc(0) == 6);

        REQUIRE(data_vec.nomiss_nbr_of_indiv_per_loc(0) == 3);

        REQUIRE(data_vec.nomiss_nbr_of_gene_per_loc_per_deme(0) == std::vector<int>{2, 2, 2});

        REQUIRE(data_vec.nomiss_nbr_of_indiv_per_loc_per_deme(0) == std::vector<int>{1, 1, 1});

        REQUIRE(data_vec.nomiss_nbr_of_deme_per_loc(0) == 3);

        REQUIRE(data_vec.nbr_allele_per_loc(0) == 3);
        REQUIRE(data_vec.allele_state_per_loc(0) == std::vector<std::array<int, 2>>{{0, 2}, {1, 2}, {2, 2}});

        REQUIRE(data_vec.get_plane_vec() == std::vector<int>{0, 0, 1, 1, 2, 2});

        REQUIRE(data_vec.same_loc_in_indiv(0, 1));
        REQUIRE(!data_vec.same_loc_in_indiv(1, 2));
        REQUIRE(data_vec.same_loc_in_indiv(2, 3));
        REQUIRE(!data_vec.same_loc_in_indiv(2, 4));

        REQUIRE(data_vec.same_deme(0, 1));
        REQUIRE(!data_vec.same_deme(1, 2));
        REQUIRE(data_vec.same_deme(2, 3));
        REQUIRE(!data_vec.same_deme(2, 4));

        REQUIRE(data_vec.nomiss_data_indiv(0)[0] == data_vec.nomiss_data_indiv(1)[0]);
        REQUIRE(data_vec.nomiss_data_indiv(1)[0] == data_vec.nomiss_data_indiv(2)[0]);

        REQUIRE(data_vec.dist_btw_deme(0, 0) == 0);
        REQUIRE(data_vec.dist_btw_deme(0, 1) == 0);
        REQUIRE(data_vec.dist_btw_deme(0, 2) == 1);
        REQUIRE(data_vec.dist_btw_deme(1, 4) == 2);
        REQUIRE(data_vec.dist_btw_deme(1, 2) == 1);

        REQUIRE(data_vec.dist_class_btw_deme(0, 0) == 0);
        REQUIRE(data_vec.dist_class_btw_deme(0, 1) == 0);
        REQUIRE(data_vec.dist_class_btw_deme(0, 2) == 1);
        REQUIRE(data_vec.dist_class_btw_deme(1, 4) == 2);
        REQUIRE(data_vec.dist_class_btw_deme(1, 2) == 1);
    }

    SECTION("update 3 diploid indiv, 3 deme, 2 independant locus")
    {
        //Simple case : 3 diploid indiv, 3 deme, 1 locus
        std::vector<std::vector<std::vector<std::pair<int, int>>>> sample_mutated_state_vec = {{{{0, 1}},
                                                                                                      {{0, 1}},
                                                                                                      {{0, 6}},
                                                                                                      {{0, 7}},
                                                                                                      {{0, 7}},
                                                                                                      {{0, 4}}},
                                                                                                     {{{0, 2}},
                                                                                                      {{0, 2}},
                                                                                                      {{0, 3}},
                                                                                                      {{0, 4}},
                                                                                                      {{0, 4}},
                                                                                                      {{0, 5}}}};
        auto &samp_param = singleton_c<samp_param_c>::instance();
        samp_param.Sample_coord_vec = {{0, 0},
                                       {0, 1},
                                       {0, 2}};
        samp_param.Ploidy = 2;
        samp_param.Sample_size_per_node = {1};
        samp_param.n_total_sample_size = 3;
        samp_param.Sequence_length = 1;
        samp_param.Chr_nbr = 2;
        samp_param.Dist_class_nbr = 3;

        std::vector<std::vector<int>> ancestry_seq = {{1}};
        data_plane_vec_c data_vec(sample_mutated_state_vec, samp_param, ancestry_seq);

        std::vector<std::vector<std::vector<std::pair<int, int>>>> sample_mutated_state_vec2 = {{{{0, 0}},
                                                                                                       {{0, 0}},
                                                                                                       {{0, 1}},
                                                                                                       {{0, 1}},
                                                                                                       {{0, 2}},
                                                                                                       {{0, 2}}},
                                                                                                      {{{0, 3}},
                                                                                                       {{0, 3}},
                                                                                                       {{0, 4}},
                                                                                                       {{0, 4}},
                                                                                                       {{0, 5}},
                                                                                                       {{0, 5}}}};

        data_vec.update_data_plane_vec(sample_mutated_state_vec2, ancestry_seq, samp_param);

        REQUIRE(data_vec.get_Ploidy() == 2);
        REQUIRE(data_vec.size() == 12);
        REQUIRE(data_vec.nbr_of_deme() == 3);
        REQUIRE(data_vec.base_nbr_locus_per_indiv() == 2); //number of locus in a specifique indiv
        REQUIRE(data_vec.nbr_of_locus_tot() == 12);
        REQUIRE(data_vec.nbr_of_indiv() == 3);

        REQUIRE(data_vec.nbr_of_indiv_per_deme(0) == 1); //deme size
        REQUIRE(data_vec.nbr_of_indiv_per_deme(1) == 1);
        REQUIRE(data_vec.nbr_of_indiv_per_deme(2) == 1);

        REQUIRE(data_vec.cumul_nbr_of_indiv_per_deme() == std::vector<int>{0, 1, 2});

        REQUIRE(data_vec.get_indiv(0) == 0);
        REQUIRE(data_vec.get_indiv(1) == 0);
        REQUIRE(data_vec.get_indiv(2) == 1);
        REQUIRE(data_vec.get_indiv(3) == 1);
        REQUIRE(data_vec.get_indiv(4) == 2);
        REQUIRE(data_vec.get_indiv(5) == 2);

        REQUIRE(data_vec.get_indiv(6) == 0);
        REQUIRE(data_vec.get_indiv(7) == 0);
        REQUIRE(data_vec.get_indiv(8) == 1);
        REQUIRE(data_vec.get_indiv(9) == 1);
        REQUIRE(data_vec.get_indiv(10) == 2);
        REQUIRE(data_vec.get_indiv(11) == 2);

        REQUIRE(data_vec.get_feature(0).Deme == 0);
        REQUIRE(data_vec.get_feature(1).Deme == 1);
        REQUIRE(data_vec.get_feature(2).Deme == 2);

        REQUIRE(data_vec.nomiss_nbr_of_gene_per_loc(0) == 6);
        REQUIRE(data_vec.nomiss_nbr_of_indiv_per_loc(0) == 3);
        REQUIRE(data_vec.nomiss_nbr_of_gene_per_loc_per_deme(0) == std::vector<int>{2, 2, 2});
        REQUIRE(data_vec.nomiss_nbr_of_indiv_per_loc_per_deme(0) == std::vector<int>{1, 1, 1});
        REQUIRE(data_vec.nomiss_nbr_of_deme_per_loc(0) == 3);
        REQUIRE(data_vec.nbr_allele_per_loc(0) == 3);
        REQUIRE(data_vec.allele_state_per_loc(0) == std::vector<std::array<int, 2>>{{0, 2}, {1, 2}, {2, 2}});

        REQUIRE(data_vec.nomiss_nbr_of_gene_per_loc(1) == 6);
        REQUIRE(data_vec.nomiss_nbr_of_indiv_per_loc(1) == 3);
        REQUIRE(data_vec.nomiss_nbr_of_gene_per_loc_per_deme(1) == std::vector<int>{2, 2, 2});
        REQUIRE(data_vec.nomiss_nbr_of_indiv_per_loc_per_deme(1) == std::vector<int>{1, 1, 1});
        REQUIRE(data_vec.nomiss_nbr_of_deme_per_loc(1) == 3);
        REQUIRE(data_vec.nbr_allele_per_loc(1) == 3);
        REQUIRE(data_vec.allele_state_per_loc(1) == std::vector<std::array<int, 2>>{{3, 2}, {4, 2}, {5, 2}});

        REQUIRE(data_vec.get_plane_vec() == std::vector<int>{0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5});

        REQUIRE(data_vec.same_loc_in_indiv(0, 1));
        REQUIRE(!data_vec.same_loc_in_indiv(0, 6));
        REQUIRE(!data_vec.same_loc_in_indiv(1, 2));
        REQUIRE(!data_vec.same_loc_in_indiv(1, 7));
        REQUIRE(data_vec.same_loc_in_indiv(2, 3));
        REQUIRE(!data_vec.same_loc_in_indiv(2, 8));
        REQUIRE(!data_vec.same_loc_in_indiv(2, 4));
        REQUIRE(!data_vec.same_loc_in_indiv(2, 9));

        REQUIRE(data_vec.same_deme(0, 1));
        REQUIRE(data_vec.same_deme(0, 6));
        REQUIRE(!data_vec.same_deme(1, 2));
        REQUIRE(!data_vec.same_deme(1, 8));
        REQUIRE(data_vec.same_deme(2, 3));
        REQUIRE(data_vec.same_deme(2, 8));
        REQUIRE(!data_vec.same_deme(2, 4));
        REQUIRE(!data_vec.same_deme(2, 10));

        REQUIRE(data_vec.nomiss_data_indiv(0)[0] == data_vec.nomiss_data_indiv(1)[0]);
        REQUIRE(data_vec.nomiss_data_indiv(1)[0] == data_vec.nomiss_data_indiv(2)[0]);

        REQUIRE(data_vec.dist_btw_deme(0, 0) == 0);
        REQUIRE(data_vec.dist_btw_deme(0, 1) == 0);
        REQUIRE(data_vec.dist_btw_deme(0, 2) == 1);
        REQUIRE(data_vec.dist_btw_deme(1, 4) == 2);
        REQUIRE(data_vec.dist_btw_deme(1, 2) == 1);

        REQUIRE(data_vec.dist_class_btw_deme(0, 0) == 0);
        REQUIRE(data_vec.dist_class_btw_deme(0, 1) == 0);
        REQUIRE(data_vec.dist_class_btw_deme(0, 2) == 1);
        REQUIRE(data_vec.dist_class_btw_deme(1, 4) == 2);
        REQUIRE(data_vec.dist_class_btw_deme(1, 2) == 1);
    }
}
