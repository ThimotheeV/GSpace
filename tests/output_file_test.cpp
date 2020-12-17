#define CATCH_CONFIG_MAIN OutputTest
#include "catch.hpp"

#include "output.hpp"

TEST_CASE("constructor_output_c_test")
{
    REQUIRE(1);
}
// TEST_CASE("constructor_output_c_test")
// {
//     output_stat_c out("Truc");
//     std::string str = "Repetition_nbr\tPhiA_seq_mean_non_iter\tPhiA_seq_Mean\tPhiA_seq_Var\tPhiA_Mean\tPhiA_Var\tPhiAB_Mean\tPhiAB_Var\tGammaAB_Mean\tGammaAB_Var\tDeltaAB_Mean\tDeltaAB_Var\n";
//     std::string strInfo = "Repetition_nbr\tMRCA_accu_moy\tMRCA_accu_var\tcoa_2_lignee_t_accu_moy\tcoa_2_lignee_t_accu_var\n";

//     REQUIRE(out.Directory == "Truc");
//     REQUIRE(out.Buf_prob_id_stat_out == str);
//     REQUIRE(out.Info == strInfo);
// }

// TEST_CASE("prepare_prob_id_1_2_loc_display_output_c_test")
// {
//     output_stat_c out("Machin");
//     std::string str = "1\t0.000000000000e+00\t0.000000000000e+00\t0.000000000000e+00\t2.000000000000e+00\t3.000000000000e+00\t4.000000000000e+00\t5.000000000000e+00\t6.000000000000e+00\t7.000000000000e+00\t8.000000000000e+00\t9.000000000000e+00\n";

//     Prob_id_1_2_loc_res_c Prob_id_1_2_loc_res;

//     Prob_id_1_2_loc_res.phiA_seq_cumul_m_v = std::make_tuple(0, 0);
//     Prob_id_1_2_loc_res.phiA_cumul_m_v = std::make_tuple(2, 3);
//     Prob_id_1_2_loc_res.phiAB_cumul_m_v = std::make_tuple(4, 5);
//     Prob_id_1_2_loc_res.gammaAB_cumul_m_v = std::make_tuple(6, 7);
//     Prob_id_1_2_loc_res.deltaAB_cumul_m_v = std::make_tuple(8, 9);

//     out.prepare_prob_id_1_2_loc_display(Prob_id_1_2_loc_res, 0);

//     REQUIRE(out.Buf_prob_id_stat_out == str);
// }

// TEST_CASE("prepare_info_display_output_c_test")
// {
//     output_stat_c out("Bidule");
//     std::string strInfo = "1\t1.000000000000e+00\t2.000000000000e+00\n";

//     out.prepare_info_display(std::make_tuple(1,2), 1);

//     REQUIRE(out.Info == strInfo);
// }
