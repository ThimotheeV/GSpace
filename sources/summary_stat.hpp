#pragma once
#include <vector>

#include "migration.hpp"
#include "calc_stat.hpp"
#include "calc_stat_dl.hpp"

using sample_mutated_state_type = std::vector<std::map<int, int>>;

std::vector<int> reconstitution_full_int_seq_indiv(std::vector<int> ancestry_seq, std::vector<std::pair<int, int>> const &indiv_mut_seq);
std::string reconstitution_full_nucl_seq_indiv(std::vector<int> ancestry_seq, std::vector<std::pair<int, int>> const &indiv_mut_seq);

//WARNING : Hacky approch, not very efficient
struct Prob_id_1_2_loc_res_c
{
    void cumul_prob_id_1_2_loc(data_plane_vec_c const &data_plane_vec, int rep);

    double PHI{0};
    std::array<double, 2> PHI_cumul_m_v = {0, 0};
};

//WARNING : Hacky approch, not very efficient
struct Prob_id_1_loc_Qr_res_c
{
    void calc_Qr(data_plane_vec_c const &data_plane_vec, int rep);
    //Qr[dist]
    std::vector<std::vector<std::vector<double>>> Qr_by_chr_by_loc;
    std::vector<std::vector<double>> Qr_by_chr;
    std::vector<double> Qr;
    std::vector<std::array<double, 2>> Qr_cumul_m_v;
};

//WARNING : Hacky approch, not very efficient
struct Prob_id_1_loc_Qwi_wd_bd_res_c
{
    //Qwi : Prob id at 1 loc within indiv
    //Qbi : Prob id at 1 loc between indiv
    //Qwd : Prob id at 1 loc within deme
    //Qbd : Prob id at 1 loc between deme

    void calc_Qwi(data_plane_vec_c const &data_plane_vec, int rep);
    void calc_Qwd_bd(data_plane_vec_c const &data_plane_vec, int rep);

    std::vector<std::vector<double>> Qwi_by_chr_by_loc;
    std::vector<double> Qwi_by_chr;
    double Qwi{0};
    std::array<double, 2> Qwi_cumul_m_v{0, 0};

    std::vector<std::vector<double>> Qwd_by_chr_by_loc;
    std::vector<double> Qwd_by_chr;
    double Qwd{0};
    std::array<double, 2> Qwd_cumul_m_v{0, 0};

    std::vector<std::vector<double>> Qbd_by_chr_by_loc;
    std::vector<double> Qbd_by_chr;
    double Qbd{0};
    std::array<double, 2> Qbd_cumul_m_v{0, 0};
};

struct freq_allele_c
{
    freq_allele_c(std::vector<int> const &ancestry_seq, int sample_size) : Ancestry_seq(ancestry_seq),
                                                                           Mut_count_by_site(ancestry_seq.size()),
                                                                           Sample_size(sample_size){};

    void metering_mut_for_frequency(sample_mutated_state_type const &sample_mutated_state);

    std::vector<int> const &Ancestry_seq;
    std::vector<std::vector<std::array<int, 2>>> Mut_count_by_site;
    int Sample_size;
};

struct coa_tree_metrics_c
{
    void calcul_coa_tree_metrics(bool approx, int ploidy, int population_size, int sample_size);

    double Theo_MRCA_mean{0};
    double Theo_MRCA_var{0};
    double Theo_2_lign_coa_time_mean{0};
    double Theo_2_lign_coa_time_var{0};
};

struct info_collector_c
{
    bool Debug = false;
    bool Check_tree_bool = false;
    bool Clock = false;
    double time_coa{0};
    double time_recomb{0};
    double time_mig{0};
    double time_construct_tree{0};
    double time_mutation{0};
    double time_simulation{0};

    bool Stats{false};
    bool Iterative_stats{false};
    bool Stats_per_loc_per_simu{false};
    bool Stats_per_chr_per_simu{false};

    bool Prob_id_1_2_loc{false};
    bool Prob_id_1_loc_Qr{false};
    bool Prob_id_1_loc_Qwi_wd_bd{false};

    bool Coa_times{false};

    bool Effective_disp{false};
    //bool Check_disp_distrib{false}; // not used anymore, may be reused if we don't want to store and output dispersal distributions/Migration matrix
    bool Print_migration_matrix_bool{false};

    std::vector<double> MRCA_time;
    std::array<double, 2> MRCA_times_cumul_mean_var{0, 0};

    std::vector<double> Gen_coa_time;
    std::array<double, 2> Gen_coa_cumul_mean_var{0, 0};
    std::vector<int> Nbr_coa;

    std::array<double, 3> Nbr_depl_0_1_tot{0, 0, 0};
    std::vector<std::vector<int>> Emp_cumul_axial_disp;
    std::vector<std::vector<double>> Emp_axial_disp;
    std::vector<double> Mean_emp_axial_disp, Sig_emp_axial_disp, Kurt_emp_axial_disp, Skew_emp_axial_disp;
    std::vector<std::vector<double>> Emp_axial_disp_mean_over_rep;
    std::vector<double> Mean_emp_axial_disp_mean_over_rep, Sig_emp_axial_disp_mean_over_rep, Kurt_emp_axial_disp_mean_over_rep, Skew_emp_axial_disp_mean_over_rep;

    std::vector<double> Fwd_axial_disp_distrib;
    double Fwd_axial_disp_mean, Fwd_axial_disp_theo_sig, Fwd_axial_disp_sig, Fwd_axial_disp_kurt, Fwd_axial_disp_skew, Ds2;
    std::vector<std::vector<double>> Cumul_fwd_disp_distrib;
    std::vector<std::vector<double>> Bcwd_disp_distrib_from_migration_matrix;
    std::vector<std::vector<double>> Cumul_bcwd_disp_distrib_from_migration_matrix;

    int rep;
};

struct output_stat_c
{
    //Speed functional test (no need to write calc values)
    bool Stat_out{true};
    Prob_id_1_2_loc_res_c Prob_id_1_2_loc_res;
    Prob_id_1_loc_Qr_res_c Prob_id_1_loc_Qr_res;
    Prob_id_1_loc_Qwi_wd_bd_res_c Prob_id_1_loc_Qwi_wd_bd_res;
    
    data_plane_vec_c Data_plane_vec;
};

//tuple(mean, var, pond)
std::array<double, 2>
calc_cumul_mean_var(std::array<double, 2> cumul_mean_var_pond, std::vector<double> const &vec_value, int pond);
