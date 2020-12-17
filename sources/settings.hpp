/**
 * @file Settings.hpp
 * @author François-David Collin <Francois-David.Collin@umontpellier.fr>
 * @brief Global settings
 * @version 0.1
 * @date 2019-01-30
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#pragma once
#include <vector>
#include <map>

#include "arg_parser.hpp"
#include "random_generator.hpp"
#include "migration.hpp"

enum class edge_effect_enum;
enum class mut_model_enum;

struct equi_base_freq_c
{
    double A{0.3};
    double T{0.2};
    double C{0.3};
    double G{0.2};
};

//Uniq instance of type, can't be copy or manualy create
template <typename Type>
class singleton_c
{
public:
    static Type &instance();

    singleton_c(const singleton_c &) = delete;
    singleton_c(singleton_c &&) = delete;
    singleton_c &operator=(const singleton_c) = delete;
    singleton_c &operator=(singleton_c &&) = delete;

    ~singleton_c();

protected:
    singleton_c() = default;
};

//Create the uniq (static) instance of Type if it doesn't exist
//If it exist just return it
template <typename Type>
Type &singleton_c<Type>::instance()
{
    static Type instance{};
    return instance;
}

/**
 * @brief properties_c contains the actual parameters used by programs
 * 
 * 
 */

struct simu_param_c
{
    std::string Setting_filename = "GSpaceSettings.txt";
    std::string Generic_data_filename = "GSpace_Simu";
    std::string Data_file_extension = ".txt";
    std::string Nodesize_matrix_filename;
    std::string Migration_matrix_filename = "";
    std::string Output_dir = ".";
    std::string Param_summary_filename = "_GSpace_param_summary.txt";
    bool Wait_for_cin_input = false;
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    bool Wait_for_final_cin_input = true;
#else
    bool Wait_for_final_cin_input = false;
#endif
    bool Genepop_output = false;
    bool Genepop_ind_file = false;
    bool Genepop_group_all_samples = false;
    bool Genepop_no_coord = false;
    bool VCF_output = false;
    bool Fasta_output = false;
    bool Fasta_single_line_seq = false;
    bool Phylip_output = false;
    bool Coordinates_output = false;
    bool Seq_char_output = false;
    bool Nodesize_matrix = false;
    bool Migration_matrix = false;
    int seed{3}; //intial random seeds

    int Repetition_nbr{0};
    bool Continuous_time_approxim{false};
};

enum class dispersal_distrib_enum
{
    none,
    uniform,
    gaussian,
    geometric,
    pareto,
    sichel
};

struct map_string_dispersal_distrib_name_c
{
    std::map<dispersal_distrib_enum, std::string> ReverseMap{
        {dispersal_distrib_enum::uniform, "uniform"},
        {dispersal_distrib_enum::gaussian, "gaussian"},
        {dispersal_distrib_enum::geometric, "geometric"},
        {dispersal_distrib_enum::pareto, "pareto"},
        {dispersal_distrib_enum::sichel, "sichel"}};
};

struct samp_param_c
{
    int Ploidy{1};
    std::vector<int> Sample_size_per_node;

    int n_total_sample_size;
    std::vector<int> Sample_coord_x;//temp storage of sample coordinates,
    std::vector<int> Sample_coord_y;//temp storage of sample coordinates, 
    std::vector<std::array<int, 2>> Sample_coord_vec;
    int Dist_class_nbr{1};

    int Chr_nbr{1};
    int Sequence_length{1};
};

struct demo_param_c
{
    int Population_size_N{1};
    int Pop_size_per_node{1}; // default value of 1 is set in check_param()
    std::vector<std::vector<int>> Nodesize_mat;
    std::vector<std::vector<double>> Migration_mat;

    std::array<int, 2> Lattice_size{1, 1};
    std::array<int, 2> Disp_dist_max{0, 0};

    int Nbr_node_sampled_x{1};
    int Nbr_node_sampled_y{1};
    int Min_sample_coord_x{1};
    int Min_sample_coord_y{1};
    int Void_sample_nodes_X{1};
    int Void_sample_nodes_Y{1};

    edge_effect_enum Edge_effects;

    dispersal_distrib_enum Dispersal_distrib;
    std::array<std::function<double(int)>, 2> Disp_func;
    double Proba_migr{0};
    double G_geo_param{0.5}; //Geometric_shape;
    double N_pareto_param{2};
    double Gamma_sichel_param{0};
    double Xi_sichel_param{0};
    double Omega_sichel_param{0};

    std::array<double, 3> Sichel_pars{1, 1, 1};
};

struct recomb_param_c
{
    double Unscaled_recomb_rate{0.005};
    double Scaled_recomb_rate_rho{0};
};

struct muta_param_c
{
    mut_model_enum Mod_mut_name;
    double Scaled_mut_rate_theta{0};
    double Unscaled_mut_rate_mu{0.005};
    std::vector<std::vector<int>> Ancestry_seq;
    //If no other information, use nucl limit
    int K_min{1};
    int K_max{10};
    double P_gsm{0.22};
    std::array<double, 2> Ratio_transi_transver{0.5, 1};
    equi_base_freq_c Equi_base_freq;
};

/**
 * @brief Read settings file
 * 
 */

std::string const read_file(std::string const &filename);

std::string read_write_cmdline(int argc, char **argv);
void parser_str(std::string const &file_str);
void check_param();
void apply_param();
void compute_sample_coord(samp_param_c &samp_param, demo_param_c const &demo_param);
void output_screen_info(std::string const &version, simu_param_c const &simu_param, samp_param_c const &samp_param, demo_param_c const &demo_param, muta_param_c const &muta_param, recomb_param_c const &recomb_param);
