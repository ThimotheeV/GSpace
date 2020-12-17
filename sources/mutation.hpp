#pragma once
#include <vector>
#include <tuple>
#include <map>
#include <iostream>
#include <random>
#include <list>

#include "random_generator.hpp"
#include "coalescence_table.hpp"
#include "settings.hpp"
struct map_string_nuc_name_c
{
    //for glib 0 is equal to missing value
    std::map<std::string, int> Map{
        {"A", 1},
        {"a", 1},
        {"T", 2},
        {"t", 2},
        {"C", 3},
        {"c", 3},
        {"G", 4},
        {"g", 4}};
    //for glib 0 is equal to missing value
    std::map<int, std::string> Inverse_map{
        {1, "A"},
        {2, "T"},
        {3, "C"},
        {4, "G"}};
};

enum class mut_model_enum
{
    iam,
    kam,
    smm,
    gsm,
    /*****/
    jcm,
    k80,
    f81,
    hky,
    tn93
};

/**
 * @brief Map from parameter string to param enum
 * 
 */
struct map_string_mod_mut_name_c
{
    std::map<std::string, mut_model_enum> Map{
        {"iam", mut_model_enum::iam},
        {"kam", mut_model_enum::kam},
        {"smm", mut_model_enum::smm},
        {"gsm", mut_model_enum::gsm},
        /*****/
        {"jcm", mut_model_enum::jcm},
        {"k80", mut_model_enum::k80},
        {"f81", mut_model_enum::f81},
        {"hky", mut_model_enum::hky},
        {"tn93", mut_model_enum::tn93}};

    std::map<mut_model_enum, std::string> ReverseMap{
        {mut_model_enum::iam, "iam"},
        {mut_model_enum::kam, "kam"},
        {mut_model_enum::smm, "smm"},
        {mut_model_enum::gsm, "gsm"},
        /*****/
        {mut_model_enum::jcm, "jcm"},
        {mut_model_enum::k80, "k80"},
        {mut_model_enum::f81, "f81"},
        {mut_model_enum::hky, "hky"},
        {mut_model_enum::tn93, "tn93"}};
};

//}

//create cumulative distribution for nucleotide relative frequency changes ( = one row of mutation matrix)
std::array<int, 4> cumul_distrib(std::array<double, 4> prob_vec);
//draws from the ditribution computed above
int nucl_distri(std::array<int, 4> const &cumul_prob_vec, double random);

// std::vector<std::tuple<childs_node, time_t, gen), nbr_leaf_for_node>>;
using top_down_tree_t = std::vector<std::tuple<std::vector<int>, double, int, int>>;
// vector of maps (one sample = one map) map of < site nbr, state >
using sample_mutated_state_type = std::vector<std::map<int, int>>;

/**************************************************************/
/*                        Mut mod                             */
/**************************************************************/

struct mut_mod_c_iam
{
    //not necessary because there is no other constructor than the empty constructor
    mut_mod_c_iam() = default;
    explicit mut_mod_c_iam(int nbr_site);

    void apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site);

private:
    std::vector<int> Count;
};

struct mut_mod_c_kam
{
    mut_mod_c_kam() = default; // equiv to mut_mod_c_kam(); because the default constructor is the empty one
    mut_mod_c_kam(int k_min, int k_max, rand_gen_c *rand_gen);

    void apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site);

    rand_gen_c *Rand_gen{nullptr};
    std::uniform_int_distribution<> Uni_int_kmin_kmax_distrib;
};

struct mut_mod_c_smm
{
    mut_mod_c_smm() = default;
    mut_mod_c_smm(int k_min, int k_max, rand_gen_c *rand_gen);

    void apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site);

    int K_min{1};
    int K_max{1};

    rand_gen_c *Rand_gen{nullptr};
};

struct mut_mod_c_gsm
{
    mut_mod_c_gsm() = default;

    mut_mod_c_gsm(int k_min, int k_max, double p_gsm, rand_gen_c *rand_gen);

    void apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site);

    int K_min{1};
    int K_max{1};

    rand_gen_c *Rand_gen{nullptr};
    std::geometric_distribution<> Geo_distrib;
};

// //Keep in memory the next "mutation number representation" (power of 2) for the site.
// struct mut_mod_c_ism
// {
//     mut_mod_c_ism() = default;
//     explicit mut_mod_c_ism(int number_site);

//     void apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site);

//     std::vector<int> Next_indivs_mut_state_t_vector;
// };

struct mut_mod_c_jcm
{
    mut_mod_c_jcm() = default;
    explicit mut_mod_c_jcm(rand_gen_c *rand_gen);

    void apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site);

    rand_gen_c *Rand_gen{nullptr};
    //Nasty but bug proof
    std::array<int, 4> Cumul_distrib_nucl_A{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
    std::array<int, 4> Cumul_distrib_nucl_C{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
    std::array<int, 4> Cumul_distrib_nucl_G{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
    std::array<int, 4> Cumul_distrib_nucl_T{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
};

struct mut_mod_c_k80
{
    mut_mod_c_k80() = default;
    mut_mod_c_k80(std::array<double, 2> ratio_transi_transver, rand_gen_c *rand_gen);

    void apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site);

    rand_gen_c *Rand_gen{nullptr};

    std::array<int, 4> Cumul_distrib_nucl_A{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
    std::array<int, 4> Cumul_distrib_nucl_C{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
    std::array<int, 4> Cumul_distrib_nucl_G{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
    std::array<int, 4> Cumul_distrib_nucl_T{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
};

struct mut_mod_c_f81
{
    mut_mod_c_f81() = default;
    mut_mod_c_f81(equi_base_freq_c equi_base_freq, rand_gen_c *rand_gen);

    void apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site);

    rand_gen_c *Rand_gen{nullptr};

    std::array<int, 4> Cumul_distrib_nucl_A{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
    std::array<int, 4> Cumul_distrib_nucl_C{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
    std::array<int, 4> Cumul_distrib_nucl_G{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
    std::array<int, 4> Cumul_distrib_nucl_T{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
};

struct mut_mod_c_hky
{
    mut_mod_c_hky() = default;
    mut_mod_c_hky(std::array<double, 2> ratio_transi_transver, equi_base_freq_c equi_base_freq, rand_gen_c *rand_gen);

    void apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site);

    rand_gen_c *Rand_gen{nullptr};

    std::array<int, 4> Cumul_distrib_nucl_A{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
    std::array<int, 4> Cumul_distrib_nucl_C{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
    std::array<int, 4> Cumul_distrib_nucl_G{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
    std::array<int, 4> Cumul_distrib_nucl_T{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
};

struct mut_mod_c_tn93
{
    mut_mod_c_tn93() = default;
    mut_mod_c_tn93(std::array<double, 2> ratio_transi_transver, equi_base_freq_c equi_base_freq, rand_gen_c *rand_gen);

    void apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site);

    rand_gen_c *Rand_gen{nullptr};

    std::array<int, 4> Cumul_distrib_nucl_A{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
    std::array<int, 4> Cumul_distrib_nucl_C{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
    std::array<int, 4> Cumul_distrib_nucl_G{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
    std::array<int, 4> Cumul_distrib_nucl_T{25 * PRECISION / 100, 50 * PRECISION / 100, 75 * PRECISION / 100, 1 * PRECISION};
};

/**************************************************************/
/*                        Mut struct_arg                      */
/**************************************************************/

class mut_model_c
{
public:
    //looks like an empty constructor but no ! (necessary because of the emptyy default constructor of mut_genenerator_c)
    //it is implemented in mutation.cpp, and fill mut_model with global variables
    mut_model_c();

    void apply_mut_to_site(mut_model_enum mut_name, std::map<int, int>::iterator indivs_mut_state_node_site);

    //needs to initialize all models for the switch in choose_mod()
    mut_mod_c_iam Mut_mod_iam;
    mut_mod_c_kam Mut_mod_kam;
    mut_mod_c_smm Mut_mod_smm;
    mut_mod_c_gsm Mut_mod_gsm;
    /*****/
    //mut_mod_c_ism Mut_mod_ism;
    mut_mod_c_jcm Mut_mod_jcm;
    mut_mod_c_k80 Mut_mod_k80;
    mut_mod_c_f81 Mut_mod_f81;
    mut_mod_c_hky Mut_mod_hky;
    mut_mod_c_tn93 Mut_mod_tn93;
};

struct mut_genenerator_c
{
    sample_mutated_state_type &mut_generation(top_down_tree_t const &gene_tree, std::vector<int> const &ancestry_seq, std::array<int, 2> const &begin_end_seq_current_tree, int MRCA_index);

    void browse_tree(bool approx, top_down_tree_t const &gene_tree, std::vector<int> &child_nodes, std::vector<int> const &ancestry_seq);

    void mut_process(bool approx, top_down_tree_t const &gene_tree, int ancester_index, int node_index, std::vector<int> const &ancestry_seq);

    void mut_events(bool approx, int node_index, double long_branch, double mut_rate, std::vector<int> const &ancestry_seq, rand_gen_c &rand_gen);

    //Membre
    mut_model_c Mut_model;
    std::array<int, 2> Begin_end_sequence{0, 0};
    sample_mutated_state_type All_indivs_mut_state; // mutated sites and states for each sample, for the current chunk
};

// make all the mutation process on all trees
// set each tree, create and apply mut_generator on each tree, and return the merged map of mutated sites and states for all chunks
std::vector<std::vector<std::pair<int, int>>> apply_mut_to_sample(coa_table_c &coa_table, int next_max_node_nbr, std::vector<int> const &ancestry_seq, muta_param_c const &muta_param, samp_param_c const &samp_param, demo_param_c const &demo_param);
