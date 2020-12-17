#include <iostream>
#include <cassert>

#include "mutation.hpp"
#include "summary_stat.hpp"

std::array<int, 4> cumul_distrib(std::array<double, 4> prob_vec)
{
    double sum = 0;
    for (auto &value : prob_vec)
    {
        value += sum;
        sum = value;
    }

    std::array<int, 4> result{0, 0, 0, 0};
    for (int i = 0; i < 4; ++i)
    {
        result.at(i) = prob_vec.at(i) / sum * PRECISION;
    }

    return result;
}

int nucl_distri(std::array<int, 4> const &cumul_prob_vec, double random)
{
    int nucl = 1;
    // <= right continuity
    while (random > cumul_prob_vec.at(nucl - 1))
    {
        ++nucl;
    }
    return nucl;
}

/**************************************************************/
/*                        Mut mod                             */
/**************************************************************/
mut_mod_c_iam::mut_mod_c_iam(int nbr_site)
{
    Count.resize(nbr_site, 1);
}

void mut_mod_c_iam::apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site)
{
    indivs_mut_state_node_site->second = ++Count.at(indivs_mut_state_node_site->first);
}

mut_mod_c_kam::mut_mod_c_kam(int k_min, int k_max, rand_gen_c *rand_gen)
{
    Rand_gen = rand_gen;

    Uni_int_kmin_kmax_distrib = std::uniform_int_distribution<>(k_min, k_max);
}

void mut_mod_c_kam::apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site)
{
    auto temp = Uni_int_kmin_kmax_distrib(Rand_gen->Seed_gen);
    while (temp == indivs_mut_state_node_site->second)
    {
        temp = Uni_int_kmin_kmax_distrib(Rand_gen->Seed_gen);
    }
    indivs_mut_state_node_site->second = temp;
}

mut_mod_c_smm::mut_mod_c_smm(int k_min, int k_max, rand_gen_c *rand_gen)
{
    Rand_gen = rand_gen;

    K_min = k_min;
    K_max = k_max;
}

void mut_mod_c_smm::apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site)
{
    if (Rand_gen->rand_bool())
    {
        if (indivs_mut_state_node_site->second != K_max)
        {
            indivs_mut_state_node_site->second += 1;
        }
    }
    else
    {
        if (indivs_mut_state_node_site->second != K_min)
        {
            indivs_mut_state_node_site->second -= 1;
        }
    }
}

mut_mod_c_gsm::mut_mod_c_gsm(int k_min, int k_max, double p_gsm, rand_gen_c *rand_gen)
{
    Rand_gen = rand_gen;

    K_min = k_min;
    K_max = k_max;

    //use q = 1-p instead off p in geometric distribution
    Geo_distrib = std::geometric_distribution<>(p_gsm); // <> = default template type = int
}

void mut_mod_c_gsm::apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site)
{
    auto repet = Geo_distrib(Rand_gen->Seed_gen) + 1;

    if (Rand_gen->rand_bool())
    {
        if (indivs_mut_state_node_site->second + repet <= K_max)
        {
            indivs_mut_state_node_site->second += repet;
        }
    }
    else
    {
        if (indivs_mut_state_node_site->second - repet >= K_min)
        {
            indivs_mut_state_node_site->second -= repet;
        }
    }
}

// mut_mod_c_ism::mut_mod_c_ism(int site_nbr)
// {
//     // value for next mutated state, ini 1
//     Next_indivs_mut_state_t_vector = std::vector<int>(site_nbr, 1);
// }

// void mut_mod_c_ism::apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site)
// {
//     //Keep in memory the next "mutation number representation" (power of 2) for the site.
//     auto next_indivs_mut_state_t = Next_indivs_mut_state_t_vector.at(indivs_mut_state_node_site->first);
//     indivs_mut_state_node_site->second += next_indivs_mut_state_t;
//     //next power of 2 (by movings bit by 1 in left)
//     Next_indivs_mut_state_t_vector.at(indivs_mut_state_node_site->first) = next_indivs_mut_state_t << 1;
// }

mut_mod_c_jcm::mut_mod_c_jcm(rand_gen_c *rand_gen)
{
    Rand_gen = rand_gen;

    Cumul_distrib_nucl_A = cumul_distrib({0, 1, 1, 1});
    Cumul_distrib_nucl_T = cumul_distrib({1, 0, 1, 1});
    Cumul_distrib_nucl_C = cumul_distrib({1, 1, 0, 1});
    Cumul_distrib_nucl_G = cumul_distrib({1, 1, 1, 0});
}

//Possibility of mutation in the same nucl
void mut_mod_c_jcm::apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site)
{
    switch (indivs_mut_state_node_site->second)
    {
    case 1:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_A, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    case 2:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_T, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    case 3:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_C, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    case 4:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_G, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    }
}

mut_mod_c_k80::mut_mod_c_k80(std::array<double, 2> ratio_transi_transver, rand_gen_c *rand_gen)
{
    Rand_gen = rand_gen;

    //TODO : regarder cette histoire de transi_transfer
    Cumul_distrib_nucl_A = cumul_distrib({0, 1, 1, ratio_transi_transver.at(0)});
    Cumul_distrib_nucl_T = cumul_distrib({1, 0, ratio_transi_transver.at(0), 1});
    Cumul_distrib_nucl_C = cumul_distrib({1, ratio_transi_transver.at(0), 0, 1});
    Cumul_distrib_nucl_G = cumul_distrib({ratio_transi_transver.at(0), 1, 1, 0});
}

void mut_mod_c_k80::apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site)
{
    switch (indivs_mut_state_node_site->second)
    {
    case 1:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_A, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    case 2:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_T, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    case 3:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_C, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    case 4:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_G, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    }
}

mut_mod_c_f81::mut_mod_c_f81(equi_base_freq_c equi_base_freq, rand_gen_c *rand_gen)
{
    Rand_gen = rand_gen;

    Cumul_distrib_nucl_A = cumul_distrib({0, equi_base_freq.T, equi_base_freq.C, equi_base_freq.G});
    Cumul_distrib_nucl_T = cumul_distrib({equi_base_freq.A, 0, equi_base_freq.C, equi_base_freq.G});
    Cumul_distrib_nucl_C = cumul_distrib({equi_base_freq.A, equi_base_freq.T, 0, equi_base_freq.G});
    Cumul_distrib_nucl_G = cumul_distrib({equi_base_freq.A, equi_base_freq.T, equi_base_freq.C, 0});
}

void mut_mod_c_f81::apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site)
{
    switch (indivs_mut_state_node_site->second)
    {
    case 1:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_A, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    case 2:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_T, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    case 3:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_C, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    case 4:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_G, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    }
}

mut_mod_c_hky::mut_mod_c_hky(std::array<double, 2> ratio_transi_transver, equi_base_freq_c equi_base_freq, rand_gen_c *rand_gen)
{
    Rand_gen = rand_gen;

    Cumul_distrib_nucl_A = cumul_distrib({0, equi_base_freq.T, equi_base_freq.C, ratio_transi_transver.at(0) * equi_base_freq.G});
    Cumul_distrib_nucl_T = cumul_distrib({equi_base_freq.A, 0, ratio_transi_transver.at(0) * equi_base_freq.C, equi_base_freq.G});
    Cumul_distrib_nucl_C = cumul_distrib({equi_base_freq.A, ratio_transi_transver.at(0) * equi_base_freq.T, 0, equi_base_freq.G});
    Cumul_distrib_nucl_G = cumul_distrib({ratio_transi_transver.at(0) * equi_base_freq.A, equi_base_freq.T, equi_base_freq.C, 0});
}

void mut_mod_c_hky::apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site)
{
    switch (indivs_mut_state_node_site->second)
    {
    case 1:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_A, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    case 2:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_T, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    case 3:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_C, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    case 4:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_G, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    }
}

//TODO rename tn933
mut_mod_c_tn93::mut_mod_c_tn93(std::array<double, 2> ratio_transi_transver, equi_base_freq_c equi_base_freq, rand_gen_c *rand_gen)
{
    Rand_gen = rand_gen;

    Cumul_distrib_nucl_A = cumul_distrib({0, equi_base_freq.T, equi_base_freq.C, ratio_transi_transver.at(0) * equi_base_freq.G});
    Cumul_distrib_nucl_T = cumul_distrib({equi_base_freq.A, 0, ratio_transi_transver.at(1) * equi_base_freq.C, equi_base_freq.G});
    Cumul_distrib_nucl_C = cumul_distrib({equi_base_freq.A, ratio_transi_transver.at(1) * equi_base_freq.T, 0, equi_base_freq.G});
    Cumul_distrib_nucl_G = cumul_distrib({ratio_transi_transver.at(0) * equi_base_freq.A, equi_base_freq.T, equi_base_freq.C, 0});
}

void mut_mod_c_tn93::apply_mut(std::map<int, int>::iterator indivs_mut_state_node_site)
{
    switch (indivs_mut_state_node_site->second)
    {
    case 1:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_A, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    case 2:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_T, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    case 3:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_C, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    case 4:
    {
        indivs_mut_state_node_site->second = nucl_distri(Cumul_distrib_nucl_G, Rand_gen->int_0_PRECISION_rand());
        break;
    }
    }
}

/**************************************************************/
/*                        Mut algo                            */
/**************************************************************/

mut_model_c::mut_model_c()
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto const &muta_param = singleton_c<muta_param_c>::instance();
    auto const &samp_param = singleton_c<samp_param_c>::instance();

    Mut_mod_iam = mut_mod_c_iam(samp_param.Sequence_length);
    Mut_mod_kam = mut_mod_c_kam(muta_param.K_min, muta_param.K_max, &rand_gen);
    Mut_mod_smm = mut_mod_c_smm(muta_param.K_min, muta_param.K_max, &rand_gen);
    Mut_mod_gsm = mut_mod_c_gsm(muta_param.K_min, muta_param.K_max, muta_param.P_gsm, &rand_gen);
    /*****/
    //Mut_mod_ism = mut_mod_c_ism(samp_param.Sequence_length);
    Mut_mod_jcm = mut_mod_c_jcm(&rand_gen);
    Mut_mod_k80 = mut_mod_c_k80(muta_param.Ratio_transi_transver, &rand_gen);
    Mut_mod_f81 = mut_mod_c_f81(muta_param.Equi_base_freq, &rand_gen);
    Mut_mod_hky = mut_mod_c_hky(muta_param.Ratio_transi_transver, muta_param.Equi_base_freq, &rand_gen);
    Mut_mod_tn93 = mut_mod_c_tn93(muta_param.Ratio_transi_transver, muta_param.Equi_base_freq, &rand_gen);
}

void mut_model_c::apply_mut_to_site(mut_model_enum mut_name, std::map<int, int>::iterator indivs_mut_state_node_site)
{
    switch (mut_name)
    {
    case mut_model_enum::iam:
    {
        Mut_mod_iam.apply_mut(indivs_mut_state_node_site);
        break;
    }
    case mut_model_enum::kam:
    {
        Mut_mod_kam.apply_mut(indivs_mut_state_node_site);
        break;
    }
    case mut_model_enum::smm:
    {
        Mut_mod_smm.apply_mut(indivs_mut_state_node_site);
        break;
    }
    case mut_model_enum::gsm:
    {
        Mut_mod_gsm.apply_mut(indivs_mut_state_node_site);
        break;
    }
    /*****/
    case mut_model_enum::jcm:
    {
        Mut_mod_jcm.apply_mut(indivs_mut_state_node_site);
        break;
    }
    case mut_model_enum::k80:
    {
        Mut_mod_k80.apply_mut(indivs_mut_state_node_site);
        break;
    }
    case mut_model_enum::f81:
    {
        Mut_mod_f81.apply_mut(indivs_mut_state_node_site);
        break;
    }
    case mut_model_enum::hky:
    {
        Mut_mod_hky.apply_mut(indivs_mut_state_node_site);
        break;
    }
    case mut_model_enum::tn93:
    {
        Mut_mod_tn93.apply_mut(indivs_mut_state_node_site);
        break;
    }
    }
}

#include "debug.hpp"
//Prepare the FIFO node_fifo and prepare the map to store the differents mut sets (one per node).
// FIFO = a queue (une fifo) = first in first out
//ancestry_seq not copy during mut process
sample_mutated_state_type &mut_genenerator_c::mut_generation(top_down_tree_t const &gene_tree, std::vector<int> const &ancestry_seq, std::array<int, 2> const &begin_end_seq_current_tree, int MRCA_index)
{
    auto const &info_collector = singleton_c<info_collector_c>::instance();
    auto const &samp_param = singleton_c<samp_param_c>::instance();

    if (info_collector.Check_tree_bool)
    {
        check_tree(gene_tree, MRCA_index, samp_param.n_total_sample_size * samp_param.Ploidy, info_collector.rep);
    }

    auto const &simu_param = singleton_c<simu_param_c>::instance(); // for Continuous_time_approxim

    Begin_end_sequence = begin_end_seq_current_tree;

    // mutated sites and states for each sample, for the current chunk
    All_indivs_mut_state = std::vector<std::map<int, int>>(MRCA_index + 1);

    std::vector<int> node_fifo;
    node_fifo.reserve(MRCA_index);

    node_fifo.push_back(MRCA_index);

    // each map will be "moved" in the vector from the parent to its last child
    // at the end, only the first 'sample_size' vector cell will content a map
    browse_tree(simu_param.Continuous_time_approxim, gene_tree, node_fifo, ancestry_seq);

    return All_indivs_mut_state;
}

//Browse tree in wide with an FIFO that allows to browse the nodes (more or less) from left to right

void mut_genenerator_c::browse_tree(bool approx, top_down_tree_t const &gene_tree, std::vector<int> &child_nodes, std::vector<int> const &ancestry_seq)
{
    auto const &samp_param = singleton_c<samp_param_c>::instance();
    auto &info_collector = singleton_c<info_collector_c>::instance();
    if (info_collector.Coa_times)
    {
        info_collector.Gen_coa_time.push_back(0); // ini for mean coal time between pairs of samples (leaves)
    }
    //While all nodes not browse
    while (!child_nodes.empty())
    {
        std::vector<int> current_nodes;
        //Max number of node in the same level
        current_nodes.reserve(samp_param.n_total_sample_size * samp_param.Ploidy);

        child_nodes.swap(current_nodes);
        //after swap,
        //current_nodes contains all nodes at this level of the tree
        // and child_nodes is empty and will later contains all children nodes from all nodes at this level
        for (int ancest_index : current_nodes)
        {
            std::vector<int> children = std::get<0>(gene_tree.at(ancest_index));

            int nbr_of_pair = 0;
            if (info_collector.Coa_times)
            {
                if (std::get<3>(gene_tree.at(ancest_index)) > 1)
                {
                    nbr_of_pair = combination(2, std::get<3>(gene_tree.at(ancest_index)));
                }
            }

            //Node's tree have severals childs or leaves!
            auto last_child = --children.end();
            if (!children.empty()) // test if it is a leaf
            {
                // for all children, except last one (for which the map is moved to below), copy the mut_state info and apply mut
                for (auto child = children.begin(); child != last_child; ++child)
                {
                    child_nodes.push_back(*child);
                    All_indivs_mut_state[*child] = All_indivs_mut_state.at(ancest_index);
                    mut_process(approx, gene_tree, ancest_index, *child, ancestry_seq);
                    //Information
                    if (info_collector.Coa_times)
                    {
                        if (std::get<3>(gene_tree.at(*child)) > 1)
                        {
                            nbr_of_pair -= combination(2, std::get<3>(gene_tree[*child]));
                        }
                    }
                }
                //last child map can be move
                child_nodes.push_back(*last_child);
                All_indivs_mut_state.at(*last_child) = std::move(All_indivs_mut_state.at(ancest_index));
                mut_process(approx, gene_tree, ancest_index, *last_child, ancestry_seq);
                //Information
                if (info_collector.Coa_times)
                {
                    if (std::get<3>(gene_tree.at(*last_child)) > 1)
                    {
                        nbr_of_pair -= combination(2, std::get<3>(gene_tree.at(*last_child)));
                    }
                }
            }
            if (info_collector.Coa_times)
            {
                if (approx)
                {
                    info_collector.Gen_coa_time.back() += (nbr_of_pair * std::get<1>(gene_tree.at(ancest_index)));
                }
                else
                {
                    info_collector.Gen_coa_time.back() += (nbr_of_pair * std::get<2>(gene_tree.at(ancest_index)));
                }
            }
        }
    }
    //Coa mean time between 2 lineage
    if (info_collector.Coa_times)
    {
        if (samp_param.n_total_sample_size * samp_param.Ploidy > 1)
        {
            info_collector.Gen_coa_time.back() /= combination(2, samp_param.n_total_sample_size * samp_param.Ploidy);
        }
    }
}

//Calculates the branch length and apply mutations
void mut_genenerator_c::mut_process(bool approx, top_down_tree_t const &gene_tree, int ancester_index, int node_index, std::vector<int> const &ancestry_seq)
{
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto const &muta_param = singleton_c<muta_param_c>::instance();

    double long_left_branch;

    //Hudson algorithm uses relative time between node
    if (approx)
    {
        long_left_branch = std::get<1>(gene_tree[ancester_index]) - std::get<1>(gene_tree[node_index]);
        mut_events(approx, node_index, long_left_branch, muta_param.Scaled_mut_rate_theta, ancestry_seq, rand_gen);
    }
    //Gen by gen algo uses nbr of gen
    else
    {
        long_left_branch = std::get<2>(gene_tree[ancester_index]) - std::get<2>(gene_tree[node_index]);
        mut_events(approx, node_index, long_left_branch, muta_param.Unscaled_mut_rate_mu, ancestry_seq, rand_gen);
    }
}

//Draws in a poisson/binomial distribution the number of muts that will occur in that branch and uniformely distributes them across sites

void mut_genenerator_c::mut_events(bool approx, int node_index, double long_branch, double mut_rate, std::vector<int> const &ancestry_seq, rand_gen_c &rand_gen)
{
    auto const &muta_param = singleton_c<muta_param_c>::instance();
    std::map<int, int> &indivs_mut_state_node = All_indivs_mut_state[node_index];

    auto site_nbr = Begin_end_sequence.at(1) - Begin_end_sequence.at(0);

    int mut_nbr;
    //uses local implementation of random distributions (cstl)
    //Binomial distribution need a number of "trial" so can only be use in gen by gen algo
    if (approx)
    {
        //Binomial distribution can be approxim by a poisson_distrib if mut_rate is small
        std::poisson_distribution<int> poisson_distrib(mut_rate * long_branch * site_nbr);
        mut_nbr = poisson_distrib(rand_gen.Seed_gen);
    }
    else
    {
        std::binomial_distribution<long int> binomial_distrib(long_branch * site_nbr, mut_rate);
        mut_nbr = binomial_distrib(rand_gen.Seed_gen);
    }

    //To choose mutated sites
    //Use -1 in end because sequence of 10 sites is implemented as [0,10[
    std::uniform_int_distribution<int> uniform_distribution(Begin_end_sequence.at(0), Begin_end_sequence.at(1) - 1);
    int site_index;

    while (mut_nbr != 0)
    {
        site_index = uniform_distribution(rand_gen.Seed_gen);
        // std::map<int, int>::iterator : look it the map if the site is present (already mutated) otherwise return the end of the map
        auto indivs_mut_state_site = indivs_mut_state_node.find(site_index);
        //if indivs_mut_state_site doesn't yet exist
        if (indivs_mut_state_site == indivs_mut_state_node.end())
        {
            //Emplace return a pair consisting of an iter to the inserted element and a bool.
            indivs_mut_state_site = indivs_mut_state_node.emplace(site_index, ancestry_seq[site_index]).first;
        }

        //TODO chercher moyen de ne pas faire le switch a chaque mut
        Mut_model.apply_mut_to_site(muta_param.Mod_mut_name, indivs_mut_state_site);
        --mut_nbr;
    }
}
#include <chrono>
//TODO : Not covered by unit test
std::vector<std::vector<std::pair<int, int>>> apply_mut_to_sample(coa_table_c &coa_table, int next_max_node_nbr, std::vector<int> const &ancestry_seq, muta_param_c const &muta_param, samp_param_c const &samp_param, demo_param_c const &demo_param)
{
    auto &info_collect = singleton_c<info_collector_c>::instance();

    std::vector<std::vector<std::pair<int, int>>> sample_mutated_state(samp_param.n_total_sample_size * samp_param.Ploidy);
    for (auto &mut_stat : sample_mutated_state)
    {
        //reserve => nbr mut / gen * Theo mean MRCA
        mut_stat.reserve(muta_param.Unscaled_mut_rate_mu * samp_param.Sequence_length * (2 * demo_param.Population_size_N * samp_param.Ploidy) * (1 - (1.0 / samp_param.n_total_sample_size * samp_param.Ploidy)));
    }

    mut_genenerator_c mut_generator;

    //[begin, end[ non recombinant chunk
    std::tuple<int, int, double> begin_end_MRCAtime_non_recomb_chunk;
    //tree_gen constructor
    tree_gen_c tree_gen(coa_table, next_max_node_nbr, samp_param.n_total_sample_size * samp_param.Ploidy);
    /*************************/
    do
    {
        std::chrono::time_point<std::chrono::high_resolution_clock> debut;
        if (info_collect.Clock)
        {
            debut = std::chrono::high_resolution_clock::now();
        }

        begin_end_MRCAtime_non_recomb_chunk = tree_gen.set_next_tree();

        if (info_collect.Clock)
        {
            auto fin = std::chrono::high_resolution_clock::now();
            info_collect.time_construct_tree += (std::chrono::duration_cast<std::chrono::nanoseconds>(fin - debut).count() * 0.000000001);
        }

        if (info_collect.Coa_times)
        {
            info_collect.MRCA_time.push_back(std::get<2>(begin_end_MRCAtime_non_recomb_chunk));
        }

        std::array<int, 2> begin_end_seq = {std::get<0>(begin_end_MRCAtime_non_recomb_chunk), std::get<1>(begin_end_MRCAtime_non_recomb_chunk)};

        if (info_collect.Clock)
        {
            debut = std::chrono::high_resolution_clock::now();
        }
        //list of chunk for each vector
        auto temp = mut_generator.mut_generation(tree_gen.get_tree(), ancestry_seq, begin_end_seq, tree_gen.MRCA_node_nbr);

        auto temp_itr = temp.begin();
        for (auto &mut_stat : sample_mutated_state)
        {
            for (auto const &mut : *temp_itr)
            {
                mut_stat.emplace_back(mut);
            }
            ++temp_itr;
        }

        if (info_collect.Clock)
        {
            auto fin = std::chrono::high_resolution_clock::now();
            info_collect.time_mutation += (std::chrono::duration_cast<std::chrono::nanoseconds>(fin - debut).count() * 0.000000001);
        }

    } while (std::get<1>(begin_end_MRCAtime_non_recomb_chunk) < samp_param.Sequence_length);

    for (auto &mut_stat : sample_mutated_state)
    {
        mut_stat.shrink_to_fit();
    }

    return sample_mutated_state;
}
