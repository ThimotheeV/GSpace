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

int mut_mod_c_iam::apply_mut(int site_index)
{
    return ++Count.at(site_index);
}

mut_mod_c_kam::mut_mod_c_kam(int k_min, int k_max, rand_gen_c *rand_gen)
{
    Rand_gen = rand_gen;

    Uni_int_kmin_kmax_distrib = std::uniform_int_distribution<>(k_min, k_max);
}

int mut_mod_c_kam::apply_mut(int state)
{
    auto temp = Uni_int_kmin_kmax_distrib(Rand_gen->Seed_gen);
    while (temp == state)
    {
        temp = Uni_int_kmin_kmax_distrib(Rand_gen->Seed_gen);
    }
    return temp;
}

mut_mod_c_smm::mut_mod_c_smm(int k_min, int k_max, rand_gen_c *rand_gen)
{
    Rand_gen = rand_gen;

    K_min = k_min;
    K_max = k_max;
}

int mut_mod_c_smm::apply_mut(int state)
{
    if (Rand_gen->rand_bool())
    {
        if (state != K_max)
        {
            state += 1;
        }
    }
    else
    {
        if (state != K_min)
        {
            state -= 1;
        }
    }
    return state;
}

mut_mod_c_gsm::mut_mod_c_gsm(int k_min, int k_max, double p_gsm, rand_gen_c *rand_gen)
{
    Rand_gen = rand_gen;

    K_min = k_min;
    K_max = k_max;

    //use q = 1-p instead off p in geometric distribution
    Geo_distrib = std::geometric_distribution<>(p_gsm); // <> = default template type = int
}

int mut_mod_c_gsm::apply_mut(int state)
{
    auto repet = Geo_distrib(Rand_gen->Seed_gen) + 1;

    if (Rand_gen->rand_bool())
    {
        if (state + repet <= K_max)
        {
            state += repet;
        }
    }
    else
    {
        if (state - repet >= K_min)
        {
            state -= repet;
        }
    }
    return state;
}

// mut_mod_c_ism::mut_mod_c_ism(int site_nbr)
// {
//     // value for next mutated state, ini 1
//     Next_indivs_mut_state_t_vector = std::vector<int>(site_nbr, 1);
// }

// int mut_mod_c_ism::apply_mut(int state)
// {
//     //Keep in memory the next "mutation number representation" (power of 2) for the site.
//     auto next_indivs_mut_state_t = Next_indivs_mut_state_t_vector.at(indivs_mut_state_node_site->first);
//     state += next_indivs_mut_state_t;
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
int mut_mod_c_jcm::apply_mut(int state)
{
    switch (state)
    {
    case 1:
    {
        return nucl_distri(Cumul_distrib_nucl_A, Rand_gen->int_0_PRECISION_rand());
    }
    case 2:
    {
        return nucl_distri(Cumul_distrib_nucl_T, Rand_gen->int_0_PRECISION_rand());
    }
    case 3:
    {
        return nucl_distri(Cumul_distrib_nucl_C, Rand_gen->int_0_PRECISION_rand());
    }
    case 4:
    {
        return nucl_distri(Cumul_distrib_nucl_G, Rand_gen->int_0_PRECISION_rand());
    }
    default:
    {
        return -1;
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

int mut_mod_c_k80::apply_mut(int state)
{
    switch (state)
    {
    case 1:
    {
        return nucl_distri(Cumul_distrib_nucl_A, Rand_gen->int_0_PRECISION_rand());
    }
    case 2:
    {
        return nucl_distri(Cumul_distrib_nucl_T, Rand_gen->int_0_PRECISION_rand());
    }
    case 3:
    {
        return nucl_distri(Cumul_distrib_nucl_C, Rand_gen->int_0_PRECISION_rand());
    }
    case 4:
    {
        return nucl_distri(Cumul_distrib_nucl_G, Rand_gen->int_0_PRECISION_rand());
    }
    default:
    {
        return -1;
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

int mut_mod_c_f81::apply_mut(int state)
{
    switch (state)
    {
    case 1:
    {
        return nucl_distri(Cumul_distrib_nucl_A, Rand_gen->int_0_PRECISION_rand());
    }
    case 2:
    {
        return nucl_distri(Cumul_distrib_nucl_T, Rand_gen->int_0_PRECISION_rand());
    }
    case 3:
    {
        return nucl_distri(Cumul_distrib_nucl_C, Rand_gen->int_0_PRECISION_rand());
    }
    case 4:
    {
        return nucl_distri(Cumul_distrib_nucl_G, Rand_gen->int_0_PRECISION_rand());
    }
    default:
    {
        return -1;
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

int mut_mod_c_hky::apply_mut(int state)
{
    switch (state)
    {
    case 1:
    {
        return nucl_distri(Cumul_distrib_nucl_A, Rand_gen->int_0_PRECISION_rand());
    }
    case 2:
    {
        return nucl_distri(Cumul_distrib_nucl_T, Rand_gen->int_0_PRECISION_rand());
    }
    case 3:
    {
        return nucl_distri(Cumul_distrib_nucl_C, Rand_gen->int_0_PRECISION_rand());
    }
    case 4:
    {
        return nucl_distri(Cumul_distrib_nucl_G, Rand_gen->int_0_PRECISION_rand());
    }
    default:
    {
        return -1;
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

int mut_mod_c_tn93::apply_mut(int state)
{
    switch (state)
    {
    case 1:
    {
        return nucl_distri(Cumul_distrib_nucl_A, Rand_gen->int_0_PRECISION_rand());
    }
    case 2:
    {
        return nucl_distri(Cumul_distrib_nucl_T, Rand_gen->int_0_PRECISION_rand());
    }
    case 3:
    {
        return nucl_distri(Cumul_distrib_nucl_C, Rand_gen->int_0_PRECISION_rand());
    }
    case 4:
    {
        return nucl_distri(Cumul_distrib_nucl_G, Rand_gen->int_0_PRECISION_rand());
    }
    default:
    {
        return -1;
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

int mut_model_c::apply_mut_to_site(mut_model_enum mut_name, int state, int site_index)
{
    switch (mut_name)
    {
    case mut_model_enum::iam:
    {
        return Mut_mod_iam.apply_mut(site_index);
    }
    case mut_model_enum::kam:
    {
        return Mut_mod_kam.apply_mut(state);
    }
    case mut_model_enum::smm:
    {
        return Mut_mod_smm.apply_mut(state);
    }
    case mut_model_enum::gsm:
    {
        return Mut_mod_gsm.apply_mut(state);
    }
    /*****/
    case mut_model_enum::jcm:
    {
        return Mut_mod_jcm.apply_mut(state);
    }
    case mut_model_enum::k80:
    {
        return Mut_mod_k80.apply_mut(state);
    }
    case mut_model_enum::f81:
    {
        return Mut_mod_f81.apply_mut(state);
    }
    case mut_model_enum::hky:
    {
        return Mut_mod_hky.apply_mut(state);
    }
    case mut_model_enum::tn93:
    {
        return Mut_mod_tn93.apply_mut(state);
    }
    default:
    {
        return -1;
    }
    }
}

#include "debug.hpp"
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

        if (info_collect.MRCA_record)
        {
            info_collect.MRCA_time.push_back(std::get<2>(begin_end_MRCAtime_non_recomb_chunk));
        }
        std::array<int, 2> begin_end_seq = {std::get<0>(begin_end_MRCAtime_non_recomb_chunk), std::get<1>(begin_end_MRCAtime_non_recomb_chunk)};
        // std::cout<<std::get<0>(begin_end_MRCAtime_non_recomb_chunk)<<" "<<std::get<1>(begin_end_MRCAtime_non_recomb_chunk)<<std::endl;
        if (info_collect.Clock)
        {
            debut = std::chrono::high_resolution_clock::now();
        }
        //list of chunk for each vector
        auto temp = mut_generator.mut_generation(tree_gen.get_tree(), ancestry_seq, begin_end_seq, muta_param, tree_gen.MRCA_node_nbr, samp_param.n_total_sample_size * samp_param.Ploidy);

        auto temp_itr = temp.begin();
        for (auto &mut_stat : sample_mutated_state)
        {
            mut_stat.reserve(mut_stat.size() + temp_itr->size());
            for (auto const &value : *temp_itr)
            {
                mut_stat.push_back(value);
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

//Browse tree in deep, create markov chain of mut for each site and sample
sample_mutated_state_type &mut_genenerator_c::mut_generation(top_down_tree_t const &gene_tree, std::vector<int> const &ancestry_seq, std::array<int, 2> const &begin_end_seq_current_tree, muta_param_c const &muta_param, int MRCA_index, int sample_size)
{
    auto const &simu_param = singleton_c<simu_param_c>::instance();
    auto const &samp_param = singleton_c<samp_param_c>::instance();
    auto &rand_gen = singleton_c<rand_gen_c>::instance();

    auto const &info_collector = singleton_c<info_collector_c>::instance();

    if (info_collector.Check_tree_bool)
    {
        check_tree(gene_tree, MRCA_index, samp_param.n_total_sample_size * samp_param.Ploidy, info_collector.rep);
    }
    //Ini
    bool descend = true;
    long int nbr_site = std::get<1>(begin_end_seq_current_tree) - std::get<0>(begin_end_seq_current_tree);
    Mod_mut_name = muta_param.Mod_mut_name;

    //Draw a site
    Site_distri = std::uniform_int_distribution<long int>(std::get<0>(begin_end_seq_current_tree), std::get<1>(begin_end_seq_current_tree) - 1);
    //For each site, the markov chain mutation
    Mut_chain = std::map<int, std::list<int>>();
    //Recap all mut emplacement by branch
    Mut_register = std::list<std::list<int>>{};

    // mutated sites and states for each sample, for the current chunk
    All_indivs_mut_state = sample_mutated_state_type(sample_size);

    int last_node = -1;

    //Recap of the browsing branch (chain of parent node)
    std::vector<int> ligneage;
    ligneage.reserve(sample_size);
    ligneage.push_back(MRCA_index);

    while (!ligneage.empty())
    {
        int parent_index = ligneage.back();
        auto const &branch = gene_tree.at(parent_index);
        if (descend)
        {
            //leaf = ligneage.back() and return to previous node next turn => ligneage.pop_back();
            //If no child
            if (std::get<0>(branch).empty())
            {
                descend = false;
                ligneage.pop_back();
                /********************/
                if (!Mut_register.empty())
                {
                    harvest_mut(parent_index);
                    remove_mut();
                }
            }
            else
            {
                //internal node go to next node => ligneage.push_back(child);
                int child = std::get<0>(branch).back();
                ligneage.push_back(child);
                /********************/
                int nbr_mut = draw_nbr_mut(muta_param, rand_gen, branch, gene_tree.at(child), nbr_site, simu_param.Continuous_time_approxim);
                //need to add all node
                Mut_register.emplace_back(std::list<int>{});
                if (nbr_mut != 0)
                {
                    add_mut(ancestry_seq, rand_gen, child, nbr_mut);
                }
            }
        }
        else
        {
            //have already visit this node, search the first child != than previous
            //if exist go to it, if not destroy and return ta parent
            bool no_more_new_child = true;
            auto rev_itr = std::get<0>(branch).rbegin();
            while (no_more_new_child && rev_itr != std::get<0>(branch).rend())
            {
                if (*rev_itr == last_node)
                {
                    no_more_new_child = false;
                }
                ++rev_itr;
            }
            /********************/

            //rev_itr = null => internal node with no more child = delete
            if (no_more_new_child || rev_itr == std::get<0>(branch).rend())
            {
                ligneage.pop_back();
                /********************/
                if (!Mut_register.empty())
                {
                    remove_mut();
                }
            }
            //rev_itr = child => internal node with child, child not mutate
            else
            {
                descend = true;
                int child = *rev_itr;

                ligneage.push_back(child);
                /********************/
                int nbr_mut = draw_nbr_mut(muta_param, rand_gen, branch, gene_tree.at(child), nbr_site, simu_param.Continuous_time_approxim);
                Mut_register.emplace_back(std::list<int>{});
                if (nbr_mut != 0)
                {
                    add_mut(ancestry_seq, rand_gen, child, nbr_mut);
                }
            }
        }
        //needed for mutation
        last_node = parent_index;
    }

    return All_indivs_mut_state;
}

void mut_genenerator_c::add_mut(std::vector<int> const &ancester_state, rand_gen_c &seed_gen, int node_index, int nbr_mut)
{
    auto &tskit = singleton_c<tskit_struct_c>::instance();

    std::vector<long int> sites_index;
    //+1 for efficient sorting
    sites_index.reserve(nbr_mut + 1);

    for (int mut = 0; mut < nbr_mut; ++mut)
    {
        sites_index.push_back(Site_distri(seed_gen.Seed_gen));
    }
    std::sort(sites_index.begin(), sites_index.end());

    sites_index.shrink_to_fit();

    std::list<int> &node_mut_register = Mut_register.back();

    int prev_index = -1;
    std::list<int> *mut_chain_site = nullptr;
    for (int index_site : sites_index)
    {
        //Sorting site => cluster multi mut
        if (index_site != prev_index)
        {
            //first mut for index_site for this node
            auto search = Mut_chain.find(index_site);
            if (search == Mut_chain.end())
            {
                //first mut for index_site for all the tree
                //emplace : Returns a pair consisting of an iterator to the inserted element and a bool denoting whether the insertion took place
                search = Mut_chain.emplace(index_site, std::list<int>{ancester_state.at(index_site)}).first;
                mut_chain_site = &search->second;
            }
            else
            {
                //first mut for index_site for this node = take previous state
                mut_chain_site = &search->second;
                mut_chain_site->push_back(mut_chain_site->back());
            }
            //index_site mut in this node
            node_mut_register.push_back(index_site);
        }
        //apply mut
        mut_chain_site->back() = Mut_model.apply_mut_to_site(Mod_mut_name, mut_chain_site->back(), index_site);
        if (tskit.Tskit_output)
        {
            tskit.add_new_mut(index_site, node_index, std::to_string(mut_chain_site->back()));
        }
        prev_index = index_site;
    }
}

void mut_genenerator_c::remove_mut()
{
    for (int index_site : Mut_register.back())
    {
        //find : Return iterator to an element with key equivalent to key
        auto chain = Mut_chain.find(index_site);
        chain->second.pop_back();
        if (chain->second.empty())
        {
            Mut_chain.erase(chain);
        }
    }
    Mut_register.pop_back();
}

//one sample
void mut_genenerator_c::harvest_mut(int num_sample)
{
    for (auto const &state : Mut_chain)
    {
        All_indivs_mut_state.at(num_sample).emplace_back(state.first, (state.second).back());
    }
}

int mut_genenerator_c::draw_nbr_mut(muta_param_c const &muta_param, rand_gen_c &rand_gen, std::tuple<std::vector<int>, double, int> const &parent_node, std::tuple<std::vector<int>, double, int> const &child_node, long int nbr_site, bool approx)
{
    int nbr_mut = 0;
    if (approx)
    {
        double long_branch = std::get<1>(parent_node) - std::get<1>(child_node);
        std::poisson_distribution<long int> poisson_distrib(muta_param.Scaled_mut_rate_theta * long_branch * nbr_site);
        nbr_mut = poisson_distrib(rand_gen.Seed_gen);
    }
    else
    {
        int long_branch = std::get<2>(parent_node) - std::get<2>(child_node);
        long long int nbr_trial = nbr_site * long_branch;
        if (nbr_trial > 10000)
        {
            std::poisson_distribution<long long int> nbr_of_mut(nbr_trial * muta_param.Unscaled_mut_rate_mu);
            nbr_mut = nbr_of_mut(rand_gen.Seed_gen);
        }
        else
        {
            std::binomial_distribution<int> nbr_of_mut(nbr_trial, muta_param.Unscaled_mut_rate_mu);
            nbr_mut = nbr_of_mut(rand_gen.Seed_gen);
        }
    }

    return nbr_mut;
}