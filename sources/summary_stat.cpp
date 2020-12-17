#include <tuple>
#include <numeric>
#include <map>
#include <cmath>
#include <iostream>

#include "summary_stat.hpp"
#include "mutation.hpp" // TODO a virer quand code ci dessous déplacé dans mutation

//TODO a déplacer dans mutation.hpp/cpp
std::vector<int> reconstitution_full_int_seq_indiv(std::vector<int> ancestry_seq, std::vector<std::pair<int, int>> const &indiv_mut_seq)
{
    for (auto const &mut : indiv_mut_seq)
    {
        ancestry_seq[mut.first] = mut.second;
    }

    return ancestry_seq;
}

//TODO a déplacer dans mutation.hpp/cpp
std::string reconstitution_full_nucl_seq_indiv(std::vector<int> ancestry_seq, std::vector<std::pair<int, int>> const &indiv_mut_seq)
{
    std::string nucleotide_sequence;

    for (auto const &mut : indiv_mut_seq)
    {
        ancestry_seq[mut.first] = mut.second;
    }
    map_string_nuc_name_c map_string_nuc_name;

    for (auto const &int_site : ancestry_seq)
    {
        nucleotide_sequence.append(map_string_nuc_name.Inverse_map.at(int_site));
    }
    return nucleotide_sequence;
}

void Prob_id_1_loc_Qr_res_c::calc_Qr(data_plane_vec_c const &data_plane_vec, int rep)
{
    auto const &samp_param = singleton_c<samp_param_c>::instance();

    Qr_by_chr_by_loc = std::vector<std::vector<std::vector<double>>>(samp_param.Chr_nbr, std::vector<std::vector<double>>(samp_param.Sequence_length, std::vector<double>(data_plane_vec.nbr_of_dist_class())));
    Qr_by_chr = std::vector<std::vector<double>>(samp_param.Chr_nbr, std::vector<double>(data_plane_vec.nbr_of_dist_class()));
    Qr = std::vector<double>(data_plane_vec.nbr_of_dist_class());

    if (Qr_cumul_m_v.empty())
    {
        Qr_cumul_m_v = std::vector<std::array<double, 2>>(data_plane_vec.nbr_of_dist_class(), {0, 0});
    }

    std::vector<std::vector<double>> temp_for_all(data_plane_vec.nbr_of_dist_class(), std::vector<double>(samp_param.Chr_nbr));
    for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
    {
        std::vector<std::vector<double>> temp_for_chr(data_plane_vec.nbr_of_dist_class(), std::vector<double>(samp_param.Sequence_length));
        for (int locus = 0; locus < samp_param.Sequence_length; ++locus)
        {
            //locus + samp_param.Sequence_length * chr to browse one locus in each chr (and not chr time the first chr)
            auto temp_qr = calc_qr_loc_by_loc(data_plane_vec, locus + samp_param.Sequence_length * chr);
            for (int dist = 0; dist < data_plane_vec.nbr_of_dist_class(); ++dist)
            {
                Qr_by_chr_by_loc[chr][locus][dist] = static_cast<double>(temp_qr[dist].at(0)) / temp_qr[dist].at(1);
                temp_for_chr[dist][locus] = Qr_by_chr_by_loc[chr][locus][dist];
            }
        }
        for (int dist = 0; dist < data_plane_vec.nbr_of_dist_class(); ++dist)
        {
            Qr_by_chr[chr][dist] = mean(temp_for_chr[dist]);
            temp_for_all[dist][chr] = Qr_by_chr[chr][dist];
        }
    }
    for (int dist = 0; dist < data_plane_vec.nbr_of_dist_class(); ++dist)
    {
        Qr[dist] = mean(temp_for_all[dist]);
        Qr_cumul_m_v[dist] = calc_cumul_mean_var(Qr_cumul_m_v[dist], {Qr[dist]}, rep);
    }
}

void Prob_id_1_loc_Qwi_wd_bd_res_c::calc_Qwi(data_plane_vec_c const &data_plane_vec, int rep)
{
    auto const &samp_param = singleton_c<samp_param_c>::instance();

    Qwi_by_chr_by_loc = std::vector<std::vector<double>>(samp_param.Chr_nbr, std::vector<double>(samp_param.Sequence_length));
    Qwi_by_chr = std::vector<double>(samp_param.Chr_nbr);

    std::vector<double> temp_for_all(samp_param.Chr_nbr);
    for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
    {
        std::vector<double> temp_for_chr(samp_param.Sequence_length);
        for (int locus = 0; locus < samp_param.Sequence_length; ++locus)
        {
            //Qwi by loc, //locus + samp_param.Sequence_length * chr to browse one locus in each chr (and not chr time the first chr)
            auto temp = calc_Q_intra_indiv_per_locus(data_plane_vec, locus + samp_param.Sequence_length * chr);
            Qwi_by_chr_by_loc[chr][locus] = static_cast<double>(temp.at(0)) / temp.at(1);
            temp_for_chr[locus] = Qwi_by_chr_by_loc[chr][locus];
        }
        //Qwi by chr
        Qwi_by_chr[chr] = mean(temp_for_chr);
        temp_for_all[chr] = Qwi_by_chr[chr];
    }
    //Qwi by chr
    Qwi = mean(temp_for_all);
    Qwi_cumul_m_v = calc_cumul_mean_var(Qwi_cumul_m_v, {Qwi}, rep);
}

void Prob_id_1_loc_Qwi_wd_bd_res_c::calc_Qwd_bd(data_plane_vec_c const &data_plane_vec, int rep)
{
    auto const &samp_param = singleton_c<samp_param_c>::instance();

    Qwd_by_chr_by_loc = std::vector<std::vector<double>>(samp_param.Chr_nbr, std::vector<double>(samp_param.Sequence_length));
    Qbd_by_chr_by_loc = std::vector<std::vector<double>>(samp_param.Chr_nbr, std::vector<double>(samp_param.Sequence_length));

    Qwd_by_chr = std::vector<double>(samp_param.Chr_nbr);
    Qbd_by_chr = std::vector<double>(samp_param.Chr_nbr);

    std::vector<double> temp_for_all_Qwd(samp_param.Chr_nbr);
    std::vector<double> temp_for_all_Qbd(samp_param.Chr_nbr);

    for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
    {
        std::vector<double> temp_for_chr_Qwd(samp_param.Sequence_length);
        std::vector<double> temp_for_chr_Qbd(samp_param.Sequence_length);
        for (int locus = 0; locus < samp_param.Sequence_length; ++locus)
        {
            //locus + samp_param.Sequence_length * chr to browse one locus in each chr (and not chr time the first chr)
            auto temp_Qwd = calc_Q_inter_indiv_per_locus(data_plane_vec, locus + samp_param.Sequence_length * chr);
            auto temp_Qbd = calc_Q_inter_deme_per_locus(data_plane_vec, locus + samp_param.Sequence_length * chr);

            Qwd_by_chr_by_loc[chr][locus] = static_cast<double>(temp_Qwd.at(0)) / temp_Qwd.at(1);
            temp_for_chr_Qwd[locus] = Qwd_by_chr_by_loc[chr][locus];

            Qbd_by_chr_by_loc[chr][locus] = static_cast<double>(temp_Qbd.at(0)) / temp_Qbd.at(1);
            temp_for_chr_Qbd[locus] = Qbd_by_chr_by_loc[chr][locus];
        }

        Qwd_by_chr[chr] = mean(temp_for_chr_Qwd);
        temp_for_all_Qwd[chr] = Qwd_by_chr[chr];

        Qbd_by_chr[chr] = mean(temp_for_chr_Qbd);
        temp_for_all_Qbd[chr] = Qbd_by_chr[chr];
    }
    //Qwd by chr
    Qwd = mean(temp_for_all_Qwd);
    Qwd_cumul_m_v = calc_cumul_mean_var(Qwd_cumul_m_v, {Qwd}, rep);

    Qbd = mean(temp_for_all_Qbd);
    Qbd_cumul_m_v = calc_cumul_mean_var(Qbd_cumul_m_v, {Qbd}, rep);
}

void Prob_id_1_2_loc_res_c::cumul_prob_id_1_2_loc(data_plane_vec_c const &data_plane_vec, int rep)
{
    auto const &samp_param = singleton_c<samp_param_c>::instance();

    auto temp = calc_phi_ij(data_plane_vec, samp_param.Ploidy);
    PHI = mean(temp);
    PHI_cumul_m_v = calc_cumul_mean_var(PHI_cumul_m_v, temp, rep);
}

#include <algorithm>

void freq_allele_c::metering_mut_for_frequency(sample_mutated_state_type const &sample_mutated_state)
{
    std::vector<std::vector<int>> muts_by_site(Ancestry_seq.size());
    //Browse each indiv_c
    for (auto const &indiv_c : sample_mutated_state)
    { //Browse each mut in indiv_c and regroup them by site
        for (auto const &mut : indiv_c)
        {
            muts_by_site[mut.first].push_back(mut.second);
        }
    }

    int num_site = -1;
    //For each site sort the mutation and count each
    for (auto &site : muts_by_site)
    {
        std::sort(site.begin(), site.end());

        //Need to know if the next mutation is the same or a new one
        int prev_mut = -1;
        ++num_site;
        //Can have more mutation in a site than the greatest number representation
        if (!site.empty())
        {
            Mut_count_by_site[num_site].reserve(site.back());
        }

        for (auto const mut : site)
        {
            //for mut is new = create a new array with the mutation number and the number of mutation
            if (mut != prev_mut)
            {
                Mut_count_by_site[num_site].push_back({mut, 1});
            }
            else
            {
                ++(Mut_count_by_site[num_site].back()).at(1);
            }
            prev_mut = mut;
        }
        //Reduce vector reserve (in memory) at vector size (number of element)
        Mut_count_by_site[num_site].shrink_to_fit();
    }
}

std::array<double, 2> calc_cumul_mean_var(std::array<double, 2> cumul_mean_var_pond, std::vector<double> const &vec_value, int pond)
{
    auto local_mean = mean(vec_value);
    auto local_mean_var = std::array<double, 2>{local_mean, var(vec_value, local_mean)};
    double ratio = (static_cast<double>(1) / (pond + 1));

    //For new variance we can calculate separetly mean and value because var = pond*(1/n)*SUM[1,(p+1)*n](val^2)-mean^2 where n is the number of value
    double cumul_var_adapt = cumul_mean_var_pond.at(1) + std::pow(cumul_mean_var_pond.at(0), 2);
    double local_var_adapt = local_mean_var.at(1) + std::pow(local_mean_var.at(0), 2);

    cumul_mean_var_pond.at(0) = ratio * (pond * cumul_mean_var_pond.at(0) + local_mean_var.at(0));
    cumul_mean_var_pond.at(1) = ratio * (pond * cumul_var_adapt + local_var_adapt) - std::pow(cumul_mean_var_pond.at(0), 2);

    return cumul_mean_var_pond;
}

void coa_tree_metrics_c::calcul_coa_tree_metrics(bool approx, int ploidy, int population_size, int sample_size)
{
    //WARNING : In this prog by construction the indiv_c in pop are always diploid
    int gene_population_size;
    if (approx)
    {
        gene_population_size = 1;
    }
    else
    {
        gene_population_size = 2 * population_size;
    }

    int gene_sample_size = ploidy * sample_size;

    Theo_MRCA_mean = (2 * gene_population_size) * (1 - (1.0 / gene_sample_size));
    Theo_MRCA_var = pow(2 * gene_population_size, 2) * (1 - (1.0 / (gene_sample_size)));

    double harmonic_number_of_sample_size = 0;
    double harmonic_number_of_square_sample_size = 0;

    for (int i = 1; i < gene_sample_size; ++i)
    {
        harmonic_number_of_sample_size += 1.0 / i;
        harmonic_number_of_square_sample_size += 1.0 / std::pow(i, 2);
    }

    Theo_2_lign_coa_time_mean = gene_population_size;
    Theo_2_lign_coa_time_var = (std::pow(2 * gene_population_size, 2) * harmonic_number_of_square_sample_size) - std::pow(2 * gene_population_size * harmonic_number_of_sample_size, 2);
}
