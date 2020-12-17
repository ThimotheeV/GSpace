#include "data_plane_vec_lib.hpp"

// std::vector<std::vector<std::vector<std::pair<int, int>>>> => chr,indiv,mut_site
data_plane_vec_c::data_plane_vec_c(std::vector<std::vector<std::vector<std::pair<int, int>>>> const &sample_mutated_state_chr, samp_param_c const &samp_param, std::vector<std::vector<int>> const &ancestry_seq)
{
    Ploidy = samp_param.Ploidy;
    Dist_class_nbr = samp_param.Dist_class_nbr;
    std::vector<std::array<int, 2>> deme_coord;
    deme_coord.reserve(samp_param.Sample_coord_vec.size());
    std::array<int, 2> old_coord = {samp_param.Sample_coord_vec[0].at(0) - 1, 0};

    for (auto const &coord : samp_param.Sample_coord_vec)
    {
        if (old_coord != coord)
        {
            deme_coord.push_back(coord);
        }
        old_coord = coord;
    }
    deme_coord.shrink_to_fit();
    Nbr_of_deme = deme_coord.size();
    double max_dist = 0;

    Dist_btw_deme.resize(Nbr_of_deme, std::vector<double>(Nbr_of_deme, 0));
    for (auto deme1 = 0; deme1 < Nbr_of_deme; ++deme1)
    {
        for (auto deme2 = 0; deme2 < Nbr_of_deme; ++deme2)
        {
            //euclidien dist
            Dist_btw_deme[deme1][deme2] = sqrt(pow(deme_coord[deme1].at(0) - deme_coord[deme2].at(0), 2) + pow(deme_coord[deme1].at(1) - deme_coord[deme2].at(1), 2));
            if (Dist_btw_deme[deme1][deme2] > max_dist)
            {
                max_dist = Dist_btw_deme[deme1][deme2];
            }
        }
    }

    double dist_btw_class = max_dist / Dist_class_nbr;
    Dist_class_btw_deme = std::vector<std::vector<int>>(Nbr_of_deme, std::vector<int>(Nbr_of_deme));
    {
        for (std::size_t deme1 = 0; deme1 < Dist_btw_deme.size(); ++deme1)
        {
            for (std::size_t deme2 = 0; deme2 < Dist_btw_deme.size(); ++deme2)
            {
                //if Dist_btw_deme  =  0 class 0 else floor(Dist_btw_deme/dist_btw_class) to be [class_limit_min, class_limit_max[
                //
                Dist_class_btw_deme[deme1][deme2] = (static_cast<bool>(Dist_btw_deme[deme1][deme2]) ? ceil(Dist_btw_deme[deme1][deme2] / dist_btw_class) - 1 : 0);
                //Handle trouble with value who are close to max_dist and who can be ceil at class + 1
                if (Dist_class_btw_deme[deme1][deme2] == Dist_class_nbr)
                {
                    --Dist_class_btw_deme[deme1][deme2];
                }
            }
        }
    }

    //WARNING : Work if sample have uniform repartition
    //TODO : A changer
    Nbr_of_indiv_per_deme.clear();
    Nbr_of_indiv_per_deme.resize(Nbr_of_deme);

    auto itr = samp_param.Sample_size_per_node.begin();
    for (auto itr2 = Nbr_of_indiv_per_deme.begin();itr2 < Nbr_of_indiv_per_deme.end(); ++itr2)
    {
        *itr2 = *itr;
        if (samp_param.Sample_size_per_node.size() != 1 ) ++itr;
    }

    Locus_nbr = samp_param.Chr_nbr * samp_param.Sequence_length;
    Cumul_nbr_of_indiv_per_deme.reserve(Nbr_of_deme);

    for (auto nbr_indiv : Nbr_of_indiv_per_deme)
    {
        Cumul_nbr_of_indiv_per_deme.push_back(Nbr_of_indiv_tot);
        Nbr_of_indiv_tot += nbr_indiv;
    }

    set_indiv_feature();

    Plane_vec.reserve(Nbr_of_indiv_tot * Locus_nbr * Ploidy);
    Allele_state_per_loc.resize(Locus_nbr);

    //Sample with no missing data
    Nomiss_nbr_of_gene_per_loc.reserve(Locus_nbr);
    Nomiss_nbr_of_indiv_per_loc.reserve(Locus_nbr);
    Nomiss_nbr_of_deme_per_loc.reserve(Locus_nbr);
    //matrix(indiv, loc) with only 1
    Nomiss_indiv_bool_per_loc.resize(Nbr_of_indiv_tot, bin_vec(Locus_nbr));

    //iterator on each indiv to find the next mut site
    std::vector<std::vector<std::pair<int, int>>::const_iterator> sample_mutated_state_indiv(Nbr_of_indiv_tot * Ploidy);
    for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
    {
        auto sample_mutated_state_indiv_itr = sample_mutated_state_indiv.begin();
        for (auto const &indiv : sample_mutated_state_chr[chr])
        {
            *sample_mutated_state_indiv_itr = indiv.cbegin();
            ++sample_mutated_state_indiv_itr;
        }

        for (int locus = 0; locus < samp_param.Sequence_length; ++locus)
        {
            std::map<int, int> temp_count_allele_state;
            for (int indiv = 0; indiv < Nbr_of_indiv_tot * Ploidy; ++indiv)
            {
                int value;
                //if nbr_site == locus
                if ((sample_mutated_state_indiv.at(indiv) != sample_mutated_state_chr[chr][indiv].cend()) && (sample_mutated_state_indiv.at(indiv)->first == locus))
                {
                    value = sample_mutated_state_indiv.at(indiv)->second;
                    ++sample_mutated_state_indiv.at(indiv);
                }
                else
                {
                    value = ancestry_seq[chr][locus];
                }

                //emplace return a pair consisting of an iterator to the inserted element (or element already in place) and a bool denoting whether the insertion took place
                auto pair = temp_count_allele_state.emplace(value, 1);
                if (!pair.second)
                {
                    pair.first->second += 1;
                }

                Plane_vec.push_back(value);
            }

            Allele_state_per_loc[chr + locus].reserve(temp_count_allele_state.size());
            for (auto const &pair : temp_count_allele_state)
            {
                Allele_state_per_loc[chr + locus].push_back({pair.first, pair.second});
            }

            //Sample have no missing data
            Nomiss_nbr_of_gene_per_loc.push_back(Nbr_of_indiv_tot * Ploidy);
            Nomiss_nbr_of_indiv_per_loc.push_back(Nbr_of_indiv_tot);
            Nomiss_nbr_of_deme_per_loc.push_back(Nbr_of_deme);
        }
    }

    //Sample have no missing data
    Nomiss_nbr_of_gene_per_loc_per_deme.resize(Locus_nbr, std::vector<int>(Nbr_of_deme));
    Nomiss_nbr_of_indiv_per_loc_per_deme.resize(Locus_nbr, std::vector<int>(Nbr_of_deme));
    for (int locus = 0; locus < Locus_nbr; ++locus)
    {
        for (int deme = 0; deme < Nbr_of_deme; ++deme)
        {
            Nomiss_nbr_of_gene_per_loc_per_deme[locus][deme] = Nbr_of_indiv_per_deme[locus] * Ploidy;
            Nomiss_nbr_of_indiv_per_loc_per_deme[locus][deme] = Nbr_of_indiv_per_deme[locus];
        }
    }

    //Genetic map

    Dist_btw_locus.resize(Locus_nbr);
    for (auto locus1 = 0; locus1 < Locus_nbr; ++locus1)
    {
        Dist_btw_locus[locus1].resize(Locus_nbr);
        for (auto locus2 = 0; locus2 < Locus_nbr; ++locus2)
        {
            Dist_btw_locus[locus1][locus2] = (locus1 != locus2);
        }
    }
}

void data_plane_vec_c::update_data_plane_vec(std::vector<std::vector<std::vector<std::pair<int, int>>>> const &sample_mutated_state_chr, std::vector<std::vector<int>> const &ancestry_seq, samp_param_c const &samp_param)
{
    Plane_vec.clear();
    Allele_state_per_loc.clear();

    Plane_vec.reserve(Nbr_of_indiv_tot * Locus_nbr * Ploidy);
    Allele_state_per_loc.resize(Locus_nbr);

    //iterator on each indiv to find the next mut site
    auto sample_mutated_state_indiv = std::vector<std::vector<std::pair<int, int>>::const_iterator>(sample_mutated_state_chr.at(0).size());
    for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
    {
        auto sample_mutated_state_indiv_itr = sample_mutated_state_indiv.begin();
        for (auto const &indiv : sample_mutated_state_chr[chr])
        {
            *sample_mutated_state_indiv_itr = indiv.cbegin();
            ++sample_mutated_state_indiv_itr;
        }

        for (int locus = 0; locus < samp_param.Sequence_length; ++locus)
        {
            std::map<int, int> temp_count_allele_state;
            for (int indiv = 0; indiv < Nbr_of_indiv_tot * Ploidy; ++indiv)
            {
                int value;
                //if nbr_site == locus
                if ((sample_mutated_state_indiv.at(indiv) != sample_mutated_state_chr[chr][indiv].cend()) && (sample_mutated_state_indiv.at(indiv)->first == locus))
                {
                    value = sample_mutated_state_indiv.at(indiv)->second;
                    ++sample_mutated_state_indiv.at(indiv);
                }
                else
                {
                    value = ancestry_seq[chr][locus];
                }

                //emplace return a pair consisting of an iterator to the inserted element (or element already in place) and a bool denoting whether the insertion took place
                auto pair = temp_count_allele_state.emplace(value, 1);
                if (!pair.second)
                {
                    pair.first->second += 1;
                }

                Plane_vec.push_back(value);
            }

            Allele_state_per_loc[chr + locus].reserve(temp_count_allele_state.size());
            for (auto const &pair : temp_count_allele_state)
            {
                Allele_state_per_loc[chr + locus].push_back({pair.first, pair.second});
            }
        }
    }
}

void data_plane_vec_c::set_indiv_feature()
{
    //Each indiv have attribut
    Indiv_feat.resize(Nbr_of_indiv_tot);
    auto indiv_feat_itr = Indiv_feat.begin();
    for (int deme = 0; deme < Nbr_of_deme; ++deme)
    {
        for (int indiv = 0; indiv < Nbr_of_indiv_per_deme[deme]; ++indiv)
        {
            indiv_feat_itr->Deme = deme;
            ++indiv_feat_itr;
        }
    }
}

// data_plane_vec_c &data_plane_vec_c::operator=(data_plane_vec_c data_plane_vec)
// {
//     std::swap(Plane_vec, data_plane_vec.Plane_vec);
//     std::swap(Indiv_feat, data_plane_vec.Indiv_feat);
//     std::swap(Dist_btw_deme, data_plane_vec.Dist_btw_deme);
//     std::swap(Dist_class_btw_deme, data_plane_vec.Dist_class_btw_deme);
//     std::swap(Dist_class_nbr, data_plane_vec.Dist_class_nbr);

//     std::swap(Nbr_of_deme, data_plane_vec.Nbr_of_deme);

//     std::swap(Nbr_of_indiv_tot, data_plane_vec.Nbr_of_indiv_tot);
//     std::swap(Nbr_of_indiv_per_deme, data_plane_vec.Nbr_of_indiv_per_deme);
//     std::swap(Cumul_nbr_of_indiv_per_deme, data_plane_vec.Cumul_nbr_of_indiv_per_deme);
//     std::swap(Locus_nbr, data_plane_vec.Locus_nbr);

//     std::swap(Nomiss_nbr_of_gene_per_loc, data_plane_vec.Nomiss_nbr_of_gene_per_loc);
//     std::swap(Nomiss_nbr_of_indiv_per_loc, data_plane_vec.Nomiss_nbr_of_indiv_per_loc);
//     std::swap(Nomiss_nbr_of_deme_per_loc, data_plane_vec.Nomiss_nbr_of_deme_per_loc);

//     std::swap(Nomiss_indiv_bool_per_loc, data_plane_vec.Nomiss_indiv_bool_per_loc);

//     std::swap(Nomiss_nbr_of_gene_per_loc_per_deme, data_plane_vec.Nomiss_nbr_of_gene_per_loc_per_deme);
//     std::swap(Nomiss_nbr_of_indiv_per_loc_per_deme, data_plane_vec.Nomiss_nbr_of_indiv_per_loc_per_deme);

//     std::swap(Allele_state_per_loc, data_plane_vec.Allele_state_per_loc);
//     std::swap(Ploidy, data_plane_vec.Ploidy);

//     return *this;
// }

int data_plane_vec_c::get_Ploidy() const
{
    return Ploidy;
}

int data_plane_vec_c::size() const
{
    return Plane_vec.size();
}
int data_plane_vec_c::base_nbr_locus_per_indiv() const
{
    return Locus_nbr;
}
int data_plane_vec_c::nbr_of_deme() const
{
    return Nbr_of_deme;
}
int data_plane_vec_c::nbr_of_locus_tot() const
{
    return Nbr_of_indiv_tot * Ploidy * Locus_nbr;
}
int data_plane_vec_c::nbr_of_indiv() const
{
    return Nbr_of_indiv_tot;
}
int data_plane_vec_c::nbr_of_indiv_per_deme(int nbr_of_deme) const
{
    return Nbr_of_indiv_per_deme[nbr_of_deme];
}
std::vector<int> const &data_plane_vec_c::cumul_nbr_of_indiv_per_deme() const
{
    return Cumul_nbr_of_indiv_per_deme;
}

int data_plane_vec_c::get_indiv(int gene_index) const
{
    int place_in_locus = gene_index % (Nbr_of_indiv_tot * Ploidy);

    return place_in_locus / Ploidy;
}
feature_c const &data_plane_vec_c::get_feature(int indiv)
{
    return Indiv_feat[indiv];
}

int data_plane_vec_c::nomiss_nbr_of_gene_per_loc(int locus) const
{
    return Nomiss_nbr_of_gene_per_loc[locus];
}
int data_plane_vec_c::nomiss_nbr_of_indiv_per_loc(int locus) const
{
    return Nomiss_nbr_of_indiv_per_loc[locus];
}
std::vector<int> const &data_plane_vec_c::nomiss_nbr_of_gene_per_loc_per_deme(int locus) const
{
    return Nomiss_nbr_of_gene_per_loc_per_deme[locus];
}
int data_plane_vec_c::nomiss_nbr_of_gene_per_loc_per_deme(int locus, int deme) const
{
    return Nomiss_nbr_of_gene_per_loc_per_deme[locus][deme];
}
std::vector<int> const &data_plane_vec_c::nomiss_nbr_of_indiv_per_loc_per_deme(int locus) const
{
    return Nomiss_nbr_of_indiv_per_loc_per_deme[locus];
}
int data_plane_vec_c::nomiss_nbr_of_indiv_per_loc_per_deme(int locus, int deme) const
{
    return Nomiss_nbr_of_indiv_per_loc_per_deme[locus][deme];
}
int data_plane_vec_c::nomiss_nbr_of_deme_per_loc(int locus) const
{
    return Nomiss_nbr_of_deme_per_loc[locus];
}

int data_plane_vec_c::nbr_allele_per_loc(int locus) const
{
    return Allele_state_per_loc[locus].size();
}

std::vector<std::array<int, 2>> const &data_plane_vec_c::allele_state_per_loc(int locus) const
{
    return Allele_state_per_loc[locus];
}

std::vector<int> const &data_plane_vec_c::get_plane_vec()
{
    return Plane_vec;
}

int data_plane_vec_c::operator[](int i) const
{
    return Plane_vec[i];
}

std::vector<int>::const_iterator data_plane_vec_c::begin() const
{
    //cbegin => const_iter
    return Plane_vec.cbegin();
}
std::vector<int>::const_iterator data_plane_vec_c::end() const
{
    return Plane_vec.cend();
}

int const &data_plane_vec_c::operator()(int locus, int deme, int indiv, int gene) const
{
    if (gene > Ploidy - 1)
    {
        ++gene;
        throw std::logic_error("( Can't show the gene " + std::to_string(gene) + " when max gene by indiv is " + std::to_string(Ploidy) + ". Contact the developpers. I exit. )"); // TODO éclaircir ce message obscur ?
    }
    return Plane_vec[(Nbr_of_indiv_tot * locus + Cumul_nbr_of_indiv_per_deme[deme] + indiv) * Ploidy + gene]; //(Ploidy - 1) if operator use in haploid return the same for gene 0 and gene 1
}

int const &data_plane_vec_c::operator()(int locus, int indiv, int gene) const
{
    if (gene > Ploidy - 1)
    {
        ++gene;
        throw std::logic_error("( Can't show the gene " + std::to_string(gene) + " when max gene by indiv is " + std::to_string(Ploidy) + ". Contact the developpers. I exit. )");
    }

    if (indiv > Nbr_of_indiv_tot - 1)
    {
        throw std::logic_error("( Only " + std::to_string(Nbr_of_indiv_tot) + " individuals for locus " + std::to_string(locus) + ". Contact the developpers. I exit. )"); // TODO éclaircir ce message obscur ?
    }

    return Plane_vec[(Nbr_of_indiv_tot * locus + indiv) * Ploidy + gene]; //(Ploidy - 1) if operator use in haploid return the same for gene 0 and gene 1
}

int data_plane_vec_c::index_begin_locus(int locus) const
{
    return Nbr_of_indiv_tot * locus * Ploidy;
}

int data_plane_vec_c::index_end_locus(int locus) const
{
    return Nbr_of_indiv_tot * (locus + 1) * Ploidy;
}
//In same locus
bool data_plane_vec_c::same_loc_in_indiv(int dpv_gene_index1, int dpv_gene_index2) const
{
    if (dpv_gene_index1 == dpv_gene_index2)
    {
        return true;
    }
    if (dpv_gene_index1 > dpv_gene_index2)
    {
        auto temp = dpv_gene_index1;
        dpv_gene_index1 = dpv_gene_index2;
        dpv_gene_index2 = temp;
    }

    bool result{false};
    if (Ploidy == 2)
    {
        result = (dpv_gene_index1 % 2 == 0) && (dpv_gene_index2 - dpv_gene_index1 == 1);
    }
    return result;
}

bool data_plane_vec_c::same_deme(int dpv_gene_index1, int dpv_gene_index2) const
{
    if (dpv_gene_index1 == dpv_gene_index2)
    {
        return true;
    }

    return Indiv_feat[get_indiv(dpv_gene_index1)].Deme == Indiv_feat[get_indiv(dpv_gene_index2)].Deme;
}

bin_vec const &data_plane_vec_c::nomiss_data_indiv(int indiv) const
{
    return Nomiss_indiv_bool_per_loc[indiv];
}

bool data_plane_vec_c::nomiss_data_indiv_per_loc(int indiv, int locus) const
{
    return Nomiss_indiv_bool_per_loc[indiv].at(locus);
}

//Passer par un tableau d'attribut des indivs
double data_plane_vec_c::dist_btw_deme(int dpv_gene_index1, int dpv_gene_index2) const
{
    if (dpv_gene_index1 == dpv_gene_index2)
    {
        return 0;
    }
    auto const deme_gen1 = Indiv_feat[get_indiv(dpv_gene_index1)].Deme;
    auto const deme_gen2 = Indiv_feat[get_indiv(dpv_gene_index2)].Deme;

    return Dist_btw_deme[deme_gen1][deme_gen2];
}

double data_plane_vec_c::dist_btw_deme_with_deme(int deme_index1, int deme_index2) const
{
    if (deme_index1 == deme_index2)
    {
        return 0;
    }

    return Dist_btw_deme[deme_index1][deme_index2];
}

int data_plane_vec_c::nbr_of_dist_class() const
{
    return Dist_class_nbr;
}

int data_plane_vec_c::dist_class_btw_deme(int dpv_gene_index1, int dpv_gene_index2) const
{
    if (dpv_gene_index1 == dpv_gene_index2)
    {
        return 0;
    }
    auto const deme_gen1 = Indiv_feat[get_indiv(dpv_gene_index1)].Deme;
    auto const deme_gen2 = Indiv_feat[get_indiv(dpv_gene_index2)].Deme;

    return Dist_class_btw_deme[deme_gen1][deme_gen2];
}

double data_plane_vec_c::dist_btw_locus(int locus_index1, int locus_index2) const
{
    if (locus_index1 == locus_index2)
    {
        return 0;
    }

    return Dist_btw_locus[locus_index1][locus_index2];
}
