#include "calc_stat_dl.hpp"
#include "calc_stat.hpp"

//<locus<chr>>
double calc_phi_ij_xy(std::array<std::array<int, 4>, 2> const &locus_chr)
{
    double minor_phi = (((locus_chr[0][0] == locus_chr[0][2]) && (locus_chr[1][0] == locus_chr[1][2])) + ((locus_chr[0][0] == locus_chr[0][3]) && (locus_chr[1][0] == locus_chr[1][3])) + ((locus_chr[0][1] == locus_chr[0][2]) && (locus_chr[1][1] == locus_chr[1][2])) + ((locus_chr[0][1] == locus_chr[0][3]) && (locus_chr[1][1] == locus_chr[1][3]))) / 16.0;
    double gamma1 = (((locus_chr[0][0] == locus_chr[0][2]) && (locus_chr[1][0] == locus_chr[1][3])) + ((locus_chr[0][0] == locus_chr[0][3]) && (locus_chr[1][0] == locus_chr[1][2])) + ((locus_chr[0][1] == locus_chr[0][2]) && (locus_chr[1][1] == locus_chr[1][3])) + ((locus_chr[0][1] == locus_chr[0][3]) && (locus_chr[1][1] == locus_chr[1][2]))) / 16.0;
    double gamma2 = (((locus_chr[0][0] == locus_chr[0][2]) && (locus_chr[1][1] == locus_chr[1][2])) + ((locus_chr[0][1] == locus_chr[0][2]) && (locus_chr[1][0] == locus_chr[1][2])) + ((locus_chr[0][0] == locus_chr[0][3]) && (locus_chr[1][1] == locus_chr[1][3])) + ((locus_chr[0][1] == locus_chr[0][3]) && (locus_chr[1][0] == locus_chr[1][3]))) / 16.0;
    double delta = (((locus_chr[0][0] == locus_chr[0][2]) && (locus_chr[1][1] == locus_chr[1][3])) + ((locus_chr[0][0] == locus_chr[0][3]) && (locus_chr[1][1] == locus_chr[1][2])) + ((locus_chr[0][1] == locus_chr[0][2]) && (locus_chr[1][0] == locus_chr[1][3])) + ((locus_chr[0][1] == locus_chr[0][3]) && (locus_chr[1][0] == locus_chr[1][2]))) / 16.0;
    return (minor_phi + gamma1 + gamma2 + delta);
}
//Phi for all indiv, calc without missing data
std::vector<double> calc_phi_ij(data_plane_vec_c const &data_plane_vec, int ploidy)
{
    int locus_pair_nbr = (data_plane_vec.base_nbr_locus_per_indiv() * (data_plane_vec.base_nbr_locus_per_indiv() - 1)) / 2;

    std::vector<double> result;
    result.reserve(locus_pair_nbr);

    for (int locus_i = 0; locus_i < data_plane_vec.base_nbr_locus_per_indiv(); ++locus_i)
    {
        for (int locus_j = locus_i + 1; locus_j < data_plane_vec.base_nbr_locus_per_indiv(); ++locus_j)
        {
            int indiv_pair_nbr = 0;
            double value = 0;
            if (ploidy == 2)
            {
                for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv(); ++indiv)
                {
                    int locus_i_indiv1_gene1 = data_plane_vec(locus_i, indiv, 0);
                    int locus_i_indiv1_gene2 = data_plane_vec(locus_i, indiv, 1);

                    int locus_j_indiv1_gene1 = data_plane_vec(locus_j, indiv, 0);
                    int locus_j_indiv1_gene2 = data_plane_vec(locus_j, indiv, 1);

                    for (int other_indiv = indiv + 1; other_indiv < data_plane_vec.nbr_of_indiv(); ++other_indiv)
                    {
                        int locus_i_indiv2_gene1 = data_plane_vec(locus_i, other_indiv, 0);
                        int locus_i_indiv2_gene2 = data_plane_vec(locus_i, other_indiv, 1);

                        int locus_j_indiv2_gene1 = data_plane_vec(locus_j, other_indiv, 0);
                        int locus_j_indiv2_gene2 = data_plane_vec(locus_j, other_indiv, 1);

                        value += calc_phi_ij_xy({{{locus_i_indiv1_gene1, locus_i_indiv1_gene2, locus_i_indiv2_gene1, locus_i_indiv2_gene2},
                                                  {locus_j_indiv1_gene1, locus_j_indiv1_gene2, locus_j_indiv2_gene1, locus_j_indiv2_gene2}}});
                        ++indiv_pair_nbr;
                    }
                }
            }
            else
            {
                for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv(); ++indiv)
                {
                    int locus_i_indiv1 = data_plane_vec(locus_i, indiv, 0);
                    int locus_j_indiv1 = data_plane_vec(locus_j, indiv, 0);

                    for (int other_indiv = indiv + 1; other_indiv < data_plane_vec.nbr_of_indiv(); ++other_indiv)
                    {
                        int locus_i_indiv2 = data_plane_vec(locus_i, other_indiv, 0);
                        int locus_j_indiv2 = data_plane_vec(locus_j, other_indiv, 0);

                        value += (locus_i_indiv1 == locus_i_indiv2) && (locus_j_indiv1 == locus_j_indiv2);
                        ++indiv_pair_nbr;
                    }
                }
            }
            result.push_back(value / indiv_pair_nbr);
        }
    }
    return result;
}

//<locus_pair_nbr<value eta>>
double calc_eta_ij_xy(data_plane_vec_c const &data_plane_vec, int locus_i, double Q2_loc_i, double Q1_loc_i, int locus_j, double Q2_loc_j, double Q1_loc_j, int deme_x, int deme_y)
{
    double phi = 0;
    double div = 0;
    if (data_plane_vec.get_Ploidy() == 2)
    {
        for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv_per_deme(deme_x); ++indiv)
        {
            int locus_i_indiv1_gene1 = data_plane_vec(locus_i, deme_x, indiv, 0);
            int locus_i_indiv1_gene2 = data_plane_vec(locus_i, deme_x, indiv, 1);

            int locus_j_indiv1_gene1 = data_plane_vec(locus_j, deme_x, indiv, 0);
            int locus_j_indiv1_gene2 = data_plane_vec(locus_j, deme_x, indiv, 1);
            //In AAaaBBbb if A&B missing => phi = 0, etc
            if (((locus_i_indiv1_gene1 != 0) || (locus_i_indiv1_gene2 != 0)) && ((locus_j_indiv1_gene1 != 0) || (locus_j_indiv1_gene2 != 0)))
            {
                for (int indiv_other_deme = 0; indiv_other_deme < data_plane_vec.nbr_of_indiv_per_deme(deme_y); ++indiv_other_deme)
                {
                    int locus_i_indiv2_gene1 = data_plane_vec(locus_i, deme_y, indiv_other_deme, 0);
                    int locus_i_indiv2_gene2 = data_plane_vec(locus_i, deme_y, indiv_other_deme, 1);

                    int locus_j_indiv2_gene1 = data_plane_vec(locus_j, deme_y, indiv_other_deme, 0);
                    int locus_j_indiv2_gene2 = data_plane_vec(locus_j, deme_y, indiv_other_deme, 1);
                    if (((locus_i_indiv2_gene1 != 0) || (locus_i_indiv2_gene2 != 0)) && ((locus_j_indiv2_gene1 != 0) || (locus_j_indiv2_gene2 != 0)))
                    {
                        int miss = 0;
                        double pond;
                        //handle missing value
                        std::array<std::array<int, 4>, 2> locus_chr = {
                            {{locus_i_indiv1_gene1, locus_i_indiv1_gene2, locus_i_indiv2_gene1, locus_i_indiv2_gene2},
                             {locus_j_indiv1_gene1, locus_j_indiv1_gene2, locus_j_indiv2_gene1, locus_j_indiv2_gene2}}};
                        for (auto &value_loc : locus_chr)
                        {
                            for (auto &value : value_loc)
                            { //To avoid 0==0 test in phi calculation
                                if (value == 0)
                                {
                                    value = --miss;
                                }
                            }
                        }
                        //pond depend of missing value number
                        switch (miss)
                        {
                        case 0:
                        {
                            pond = 1;
                            break;
                        }
                        case -1:
                        {
                            pond = 0.5;
                            break;
                        }
                        case -2:
                        {
                            pond = 0.25;
                            break;
                        }
                        case -3:
                        {
                            pond = 0.125;
                            break;
                        }
                        default:
                        {
                            pond = 0.0625;
                            break;
                        }
                        }

                        phi += pond * calc_phi_ij_xy(locus_chr);
                        div += pond;
                    }
                }
            }
        }
    }
    else
    {
        for (int indiv = 0; indiv < data_plane_vec.nbr_of_indiv_per_deme(deme_x); ++indiv)
        {
            int locus_i_indiv1 = data_plane_vec(locus_i, deme_x, indiv, 0);
            int locus_j_indiv1 = data_plane_vec(locus_j, deme_x, indiv, 0);
            if ((locus_i_indiv1 != 0) && (locus_j_indiv1 != 0))
            {
                for (int indiv_other_deme = 0; indiv_other_deme < data_plane_vec.nbr_of_indiv_per_deme(deme_y); ++indiv_other_deme)
                {
                    int locus_i_indiv2 = data_plane_vec(locus_i, deme_y, indiv_other_deme, 0);
                    int locus_j_indiv2 = data_plane_vec(locus_j, deme_y, indiv_other_deme, 0);
                    if ((locus_i_indiv2 != 0) && (locus_j_indiv2 != 0))
                    {
                        phi += (locus_i_indiv1 == locus_i_indiv2) && (locus_j_indiv1 == locus_j_indiv2);
                        ++div;
                    }
                }
            }
        }
    }

    return (phi / div - Q1_loc_i * Q1_loc_j) / ((1 - Q2_loc_i) * (1 - Q2_loc_j));
}
//<deme_pair_nbr,<dist-deme, dist-locus, value eta>>
std::vector<std::array<double, 3>> calc_eta_ij(data_plane_vec_c const &data_plane_vec, int locus_i, double Q2_loc_i, int locus_j, double Q2_loc_j)
{
    int deme_pair_nbr = (data_plane_vec.nbr_of_deme() * (data_plane_vec.nbr_of_deme() - 1)) / 2;
    std::vector<std::array<double, 3>> result;
    result.reserve(deme_pair_nbr);

    std::vector<std::array<int, 2>> Q1_locus_i_deme;
    Q1_locus_i_deme.reserve(data_plane_vec.nbr_of_deme());

    std::vector<std::array<int, 2>> Q1_locus_j_deme;
    Q1_locus_j_deme.reserve(data_plane_vec.nbr_of_deme());

    for (int deme = 0; deme < data_plane_vec.nbr_of_deme(); ++deme)
    {
        Q1_locus_i_deme.push_back(calc_Q_inter_indiv_per_locus_per_deme(data_plane_vec, locus_i, deme));

        Q1_locus_j_deme.push_back(calc_Q_inter_indiv_per_locus_per_deme(data_plane_vec, locus_j, deme));
    }

    double dist_locus = data_plane_vec.dist_btw_locus(locus_i, locus_j);
    for (int deme_x = 0; deme_x < data_plane_vec.nbr_of_deme(); ++deme_x)
    {
        for (int deme_y = deme_x + 1; deme_y < data_plane_vec.nbr_of_deme(); ++deme_y)
        {
            auto Q1_loc_i_xy = static_cast<double>(Q1_locus_i_deme[deme_x].at(0) + Q1_locus_i_deme[deme_y].at(0)) / (Q1_locus_i_deme[deme_x].at(1) + Q1_locus_i_deme[deme_y].at(1));
            auto Q1_loc_j_xy = static_cast<double>(Q1_locus_j_deme[deme_x].at(0) + Q1_locus_j_deme[deme_y].at(0)) / (Q1_locus_j_deme[deme_x].at(1) + Q1_locus_j_deme[deme_y].at(1));
            auto eta = calc_eta_ij_xy(data_plane_vec, locus_i, Q2_loc_i, Q1_loc_i_xy, locus_j, Q2_loc_j, Q1_loc_j_xy, deme_x, deme_y);
            //dist-deme, dist-locus, value eta
            result.push_back({{data_plane_vec.dist_btw_deme_with_deme(deme_x, deme_y), dist_locus, eta}});
        }
    }

    return result;
}

//<pair of deme * pair of locus,<dist-locus, dist-deme, value eta>>
std::vector<std::array<double, 3>> calc_eta(data_plane_vec_c const &data_plane_vec)
{
    int locus_pair_nbr = (data_plane_vec.base_nbr_locus_per_indiv() * (data_plane_vec.base_nbr_locus_per_indiv() - 1)) / 2;
    int deme_pair_nbr = (data_plane_vec.nbr_of_deme() * (data_plane_vec.nbr_of_deme() - 1)) / 2;

    std::vector<std::array<double, 3>> result;
    result.reserve(locus_pair_nbr * deme_pair_nbr);

    std::vector<double> Q2_locus;
    Q2_locus.reserve(data_plane_vec.base_nbr_locus_per_indiv());
    for (int locus = 0; locus < data_plane_vec.base_nbr_locus_per_indiv(); ++locus)
    {
        auto temp = calc_Q_inter_deme_per_locus(data_plane_vec, locus);
        Q2_locus.push_back(static_cast<double>(temp.at(0)) / temp.at(1));
    }

    for (int locus_i = 0; locus_i < data_plane_vec.base_nbr_locus_per_indiv(); ++locus_i)
    {
        for (int locus_j = locus_i + 1; locus_j < data_plane_vec.base_nbr_locus_per_indiv(); ++locus_j)
        {
            auto temp_vec = calc_eta_ij(data_plane_vec, locus_i, Q2_locus[locus_i], locus_j, Q2_locus[locus_j]);
            result.insert(result.end(), temp_vec.begin(), temp_vec.end());
        }
    }

    return result;
}

//Version diploide and +1 indiv/deme
std::vector<std::array<double, 3>> calc_eta_q1_version(data_plane_vec_c const &data_plane_vec)
{
    int locus_pair_nbr = (data_plane_vec.base_nbr_locus_per_indiv() * (data_plane_vec.base_nbr_locus_per_indiv() - 1)) / 2;
    int deme_pair_nbr = (data_plane_vec.nbr_of_deme() * (data_plane_vec.nbr_of_deme() - 1)) / 2;

    std::vector<std::array<double, 3>> result;
    result.reserve(locus_pair_nbr * deme_pair_nbr);

    for (int locus_i = 0; locus_i < data_plane_vec.base_nbr_locus_per_indiv(); ++locus_i)
    {
        for (int locus_j = locus_i + 1; locus_j < data_plane_vec.base_nbr_locus_per_indiv(); ++locus_j)
        {
            std::vector<std::array<int, 2>> Q1_locus_i_deme;
            Q1_locus_i_deme.reserve(data_plane_vec.nbr_of_deme());

            std::vector<std::array<int, 2>> Q1_locus_j_deme;
            Q1_locus_j_deme.reserve(data_plane_vec.nbr_of_deme());

            for (int deme = 0; deme < data_plane_vec.nbr_of_deme(); ++deme)
            {
                Q1_locus_i_deme.push_back(calc_Q_inter_indiv_per_locus_per_deme(data_plane_vec, locus_i, deme));

                Q1_locus_j_deme.push_back(calc_Q_inter_indiv_per_locus_per_deme(data_plane_vec, locus_j, deme));
            }

            double dist_locus = data_plane_vec.dist_btw_locus(locus_i, locus_j);

            for (int deme_x = 0; deme_x < data_plane_vec.nbr_of_deme(); ++deme_x)
            {
                for (int deme_y = deme_x + 1; deme_y < data_plane_vec.nbr_of_deme(); ++deme_y)
                {
                    auto Q1_loc_i_xy = static_cast<double>(Q1_locus_i_deme[deme_x].at(0) + Q1_locus_i_deme[deme_y].at(0)) / (Q1_locus_i_deme[deme_x].at(1) + Q1_locus_i_deme[deme_y].at(1));
                    auto Q1_loc_j_xy = static_cast<double>(Q1_locus_j_deme[deme_x].at(0) + Q1_locus_j_deme[deme_y].at(0)) / (Q1_locus_j_deme[deme_x].at(1) + Q1_locus_j_deme[deme_y].at(1));
                    auto eta = calc_eta_ij_xy(data_plane_vec, locus_i, Q1_loc_i_xy, Q1_loc_i_xy, locus_j, Q1_loc_j_xy, Q1_loc_j_xy, deme_x, deme_y);
                    //dist-deme, dist-locus, value eta
                    result.push_back({{data_plane_vec.dist_btw_deme_with_deme(deme_x, deme_y), dist_locus, eta}});
                }
            }
        }
    }
    return result;
}
//Continous habitat isolation by distance
std::vector<std::array<double, 3>> calc_eta_1_indiv_deme_v(data_plane_vec_c const &data_plane_vec)
{
    int locus_pair_nbr = (data_plane_vec.base_nbr_locus_per_indiv() * (data_plane_vec.base_nbr_locus_per_indiv() - 1)) / 2;
    int deme_pair_nbr = (data_plane_vec.nbr_of_deme() * (data_plane_vec.nbr_of_deme() - 1)) / 2;

    std::vector<std::array<double, 3>> result;
    result.reserve(locus_pair_nbr * deme_pair_nbr);

    std::vector<double> Q2_locus;
    Q2_locus.reserve(data_plane_vec.base_nbr_locus_per_indiv());
    std::vector<double> Q0_locus;
    Q0_locus.reserve(data_plane_vec.base_nbr_locus_per_indiv());
    for (int locus = 0; locus < data_plane_vec.base_nbr_locus_per_indiv(); ++locus)
    {
        auto temp = calc_Q_inter_deme_per_locus(data_plane_vec, locus);
        Q2_locus.push_back(static_cast<double>(temp.at(0)) / temp.at(1));
        if (data_plane_vec.get_Ploidy() == 2)
        {
            temp = calc_Q_intra_indiv_per_locus(data_plane_vec, locus);
            Q0_locus.push_back(static_cast<double>(temp.at(0)) / temp.at(1));
        }
    }

    for (int locus_i = 0; locus_i < data_plane_vec.base_nbr_locus_per_indiv(); ++locus_i)
    {
        for (int locus_j = locus_i + 1; locus_j < data_plane_vec.base_nbr_locus_per_indiv(); ++locus_j)
        {
            double dist_locus = data_plane_vec.dist_btw_locus(locus_i, locus_j);

            for (int deme_x = 0; deme_x < data_plane_vec.nbr_of_deme(); ++deme_x)
            {
                for (int deme_y = deme_x + 1; deme_y < data_plane_vec.nbr_of_deme(); ++deme_y)
                {
                    double eta;
                    if (data_plane_vec.get_Ploidy() == 2)
                    {
                        eta = calc_eta_ij_xy(data_plane_vec, locus_i, Q2_locus[locus_i], Q0_locus[locus_i], locus_j, Q2_locus[locus_j], Q0_locus[locus_j], deme_x, deme_y);
                    }
                    if (data_plane_vec.get_Ploidy() == 1)
                    {
                        eta = calc_eta_ij_xy(data_plane_vec, locus_i, Q2_locus[locus_i], Q2_locus[locus_i], locus_j, Q2_locus[locus_j], Q2_locus[locus_j], deme_x, deme_y);
                    }
                    //dist-deme, dist-locus, value eta
                    result.push_back({{data_plane_vec.dist_btw_deme_with_deme(deme_x, deme_y), dist_locus, eta}});
                }
            }
        }
    }

    return result;
}