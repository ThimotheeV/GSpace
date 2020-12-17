#include "migration.hpp"
#include "common_tools.hpp"
#include "charfunc.hpp"
#include "summary_stat.hpp"

/*************************\
|       STATISTIC         |
|       DISPERSAL         |
|      DISTRIBUTION       |
\*************************/

/********fwd_disp_distrib_c********/

std::function<double(int)> fwd_disp_distrib_c::uniform()
{
    return [](int i) {
        return 1;
    };
}

std::function<double(int)> fwd_disp_distrib_c::geometric(double g_geo_param)
{
    return [g_geo_param](int i) {
        return g_geo_param * (std::pow((1 - g_geo_param), i - 1));
    };
}

//Not a pareto => Zeta distribution (discret equivalent of pareto)
std::function<double(int)> fwd_disp_distrib_c::pareto(double n_pareto_param)
{
    //calc value for Riemann zeta function
    //calc for a precision at 1e-8
    double zeta_sum = 0;
    double value = 1;
    for (int denom = 2; value > 1e-8; ++denom)
    {
        zeta_sum += value;
        value = 1.0 / std::pow(denom, n_pareto_param);
    }
    //Not a pareto => Zeta distribution (discret equivalent of pareto)
    return [n_pareto_param, zeta_sum](int i) {
        return (1 / std::pow(i, n_pareto_param)) / zeta_sum;
    };
}

//base on annexe of Rousset 1997
//no d because sum = d*sum(proba) so d*proba/sum = proba/sum(proba)
std::function<double(int)> fwd_disp_distrib_c::gaussian(int dist_max)
{
    int nbr_pos_case = 2 * dist_max;
    return [nbr_pos_case](int distance) {
        return combination(distance + nbr_pos_case / 2, nbr_pos_case) * pow(2, -nbr_pos_case);
    };
}

std::function<double(int)> fwd_disp_distrib_c::gamma_or_sichel(demo_param_c const &demo_p, int dim)
{
    //WARNING : sichel_pars[3] is a proba of axial uniform dispersal, 0 now but can be change later
    std::array<double, 4> sichel_pars{demo_p.Gamma_sichel_param, demo_p.Xi_sichel_param, demo_p.Omega_sichel_param, 0};
    int lpson = 0;
    double pi = 3.141592;
    long double d_max = demo_p.Disp_dist_max.at(dim);

    //sichel_pars[0]: Sichel gamma
    //sichel_pars[1]: Sichel xi
    //sichel_pars[2]: Sichel omega
    //sichel_pars[2]<0 is used as indicator that the distribution is the reciprocal Gamma fn.
    //sichel_pars[3] non charfunc fraction
    //sichel_pars[3] is a proba of axial uniform dispersal

    //what follows somehow assumes dim_reseau1=dim_reseau2 ?
    if (demo_p.Edge_effects == edge_effect_enum::circular)
    {
        lpson = demo_p.Lattice_size.at(dim);
    } /*ca donne des transformÈes sur un habitat fermÈ
                                                  => Sicheltable[->lpson]->Sicheltable[1]... pas terrible quand on veut calculer var et kurtosis;
                                                  ou bien il faut que lpson>2 dxmax...*/
    else
    {
        lpson = 100 * d_max;
    } // pour Èviter effet repliement
    long double DeuxPiSurNx = 2 * pi / lpson;

    std::vector<long double> costable(lpson);
    std::vector<double> migra(d_max + 1);
    //sichel_table=ldvector(lpson/2+1);
    std::vector<long double> sichel_table(lpson / 2 + 1);
    double temp = 0.0;
    double norm = 0.0;

    for (int i = 0; i < lpson; i++)
    {
        costable[i] = static_cast<long double>(cos(DeuxPiSurNx * i));
    }

    for (int i = 0; i < lpson / 2 + 1; i++)
    { //la valeur en lpson/2+1 est utilisÈe
        sichel_table[i] = Sichelcharfunc((long double)(2 * pi * i / lpson), sichel_pars);
    }

    for (int dx = 1; dx <= d_max; dx++)
    {
        temp = 1.0; // c'est charfn[0]
        for (int ii = 1; ii < lpson / 2; ii++)
        {
            //2 parce que les 2 termes sym de la transformÈe. mais c'est la valeur en |dx|
            temp += 2 * sichel_table[ii] * costable[(dx * ii) % lpson];
            //if ((dx==3500) && (ii%500==0)) {cout<<ii<<" "<<temp<<endl;getchar();}
        }
        if (lpson % 2 == 0)
        {
            temp += sichel_table[lpson / 2] * costable[(dx % 2) * lpson / 2];
        }
        temp /= lpson;
        migra[dx] = temp; // migration selon Sichel
    }
    // sichel_pars[3] is a proba of axial uniform dispersal, not Mig !
    for (int dist = 1; dist <= d_max; dist++)
    {
        temp = migra[dist] * (1. - sichel_pars[3]);
        if (dist > 0)
        {
            temp += sichel_pars[3] / lpson;
        }
        norm += ((dist == 0 ? 1. : 2.) * temp);
        migra[dist] = temp;
    } // for dist ...d_max
    for (int dist = 1; dist <= d_max; dist++)
    {
        migra[dist] /= norm;
    }

    return [migra](int dx) {
        return migra[dx] / 2;
    };
}

fwd_disp_distrib_c::fwd_disp_distrib_c(demo_param_c const &demo_param)
{

    Fwd_distrib[0] = construct_fwd_disp_distrib(demo_param.Disp_dist_max[0], demo_param.Proba_migr, demo_param.Disp_func.at(0));

    auto &info_collect = singleton_c<info_collector_c>::instance();

    if (info_collect.Check_disp_distrib)
    {
        if (demo_param.Dispersal_distrib == dispersal_distrib_enum::uniform)
        {
            info_collect.Fwd_axial_disp_theo_sig = INFINITY;
        }
        if (demo_param.Dispersal_distrib == dispersal_distrib_enum::geometric)
        {
            //        if(!Simple1DProductDispBool) {//RL072014 wascmp_nocase(twoD_disp,"1DProductWithoutm0")==0) {
            //            long double condsig2;
            //            if (dimRes2/vide>1) { //2D
            //              /*cf code R
            //                condaxialS2fromg<-function(gv) {
            //                    return((1+gv)/((2-gv)*(1-gv)^2))
            //                }\n\                      */
            //                condsig2=(1.+GeoG)/((2-GeoG)*(1.-GeoG)*(1.-GeoG));
            //            } else {
            //                condsig2=(1.+GeoG)/((1.-GeoG)*(1.-GeoG));
            //            }
            //            Nhm=ploidy*TVpars[0].initialDens*Mig;
            //            Dhs2=Nhm*condsig2;
            //            sig2=condsig2*Mig;
            //        } else {
            info_collect.Fwd_axial_disp_theo_sig = demo_param.Proba_migr * (1. + demo_param.G_geo_param) / ((1. - demo_param.G_geo_param) * (1. - demo_param.G_geo_param));
        }
        if (demo_param.Dispersal_distrib == dispersal_distrib_enum::pareto)
        {
            info_collect.Fwd_axial_disp_theo_sig = NAN; // TODO il doit bien exister une formule.
        }
        //    if (model_mig=='S') {
        ////cout<<"non charfunc fraction: "<<SichelDispPars[3]<<endl;
        //        if (SichelDispPars[2]<0) {
        //            sig2=-(1.-SichelDispPars[3])*SichelDispPars[1]/(2.*(1.+SichelDispPars[0]));
        //        } else {
        //            sig2=-(1.-SichelDispPars[3])*SichelDispPars[1]*bessel_k(1.+SichelDispPars[0],SichelDispPars[2],1)
        //                                            /(2.*bessel_k(SichelDispPars[0],SichelDispPars[2],1));
        //        }
    }
    if (demo_param.Disp_dist_max[1] == demo_param.Disp_dist_max[0])
    {
        Fwd_distrib[1] = Fwd_distrib[0];
    }
    else
    {
        Fwd_distrib[1] = construct_fwd_disp_distrib(demo_param.Disp_dist_max[1], demo_param.Proba_migr, demo_param.Disp_func.at(1));
    }
}

std::shared_ptr<std::vector<double>> fwd_disp_distrib_c::construct_fwd_disp_distrib(int dist_max, double proba_mig, std::function<double(int)> const &calc_proba)
{
    auto fwd_distrib = std::make_shared<std::vector<double>>();
    fwd_distrib->reserve(dist_max + 1);

    double cumul_sum = 0;
    //Create non adjust forward distriution

    fwd_distrib->push_back(1 - proba_mig);

    for (int i = 1; i <= dist_max; ++i)
    {
        double proba = calc_proba(i);
        cumul_sum += proba;
        fwd_distrib->push_back(proba);
    }
    cumul_sum /= proba_mig; //Fwd_distrib[i] with i € [1,n] be between 0 and 1
    //Fwd_distrib[i] *= proba_mig be between 0 and proba_mig

    //Adjust fw distri
    //one direction => create a semi-distri
    for (int i = 1; i <= dist_max; ++i)
    {
        fwd_distrib->at(i) /= cumul_sum; // Fwd_distrib[i] /= cumul_sum , RL I removed Fwd_distrib[i] /= 2
    }

    auto &info_collect = singleton_c<info_collector_c>::instance();
    if (info_collect.Check_disp_distrib)
    {
        if (info_collect.Fwd_axial_disp_distrib.empty())
        {
            info_collect.Fwd_axial_disp_distrib.reserve(dist_max + 1);
            info_collect.Fwd_axial_disp_mean = 0.0;
            info_collect.Fwd_axial_disp_sig = 0.0;
            info_collect.Fwd_axial_disp_kurt = 0.0;

            for (int i = 0; i <= dist_max; ++i)
            {
                info_collect.Fwd_axial_disp_distrib.push_back(fwd_distrib->at(i));
                info_collect.Fwd_axial_disp_mean += info_collect.Fwd_axial_disp_distrib.at(i) * i;
                info_collect.Fwd_axial_disp_sig += info_collect.Fwd_axial_disp_distrib.at(i) * i * i;
                info_collect.Fwd_axial_disp_kurt += info_collect.Fwd_axial_disp_distrib.at(i) * i * i * i * i;
            }
            info_collect.Fwd_axial_disp_kurt = info_collect.Fwd_axial_disp_kurt / (info_collect.Fwd_axial_disp_sig * info_collect.Fwd_axial_disp_sig) - 3.0;
        }
    }

    return fwd_distrib;
}

/********cumul_fwd_disp_distrib_c********/

cumul_fwd_disp_distrib_c::cumul_fwd_disp_distrib_c(rand_gen_c *rand_gen, fwd_disp_distrib_c const &fwd_distrib)
{
    Rand_gen = rand_gen; // for the draw() function
    // assign pointer for dim 0
    Cumul_fwd_distrib[0] = construct_cumul_fwd_disp_distrib(fwd_distrib[0]);

    // assign pointer for dim1
    if (fwd_distrib[0] == fwd_distrib[1])
    {
        Cumul_fwd_distrib[1] = Cumul_fwd_distrib[0];
    }
    else
    {
        Cumul_fwd_distrib[1] = construct_cumul_fwd_disp_distrib(fwd_distrib[1]);
    }

    auto &info_collect = singleton_c<info_collector_c>::instance();
    if (info_collect.Check_disp_distrib)
    {
        if (fwd_distrib[0] == fwd_distrib[1])
        {
            info_collect.Cumul_fwd_disp_distrib.resize(1);
        }
        else
        {
            info_collect.Cumul_fwd_disp_distrib.resize(2);
        }

        for (std::size_t i = 0; i < info_collect.Cumul_fwd_disp_distrib.size(); ++i)
        {
            info_collect.Cumul_fwd_disp_distrib.at(i).resize(Cumul_fwd_distrib.at(i)->size());
            for (std::size_t j = 0; j < info_collect.Cumul_fwd_disp_distrib.at(i).size(); ++j)
            {
                info_collect.Cumul_fwd_disp_distrib.at(i).at(std::get<0>(Cumul_fwd_distrib.at(i)->at(j)) + fwd_distrib[i].size() - 1) = 1. * std::get<1>(Cumul_fwd_distrib.at(i)->at(j)) / static_cast<double>(PRECISION);
            }
        }
    }
}

// this is the function to modify to implement asymetric migration (either two fwd distrib or a factor of asymetry)
std::shared_ptr<std::vector<std::tuple<int, int>>> cumul_fwd_disp_distrib_c::construct_cumul_fwd_disp_distrib(std::vector<double> const &fwd_distrib)
{
    int dist_max = fwd_distrib.size() - 1;
    //Bi_directionnal distri
    auto cumul_fwd_distrib = std::make_shared<std::vector<std::tuple<int, int>>>();
    cumul_fwd_distrib->reserve(2 * dist_max + 1);

    int cumul_proba = 0;

    //Prepare forward cumul distri with fw adjust
    // asymetry should be implemented here
    for (int i = 0; i <= 2 * dist_max; ++i)
    {
        //Symetric distri : dist_max ... 1 0 1 ... dist_max
        // construct the tuple directly in the vector
        //WARNING *10-6 for translate in int => optimization purpose, more light to drawn in int and compare int vs int
        if (i != dist_max)
        {
            cumul_fwd_distrib->emplace_back(i - dist_max, (fwd_distrib[std::abs(i - dist_max)] / 2) * PRECISION);
        }
        else
        {
            cumul_fwd_distrib->emplace_back(i - dist_max, (fwd_distrib[std::abs(i - dist_max)]) * PRECISION);
        }
    }

    //Sort distriution forward cumul by value
    std::sort(cumul_fwd_distrib->begin(), cumul_fwd_distrib->end(),
              [](auto a, auto b) {
                  return std::get<1>(a) > std::get<1>(b);
              });

    //Cumul value in cumul distri fw
    for (int i = 0; i <= 2 * dist_max; ++i)
    {
        cumul_proba += std::get<1>(cumul_fwd_distrib->at(i));
        std::get<1>(cumul_fwd_distrib->at(i)) = cumul_proba;
    }

    return cumul_fwd_distrib;
}

//draw one value (distance) from cumulated fwd distri
int cumul_fwd_disp_distrib_c::draw_value(int dim)
{
    int value = 0;

    // at() = prescription for all vector access using a variable (ie not a number)
    if ((Cumul_fwd_distrib.at(dim))->size() > 1)
    {
        int random_numb = Rand_gen->int_0_PRECISION_rand();
        // the vector of tuple is ordered (decreasing prob) for optimization
        // from -dx to +dx, and can be asymetric (asymetry should be set during the computation of the cumul forward distrib)
        for (std::tuple<int, int> const &dist_and_value : *(Cumul_fwd_distrib.at(dim)))
        {
            if (random_numb <= std::get<1>(dist_and_value))
            {
                value = std::get<0>(dist_and_value);
                break;
            }
        }
    }

    return value;
}

//draw two values (distances) from cumulated fwd distri
std::array<int, 2> cumul_fwd_disp_distrib_c::draw_values()
{
    std::array<int, 2> values{0, 0};

    //WARNING : if one day more than 2 dim was implement need to change here
    for (int dim = 0; dim < 2; ++dim)
    {
        // at() = prescription for all vector access using a variable (ie not a number)
        if ((Cumul_fwd_distrib.at(dim))->size() > 1)
        {
            int random_numb = Rand_gen->int_0_PRECISION_rand();
            // the vector of tuple is ordered (decreasing prob) for optimization
            // from -dx to +dx, and can be asymetric (asymetry should be set during the computation of the cumul forward distrib)
            for (std::tuple<int, int> const &dist_and_value : *(Cumul_fwd_distrib.at(dim)))
            {
                if (random_numb <= std::get<1>(dist_and_value))
                {
                    values.at(dim) = std::get<0>(dist_and_value);
                    break;
                }
            }
        }
    }

    return values;
}

/********cumul_bcwd_disp_distrib_c********/

cumul_bcwd_disp_distrib_c::cumul_bcwd_disp_distrib_c(rand_gen_c *rand_gen, node_lattice_c const &node_lati) : Node_lati(node_lati)
{
    Rand_gen = rand_gen;
}

void cumul_bcwd_disp_distrib_c::calc_bcwd_distrib(extend_lattice_c &rmap)
{
    Cumul_bcwd_distrib = rmap.cumul_normalized_bcwd(rmap.compute_migr_nb_reaching_focal_node(Node_lati.Lat.Boundary_effect, Node_lati.Coord));
}

//pull one value (node) from cumulated bcwd distri
node_lattice_c *cumul_bcwd_disp_distrib_c::draw_values()
{
    int random_numb = Rand_gen->int_0_PRECISION_rand();

    for (std::tuple<node_lattice_c *, int> const &node_and_value : Cumul_bcwd_distrib)
    {
        if (random_numb <= std::get<1>(node_and_value))
        {
            return std::get<0>(node_and_value);
        }
    }
    return nullptr;
}

/*************************\
|          NODE           |
\*************************/

//Constructeur
//Node in homogeneous Lattice
// lat is reference to a lattice object, the node should thus be constructed using ' : '
// so that the Lat is filled with the reference during its construction and not by copy after its construction because a ref can not be empty.
node_lattice_c::node_lattice_c(lattice_c &lat, std::array<int, 2> coord, int subpop_size) : Coord(coord), Bcwd_distrib(*this), Lat(lat) //Footnote : Array Coord initialized by copy with this method (avoid non initialize and copy)
{
    Subpop_size = subpop_size;
    Indivs_in_pop.reserve(subpop_size);
}
//Node in heterogeneous Lattice
node_lattice_c::node_lattice_c(lattice_c &lat, std::array<int, 2> coord, int subpop_size, rand_gen_c *rand_gen) : Coord(coord), Bcwd_distrib(rand_gen, *this), Lat(lat)
{
    Subpop_size = subpop_size;
    Indivs_in_pop.reserve(subpop_size);
}

void node_lattice_c::add_indiv(indiv_c *ind)
{
    if (Indivs_in_pop.empty())
    {
        Lat.Nodes_with_lineage.push_back(this);
    }

    if (Indivs_in_pop.capacity() < (Indivs_in_pop.size() + 1))
    {
        Indivs_in_pop.reserve((Indivs_in_pop.size() + 1) * 2);
    }

    ind->Node_lat = this;
    // emplace_back() directly construct the tuple in the vector
    Indivs_in_pop.emplace_back(-1, ind);
}

/*************************\
|         LATTICE         |
\*************************/

//For homogeneous lattice, Homogeneous = true
lattice_c::lattice_c(rand_gen_c &rand_gen, fwd_disp_distrib_c &fwd_distrib, demo_param_c const &demo_param)
{
    Fwd_distrib = std::move(fwd_distrib);
    Cumul_fwd_distrib = cumul_fwd_disp_distrib_c(&rand_gen, Fwd_distrib);

    Boundary_effect = demo_param.Edge_effects;
    Lat_size = demo_param.Lattice_size;
    Disp_dist_max = demo_param.Disp_dist_max;

    assign_homogeneous_subpopsizes(demo_param.Pop_size_per_node, rand_gen);
}

//For heterogeneous lattice with custom subpopsize matrix but with homogeneous dispersal
lattice_c::lattice_c(std::vector<std::vector<int>> const &subpopsize_mat, rand_gen_c &rand_gen, fwd_disp_distrib_c &fwd_distrib, demo_param_c const &demo_param) : Homogeneous(false)
{
    Fwd_distrib = std::move(fwd_distrib);
    Cumul_fwd_distrib = cumul_fwd_disp_distrib_c(&rand_gen, Fwd_distrib);

    Boundary_effect = demo_param.Edge_effects;
    Lat_size = demo_param.Lattice_size;
    Disp_dist_max = demo_param.Disp_dist_max;

    assign_custom_subpopsizes(subpopsize_mat, rand_gen);
}

//For heterogeneous lattice with custom fwd dispersal distributions but homogeneous subpopsizes
lattice_c::lattice_c(rand_gen_c &rand_gen, std::vector<fwd_disp_distrib_c> &fwd_disp_matx, demo_param_c const &demo_param) : Homogeneous(false)
{
    Boundary_effect = demo_param.Edge_effects;
    Lat_size = demo_param.Lattice_size;
    Disp_dist_max = demo_param.Disp_dist_max; // max of max_dists from the different distributions

    assign_homogeneous_subpopsizes(demo_param.Pop_size_per_node, rand_gen);

    auto fwd_disp_matx_itr = fwd_disp_matx.begin();

    for (int x = 0; x <= Lat_size[0]; ++x)
    {
        for (int y = 0; y <= Lat_size[1]; ++y)
        {
            Lattice.back().Fwd_distrib = std::move(*fwd_disp_matx_itr);
            ++fwd_disp_matx_itr;
        }
    }
}

//For heterogeneous lattice with custom migration rate matrix but homogeneous subpopsizes
lattice_c::lattice_c(rand_gen_c &rand_gen, std::vector<std::vector<double>> const &mig_prob_matx, demo_param_c const &demo_param) : Homogeneous(false)
{
    Boundary_effect = demo_param.Edge_effects;
    Lat_size = demo_param.Lattice_size;

    assign_homogeneous_subpopsizes(demo_param.Pop_size_per_node, rand_gen);
    compute_bcwd_disp_distrib(mig_prob_matx);
}

//For heterogeneous lattice with custom migration rate matrix and custom subpopsize matrix
lattice_c::lattice_c(std::vector<std::vector<int>> const &subpopsize_mat, rand_gen_c &rand_gen, std::vector<std::vector<double>> const &mig_prob_matx, demo_param_c const &demo_param) : Homogeneous(false)
{
    Boundary_effect = demo_param.Edge_effects;
    Lat_size = demo_param.Lattice_size;

    assign_custom_subpopsizes(subpopsize_mat, rand_gen);
    compute_bcwd_disp_distrib(mig_prob_matx);
}

void lattice_c::assign_homogeneous_subpopsizes(int const &subpopsize, rand_gen_c &rand_gen)
{
    Lattice.reserve((Lat_size[0] + 1) * (Lat_size[1] + 1));

    for (int x = 0; x <= Lat_size[0]; ++x)
    {
        for (int y = 0; y <= Lat_size[1]; ++y)
        {
            //Bcwd_distrib of each node needs a reference to rand_gen, and a reference can not be empty
            Lattice.emplace_back(*this, std::array<int, 2>{x, y}, subpopsize, &rand_gen);
        }
    }
}

void lattice_c::assign_custom_subpopsizes(std::vector<std::vector<int>> const &subpopsize_mat, rand_gen_c &rand_gen)
{
    Lattice.reserve((Lat_size[0] + 1) * (Lat_size[1] + 1));

    for (int x = 0; x <= Lat_size[0]; ++x)
    {
        for (int y = 0; y <= Lat_size[1]; ++y)
        {
            Lattice.emplace_back(*this, std::array<int, 2>{x, y}, subpopsize_mat[x][y], &rand_gen);
        }
    }
}

void lattice_c::compute_bcwd_disp_distrib(std::vector<std::vector<double>> const &mig_prob_matx)
{
    auto &info_collect = singleton_c<info_collector_c>::instance();
    if (info_collect.Print_migration_matrix_bool)
    {
        info_collect.Bcwd_disp_distrib_from_migration_matrix.resize((Lat_size[0] + 1) * (Lat_size[1] + 1)); // TDODO ? remplacer (Lat_size[0] + 1) * (Lat_size[1] + 1) par Lattice.size() ?
        for (auto it = info_collect.Bcwd_disp_distrib_from_migration_matrix.begin(); it < info_collect.Bcwd_disp_distrib_from_migration_matrix.end(); ++it)
        {
            it->resize((Lat_size[0] + 1) * (Lat_size[1] + 1)); // TDODO ? remplacer (Lat_size[0] + 1) * (Lat_size[1] + 1) par Lattice.size() ?
        }
    }

    //Calc of backward
    int num_focal_node = 0;
    for (int x_node = 0; x_node <= Lat_size[0]; ++x_node)
    {
        for (int y_node = 0; y_node <= Lat_size[1]; ++y_node)
        {
            node_lattice_c *node = this->node({x_node, y_node});
            node->Bcwd_distrib.Cumul_bcwd_distrib.reserve(Lattice.size());

            std::vector<std::tuple<node_lattice_c *, double>> distrib;
            distrib.reserve(Lattice.size());

            int num_migr_node = 0;
            double Immigrant_nbr_sum = 0;
            double Immigrant_nbr;
            for (int x = 0; x <= Lat_size[0]; ++x)
            {
                for (int y = 0; y <= Lat_size[1]; ++y)
                {
                    Immigrant_nbr = this->node({x, y})->Subpop_size * mig_prob_matx[num_migr_node][num_focal_node];
                    distrib.emplace_back(this->node({x, y}), Immigrant_nbr); // TODO : à vérifier avec Tim, c'était mig_prob_matx[num_focal_node][num_migr_node]);
                    Immigrant_nbr_sum += Immigrant_nbr;
                    ++num_migr_node;
                }
            }

            for (auto it1 = distrib.begin(); it1 < distrib.end(); ++it1)
            {
                std::get<1>(*it1) /= Immigrant_nbr_sum;
            }

            if (info_collect.Print_migration_matrix_bool)
            {
                auto it2 = info_collect.Bcwd_disp_distrib_from_migration_matrix.at(num_focal_node).begin();
                for (auto it1 = distrib.begin(); it1 < distrib.end(); ++it1, ++it2)
                {
                    (*it2) = std::get<1>(*it1);
                }
            }

            std::sort(distrib.begin(), distrib.end(),
                      [](auto a, auto b) {
                          return std::get<1>(a) > std::get<1>(b);
                      });

            double cumul = 0;

            for (auto const &node_migr_nbr : distrib)
            {
                cumul += std::get<1>(node_migr_nbr);
                node->Bcwd_distrib.Cumul_bcwd_distrib.emplace_back(std::get<0>(node_migr_nbr), cumul * PRECISION);
            }
            ++num_focal_node;
        }
    }
}

//Nodes have number between 0 and n with (0,0) -> 0 and (x,y) -> n (browse by column)
// the lattice is build/explored in columns bottom up
int lattice_c::hash(std::array<int, 2> const &coord)
{
    return coord[0] * (Lat_size[1] + 1) + coord[1];
}

//Add indiv_c by coord
void lattice_c::add_indiv(indiv_c *ind, std::array<int, 2> const &coord)
{
    auto node_lati = node(coord);
    node_lati->add_indiv(ind);
}

//Access to node by coord
node_lattice_c *lattice_c::node(std::array<int, 2> const &coord)
{
    auto hsh = hash(coord);
    return &(Lattice.at(hsh));
}

/********remap********/

extend_lattice_c::extend_lattice_c(lattice_c &lat) : Lat(lat)
{
    //Map all the "real" coord for each poosible coord. Coord in the lattice stay the same bu coord out will be translate (depend on apply_edge_effect)
    Remap_extend_lattice = std::vector<std::vector<std::array<int, 2>>>(Lat.Lat_size[0] + (Lat.Disp_dist_max[0] * 2) + 1, std::vector<std::array<int, 2>>(Lat.Lat_size[1] + (Lat.Disp_dist_max[1] * 2) + 1, std::array<int, 2>()));

    std::array<int, 2> real_values{-1, -1};

    for (int x = 0; x <= Lat.Lat_size[0] + (Lat.Disp_dist_max[0] * 2); ++x)
    {
        real_values[0] = x - Lat.Disp_dist_max[0];

        for (int y = 0; y <= Lat.Lat_size[1] + (Lat.Disp_dist_max[1] * 2); ++y)
        {
            real_values[1] = y - Lat.Disp_dist_max[1];

            for (int dim = 0; dim < 2; ++dim)
            {
                //Out lattice
                switch (Lat.Boundary_effect)
                {
                case edge_effect_enum::reflecting:
                {
                    while ((real_values.at(dim) < 0) || (real_values.at(dim) > Lat.Lat_size.at(dim)))
                    {
                        if (real_values.at(dim) < 0)
                        {
                            real_values.at(dim) = -real_values.at(dim);
                        }
                        else
                        { // Lat_size.at(dim) - ((real_values.at(dim) - Lat_size.at(dim))
                            real_values.at(dim) = 2 * Lat.Lat_size.at(dim) - real_values.at(dim);
                        }
                    }
                    break;
                }

                case edge_effect_enum::circular:
                {
                    while ((real_values.at(dim) < 0) || (real_values.at(dim) > Lat.Lat_size.at(dim)))
                    {
                        if (real_values.at(dim) < 0)
                        {
                            real_values.at(dim) = Lat.Lat_size.at(dim) + real_values.at(dim) + 1;
                        }
                        else
                        {
                            real_values.at(dim) = -Lat.Lat_size.at(dim) + real_values.at(dim) - 1;
                        }
                    }
                    break;
                }
                //Need to have a value in the lattice (for Neig_migr_nbr)
                case edge_effect_enum::absorbing:
                {
                    if ((real_values.at(dim) < 0) || (real_values.at(dim) > Lat.Lat_size.at(0)))
                    {
                        real_values.at(dim) = -1;
                    }
                    break;
                }
                }
            }
            Remap_extend_lattice[x][y] = real_values;
        }
    }
    Neig_migr_nbr = std::vector<std::vector<double>>(Lat.Lat_size[0] + 1, std::vector<double>(Lat.Lat_size[1] + 1, 0));
}

//Acces to "map"
std::array<int, 2> extend_lattice_c::apply_edge_effect(std::array<int, 2> const &coord)
{
    return Remap_extend_lattice[coord[0] + Lat.Disp_dist_max[0]][coord[1] + Lat.Disp_dist_max[1]];
}

std::vector<std::tuple<node_lattice_c *, double>> extend_lattice_c::compute_migr_nb_reaching_focal_node(edge_effect_enum edg_effect, std::array<int, 2> const &coord_focal_node)
{
    if (edg_effect != edge_effect_enum::absorbing)
    {
        //Browse a square of Disp_dist_max edge size around coord of target node
        for (int ext_lat_x = coord_focal_node[0] - Lat.Disp_dist_max[0]; ext_lat_x <= coord_focal_node[0] + Lat.Disp_dist_max[0]; ++ext_lat_x)
        {
            for (int ext_lat_y = coord_focal_node[1] - Lat.Disp_dist_max[1]; ext_lat_y <= coord_focal_node[1] + Lat.Disp_dist_max[1]; ++ext_lat_y)
            {
                //WARNING = Maybe initialisation will provide problems
                std::array<int, 2> true_lat_coord{0, 0};
                //coord in the real lattice
                true_lat_coord = apply_edge_effect({ext_lat_x, ext_lat_y});
                //Search the "real" node
                auto node_ptr = Lat.node(true_lat_coord);
                fwd_disp_distrib_c *node_fwd_distrib = &(node_ptr->Fwd_distrib);
                //If "real" node's Fwd_distrib not available take lattice's Fwd_distrib
                if (node_fwd_distrib->empty())
                {
                    node_fwd_distrib = &(Lat.Fwd_distrib);
                }

                //See if "real" node can reach target node
                //node_fwd_distrib->size() = dist_max for node's fw distri
                //WARNING only works for symetric forward, need to be adapted for asymetric fwrd
                if (((*node_fwd_distrib).size(0) > static_cast<std::size_t>(abs(coord_focal_node[0] - ext_lat_x))) && ((*node_fwd_distrib).size(1) > static_cast<std::size_t>(abs(coord_focal_node[1] - ext_lat_y))))
                {
                    //P(dist between "false" x and target node x) * P(dist between "false" y and target node y)
                    double prob;

                    if ((abs(coord_focal_node[0] - ext_lat_x) != 0) && (abs(coord_focal_node[1] - ext_lat_y) != 0))
                    {
                        prob = (*node_fwd_distrib)[0][abs(coord_focal_node[0] - ext_lat_x)] * (*node_fwd_distrib)[1][abs(coord_focal_node[1] - ext_lat_y)] / 4;
                    }
                    else
                    {
                        if (abs(coord_focal_node[0] - ext_lat_x) != 0)
                        {
                            prob = (*node_fwd_distrib)[0][abs(coord_focal_node[0] - ext_lat_x)] * (*node_fwd_distrib)[1][abs(coord_focal_node[1] - ext_lat_y)] / 2;
                        }
                        else if (abs(coord_focal_node[1] - ext_lat_y) != 0)
                        {
                            prob = (*node_fwd_distrib)[0][abs(coord_focal_node[0] - ext_lat_x)] * (*node_fwd_distrib)[1][abs(coord_focal_node[1] - ext_lat_y)] / 2;
                        }
                        else
                        {
                            prob = (*node_fwd_distrib)[0][abs(coord_focal_node[0] - ext_lat_x)] * (*node_fwd_distrib)[1][abs(coord_focal_node[1] - ext_lat_y)];
                        }
                    }
                    //create the neighbord's table
                    //WARNING *10-6 for translate in int => optimization purpose, more light to drawn in int and compare int vs int
                    Neig_migr_nbr[true_lat_coord[0]][true_lat_coord[1]] += (prob * node_ptr->Subpop_size);
                }
            }
        }
    }
    else
    {
        //Browse a square of Disp_dist_max edge size around coord of target node
        //if absorbing effect then need to consdier only nodes inside the lattice
        for (int ext_lat_x = coord_focal_node[0] - Lat.Disp_dist_max[0]; ext_lat_x <= coord_focal_node[0] + Lat.Disp_dist_max[0]; ++ext_lat_x)
        {
            if ((ext_lat_x >= 0) && (ext_lat_x <= Lat.Lat_size[0]))
            {
                for (int ext_lat_y = coord_focal_node[1] - Lat.Disp_dist_max[1]; ext_lat_y <= coord_focal_node[1] + Lat.Disp_dist_max[1]; ++ext_lat_y)
                {
                    if ((ext_lat_y >= 0) && (ext_lat_y <= Lat.Lat_size[1]))
                    {
                        auto node_ptr = Lat.node({ext_lat_x, ext_lat_y});
                        fwd_disp_distrib_c *node_fwd_distrib = &(node_ptr->Fwd_distrib);
                        //If node's Fwd_distrib not available take lattice's Fwd_distrib
                        if (node_fwd_distrib->empty())
                        {
                            node_fwd_distrib = &(Lat.Fwd_distrib);
                        }

                        //See if node can reach target node
                        //node_fwd_distrib->size() = dist_max for node's fw distri
                        //WARNING only works for symetric forward, need to be adapted for asymetric fwrd
                        if (((*node_fwd_distrib).size(0) > static_cast<std::size_t>(abs(coord_focal_node[0] - ext_lat_x))) && ((*node_fwd_distrib).size(1) > static_cast<std::size_t>(abs(coord_focal_node[1] - ext_lat_y))))
                        {
                            //P(dist between "false" x and target node x) * P(dist between "false" y and target node y)
                            double prob;
                            if ((abs(coord_focal_node[0] - ext_lat_x) != 0) && (abs(coord_focal_node[1] - ext_lat_y) != 0))
                            {
                                prob = (*node_fwd_distrib)[0][abs(coord_focal_node[0] - ext_lat_x)] * (*node_fwd_distrib)[1][abs(coord_focal_node[1] - ext_lat_y)] / 4;
                            }
                            else
                            {
                                if (abs(coord_focal_node[0] - ext_lat_x) != 0)
                                {
                                    prob = (*node_fwd_distrib)[0][abs(coord_focal_node[0] - ext_lat_x)] * (*node_fwd_distrib)[1][abs(coord_focal_node[1] - ext_lat_y)] / 2;
                                }
                                else if (abs(coord_focal_node[1] - ext_lat_y) != 0)
                                {
                                    prob = (*node_fwd_distrib)[0][abs(coord_focal_node[0] - ext_lat_x)] * (*node_fwd_distrib)[1][abs(coord_focal_node[1] - ext_lat_y)] / 2;
                                }
                                else
                                {
                                    prob = (*node_fwd_distrib)[0][abs(coord_focal_node[0] - ext_lat_x)] * (*node_fwd_distrib)[1][abs(coord_focal_node[1] - ext_lat_y)];
                                }
                            }
                            //create the neighbord's table
                            Neig_migr_nbr[ext_lat_x][ext_lat_y] += (prob * node_ptr->Subpop_size);
                        }
                    }
                }
            }
        }
    }

    //flatten the tab for efficiency
    //NOTE : can be done directly on the flatten tab
    std::vector<std::tuple<node_lattice_c *, double>> flattened_tab;
    //+1 in the end is to permit efficient sorting
    //Max = lattice size
    flattened_tab.reserve(((2 * Lat.Disp_dist_max[0]) + 1) * ((2 * Lat.Disp_dist_max[1]) + 1) + 1);

    if (edg_effect != edge_effect_enum::absorbing)
    {
        for (int ext_lat_x = coord_focal_node[0] - Lat.Disp_dist_max[0]; ext_lat_x <= coord_focal_node[0] + Lat.Disp_dist_max[0]; ++ext_lat_x)
        {
            for (int ext_lat_y = coord_focal_node[1] - Lat.Disp_dist_max[1]; ext_lat_y <= coord_focal_node[1] + Lat.Disp_dist_max[1]; ++ext_lat_y)
            {
                //WARNING = Maybe initialisation will provide problems
                std::array<int, 2> true_lat{0, 0};
                true_lat = apply_edge_effect({ext_lat_x, ext_lat_y});
                auto node_ptr = Lat.node(true_lat);
                //flattening freq neig node

                if (Neig_migr_nbr[true_lat[0]][true_lat[1]] != 0)
                {
                    flattened_tab.emplace_back(node_ptr, Neig_migr_nbr[true_lat[0]][true_lat[1]]);
                    //Permit to scan multiple time the same "real" node without mess in freq AND to reuse Neig_migr_nbr
                    Neig_migr_nbr[true_lat[0]][true_lat[1]] = 0;
                }
            }
        }
    }
    else
    {
        for (int ext_lat_x = coord_focal_node[0] - Lat.Disp_dist_max[0]; ext_lat_x <= coord_focal_node[0] + Lat.Disp_dist_max[0]; ++ext_lat_x)
        {
            if ((ext_lat_x >= 0) && (ext_lat_x <= Lat.Lat_size[0]))
            {
                for (int ext_lat_y = coord_focal_node[1] - Lat.Disp_dist_max[1]; ext_lat_y <= coord_focal_node[1] + Lat.Disp_dist_max[1]; ++ext_lat_y)
                {
                    if ((ext_lat_y >= 0) && (ext_lat_y <= Lat.Lat_size[1]))
                    {
                        auto node_ptr = Lat.node({ext_lat_x, ext_lat_y});
                        //flattening freq neig node

                        if (Neig_migr_nbr[ext_lat_x][ext_lat_y] != 0)
                        {
                            flattened_tab.emplace_back(node_ptr, Neig_migr_nbr[ext_lat_x][ext_lat_y]);
                            //Permit to scan multiple time the same "real" node without mess in freq AND to reuse Neig_migr_nbr
                            Neig_migr_nbr[ext_lat_x][ext_lat_y] = 0;
                        }
                    }
                }
            }
        }
    }

    //Sort A > B
    std::sort(flattened_tab.begin(), flattened_tab.end(),
              [](auto a, auto b) {
                  return std::get<1>(a) > std::get<1>(b);
              });

    //Cut all the 0
    flattened_tab.shrink_to_fit();

    return flattened_tab;
}

std::vector<std::tuple<node_lattice_c *, int>> extend_lattice_c::cumul_normalized_bcwd(std::vector<std::tuple<node_lattice_c *, double>> non_cumul_bcwd)
{
    std::vector<std::tuple<node_lattice_c *, int>> int_non_cumul_bcwd(non_cumul_bcwd.size());
    //Cumul
    for (std::size_t i = 1; i < non_cumul_bcwd.size(); ++i)
    {
        std::get<1>(non_cumul_bcwd[i]) += std::get<1>(non_cumul_bcwd[i - 1]);
    }

    double cumul = std::get<1>(non_cumul_bcwd[non_cumul_bcwd.size() - 1]);
    //Normalize between 0 and 1
    auto int_non_cumul_bcwd_itr = int_non_cumul_bcwd.begin();
    for (auto &node_and_value : non_cumul_bcwd)
    {
        std::get<0>(*int_non_cumul_bcwd_itr) = std::get<0>(node_and_value);
        //WARNING *10-6 for translate in int => optimization purpose, more light to drawn in int and compare int vs int
        std::get<1>(*int_non_cumul_bcwd_itr) = (std::get<1>(node_and_value) / cumul) * PRECISION;
        ++int_non_cumul_bcwd_itr;
    }

    return int_non_cumul_bcwd;
}

void add_indiv_at_node(indiv_c *indiv, node_lattice_c *node)
{
    node->add_indiv(indiv);
    indiv->Node_lat = node;
}

/*************************\
|        MIGRATION        |
\*************************/

void migration(lattice_c &lat, extend_lattice_c *rmap, indiv_stock_c &indiv_stock)
{
    if (lat.Homogeneous)
    {
        homogene_migration(lat, rmap, indiv_stock);
    }
    else
    {
        heterogene_migration(lat, rmap, indiv_stock);
    }
}
//TODO : temporary
#include "summary_stat.hpp"
void homogene_migration(lattice_c &lat, extend_lattice_c *rmap, indiv_stock_c &indiv_stock)
{
    for (auto indiv : indiv_stock)
    {
        std::array<int, 2> dist{0, 0};
        dist = lat.Cumul_fwd_distrib.draw_values();
        // Compute the empirical immigration rate and empirical dispersal distribution
        // singleton (objet global) pour collecter des info
        auto &info_collect = singleton_c<info_collector_c>::instance();
        auto const &demo_param = singleton_c<demo_param_c>::instance();
        if (info_collect.Effective_disp)
        {
            for (size_t dim = 0; dim < (!(lat.Lat_size[1]) ? 1 : 2); ++dim)
            {
                if (dist[dim] == 0)
                {
                    ++(info_collect.Nbr_depl_0_1_tot[0]);
                }
                else
                {
                    ++(info_collect.Nbr_depl_0_1_tot[1]);
                }
                ++(info_collect.Nbr_depl_0_1_tot[2]);
                ++(info_collect.Emp_cumul_axial_disp[dim].at(dist[dim] + demo_param.Disp_dist_max[dim]));
            }
        }
        //needed because the lattice(s) have been emptied before the migration step
        node_lattice_c *destination_node = indiv->Node_lat;
        if ((dist[0] != 0) || (dist[1] != 0))
        {
            destination_node = apply_movment(lat, *rmap, indiv, dist);
        }
        destination_node->add_indiv(indiv);
    }
}

node_lattice_c *apply_movment(lattice_c &lat, extend_lattice_c &rmap, indiv_c const *indiv, std::array<int, 2> const &dist)
{
    auto coord = (indiv->Node_lat)->Coord;

    std::array<int, 2> new_coord = {coord[0] + dist[0], coord[1] + dist[1]};

    if (lat.Boundary_effect != edge_effect_enum::absorbing)
    {
        new_coord = rmap.apply_edge_effect(new_coord);
    }
    else
    { // redraw for the absorbing case
        for (int dim = 0; dim < 2; ++dim)
        {
            while ((new_coord.at(dim) < 0) || (new_coord.at(dim) > lat.Lat_size.at(dim)))
            {
                new_coord.at(dim) = coord.at(dim) + lat.Cumul_fwd_distrib.draw_value(dim);
            }
        }
    }

    return lat.node(new_coord);
}

// for heterogenous migration, directly draw a node and not a distance
void heterogene_migration(lattice_c &lat, extend_lattice_c *rmap, indiv_stock_c &indiv_stock)
{
    for (auto &indiv : indiv_stock)
    {
        auto &bcwd_distrib = (indiv->Node_lat)->Bcwd_distrib;
        if (bcwd_distrib.empty())
        {
            bcwd_distrib.calc_bcwd_distrib(*rmap);
        }
        node_lattice_c *destination_node = bcwd_distrib.draw_values();
        // Compute the empirical immigration rate and empirical dispersal distribution
        // singleton (objet global) pour collecter des info
        auto &info_collect = singleton_c<info_collector_c>::instance();
        auto const &demo_param = singleton_c<demo_param_c>::instance();
        if (info_collect.Effective_disp)
        {
            std::array<int, 2> from_coord = indiv->Node_lat->Coord, to_coord = destination_node->Coord;

            for (size_t dim = 0; dim <= static_cast<size_t>(lat.Lat_size[1]); ++dim)
            {
                int dist = from_coord[dim] - to_coord[dim];
                if (dist == 0)
                {
                    ++(info_collect.Nbr_depl_0_1_tot[0]);
                }
                else
                {
                    ++(info_collect.Nbr_depl_0_1_tot[1]);
                }
                ++(info_collect.Nbr_depl_0_1_tot[2]);
                ++(info_collect.Emp_cumul_axial_disp[dim].at(dist + demo_param.Disp_dist_max[dim]));
            }
        }

        destination_node->add_indiv(indiv);
    }
}
