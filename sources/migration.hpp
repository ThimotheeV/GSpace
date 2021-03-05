#pragma once
#include <array>
#include <memory>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <iostream>
#include <functional>
#include <tuple>

#include "segment.hpp"
#include "random_generator.hpp"
#include "settings.hpp"

enum class edge_effect_enum
{
    reflecting,
    absorbing,
    circular
};

struct map_string_edge_effect_name_c
{
    std::map<edge_effect_enum, std::string> ReverseMap{
        {edge_effect_enum::reflecting, "reflecting"},
        {edge_effect_enum::absorbing, "absorbing"},
        {edge_effect_enum::circular, "circular"}};
};

//Pre declaration du to the include mess with compilator
// so that the compiler knows the object name even if the object is not defined before being referenced somewhere
// the object may then be define later in this file or in another (included later)
struct indiv_c;
struct indiv_stock_c;
struct lattice_c;
struct node_lattice_c;
struct extend_lattice_c;
struct demo_param_c;

/*************************\
|       DISPERSAL         |
|      DISTRIBUTION       |
\*************************/
struct fwd_disp_distrib_c
{
    static std::function<double(int)> uniform();
    static std::function<double(int)> geometric(double g_geo_param);
    //Not a pareto => Zeta distribution (discret equivalent of pareto)
    static std::function<double(int)> pareto(double n_pareto_param);
    static std::function<double(int)> gaussian(int dist_max);
    static std::function<double(int)> gamma_or_sichel(demo_param_c const &demo_p, int dim);

    fwd_disp_distrib_c() = default;
    //explicit = disable the copy constructor.
    //It is a good practice when declaring a constructor with a single param
    explicit fwd_disp_distrib_c(demo_param_c const &demo_param);
    //build the dispersal distribution
    //the shared_ptr is automatically freed when it is not used = not referenced anywhere (when there is no more variable with its adress)
    std::shared_ptr<std::vector<double>> construct_fwd_disp_distrib(int dist_max, double proba_mig, std::function<double(int)> const &calc_proba);

    bool empty()
    {
        return Fwd_distrib.at(0) == nullptr;
    }
    //disp_dist_max can be different in both dimensions
    std::size_t size(std::size_t i)
    {
        return Fwd_distrib.at(i)->size();
    }
    //returns a reference, thus the object can be modified =  not constant method
    std::vector<double> &operator[](std::size_t i)
    {
        return *(Fwd_distrib.at(i));
    }
    //returns a copy, thus only the copy can be modified, not the original object = constant method
    std::vector<double> const &operator[](std::size_t i) const
    {
        return *(Fwd_distrib.at(i));
    }

private:
    //ptr to the vectors containing the dispersal distribution, could be private because generally not used ([] is used instead)
    std::array<std::shared_ptr<std::vector<double>>, 2> Fwd_distrib{nullptr, nullptr};
};

//for homogeneous lattice or homogeneous area
//computes cumulative forward distribution from a fwd_disp_distrib object
//fw can be an attribute of a lattice or of a node
struct cumul_fwd_disp_distrib_c
{
    cumul_fwd_disp_distrib_c() = default;
    cumul_fwd_disp_distrib_c(rand_gen_c *rand_gen, fwd_disp_distrib_c const &fwd_distrib);
    std::shared_ptr<std::vector<std::tuple<int, int>>> construct_cumul_fwd_disp_distrib(std::vector<double> const &fwd_distrib);
    int draw_value(int dim);          // draw a distance from the cumul distrib
    std::array<int, 2> draw_values(); // draw two distances from cumul distrib (for optimisation reasons)

    //shared pointer because can be the same for the two dimensions
    // shared -> avoid to manually free the pointer(s)
    // the vector of tuple is ordered (decreasing prob) for optimization
    // from -dx to +dx, and can be asymetric (asymetry should be set during the computation of the cumul forward distrib)
    std::array<std::shared_ptr<std::vector<std::tuple<int, int>>>, 2> Cumul_fwd_distrib{nullptr, nullptr};
    rand_gen_c *Rand_gen{nullptr};
};

// draw a pointer to a node, not a distance
struct cumul_bcwd_disp_distrib_c
{
    // no default constructor because the class has a ref and a ref can not be empty

    // explicit = prescription for a constructor with a single arg, tells that it is not a constructor by copy.
    explicit cumul_bcwd_disp_distrib_c(node_lattice_c const &node_lati) : Node_lati(node_lati){};

    cumul_bcwd_disp_distrib_c(rand_gen_c *rand_gen, node_lattice_c const &node_lati);

    //reprendre ici
    void calc_bcwd_distrib(extend_lattice_c &rmap);

    // draw a pointer to a node, not a distance
    node_lattice_c *draw_values();

    bool empty()
    {
        return Cumul_bcwd_distrib.empty();
    }

    std::vector<std::tuple<node_lattice_c *, int>> Cumul_bcwd_distrib;
    rand_gen_c *Rand_gen{nullptr};
    // ref to the node to whch the bcwd is attached
    // needed to construct the backwd only when necessary
    node_lattice_c const &Node_lati;
};

/*************************\
|          NODE           |
\*************************/
struct node_lattice_c
{
    //Constructor
    //Homogeneous lattice
    node_lattice_c(lattice_c &lat, std::array<int, 2> coord, int subpop_size);
    //Heterogeneous lattice
    node_lattice_c(lattice_c &lat, std::array<int, 2> coord, int subpop_size, rand_gen_c *rand_gen);

    void add_indiv(indiv_c *ind);

    std::array<int, 2> Coord;
    //Uniq_ptr pas possible, nécéssite un transfert de l'ownership sur Indivs_in_pop (pas sûr que ce soit intérressant !)
    std::vector<std::tuple<int, indiv_c *>> Indivs_in_pop; // mostly for coalescence, the int is the nbr of its parent at the current generation
    int Subpop_size{0};
    fwd_disp_distrib_c Fwd_distrib;         // empty when heterogeneous lattice
    cumul_bcwd_disp_distrib_c Bcwd_distrib; // empty when homogeneous lattice
    lattice_c &Lat;                         // ref to its lattice
};

/*************************\
|         LATTICE         |
\*************************/
struct lattice_c
{
    //empty constructor
    lattice_c() = default ;
    //Constructor fo Homogeneous (subpopsize & disp) lattice
    lattice_c(rand_gen_c &rand_gen, fwd_disp_distrib_c &fwd_distrib, demo_param_c const &demo_param);
    //Constructor for heterogeneouslattice with custom subpopsize matrix but homogeneous dispersal
    lattice_c(std::vector<std::vector<int>> const &subpopsize_mat, rand_gen_c &rand_gen, fwd_disp_distrib_c &fwd_distrib, demo_param_c const &demo_param);
    //Constructor for lattice with heterogeneous dispersal distributions
    lattice_c(rand_gen_c &rand_gen, std::vector<fwd_disp_distrib_c> &fwd_disp_matx, demo_param_c const &demo_param);
    //Constructor for lattice with custom migration (m) matrix
    // ??? mig_prob_matx not cumul, not sorted but normalized
    lattice_c(rand_gen_c &rand_gen, std::vector<std::vector<double>> const &mig_prob_matx, demo_param_c const &demo_param);
    //Constructor for lattice with custom migration (m) matrix and custom subpopsize matrix
    lattice_c(std::vector<std::vector<int>> const &subpopsize_mat, rand_gen_c &rand_gen, std::vector<std::vector<double>> const &mig_prob_matx, demo_param_c const &demo_param);

    //2 functions used in constructors to fill subpopsizes from custom matrix and backward dispersal distribution from migration matrix
    void assign_homogeneous_subpopsizes(int const &subpopsize, rand_gen_c &rand_gen);
    void assign_custom_subpopsizes(std::vector<std::vector<int>> const &subpopsize_mat, rand_gen_c &rand_gen);
    void compute_bcwd_disp_distrib(std::vector<std::vector<double>> const &mig_prob_matx);
    //transform 2 coords into one unique int value
    int hash(std::array<int, 2> const &coord);
    // add an individual on a coord
    void add_indiv(indiv_c *ind, std::array<int, 2> const &coord);
    //access to node by coord
    node_lattice_c *node(std::array<int, 2> const &coord);

    bool Homogeneous{true}; // reserve the mem zone and put true during declaration (i.e. before construction)
    std::array<int, 2> Lat_size{1, 1};
    std::array<int, 2> Disp_dist_max{0, 0}; // max of max_dists if heterogeneous dispersal

    std::vector<node_lattice_c> Lattice;
    std::vector<node_lattice_c *> Nodes_with_lineage; //only nodes with lineages (for coalescence process)

    edge_effect_enum Boundary_effect{edge_effect_enum::reflecting};
    fwd_disp_distrib_c Fwd_distrib;
    cumul_fwd_disp_distrib_c Cumul_fwd_distrib;
};

struct extend_lattice_c
{
    //Used to apply edge effects on fw disp for homogeneous lattice or to calc bwk dispers distri for heterogeneous lattice
    explicit extend_lattice_c(lattice_c &lat);

    //lattice_to_lattice or neig_to_lat = all nodes than can reach the focal node
    std::array<int, 2> apply_edge_effect(std::array<int, 2> const &coord);

    //all nodes than can reach the focal node
    std::vector<std::tuple<node_lattice_c *, double>> compute_migr_nb_reaching_focal_node(edge_effect_enum edg_effect, std::array<int, 2> const &coord_focal_node);
    std::vector<std::tuple<node_lattice_c *, int>> cumul_normalized_bcwd(std::vector<std::tuple<node_lattice_c *, double>> non_cumul_bcwd);

    //Need to browse by y and not by x (cache direction) !
    std::vector<std::vector<std::array<int, 2>>> Remap_extend_lattice;
    //Representation of lattice (same size)
    // squarre of minimal size with nodes that can reach the focal node
    std::vector<std::vector<double>> Neig_migr_nbr;
    //ref to the the real lattice, needed to reach the different nodes
    //should allow the consideration of multiple lattices
    lattice_c &Lat;
};

//"outside" the lattice
void add_indiv_at_node(indiv_c *indiv, node_lattice_c *node);

/*************************\
|        MIGRATION        |
\*************************/
void migration(lattice_c &lat, extend_lattice_c *rmap, indiv_stock_c &indiv_stock);

void homogene_migration(lattice_c &lat, extend_lattice_c *rmap, indiv_stock_c &indiv_stock);
node_lattice_c *apply_movment(lattice_c &lat, extend_lattice_c &rmap, indiv_c const *indiv, std::array<int, 2> const &dist);

void heterogene_migration(lattice_c &lat, extend_lattice_c *rmap, indiv_stock_c &indiv_stock);
