#include <cassert>

#include "kelleher_algo.hpp"

//This implementation come from 2016 Kelleher paper "Efficient Coalescence Simulation and genealogical Analysis for Large Sample Sizes"
// See "simulation.cpp" for the global algorithm"

void struct_arg_c::ini(indiv_stock_c const &indiv_vec, samp_param_c const &samp_param, int chr_index)
{
    W_next_node_nbr = 0;

    //Change nbr dependant of Ploidy
    auto chr_nbr_per_loc = samp_param.Ploidy * samp_param.n_total_sample_size;

    //L_cumul_nbr_brkpt_per_seg contains number of brkpt_site per seg, see "fenwick_tree.hpp"
    L_cumul_nbr_brkpt_per_seg.create_fenwick_tree(chr_nbr_per_loc);

    for (auto &indiv_ptr : indiv_vec)
    {
        // Here "i" is the node number and "index" the seg identifier
        auto seg_ptr = Pool_seg.get_new_obj();
        auto index = Pool_seg.get_index_from_obj(seg_ptr);
        L_cumul_nbr_brkpt_per_seg.set_value_at_index(index, samp_param.Sequence_length - 1);

        *seg_ptr = seg(0, samp_param.Sequence_length, W_next_node_nbr++, nullptr, nullptr, indiv_ptr);
        seg *seg_ptr2 = nullptr;

        if (samp_param.Ploidy == 2)
        {
            seg_ptr2 = Pool_seg.get_new_obj();
            auto index2 = Pool_seg.get_index_from_obj(seg_ptr2);
            L_cumul_nbr_brkpt_per_seg.set_value_at_index(index2, samp_param.Sequence_length - 1);

            *seg_ptr2 = seg(0, samp_param.Sequence_length, W_next_node_nbr++, nullptr, nullptr, indiv_ptr);
        }
        indiv_ptr->update_indiv(chr_index, seg_ptr, seg_ptr2);
    }
    //S_intersection_count represents the number of seg present in [recomb_site, next_key_in_map[
    S_intersection_count.emplace(0, chr_nbr_per_loc);
    S_intersection_count.emplace(samp_param.Sequence_length, -1);

    T_time = 0;
}

//If true recomb event, else coa event
bool struct_arg_c::choose_event(double scaled_recomb_rate_rho, indiv_stock_c &indiv_vec, rand_gen_c &rand_gen)
{
    double lambda_r = scaled_recomb_rate_rho * L_cumul_nbr_brkpt_per_seg.get_tot_cumul_freq();
    double lambda = lambda_r + indiv_vec.size() * (indiv_vec.size() - 1);

    std::exponential_distribution<double> exponential_distribution(lambda);

    T_time += exponential_distribution(rand_gen.Seed_gen);

    return rand_gen.int_0_PRECISION_rand() < (lambda_r / lambda) * PRECISION;
}

/**************************************************************/
/*                      Reco algorithm                        */
/**************************************************************/

bool alg_recomb_c::ini(seg *y_choosen_seg, int k_choose_brkpt)
{
    Y_choosen_seg = y_choosen_seg;
    K_choose_brkpt = k_choose_brkpt;
    X_prev_seg = Y_choosen_seg->Previous_seg;

    //If Y_choosen_seg->L_left_brkpt => K_choose_brkpt break is between Y_choosen_seg and X_prev_seg
    return (Y_choosen_seg->L_left_brkpt < K_choose_brkpt);
}

void alg_recomb_c::break_between_segment()
{
    if (X_prev_seg != nullptr)
    {
        X_prev_seg->Next_seg = nullptr;
    }
    Y_choosen_seg->Previous_seg = nullptr;
    Z_new_seg = Y_choosen_seg;
}

void alg_recomb_c::break_within_segment(obj_mem_manag_c<seg> &pool_seg, fenwick_tree_c &cumul_nbr_brkpt_per_seg)
{
    Z_new_seg = pool_seg.get_new_obj();
    seg *y_next_seg = Y_choosen_seg->Next_seg;
    *Z_new_seg = seg(K_choose_brkpt, Y_choosen_seg->R_right_brkpt, Y_choosen_seg->U_ident_node, nullptr, y_next_seg);

    if (y_next_seg != nullptr)
    {
        y_next_seg->Previous_seg = Z_new_seg;
    }
    //break chain of seg for Y_choosen_seg
    Y_choosen_seg->Next_seg = nullptr;
    Y_choosen_seg->R_right_brkpt = K_choose_brkpt;
    X_prev_seg = Y_choosen_seg;

    cumul_nbr_brkpt_per_seg.up_at_index_with_value(pool_seg.get_index_from_obj(Y_choosen_seg), K_choose_brkpt - Z_new_seg->R_right_brkpt);
}

std::array<seg *, 2> alg_recomb_c::update_population(obj_mem_manag_c<seg> &pool_seg, fenwick_tree_c &cumul_nbr_brkpt_per_seg)
{
    auto index = pool_seg.get_index_from_obj(Z_new_seg);
    //Don't keeps track of left brkpt because they have no interest
    cumul_nbr_brkpt_per_seg.set_value_at_index(index, Z_new_seg->R_right_brkpt - Z_new_seg->L_left_brkpt - 1);
    return {X_prev_seg, Z_new_seg};
}

/**************************************************************/
/*                 Coalescence algorithm                      */
/**************************************************************/

void alg_coa_c::choose_ancestors(seg *seg1, seg *seg2)
{
    X_prev_seg = seg1;
    Y_choosen_seg = seg2;
}

//True means X_prev_seg and Y_choosen_seg overlap and have equal L_left_brkpt
bool alg_coa_c::choose_case(obj_mem_manag_c<seg> &pool_seg)
{
    bool overlap;
    //X_prev_seg->L_left_brkpt need to be smaller than Y_choosen_seg->L_left_brkpt for the rest of algorithm
    if (Y_choosen_seg->L_left_brkpt < X_prev_seg->L_left_brkpt)
    {
        seg *temp_seg;
        temp_seg = X_prev_seg;
        X_prev_seg = Y_choosen_seg;
        Y_choosen_seg = temp_seg;
    }

    //X_prev_seg and Y_choosen_seg don't overlap.
    if (X_prev_seg->R_right_brkpt <= Y_choosen_seg->L_left_brkpt)
    {
        Alpha_new_seg = X_prev_seg;
        X_prev_seg = X_prev_seg->Next_seg;
        Alpha_new_seg->Next_seg = nullptr;
        overlap = false;
    }
    //First part of X_prev_seg and Y_choosen_seg don't overlap. This part can be look like a seg who don't overlap at all.
    else if (X_prev_seg->L_left_brkpt != Y_choosen_seg->L_left_brkpt)
    {
        Alpha_new_seg = pool_seg.get_new_obj();
        *Alpha_new_seg = seg(X_prev_seg->L_left_brkpt, Y_choosen_seg->L_left_brkpt, X_prev_seg->U_ident_node);
        X_prev_seg->L_left_brkpt = Y_choosen_seg->L_left_brkpt;
        overlap = false;
    }
    else
    {
        overlap = true;
    }

    return overlap;
}

//True means they are two segs in this overlapping area
bool alg_coa_c::coa(std::map<int, int> &intersection_count, int &next_node_ident)
{
    if (!(Bool_coa))
    {
        Bool_coa = true;
        next_node_ident += 1;
    }

    U_ident_node = next_node_ident - 1;
    L_left_brkpt = X_prev_seg->L_left_brkpt;
    //Search for overlapping area
    Rp_min_right_brkpt = std::min(X_prev_seg->R_right_brkpt, Y_choosen_seg->R_right_brkpt);

    //If intersection_count's node for L_left_brkpt and Rp_min_right_brkpt aren't up
    //  create them with the value of their previous nodes (value == Intersection_counter for their areas)
    auto node_j = intersection_count.find(L_left_brkpt);
    if (node_j == intersection_count.end())
    {
        node_j = std::prev(intersection_count.upper_bound(L_left_brkpt));
        intersection_count.emplace(L_left_brkpt, node_j->second);
    }

    node_j = intersection_count.find(Rp_min_right_brkpt);
    if (node_j == intersection_count.end())
    {
        node_j = std::prev(intersection_count.upper_bound(Rp_min_right_brkpt));
        intersection_count.emplace(Rp_min_right_brkpt, node_j->second);
    }

    return (intersection_count.at(L_left_brkpt) != 2);
}

//Allow to skip decrement_overlaps and thus create a Alpha_new_seg for this MRCA area
//Include all area with R_right_brkpt for set coalescence event
//Can happen only once because by definition there is no more overlap in the area
void alg_coa_c::seg_mrca(std::map<int, int> &intersection_count)
{
    intersection_count.at(L_left_brkpt) = 0;
    R_right_brkpt = intersection_count.upper_bound(L_left_brkpt)->first;
}

//Reduced overlap for all nodes in this area and create a new segment representing the fusion of the two previous segs
void alg_coa_c::decrement_overlaps(obj_mem_manag_c<seg> &pool_seg, std::map<int, int> &intersection_count)
{
    R_right_brkpt = L_left_brkpt;
    auto Node_right_r = intersection_count.find(R_right_brkpt);
    while ((Node_right_r != intersection_count.end()) && (Node_right_r->second != 2) && (R_right_brkpt < Rp_min_right_brkpt))
    {
        Node_right_r->second -= 1;
        R_right_brkpt = intersection_count.upper_bound(R_right_brkpt)->first;
        Node_right_r = intersection_count.find(R_right_brkpt);
    }
    Alpha_new_seg = pool_seg.get_new_obj();
    *Alpha_new_seg = seg(L_left_brkpt, R_right_brkpt, U_ident_node);
}

void alg_coa_c::update_segs(obj_mem_manag_c<seg> &pool_seg, fenwick_tree_c &cumul_nbr_brkpt_per_seg, coa_table_c &coa_table, double time_t, int gen)
{
    //parent_node X_prev_seg or Y_choosen_seg is already a parent just need to update coalescence table.
    // if(U_ident_node >= 193)
    // {
    //     std::cout<<L_left_brkpt<<", "<<R_right_brkpt<<" "<<U_ident_node<<" "<<X_prev_seg->U_ident_node<<" "<<Y_choosen_seg->U_ident_node<<" "<<gen<<std::endl;
    // }
    coa_table.set_coa_event(L_left_brkpt, R_right_brkpt, U_ident_node, {X_prev_seg->U_ident_node, Y_choosen_seg->U_ident_node}, time_t, gen);

    //If area go to the end of the X_prev_seg replace it by its next one. Delete X_prev_seg of simulation.
    if (X_prev_seg->R_right_brkpt == R_right_brkpt)
    {
        seg *Next_seg_x = X_prev_seg->Next_seg;
        auto index = pool_seg.get_index_from_obj(X_prev_seg);

        cumul_nbr_brkpt_per_seg.set_value_at_index(index, 0);
        pool_seg.free_obj(X_prev_seg);
        X_prev_seg = Next_seg_x;
    }
    //else shift begin of X_prev_seg at the end of the area
    else
    {
        X_prev_seg->L_left_brkpt = R_right_brkpt;
    }

    //Idem for Y_choosen_seg
    if (Y_choosen_seg->R_right_brkpt == R_right_brkpt)
    {
        seg *Next_seg_y = Y_choosen_seg->Next_seg;
        auto index = pool_seg.get_index_from_obj(Y_choosen_seg);

        cumul_nbr_brkpt_per_seg.set_value_at_index(index, 0);
        pool_seg.free_obj(Y_choosen_seg);
        Y_choosen_seg = Next_seg_y;
    }
    else
    {
        Y_choosen_seg->L_left_brkpt = R_right_brkpt;
    }
}

void alg_coa_c::update_brkpts(obj_mem_manag_c<seg> &pool_seg, fenwick_tree_c &cumul_nbr_brkpt_per_seg)
{
    //Alpha doesn't exist for MRCA case, this area disappears from the chain of the parent
    if (Alpha_new_seg == nullptr)
    {
        return;
    }

    std::size_t index;
    //Alpha represents the first seg of the new individual
    if (Z_new_seg == nullptr)
    {
        index = pool_seg.get_index_from_obj(Alpha_new_seg);
        New_first_seg = Alpha_new_seg;
        cumul_nbr_brkpt_per_seg.set_value_at_index(index, Alpha_new_seg->R_right_brkpt - Alpha_new_seg->L_left_brkpt - 1);
    }
    else
    {
        //If bool true two segs can be merge into one
        Bool_defrag |= ((Z_new_seg->R_right_brkpt == Alpha_new_seg->L_left_brkpt) && (Z_new_seg->U_ident_node == Alpha_new_seg->U_ident_node));
        Z_new_seg->Next_seg = Alpha_new_seg;
        index = pool_seg.get_index_from_obj(Alpha_new_seg);
        //With Alpha_new_seg->R_right_brkpt - Z_new_seg->R_right_brkpt, FW(Seg.index) keeps track of all the brkpt between Prev(Seg).right and Seg.left
        cumul_nbr_brkpt_per_seg.set_value_at_index(index, Alpha_new_seg->R_right_brkpt - Z_new_seg->R_right_brkpt);
    }
    //Store Alpha_new_seg has previous seg
    Alpha_new_seg->Previous_seg = Z_new_seg;
    Z_new_seg = Alpha_new_seg;
}

//To limit the time of node access if two consecutive nodes have the same  value, one of them was deleted
void alg_coa_c::prune_tree(std::map<int, int> &intersection_count)
{
    auto node_j = intersection_count.begin();

    while (node_j != intersection_count.end())
    {
        auto node_k = std::next(node_j);
        if (node_j->second == node_k->second)
        {
            intersection_count.erase(node_k);
        }
        node_j = std::next(node_j);
    }
}

//To limit the time of seg browse if two consecutive segs (x,y) in the same invidual have x.end == y.begin, they are merged
void alg_coa_c::reshape_segs(obj_mem_manag_c<seg> &pool_seg, fenwick_tree_c &cumul_nbr_brkpt_per_seg)
{
    Y_choosen_seg = Z_new_seg;
    if (Y_choosen_seg == nullptr)
    {
        return;
    }
    while (Y_choosen_seg->Previous_seg != nullptr)
    {
        X_prev_seg = Y_choosen_seg->Previous_seg;
        if ((X_prev_seg->R_right_brkpt == Y_choosen_seg->L_left_brkpt) && (X_prev_seg->U_ident_node == Y_choosen_seg->U_ident_node))
        {
            X_prev_seg->R_right_brkpt = Y_choosen_seg->R_right_brkpt;
            X_prev_seg->Next_seg = Y_choosen_seg->Next_seg;
            if (Y_choosen_seg->Next_seg != nullptr)
            {
                Y_choosen_seg->Next_seg->Previous_seg = X_prev_seg;
            }

            auto index_seg_x = pool_seg.get_index_from_obj(X_prev_seg);
            auto index_seg_y = pool_seg.get_index_from_obj(Y_choosen_seg);
            cumul_nbr_brkpt_per_seg.up_at_index_with_value(index_seg_x, Y_choosen_seg->R_right_brkpt - Y_choosen_seg->L_left_brkpt);
            cumul_nbr_brkpt_per_seg.set_value_at_index(index_seg_y, 0);

            pool_seg.free_obj(Y_choosen_seg);
        }
        Y_choosen_seg = X_prev_seg;
    }
}

//For swap_chain in simulator.cpp
//have a order last_chr1 <-> first_chr2
void simple_coa(struct_arg_c &struct_arg, seg *chr1, seg *first_chr2)
{
    while (chr1->Next_seg != nullptr)
    {
        chr1 = chr1->Next_seg;
    }

    if ((chr1->R_right_brkpt > first_chr2->L_left_brkpt) || (first_chr2->Previous_seg != nullptr))
    {
        throw std::logic_error("( In simple_coa() : Problem in swap chain because of non coherant segments. Contact the developpers. I exit. )");
    }

    chr1->Next_seg = first_chr2;
    first_chr2->Previous_seg = chr1;

    auto index = struct_arg.Pool_seg.get_index_from_obj(first_chr2);
    struct_arg.L_cumul_nbr_brkpt_per_seg.set_value_at_index(index, first_chr2->R_right_brkpt - chr1->R_right_brkpt);
}
