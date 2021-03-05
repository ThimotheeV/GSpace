#include <numeric>
#include <algorithm>
#include <iostream>

#include "coalescence_table.hpp"
#include "settings.hpp"

/**************************************************************/
/*                 Coalescence table                          */
/**************************************************************/

void coa_table_c::set_coa_event(int left_brkpt_l, int right_brkpt_r, int num_node_u, std::vector<int> childs_node, double T_time, int gen)
{
    auto event = std::make_tuple(left_brkpt_l, right_brkpt_r, num_node_u, childs_node, T_time, gen);
    if (Coalescence_table.capacity() <= Coalescence_table.size())
    {
        Coalescence_table.reserve(Coalescence_table.capacity() * 2);
    }
    Coalescence_table.push_back(event);
}

void coa_table_c::set_coa_event(coa_event_t event)
{
    if (Coalescence_table.capacity() <= Coalescence_table.size())
    {
        Coalescence_table.reserve(Coalescence_table.capacity() * 2);
    }
    Coalescence_table.push_back(event);
}

std::vector<coa_event_t> const &coa_table_c::get_coa_table()
{
    return Coalescence_table;
}

int &coa_table_c::get_left_brkpt_l(coa_event_t &event)
{
    return std::get<0>(event);
}

int const &coa_table_c::get_left_brkpt_l(coa_event_t const &event)
{
    return std::get<0>(event);
}

int &coa_table_c::get_right_brkpt_r(coa_event_t &event)
{
    return std::get<1>(event);
}

int const &coa_table_c::get_right_brkpt_r(coa_event_t const &event)
{
    return std::get<1>(event);
}

int const &coa_table_c::get_num_node_u(coa_event_t const &event)
{
    return std::get<2>(event);
}

std::vector<int> const &coa_table_c::get_childs_node_c(coa_event_t const &event)
{
    return std::get<3>(event);
}

double const &coa_table_c::get_time_t(coa_event_t const &event)
{
    return std::get<4>(event);
}

int const &coa_table_c::get_gen(coa_event_t const &event)
{
    return std::get<5>(event);
}

coa_event_t const &coa_table_c::operator[](int index) const
{
    return Coalescence_table[index];
}

//Some Accessor

std::vector<coa_event_t>::iterator coa_table_c::begin()
{
    return Coalescence_table.begin();
}

std::vector<coa_event_t>::iterator coa_table_c::end()
{
    return Coalescence_table.end();
}

std::vector<coa_event_t>::const_iterator coa_table_c::begin() const
{
    return Coalescence_table.cbegin();
}

std::vector<coa_event_t>::const_iterator coa_table_c::end() const
{
    return Coalescence_table.cend();
}

std::size_t coa_table_c::size()
{
    return Coalescence_table.size();
}

void coa_table_c::sort_by_num_u(std::size_t index_begin_zone_at_group)
{
    //sort by left_brkpt_l, right_brkpt_l and node_parent
    // lambda function stored in a variable, the variable is a function !
    auto p_l_r_sort_func = [](coa_event_t &event1, coa_event_t &event2) {
        if (get_num_node_u(event1) < get_num_node_u(event2)) //Leaft_brkpt
        {
            return true;
        }
        else if (get_num_node_u(event1) == get_num_node_u(event2))
        {
            if (get_left_brkpt_l(event1) < get_left_brkpt_l(event2))
            {
                return true;
            }
            else if (get_left_brkpt_l(event1) == get_left_brkpt_l(event2))
            {
                return get_right_brkpt_r(event1) < get_right_brkpt_r(event2);
            }
        }
        return false;
    };

    sort(Coalescence_table.begin() + index_begin_zone_at_group, Coalescence_table.end(), p_l_r_sort_func);
}

//Use after a multi-coa event
//Work if gen/gen algo (^^)
//WARNING : Use after a potential multicoa event, node num need to be num and num+1
void coa_table_c::group_multi_coa(std::size_t index_begin_zone_at_group)
{
    //Special case du to reshape_segs in multi coa process.
    //p parent, c child => need to cut parent before other algo
    // |---|---| c
    // |-------| p

    //sort by node_parent (parent first), left_brkpt_l, right_brkpt_l
    // lambda function stored in a variable, the variable is a function !
    auto p_l_r_sort_func = [](coa_event_t &event1, coa_event_t &event2) {
        if (get_num_node_u(event1) > get_num_node_u(event2)) //Leaft_brkpt
        {
            return true;
        }
        else if (get_num_node_u(event1) == get_num_node_u(event2))
        {
            if (get_left_brkpt_l(event1) > get_left_brkpt_l(event2))
            {
                return true;
            }
            else if (get_left_brkpt_l(event1) == get_left_brkpt_l(event2))
            {
                return get_right_brkpt_r(event1) > get_right_brkpt_r(event2);
            }
        }
        return false;
    };

    std::size_t parent_index = index_begin_zone_at_group;
    std::size_t child_zone = index_begin_zone_at_group;

    while (parent_index < Coalescence_table.size() - 1)
    {
        std::size_t end = Coalescence_table.size();
        sort(Coalescence_table.begin() + index_begin_zone_at_group, Coalescence_table.end(), p_l_r_sort_func);

        while ((get_num_node_u(Coalescence_table[parent_index]) == get_num_node_u(Coalescence_table[child_zone])) && (child_zone < end))
        {
            ++child_zone;
        }

        while (parent_index < child_zone)
        {
            std::size_t child_index = child_zone;
            auto *parent = &Coalescence_table[parent_index];
            while (child_index < end)
            {
                auto *child = &Coalescence_table[child_index];
                std::size_t chi = 0;
                while (chi < get_childs_node_c(*parent).size())
                {
                    int cut = false;
                    if (get_childs_node_c(*parent).at(chi) == get_num_node_u(*child))
                    {
                        //Lp left parent, Rp right parent
                        //Lc left child, Rc right child
                        // |-------| p  OR // |-------| p   OR //    |-----| p
                        // |---| c         //   |---| c        //  |----| c
                        //parent need to be cut between Rp and Rc, the original will be [Lp : Rc] and the copy will be [Lc : Rp]
                        if ((get_right_brkpt_r(*child) > get_left_brkpt_l(*parent)) && (get_right_brkpt_r(*child) < get_right_brkpt_r(*parent)))
                        {
                            //split coa_event in two : one zone overlapping and the rest
                            //push back the copy
                            set_coa_event(*parent);
                            //WARNING potentialy bad ref du to resize of vector
                            child = &Coalescence_table[child_index];
                            parent = &Coalescence_table[parent_index];
                            //Original
                            get_right_brkpt_r(*parent) = get_right_brkpt_r(*child);
                            //Copy
                            get_left_brkpt_l(Coalescence_table[Coalescence_table.size() - 1]) = get_right_brkpt_r(*child);
                            cut = true;
                        }
                        // |-------| p
                        //     |---| c
                        if ((get_left_brkpt_l(*child) > get_left_brkpt_l(*parent)) && (get_right_brkpt_r(*child) == get_right_brkpt_r(*parent)))
                        {
                            //split coa_event in two : one zone overlapping and the rest
                            //push back the copy
                            set_coa_event(*parent);
                            //WARNING potentialy bad ref du to resize of vector
                            child = &Coalescence_table[child_index];
                            parent = &Coalescence_table[parent_index];
                            //Original
                            get_right_brkpt_r(*parent) = get_left_brkpt_l(*child);
                            //Copy
                            get_left_brkpt_l(Coalescence_table[Coalescence_table.size() - 1]) = get_left_brkpt_l(*child);
                            cut = true;
                        }

                        // |---| p                        //   |---| p                       // |---| p
                        // |-------| c                    // |-------| c                     //   |-----| c
                        //child need to be cut between Rp and Rc, the original will be [Lc : Rp] and the copy will be [Lp : Rc]
                        if ((get_right_brkpt_r(*child) > get_right_brkpt_r(*parent)) && (get_left_brkpt_l(*child) < get_right_brkpt_r(*parent)))
                        {
                            set_coa_event(*child);

                            child = &Coalescence_table[child_index];
                            parent = &Coalescence_table[parent_index];

                            get_right_brkpt_r(*child) = get_right_brkpt_r(*parent);
                            get_left_brkpt_l(Coalescence_table[Coalescence_table.size() - 1]) = get_right_brkpt_r(*parent);
                            cut = true;
                        }
                        //     |---| p
                        // |-------| c
                        //child need to be cut between Lc and Lp
                        if ((get_left_brkpt_l(*child) < get_left_brkpt_l(*parent)) && (get_right_brkpt_r(*child) == get_right_brkpt_r(*parent)))
                        {
                            set_coa_event(*child);

                            child = &Coalescence_table[child_index];
                            parent = &Coalescence_table[parent_index];

                            get_right_brkpt_r(*child) = get_left_brkpt_l(*parent);
                            get_left_brkpt_l(Coalescence_table[Coalescence_table.size() - 1]) = get_left_brkpt_l(*parent);
                            cut = true;
                        }
                    }
                    if (!cut)
                    {
                        ++chi;
                    }
                }
                ++child_index;
            }
            ++parent_index;
        }
    }
    //sort by left_brkpt_l, right_brkpt_l and node_parent
    // lambda function stored in a variable, the variable is a function !
    auto l_r_sort_func = [](coa_event_t &event1, coa_event_t &event2) {
        if (get_left_brkpt_l(event1) < get_left_brkpt_l(event2)) //Leaft_brkpt
        {
            return true;
        }
        else if (get_left_brkpt_l(event1) == get_left_brkpt_l(event2))
        {
            if (get_right_brkpt_r(event1) > get_right_brkpt_r(event2))
            {
                return true;
            }
            else if (get_right_brkpt_r(event1) == get_right_brkpt_r(event2))
            {
                return get_num_node_u(event1) < get_num_node_u(event2);
            }
        }
        return false;
    };

    std::size_t i = Coalescence_table.size() - 1;

    //Safety first
    sort(Coalescence_table.begin() + index_begin_zone_at_group, Coalescence_table.end(), l_r_sort_func);

    while (i > index_begin_zone_at_group)
    {
        //in certain case some coa event of the multi event can be separe by coa event in the homolog chr
        for (int j = i - 1; j >= static_cast<int>(index_begin_zone_at_group); --j)
        {
            //overlap here
            if (get_left_brkpt_l(Coalescence_table[i]) == get_left_brkpt_l(Coalescence_table[j]))
            {
                if (get_right_brkpt_r(Coalescence_table[i]) == get_right_brkpt_r(Coalescence_table[j]))
                {
                    bool intermediate_node = false;
                    //Rassemble the two in the same event
                    //by adding the childs of the intermediate event j to the kept event i
                    //and remove event j
                    for (std::size_t k = 0; k < get_childs_node_c(Coalescence_table[i]).size(); ++k)
                    {
                        //Check if one of the child nodes is the intermediate node
                        //Delete the intermediate node
                        if (get_childs_node_c(Coalescence_table[i])[k] == get_num_node_u(Coalescence_table[j]))
                        {
                            std::get<3>(Coalescence_table[i]).erase(get_childs_node_c(Coalescence_table[i]).begin() + k);
                            std::get<3>(Coalescence_table[i]).insert(get_childs_node_c(Coalescence_table[i]).end(), get_childs_node_c(Coalescence_table[j]).begin(), get_childs_node_c(Coalescence_table[j]).end());
                            intermediate_node = true;
                            break;
                        }
                    }
                    //Delete the intermediate event
                    if (intermediate_node)
                    {
                        Coalescence_table.erase(Coalescence_table.begin() + j);
                    }
                }
            }
        }
        --i;
    }
}

/**************************************************************/
/*                 Tree Generator                             */
/**************************************************************/

tree_gen_c::tree_gen_c(coa_table_c const &coa_table, int next_node_ident, int n_sample_size)
{
    I_insert_order_recomb = sort_table_I(coa_table);
    R_remove_order_recomb = sort_table_R(coa_table);
    Max_u = next_node_ident;
    Leaf_number = n_sample_size;

    // create first empty Ancestry_tree
    // std::vector<std::tuple<childs_node, time_t, gen, nbr_leaf_for_node>>;
    Ancestry_tree = top_down_tree_t(Max_u, std::make_tuple(std::vector<int>(), 0, 0));
}

top_down_tree_t const &tree_gen_c::get_tree()
{
    return Ancestry_tree;
}

//TODO : Peut etre à changer
//The last non null node is the MRCA
int find_MRCA(top_down_tree_t const &ancestry_tree)
{
    int number_u_MRCA = ancestry_tree.size() - 1;
    while ((std::get<0>(ancestry_tree[number_u_MRCA])).empty())
    {
        --number_u_MRCA;
    }
    return number_u_MRCA;
}

//Same than previous function but in a different order due to the need to gennerate one tree per tree with a internal Ancestry_tree
std::tuple<int, int, double> tree_gen_c::set_next_tree()
{
    std::tuple<int, int, double> begin_end_MRCAtime_non_recomb_chunk;
    double time_MRCA;

    //T4. Remove recomb
    //Not use at first call because Begin_sequence != get_right_brkpt_r(R_remove_order_recomb[index_table_R])
    // check Kelleher Supp Mat 2  and exemple of I & R tables  in XXX to understand this part.
    Begin_sequence = get_left_brkpt_l(I_insert_order_recomb[index_table_I]);
    while (Begin_sequence == get_right_brkpt_r(R_remove_order_recomb[index_table_R])) // never happened for first tree
    {
        Ancestry_tree[get_num_node_u(R_remove_order_recomb[index_table_R])] = std::make_tuple(std::vector<int>(), 0, 0);
        index_table_R += 1;
    }

    //T2. Insert recomb
    // while on the same chunk
    while ((Begin_sequence == get_left_brkpt_l(I_insert_order_recomb[index_table_I])) && (index_table_I < I_insert_order_recomb.size()))
    {
        coa_event_t event = I_insert_order_recomb[index_table_I];
        Ancestry_tree[get_num_node_u(event)] = std::make_tuple(get_childs_node_c(event), get_time_t(event), get_gen(event));
        index_table_I += 1;
    }

    //T3bis. Calculate MRCA and begin_end sequence
    // Find the MRCA
    MRCA_node_nbr = find_MRCA(Ancestry_tree);
    //Find MRCA time
    time_MRCA = std::get<1>(Ancestry_tree[MRCA_node_nbr]); // continuous time
    if (time_MRCA == 0)
    {
        time_MRCA = std::get<2>(Ancestry_tree[MRCA_node_nbr]); // in generations
    }
    //By construction R_remove_order_recomb[index_table_R] is the event that have the right_brkpt of the current chunk
    begin_end_MRCAtime_non_recomb_chunk = std::make_tuple(Begin_sequence, get_right_brkpt_r(R_remove_order_recomb[index_table_R]), time_MRCA);

    return begin_end_MRCAtime_non_recomb_chunk;
}

//Function for sort table I in purpose to generate tree (see Kelleher 2016, p.13)
coa_table_c tree_gen_c::sort_table_I(coa_table_c i_insert_order_recomb)
{
    //Use anonymous function as classifier to sort function
    std::sort(i_insert_order_recomb.begin(), i_insert_order_recomb.end(),
              [](coa_event_t &event1, coa_event_t &event2) {
                  if (coa_table_c::get_left_brkpt_l(event1) < coa_table_c::get_left_brkpt_l(event2)) //Left_brkpt
                  {
                      return true;
                  }
                  else if (coa_table_c::get_left_brkpt_l(event1) == coa_table_c::get_left_brkpt_l(event2)) //Left_brkpt and time
                  {
                      if (coa_table_c::get_gen(event1) != -1)
                      {
                          return coa_table_c::get_gen(event1) < coa_table_c::get_gen(event2);
                      }
                      else
                      {
                          return coa_table_c::get_time_t(event1) < coa_table_c::get_time_t(event2);
                      }
                  }
                  else
                  {
                      return false;
                  }
              });
    return i_insert_order_recomb;
}

// Function for sort table O in purpose to generate tree (see Kelleher 2016, p.13)
coa_table_c tree_gen_c::sort_table_R(coa_table_c r_remove_order_recomb)
{
    //Use anonymous function as classifier to sort function
    std::sort(r_remove_order_recomb.begin(), r_remove_order_recomb.end(),
              [](coa_event_t &event1, coa_event_t &event2) {
                  if (coa_table_c::get_right_brkpt_r(event1) < coa_table_c::get_right_brkpt_r(event2)) //Right_brkpt
                  {
                      return true;
                  }
                  else if (coa_table_c::get_right_brkpt_r(event1) == coa_table_c::get_right_brkpt_r(event2)) //Right_brkpt and time
                  {
                      if (coa_table_c::get_gen(event1) != -1)
                      {
                          return coa_table_c::get_gen(event1) > coa_table_c::get_gen(event2);
                      }
                      else
                      {
                          return coa_table_c::get_time_t(event1) > coa_table_c::get_time_t(event2);
                      }
                  }
                  else
                  {
                      return false;
                  }
              });
    return r_remove_order_recomb;
}