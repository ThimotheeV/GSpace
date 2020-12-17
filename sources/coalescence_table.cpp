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
    if (Coalescence_table.capacity() < Coalescence_table.size() + 1)
    {
        Coalescence_table.reserve(Coalescence_table.capacity() * 2);
    }
    Coalescence_table.push_back(event);
}

void coa_table_c::set_coa_event(coa_event_t event)
{
    if (Coalescence_table.capacity() < Coalescence_table.size() + 1)
    {
        Coalescence_table.reserve(Coalescence_table.capacity() * 2);
    }
    Coalescence_table.push_back(std::move(event));
}

std::vector<coa_event_t> const &coa_table_c::get_coa_table()
{
    return Coalescence_table;
}

int const &coa_table_c::get_left_brkpt_l(coa_event_t const &event)
{
    return std::get<0>(event);
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

std::size_t coa_table_c::size()
{
    return Coalescence_table.size();
}

//Use after a multi-coa event
// TODO WARNING : probleme si multicoalescence sur homologues sur même zone génomique?
void coa_table_c::group_multi_coa(std::size_t index_begin_zone_at_group)
{
    if ((!Coalescence_table.empty()) && (get_gen(Coalescence_table[0]) == -1))
    {
        throw std::logic_error("( Use this function in a gen_by_gen coa_table. Contact the developpers. I exit. )"); // TODO eclaircir ce mesage obscur
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
            if (get_right_brkpt_r(event1) < get_right_brkpt_r(event2))
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

    //begin sorting
    //just handle the new coalescence entry not the old one (already handle)
    sort(Coalescence_table.begin() + index_begin_zone_at_group, Coalescence_table.end(), l_r_sort_func);

    std::size_t i = index_begin_zone_at_group;

    //Separate overlapping coalescence and the rest
    while (i < Coalescence_table.size() - 1)
    {
        //overlap here
        if (get_left_brkpt_l(Coalescence_table[i]) == get_left_brkpt_l(Coalescence_table[i + 1]))
        {
            //Need to resize it to make obvious overlap
            if (get_right_brkpt_r(Coalescence_table[i]) != get_right_brkpt_r(Coalescence_table[i + 1]))
            {
                //split coa_event in two : one zone overlapping and the rest
                set_coa_event(Coalescence_table[i + 1]);
                //TODO : Changer les get en accesseurs
                std::get<1>(Coalescence_table[i + 1]) = get_right_brkpt_r(Coalescence_table[i]);
                std::get<0>(Coalescence_table[Coalescence_table.size() - 1]) = get_right_brkpt_r(Coalescence_table[i]);
                //sort and continue
                sort(Coalescence_table.begin() + index_begin_zone_at_group, Coalescence_table.end(), l_r_sort_func);
                //begin again
                i = index_begin_zone_at_group - 1;
            }
        }
        else
        {
            if (get_right_brkpt_r(Coalescence_table[i]) == get_right_brkpt_r(Coalescence_table[i + 1]))
            {
                //split coa_event in two : one zone overlapping and the rest
                set_coa_event(Coalescence_table[i]);
                std::get<1>(Coalescence_table[i]) = get_left_brkpt_l(Coalescence_table[i + 1]);
                std::get<0>(Coalescence_table[Coalescence_table.size() - 1]) = get_left_brkpt_l(Coalescence_table[i + 1]);
                //sort and continue
                sort(Coalescence_table.begin() + index_begin_zone_at_group, Coalescence_table.end(), l_r_sort_func);
                //begin again
                i = index_begin_zone_at_group - 1;
            }
        }
        ++i;
    }
    i = Coalescence_table.size() - 1;

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
    Ancestry_tree = top_down_tree_t(Max_u, std::make_tuple(std::vector<int>(), 0, 0, 1));
}

top_down_tree_t const &tree_gen_c::get_tree()
{
    return Ancestry_tree;
}

int tree_gen_c::nbr_leaf_for_node(std::vector<int> const &node_childs)
{
    int sum_leaf = 0;
    for (auto const &child : node_childs)
    {
        sum_leaf += std::get<3>(Ancestry_tree[child]);
    }

    return sum_leaf;
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
        Ancestry_tree[get_num_node_u(R_remove_order_recomb[index_table_R])] = std::make_tuple(std::vector<int>(), 0, 0, 0);
        index_table_R += 1;
    }

    //T2. Insert recomb
    // while on the same chunk
    while ((Begin_sequence == get_left_brkpt_l(I_insert_order_recomb[index_table_I])) && (index_table_I < I_insert_order_recomb.size()))
    {
        coa_event_t event = I_insert_order_recomb[index_table_I];
        Ancestry_tree[get_num_node_u(event)] = std::make_tuple(get_childs_node_c(event), get_time_t(event), get_gen(event), nbr_leaf_for_node(get_childs_node_c(event)));
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
