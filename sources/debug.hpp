#include <vector>
#include <string>

static void check_tree(std::vector<std::tuple<std::vector<int>, double, int>> const &gene_tree, int MRCA_index, int nbr_of_sample, int rep)
{
    std::vector<int> node_file;
    node_file.reserve(MRCA_index);

    node_file.push_back(MRCA_index);

    for (std::size_t index = 0; index < node_file.size(); ++index)
    {
        std::vector<int> children_s_node = std::get<0>(gene_tree[node_file[index]]);

        if (!children_s_node.empty()) // test if it is a leaf
        {
            // for all children, except last one (for which the map is moved to below), copy the mut_state info and apply mut
            //TODO change iter to child
            for (int child : children_s_node)
            {
                node_file.push_back(child);
            }
            node_file[index] = -1;
        }
    }
    auto verif = 0;
    for (int value : node_file)
    {
        if (value != -1)
        {
            ++verif;
        }
    }

    if (verif != nbr_of_sample)
    {
        std::string str = "( In check_tree() : Tree not valid at rep " + std::to_string(rep) + ". Contact the developpers. I exit. )";
        throw std::logic_error(str);
    }
}
