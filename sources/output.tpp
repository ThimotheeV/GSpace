#include <vector>
#include <string>
#include <fstream>

template <typename values>
void print_output(std::string path_to_file, std::vector<values> vec_value, std::string open_file_mode)
{
    std::ofstream file;
    if (open_file_mode == "over")
    {
        file.open(path_to_file, std::ios::out);
    } else {
        file.open(path_to_file, std::ios::app);
    }
    auto vec_value_itr = vec_value.begin();
    if (file.is_open())
    {
        for (std::size_t i = 0; i < vec_value.size() - 1; ++i)
        {
            file << *vec_value_itr << "\t";
            ++vec_value_itr;
        }
        file << *vec_value_itr << "\n";
        file.close();
    }
    else
    {
        throw std::invalid_argument("( Unable to open " + path_to_file + ". I exit. )");
    }
}
