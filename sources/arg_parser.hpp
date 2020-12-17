#pragma once

#include <string>
#include <vector>
#include <sstream>

std::string remove_underscores(std::string str);
std::string remove_spaces_tab_underscores(std::string str);
std::string remove_spaces_tab_in_range(std::string str, int pos_beg, int pos_end);
std::string remove_comma_semicolons_in_range(std::string str, int pos_beg, int pos_end);
std::string replace_spaces_tab_in_range_by_comma(std::string str, int pos_beg, int pos_end);
std::string replace_spaces_tab_underscores_by_comma(std::string str);
std::string str_tolower(std::string str);
std::vector<std::string> slice_by_char(std::string const &str, char sep);
std::vector<std::string> slice_by_comma_semicolon(std::string const &str);
std::vector<std::string> slice_unix_windows_file_by_line(std::string const &str);
std::vector<std::vector<std::string>> slice_str_vec_by_string(std::vector<std::string> const &vec_str, std::string const &sep);

std::vector<int> int_vector_parse_by_comma(std::string const &str);
std::vector<double> double_vector_parse_by_comma(std::string const &str);
std::vector<std::string> str_vector_parse_by_comma(std::string const &str);
std::vector<std::string> str_vector_parse_by_semicolon(std::string const &str);
std::vector<std::string> str_vector_parse_by_comma_semicolon(std::string const &str);

//#include <experimental/charconv>

// template <typename typen>
// std::vector<typen> vector_parse(std::string const &str, char c)
// {
//     std::vector<std::string> temp_result = slice_by_char(str, c);

//     std::vector<typen> result(temp_result.size());
//     auto result_itr = result.begin();
//     for (auto const &value : temp_result)
//     {
//         std::from_chars(value.data(), value.data()+value.size(), *result_itr);
//         ++result_itr;
//     }
//     return result;
// }

//TODO : temp place
bool convert_str_bool(std::string const &str);
bool str_has_only_digits(std::string const &str);
