#pragma once

#include <filesystem>

#include "summary_stat.hpp"
#include "mutation.hpp"

double const PI = 3.141592653589793;

template <typename values>
void print_output(std::string path_to_file, std::vector<values> vec_value, std::string open_file_mode);

void output_stat_files(simu_param_c const &simu_p, info_collector_c const &info, output_stat_c const &output, int rep, int ploidy);

void genepop_output(std::string const &path_to_file, std::vector<std::array<int, 2>> const &coord_vec, std::vector<std::vector<int>> const &ancestry_seq, std::vector<std::vector<std::vector<std::pair<int, int>>>> const &sample_mutated_state_chr, std::string const &mut_mod, int ploidy);
void vcf_output(std::string const &path_to_file, std::vector<std::array<int, 2>> const &coord_vec, std::vector<std::vector<int>> const &ancestry_seq, std::vector<std::vector<std::vector<std::pair<int, int>>>> const &sample_mutated_state_chr, int ploidy);
void fasta_output(std::string const &path_to_file, std::string const &generic_data_filename, std::vector<std::array<int, 2>> const &coord_vec, std::vector<std::vector<int>> const &ancestry_seq, std::vector<std::vector<std::vector<std::pair<int, int>>>> const &sample_mutated_state_chr, int ploidy, int rep);
void phylip_output(std::string const &path_to_file, std::string const &generic_data_filename, std::vector<std::array<int, 2>> const &coord_vec, std::vector<std::vector<int>> const &ancestry_seq, std::vector<std::vector<std::pair<int, int>>> const &sample_mutated_state_ind, int ploidy, int chrom);
void seq_char_output(std::string const &path_to_file, std::vector<std::array<int, 2>> const &coord_vec, std::vector<std::vector<std::vector<std::pair<int, int>>>> const &sample_mutated_state_chr, int ploidy);
void coord_output(std::string const &path_to_file, std::vector<std::array<int, 2>> const &coord_vec);
void write_beforerun_param_settings_summary(std::string const &path_to_file, std::string const &version, simu_param_c const &simu_param, samp_param_c const &samp_param, demo_param_c const &demo_param, muta_param_c const &muta_param, recomb_param_c const &recomb_param);

void write_afterrun_param_settings_summary(std::string const &path_to_file, simu_param_c const &simu_param, samp_param_c const &samp_param, demo_param_c const &demo_param);

#include "output.tpp"
