#include <sstream>
#include <fstream>
#include <iterator>
#include <functional>
#include <iostream>

#include "settings.hpp"
#include "mutation.hpp"
#include "output.hpp"
#include "simulator.hpp"

extern double const DBL_EPSILON = 1e-6;

//TODO : If file doesn't exist, don't throw
std::string const read_file(std::string const &filename)
{
    //Put filename in a stream
    //ifstream constructor
    std::ifstream file_str(filename);
    if (!file_str)
    {
        throw std::logic_error("( Could not open file : " + filename + ". I exit. )");
    }
    else
    {
        //Test if the stream is good
        if (file_str.bad())
        {
            throw std::logic_error("( Could not open file : " + filename + ". I exit. )");
        }
        //Don't skip space
        file_str >> std::noskipws;
        //Convert the stream in a string
        const std::string &fstr{std::istream_iterator<char>{file_str}, {}};
        return fstr;
    }
}

/**
 * @brief Parser file to param 
 * 
 */

std::string read_write_cmdline(int argc, char **argv)
{
    std::string result;
    if (argc > 1)
    {
        for (int arg = 1; arg < argc - 1; ++arg)
        {
            result += argv[arg];
            result += "\n";
        }
        result += argv[argc - 1];
    }
    //Put filename in a stream
    std::ofstream file_str("cmdline_settings.txt");
    if (!file_str)
    {
        throw std::logic_error("( Could not open file 'cmdline_settings.txt'. I exit. )");
    }
    else
    {
        file_str << argv[0] << "\n";
        file_str << result;
        file_str.close();
    }

    return result;
}

void parser_str(std::string const &file_str)
{
    const auto &str_vec = slice_unix_windows_file_by_line(file_str);

    //Call uniq instance of rand_gen, simu_param_c, demo_param_c, demo_param_c, recomb_param_c, muta_param_c
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();

    for (auto const &line : str_vec)
    {
        //Clean and parse key=value in a vector[key, value]
        const auto &clean_line = remove_spaces_tab_in_range(line, 0, static_cast<int>(line.size()));
        auto line_pair_key_value = slice_by_char(clean_line, '=');
        const auto &param_name = remove_underscores(str_tolower(line_pair_key_value[0]));

        if (line_pair_key_value[1].empty())
        {
            throw std::invalid_argument("( Setting " + param_name + " seems to have no value. Check settings. I exit. )");
        }

        ////SIMULATION PARAMETERS
        if (param_name == "settingfilename" || param_name == "settingsfilename" || param_name == "settingfile" || param_name == "settingsfile")
        {
            simu_param.Setting_filename = line_pair_key_value[1];
            continue;
        }
        if (param_name == "randomseeds")
        {
            int const seed = std::stoi(line_pair_key_value[1]);
            if (seed < 1 || seed > std::numeric_limits<int>::max())
            {
                throw std::invalid_argument("( Random_Seeds should always be between 1 and max(int)=" +
                                            std::to_string(std::numeric_limits<int>::max()) + ". Check settings. I exit. )");
            }

            simu_param.seed = seed;
            rand_gen.put_seed(seed);
            continue;
        }
        if (param_name == "runnumber" || param_name == "runnbr")
        {
            simu_param.Repetition_nbr = std::stoi(line_pair_key_value[1]);
            if (simu_param.Repetition_nbr < 1 || simu_param.Repetition_nbr > std::numeric_limits<int>::max())
            {
                throw std::invalid_argument("( Run_Number should always be between 1 and max(int)=" +
                                            std::to_string(std::numeric_limits<int>::max()) + ". Check settings. I exit. )");
            }

            continue;
        }
        if (param_name == "approximatetime")
        {
            try
            {
                simu_param.Continuous_time_approxim = convert_str_bool(line_pair_key_value[1]);
            }
            catch (const std::logic_error &e)
            {
                //ecritre un truc ici
            }
            continue;
        }
        if (param_name == "pause")
        {
            const std::string value = str_tolower(line_pair_key_value[1]);

            if (value == "onerror")
            {
                simu_param.Wait_for_cin_input = true;
                simu_param.Wait_for_final_cin_input = true;
            }
            else if (value == "final")
            {
                simu_param.Wait_for_cin_input = false;
                simu_param.Wait_for_final_cin_input = true;
            }
            else if (value == "never")
            {
                simu_param.Wait_for_cin_input = false;
                simu_param.Wait_for_final_cin_input = false;
            }
            else if (value == "default")
            {
                simu_param.Wait_for_cin_input = false;
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
                simu_param.Wait_for_final_cin_input = true;
#else
                simu_param.Wait_for_final_cin_input = false;
#endif
            }
            else
            {
                throw std::invalid_argument("( Unknown arguments for Pause= " + value + ". Only OnError, Final, Never (Default=Final for windows, Never for Linux) are currently implemented. Check settings. I exit. )");
            }
            continue;
        }

        ////OUTPUT FILE FORMAT OPTIONS
        if (param_name == "outputdir")
        {
            simu_param.Output_dir = line_pair_key_value[1];
            continue;
        }
        if (param_name == "datafilename" || param_name == "outputfilename")
        {
            simu_param.Generic_data_filename = line_pair_key_value[1];
            continue;
        }
        if (param_name == "datafileextension" || param_name == "outputfileextension")
        {
            simu_param.Data_file_extension = line_pair_key_value[1];
            continue;
        }
        if (param_name == "genepop")
        {
            simu_param.Genepop_output = convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (param_name == "genepopindfile")
        {
            simu_param.Genepop_ind_file = convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (param_name == "genepopgroupallsamples" || param_name == "genepopgroupallsample")
        {
            simu_param.Genepop_group_all_samples = convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (param_name == "genepopnocoord")
        {
            simu_param.Genepop_no_coord = convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (param_name == "vcf")
        {
            simu_param.VCF_output = convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (param_name == "coordfile" || param_name == "coordinatefile")
        {
            simu_param.Coordinates_output = convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (param_name == "seqcharfile" || param_name == "sequencecharacteristicsfile" || param_name == "sequencecharacteristicfile")
        {
            simu_param.Seq_char_output = convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (param_name == "fasta")
        {
            simu_param.Fasta_output = convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (param_name == "fastasinglelineseq")
        {
            simu_param.Fasta_single_line_seq = convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (param_name == "phylip")
        {
            simu_param.Phylip_output = convert_str_bool(line_pair_key_value[1]);
            if (simu_param.Phylip_output)
                simu_param.Seq_char_output = true;
            continue;
        }

        ////MARKERS PARAMETERS
        if ((param_name == "ploidy") || (param_name == "ploïdy"))
        {
            auto const temp = str_tolower(line_pair_key_value[1]);

            if (temp == "diploid" || temp == "diploïd")
            {
                samp_param.Ploidy = 2;
            }
            else
            {
                if (temp != "haploid" && temp != "haploïd")
                {
                    throw std::invalid_argument("( Ploidy argument should be Haploid or Diploid and not " + temp + ". Check settings. I exit. )");
                }
                samp_param.Ploidy = 1;
            }
            continue;
        }
        // Chromosome_Number
        if (param_name == "locusnumber" || param_name == "locusnbr" || param_name == "chromosomenumber" || param_name == "chromosomenbr")
        {
            samp_param.Chr_nbr = std::stoi(line_pair_key_value[1]);
            if (samp_param.Chr_nbr < 1 || samp_param.Chr_nbr > std::numeric_limits<int>::max())
            {
                throw std::invalid_argument("( Locus_Number or Chromosome_Number should always be between 1 and max(int)=" +
                                            std::to_string(std::numeric_limits<int>::max()) + ". Check settings. )");
            }
            continue;
        }
        //muta_param_c
        if (param_name == "mutationmodel")
        {
            auto const &temp = str_tolower(line_pair_key_value[1]);
            map_string_mod_mut_name_c m;
            const auto &mut_name = m.Map.find(temp);

            if (mut_name != m.Map.end())
            {
                muta_param.Mod_mut_name = mut_name->second;
            }
            else
            {
                throw std::logic_error("( Mutation_Model=" + temp + " specification is not valid. Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "mutationrate")
        {
            muta_param.Unscaled_mut_rate_mu = std::stod(line_pair_key_value[1]);
            if (muta_param.Unscaled_mut_rate_mu < 0.0 || muta_param.Unscaled_mut_rate_mu > 1.0)
            {
                throw std::invalid_argument("( The Mutation_Rate should always be between 0 and 1. Check settings. I exit. )");
            }
            continue;
        }
        //Seq of ATCG with comma separation
        if (param_name == "mrcasequence" || param_name == "mrcanucleotidicsequence" || param_name == "mrcanucleotidicseq")
        {
            auto const &temp = str_vector_parse_by_comma_semicolon(line_pair_key_value[1]);
            //Nbr of chromosome
            muta_param.Ancestry_seq.resize(temp.size());
            //iterator -> vector<int>
            auto ancestry_seq_itr = muta_param.Ancestry_seq.begin();

            for (auto const &seq : temp)
            {
                //Translate alphabetic nucl in numeric value
                map_string_nuc_name_c nuc_name;

                std::size_t index = seq.find_first_not_of("ATGCatgc0");
                if (index != std::string::npos)
                {
                    std::string outputStr(1, seq[index]);
                    throw std::invalid_argument("( Ancestry sequence " + seq + " contains (at least) one character " + outputStr + " different from A/T/G/C or 0. Check MRCA_sequence settings. I exit. )");
                }
                else
                {
                    ancestry_seq_itr->clear();
                    ancestry_seq_itr->reserve(seq.size());
                    for (auto const car : seq)
                    {
                        std::string carStr(1, car);
                        if (carStr == "0")
                        {
                            ancestry_seq_itr->push_back(0);
                        }
                        else
                        {
                            ancestry_seq_itr->push_back(nuc_name.Map.at(carStr));
                        }
                    }
                }
                ++ancestry_seq_itr;
            }

            // size of the MRCA sequence will be chekced later in check_param()
            continue;
        }
        //Sequence of Allelic value with comma separation and semicolon between chr
        if (param_name == "mrcaallelicstate" || param_name == "mrcaallelicstate" || param_name == "allelicstateMRCA" || param_name == "allelicstatesMRCA" || param_name == "mrcaallelicsequence" || param_name == "mrcaallelicseq")
        {
            auto const &temp = slice_by_char(line_pair_key_value[1], ';');
            //Nbr of chromosome
            muta_param.Ancestry_seq.resize(temp.size());
            //iterator -> vector<int>
            auto ancestry_seq_itr = muta_param.Ancestry_seq.begin();
            int chr = 1;
            for (auto const &seq : temp)
            {
                int allele = 1;
                auto const &allelic_seq = slice_by_char(seq, ',');
                ancestry_seq_itr->resize(allelic_seq.size());
                auto allelic_ancestry_seq_itr = ancestry_seq_itr->begin();
                for (auto const &state : allelic_seq)
                {
                    if (str_has_only_digits(state))
                    {
                        int intvalue = stoi(state);
                        if ((intvalue >= 0) && (intvalue <= 999))
                        {
                            *allelic_ancestry_seq_itr = intvalue;
                        }
                        else
                        {
                            throw std::invalid_argument("( Ancestry allelic state found to be " + std::to_string(intvalue) + " an integer larger than 999. Check MRCA_sequence settings. I exit. )");
                        }
                    }
                    else
                    {
                        throw std::invalid_argument("( In chromosome number " + std::to_string(chr) + " ancestry allelic state number " + std::to_string(allele) + " is not a single integer number (for allelic data) or a single string of A/T/G/C or 0 but a vector with , or ;. Check MRCA_sequence settings. I exit. )");
                    }
                    ++allelic_ancestry_seq_itr;
                    ++allele;
                }
                ++ancestry_seq_itr;
                ++chr;
            }
            // size of the MRCA sequence will be chekced later in check_param()
            continue;
        }
        if (param_name == "alleliclowerbound")
        {
            muta_param.K_min = std::stoi(line_pair_key_value[1]);
            if (muta_param.K_min < 1 || muta_param.K_min > 999)
            {
                throw std::invalid_argument("( Allelic_Lower_Bound should always be between 1 and 999. Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "allelicupperbound")
        {
            muta_param.K_max = std::stoi(line_pair_key_value[1]);
            if (muta_param.K_max < 1 || muta_param.K_max > 999)
            {
                throw std::invalid_argument("( Allelic_Upper_Bound should always be between 1 and 999. Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "pgsm")
        {
            muta_param.P_gsm = std::stod(line_pair_key_value[1]);
            if (muta_param.P_gsm < 0 || muta_param.P_gsm > 1)
            {
                throw std::invalid_argument("( P_GSM should always be between 0 and 1. Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "transitiontransversionratio")
        {
            std::get<0>(muta_param.Ratio_transi_transver) = std::stod(line_pair_key_value[1]);
            if (std::get<0>(muta_param.Ratio_transi_transver) < 0)
            {
                throw std::invalid_argument("( Transition_Transversion_ratio should always be positive. Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "transition1transition2ratio")
        {
            std::get<1>(muta_param.Ratio_transi_transver) = std::stod(line_pair_key_value[1]);
            if (std::get<1>(muta_param.Ratio_transi_transver) < 0)
            {
                throw std::invalid_argument("( Transition1_Transition2_ratio should always be positive. Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "equilibriumfrequencies")
        {
            auto temp_vec = double_vector_parse_by_comma(line_pair_key_value[1]);
            muta_param.Equi_base_freq.A = temp_vec[0];
            muta_param.Equi_base_freq.G = temp_vec[1];
            muta_param.Equi_base_freq.C = temp_vec[2];
            muta_param.Equi_base_freq.T = temp_vec[3];
            if (muta_param.Equi_base_freq.A < 0 || muta_param.Equi_base_freq.A > 1 ||
                muta_param.Equi_base_freq.G < 0 || muta_param.Equi_base_freq.G > 1 ||
                muta_param.Equi_base_freq.C < 0 || muta_param.Equi_base_freq.C > 1 ||
                muta_param.Equi_base_freq.T < 0 || muta_param.Equi_base_freq.T > 1)
            {
                throw std::invalid_argument("( Equilibrium_Frequencies should always be between 0 and 1. Check settings. I exit. )");
            }
            continue;
        }

        ////SEQUENCE SPECIFIC SETTINGS
        // sequencelength (user version) ; sequencesize (IBDSim compatibility)
        if (param_name == "sequencelength" || param_name == "sequencesize")
        {
            samp_param.Sequence_length = std::stoi(line_pair_key_value[1]);
            if (samp_param.Sequence_length < 1 || samp_param.Sequence_length > std::numeric_limits<int>::max())
            {
                throw std::invalid_argument("( Sequence_Length/Sequence_Size should always be between 1 and max(int)=" +
                                            std::to_string(std::numeric_limits<int>::max()) + ". Check settings. I exit. )");
            }
            continue;
        }

        ////RECOMBINATION PARAMETERS
        //recomb_param_c
        if (param_name == "recombinationrate")
        {
            recomb_param.Unscaled_recomb_rate = std::stod(line_pair_key_value[1]);
            if (recomb_param.Unscaled_recomb_rate < 0.0 || recomb_param.Unscaled_recomb_rate > 1.0)
            {
                throw std::invalid_argument("( Recombination_Rate should always be between 0 and 1. Check settings. I exit. )");
            }
            continue;
        }

        ////DEMOGRAPHIC SETTINGS
        ///LATTICE
        // User-defined matrix of subpop sizes /  densities
        if (param_name == "nodesizematrix" || param_name == "subpopsizematrix" || param_name == "densitymatrix")
        {
            simu_param.Nodesize_matrix = convert_str_bool(line_pair_key_value[1]);
            continue;
        }
        if (param_name == "nodesizematrixfilename" ||
            param_name == "subpopsizematrixfilename" || param_name == "indpernodematrixfilename" || param_name == "densitymatrixfilename")
        {
            auto const temp = str_tolower(line_pair_key_value[1]);
            auto const temp_vec = slice_by_char(temp, '.');

            if (temp_vec[temp_vec.size() - 1] != "txt")
            {
                throw std::invalid_argument("( Node_Size_Matrix_Filename arguments does not seems to be a text file name (i.e. no .txt extension). Check settings. I exit. )");
            }
            else
            {
                simu_param.Nodesize_matrix_filename = temp;
            }
            continue;
        }
        if (param_name == "latticesizex")
        {
            demo_param.Lattice_size[0] = std::stoi(line_pair_key_value[1]) - 1;
            if (demo_param.Lattice_size[0] < 0 || demo_param.Lattice_size[0] > std::numeric_limits<int>::max())
            {
                throw std::invalid_argument("( Lattice_Size_X should always be between 1 and max(int)=" +
                                            std::to_string(std::numeric_limits<int>::max()) + ". Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "latticesizey")
        {
            demo_param.Lattice_size[1] = std::stoi(line_pair_key_value[1]) - 1;
            if (demo_param.Lattice_size[1] < 0 || demo_param.Lattice_size[1] > std::numeric_limits<int>::max())
            {
                throw std::invalid_argument("( Lattice_Size_Y should always be between 1 and max(int)=" +
                                            std::to_string(std::numeric_limits<int>::max()) + ". Check settings. I exit. )");
            }

            continue;
        }
        // indpernode (user version) ; indperpop (IBDSim compatibility)
        if (param_name == "indpernode" || param_name == "indperpop")
        {
            demo_param.Pop_size_per_node = std::stoi(line_pair_key_value[1]);
            if (demo_param.Pop_size_per_node < 1 || demo_param.Pop_size_per_node > std::numeric_limits<int>::max())
            {
                throw std::invalid_argument("( Ind_Per_Node/Ind_Per_Pop should always be between 1 and max(int)=" +
                                            std::to_string(std::numeric_limits<int>::max()) + ". Check settings. I exit. )");
            }
            continue;
        }

        ///SAMPLES
        //sample_disposition
        if (param_name == "samplecoordinatey" || param_name == "samplecoordy")
        {
            samp_param.Sample_coord_y = int_vector_parse_by_comma(line_pair_key_value[1]);
            continue;
        }
        if (param_name == "samplecoordinatex" || param_name == "samplecoordx")
        {
            samp_param.Sample_coord_x = int_vector_parse_by_comma(line_pair_key_value[1]);
            continue;
        }
        if (param_name == "samplesizex")
        {
            demo_param.Nbr_node_sampled_x = std::stoi(line_pair_key_value[1]);
            if (demo_param.Nbr_node_sampled_x < 1 || demo_param.Nbr_node_sampled_x > std::numeric_limits<int>::max())
            {
                throw std::invalid_argument("( Sample_Size_X should always be between 1 and max(int)=" +
                                            std::to_string(std::numeric_limits<int>::max()) + ". Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "samplesizey")
        {
            demo_param.Nbr_node_sampled_y = std::stoi(line_pair_key_value[1]);
            if (demo_param.Nbr_node_sampled_y < 1 || demo_param.Nbr_node_sampled_y > std::numeric_limits<int>::max())
            {
                throw std::invalid_argument("( Sample_Size_Y should always be between 1 and max(int)=" +
                                            std::to_string(std::numeric_limits<int>::max()) + ". Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "minsamplecoordinatex" || param_name == "minsamplecoordx")
        {
            demo_param.Min_sample_coord_x = std::stoi(line_pair_key_value[1]) - 1;
            if (demo_param.Min_sample_coord_x < 0 || demo_param.Min_sample_coord_x > std::numeric_limits<int>::max())
            {
                throw std::invalid_argument("( Min_Sample_Coordinate_X should always be between 1 and max(int)=" +
                                            std::to_string(std::numeric_limits<int>::max()) + ". Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "minsamplecoordinatey" || param_name == "minsamplecoordy")
        {
            demo_param.Min_sample_coord_y = std::stoi(line_pair_key_value[1]) - 1;
            if (demo_param.Min_sample_coord_y < 0 || demo_param.Min_sample_coord_y > std::numeric_limits<int>::max())
            {
                throw std::invalid_argument("( Min_Sample_Coordinate_Y should always be between 1 and max(int)=" +
                                            std::to_string(std::numeric_limits<int>::max()) + ". Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "voidsamplenodesx" || param_name == "voidsamplenodex" || param_name == "voidsamplepopx")
        {
            demo_param.Void_sample_nodes_X = std::stoi(line_pair_key_value[1]);

            if (demo_param.Void_sample_nodes_X < 0 || demo_param.Void_sample_nodes_X > std::numeric_limits<int>::max())
            {
                throw std::invalid_argument("( Void_Sample_Nodes_X should always be between 1 and max(int)=" +
                                            std::to_string(std::numeric_limits<int>::max()) + ". Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "voidsamplenodesy" || param_name == "voidsamplenodey" || param_name == "voidsamplepopy")
        {
            demo_param.Void_sample_nodes_Y = std::stoi(line_pair_key_value[1]);

            if (demo_param.Void_sample_nodes_Y < 0 || demo_param.Void_sample_nodes_Y > std::numeric_limits<int>::max())
            {
                throw std::invalid_argument("( Void_Sample_Nodes_Y should always be between 1 and max(int)=" +
                                            std::to_string(std::numeric_limits<int>::max()) + ". Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "voidsamplenodes" || param_name == "voidsamplenode" || param_name == "voidsamplepop")
        {
            demo_param.Void_sample_nodes_X = std::stoi(line_pair_key_value[1]);
            demo_param.Void_sample_nodes_Y = std::stoi(line_pair_key_value[1]);

            if (demo_param.Void_sample_nodes_X < 0 || demo_param.Void_sample_nodes_X > std::numeric_limits<int>::max())
            {
                throw std::invalid_argument("( Void_Sample_Nodes should always be between 1 and max(int)=" +
                                            std::to_string(std::numeric_limits<int>::max()) + ". Check settings. I exit. )");
            }
            continue;
        }
        // indpernodesampled (user version) ; indperpopsampled (IBDSim compatibility)
        if (param_name == "indpernodesampled" || param_name == "indperpopsampled")
        {
            auto const &temp = str_vector_parse_by_comma_semicolon(line_pair_key_value[1]);
            //size should be the Nbr of sampled nodes = demo_param.Nbr_node_sampled_x * demo_param.Nbr_node_sampled_y, checked in Check_settings()
            samp_param.Sample_size_per_node.reserve(temp.size());
            for (auto const &val_str : temp)
            {
                if (str_has_only_digits(val_str))
                {
                    int intvalue = stoi(val_str);
                    if ((intvalue >= 0) && (intvalue <= std::numeric_limits<int>::max()))
                    {
                        samp_param.Sample_size_per_node.push_back(intvalue);
                    }
                    else
                    {
                        throw std::invalid_argument("( At least one value in Ind_Per_Node_Sampled found to be " + std::to_string(intvalue) + " an integer larger than max(int)=" + std::to_string(std::numeric_limits<int>::max()) + " or <0. Check Ind_Per_Node_Sampled settings. I exit. )");
                    }
                }
                else
                {
                    throw std::invalid_argument("( At least one value in Ind_Per_Node_Sampled found to be " + val_str + " but only integers are allowed. Check Ind_Per_Node_Sampled settings. I exit. )");
                }
            }
            continue;
        }

        ///DISPERSAL
        // User-defined matrix of forward migration rates
        if (param_name == "migrationmatrix")
        {
            simu_param.Migration_matrix = convert_str_bool(line_pair_key_value[1]);
            if (simu_param.Migration_matrix_filename.empty() || simu_param.Migration_matrix_filename == "none")
            {
                simu_param.Migration_matrix_filename = "MigrationMatrix.txt";
            }
            continue;
        }
        if (param_name == "migrationmatrixfilename")
        {
            auto const temp = str_tolower(line_pair_key_value[1]);
            auto const temp_vec = slice_by_char(temp, '.');

            if (temp_vec[temp_vec.size() - 1] != "txt")
            {
                throw std::invalid_argument("( Migration_Matrix arguments does not seems to be a text file name (i.e. no '.txt' extension). Check settings. I exit. )");
            }
            else
            {
                simu_param.Migration_matrix_filename = temp;
            }
            continue;
        }
        if (param_name == "dispersaldistribution")
        {
            auto const temp = str_tolower(line_pair_key_value[1]);

            if ((temp == "u") || (temp == "uniform"))
            {
                demo_param.Dispersal_distrib = dispersal_distrib_enum::uniform;
            }
            else if ((temp == "n") || (temp == "gaussian"))
            {
                demo_param.Dispersal_distrib = dispersal_distrib_enum::gaussian;
            }
            else if ((temp == "g") || (temp == "geometric"))
            {
                demo_param.Dispersal_distrib = dispersal_distrib_enum::geometric;
            }
            else if ((temp == "p") || (temp == "pareto"))
            {
                demo_param.Dispersal_distrib = dispersal_distrib_enum::pareto;
            }
            else if ((temp == "s") || (temp == "sichel"))
            {
                demo_param.Dispersal_distrib = dispersal_distrib_enum::sichel;
            }
            else
            {
                throw std::invalid_argument("( Unknown arguments for Disperal_Distribution=" + temp + ". Only  Uniform/u, Gaussian/n, Geometric/g, Pareto/p, Sichel/s are currently implemented. Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "geometricshape")
        {
            demo_param.G_geo_param = std::stod(line_pair_key_value[1]);
            if (demo_param.G_geo_param < 0.0 || demo_param.G_geo_param > std::numeric_limits<double>::max())
            {
                throw std::invalid_argument("( Geometric_Shape should always be between 0.0 and max(double)=" +
                                            std::to_string(std::numeric_limits<double>::max()) + ". Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "paretoshape")
        {
            demo_param.N_pareto_param = std::stod(line_pair_key_value[1]);
            if (demo_param.N_pareto_param < 0.0 || demo_param.N_pareto_param > std::numeric_limits<double>::max())
            {
                throw std::invalid_argument("( Pareto_Shape should always be between 0.0 and max(double)=" +
                                            std::to_string(std::numeric_limits<double>::max()) + ". Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "sichelgamma")
        {
            demo_param.Gamma_sichel_param = std::stod(line_pair_key_value[1]);
            if (demo_param.Gamma_sichel_param > 0.0 || demo_param.Gamma_sichel_param > std::numeric_limits<double>::min())
            {
                throw std::invalid_argument("( Sichel_Gamma should always be a negative real between 0.0 and min(double)=" +
                                            std::to_string(std::numeric_limits<double>::min()) + ". Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "sichelxi")
        {
            demo_param.Xi_sichel_param = std::stod(line_pair_key_value[1]);
            if (demo_param.Xi_sichel_param < std::numeric_limits<double>::min() || demo_param.Xi_sichel_param > std::numeric_limits<double>::max())
            {
                throw std::invalid_argument("( Sichel_Xi should always be between min(double)=" +
                                            std::to_string(std::numeric_limits<double>::min()) + " and max(double)=" + std::to_string(std::numeric_limits<double>::max()) + ". Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "sichelomega")
        {
            demo_param.Omega_sichel_param = std::stod(line_pair_key_value[1]);
            if (demo_param.Xi_sichel_param < std::numeric_limits<double>::min() || demo_param.Xi_sichel_param > std::numeric_limits<double>::max())
            {
                throw std::invalid_argument("( Sichel_Omega should always be between min(double)=" +
                                            std::to_string(std::numeric_limits<double>::min()) + " and max(double)=" + std::to_string(std::numeric_limits<double>::max()) + ". Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "latticeboundaries" || param_name == "edgeeffects")
        {
            auto const temp = str_tolower(line_pair_key_value[1]);
            if (temp == "reflecting")
            {
                demo_param.Edge_effects = edge_effect_enum::reflecting;
            }
            else if (temp == "absorbing")
            {
                demo_param.Edge_effects = edge_effect_enum::absorbing;
            }
            else if (temp == "circular")
            {
                demo_param.Edge_effects = edge_effect_enum::circular;
            }
            else
            {
                throw std::invalid_argument("( Unknown arguments for Lattice_Boundaries/Edge_Effects=" + temp + ". Only Reflecting/u, Absorbing/g, Circular/p are currently implemented. Check settings. I exit. )");
            }
            continue;
        }
        if (param_name == "totalemigrationrate" || param_name == "emigrationrate" || param_name == "totalemmigrationrate")
        {
            demo_param.Proba_migr = std::stod(line_pair_key_value[1]);
            if (demo_param.Proba_migr < 0.0 || demo_param.Proba_migr > 1.0)
            {
                throw std::invalid_argument("( Total_emigration_Rate should always be between 0.0 and 1.0. Check settings. I exit. )");
            }
            continue;
        }
        // dispdistmax (user version) ; distmax (IBDSim compatibility)
        if (param_name == "dispdistmax" || param_name == "distmax")
        {
            auto temp = int_vector_parse_by_comma(line_pair_key_value[1]);
            demo_param.Disp_dist_max[0] = temp[0];
            if (demo_param.Disp_dist_max[0] < 0 || demo_param.Disp_dist_max[0] > std::numeric_limits<int>::max())
            {
                throw std::invalid_argument("( Disp_Dist_Max should always be between 0 and max(int)=" +
                                            std::to_string(std::numeric_limits<int>::max()) + ". Check settings. I exit. )");
            }
            if (temp.size() == 1)
            {
                demo_param.Disp_dist_max[1] = temp[0];
            }
            else
            {
                demo_param.Disp_dist_max[1] = temp[1];
                if (demo_param.Disp_dist_max[1] < 0 || demo_param.Disp_dist_max[1] > std::numeric_limits<int>::max())
                {
                    throw std::invalid_argument("( Disp_Dist_Max should always be between 0 and max(int)=" +
                                                std::to_string(std::numeric_limits<int>::max()) + ". Check settings. I exit. )");
                }
            }
            continue;
        }
        //VARIOUS COMPUTATION OPTION FOR DEV
        // postsimcomputation (user version) ; diagnostictables (IBDSim compatibility)
        if (param_name == "postsimcomputations" || param_name == "postsimcomputation" || param_name == "diagnostictables")
        {
            info_collect.Stats = true;
            auto const temp = str_vector_parse_by_comma_semicolon(str_tolower(line_pair_key_value[1]));
            for (auto const &temp2 : temp)
            {
                auto const &value = remove_underscores(temp2);
                if (value == "probid2loc" || value == "probid12loc" || value == "identityprobability2locus")
                {
                    info_collect.Prob_id_1_2_loc = true;
                    continue;
                }
                else if (value == "probid1loc" || value == "identityprobability1locus" || value == "identityprobability" || value == "qstat")
                {
                    //TODO : Les deux s'overlap dans pas mal de situation, ce serait dommage de ralentir les simu en calculant plusieurs fois les mêmes valeurs
                    info_collect.Prob_id_1_loc_Qr = true;
                    info_collect.Prob_id_1_loc_Qwi_wd_bd = true;
                    continue;
                }
                else if (value == "districoa" || value == "coalescencetimes" || value == "coatimes")
                {
                    info_collect.Coa_times = true;
                    continue;
                }
                else if (value == "effectivedispersal" || value == "empiricaldispersal" || value == "dispersal")
                {
                    info_collect.Effective_disp = true;
                    continue;
                }
                else if (value == "iterative" || value == "iterativestatistics")
                {
                    info_collect.Iterative_stats = true;
                    continue;
                }
                else if (value == "perloc" || value == "perlocstatistics")
                {
                    info_collect.Stats_per_loc_per_simu = true;
                    continue;
                }
                else if (value == "perchrom" || value == "perchromstatistics")
                {
                    info_collect.Stats_per_chr_per_simu = true;
                    continue;
                }
                else
                {
                    if (value.size() != 0)
                    {
                        throw std::invalid_argument("( Unknown arguments for Post_Sim_Computation/Diagnostic_Tables=" + value + ". Only probid2loc, qstat, districoa, effectivedispersal, iterative, perloc, perchrom are currently implemented. Check settings. I exit. )");
                    }
                }
            }
            continue;
        }
        if (param_name == "distclassnbr")
        {
            samp_param.Dist_class_nbr = std::stoi(line_pair_key_value[1]);
            continue;
        }
        throw std::invalid_argument("( Unknown keyworld : " + line + " from cmdline_settings.txt or GSpaceSettings.txt. I exit. )");
    }
}

void check_param()
{
    auto &simu_param = singleton_c<simu_param_c>::instance();
    auto const &samp_param = singleton_c<samp_param_c>::instance();
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();
    auto &demo_param = singleton_c<demo_param_c>::instance();

    ////       CHECKING OUTPUT PARAMETERS
    if (simu_param.Genepop_ind_file && simu_param.Genepop_group_all_samples)
    {
        throw std::invalid_argument("( Genepop_Ind_File and Genepop_Group_All_Samples can't be set at True simultaneously. Check settings. I exit. )");
    }

    ////        CHECKING DEMOGRAPHIQUE PARAMETERS
    if (simu_param.Nodesize_matrix)
    {
        if (simu_param.Nodesize_matrix_filename.empty())
        {
            simu_param.Nodesize_matrix_filename = std::string("Node_Size_Matrix.txt");
        }
    }
    else
    {
        if (!simu_param.Nodesize_matrix_filename.empty())
        {
            throw std::invalid_argument("( Node_Size_Matrix_Filename=" + simu_param.Nodesize_matrix_filename + " can not be considered because Node_Size_Matrix is not set to true . Check settings. I exit. )");
        }
    }

    ////       CHECKING MUTATIONAL PARAMETERS

    if (muta_param.K_min > muta_param.K_max)
    {
        throw std::invalid_argument("( Allelic_Lower_Bound should always be smaller than Allelic_Upper_Bound. Check settings. I exit. )");
    }

    if (fabs(1.0 - (muta_param.Equi_base_freq.A + muta_param.Equi_base_freq.G + muta_param.Equi_base_freq.C + muta_param.Equi_base_freq.T)) > DBL_EPSILON)
    {
        throw std::invalid_argument("( Equilibrium_Frequencies should sum up to 1. Check settings. I exit. )");
    }

    if (!((muta_param.Mod_mut_name == mut_model_enum::iam) || (muta_param.Mod_mut_name == mut_model_enum::kam) || (muta_param.Mod_mut_name == mut_model_enum::smm) || (muta_param.Mod_mut_name == mut_model_enum::gsm)))
    {
        if (simu_param.Genepop_output)
        {
            throw std::invalid_argument("( Genepop format can't deal with sequence data but only wityh allelic data. Check settings. I exit. )");
        }
    }
    else
    {
        if (simu_param.VCF_output)
        {
            throw std::invalid_argument("( VCF format can't deal with allelic data but only sequence data. Check settings. I exit. )");
        }
        if (simu_param.Fasta_output)
        {
            throw std::invalid_argument("( Fasta format can't deal with allelic data but only sequence data. Check settings. I exit. )");
        }
        if (simu_param.Phylip_output)
        {
            throw std::invalid_argument("( Phylip format can't deal with allelic data but only sequence data. Check settings. I exit. )");
        }
    }
    // set equilibrium distribution of allelic/nucleotides ancestral states
    std::uniform_int_distribution<int> distrib;
    if (muta_param.Mod_mut_name == mut_model_enum::iam)
    {
        distrib = std::uniform_int_distribution<int>(1, 1);
    }
    else if ((muta_param.Mod_mut_name == mut_model_enum::kam) || (muta_param.Mod_mut_name == mut_model_enum::smm) || (muta_param.Mod_mut_name == mut_model_enum::gsm))
    {
        distrib = std::uniform_int_distribution<int>(muta_param.K_min, muta_param.K_max);
    }
    else
    {
        distrib = std::uniform_int_distribution<int>(1, 4);
    }

    ////       CHECKING Ancestry sequence settings
    // IAM special case
    if (muta_param.Mod_mut_name == mut_model_enum::iam)
    {
        if (!muta_param.Ancestry_seq.empty())
        {
            throw std::logic_error("( No MRCA allelic state is needed for IAM model. Check settings. I exit. )");
        }
        muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(samp_param.Sequence_length, 1));
    }
    else
    {
        //If MRCA_Sequence not provide or just one zero
        if (muta_param.Ancestry_seq.empty() || ((muta_param.Ancestry_seq.size() == 1) && (muta_param.Ancestry_seq[0].size() == 1) && (muta_param.Ancestry_seq[0][0] == 0)))
        {
            muta_param.Ancestry_seq = std::vector<std::vector<int>>(samp_param.Chr_nbr, std::vector<int>(1, 0));
        }
        //If MRCA_Sequence number and Chromosome number don't match

        if (static_cast<int>(muta_param.Ancestry_seq.size()) != samp_param.Chr_nbr)
        {
            throw std::invalid_argument("( " + std::to_string(muta_param.Ancestry_seq.size()) + " MRCA_Allelic/Nucleotidic_Sequence have been provide but " + std::to_string(samp_param.Chr_nbr) + " chromosomes are needed. Provide all (complete sequence or with a 0) or provide none. Check MRCA_sequence settings. I exit. )");
        }
        //Checks if ancestry size for each chromosome are equal to sequence size parameter and if the specified allelic bounds contain the state of the MRCA
        else if ((muta_param.Mod_mut_name == mut_model_enum::kam) || (muta_param.Mod_mut_name == mut_model_enum::smm) || (muta_param.Mod_mut_name == mut_model_enum::gsm))
        {
            int chr = 1;
            for (auto &ancestry : muta_param.Ancestry_seq)
            {
                if ((ancestry.size() == 1) && (ancestry[0] == 0)) // random assignation of MRCA states for each locus/site
                {
                    std::cout << "Random assignation of Ancestral (MRCA) allelic states for the " + std::to_string(samp_param.Sequence_length) + " loci for chromosome " + std::to_string(chr) + "." << std::endl;
                    if (simu_param.Wait_for_cin_input)
                    {
                        std::cout << "Press any key to resume." << std::endl;
                        std::cin.get();
                    }
                    ancestry.resize(samp_param.Sequence_length);
                    for (int i = 0; i < samp_param.Sequence_length; ++i)
                    {
                        ancestry.at(i) = distrib(rand_gen.Seed_gen);
                    }
                }
                else if (ancestry.size() != static_cast<std::size_t>(samp_param.Sequence_length))
                {
                    throw std::logic_error("( Ancestral (MRCA) allelic state should have the same lenght than the value provide in Sequence_size. Check " +
                                           std::to_string(chr) + "th chromosome MRCA_Sequence settings. I exit. )");
                }
                else
                {
                    for (auto &state : ancestry)
                    {
                        if (state == 0)
                        {
                            state = distrib(rand_gen.Seed_gen);
                        }
                        else if (state < muta_param.K_min || state > muta_param.K_max)
                        {
                            throw std::logic_error("( Ancestral (MRCA) allelic state should always be inbetween Allelic_Lower_Bound and Allelic_Upper_Bound values. Check " + std::to_string(chr) + "th chromosome MRCA_Sequence settings. I exit. )");
                        }
                    }
                }
                ++chr;
            }
        }
        // Sequence data
        else
        {
            int chr = 1;
            for (auto &ancestry : muta_param.Ancestry_seq)
            {
                if ((ancestry.size() == 1) && (ancestry[0] == 0)) // random assignation of MRCA states for each site
                {
                    std::cout << "Random assignation " + std::to_string(chr) + " chromosome MRCA nucleotidic states. Press any key to resume." << std::endl;
                    if (simu_param.Wait_for_cin_input)
                    {
                        std::cin.get();
                    }
                    ancestry.resize(samp_param.Sequence_length);
                    for (int i = 0; i < samp_param.Sequence_length; ++i)
                    {
                        ancestry.at(i) = distrib(rand_gen.Seed_gen);
                    }
                }
                else if (ancestry.size() != static_cast<std::size_t>(samp_param.Sequence_length))
                {
                    throw std::invalid_argument("( Ancestral " + std::to_string(chr) + " chromosome MRCA_Sequence and simulated sequences (Sequence_size) don't have the same size : " + std::to_string(ancestry.size()) + " vs. " + std::to_string(samp_param.Sequence_length) + ", and MRCA_sequence is not 0 (setting for random assignation of MRCA nucleotide states). Check MRCA_sequence and Sequence_size settings. I exit. )");
                }
                //Check if some zero have been put in ancestral seq
                else
                {
                    for (auto &anc_mut_stat : ancestry)
                    {
                        if (anc_mut_stat == 0)
                        {
                            anc_mut_stat = distrib(rand_gen.Seed_gen);
                        }
                    }
                }
            }
            ++chr;
        }
    }

    ////       CHECKING Nodesize and Migration matrix settings

    if (simu_param.Migration_matrix)
    {
        if (simu_param.Migration_matrix_filename.empty())
        {
            simu_param.Migration_matrix_filename = std::string("Migration_Matrix.txt");
        }
        else
        {
            if (!simu_param.Migration_matrix_filename.empty())
            {
                throw std::invalid_argument("( Migration_Matrix_Filename=" + simu_param.Migration_matrix_filename + " will not be considered because Migration_Matrix is not set to true . Check settings. I exit. )");
            }
            if ((demo_param.Dispersal_distrib != dispersal_distrib_enum::none) || (demo_param.Proba_migr != 0.0))
            {
                map_string_dispersal_distrib_name_c d;
                throw std::logic_error("( Both [ Disp_Dist_Max=" + std::to_string(demo_param.Disp_dist_max[0]) +
                                       "or Dispersal_Distribution=" + d.ReverseMap[demo_param.Dispersal_distrib] +
                                       "or Total_emigration_Rate=" + std::to_string(demo_param.Proba_migr) +
                                       " ] and [ Migration_matrix_filename=" + simu_param.Migration_matrix_filename +
                                       " and Migration_Matrix=" + std::to_string(simu_param.Migration_matrix) +
                                       " ] are defined in the settings. Those options are mutually exclusive. Check dispersal settings. I exit. )");
            }
            if (simu_param.Migration_matrix)
            {
                auto file_str = read_file(simu_param.Migration_matrix_filename);
                std::vector<std::string> file_line_str_vec = slice_unix_windows_file_by_line(file_str);

                if (file_line_str_vec.size() != static_cast<std::size_t>((demo_param.Lattice_size[0] + 1) * (demo_param.Lattice_size[1] + 1)))
                {
                    throw std::logic_error("( The Migration_Matrix file (" + simu_param.Migration_matrix_filename + ") does not have a number of lines corresponding to the number of subpopulations on the lattice = Lattice_Size_X *Lattice_Size_Y = " + std::to_string((demo_param.Lattice_size[0] + 1) * (demo_param.Lattice_size[1] + 1)) + ".Check Migration Matrix settings.I exit. )");
                }

                demo_param.Migration_mat.resize((demo_param.Lattice_size[0] + 1) * (demo_param.Lattice_size[1] + 1));

                for (std::size_t i = 0; i < file_line_str_vec.size(); ++i)
                {
                    std::string line_str = replace_spaces_tab_underscores_by_comma(file_line_str_vec.at(i));
                    auto line_double_vec = double_vector_parse_by_comma(line_str);
                    if (line_double_vec.size() != static_cast<std::size_t>((demo_param.Lattice_size[0] + 1) * (demo_param.Lattice_size[1] + 1)))
                    {
                        throw std::logic_error("( In the Migration_Matrix file (" + simu_param.Migration_matrix_filename + "), line " + std::to_string(i + 1) + " does not have a number of elements corresponding to the number of subpopulations on the lattice = Lattice_Size_X * Lattice_Size_Y = " + std::to_string((demo_param.Lattice_size[0] + 1) * (demo_param.Lattice_size[1] + 1)) + ". Check settings. I exit. )");
                    }
                    double sum{0.0};
                    for (std::size_t j = 0; j < line_double_vec.size(); ++j)
                    {
                        sum += line_double_vec.at(j);
                    }
                    if (fabs(sum - 1.0) > DBL_EPSILON)
                    {
                        throw std::logic_error("( In the Migration_Matrix file (" + simu_param.Migration_matrix_filename + "), elements of line " + std::to_string(i + 1) + " does not sum to 1.0\n" + line_str + "\nCheck Migration Matrix settings. I exit. )");
                    }

                    for (std::size_t j = 0; j < line_double_vec.size(); ++j)
                    {
                        demo_param.Migration_mat.at(i).push_back(line_double_vec.at(j));
                    }
                }
            }
            if ((simu_param.Nodesize_matrix) && (demo_param.Pop_size_per_node != 0))
            {
                throw std::logic_error("( Both Ind_Per_Node=" + std::to_string(demo_param.Pop_size_per_node) + " and Nodesize_matrix_filename=" + simu_param.Nodesize_matrix_filename + " Nodesize_Matrix=" + std::to_string(simu_param.Nodesize_matrix) + " are defined in the settings. Those options are mutually exclusive. Check settings. I exit. )");
            }
            if (simu_param.Nodesize_matrix)
            {
                auto file_str = read_file(simu_param.Nodesize_matrix_filename);
                std::vector<std::string> file_line_str_vec = slice_unix_windows_file_by_line(file_str);

                if (file_line_str_vec.size() != static_cast<std::size_t>(demo_param.Lattice_size[1] + 1))
                {
                    throw std::logic_error("( The subpopsize_matrix file (" + simu_param.Nodesize_matrix_filename + ") does not have a number of lines corresponding to Lattice_Size_Y=" + std::to_string(demo_param.Lattice_size[1] + 1) + ". Check settings. I exit. )");
                }

                demo_param.Nodesize_mat.resize(demo_param.Lattice_size[1] + 1);

                for (int i = file_line_str_vec.size() - 1; i >= 0; i--)
                {
                    std::string line_str = replace_spaces_tab_underscores_by_comma(file_line_str_vec.at(i));
                    auto line_double_vec = double_vector_parse_by_comma(line_str);
                    if (line_double_vec.size() != static_cast<std::size_t>(demo_param.Lattice_size[0] + 1))
                    {
                        throw std::logic_error("( In the subpopsize_matrix file (" + simu_param.Nodesize_matrix_filename + "), line " + std::to_string(i + 1) + " does not have a number of elements corresponding to Lattice_Size_X=" + std::to_string(demo_param.Lattice_size[1] + 1) + ". Check settings. I exit. )");
                    }
                    double intpart;
                    for (std::size_t j = 0; j < line_double_vec.size(); ++j)
                    {
                        if (std::modf(line_double_vec.at(j), &intpart) != 0.0)
                        {
                            throw std::logic_error("( In the subpopsize_matrix file (" + simu_param.Nodesize_matrix_filename + "), line " + std::to_string(i + 1) + " contains non-integer values:\n" + line_str + "\nCheck Nodesize_Matrix settings. I exit. )");
                        }
                    }

                    auto line_int_vec = int_vector_parse_by_comma(line_str);
                    if (line_int_vec.size() != static_cast<std::size_t>(demo_param.Lattice_size[0] + 1))
                    {
                        throw std::logic_error("( In the subpopsize_matrix file (" + simu_param.Nodesize_matrix_filename + "), line " + std::to_string(i + 1) + " does not have a number of elements corresponding to Lattice_Size_X = " +
                                               std::to_string(demo_param.Lattice_size[1] + 1) + ". Check settings. I exit. )");
                    }
                    for (std::size_t j = 0; j < line_int_vec.size(); ++j)
                    {
                        demo_param.Nodesize_mat.at(i).push_back(line_int_vec.at(j));
                    }
                }
            }
        }

        ////       CHECKING LATTICE AND SAMPLE  PARAMETERS
        if ((samp_param.Sample_size_per_node.size() > 1) && (static_cast<int>(samp_param.Sample_size_per_node.size()) != (demo_param.Nbr_node_sampled_x * demo_param.Nbr_node_sampled_y)))
        {
            throw std::logic_error("( The size of the vector Ind_Per_Node_Sampled is neither 1 (same sample size for each node) nor not equal to Sample_Size_X * Sample_Size_Y (specific sample szie for each node). Check sample settings. I exit. )");
        }
        if (demo_param.Min_sample_coord_x < 0 || demo_param.Min_sample_coord_x > demo_param.Lattice_size[0] + 1 ||
            demo_param.Min_sample_coord_y < 0 || demo_param.Min_sample_coord_y > demo_param.Lattice_size[1] + 1)
        {
            throw std::logic_error("( Problem with sample position on the lattice: Min_Sample_Coordinate_X/Y < 1 or > lattice dimensions. Check sample settings. I exit. )");
        }

        if (demo_param.Lattice_size[0] < ((demo_param.Nbr_node_sampled_x * demo_param.Void_sample_nodes_X) + demo_param.Min_sample_coord_x - 1))
        {
            throw std::logic_error("( Habitat dimension Lattice_Size_X < Sample_Size_X*Void_Sample_Nodes_X + Min_Sample_Coord_X. Check sample settings. I exit. )");
        }
        if (demo_param.Lattice_size[1] < ((demo_param.Nbr_node_sampled_y * demo_param.Void_sample_nodes_Y) + demo_param.Min_sample_coord_y - 1))
        {
            throw std::logic_error("( Habitat dimension Lattice_Size_Y < Sample_Size_Y*Void_Sample_Nodes_Y + Min_Sample_Coord_Y. Check sample settings. I exit. )");
        }
        if (samp_param.Sample_coord_x.size() != samp_param.Sample_coord_x.size())
        {
            throw std::logic_error("( Sample_Coordinates_X size (" + std::to_string(samp_param.Sample_coord_x.size()) + ") != Sample_coordinates_y size (" + std::to_string(samp_param.Sample_coord_y.size()) + "). Check sample settings. I exit. )");
        }
    }
}

void apply_param()
{
    static bool executed = false;
    //TODO @Raph : Je ne comprend pas l'utilité ?
    if (executed)
    {
        throw std::logic_error("( apply_param() has been called twice... I exit. )");
    }
    executed = true;

    //Call uniq instance of simu_param_c, demo_param_c, demo_param_c, recomb_param_c, muta_param_c
    auto &demo_param = singleton_c<demo_param_c>::instance();
    auto &samp_param = singleton_c<samp_param_c>::instance();
    auto &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &muta_param = singleton_c<muta_param_c>::instance();

    //No dispertion if lattice_size < 2
    for (int dim = 0; dim < 2; ++dim)
    {
        if (demo_param.Lattice_size.at(dim) < 2)
        {
            demo_param.Disp_dist_max.at(dim) = 0;
        }
    }

    if (demo_param.Dispersal_distrib == dispersal_distrib_enum::uniform)
    {
        demo_param.Disp_func.at(0) = fwd_disp_distrib_c::uniform();
        demo_param.Disp_func.at(1) = fwd_disp_distrib_c::uniform();
    }
    else if (demo_param.Dispersal_distrib == dispersal_distrib_enum::gaussian)
    {
        demo_param.Disp_func.at(0) = fwd_disp_distrib_c::gaussian(demo_param.Disp_dist_max.at(0));
        demo_param.Disp_func.at(1) = fwd_disp_distrib_c::gaussian(demo_param.Disp_dist_max.at(1));
    }
    else if (demo_param.Dispersal_distrib == dispersal_distrib_enum::geometric)
    {
        demo_param.Disp_func.at(0) = fwd_disp_distrib_c::geometric(demo_param.G_geo_param);
        demo_param.Disp_func.at(1) = fwd_disp_distrib_c::geometric(demo_param.G_geo_param);
    }
    else if (demo_param.Dispersal_distrib == dispersal_distrib_enum::pareto)
    {
        demo_param.Disp_func.at(0) = fwd_disp_distrib_c::pareto(demo_param.N_pareto_param);
        demo_param.Disp_func.at(1) = fwd_disp_distrib_c::pareto(demo_param.N_pareto_param);
    }
    else
    {
        demo_param.Disp_func.at(0) = fwd_disp_distrib_c::gamma_or_sichel(demo_param, 0);
    }

    demo_param.Population_size_N = demo_param.Pop_size_per_node * (demo_param.Lattice_size[0] + 1) * (demo_param.Lattice_size[1] + 1);

    //2.N.mu is for haploid pop
    recomb_param.Scaled_recomb_rate_rho = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * recomb_param.Unscaled_recomb_rate;
    muta_param.Scaled_mut_rate_theta = samp_param.Ploidy * 2.0 * demo_param.Population_size_N * muta_param.Unscaled_mut_rate_mu;

    compute_sample_coord(samp_param, demo_param);
}

void compute_sample_coord(samp_param_c &samp_param, demo_param_c const &demo_param)
{
    //If Sample_coord are use
    if (!samp_param.Sample_coord_y.empty())
    {
        samp_param.n_total_sample_size = samp_param.Sample_coord_y.size();
        samp_param.Sample_coord_vec.resize(samp_param.Sample_coord_y.size());

        for (std::size_t i = 0; i < samp_param.Sample_coord_y.size(); ++i)
        {
            samp_param.Sample_coord_vec[i] = std::array<int, 2>{samp_param.Sample_coord_x[i] - 1, samp_param.Sample_coord_y[i] - 1};
        }
        std::sort(samp_param.Sample_coord_vec.begin(), samp_param.Sample_coord_vec.end(),
                  [](std::array<int, 2> const &a, std::array<int, 2> const &b) -> bool {
                      if (a.at(0) < b.at(0))
                      {
                          return true;
                      }
                      else
                      {
                          if (a.at(0) == b.at(0))
                          {
                              return a.at(1) < b.at(1);
                          }
                          else
                          {
                              return false;
                          }
                      }
                  });
        samp_param.Sample_coord_x.clear();
        samp_param.Sample_coord_y.clear();
    }
    else
    {
        int nbr_node = demo_param.Nbr_node_sampled_x * demo_param.Nbr_node_sampled_y;
        //WARNING : Sample number (haploide or diploid) in nbr of individuals not genes !

        if (samp_param.Sample_size_per_node.size() == 1) // same sample size for each sampled node
        {
            samp_param.n_total_sample_size = samp_param.Sample_size_per_node[0] * nbr_node;
        }
        else
        {
            samp_param.n_total_sample_size = 0;
            for (auto node_sample_size : samp_param.Sample_size_per_node)
            {
                samp_param.n_total_sample_size += node_sample_size;
            }
        }

        samp_param.Sample_coord_vec.reserve(samp_param.n_total_sample_size);
        int x = demo_param.Min_sample_coord_x;
        int y = demo_param.Min_sample_coord_y;
        auto itr = samp_param.Sample_size_per_node.begin();

        while (nbr_node > 0)
        {
            int nbr_sample = *itr;

            while (nbr_sample > 0)
            {
                samp_param.Sample_coord_vec.push_back({x, y});
                --nbr_sample;
            }

            y += demo_param.Void_sample_nodes_Y;
            //
            if (y == demo_param.Nbr_node_sampled_y * demo_param.Void_sample_nodes_Y + demo_param.Min_sample_coord_y)
            {
                y = demo_param.Min_sample_coord_y;
                x += demo_param.Void_sample_nodes_X;
            }
            --nbr_node;

            if (samp_param.Sample_size_per_node.size() != 1) // same sample size for each sampled node
            {
                ++itr;
            }
        }
    }
}

void output_screen_info(std::string const &version, simu_param_c const &simu_param, samp_param_c const &samp_param, demo_param_c const &demo_param, muta_param_c const &muta_param, recomb_param_c const &recomb_param)
{
    map_string_mod_mut_name_c m;
    map_string_edge_effect_name_c e;

    std::cout << "\n\n==========================================================================" << std::endl;
    std::cout << "         This is  GSpace  v" << version << "    " << std::endl;
    std::cout << "               (Virgoulay et al. 2020 Bioinformatics)                       " << std::endl;
    std::cout << "         an exact coalescent simulator of genetic /  genomic data           " << std::endl;
    std::cout << "            under generalized models of isolation by distance               " << std::endl;
    std::cout << "============================================================================" << std::endl;
    std::cout << "Settings summary : Generic output filename is " << simu_param.Generic_data_filename << std::endl;
    std::cout << " Simulation of " << simu_param.Repetition_nbr << " data sets" << std::endl;
    std::cout << "   with " << samp_param.Chr_nbr << " chromosomes / independant loci";
    std::cout << " with " << samp_param.Sequence_length << " linked sites /  loci each. " << std::endl;
    std::cout << " Mutation model is " << m.ReverseMap[muta_param.Mod_mut_name] << std::endl;
    std::cout << "   with a mutation rate of " << muta_param.Unscaled_mut_rate_mu << " mutations per site per generation." << std::endl;
    std::cout << "   and a recombination rate of " << recomb_param.Unscaled_recomb_rate << " between adjacent sites per generation." << std::endl;
    if (samp_param.Sample_size_per_node.size() == 1)
    {
        std::cout << "Homogeneous sample of size (" << demo_param.Nbr_node_sampled_x << "x";
        std::cout << demo_param.Nbr_node_sampled_y << ")*" << samp_param.Sample_size_per_node[0] << " = " << demo_param.Nbr_node_sampled_x * demo_param.Nbr_node_sampled_y * samp_param.Sample_size_per_node[0];
    }
    else
    {
        std::cout << "User specified sample size for each sampled node on a squarre (" << demo_param.Nbr_node_sampled_x << " x " << demo_param.Nbr_node_sampled_y << ") :" << std::endl;
        size_t node = 0;
        for (auto node_sample_size : samp_param.Sample_size_per_node)
        {
            std::cout << " " << node_sample_size << "(" << samp_param.Sample_coord_vec.at(node)[0] << "," << samp_param.Sample_coord_vec.at(node)[1] << ")";
            node += node_sample_size;
        }
        std::cout << std::endl;
        std::cout << " with a total sample size of : " << samp_param.n_total_sample_size;
    }
    if (samp_param.Ploidy == 1)
    {
        std::cout << " haploid individuals " << std::endl;
    }
    if (samp_param.Ploidy == 2)
    {
        std::cout << " diploid individuals " << std::endl;
    }
    std::cout << "evolving on a " << demo_param.Lattice_size[0] + 1 << " x " << demo_param.Lattice_size[1] + 1 << " lattice";
    if (simu_param.Nodesize_matrix)
    {
        std::cout << std::endl
                  << " where the size of each subpopulation is given in the file " << simu_param.Nodesize_matrix_filename << "." << std::endl;
    }
    else
    {
        std::cout << " with " << e.ReverseMap[demo_param.Edge_effects] << " boundaries" << std::endl;
        std::cout << "  where each node carries " << demo_param.Pop_size_per_node << " individuals." << std::endl;
    }
    if (simu_param.Migration_matrix)
    {
        std::cout << "Custom migration matrix given in the file " << simu_param.Migration_matrix_filename << "." << std::endl;
    }
    std::cout << "Dispersal settings are summarized in the "<< simu_param.Generic_data_filename + simu_param.Param_summary_filename <<" file. " << std::endl;
    std::cout << "============================================================================\n"
              << std::endl;
}
