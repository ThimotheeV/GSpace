#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream> //to convert double to String with specified precision
#include <string>
#include <iostream>
#include <fstream>

#include "output.hpp"
#include "mutation.hpp"
#include "simulator.hpp"

extern double const DBL_EPSILON;

void output_stat_files(simu_param_c const &simu_p, info_collector_c const &info, output_stat_c const &output, int rep, int ploidy)
{
    std::vector<std::string> head;
    int nbr_class = output.Prob_id_1_loc_Qr_res.Qr.size();
    head.reserve(13 + 3 * nbr_class);

    head.emplace_back("Sim");

    if (info.Prob_id_1_2_loc)
    {
        head.emplace_back("Prob_id_1_2_loc");
        if (info.Iterative_stats)
        {
            head.emplace_back("Prob_id_1_2_loc_cumul_mean");
            head.emplace_back("Prob_id_1_2_loc_cumul_var");
        }
    }
    if (info.Prob_id_1_loc_Qr)
    {
        for (int dist = 0; dist < nbr_class; ++dist)
        {
            std::string temp = "Q" + std::to_string(dist);
            head.emplace_back(temp);

            if (info.Iterative_stats)
            {
                std::string temp_m = temp + "_cumul_mean";
                head.emplace_back(temp_m);
                std::string temp_v = temp + "_cumul_var";
                head.emplace_back(temp_v);
            }
        }
    }
    if (info.Prob_id_1_loc_Qwi_wd_bd)
    {
        if (ploidy > 1)
        {
            head.emplace_back("Qwi");
            if (info.Iterative_stats)
            {
                head.emplace_back("Qwi_cumul_mean");
                head.emplace_back("Qwi_cumul_var");
            }
        }

        head.emplace_back("Qwd");
        if (info.Iterative_stats)
        {
            head.emplace_back("Qwd_cumul_mean");
            head.emplace_back("Qwd_cumul_var");
        }

        head.emplace_back("Qbd");
        if (info.Iterative_stats)
        {
            head.emplace_back("Qbd_cumul_mean");
            head.emplace_back("Qbd_cumul_var");
        }
    }

    //+++++++++++++++++++++++++++++++++++++++++++++++//
    head.shrink_to_fit();
    if (rep == 0)
    {
        print_output(simu_p.Generic_data_filename+"_Global_Stats.txt", head,"over");
    }
    //for output
    rep += 1; // RL risqué de changer ca ici, plutot faire une variable locale.
    std::vector<double> stats_run;
    stats_run.reserve(13 + 3 * nbr_class);

    stats_run.push_back(rep);

    if (info.Prob_id_1_2_loc)
    {
        stats_run.push_back(output.Prob_id_1_2_loc_res.PHI);
        if (info.Iterative_stats)
        {
            stats_run.push_back(output.Prob_id_1_2_loc_res.PHI_cumul_m_v.at(0));
            stats_run.push_back(output.Prob_id_1_2_loc_res.PHI_cumul_m_v.at(1));
        }
    }
    if (info.Prob_id_1_loc_Qr)
    {
        for (int dist = 0; dist < nbr_class; ++dist)
        {
            stats_run.push_back(output.Prob_id_1_loc_Qr_res.Qr[dist]);

            if (info.Iterative_stats)
            {
                stats_run.push_back(output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[dist].at(0));
                stats_run.push_back(output.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[dist].at(1));
            }
        }
    }
    if (info.Prob_id_1_loc_Qwi_wd_bd)
    {
        if (ploidy > 1)
        {
            stats_run.push_back(output.Prob_id_1_loc_Qwi_wd_bd_res.Qwi);
            if (info.Iterative_stats)
            {
                stats_run.push_back(output.Prob_id_1_loc_Qwi_wd_bd_res.Qwi_cumul_m_v.at(0));
                stats_run.push_back(output.Prob_id_1_loc_Qwi_wd_bd_res.Qwi_cumul_m_v.at(1));
            }
        }

        stats_run.push_back(output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd);
        if (info.Iterative_stats)
        {
            stats_run.push_back(output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(0));
            stats_run.push_back(output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v.at(1));
        }

        stats_run.push_back(output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd);
        if (info.Iterative_stats)
        {
            stats_run.push_back(output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(0));
            stats_run.push_back(output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v.at(1));
        }
    }
    stats_run.shrink_to_fit();
    print_output(simu_p.Generic_data_filename + "_Global_Stats.txt", stats_run, "app");

    //+++++++++++++++++++++++++++++++++++++++++++++++//
    head.clear();
    head.reserve(4 + nbr_class);

    if (info.Prob_id_1_2_loc)
    {
        head.emplace_back("Prob_id_1_2_loc");
    }
    if (info.Prob_id_1_loc_Qr)
    {
        for (int dist = 0; dist < nbr_class; ++dist)
        {
            std::string temp = "Q" + std::to_string(dist);
            head.push_back(temp);
        }
    }
    if (info.Prob_id_1_loc_Qwi_wd_bd)
    {
        if (ploidy > 1)
        {
            head.emplace_back("Qwi");
        }
        head.emplace_back("Qwd");
        head.emplace_back("Qbd");
    }

    if (info.Stats_per_loc_per_simu)
    {
        int chr_nbr = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_by_chr_by_loc.size();
        int loc_nbr = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_by_chr_by_loc[0].size();
        for (int chr = 0; chr < chr_nbr; ++chr)
        {
            //print head
            std::string temp = simu_p.Generic_data_filename + "_" + std::to_string(rep) + "_Chr_" + std::to_string(chr + 1) + "_Stats_per_loc.txt";
            head[0] = "Locus";
            print_output(temp, head,"over");

            for (int loc = 0; loc < loc_nbr; ++loc)
            {
                std::vector<double> stats_loc;
                stats_loc.reserve(4 + nbr_class);

                stats_loc.push_back(loc + 1);

                if (info.Prob_id_1_loc_Qr)
                {
                    for (int dist = 0; dist < nbr_class; ++dist)
                    {
                        stats_loc.push_back(output.Prob_id_1_loc_Qr_res.Qr_by_chr_by_loc[chr][loc][dist]);
                    }
                }
                if (info.Prob_id_1_loc_Qwi_wd_bd)
                {
                    if (ploidy > 1)
                    {
                        stats_loc.push_back(output.Prob_id_1_loc_Qwi_wd_bd_res.Qwi_by_chr_by_loc[chr][loc]);
                    }
                    stats_loc.push_back(output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_by_chr_by_loc[chr][loc]);
                    stats_loc.push_back(output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_by_chr_by_loc[chr][loc]);
                }

                stats_loc.shrink_to_fit();
                print_output(temp, stats_loc,"app");
            }
        }
    }

    //+++++++++++++++++++++++++++++++++++++++++++++++//
    head.clear();
    head.reserve(4 + nbr_class);

    if (info.Prob_id_1_2_loc)
    {
        head.emplace_back("Prob_id_1_2_loc_mean");
    }
    if (info.Prob_id_1_loc_Qr)
    {
        for (int dist = 0; dist < nbr_class; ++dist)
        {
            std::string temp = "Q" + std::to_string(dist) + "_mean";
            head.push_back(temp);
        }
    }
    if (info.Prob_id_1_loc_Qwi_wd_bd)
    {
        if (ploidy > 1)
        {
            head.emplace_back("Qwi_mean");
        }
        head.emplace_back("Qwd_mean");
        head.emplace_back("Qbd_mean");
    }
    if (info.Stats_per_chr_per_simu)
    {
        int chr_nbr = output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_by_chr_by_loc.size();
        //print head
        std::string temp = simu_p.Generic_data_filename + "_" + std::to_string(rep) + "_Stats_per_chr.txt";
        head[0] = "Chr";
        print_output(temp, head, "over");
        for (int chr = 0; chr < chr_nbr; ++chr)
        {
            std::vector<double> stats_loc;
            stats_loc.reserve(4 + nbr_class);

            stats_loc.push_back(chr + 1);

            if (info.Prob_id_1_loc_Qr)
            {
                for (int dist = 0; dist < nbr_class; ++dist)
                {
                    stats_loc.push_back(output.Prob_id_1_loc_Qr_res.Qr_by_chr[chr][dist]);
                }
            }
            if (info.Prob_id_1_loc_Qwi_wd_bd)
            {
                if (ploidy > 1)
                {
                    stats_loc.push_back(output.Prob_id_1_loc_Qwi_wd_bd_res.Qwi_by_chr[chr]);
                }

                stats_loc.push_back(output.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_by_chr[chr]);
                stats_loc.push_back(output.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_by_chr[chr]);
            }

            stats_loc.shrink_to_fit();
            print_output(temp, stats_loc, "app");
        }
    }
    if (info.Effective_disp)
    {
        //print head
        head.clear();
        head.resize(3);
        std::string filename = simu_p.Generic_data_filename + "_" + std::to_string(rep) + "_Emp_disp.txt";
        head[0] = "Step";
        head[1] = "X";
        if( info.Emp_axial_disp.size() == 2 )
        {
            head[2] = "Y";
        }
        print_output(filename, head, "over");
        
        std::ostringstream streamObj;
        // Set precision to 2 digits
        streamObj << std::setprecision(15) << std::setw(17);
        
        std::vector<std::string> stats_run(3,"");
        std::size_t max_size;
        if (info.Emp_axial_disp.size() == 2)
        {
            max_size = std::max(info.Emp_axial_disp.at(0).size(),info.Emp_axial_disp.at(1).size());
        }
        else{
            max_size = info.Emp_axial_disp.at(0).size();
        }
        for (std::size_t i = 0; i < max_size; ++i)
        {
            stats_run[0] += std::to_string(static_cast<int>(i - (max_size - 1)/2));
            if (fabs(static_cast<int>(i - (max_size - 1)/2)) <= ((info.Emp_axial_disp.at(0).size() - 1) / 2))
            {
                //Add double to stream to set the precision & width
                streamObj.str(""); streamObj.clear();
                if (info.Emp_axial_disp.at(0).size() < max_size)
                {
                    streamObj << info.Emp_axial_disp.at(0).at(i - ((max_size - 1)/2 - (info.Emp_axial_disp.at(0).size() - 1)/2));
                }
                else
                {
                    streamObj << info.Emp_axial_disp.at(0).at(i);
                }
                stats_run[0] += "\t" + streamObj.str();
            }
            else
            {
                stats_run[0] += "\t0.00000000000000000";
            }
            if (info.Emp_axial_disp.size() == 2)
            {
                if (fabs(static_cast<int>(i - (max_size - 1)/2)) <= ((info.Emp_axial_disp.at(1).size() - 1) /  2))
                {
                    //Add double to stream to set the precision & width
                    streamObj.str(""); streamObj.clear();
                    if (info.Emp_axial_disp.at(1).size() < max_size)
                    {
                        streamObj << info.Emp_axial_disp.at(1).at(i - ((max_size - 1)/2 - (info.Emp_axial_disp.at(1).size() - 1)/2));
                    }
                    else
                    {
                        streamObj << info.Emp_axial_disp.at(1).at(i);
                    }
                    stats_run[0] += "\t" + streamObj.str();
                }
                else
                {
                    //Add double to stream to set the precision & width
                    stats_run[0] += "\t0.00000000000000000";
                }
            }
//            if (i != (max_size - 1))
//            {
                stats_run[0] += "\n";
//            }
        }
        //Add double to stream to set the precision & width
        streamObj.str(""); streamObj.clear();
        streamObj << info.Mean_emp_axial_disp.at(0);
        stats_run[2] += "\nEmpirical dispersal mean (X [, Y]) :\t" + streamObj.str();
        if (info.Mean_emp_axial_disp.size() == 2)
        {
            //Add double to stream to set the precision & width
            streamObj.str(""); streamObj.clear();
            streamObj << info.Mean_emp_axial_disp.at(1);
            stats_run[2] += "\t" + streamObj.str();
        }
        stats_run[2] += "\n";
        //Add double to stream to set the precision & width
        streamObj.str(""); streamObj.clear();
        streamObj << info.Sig_emp_axial_disp.at(0);
        stats_run[2] += "Empirical dispersal sigma2 (X [, Y]) :\t" + streamObj.str();
        if (info.Sig_emp_axial_disp.size() == 2)
        {
         //Add double to stream to set the precision & width
         streamObj.str(""); streamObj.clear();
         streamObj << info.Sig_emp_axial_disp.at(1);
         stats_run[2] += "\t" + streamObj.str();
        }
        stats_run[2] += "\n";
        //Add double to stream to set the precision & width
        streamObj.str(""); streamObj.clear();
        streamObj << info.Kurt_emp_axial_disp.at(0);
        stats_run[2] += "Empirical dispersal kurtosis (X [, Y]) :\t" + streamObj.str();
         if (info.Kurt_emp_axial_disp.size() == 2)
         {
             //Add double to stream to set the precision & width
             streamObj.str(""); streamObj.clear();
             streamObj << info.Kurt_emp_axial_disp.at(1);
             stats_run[2] += "\t" + streamObj.str();
         }

        print_output(filename, stats_run, "app");
    }
}

//sample_mutated_state_chr => chr, indiv, mut {site, state}
void genepop_output(std::string const &path_to_file, std::vector<std::array<int, 2>> const &coord_vec, std::vector<std::vector<int>> const &ancestry_seq, std::vector<std::vector<std::vector<std::pair<int, int>>>> const &sample_mutated_state_chr, std::string const &mut_mod, int ploidy)
{
    std::ofstream genepop(path_to_file);
    auto const &simu_param = singleton_c<simu_param_c>::instance();

    if (genepop.is_open())
    {
        if (coord_vec.size() * ploidy != sample_mutated_state_chr.at(0).size())
        {
            throw std::logic_error("( In genepop_output() : coord_vec.size() * ploidy != sample_mutated_state.size() (nbr of indiv). Contact the developpers. I exit. )");
        }

        genepop << "This file has been generated by the GSpace program.";

        std::size_t tot_nbr_of_loc = sample_mutated_state_chr.size() * ancestry_seq[0].size();
        for (std::size_t lign = 0; lign < tot_nbr_of_loc; ++lign)
        {
            genepop << std::endl;
            genepop << "locus" << std::to_string(lign + 1) << "_" << mut_mod;
        }
        genepop << std::endl;
        genepop << "pop";

        auto coord_itr = coord_vec.begin();
        std::array<int, 2> previous_coord = *coord_itr;

        // loop over individuals
        for (std::size_t indiv = 0; indiv < sample_mutated_state_chr.at(0).size(); indiv += ploidy)
        {
            //If coord change, pop change, or if ind genepop file format
            if ((!simu_param.Genepop_group_all_samples) && (simu_param.Genepop_ind_file || ((previous_coord.at(1) != coord_itr->at(1)) || (previous_coord.at(0) != coord_itr->at(0)))))
            {
                genepop << std::endl;
                genepop << "pop";
            }

            genepop << std::endl;
            if (simu_param.Genepop_no_coord)
            {
                genepop << indiv / ploidy + 1 << " ,";
            }
            else
            {
                genepop << std::to_string(coord_itr->at(0) + 1) << " " << std::to_string(coord_itr->at(1) + 1) << " ,";
            }

            //Loop over chromosomes
            for (std::size_t chr = 0; chr < sample_mutated_state_chr.size(); ++chr)
            {
                //Full int sequence reconstruction
                std::vector<int> full_seq_indiv = reconstitution_full_int_seq_indiv(ancestry_seq[chr], sample_mutated_state_chr[chr][indiv]);
                std::vector<int> full_seq2_indiv;

                if (ploidy == 2)
                {
                    full_seq2_indiv = reconstitution_full_int_seq_indiv(ancestry_seq[chr], sample_mutated_state_chr[chr][indiv + 1]);
                }

                for (std::size_t locus = 0; locus < ancestry_seq[chr].size(); ++locus)
                {
                    genepop << " ";
                    if (full_seq_indiv[locus] > 999)
                    {
                        throw std::logic_error("( An allele was found to be " + std::to_string(full_seq_indiv[locus]) + ", but Genepop format can't handle alleles larger than 999.  I exit. )");
                    }
                    if (full_seq_indiv[locus] < 10)
                    {
                        genepop << "0";
                    }
                    if (full_seq_indiv[locus] < 100)
                    {
                        genepop << "0";
                    }

                    genepop << std::to_string(full_seq_indiv[locus]);

                    if (ploidy == 2)
                    {
                        if (full_seq2_indiv[locus] > 999)
                        {
                            throw std::logic_error("( An allele was found to be " + std::to_string(full_seq2_indiv[locus]) + ", but Genepop format can't handle alleles larger than 999.  I exit. )");
                        }
                        if (full_seq2_indiv[locus] < 10)
                        {
                            genepop << "0";
                        }
                        if (full_seq2_indiv[locus] < 100)
                        {
                            genepop << "0";
                        }

                        genepop << std::to_string(full_seq2_indiv[locus]);
                    }
                }
            }

            previous_coord = *coord_itr;
            ++coord_itr;
        }
        genepop << std::endl;
        genepop.close();
    }
    else
    {
        throw std::invalid_argument("( Unable to open " + path_to_file + ". I exit. )");
    }
}

// std::vector<std::vector<std::vector<std::pair<int, int>>>> => chr,indiv,mut_site
void vcf_output(std::string const &path_to_file, std::vector<std::array<int, 2>> const &coord_vec, std::vector<std::vector<int>> const &ancestry_seq, std::vector<std::vector<std::vector<std::pair<int, int>>>> const &sample_mutated_state_chr, int ploidy)
{
    //     Chromosome Name
    //     Chromosome Position
    //     ID
    //         This is generally used to reference an annotated variant in dbSNP or other curate variant database.
    //     Reference base(s)
    //         What is the reference’s base at this position
    //     Alternate base(s)
    //         The variants found in your dataset that differ from the reference
    //     Variant Quality
    //         Phred-scaled quality for the observed ALT
    //     Filter
    //         Whether or not this has passed all filters – generally a QC measure in variant calling algorithms
    //     Info
    //         This is for additional information, generally describing the nature of the position/variants with respect to other data.
    //     Format
    //         This is for additional information, describing the nature of information display for each sample.

    std::ofstream vcf(path_to_file);
    if (vcf.is_open())
    {
        map_string_nuc_name_c map_string_nuc_name;
        vcf << "##fileformat=VCFv4.3" << std::endl;
        vcf << "##FORMAT=<ID=GT, Number=1,Type=String,Description=\"Genotype\">" << std::endl;
        for (std::size_t i = 0; i < coord_vec.size(); ++i)
        {
            vcf << "##SAMPLE=<ID=Indiv" << i + 1 << ",Assay=WholeGenome,Description=\"Individual from pop[" << coord_vec[i].at(0) + 1 << ";" << coord_vec[i].at(1) + 1 << "]\">" << std::endl;
        }
        vcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        for (std::size_t i = 0; i < coord_vec.size(); ++i)
        {
            vcf << "\tINDIV" << std::to_string(i + 1);
        }
        vcf << std::endl;
        //Per chromosome
        for (std::size_t chr = 0; chr < sample_mutated_state_chr.size(); ++chr)
        {
            std::size_t nbr_haplotype = sample_mutated_state_chr[chr].size();
            std::vector<std::vector<int>> sample_mutated_state(nbr_haplotype);
            for (std::size_t haplotype = 0; haplotype < nbr_haplotype; ++haplotype)
            {
                sample_mutated_state[haplotype] = reconstitution_full_int_seq_indiv(ancestry_seq[chr], sample_mutated_state_chr[chr][haplotype]);
            }
            for (std::size_t site = 0; site < ancestry_seq[chr].size(); ++site)
            {
                std::vector<int> alt_value;
                alt_value.reserve(ancestry_seq[chr].size());
                alt_value.push_back(ancestry_seq[chr][site]);
                vcf << (chr + 1) << "\t" << (site + 1) << "\t"
                    << "."
                    << "\t"
                    << map_string_nuc_name.Inverse_map.at(ancestry_seq[chr][site]) << "\t";

                for (std::size_t haplotype = 0; haplotype < nbr_haplotype; ++haplotype)
                {
                    bool value_not_know = true;
                    std::size_t index = 0;
                    while (value_not_know && (index < alt_value.size()))
                    {
                        if (sample_mutated_state[haplotype][site] == alt_value[index])
                        {
                            value_not_know = false;
                        }
                        ++index;
                    }
                    if (value_not_know)
                    {
                        alt_value.push_back(sample_mutated_state[haplotype][site]);
                    }
                }

                if (alt_value.size() > 1)
                {
                    for (std::size_t index = 1; index < alt_value.size(); ++index)
                    {
                        if (index < (alt_value.size() - 1))
                        {
                            vcf << map_string_nuc_name.Inverse_map.at(alt_value[index]) << ",";
                        }
                        else
                        {
                            vcf << map_string_nuc_name.Inverse_map.at(alt_value[index]);
                        }
                    }
                }
                else
                {
                    vcf << ".";
                }

                vcf << "\t"
                    << "."
                    << "\t"
                    << "."
                    << "\t"
                    << "."
                    << "\t"
                    << "GT";
                for (std::size_t indiv = 0; indiv < nbr_haplotype; indiv += ploidy)
                {
                    std::size_t index1 = 0;
                    std::size_t index2 = 0;
                    for (std::size_t i = 0; i < alt_value.size(); ++i)
                    {
                        if (sample_mutated_state[indiv][site] == alt_value[i])
                        {
                            index1 = i;
                        }
                        if (ploidy == 2)
                        {
                            if (sample_mutated_state[indiv + 1][site] == alt_value[i])
                            {
                                index2 = i;
                            }
                        }
                    }
                    vcf << "\t" << index1;
                    if (ploidy == 2)
                    {
                        vcf << "|" << index2;
                    }
                }
                vcf << std::endl;
            }
        }
        vcf.close();
    }
    else
    {
        throw std::invalid_argument("( Unable to open " + path_to_file + ". I exit. )");
    }
}

void fasta_output(std::string const &path_to_file, std::string const &generic_data_filename, std::vector<std::array<int, 2>> const &coord_vec, std::vector<std::vector<int>> const &ancestry_seq, std::vector<std::vector<std::vector<std::pair<int, int>>>> const &sample_mutated_state_chr, int ploidy, int rep)
{
    std::ofstream fasta(path_to_file);
    auto const &simu_param = singleton_c<simu_param_c>::instance();

    if (fasta.is_open())
    {
        if (coord_vec.size() * ploidy != sample_mutated_state_chr.at(0).size())
        {
            throw std::logic_error("( In fasta_output() : coord_vec.size() * ploidy != sample_mutated_state.size(). I exit. )");
        }

        auto coord_itr = coord_vec.begin();

        //write the ancestral sequence for each chromosome)
        for (std::size_t chr = 0; chr < sample_mutated_state_chr.size(); ++chr)
        {
            fasta << ">Gspace_" << generic_data_filename << "_" << rep + 1 << "_ancestral_sequence_chrom_" << chr + 1 << "_mutNbr_0" << std::endl;

            //Full nucleotide sequence reconstruction
            std::vector<std::pair<int, int>> anc_seq_mut_seq{};
            std::string full_anc_seq = reconstitution_full_nucl_seq_indiv(ancestry_seq[chr], anc_seq_mut_seq);

            if (simu_param.Fasta_single_line_seq)
            {
                fasta << full_anc_seq << std::endl;
            }
            else
            {
                size_t start_str = 0, site_nbr = 80;
                while (start_str < (full_anc_seq.size()))
                {
                    fasta << full_anc_seq.substr(start_str, site_nbr) << std::endl;
                    start_str += site_nbr;
                }
            }
        }
        //loop over individuals
        for (std::size_t ind = 0; ind < sample_mutated_state_chr.at(0).size(); ind += ploidy)
        {

            //Loop over chromosomes
            for (std::size_t chr = 0; chr < sample_mutated_state_chr.size(); ++chr)
            {
                //write the single line description of the first sequence
                fasta << ">Gspace_" << generic_data_filename << "_" << rep + 1 << "_ind_" << ind / ploidy + 1 << "_coord_" << std::to_string(coord_itr->at(0) + 1) << "_" << std::to_string(coord_itr->at(1) + 1) << "_chrom_" << chr + 1 << "_1_mutNbr_" << sample_mutated_state_chr[chr][ind].size() << std::endl;

                //Full nucleotide sequence reconstruction
                std::string full_seq_indiv = reconstitution_full_nucl_seq_indiv(ancestry_seq[chr], sample_mutated_state_chr[chr][ind]);

                if (simu_param.Fasta_single_line_seq)
                {
                    fasta << full_seq_indiv << std::endl;
                }
                else
                {
                    size_t start_str = 0, site_nbr = 80;
                    while (start_str < (full_seq_indiv.size()))
                    {
                        fasta << full_seq_indiv.substr(start_str, site_nbr) << std::endl;
                        start_str += site_nbr;
                    }
                }

                std::string full_seq2_indiv;
                if (ploidy == 2)
                {
                    //write the single line description of the first sequence
                    fasta << ">Gspace_" << generic_data_filename << "_" << rep + 1 << "_ind_" << ind / ploidy + 1 << "_coord_" << std::to_string(coord_itr->at(0) + 1) << "_" << std::to_string(coord_itr->at(1) + 1) << "_chrom_" << chr + 1 << "_2_mutNbr_" << sample_mutated_state_chr[chr][ind + 1].size() << std::endl;

                    //Full nucleotide sequence reconstruction
                    full_seq2_indiv = reconstitution_full_nucl_seq_indiv(ancestry_seq[chr], sample_mutated_state_chr[chr][ind + 1]);

                    if (simu_param.Fasta_single_line_seq)
                    {
                        fasta << full_seq2_indiv << std::endl;
                    }
                    else
                    {
                        size_t start_str = 0, site_nbr = 80;
                        while (start_str < (full_seq2_indiv.size()))
                        {
                            fasta << full_seq2_indiv.substr(start_str, site_nbr) << std::endl;
                            start_str += site_nbr;
                        }
                    }
                }
            }
            ++coord_itr;
        }
        fasta.close();
    }
    else
    {
        throw std::invalid_argument("( Unable to open " + path_to_file + ". I exit. )");
    }
}

void phylip_output(std::string const &path_to_file, std::string const &generic_data_filename, std::vector<std::array<int, 2>> const &coord_vec, std::vector<std::vector<int>> const &ancestry_seq, std::vector<std::vector<std::pair<int, int>>> const &sample_mutated_state_ind, int ploidy, int chrom)
{
     std::ofstream phylip(path_to_file);

     if (phylip.is_open())
     {
         if (coord_vec.size() * ploidy != sample_mutated_state_ind.size())
         {
             throw std::logic_error("( In phylip_output() : coord_vec.size() * ploidy != sample_mutated_state.size(). I exit. )");
         }

         // + 1 for the ancestral sequence
         phylip << sample_mutated_state_ind.size() + 1 << " " << ancestry_seq[0].size() << std::endl;

         //write the ancestral sequence for each chromosome
         //Full nucleotide sequence reconstruction
         std::vector<std::pair<int, int>> anc_seq_mut_seq{};
         std::string full_anc_seq = reconstitution_full_nucl_seq_indiv(ancestry_seq[chrom], anc_seq_mut_seq);

         phylip.fill(' ');
         phylip.width(10);
         phylip << std::left << "Anc_" + std::to_string(chrom + 1) + " " << full_anc_seq << std::endl;

         //loop over individuals
         for (std::size_t ind = 0; ind < sample_mutated_state_ind.size(); ind += ploidy)
         {
             //Full nucleotide sequence reconstruction
             std::string full_seq_indiv = reconstitution_full_nucl_seq_indiv(ancestry_seq[chrom], sample_mutated_state_ind[ind]);

             phylip.fill(' ');
             phylip.width(10);
             phylip << std::left << std::to_string( (int) ind / ploidy + 1) + "_" + std::to_string(chrom + 1) + "_1 " << full_seq_indiv << std::endl;

             std::string full_seq2_indiv;
             if (ploidy == 2)
             {
                 //Full nucleotide sequence reconstruction
                 full_seq2_indiv = reconstitution_full_nucl_seq_indiv(ancestry_seq[chrom], sample_mutated_state_ind[ind + 1]);

                 phylip.fill(' ');
                 phylip.width(10);
                 phylip << std::left << std::to_string( (int) ind / ploidy + 1) + "_" + std::to_string(chrom + 1) + "_2 " << full_seq2_indiv << std::endl;
             }
         }
         phylip.close();
     }
     else
     {
         throw std::invalid_argument("( Unable to open " + path_to_file + ". I exit. )");
     }
}

void seq_char_output(std::string const &path_to_file, std::vector<std::array<int, 2>> const &coord_vec, std::vector<std::vector<std::vector<std::pair<int, int>>>> const &sample_mutated_state_chr, int ploidy)
{
    std::ofstream seq_char(path_to_file);
    if (seq_char.is_open())
    {
        auto coord_itr = coord_vec.begin();

        seq_char << "ind\tchrom\tphase\tX\tY\tmut_nb" << std::endl;

        for (std::size_t ind = 0; ind < sample_mutated_state_chr.at(0).size(); ind += ploidy)
        {

            //Loop over chromosomes
            for (std::size_t chr = 0; chr < sample_mutated_state_chr.size(); ++chr)
            {
                seq_char << ind / ploidy + 1 << "\t" << chr + 1 << "\t1\t" << std::to_string(coord_itr->at(0) + 1) << "\t" << std::to_string(coord_itr->at(1) + 1) << "\t" << sample_mutated_state_chr[chr][ind].size() << std::endl;

                if (ploidy == 2)
                {
                    seq_char << ind / ploidy + 1 << "\t" << chr + 1 << "\t2\t" << std::to_string(coord_itr->at(0) + 1) << "\t" << std::to_string(coord_itr->at(1) + 1) << "\t" << sample_mutated_state_chr[chr][ind + 1].size() << std::endl;
                }
            }
            ++coord_itr;
        }
        seq_char.close();
    }
    else
    {
        throw std::invalid_argument("( Unable to open " + path_to_file + ". I exit. )");
    }
}

void coord_output(std::string const &path_to_file, std::vector<std::array<int, 2>> const &coord_vec)
{
    std::ofstream coord(path_to_file);
    if (coord.is_open())
    {

        coord << "X Y" << std::endl;
        for (auto coord_itr = coord_vec.begin(); coord_itr < coord_vec.end(); ++coord_itr)
        {
            coord << std::to_string(coord_itr->at(0) + 1) << " " << std::to_string(coord_itr->at(1) + 1) << std::endl;
        }
        coord.close();
    }
    else
    {
        throw std::invalid_argument("( Unable to open " + path_to_file + ". I exit. )");
    }
}

void write_beforerun_param_settings_summary(std::string const &path_to_file, std::string const &version, simu_param_c const &simu_param, samp_param_c const &samp_param, demo_param_c const &demo_param, muta_param_c const &muta_param, recomb_param_c const &recomb_param)
{
    std::ofstream param_summary_f(path_to_file);
    if (param_summary_f.is_open())
    {
        map_string_mod_mut_name_c m;
        map_string_edge_effect_name_c e;
        map_string_dispersal_distrib_name_c d;

        param_summary_f << "GSpace v " << version << std::endl
                        << std::endl;
        param_summary_f << "Random Seeds : " << simu_param.seed << std::endl;
        param_summary_f << "Generic output filename : " << simu_param.Generic_data_filename << std::endl;
        param_summary_f << "Number of simulated data sets : " << simu_param.Repetition_nbr << std::endl;
        param_summary_f << "Number of chromosomes / independant loci : " << samp_param.Chr_nbr << std::endl;
        param_summary_f << "Number of linked sites / loci per chromosome : " << samp_param.Sequence_length << std::endl;
        param_summary_f << "Mutation model : " << m.ReverseMap[muta_param.Mod_mut_name] << std::endl;
        param_summary_f << "Mutation rate : " << muta_param.Unscaled_mut_rate_mu << std::endl;
        param_summary_f << "Recombination rate : " << recomb_param.Unscaled_recomb_rate << std::endl;
        param_summary_f << "Sample size X : " << demo_param.Nbr_node_sampled_x << std::endl;
        param_summary_f << "Sample size Y : " << demo_param.Nbr_node_sampled_y << std::endl;
        param_summary_f << "Number of sampled individuals per node :" << std::endl;
        size_t node=0;
        for (auto node_sample_size : samp_param.Sample_size_per_node)
        {
            param_summary_f <<  " " << node_sample_size << "(" << samp_param.Sample_coord_vec.at(node)[0] << "," << samp_param.Sample_coord_vec.at(node)[1] << ")";
            node+=node_sample_size;

        }
        param_summary_f << std::endl;
        param_summary_f << "Total number of sampled individuals : " << samp_param.n_total_sample_size << std::endl;
        param_summary_f << "Ploidy : " << samp_param.Ploidy << std::endl;
        param_summary_f << "Lattice size X : " << demo_param.Lattice_size[0] + 1 << std::endl;
        param_summary_f << "Lattice size Y : " << demo_param.Lattice_size[1] + 1 << std::endl;
        param_summary_f << "Edge effects : " << e.ReverseMap[demo_param.Edge_effects] << std::endl;
        param_summary_f << "Number of individuals per lattice node : " << demo_param.Pop_size_per_node << std::endl;
        param_summary_f << "Dispersal distribution : " << d.ReverseMap[demo_param.Dispersal_distrib] << std::endl;
        param_summary_f << "Total emigration rate : " << demo_param.Proba_migr << std::endl;
        switch (demo_param.Dispersal_distrib)
        {
        case dispersal_distrib_enum::none:
        {
            break;
        }
        case dispersal_distrib_enum::uniform:
        {
            break;
        }
        case dispersal_distrib_enum::gaussian:
        {
            break;
        }
        case dispersal_distrib_enum::geometric:
        {
            param_summary_f << "Geometric shape : " << demo_param.G_geo_param << std::endl;
            break;
        }
        case dispersal_distrib_enum::pareto:
        {
            param_summary_f << "Pareto shape : " << demo_param.N_pareto_param << std::endl;
            break;
        }
        case dispersal_distrib_enum::sichel:
        {
            param_summary_f << "Sichel shape : Gamma " << demo_param.Sichel_pars[0] << "; Ksi " << demo_param.Sichel_pars[1] << "; Omega " << demo_param.Sichel_pars[2] << std::endl;
            break;
        }
            // TODO sichel
            //        case 'S':
            //{
            //            if (SichelDispPars[2]<0) {
            //                simulpars<<"Inverse Gamma mixture (Chesson & Lee 2005) "<<endl ;
            //                if (SichelDispPars[3]>0) simulpars<<"also including a fraction "<<SichelDispPars[3]<<" of random dispersal [excluding focal patch!]: "<<endl;
            //                simulpars<<"gamma="<<SichelDispPars[0]<<endl ;
            //                simulpars<<"kappa (i.e. lim (omega xi))= "<<SichelDispPars[1]<<endl ;
            //            } else {
            //                simulpars<<"Sichel dispersal distribution model (Chesson & Lee 2005)"<<endl;
            //                simulpars<<"gamma="<<SichelDispPars[0]<<endl ;
            //                simulpars<<"xi="<<SichelDispPars[1]<<endl ;
            //                simulpars<<"omega="<<SichelDispPars[2]<<endl ;
            //            }
            //          break;
            //}
        }
        param_summary_f << "Max dispersal distance : " << demo_param.Disp_dist_max[0] << ", " << demo_param.Disp_dist_max[1] << std::endl;
        param_summary_f << "Maximal number of immigrant between demes 2Nm : " << 2 * demo_param.Pop_size_per_node * demo_param.Proba_migr << std::endl;

        param_summary_f.close();
    }
    else
    {
        throw std::invalid_argument("( Unable to open " + path_to_file + ". I exit. )");
    }
}

void write_afterrun_param_settings_summary(std::string const &path_to_file, simu_param_c const &simu_param, samp_param_c const &samp_param, demo_param_c const &demo_param)
{
    std::ofstream param_summary_f(path_to_file, std::ios_base::app);
    if (param_summary_f.is_open())
    {
        auto &info_collect = singleton_c<info_collector_c>::instance();

        if (info_collect.Check_disp_distrib)
        {
            if (!simu_param.Migration_matrix)
            {
                param_summary_f << std::endl
                                << "Forward axial dispersal distribution (step, prob) : " << std::endl;
                int step = 0;
                for (auto const &info : info_collect.Fwd_axial_disp_distrib)
                {
                    ++step;
                    param_summary_f << step << "\t" << std::setprecision(15) << std::setw(17) << (info) << std::endl;
                }

                param_summary_f << std::endl
                                << "Forward axial dispersal mean (within dist_max) : " << info_collect.Fwd_axial_disp_mean << std::endl;
                param_summary_f << "Theoretical sigma2 without bounds : " << info_collect.Fwd_axial_disp_theo_sig << std::endl;
                param_summary_f << "Forward axial dispersal sigma2 (within dist_max) : " << info_collect.Fwd_axial_disp_sig << std::endl;
                param_summary_f << "Forward axial dispersal kurtosis (within dist_max) : " << info_collect.Fwd_axial_disp_kurt << std::endl;
                param_summary_f << "\nIf kurtosis appears wrong there may be many reasons" << std::endl;
                param_summary_f << "   It may be because the dispersal distribution is not" << std::endl;
                param_summary_f << "   computed with enough precision." << std::endl;
            }
            if (!(simu_param.Migration_matrix || simu_param.Nodesize_matrix))
            {
                info_collect.Ds2 = samp_param.Ploidy * demo_param.Pop_size_per_node * info_collect.Fwd_axial_disp_sig;
                param_summary_f << std::endl
                                << "Dsigma2 (within bounds) : " << info_collect.Ds2 << std::endl;
                param_summary_f << "Neighborhood size (within bounds) : " << (2.0 * info_collect.Ds2) << " (1D) or " << (2.0 * PI * info_collect.Ds2) << " (2D)." << std::endl;
                param_summary_f << "Expected slope (within bounds) : " << 1. / (2.0 * info_collect.Ds2) << " (1D) or " << 1.0 / (2.0 * PI * info_collect.Ds2) << " (2D)." << std::endl;

                param_summary_f << std::endl
                                << "Cumulative forward axial dispersal distribution";
                if (info_collect.Cumul_fwd_disp_distrib.size() == 2)
                {
                    param_summary_f << " on dim1";
                }
                else
                {
                    param_summary_f << " on both dimensions";
                }
                param_summary_f << " (move, prob) : " << std::endl;
                for (auto move = -demo_param.Disp_dist_max[0]; move <= demo_param.Disp_dist_max[0]; ++move)
                {
                    param_summary_f << move << "\t" << std::setprecision(15) << std::setw(17) << info_collect.Cumul_fwd_disp_distrib[0][move + demo_param.Disp_dist_max[0]] << std::endl;
                }
                if (info_collect.Cumul_fwd_disp_distrib.size() == 2)
                {
                    param_summary_f << std::endl << "Cumulative forward dispersal distribution on dim2 (move, prob) : " << std::endl;
                    for (auto move = -demo_param.Disp_dist_max[1]; move <= demo_param.Disp_dist_max[1]; ++move)
                    {
                        param_summary_f << move << "\t" << std::setprecision(15) << std::setw(17) << info_collect.Cumul_fwd_disp_distrib[1][move + demo_param.Disp_dist_max[1]] << std::endl;
                    }
                }
            }
            else
            {
                param_summary_f << std::endl
                                << "Neither Dsigma2, slope can be computed because a custom Nodesize matrix and/or a custom Migration matrix was considered." << std::endl;
            }
        }

        if (simu_param.Migration_matrix && info_collect.Print_migration_matrix_bool)
        {
            param_summary_f << std::endl
                            << "Forward migration matrix (forward migration rate m_LiCj from node i (xf, yf) -> node y (xt,yt) with i = 1 + (xf - 1)*(LatSizeX + 1) + (yj - 1) and j = 1 + (xt - 1)*(LatSizeX + 1) + (yt - 1) ) : " << std::endl;
            int num_from_node = 0;
            for (auto from_node_x = 0; from_node_x < demo_param.Lattice_size[0]; from_node_x++)
            {
                for (auto from_node_y = 0; from_node_y < demo_param.Lattice_size[1]; from_node_y++)
                {
                    int num_to_node = 0;
                    for (auto to_node_x = 0; to_node_x < demo_param.Lattice_size[0]; to_node_x++)
                    {
                        for (auto to_node_y = 0; to_node_y < demo_param.Lattice_size[1]; to_node_y++)
                        {
                            param_summary_f << std::setprecision(15) << std::setw(17) << demo_param.Migration_mat.at(num_from_node).at(num_to_node) << " ";
                            num_to_node++;
                        }
                    }
                    param_summary_f << std::endl;
                    num_from_node++;
                }
            }
            param_summary_f << std::endl
                            << "Backward migration matrix (backward migrant number Nm_LiCj i -> j with i and j = (x-1)*LatSizeX + (y-1) ) : " << std::endl;
            for (auto it1 = info_collect.Bcwd_disp_distrib_from_migration_matrix.begin(); it1 < info_collect.Bcwd_disp_distrib_from_migration_matrix.end(); ++it1)
            {
                for (auto it2 = it1->begin(); it2 < it1->end(); ++it2)
                {
                    param_summary_f << std::setprecision(15) << std::setw(17) << (*it2) << " ";
                }
                param_summary_f << std::endl;
            }
        }
        if (info_collect.Effective_disp)
        {
            size_t max_size;
            
            param_summary_f << std::endl;
            param_summary_f << "Empirical dispersal distribution ( step, freqX";
            if( info_collect.Emp_axial_disp_mean_over_rep.size() == 2 )
            {
                param_summary_f << ", freqY ";
            }
            param_summary_f << ") : " << std::endl;
            
            if (info_collect.Emp_axial_disp_mean_over_rep.size() == 2)
            {
                max_size = std::max(info_collect.Emp_axial_disp_mean_over_rep.at(0).size(),info_collect.Emp_axial_disp_mean_over_rep.at(1).size());
            }
            else{
                max_size = info_collect.Emp_axial_disp_mean_over_rep.at(0).size();
            }
            for (std::size_t i = 0; i < max_size; ++i)
            {
                            param_summary_f << static_cast<int>(i - (max_size - 1)/2);
                if (fabs(static_cast<int>(i - (max_size - 1)/2)) <= ((info_collect.Emp_axial_disp_mean_over_rep.at(0).size() - 1) / 2) )
                {
                    param_summary_f << "\t" << std::setprecision(15) << std::setw(17)
                                    << info_collect.Emp_axial_disp_mean_over_rep.at(0).at(i - ((max_size - 1)/2 - (info_collect.Emp_axial_disp_mean_over_rep.at(0).size() - 1)/2));
                }
                else
                {
                    param_summary_f << "\t  0.000000        ";
                }
                if (info_collect.Emp_axial_disp_mean_over_rep.size() == 2)
                {
                    if (fabs(static_cast<int>(i - (max_size - 1)/2)) <= ((info_collect.Emp_axial_disp_mean_over_rep.at(1).size() - 1) /  2))
                    {
                        param_summary_f << "\t" << std::setprecision(15) << std::setw(17) << info_collect.Emp_axial_disp_mean_over_rep.at(1).at(i - ((max_size - 1)/2 - (info_collect.Emp_axial_disp_mean_over_rep.at(1).size() - 1)/2));
                    }
                    else
                    {
                        param_summary_f << "\t  0.000000        ";
                    }
                }
                param_summary_f << std::endl;
            }
            
            size_t dim = 0;
            for (auto itr = info_collect.Emp_axial_disp_mean_over_rep.begin(); itr < info_collect.Emp_axial_disp_mean_over_rep.end(); ++(itr), ++dim)
            {
                double check_sum=0.0;
                int dist;
                size_t step=0;
                for (auto itr2 = itr->begin(); itr2 < itr->end(); ++(itr2), ++step)
                {
                    dist = static_cast<int>(step) - static_cast<int>(itr->size() - 1)/2;
                    info_collect.Mean_emp_axial_disp_mean_over_rep[dim] += (*itr2) * fabs(dist);
                    info_collect.Sig_emp_axial_disp_mean_over_rep[dim] += (*itr2) * dist * dist;
                    info_collect.Kurt_emp_axial_disp_mean_over_rep[dim] += (*itr2) * dist * dist * dist * dist;
                    check_sum += (*itr2);
                }
                if ( fabs(1.0 - check_sum) > DBL_EPSILON)
                {
                    throw std::logic_error("( During computation of the mean empirical dispersal distribution among runs, GSpace found that the sum of distribution terms is not 1.0 but " + std::to_string(check_sum) + "\nContact the developpers. I exit. )");
                }
                info_collect.Kurt_emp_axial_disp_mean_over_rep[dim] = info_collect.Kurt_emp_axial_disp_mean_over_rep[dim] / (info_collect.Sig_emp_axial_disp_mean_over_rep[dim] * info_collect.Sig_emp_axial_disp_mean_over_rep[dim]) - 3.0;
            }

            param_summary_f << std::endl
            << "Empirical dispersal mean (X [, Y]) : " << std::setprecision(15) << std::setw(17) << info_collect.Mean_emp_axial_disp_mean_over_rep[0];
            if (info_collect.Mean_emp_axial_disp_mean_over_rep.size() == 2)
            {
                param_summary_f << "\t" << std::setprecision(15) << std::setw(17) << info_collect.Mean_emp_axial_disp_mean_over_rep[1];
            }
            param_summary_f << std::endl;
            param_summary_f << "Empirical dispersal sigma2 (X [, Y]) : " << std::setprecision(15) << std::setw(17) << info_collect.Sig_emp_axial_disp_mean_over_rep[0];
            if (info_collect.Sig_emp_axial_disp_mean_over_rep.size() == 2)
            {
                param_summary_f << "\t" << std::setprecision(15) << std::setw(17) << info_collect.Sig_emp_axial_disp_mean_over_rep[1];
            }
            param_summary_f << std::endl;
            param_summary_f << "Empirical dispersal kurtosis (X [, Y]) : " << std::setprecision(15) << std::setw(17) << info_collect.Kurt_emp_axial_disp_mean_over_rep[0];
            if (info_collect.Kurt_emp_axial_disp_mean_over_rep.size() == 2)
            {
                param_summary_f << "\t" << std::setprecision(15) << std::setw(17) << info_collect.Kurt_emp_axial_disp_mean_over_rep[1];
            }
            param_summary_f << std::endl;
        }
        param_summary_f.close();
    }
    else
    {
        throw std::invalid_argument("( Unable to open " + path_to_file + ". I exit. )");
    }
}
