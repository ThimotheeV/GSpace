#include <chrono>

#include "simulator.hpp"

extern double const DBL_EPSILON;

void sample_simulator(output_stat_c &output_stat, int rep)
{
    //Global variable not to be modifie
    auto &rand_gen = singleton_c<rand_gen_c>::instance();
    auto const &simu_param = singleton_c<simu_param_c>::instance();
    auto const &demo_param = singleton_c<demo_param_c>::instance();
    auto const &samp_param = singleton_c<samp_param_c>::instance();
    auto const &recomb_param = singleton_c<recomb_param_c>::instance();
    auto const &muta_param = singleton_c<muta_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();
    auto &tskit = singleton_c<tskit_struct_c>::instance();

    std::chrono::time_point<std::chrono::high_resolution_clock> debut;
    if (info_collect.Clock)
    {
        debut = std::chrono::high_resolution_clock::now();
    }
    
    info_collect.rep = rep;

    std::vector<coa_table_c> coa_table_vec(samp_param.Chr_nbr);

    std::vector<int> next_max_node_nbr_vec;
    //ini migration lat and remap

    fwd_disp_distrib_c fwd_distrib(demo_param);
    //lattice_c lat(rand_gen, fwd_distrib, demo_param);

    lattice_c *lat;
    if (simu_param.Migration_matrix && simu_param.Nodesize_matrix)
    { //heterogeneous lattice with migration matrix but homogeneous subpopsizes
        lat = new lattice_c(demo_param.Nodesize_mat, rand_gen, demo_param.Migration_mat, demo_param);
    }
    else if (simu_param.Migration_matrix)
    { //heterogeneous lattice with migration matrix but homogeneous subpopsizes
        lat = new lattice_c(rand_gen, demo_param.Migration_mat, demo_param);
    }
    else if (simu_param.Nodesize_matrix)
    { //heterogeneous lattice with subpopsize matrix but homogeneous dispersal
        lat = new lattice_c(demo_param.Nodesize_mat, rand_gen, fwd_distrib, demo_param);
    }
    else
    { // homogeneous lattice
        lat = new lattice_c(rand_gen, fwd_distrib, demo_param);
    }
    extend_lattice_c rmap(*lat);

    //TODO faire une foncytion d'intialisation de effective disp.
    if (info_collect.Effective_disp)
    {
        size_t dim = lat->Lat_size.size();
        info_collect.Emp_axial_disp.clear();
        info_collect.Emp_cumul_axial_disp.clear();
        info_collect.Mean_emp_axial_disp.clear();
        info_collect.Sig_emp_axial_disp.clear();
        info_collect.Kurt_emp_axial_disp.clear();
        info_collect.Skew_emp_axial_disp.clear();
        info_collect.Emp_axial_disp.resize(dim);
        info_collect.Emp_cumul_axial_disp.resize(dim);
        info_collect.Mean_emp_axial_disp.resize(dim, 0.0);
        info_collect.Sig_emp_axial_disp.resize(dim, 0.0);
        info_collect.Kurt_emp_axial_disp.resize(dim, 0.0);
        info_collect.Skew_emp_axial_disp.resize(dim, 0.0);
        if (rep == 0)
        {
            info_collect.Emp_axial_disp_mean_over_rep.resize(dim);
            info_collect.Mean_emp_axial_disp_mean_over_rep.resize(dim, 0.0);
            info_collect.Sig_emp_axial_disp_mean_over_rep.resize(dim, 0.0);
            info_collect.Kurt_emp_axial_disp_mean_over_rep.resize(dim, 0.0);
            info_collect.Skew_emp_axial_disp_mean_over_rep.resize(dim, 0.0);
        }

        auto itr2 = info_collect.Emp_cumul_axial_disp.begin();
        auto itr3 = info_collect.Emp_axial_disp_mean_over_rep.begin();
        dim = 0;
        for (auto itr = info_collect.Emp_axial_disp.begin(); itr < info_collect.Emp_axial_disp.end(); ++(itr), ++dim, ++(itr2), ++(itr3))
        {
            itr->resize(2 * demo_param.Disp_dist_max[dim] + 1, 0.0);
            itr2->resize(2 * demo_param.Disp_dist_max[dim] + 1, 0);
            if (rep == 0)
            {
                itr3->resize(2 * demo_param.Disp_dist_max[dim] + 1, 0.0);
            }
        }
    }

    if (simu_param.Continuous_time_approxim)
    {
        // auto max_node_nbr_per_chr = ARG_simulator_continuous_time(lat, coa_table_vec, samp_param, recomb_param, rand_gen);
        next_max_node_nbr_vec = ARG_simulator_continuous_time(*lat, coa_table_vec, samp_param, recomb_param, rand_gen);
    }
    else
    {
        next_max_node_nbr_vec = ARG_simulator_gen_by_gen(*lat, rmap, coa_table_vec, rand_gen);
    }

    delete lat;

    if (info_collect.Clock)
    {
        auto fin = std::chrono::high_resolution_clock::now();
        info_collect.time_simulation += (std::chrono::duration_cast<std::chrono::nanoseconds>(fin - debut).count() * 0.000000001);
    }

    if (next_max_node_nbr_vec.size() != coa_table_vec.size())
    {
        throw std::logic_error("( Something went horribly wrong in this run. Safety first : disinfect your computer with holy water then contact the IT slave who is responsible for this pile of code... I exit. )");
    }

    // if (info_collect.Debug_m)
    // {
    //     info_collect.Debug_mut.resize(samp_param.Sequence_length);
    // }

    if ((muta_param.Unscaled_mut_rate_mu != 0) || (muta_param.Scaled_mut_rate_theta != 0))
    {
        //chr, indiv, mut{pos, state}
        std::vector<std::vector<std::vector<std::pair<int, int>>>> sample_mutated_state_chr;
        sample_mutated_state_chr.reserve(samp_param.Chr_nbr);
        if (info_collect.MRCA_record)
        {
            info_collect.MRCA_time.clear();
        }

        for (int chr = 0; chr < samp_param.Chr_nbr; ++chr)
        {
            if (tskit.Tskit_output)
            {
                tskit.add_new_chrom(coa_table_vec.at(chr));
            }
            sample_mutated_state_chr.push_back(apply_mut_to_sample(coa_table_vec.at(chr), next_max_node_nbr_vec.at(chr), muta_param.Ancestry_seq.at(chr), muta_param, samp_param, demo_param));
            ++tskit.Chr_num;
        }

        if (info_collect.Stats)
        {
            if (rep == 0)
            {
                output_stat.Data_plane_vec = data_plane_vec_c(sample_mutated_state_chr, samp_param, muta_param.Ancestry_seq);
            }
            else
            {
                output_stat.Data_plane_vec.update_data_plane_vec(sample_mutated_state_chr, muta_param.Ancestry_seq, samp_param);
            }

            if (info_collect.Prob_id_1_loc_Qr)
            {
                output_stat.Prob_id_1_loc_Qr_res.calc_Qr(output_stat.Data_plane_vec, rep);
            }

            if (info_collect.Prob_id_1_loc_Qwi_wd_bd)
            {
                if (samp_param.Ploidy > 1)
                {
                    output_stat.Prob_id_1_loc_Qwi_wd_bd_res.calc_Qwi(output_stat.Data_plane_vec, rep);
                }
                output_stat.Prob_id_1_loc_Qwi_wd_bd_res.calc_Qwd_bd(output_stat.Data_plane_vec, rep);
            }

            if (info_collect.Prob_id_1_2_loc)
            {
                output_stat.Prob_id_1_2_loc_res.cumul_prob_id_1_2_loc(output_stat.Data_plane_vec, rep);
            }

            if (info_collect.Effective_disp)
            {
                auto itr2 = info_collect.Emp_axial_disp.begin();
                auto itr3 = info_collect.Emp_axial_disp_mean_over_rep.begin();
                size_t dim = 0;
                for (auto itr1 = info_collect.Emp_cumul_axial_disp.begin(); itr1 < info_collect.Emp_cumul_axial_disp.end(); ++(itr1), ++(itr2), ++(itr3))
                {

                    int total_disp_draw_nbr = 0;
                    for (auto itr12 = itr1->begin(); itr12 < itr1->end(); ++(itr12))
                    {
                        total_disp_draw_nbr += *(itr12);
                    }
                    double check_sum = 0.0;
                    int dist;
                    size_t step = 0;
                    auto itr22 = itr2->begin();
                    auto itr32 = itr3->begin();
                    for (auto itr12 = itr1->begin(); itr12 < itr1->end(); ++(itr12), ++(itr22), ++(itr32), ++step)
                    {
                        dist = static_cast<int>(step) - static_cast<int>(itr1->size() - 1) / 2;
                        (*itr22) = 1. * (*itr12) / (1. * total_disp_draw_nbr);
                        info_collect.Mean_emp_axial_disp[dim] += (*itr22) * fabs(dist);
                        info_collect.Sig_emp_axial_disp[dim] += (*itr22) * dist * dist;
                        info_collect.Kurt_emp_axial_disp[dim] += (*itr22) * dist * dist * dist * dist;
                        info_collect.Skew_emp_axial_disp[dim] += (*itr22) * dist * dist * dist;
                        check_sum += (*itr22);
                        (*itr32) += (*itr22) / simu_param.Repetition_nbr;
                        //                        std::cout << "dim=" << dim << " *itr22 " << *itr22 << " = *itr12 " << *itr12 << " / total_disp_draw_nbr " << total_disp_draw_nbr << std::endl;
                        //                        std::cout << "      *itr32 " << *itr32 << " += *itr22 " << *itr22 << " / simu_param.Repetition_nbr " << simu_param.Repetition_nbr << std::endl;
                    }
                    if (fabs(1.0 - check_sum) > DBL_EPSILON)
                    {
                        throw std::logic_error("( During computation of empirical dispersal distribution, GSpace found that the sum of distribution terms is not 1.0 but " + std::to_string(check_sum) + "\nContact the developpers. I exit. )");
                    }
                    info_collect.Kurt_emp_axial_disp[dim] = info_collect.Kurt_emp_axial_disp[dim] / (info_collect.Sig_emp_axial_disp[dim] * info_collect.Sig_emp_axial_disp[dim]) - 3.0;
                    info_collect.Skew_emp_axial_disp[dim] = info_collect.Skew_emp_axial_disp[dim] / pow(info_collect.Sig_emp_axial_disp[dim], double(3 / 2));

                    ++dim;
                }
            }
            if (output_stat.Stat_out)
            {
                output_stat_files(simu_param, info_collect, output_stat, rep, samp_param.Ploidy);
            }

            //Calcul MRCA_time and 2_Gen_coa_time
            if (info_collect.MRCA_record)
            {
                info_collect.MRCA_record_cumul_mean_var = calc_cumul_mean_var(info_collect.MRCA_record_cumul_mean_var, info_collect.MRCA_time, rep);
                info_collect.MRCA_time.clear();
            }
        }
        if (simu_param.Genepop_output)
        {
            map_string_mod_mut_name_c map_string_mod_mut_name;
            auto genepop_filename = simu_param.Generic_data_filename;
            genepop_filename += "_GP_" + std::to_string(rep + 1) + simu_param.Data_file_extension;
            genepop_output(genepop_filename, samp_param.Sample_coord_vec, muta_param.Ancestry_seq, sample_mutated_state_chr, map_string_mod_mut_name.ReverseMap[muta_param.Mod_mut_name], samp_param.Ploidy);
        }
        if (simu_param.VCF_output)
        {
            auto vcf_filename = simu_param.Generic_data_filename;
            vcf_filename += "_VCF_" + std::to_string(rep + 1) + ".vcf";
            vcf_output(vcf_filename, samp_param.Sample_coord_vec, muta_param.Ancestry_seq, sample_mutated_state_chr, samp_param.Ploidy);
        }
        if (simu_param.Fasta_output)
        {
            map_string_mod_mut_name_c map_string_mod_mut_name;
            auto fasta_filename = simu_param.Generic_data_filename;
            fasta_filename += "_Fasta_" + std::to_string(rep + 1) + ".fa";
            fasta_output(fasta_filename, simu_param.Generic_data_filename, samp_param.Sample_coord_vec, muta_param.Ancestry_seq, sample_mutated_state_chr, samp_param.Ploidy, rep);
        }
        if (simu_param.Phylip_output)
        {
            for (std::size_t chrom = 0; chrom < sample_mutated_state_chr.size(); chrom++)
            {
                auto phylip_filename = simu_param.Generic_data_filename;
                phylip_filename += "_Phylip_" + std::to_string(rep + 1) + "_chrom" + std::to_string(chrom + 1) + ".phy";
                phylip_output(phylip_filename, simu_param.Generic_data_filename, samp_param.Sample_coord_vec, muta_param.Ancestry_seq, sample_mutated_state_chr.at(chrom), samp_param.Ploidy, chrom);
            }
        }
        if (simu_param.Seq_char_output)
        {
            map_string_mod_mut_name_c map_string_mod_mut_name;
            auto seq_char_filename = simu_param.Generic_data_filename;
            seq_char_filename += "_seq_char_" + std::to_string(rep + 1) + ".txt";
            seq_char_output(seq_char_filename, samp_param.Sample_coord_vec, sample_mutated_state_chr, samp_param.Ploidy);
        }
        if (simu_param.Coordinates_output)
        {
            auto coord_filename = simu_param.Generic_data_filename;
            coord_filename += "_coord_" + std::to_string(rep + 1) + ".txt";
            coord_output(coord_filename, samp_param.Sample_coord_vec);
        }
        if (simu_param.Tskit_output)
        {
            auto coord_filename = simu_param.Generic_data_filename;
            coord_filename += "_tskit_" + std::to_string(rep + 1) + ".tree";
            tskit.sort_and_output(coord_filename);
        }
    }
}

std::vector<int> ARG_simulator_continuous_time(lattice_c &lat, std::vector<coa_table_c> &coa_table_vec, samp_param_c const &samp_param, recomb_param_c const &recomb_param, rand_gen_c &rand_gen)
{

    auto &info_collect = singleton_c<info_collector_c>::instance();

    std::vector<int> next_node_nbr_per_chr;
    next_node_nbr_per_chr.reserve(samp_param.Chr_nbr);

    // std::cout << "recomb_param.Scaled_recomb_rate_rho : " << recomb_param.Scaled_recomb_rate_rho << std::endl;

    for (auto &coa_table_chr : coa_table_vec)
    {
        //empty constructor
        struct_arg_c algo_hudson;
        indiv_stock_c indiv_vec(1);
        //TMRCA tends to 2 for large sample sizes
        indiv_vec.ini(lat, samp_param);
        algo_hudson.ini(indiv_vec, samp_param, 0);

        int count = 0;

        if (samp_param.Ploidy == 2)
        {
            segregate_chr_btw_indiv(lat, indiv_vec, 1, rand_gen);
        }

        while (indiv_vec.size() != 0)
        {
            count += 1;
            //std::cout << "Nombre d'individus en recomb: " << indiv_vec.size() << std::endl;

            std::chrono::time_point<std::chrono::high_resolution_clock> debut;
            if (info_collect.Clock)
            {
                debut = std::chrono::high_resolution_clock::now();
            }

            if (algo_hudson.choose_event(recomb_param.Scaled_recomb_rate_rho, indiv_vec, rand_gen))
            /***************************************************/
            /*                Recombination                    */
            /***************************************************/
            {
                //Choose a brkpt with is abbsolute number (between 0 and sum(brkpt/segment))
                std::uniform_int_distribution<long int> uniform_distribution(1, algo_hudson.L_cumul_nbr_brkpt_per_seg.get_tot_cumul_freq());
                long int h_choosen_brkpt = uniform_distribution(rand_gen.Seed_gen);

                int choosen_seg_index = algo_hudson.L_cumul_nbr_brkpt_per_seg.seg_index(h_choosen_brkpt);
                auto y_choosen_seg = algo_hudson.Pool_seg.get_obj_from_index(choosen_seg_index);
                //Transform absolute position h in all seg on relative position k in y_choosen_seg
                auto k_choose_brkpt = y_choosen_seg->R_right_brkpt - algo_hudson.L_cumul_nbr_brkpt_per_seg.get_cumul_freq_at_index(choosen_seg_index) + h_choosen_brkpt - 1;

                std::array<seg *, 2> new_seg = apply_first_recomb(algo_hudson.Pool_seg, algo_hudson.L_cumul_nbr_brkpt_per_seg, y_choosen_seg, k_choose_brkpt);
                //WARNING : Create a new indiv, no need to update the old because first part of the seg is still in him
                indiv_c *new_indiv = indiv_vec.new_indiv(lat.node({0, 0}));
                new_indiv->update_indiv(0, new_seg.at(1), nullptr);
                lat.add_indiv(new_indiv, {0, 0});

                if (info_collect.Clock)
                {
                    auto fin = std::chrono::high_resolution_clock::now();
                    info_collect.time_recomb += (std::chrono::duration_cast<std::chrono::nanoseconds>(fin - debut).count() * 0.000000001);
                }
            }
            else
            /***************************************************/
            /*                 Coalescence                     */
            /***************************************************/
            {
                int num_indiv_1, num_indiv_2;
                std::uniform_int_distribution<long int> uniform_distribution(0, indiv_vec.size() - 1);
                num_indiv_1 = uniform_distribution(rand_gen.Seed_gen);
                num_indiv_2 = uniform_distribution(rand_gen.Seed_gen);

                //Choose two differents seg
                while (num_indiv_1 == num_indiv_2)
                {
                    num_indiv_2 = uniform_distribution(rand_gen.Seed_gen);
                }

                auto indiv_1 = indiv_vec[num_indiv_1];
                auto indiv_2 = indiv_vec[num_indiv_2];

                std::array<seg *, 2> two_lign = {indiv_1->get_chr(0, 0), indiv_2->get_chr(0, 0)};

                //Destructor for indiv 2 and clean indiv_1 for reuse
                indiv_vec.clean_indiv_at_chr(indiv_1->Ident, 0);
                indiv_vec.erase(indiv_2->Ident);

                //the -1 is a flag for nbr of gen (not useable in hudson algo)
                seg *new_seg = apply_coa_two_lineage(coa_table_chr, algo_hudson.Pool_seg, algo_hudson.L_cumul_nbr_brkpt_per_seg, algo_hudson.S_intersection_count, algo_hudson.W_next_node_nbr, algo_hudson.T_time, -1, two_lign.at(0), two_lign.at(1));
                if (new_seg != nullptr)
                {
                    indiv_1->update_indiv(0, new_seg, nullptr);
                }
                else
                {
                    indiv_vec.erase(indiv_1->Ident);
                }
                if (info_collect.Clock)
                {
                    auto fin = std::chrono::high_resolution_clock::now();
                    info_collect.time_coa += (std::chrono::duration_cast<std::chrono::nanoseconds>(fin - debut).count() * 0.000000001);
                }
            }
            /***************************************************/
            /*                   Migration                     */
            /***************************************************/
        }
        next_node_nbr_per_chr.push_back(algo_hudson.W_next_node_nbr);
    }
    return next_node_nbr_per_chr;
}

std::vector<int> ARG_simulator_gen_by_gen(lattice_c &lat, extend_lattice_c &rmap, std::vector<coa_table_c> &coa_table_vec, rand_gen_c &rand_gen)
{
    auto const &samp_param = singleton_c<samp_param_c>::instance();
    auto const &recomb_param = singleton_c<recomb_param_c>::instance();
    auto &info_collect = singleton_c<info_collector_c>::instance();
    //Empty Constructor
    std::vector<int> next_node_nbr_per_chr;
    next_node_nbr_per_chr.reserve(samp_param.Chr_nbr);

    std::vector<struct_arg_c> struct_gen_by_gen_vec(samp_param.Chr_nbr);
    indiv_stock_c indiv_vec(samp_param.Chr_nbr);
    indiv_vec.ini(lat, samp_param);

    for (int chr_index = 0; chr_index < samp_param.Chr_nbr; ++chr_index)
    {
        struct_gen_by_gen_vec.at(chr_index).ini(indiv_vec, samp_param, chr_index);
    }

    int gen = 0;

    while (indiv_vec.size() != 0)
    {
        gen += 1;

        std::chrono::time_point<std::chrono::high_resolution_clock> debut;
        if (info_collect.Clock)
        {
            debut = std::chrono::high_resolution_clock::now();
        }

        //std::cout << indiv_vec.size() << std::endl;

        /***************************************************/
        /*                   Migration                     */
        /***************************************************/
        //Clean node before migration to avoid indiv_c double in node_c
        //Position redondante avec la position des indiv dans le lattice donc on peut nettoyer sans risque
        //TODO : Rentrer les vclean dans le if et enlever le else
        for (auto node : lat.Nodes_with_lineage)
        {
            node->Indivs_in_pop.clear();
        }
        //Clean to avoid empty node to stay here after cleaning of indiv_c
        lat.Nodes_with_lineage.clear();

        if (lat.Lattice.size() > 1)
        {
            migration(lat, &rmap, indiv_vec);
        }
        else
        { //Special case, Wright-Fisher model
            auto node = lat.node({0, 0});
            for (auto &indiv : indiv_vec)
            {
                add_indiv_at_node(indiv, node);
            }
        }

        if (info_collect.Clock)
        {
            auto fin = std::chrono::high_resolution_clock::now();
            info_collect.time_mig += (std::chrono::duration_cast<std::chrono::nanoseconds>(fin - debut).count() * 0.000000001);
        }

        if (samp_param.Ploidy == 2)
        {
            segregate_chr_btw_indiv(lat, indiv_vec, samp_param.Chr_nbr, rand_gen);
        }
        /***************************************************/
        /*                Recombination                    */
        /***************************************************/

        if ((recomb_param.Unscaled_recomb_rate > 0))
        {
            for (int chr_index = 0; chr_index < samp_param.Chr_nbr; ++chr_index)
            {
                recomb_gen_by_gen(struct_gen_by_gen_vec.at(chr_index), recomb_param.Unscaled_recomb_rate, chr_index, rand_gen);
            }
        }

        if (samp_param.Ploidy == 1)
        {
            segregate_chr_btw_indiv(lat, indiv_vec, samp_param.Chr_nbr, rand_gen);
        }

        if (info_collect.Clock)
        {
            auto fin = std::chrono::high_resolution_clock::now();
            info_collect.time_recomb += (std::chrono::duration_cast<std::chrono::nanoseconds>(fin - debut).count() * 0.000000001);
        }

        /***************************************************/
        /*                 Coalescence                     */
        /***************************************************/

        coa_gen_by_gen(indiv_vec, struct_gen_by_gen_vec, lat, coa_table_vec, gen, samp_param, rand_gen);

        if (info_collect.Clock)
        {
            auto fin = std::chrono::high_resolution_clock::now();
            info_collect.time_coa += (std::chrono::duration_cast<std::chrono::nanoseconds>(fin - debut).count() * 0.000000001);
        }
    }
    for (auto const &struct_gen_by_gen : struct_gen_by_gen_vec)
    {
        next_node_nbr_per_chr.push_back(struct_gen_by_gen.W_next_node_nbr);
    }

    return next_node_nbr_per_chr;
}

void recomb_gen_by_gen(struct_arg_c &struct_gen_by_gen, double unscaled_recomb_rate, int chr_index, rand_gen_c &rand_gen)
{
    //L_cumul_nbr_brkpt_per_seg = Fenwick tree
    long int cumul_size = struct_gen_by_gen.L_cumul_nbr_brkpt_per_seg.get_tot_cumul_freq();
    int nbr_recomb_event;
    if (cumul_size > 10000)
    {
        std::poisson_distribution<long int> distri(cumul_size * unscaled_recomb_rate);
        nbr_recomb_event = distri(rand_gen.Seed_gen);
    }
    else
    {
        std::binomial_distribution<int> distri(cumul_size, unscaled_recomb_rate);
        nbr_recomb_event = distri(rand_gen.Seed_gen);
    }

    //For each recomb event we draw a new brkpt
    while ((nbr_recomb_event > 0) && (struct_gen_by_gen.L_cumul_nbr_brkpt_per_seg.get_tot_cumul_freq() != 0))
    {
        /* PREPARATION */

        //need to use the update fenwick tree at each round (each recomb event withing a seg reduce the total number of brkpt of 2)
        std::uniform_int_distribution<long int> uni_distri(1, struct_gen_by_gen.L_cumul_nbr_brkpt_per_seg.get_tot_cumul_freq());

        long int h_choosen_brkpt = uni_distri(rand_gen.Seed_gen);
        //Index of seg who carry h_choosen_brkpt
        int choosen_seg_index = struct_gen_by_gen.L_cumul_nbr_brkpt_per_seg.seg_index(h_choosen_brkpt);
        //Find seg with index
        auto y_choosen_seg = struct_gen_by_gen.Pool_seg.get_obj_from_index(choosen_seg_index);
        //Transform absolute position h in all seg on relative position k in y_choosen_seg
        auto k_choose_brkpt = y_choosen_seg->R_right_brkpt - struct_gen_by_gen.L_cumul_nbr_brkpt_per_seg.get_cumul_freq_at_index(choosen_seg_index) + h_choosen_brkpt - 1;

        //Find homolog_seg : find the indiv who carry y_choosen_seg (carried by the first seg in the chromosome)
        seg *homolog_seg = nullptr;
        seg *first_y_choosen_seg = find_first_seg(y_choosen_seg);
        auto indiv = first_y_choosen_seg->get_indiv();

        //If first_y_choosen_seg is in the first place, homolog is in the secund ->nullptr if no homolog
        if (indiv->get_chr(chr_index, 0) == first_y_choosen_seg)
        {
            homolog_seg = indiv->get_chr(chr_index, 1);
        }
        else
        {
            homolog_seg = indiv->get_chr(chr_index, 0);
        }

        /* ACTION */

        //If y_choosen_seg has not yet recomb this gen : create a homolog chr
        if ((indiv->get_chr(chr_index, 0) == nullptr) || (indiv->get_chr(chr_index, 1) == nullptr))
        {
            auto new_seg = apply_first_recomb(struct_gen_by_gen.Pool_seg, struct_gen_by_gen.L_cumul_nbr_brkpt_per_seg, y_choosen_seg, k_choose_brkpt);
            indiv->update_homolog_seg(chr_index, new_seg.at(1));
        }
        else
        //If chr1 has an homolog, we swap chain seg after the choosen breakpoint
        //For that we apply a recomb event at the same place in the homolog chr2. We coa the new_seg1 with chr2 and the new_seg2 with new_seg1 !
        {
            crossing_over(struct_gen_by_gen, y_choosen_seg, k_choose_brkpt, first_y_choosen_seg, homolog_seg, indiv, chr_index);
        }

        --nbr_recomb_event;
    }
}

void crossing_over(struct_arg_c &struct_gen_by_gen, seg *y_choosen_seg, int choose_brkpt_chr1, seg *first_seg_chr1, seg *first_seg_chr2, indiv_c *indiv, int chr_index)
{
    //Create a new_seg1 who are not attach yet to homolog chr2
    auto prev_and_recomb_seg_chr1 = apply_first_recomb(struct_gen_by_gen.Pool_seg, struct_gen_by_gen.L_cumul_nbr_brkpt_per_seg, y_choosen_seg, choose_brkpt_chr1);

    //Find the segment who carry the choosen breakpoint in the homolog chr (chr2)
    auto seg_browser_chr2 = first_seg_chr2;
    while (seg_browser_chr2 != nullptr)
    {
        if (seg_browser_chr2->R_right_brkpt > choose_brkpt_chr1)
        {
            break;
        }
        seg_browser_chr2 = seg_browser_chr2->Next_seg;
    }

    std::array<seg *, 2> prev_and_recomb_seg_chr2 = {first_seg_chr2, nullptr};
    //Found the seg who carry or the seg who was after the choosen breakpoint
    if (seg_browser_chr2 != nullptr)
    {
        prev_and_recomb_seg_chr2 = apply_first_recomb(struct_gen_by_gen.Pool_seg, struct_gen_by_gen.L_cumul_nbr_brkpt_per_seg, seg_browser_chr2, choose_brkpt_chr1);

        //If the choosen brkp is before the first seg of chr1 (possible if chr1 is a homolog seg of the original seg in indiv)
        if (prev_and_recomb_seg_chr1.at(0) != nullptr)
        {
            simple_coa(struct_gen_by_gen, prev_and_recomb_seg_chr1.at(0), prev_and_recomb_seg_chr2.at(1));
        }
        else
        {
            first_seg_chr1 = prev_and_recomb_seg_chr2.at(1);
        }
    }
    //if new seg chr2 was the first segment of chr2 (prev_and_recomb_seg_chr2.at(0) == nullptr because in prev apply_recomb seg_browser_chr2 haven't prev seg) : new seg chr1 become the new homolog
    if (prev_and_recomb_seg_chr2.at(0) != nullptr)
    {
        simple_coa(struct_gen_by_gen, prev_and_recomb_seg_chr2.at(0), prev_and_recomb_seg_chr1.at(1));
    }
    else
    {
        first_seg_chr2 = prev_and_recomb_seg_chr1.at(1);
    }

    indiv->update_indiv(chr_index, first_seg_chr1, first_seg_chr2);
}

//Come from Kellerher
std::array<seg *, 2> apply_first_recomb(obj_mem_manag_c<seg> &pool_seg, fenwick_tree_c &cumul_nbr_brkpt_per_seg, seg *y_choosen_seg, int k_choose_brkpt)
{
    alg_recomb_c algo_recomb;

    if (!(algo_recomb.ini(y_choosen_seg, k_choose_brkpt)))
    {
        algo_recomb.break_between_segment();
    }
    else
    {
        algo_recomb.break_within_segment(pool_seg, cumul_nbr_brkpt_per_seg);
    }
    //return new seg ptr
    return algo_recomb.update_population(pool_seg, cumul_nbr_brkpt_per_seg);
}

void coa_gen_by_gen(indiv_stock_c &indiv_vec, std::vector<struct_arg_c> &struct_gen_by_gen_vec, lattice_c &lat, std::vector<coa_table_c> &coa_table_vec, int gen, samp_param_c const &samp_param, rand_gen_c &rand_gen)
{
    std::vector<std::size_t> coa_table_size(coa_table_vec.size());
    //to know where multicoa arise
    for (std::size_t chr = 0; chr < coa_table_vec.size(); ++chr)
    {
        coa_table_size.at(chr) = coa_table_vec.at(chr).size();
    }
    //Using lattice with lineage for visit only lattice node with potential coalescence event
    for (auto &node_lat : lat.Nodes_with_lineage)
    {
        std::uniform_int_distribution<int> uni_distrib(1, node_lat->Subpop_size);
        //Vector of indiv
        auto &indiv_in_pop = node_lat->Indivs_in_pop;
        //Each indiv belong to an ancester (pull in uni_distrib)
        for (auto &indiv : indiv_in_pop)
        {
            std::get<0>(indiv) = uni_distrib(rand_gen.Seed_gen);
        }
        //Sort by random nbr
        std::sort(indiv_in_pop.begin(), indiv_in_pop.end(),
                  [](auto const &tuple1, auto const &tuple2) {
                      return (std::get<0>(tuple1) < std::get<0>(tuple2));
                  });

        for (auto indiv_itr = indiv_in_pop.begin(); indiv_itr != indiv_in_pop.end(); ++indiv_itr)
        {
            //take  first indiv
            indiv_c *cell_lig1 = std::get<1>(*indiv_itr);

            //If two or more cell_lig have the same ancester, they coa
            while ((indiv_itr + 1 != indiv_in_pop.end()) && (std::get<0>(*indiv_itr) == std::get<0>(*(indiv_itr + 1))))
            {
                ++indiv_itr;
                //Choose two cell_lig
                indiv_c *cell_lig2 = std::get<1>(*indiv_itr);
                //ditrib_lign_btw_2_chr => intra_indiv_coa
                if (samp_param.Ploidy == 1)
                {
                    cell_lig1 = ditrib_and_coa_lign_btw_1_chr(coa_table_vec, indiv_vec, struct_gen_by_gen_vec, gen, cell_lig1, cell_lig2);
                }
                else
                {
                    if (samp_param.Ploidy == 2)
                    {
                        cell_lig1 = ditrib_and_coa_lign_btw_2_chr(coa_table_vec, indiv_vec, struct_gen_by_gen_vec, gen, cell_lig1, cell_lig2);
                    }
                }
            }
        }
    }
    //Search for multi_coalescence
    for (std::size_t chr = 0; chr < coa_table_vec.size(); ++chr)
    {
        if (coa_table_size.at(chr) != coa_table_vec.at(chr).size())
        {
            coa_table_vec.at(chr).group_multi_coa(coa_table_size.at(chr));
        }
    }
}

//Haploid coalescence (reverse mitose)
indiv_c *ditrib_and_coa_lign_btw_1_chr(std::vector<coa_table_c> &coa_table_vec, indiv_stock_c &indiv_vec, std::vector<struct_arg_c> &struct_arg_vec, int gen, indiv_c *indiv_1, indiv_c *indiv_2)
{
    //in the case when a partial mrca take place at the step before
    if (indiv_1 == nullptr)
    {
        return indiv_2;
    }

    seg *new_seg;

    bool empty_indiv = true;

    for (std::size_t chr_index = 0; chr_index < struct_arg_vec.size(); ++chr_index)
    {
        std::array<seg *, 2> two_lign = {indiv_1->get_chr(chr_index, 0), indiv_2->get_chr(chr_index, 0)};

        //Destructor for indiv 2 and clean indiv_1 for reuse
        indiv_vec.clean_indiv_at_chr(indiv_1->Ident, chr_index);
        indiv_vec.clean_indiv_at_chr(indiv_2->Ident, chr_index);

        auto &struct_arg = struct_arg_vec.at(chr_index);
        //Cell1 and cell2 coa
        new_seg = apply_coa_two_lineage(coa_table_vec.at(chr_index), struct_arg.Pool_seg, struct_arg.L_cumul_nbr_brkpt_per_seg, struct_arg.S_intersection_count, struct_arg.W_next_node_nbr, struct_arg.T_time, gen, two_lign.at(0), two_lign.at(1));

        //Do nothing because reach a local mrca
        if (new_seg != nullptr)
        {
            empty_indiv = false;
        }

        indiv_1->update_indiv(chr_index, new_seg, nullptr);
    }
    //WARNING : Indiv2 is still in node lat but will be clean before mig (to difficult to retreive indiv2 in node_lat to clean node_lat here)
    indiv_vec.erase(indiv_2->Ident);

    if (empty_indiv)
    {
        indiv_vec.erase(indiv_1->Ident);
        return nullptr;
    }
    else
    {
        return indiv_1;
    }
}

//Diploid coalescence (reverse meiose)
indiv_c *ditrib_and_coa_lign_btw_2_chr(std::vector<coa_table_c> &coa_table_vec, indiv_stock_c &indiv_vec, std::vector<struct_arg_c> &struct_arg_vec, int gen, indiv_c *indiv_1, indiv_c *indiv_2)
{
    //in the case when a partial mrca take place at the step before
    if (indiv_1 == nullptr)
    {
        return indiv_2;
    }

    seg *new_seg1, *new_seg2;
    auto &rand_gen = singleton_c<rand_gen_c>::instance();

    bool empty_indiv = true;

    for (std::size_t chr_index = 0; chr_index < struct_arg_vec.size(); ++chr_index)
    {
        std::array<seg *, 4> four_lign = {indiv_1->get_chr(chr_index, 0), indiv_1->get_chr(chr_index, 1), indiv_2->get_chr(chr_index, 0), indiv_2->get_chr(chr_index, 1)};

        //Destructor for indiv 2 and clean indiv_1 for reuse
        indiv_vec.clean_indiv_at_chr(indiv_1->Ident, chr_index);
        indiv_vec.clean_indiv_at_chr(indiv_2->Ident, chr_index);

        //Without autofecondation for the moment
        four_lign = random_coa_process_btw_seg(four_lign, rand_gen);
        auto &struct_arg = struct_arg_vec.at(chr_index);
        //Cell1 and cell2 have one chr in placeholder 1 they coa
        if (four_lign[1] != nullptr)
        {
            new_seg1 = apply_coa_two_lineage(coa_table_vec.at(chr_index), struct_arg.Pool_seg, struct_arg.L_cumul_nbr_brkpt_per_seg, struct_arg.S_intersection_count, struct_arg.W_next_node_nbr, struct_arg.T_time, gen, four_lign[0], four_lign[1]);
        }
        else
        {
            new_seg1 = four_lign[0];
        }

        //Cell1 and cell2 have one chr in placeholder 1 they coa
        if ((four_lign[2] != nullptr) && (four_lign[3] != nullptr))
        {
            new_seg2 = apply_coa_two_lineage(coa_table_vec.at(chr_index), struct_arg.Pool_seg, struct_arg.L_cumul_nbr_brkpt_per_seg, struct_arg.S_intersection_count, struct_arg.W_next_node_nbr, struct_arg.T_time, gen, four_lign[2], four_lign[3]);
        }
        else if (four_lign[2] != nullptr)
        {
            new_seg2 = four_lign[2];
        }
        else if (four_lign[3] != nullptr)
        {
            new_seg2 = four_lign[3];
        }
        else
        {
            new_seg2 = nullptr;
        }

        if (new_seg1 != nullptr || new_seg2 != nullptr)
        {
            empty_indiv = false;
        }

        indiv_1->update_indiv(chr_index, new_seg1, new_seg2);
    }

    indiv_vec.erase(indiv_2->Ident);

    //Do nothing because reach a local mrca
    if (empty_indiv)
    {
        indiv_vec.erase(indiv_1->Ident);
        return nullptr;
    }
    else
    {
        if (indiv_1->empty())
        {
            throw std::logic_error("( In ditrib_and_coa_lign_btw_2_chr() : indiv is return in a inconsistant state. Contact the developpers. I exit. )");
        }
        return indiv_1;
    }
}

seg *apply_coa_two_lineage(coa_table_c &coa_table, obj_mem_manag_c<seg> &pool_seg, fenwick_tree_c &cumul_nbr_brkpt_per_seg, std::map<int, int> &intersection_count, int &next_max_node_nbr, double &t_time, int gen, seg *first_seg_lineage_1, seg *first_seg_lineage_2)
{
    alg_coa_c apply_coa_two_lineage;

    //This implementation come from 2016 Kelleher paper "Efficient Coalescence Simulation and genealogical Analysis for Large Sample Sizes"
    apply_coa_two_lineage.choose_ancestors(first_seg_lineage_1, first_seg_lineage_2);

    while ((apply_coa_two_lineage.X_prev_seg != nullptr) || (apply_coa_two_lineage.Y_choosen_seg != nullptr))
    {
        apply_coa_two_lineage.Alpha_new_seg = nullptr;
        if ((apply_coa_two_lineage.X_prev_seg != nullptr) && (apply_coa_two_lineage.Y_choosen_seg != nullptr))
        {
            if (apply_coa_two_lineage.choose_case(pool_seg))
            {
                if (!(apply_coa_two_lineage.coa(intersection_count, next_max_node_nbr)))
                {
                    apply_coa_two_lineage.seg_mrca(intersection_count);
                }
                else
                {
                    apply_coa_two_lineage.decrement_overlaps(pool_seg, intersection_count);
                }
                apply_coa_two_lineage.update_segs(pool_seg, cumul_nbr_brkpt_per_seg, coa_table, t_time, gen);
            }
        }
        else
        {
            if (apply_coa_two_lineage.X_prev_seg != nullptr)
            {
                apply_coa_two_lineage.Alpha_new_seg = apply_coa_two_lineage.X_prev_seg;
                apply_coa_two_lineage.X_prev_seg = nullptr;
            }

            if (apply_coa_two_lineage.Y_choosen_seg != nullptr)
            {
                apply_coa_two_lineage.Alpha_new_seg = apply_coa_two_lineage.Y_choosen_seg;
                apply_coa_two_lineage.Y_choosen_seg = nullptr;
            }
        }
        apply_coa_two_lineage.update_brkpts(pool_seg, cumul_nbr_brkpt_per_seg);
    }
    if (apply_coa_two_lineage.Bool_coa)
    {
        apply_coa_two_lineage.prune_tree(intersection_count);
    }
    if (apply_coa_two_lineage.Bool_defrag)
    {
        apply_coa_two_lineage.reshape_segs(pool_seg, cumul_nbr_brkpt_per_seg);
    }

    return apply_coa_two_lineage.New_first_seg;
}

//Used for transform diploid indiv in 2 haploid indiv
void segregate_chr_btw_indiv(lattice_c &lat, indiv_stock_c &indiv_vec, std::size_t chr_nbr, rand_gen_c &rand_gen)
{
    auto size = indiv_vec.size();
    //Browse the already_existing_indiv
    for (std::size_t indiv = 0; indiv < size; ++indiv)
    {
        auto coord_node = indiv_vec[indiv]->Node_lat->Coord;
        indiv_c *new_indiv_ptr = nullptr;
        //not see ligneage in indiv
        bool possibly_empty = true;
        bool no_new_indiv = true;
        //first chr

        for (std::size_t chr_index = 0; chr_index < chr_nbr; ++chr_index)
        {
            auto chr0_ptr = indiv_vec[indiv]->get_chr(chr_index, 0);
            //if chr0 empty chr1 empty to, nothing happend
            if (chr0_ptr != nullptr)
            {
                //don't have chr0 not empty before in indiv[i]
                if (possibly_empty)
                {
                    auto chr1_ptr = indiv_vec[indiv]->get_chr(chr_index, 1);
                    //no ramdom here because they are the first pair of chr to be separeted
                    if (chr1_ptr != nullptr)
                    {
                        indiv_vec[indiv]->untie_chr(chr_index, 1);
                        if (no_new_indiv)
                        {
                            new_indiv_ptr = indiv_vec.new_indiv(lat.node(coord_node));
                            no_new_indiv = false;
                        }
                        new_indiv_ptr->update_indiv(chr_index, chr1_ptr, nullptr);
                    }
                    possibly_empty = false;
                }
                //have a non empty chr0 before in indiv[i]
                else
                {
                    if (rand_gen.rand_bool())
                    {
                        //Exchange between hapologenome
                        indiv_vec[indiv]->untie_chr(chr_index, 0);
                        if (no_new_indiv)
                        {
                            new_indiv_ptr = indiv_vec.new_indiv(lat.node(coord_node));
                            no_new_indiv = false;
                        }
                        new_indiv_ptr->update_indiv(chr_index, chr0_ptr, nullptr);
                        //push chr1 to chr0 placeholder
                        auto chr1_ptr = indiv_vec[indiv]->get_chr(chr_index, 1);
                        if (chr1_ptr != nullptr)
                        {
                            indiv_vec[indiv]->update_indiv(chr_index, chr1_ptr, nullptr);
                        }
                    }
                    else
                    {
                        //p=0.5 to stay in the same haplogenome
                        auto chr1_ptr = indiv_vec[indiv]->get_chr(chr_index, 1);
                        if (chr1_ptr != nullptr)
                        {
                            indiv_vec[indiv]->untie_chr(chr_index, 1);
                            if (no_new_indiv)
                            {
                                new_indiv_ptr = indiv_vec.new_indiv(lat.node(coord_node));
                                no_new_indiv = false;
                            }
                            new_indiv_ptr->update_indiv(chr_index, chr1_ptr, nullptr);
                        }
                    }
                }
            }
        }

        if (!no_new_indiv)
        {
            lat.add_indiv(new_indiv_ptr, coord_node);
        }
    }
}
