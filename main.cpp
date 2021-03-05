#include <chrono>

#include "simulator.hpp"
#include "output.hpp"

/*******Writing Code Rules********/
// - Variable name in lowercase with underscore (lower_snake_case)
// - MAJ : template variable name first letter + objet data variable name first letter
//
//
//
//
//
//
//
/*********************************/

std::string const datestring = __DATE__;
std::string const timestring = __TIME__;
std::string const version = "0.1 (Built on " + datestring + " at " + timestring + ")";

namespace fs = std::filesystem;
int main(int argc, char *argv[])
{
    //Allow to approximately calculate run-time duration
    auto debut = std::chrono::high_resolution_clock::now();

    //Search for settings name file
    std::string cmd_str = read_write_cmdline(argc, argv);
    parser_str(cmd_str);

    //Construct a instance of param_c called param
    auto const &simu_param = singleton_c<simu_param_c>::instance();

    std::cout << "reading settings file : " << simu_param.Setting_filename << "\n"
              << std::endl;

    auto file_str = read_file(simu_param.Setting_filename);
    parser_str(file_str);
    //clean str for empty memory
    file_str.erase();

    //Arg priority
    parser_str(cmd_str);
    cmd_str.erase();

    check_param();
    apply_param();

    //Construct a instance of param_c called param
    auto const &samp_param = singleton_c<samp_param_c>::instance();
    auto const &demo_param = singleton_c<demo_param_c>::instance();
    auto const &muta_param = singleton_c<muta_param_c>::instance();
    auto const &recomb_param = singleton_c<recomb_param_c>::instance();

    output_screen_info(version, simu_param, samp_param, demo_param, muta_param, recomb_param);

    auto &info_collect = singleton_c<info_collector_c>::instance();
    write_beforerun_param_settings_summary(simu_param.Generic_data_filename + simu_param.Param_summary_filename, version, simu_param, samp_param, demo_param, muta_param, recomb_param);

    output_stat_c output_stat;
    
    if (info_collect.Clock)
    {
        auto fin = std::chrono::high_resolution_clock::now();
        info_collect.time_before_simulation += (std::chrono::duration_cast<std::chrono::nanoseconds>(fin - debut).count() * 0.000000001);
    }


    for (int rep = 0; rep < simu_param.Repetition_nbr; ++rep)
    {
        std::cout << "Simulation " << rep + 1 << std::endl;
        sample_simulator(output_stat, rep);
    }

    std::chrono::time_point<std::chrono::high_resolution_clock> debut_after_simulation;
    if (info_collect.Clock)
    {
        debut_after_simulation = std::chrono::high_resolution_clock::now();
    }

    
    write_afterrun_param_settings_summary(simu_param.Generic_data_filename + simu_param.Param_summary_filename, simu_param, samp_param, demo_param);

    std::cout << std::endl;

    if (info_collect.MRCA_record)
    {
        coa_tree_metrics_c coa_tree_metric;
        coa_tree_metric.calcul_coa_tree_metrics(simu_param.Continuous_time_approxim, samp_param.Ploidy, demo_param.Population_size_N, samp_param.n_total_sample_size);

        std::cout << "    ************    " << std::endl;
        std::cout << "Theo MRCA mean : " << coa_tree_metric.Theo_MRCA_mean << std::endl;
        std::cout << "MRCA mean : " << info_collect.MRCA_record_cumul_mean_var[0] << " (" << info_collect.MRCA_record_cumul_mean_var[1] << ")" << std::endl;
        std::cout << "        " << std::endl;
        std::cout << "Theo 2 lineage coa time mean : " << coa_tree_metric.Theo_2_lign_coa_time_mean << std::endl;
        std::cout << "2 lineage coa time mean : Implement soon" << std::endl;
        std::cout << "        " << std::endl;
    }
    if (info_collect.Prob_id_1_loc_Qr)
    {
        std::cout << "    ************    " << std::endl;
        std::cout << "Qr0 mean : " << output_stat.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(0) << std::endl;
        std::cout << "Qr0 var : " << output_stat.Prob_id_1_loc_Qr_res.Qr_cumul_m_v[0].at(1) << std::endl;
        std::cout << "        " << std::endl;
    }
    if (info_collect.Prob_id_1_2_loc)
    {
        std::cout << "    ************    " << std::endl;
        std::cout << "Phi mean : " << std::get<0>(output_stat.Prob_id_1_2_loc_res.PHI_cumul_m_v) << std::endl;
        std::cout << "Phi var : " << std::get<1>(output_stat.Prob_id_1_2_loc_res.PHI_cumul_m_v) << std::endl;
        std::cout << "        " << std::endl;
    }
    if (info_collect.Prob_id_1_loc_Qwi_wd_bd)
    {
        std::cout << "    ************    " << std::endl;
        if (samp_param.Ploidy == 2)
        {
            std::cout << "Qwi mean : " << std::get<0>(output_stat.Prob_id_1_loc_Qwi_wd_bd_res.Qwi_cumul_m_v) << std::endl;
            std::cout << "Qwi var : " << std::get<1>(output_stat.Prob_id_1_loc_Qwi_wd_bd_res.Qwi_cumul_m_v) << std::endl;
            std::cout << "        " << std::endl;
        }
        std::cout << "Qwd mean : " << std::get<0>(output_stat.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v) << std::endl;
        std::cout << "Qwd var : " << std::get<1>(output_stat.Prob_id_1_loc_Qwi_wd_bd_res.Qwd_cumul_m_v) << std::endl;
        std::cout << "        " << std::endl;

        std::cout << "Qbd mean : " << std::get<0>(output_stat.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v) << std::endl;
        std::cout << "Qbd var : " << std::get<1>(output_stat.Prob_id_1_loc_Qwi_wd_bd_res.Qbd_cumul_m_v) << std::endl;
        std::cout << "        " << std::endl;
    }

    if (info_collect.Clock)
    {
        auto fin = std::chrono::high_resolution_clock::now();
        info_collect.time_after_simulation += (std::chrono::duration_cast<std::chrono::nanoseconds>(fin - debut_after_simulation).count() * 0.000000001);
    }

    
    std::cout << "Total execution time  is ";
    auto fin = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(fin - debut).count() * 0.000000001;
    std::cout << time << " seconds" << std::endl;

    if (info_collect.Clock)
    {
        std::cout << "***************" << std::endl;
        std::cout << "Intialization time" << std::endl;
        std::cout << (info_collect.time_before_simulation / time) * 100 << " % " << std::endl;
        std::cout << "***************" << std::endl;
        std::cout << "Tree simulation time" << std::endl;
        std::cout << (info_collect.time_simulation / time) * 100 << " % " << std::endl;
        std::cout << "############" << std::endl;
        std::cout << "Time for migration" << std::endl;
        std::cout << (info_collect.time_mig / time) * 100 << " %" << std::endl;
        std::cout << "Time for coalescence" << std::endl;
        std::cout << (info_collect.time_coa / time) * 100 << " %" << std::endl;
        std::cout << "Time for recombination" << std::endl;
        std::cout << (info_collect.time_recomb / time) * 100 << " %" << std::endl;
        std::cout << "***************" << std::endl;
        std::cout << "Time for tree reconstruction before mutation" << std::endl;
        std::cout << (info_collect.time_construct_tree / time) * 100 << " %" << std::endl;
        std::cout << "Time for mutation" << std::endl;
        std::cout << (info_collect.time_mutation / time) * 100 << " %" << std::endl;
        std::cout << "***************" << std::endl;
        std::cout << "After simulation time" << std::endl;
        std::cout << (info_collect.time_after_simulation / time) * 100 << " % " << std::endl;
    }

    if (simu_param.Wait_for_cin_input || simu_param.Wait_for_final_cin_input)
    {
        std::cout << "\n...Press any key to stop the program...\n" << std::endl;
        std::cin.get();
    }

    std::cout << "\nNormal ending of GSpace.\n" << std::endl;


    return EXIT_SUCCESS;
}
