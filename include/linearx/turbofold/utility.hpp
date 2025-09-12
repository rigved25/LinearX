// linearx/turbofold/utility.hpp
#pragma once

#include <filesystem>
#include <fstream>
#include <linearx/alignment/utility.hpp>
#include <linearx/partition/utility.hpp>
#include <stdexcept>

struct TurboFoldLog {
    int num_iterations;
    bool use_lazy_outside;
    bool use_prev_itr_beta;
    bool restrict_search;
    float alignment_pruning_threshold;
    float folding_pruning_threshold;
    std::vector<std::vector<AlignmentLog>> aln_logs;  // aln_logs[iteration][pair_index]
    std::vector<std::vector<PartitionLog>> pf_logs;   // pf_logs[iteration][sequence_index]
    std::vector<value_type> aln_itr_exec_times;
    std::vector<value_type> pf_itr_exec_times;
    value_type itrs_exec_time;

    TurboFoldLog(const int num_iterations, bool use_lazy_outside, bool use_prev_itr_beta, bool restrict_search,
                 float alignment_pruning_threshold, float folding_pruning_threshold)
        : num_iterations(num_iterations),
          use_lazy_outside(use_lazy_outside),
          use_prev_itr_beta(use_prev_itr_beta),
          restrict_search(restrict_search),
          alignment_pruning_threshold(alignment_pruning_threshold),
          folding_pruning_threshold(folding_pruning_threshold) {
        aln_logs.resize(num_iterations + 1);
        pf_logs.resize(num_iterations + 1);
        aln_itr_exec_times.resize(num_iterations + 1);
        pf_itr_exec_times.resize(num_iterations + 1);
        itrs_exec_time = 0.0;
    }

    inline void save_logs(const std::string &out_dir) {
        std::filesystem::create_directories(out_dir);

        itrs_exec_time = 0.0;
        for (const auto &t : aln_itr_exec_times) itrs_exec_time += t;
        for (const auto &t : pf_itr_exec_times) itrs_exec_time += t;

        // --- write top-level summary file
        std::ofstream summary(out_dir + "/turbofold_log.txt");
        if (!summary) {
            throw std::runtime_error("Failed to open summary log file for writing: " + out_dir + "/turbofold_log.txt");
        }

        summary << "num_iterations: " << num_iterations << "\n";
        summary << "use_lazy_outside: " << use_lazy_outside << "\n";
        summary << "use_prev_itr_beta: " << use_prev_itr_beta << "\n";
        summary << "restrict_search: " << restrict_search << "\n";
        summary << "alignment_pruning_threshold: " << alignment_pruning_threshold << "\n";
        summary << "folding_pruning_threshold: " << folding_pruning_threshold << "\n";

        summary << "aln_itr_exec_times (ms): ";
        for (const auto &t : aln_itr_exec_times) summary << t << " ";
        summary << "\n";

        summary << "pf_itr_exec_times (ms): ";
        for (const auto &t : pf_itr_exec_times) summary << t << " ";
        summary << "\n";

        summary << "itrs_exec_time (ms): " << itrs_exec_time << "\n";
        summary.close();

        // --- log directories
        std::string aln_log_dir = out_dir + "/logs/alignment";
        std::string pf_log_dir = out_dir + "/logs/partition";
        std::filesystem::create_directories(aln_log_dir);
        std::filesystem::create_directories(pf_log_dir);

        for (int itr = 0; itr <= num_iterations; ++itr) {
            // -- Alignment logs (skip itr 0)
            if (itr > 0) {
                std::ofstream aln_log(aln_log_dir + "/itr_" + std::to_string(itr) + ".txt");
                if (!aln_log)
                    throw std::runtime_error("Failed to open alignment log file for iteration " + std::to_string(itr));

                aln_log << "=== Alignment Logs: Iteration " << itr << " ===\n\n";
                int pair_idx = 0;
                for (const auto &log : aln_logs[itr]) {
                    aln_log << "-- Pair Index: " << pair_idx++ << " --\n";
                    aln_log << "  seq_identity: " << log.seq_identity << "\n";
                    aln_log << "  total_inside_score: " << log.total_inside_score << "\n";
                    aln_log << "  total_outside_score: " << log.total_outside_score << "\n";
                    aln_log << "  inside_exec_time: " << log.inside_exec_time << " ms\n";
                    double outside_pct =
                        (log.inside_exec_time > 0.0) ? (100.0 * log.outside_exec_time / log.inside_exec_time) : 0.0;
                    aln_log << "  outside_exec_time: " << log.outside_exec_time << " ms (" << outside_pct
                            << "% of inside)\n";
                    aln_log << "  effective_beam_size: " << log.effective_beam_size << "\n";
                }
                aln_log.close();
            }

            // -- Partition logs (always present)
            std::ofstream pf_log(pf_log_dir + "/itr_" + std::to_string(itr) + ".txt");
            if (!pf_log)
                throw std::runtime_error("Failed to open partition log file for iteration " + std::to_string(itr));

            pf_log << "=== Partition Logs: Iteration " << itr << " ===\n\n";
            int seq_idx = 0;
            for (const auto &log : pf_logs[itr]) {
                pf_log << "-- Sequence Index: " << seq_idx++ << " --\n";
                pf_log << "  free_energy_of_ensemble: " << log.free_energy_of_ensemble << "\n";
                pf_log << "  total_inside_energy: " << log.total_inside_energy << "\n";
                pf_log << "  total_outside_energy: " << log.total_outside_energy << "\n";
                pf_log << "  inside_exec_time: " << log.inside_exec_time << " ms\n";
                double outside_pct =
                    (log.inside_exec_time > 0.0) ? (100.0 * log.outside_exec_time / log.inside_exec_time) : 0.0;
                pf_log << "  outside_exec_time: " << log.outside_exec_time << " ms (" << outside_pct
                       << "% of inside)\n";
                pf_log << "  effective_beam_size: " << log.effective_beam_size << "\n";
            }
            pf_log.close();
        }
    }
};
