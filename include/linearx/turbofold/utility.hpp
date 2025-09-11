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
    std::vector<std::vector<std::pair<AlignmentInsideLog, AlignmentOutsideLog>>>
        aln_logs;  // aln_logs[iteration][pair_index]
    std::vector<std::vector<std::pair<PartitionInsideLog, PartitionOutsideLog>>>
        pf_logs;  // pf_logs[iteration][sequence_index]
    std::vector<value_type> aln_itr_exec_times;
    std::vector<value_type> pf_itr_exec_times;
    value_type itrs_exec_time;

    TurboFoldLog(const int num_iterations, bool use_lazy_outside, bool use_prev_itr_beta, bool restrict_search)
        : num_iterations(num_iterations),
          use_lazy_outside(use_lazy_outside),
          use_prev_itr_beta(use_prev_itr_beta),
          restrict_search(restrict_search) {
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

        // write top-level summary file
        std::ofstream summary(out_dir + "/turbofold_log.txt");
        if (!summary) {
            throw std::runtime_error("Failed to open summary log file for writing: " + out_dir + "/turbofold_log.txt");
        }

        summary << "num_iterations: " << num_iterations << "\n";
        summary << "use_lazy_outside: " << use_lazy_outside << "\n";
        summary << "use_prev_itr_beta: " << use_prev_itr_beta << "\n";
        summary << "restrict_search: " << restrict_search << "\n";

        summary << "aln_itr_exec_times (ms): ";
        for (const auto &t : aln_itr_exec_times) summary << t << " ";
        summary << "\n";

        summary << "pf_itr_exec_times (ms): ";
        for (const auto &t : pf_itr_exec_times) summary << t << " ";
        summary << "\n";

        summary << "itrs_exec_time (ms): " << itrs_exec_time << "\n";
        summary.close();

        // log directories
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
                for (const auto &[inside, outside] : aln_logs[itr]) {
                    aln_log << "-- Pair Index: " << pair_idx++ << " --\n";
                    aln_log << "[Alignment Inside Log]\n";
                    aln_log << "  score: " << inside.score << "\n";
                    aln_log << "  execution_time: " << inside.execution_time << " ms\n";
                    aln_log << "  beam_size: " << inside.beam_size << "\n";
                    aln_log << "  nodes_pruned: " << inside.nodes_pruned << "\n";

                    aln_log << "[Alignment Outside Log]\n";
                    aln_log << "  score: " << outside.score << "\n";
                    aln_log << "  execution_time: " << outside.execution_time << " ms\n";
                    aln_log << "  deviation_threshold: " << outside.deviation_threshold << "\n";
                    aln_log << "  effective_beam_size: " << outside.effective_beam_size << "\n";
                    aln_log << "  nodes_visited: " << outside.nodes_visited << "\n";
                    aln_log << "  nodes_pruned: " << outside.nodes_pruned << "\n";
                    aln_log << "  edges_saved: " << outside.edges_saved << "\n";
                    aln_log << "  edges_pruned: " << outside.edges_pruned << "\n\n";
                }
                aln_log.close();
            }

            // -- Partition logs (always present)
            std::ofstream pf_log(pf_log_dir + "/itr_" + std::to_string(itr) + ".txt");
            if (!pf_log)
                throw std::runtime_error("Failed to open partition log file for iteration " + std::to_string(itr));

            pf_log << "=== Partition Logs: Iteration " << itr << " ===\n\n";
            int seq_idx = 0;
            for (const auto &[inside, outside] : pf_logs[itr]) {
                pf_log << "-- Sequence Index: " << seq_idx++ << " --\n";
                pf_log << "[Partition Inside Log]\n";
                pf_log << "  energy: " << inside.energy << "\n";
                pf_log << "  execution_time: " << inside.execution_time << " ms\n";
                pf_log << "  beam_size: " << inside.beam_size << "\n";
                pf_log << "  nodes_pruned: " << inside.nodes_pruned << "\n";

                pf_log << "[Partition Outside Log]\n";
                pf_log << "  total_inside_energy: " << outside.total_inside_energy << "\n";
                pf_log << "  total_outside_energy: " << outside.total_outside_energy << "\n";
                pf_log << "  execution_time: " << outside.execution_time << " ms\n";
                pf_log << "  deviation_threshold: " << outside.deviation_threshold << "\n";
                pf_log << "  effective_beam_size: " << outside.effective_beam_size << "\n";
                pf_log << "  nodes_visited: " << outside.nodes_visited << "\n";
                pf_log << "  nodes_pruned: " << outside.nodes_pruned << "\n";
                pf_log << "  edges_saved: " << outside.edges_saved << "\n";
                pf_log << "  edges_pruned: " << outside.edges_pruned << "\n\n";
            }
            pf_log.close();
        }
    }
};
