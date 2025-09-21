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
    std::vector<std::pair<const Sequence*, const Structure*>> seqs_strucs;
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

    inline void save_threshknot_db(const std::string& out_dir) const {
        std::filesystem::create_directories(out_dir);

        std::ofstream out(out_dir + "/threshknot_db.fasta");
        if (!out) {
            throw std::runtime_error("Failed to open threshknot_db.fasta for writing in: " + out_dir);
        }

        for (const auto& [seq_ptr, struc_ptr] : seqs_strucs) {
            if (!seq_ptr || !struc_ptr) continue;  // safety

            out << ">" << seq_ptr->name << "\n";
            out << seq_ptr->seq << "\n";
            out << struc_ptr->getDotBracket() << "\n";  // assuming Structure has a getDotBracket() method
        }

        out.close();
    }

    inline void save_logs(const std::string& out_dir) {
        std::filesystem::create_directories(out_dir);

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
                for (const auto& log : aln_logs[itr]) {
                    aln_log << "----- Pair Index: " << pair_idx++ << " -----\n";
                    aln_log << "  Sequence Identity: " << log.seq_identity << "\n";
                    aln_log << "  Best Execution Time (ms): " << log.best_exec_time << " ms\n";
                    aln_log << "  Inside Execution Time (ms): " << log.exec_time.first << " ms\n";
                    double outside_pct =
                        (log.exec_time.first > 0.0) ? (100.0 * log.exec_time.second / log.exec_time.first) : 0.0;
                    aln_log << "  Outside Execution Time (ms): " << log.exec_time.second << " ms (" << outside_pct
                            << "% of inside)\n";
                    aln_log << "  Coincidence Probs Execution Time: " << log.cp_exec_time << " ms\n";
                    aln_log << "  Scores (inside, outside): " << log.total_score.first << ", " << log.total_score.second
                            << "\n";
                    aln_log << "  Effective Beam Size (inside, outside): " << log.effective_beam_size.first << ", "
                            << log.effective_beam_size.second << "\n";
                }
                aln_log.close();
            }

            // -- Partition logs (always present)
            std::ofstream pf_log(pf_log_dir + "/itr_" + std::to_string(itr) + ".txt");
            if (!pf_log)
                throw std::runtime_error("Failed to open partition log file for iteration " + std::to_string(itr));

            pf_log << "=== Partition Logs: Iteration " << itr << " ===\n\n";
            int seq_idx = 0;
            for (const auto& log : pf_logs[itr]) {
                pf_log << "----- Sequence Index: " << seq_idx++ << " -----\n";
                pf_log << "  Free Energy of Ensemble: " << log.free_energy_of_ensemble << "\n";
                pf_log << "  Inside Execution Time (ms): " << log.exec_time.first << " ms\n";
                double outside_pct =
                    (log.exec_time.first > 0.0) ? (100.0 * log.exec_time.second / log.exec_time.first) : 0.0;
                pf_log << "  Outside Execution Time (ms): " << log.exec_time.second << " ms (" << outside_pct
                       << "% of inside)\n";
                pf_log << "  BPP Execution Time (ms): " << log.bpp_exec_time << " ms\n";
                pf_log << "  Total Energy (inside, outside): " << log.total_energy.first << ", "
                       << log.total_energy.second << "\n";
                pf_log << "  Effective Beam Size (inside, outside): " << log.effective_beam_size.first << ", "
                       << log.effective_beam_size.second << "\n";
            }
            pf_log.close();
        }

        itrs_exec_time = 0.0;
        for (const auto& t : aln_itr_exec_times) itrs_exec_time += t;
        for (const auto& t : pf_itr_exec_times) itrs_exec_time += t;

        // --- write top-level summary file
        std::ofstream summary(out_dir + "/logs/log.txt");
        if (!summary) {
            throw std::runtime_error("Failed to open summary log file for writing: " + out_dir + "/logs/log.txt");
        }

        summary << "num_iterations: " << num_iterations << "\n";
        summary << "use_lazy_outside: " << use_lazy_outside << "\n";
        summary << "use_prev_itr_beta: " << use_prev_itr_beta << "\n";
        summary << "restrict_search: " << restrict_search << "\n";
        summary << "alignment_pruning_threshold: " << alignment_pruning_threshold << "\n";
        summary << "folding_pruning_threshold: " << folding_pruning_threshold << "\n";

        summary << "aln_itr_exec_times (ms): ";
        for (const auto& t : aln_itr_exec_times) summary << t << " ";
        summary << "\n";

        summary << "pf_itr_exec_times (ms): ";
        for (const auto& t : pf_itr_exec_times) summary << t << " ";
        summary << "\n";

        summary << "itrs_exec_time (ms): " << itrs_exec_time << "\n";
        summary.close();
    }
};
