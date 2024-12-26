#include "./partition.hpp"

double Partition::beam_prune(std::unordered_map<int, State> &beamstep, int beam_size) {
    if (beam_size == 0 || beamstep.size() <= beam_size) {
        return xlog(0.0);
    }
    std::vector<std::pair<double, int>> scores;
    for (auto &item : beamstep) {
        int i = item.first;
        State &cand = item.second;
        int k = i - 1;
        double newalpha = ((k >= 0 ? bestC[k].alpha : double(0.0)) + cand.alpha);
        scores.push_back(std::make_pair(newalpha, i));
    }
    if (scores.size() <= beam_size) {
        return xlog(0.0);
    }
    double threshold = Utility::quickselect(scores, 0, scores.size() - 1, scores.size() - beam_size);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }
    return threshold;
}

bool Partition::check_state(StateType type, int i, int j) const { return true; }

void Partition::update_score(State &state, const int new_score, const double prev_score) {
    if (mode == MFE) {
        state.alpha = std::max(state.alpha, prev_score + new_score);
    } else {
        Fast_LogPlusEquals(state.alpha, prev_score + new_score * INV_KT);
    }
}

void Partition::compute_inside(int beam_size) {
    auto start_time = std::chrono::high_resolution_clock::now();
    if (verbose_output) {
        fprintf(stderr, "[LinearPartition] (Seq k_id: %d) Running Inside Algorithm:\n", sequence->k_id);
    }
    for (int j = 0; j < seq->size(); j++) {
        if (verbose_output) {
            Utility::showProgressBar(j, seq->size() - 1);
        }

        // beam of H
        beam_prune(bestH[j], beam_size);
        beamstep_H(j, next_pair);
        if (j == 0) {
            continue;
        }

        // beam of Multi
        beam_prune(bestMulti[j], beam_size);
        beamstep_Multi(j, next_pair);

        // beam of P
        beam_prune(bestP[j], beam_size);
        beamstep_P(j, next_pair);

        // beam of M2
        beam_prune(bestM2[j], beam_size);
        beamstep_M2(j, next_pair);

        // beam of M
        beam_prune(bestM[j], beam_size);
        beamstep_M(j);

        beamstep_C(j);  // beam of C
    }

    // set viterbi pointer
    this->viterbi = &bestC[seq->size() - 1];

    // update/print time stats
    auto end_time = std::chrono::high_resolution_clock::now();
    inside_execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    if (verbose_output) {
        std::cerr << "  - Execution Time: " << inside_execution_time << " ms" << std::endl;
        if (mode == InsideMode::MFE) {
            std::cerr << "  - MFE (Minimum Free Energy): " << bestC[seq->size() - 1].alpha / -100.0 << " kcal/mol"
                      << std::endl;
            printf("\n%s (%.2f)\n", get_mfe_structure().c_str(), bestC[seq->size() - 1].alpha / -100.0);
        } else {
            std::cerr << "  - Free Energy of the Ensemble: " << get_ensemble_energy() << " kcal/mol" << std::endl;
        }
        std::cerr << std::endl;
    }
}

void Partition::beamstep_H(const int j, const std::vector<int> *next_pair) {
    int jnext = next_pair[seq->at(j)][j];
    while (!allow_sharp_turn && jnext < seq->size() && (jnext - j) < 4) {
        jnext = next_pair[seq->at(j)][jnext];
    }

    if (jnext < seq->size() && check_state(StateType::H, j, jnext)) {
        int tetra_hex_tri = -1;
        if (jnext - j - 1 == 4)  // 6:tetra
            tetra_hex_tri = if_tetraloops.at(j);
        else if (jnext - j - 1 == 6)  // 8:hexa
            tetra_hex_tri = if_hexaloops.at(j);
        else if (jnext - j - 1 == 3)  // 5:tri
            tetra_hex_tri = if_triloops.at(j);

        int score = -energy_model->score_hairpin(j, jnext, seq->at(j), seq->at(j + 1), seq->at(jnext - 1),
                                                 seq->at(jnext), tetra_hex_tri);
        update_score(bestH[jnext][j], score);
    }

    // for every state h in H[j]
    //   1. extend H(i, j) to H(i, jnext)
    //   2. generate P(i, j)
    for (auto &item : bestH[j]) {
        int i = item.first;
        State &state = item.second;

        // 1. extend H(i, j) to H(i, jnext)
        int jnext = next_pair[seq->at(i)][j];
        if (jnext < seq->size() && check_state(StateType::H, i, jnext)) {
            int tetra_hex_tri = -1;
            if (jnext - i - 1 == 4)  // 6:tetra
                tetra_hex_tri = if_tetraloops.at(i);
            else if (jnext - i - 1 == 6)  // 8:hexa
                tetra_hex_tri = if_hexaloops.at(i);
            else if (jnext - i - 1 == 3)  // 5:tri
                tetra_hex_tri = if_triloops.at(i);

            int score = -energy_model->score_hairpin(i, jnext, seq->at(i), seq->at(i + 1), seq->at(jnext - 1),
                                                     seq->at(jnext), tetra_hex_tri);
            update_score(bestH[jnext][i], score);
        }

        // 2. generate P(i, j)
        if (check_state(StateType::P, i, j)) {
            update_score(bestP[j][i], 0, state.alpha);
        }
    }
}

void Partition::beamstep_Multi(const int j, const std::vector<int> *next_pair) {
    for (auto &item : bestMulti[j]) {
        int i = item.first;
        State &state = item.second;

        // 1. extend Multi(i, j) to Multi(i, jnext)
        int jnext = next_pair[seq->at(i)][j];
        if (jnext < seq->size() && check_state(StateType::Multi, i, jnext)) {
            int new_score = -energy_model->score_multi_unpaired(j, jnext);
            update_score(bestMulti[jnext][i], new_score, state.alpha);
        }

        // 2. generate P(i, j)
        if (check_state(StateType::P, i, j)) {
            int new_score =
                -energy_model->score_multi(i, j, seq->at(i), seq->at(i + 1), seq->at(j - 1), seq->at(j), seq->size());
            update_score(bestP[j][i], new_score, state.alpha);
        }
    }
}

void Partition::beamstep_P(const int j, const std::vector<int> *next_pair) {
    // for every state in P[j]
    //   1. generate new P (helix/bulge)
    //   2. M = P
    //   3. M2 = M + P
    //   4. C = C + P
    for (auto &item : bestP[j]) {
        int i = item.first;
        State &state = item.second;

        // 1. generate new P (helix/bulge)
        for (int p = i - 1; p >= std::max(0, i - MAXLOOPSIZE); --p) {
            int q = next_pair[seq->at(p)][j];
            while (q < seq->size() && ((i - p) + (q - j) - 2) <= MAXLOOPSIZE) {
                // current shape is: p...i (pair) j...q
                if (check_state(StateType::P, p, q)) {
                    int new_score = -energy_model->score_single_loop(p, q, i, j, seq->at(p), seq->at(p + 1),
                                                                     seq->at(q - 1), seq->at(q), seq->at(i - 1),
                                                                     seq->at(i), seq->at(j), seq->at(j + 1));
                    update_score(bestP[q][p], new_score, state.alpha);
                }

                q = next_pair[seq->at(p)][q];
            }
        }

        // 2. M = P
        if (i > 0 && j < seq->size() - 1 && check_state(StateType::M, i, j)) {
            int new_score = -energy_model->score_M1(i, j, j, seq->at(i - 1), seq->at(i), seq->at(j),
                                                    (j + 1 < seq->size() ? seq->at(j + 1) : -1), seq->size());
            update_score(bestM[j][i], new_score, state.alpha);
        }

        // 3. M2 = M + P
        int h = i - 1;
        if (h > 0 && !bestM[h].empty()) {
            int new_score = -energy_model->score_M1(i, j, j, seq->at(i - 1), seq->at(i), seq->at(j),
                                                    (j + 1 < seq->size() ? seq->at(j + 1) : -1), seq->size());
            for (auto &m_item : bestM[h]) {
                int g = m_item.first;

                if (check_state(StateType::M2, g, j)) {
                    State &m_state = m_item.second;
                    update_score(bestM2[j][g], new_score, m_state.alpha + state.alpha);
                }
            }
        }

        // 4. C = C + P
        h = i - 1;
        if (h >= 0) {
            State &prefix_C = bestC[h];
            int new_score = -energy_model->score_external_paired(
                i, j, seq->at(h), seq->at(h + 1), seq->at(j), (j + 1 < seq->size() ? seq->at(j + 1) : -1), seq->size());
            update_score(bestC[j], new_score, prefix_C.alpha + state.alpha);
        } else {
            int new_score = -energy_model->score_external_paired(
                0, j, -1, seq->at(0), seq->at(j), (j + 1 < seq->size() ? seq->at(j + 1) : -1), seq->size());
            update_score(bestC[j], new_score, state.alpha);
        }
    }
}

void Partition::beamstep_M2(const int j, const std::vector<int> *next_pair) {
    for (auto &item : bestM2[j]) {
        int i = item.first;
        State &state = item.second;

        // 1. multi-loop = M2
        for (int p = i - 1; p >= std::max(0, i - MAXLOOPSIZE); --p) {
            int q = next_pair[seq->at(p)][j];
            if (q < seq->size() && check_state(StateType::Multi, p, q)) {
                int new_score =
                    -(energy_model->score_multi_unpaired(p, i - 1) + energy_model->score_multi_unpaired(j, q - 1));
                update_score(bestMulti[q][p], new_score, state.alpha);
            }
        }

        // 2. M = M2
        if (check_state(StateType::M, i, j)) {
            update_score(bestM[j][i], 0, state.alpha);
        }
    }
}

void Partition::beamstep_M(const int j) {
    for (auto &item : bestM[j]) {
        int i = item.first;
        State &state = item.second;

        if (j < seq->size() - 1 && check_state(StateType::M, i, j + 1)) {
            int new_score = -energy_model->score_multi_unpaired(j, j + 1);
            update_score(bestM[j + 1][i], new_score, state.alpha);
        }
    }
}

void Partition::beamstep_C(const int j) {
    if (j < seq->size() - 1) {
        update_score(bestC[j + 1], 0, bestC[j].alpha);
    }
}
