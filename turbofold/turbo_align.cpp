#include "turbofold.hpp"

double TurboAlign::beam_prune(std::unordered_map<std::pair<int, int>, HState, PairHash> &beamstep, HStateType h,
                              int beam_size) {
    std::vector<std::pair<double, std::pair<int, int>>> scores;
    scores.reserve(beamstep.size());
    for (auto &item : beamstep) {
        std::pair<int, int> aij = item.first;
        HState &cand = item.second;

        double offset = 0;
        if (use_prev_outside_score && turbofold->itr > 1) {
            int s = aij.first + aij.second;
            std::unordered_map<std::pair<int, int>, HState, PairHash> *prev_beamstep;
            switch (h) {
                case INS1:
                    prev_beamstep = &ab.bestINS1[s];
                    break;
                case INS2:
                    prev_beamstep = &ab.bestINS2[s];
                    break;
                case ALN:
                    prev_beamstep = &ab.bestALN[s];
                    break;
            }
            if (prev_beamstep->find(aij) != prev_beamstep->end()) {
                offset = prev_beamstep->at(aij).beta;
            } else {
                offset = NEG_INF;
            }
        }
        scores.push_back(std::make_pair(cand.alpha + offset, aij));
    }
    if (scores.size() <= beam_size) return NEG_INF;
    double threshold = Utility::quickselect(scores, 0, scores.size() - 1, scores.size() - beam_size);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }
    return threshold;
}
