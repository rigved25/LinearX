#include "turbofold.hpp"

LinearTurboFold::LinearTurboFold(MultiSeq *multi_seq, const int num_itr)
    : multi_seq(multi_seq), energy_model(EnergyParamsType::VIENNA), num_itr(num_itr) {

    for (int i = 0; i < multi_seq->size(); i++) {
        multi_seq->at(i).k_id = i;                                 // set k_id for each sequence
        pfs.emplace_back(this, &(multi_seq->at(i)), energy_model); // better than pfs.push_back(),
                                                                   // creates Partition object directly
                                                                   // inside the container

        // enumerate all possible k^2 sequence pairs and create LinearAlign objects
        for (int j = i + 1; j < multi_seq->size(); j++) {
            alns.emplace_back(&(multi_seq->at(i)), &(multi_seq->at(j)));
        }
    }
}

LinearAlign *LinearTurboFold::get_alignment_obj(int k1, int k2) {
    // assert k1 and k2 less than multi_seq->size()
    if (k1 >= multi_seq->size() || k2 >= multi_seq->size()) {
        throw std::out_of_range("Seq index k out of range");
    }
    if (k1 > k2) {
        return this->get_alignment_obj(k2, k1);
    }
    int index = k1 * (multi_seq->size() - 1) - (k1 * (k1 + 1)) / 2 + (k2 - k1 - 1);
    return &(alns[index]);
}

float LinearTurboFold::get_extrinsic_info(const Seq *x, const int i, const int j) {
    float output = 0.0;
    for (const Seq &y : *multi_seq) {
        if (x->k_id == y.k_id)
            continue;
        if (pfs[y.k_id].bpp == nullptr) {
            // std::cout << "[WARNING] BPP matrix for seq k: " << y.k_id << " is not computed yet." << std::endl;
            return 1.0; // remove later
        }

        LinearAlign *aln_obj = get_alignment_obj(x->k_id, y.k_id);
        for (int q = 0; q < y.length(); ++q) {
            for (const auto &item : pfs[y.k_id].bpp[q]) {
                const int p = item.first;
                const float prob = item.second;
                output += prob * aln_obj->get_bpp(i, p) * aln_obj->get_bpp(j, q);
            }
        }
    }
    return output;
}

void LinearTurboFold::run() {
    // step 0: initially compute partition function for all k sequences
    for (TurboPartition &pf : pfs) {
        pf.compute_inside();
        pf.compute_outside();
        pf.calc_prob_accm();
        pf.print_alpha_beta();
    }

    for (int itr = 1; itr < 2; ++itr) {
        // align step
        for (LinearAlign &aln : alns) {
            int k1 = aln.sequence1->k_id;
            int k2 = aln.sequence2->k_id;
            aln.set_prob_accm(pfs[k1].prob_accm, pfs[k2].prob_accm);
            aln.compute_inside();
            aln.run_backward_phase(false);
            // aln.compute_outside();
        }

        // fold step
        for (TurboPartition &pf : pfs) {
            pf.compute_inside();
            pf.compute_outside();
            pf.calc_prob_accm();
        }
    }

    // p1->compute_inside();
    // p1->compute_outside();
    // p2->compute_inside();
    // p2->compute_outside();
    // p1->calc_prob_accm();
    // p2->calc_prob_accm();

    // linAln->set_prob_accm(p1->prob_accm, p2->prob_accm);
    // linAln->compute_inside();

    // std::cout << "Traceback:" << std::endl;
    // linAln->traceback();
}