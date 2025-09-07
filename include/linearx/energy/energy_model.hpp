// linearx/energy/energy_model.hpp
#pragma once

#include <iostream>
#include <linearx/energy/params/EnergyParams.hpp>
#include <vector>

#define NUC_TO_PAIR(x, y)                \
    (x == 1 ? (y == 4 ? 5 : 0)           \
            : (x == 2 ? (y == 3 ? 1 : 0) \
                      : (x == 3 ? (y == 2 ? 2 : (y == 4 ? 3 : 0)) : (x == 4 ? (y == 3 ? 4 : (y == 1 ? 6 : 0)) : 0))))

class EnergyModel {
   private:
    inline const static char Triloops[241] =
        "CAACG "
        "GUUAC ";

    inline const static char Tetraloops[281] =
        "CAACGG "
        "CCAAGG "
        "CCACGG "
        "CCCAGG "
        "CCGAGG "
        "CCGCGG "
        "CCUAGG "
        "CCUCGG "
        "CUAAGG "
        "CUACGG "
        "CUCAGG "
        "CUCCGG "
        "CUGCGG "
        "CUUAGG "
        "CUUCGG "
        "CUUUGG ";

    inline const static char Hexaloops[361] =
        "ACAGUACU "
        "ACAGUGAU "
        "ACAGUGCU "
        "ACAGUGUU ";

   public:
    const EnergyParams epm;
    const bool use_special_hairpin;

    EnergyModel(EnergyParamsType epm_type) : epm(epm_type), use_special_hairpin(true) {}

    inline void init_tetra_hex_tri(const std::string &seq, const int seq_length, std::vector<int> &if_tetraloops,
                                   std::vector<int> &if_hexaloops, std::vector<int> &if_triloops) const {
        // TetraLoops
        if_tetraloops.resize(seq_length - 5 < 0 ? 0 : seq_length - 5, -1);
        for (int i = 0; i < seq_length - 5; ++i) {
            if (!(seq[i] == 'C' && seq[i + 5] == 'G')) continue;
            const char *ts;
            if ((ts = strstr(Tetraloops, seq.substr(i, 6).c_str()))) if_tetraloops[i] = (ts - Tetraloops) / 7;
        }

        // Triloops
        if_triloops.resize(seq_length - 4 < 0 ? 0 : seq_length - 4, -1);
        for (int i = 0; i < seq_length - 4; ++i) {
            if (!((seq[i] == 'C' && seq[i + 4] == 'G') || (seq[i] == 'G' && seq[i + 4] == 'C'))) continue;
            const char *ts;
            if ((ts = strstr(Triloops, seq.substr(i, 5).c_str()))) if_triloops[i] = (ts - Triloops) / 6;
        }

        // Hexaloops
        if_hexaloops.resize(seq_length - 7 < 0 ? 0 : seq_length - 7, -1);
        for (int i = 0; i < seq_length - 7; ++i) {
            if (!(seq[i] == 'A' && seq[i + 7] == 'U')) continue;
            const char *ts;
            if ((ts = strstr(Hexaloops, seq.substr(i, 8).c_str()))) if_hexaloops[i] = (ts - Hexaloops) / 9;
        }
    }

    inline int score_hairpin(const int i, const int j, const int nuci, const int nuci1, const int nucj_1,
                             const int nucj, const int tetra_hex_tri_index = -1) const {
        const int size = j - i - 1;
        const int type = NUC_TO_PAIR(nuci, nucj);

        int energy = 0;
        if (size <= 30) {
            energy += (*epm.hairpin37)[size];
        } else {
            energy += (*epm.hairpin37)[30] + (int)(linearx::constants::energy::LXC37 * log((size) / 30.));
        }

        bool specialHP_exists = (*epm.Tetraloop37 && *epm.Hexaloop37 && *epm.Triloop37);

        // [TODO] uncomment
        // if (specialHP_exists && size < 3) return energy; /* should only be the case when folding alignments */
        // if (!specialHP_exists && size == 3) return energy + ((type > 2 || type == 0) ? epm.TerminalAU37 : 0);

        if (specialHP_exists && use_special_hairpin) {
            if (size == 4 && tetra_hex_tri_index > -1) {
                return (*epm.Tetraloop37)[tetra_hex_tri_index];
            } else if (size == 6 && tetra_hex_tri_index > -1) {
                return (*epm.Hexaloop37)[tetra_hex_tri_index];
            } else if (size == 3) {
                if (tetra_hex_tri_index > -1) {
                    return (*epm.Triloop37)[tetra_hex_tri_index];
                }
                return (energy + (type > 2 ? epm.TerminalAU37 : 0));
            }
        }        

        energy += (*epm.mismatchH37)[type][nuci1][nucj_1];
        return energy;
    };

    inline int score_single_loop(const int i, const int j, const int p, const int q, const int nuci, const int nuci1,
                                 const int nucj_1, const int nucj, const int nucp_1, const int nucp, const int nucq,
                                 const int nucq1) const {
        const int type = NUC_TO_PAIR(nuci, nucj);
        const int type_2 = NUC_TO_PAIR(nucq, nucp);
        const int n1 = p - i - 1;
        const int n2 = j - q - 1;
        const int nl = (n1 > n2) ? n1 : n2;
        const int ns = (n1 > n2) ? n2 : n1;
        int energy = 0;

        // return v_score_single(i, j, p, q, nuci, nuci1, nucj_1, nucj, nucp_1, nucp, nucq, nucq1);

        if (nl == 0) return (*epm.stack37)[type][type_2]; /* stack */

        if (ns == 0) { /* bulge */
            energy = (nl <= linearx::constants::energy::MAXLOOPSIZE)
                         ? (*epm.bulge37)[nl]
                         : ((*epm.bulge37)[30] + (int)(linearx::constants::energy::LXC37 * log(nl / 30.)));
            if (nl == 1)
                energy += (*epm.stack37)[type][type_2];
            else {
                if (type > 2 || type == 0) energy += epm.TerminalAU37;
                if (type_2 > 2 || type_2 == 0) energy += epm.TerminalAU37;
            }
            return energy;
        } else { /* interior loop */
            if (ns == 1) {
                if (nl == 1) /* 1x1 loop */
                    return (*epm.int11_37)[type][type_2][nuci1][nucj_1];
                if (nl == 2) { /* 2x1 loop */
                    if (n1 == 1)
                        energy = (*epm.int21_37)[type][type_2][nuci1][nucq1][nucj_1];
                    else
                        energy = (*epm.int21_37)[type_2][type][nucq1][nuci1][nucp_1];
                    return energy;
                } else if ((*epm.mismatch1nI37)) { /* 1xn loop */
                    energy = (nl + 1 <= linearx::constants::energy::MAXLOOPSIZE)
                                 ? ((*epm.internal_loop37)[nl + 1])
                                 : ((*epm.internal_loop37)[30] +
                                    (int)(linearx::constants::energy::LXC37 * log((nl + 1) / 30.)));
                    energy += std::min(epm.MAX_NINIO, (nl - ns) * epm.ninio37);
                    energy += (*epm.mismatch1nI37)[type][nuci1][nucj_1] + (*epm.mismatch1nI37)[type_2][nucq1][nucp_1];
                    return energy;
                } else {
                    goto generic;
                }
            } else if (ns == 2) {
                if (nl == 2) { /* 2x2 loop */
                    return (*epm.int22_37)[type][type_2][nuci1][nucp_1][nucq1][nucj_1];
                } else if (nl == 3 && (*epm.mismatch23I37)) { /* 2x3 loop */
                    energy = (*epm.internal_loop37)[5] + epm.ninio37;
                    energy += (*epm.mismatch23I37)[type][nuci1][nucj_1] + (*epm.mismatch23I37)[type_2][nucq1][nucp_1];
                    return energy;
                } else {
                    goto generic;
                }
            }

        generic: { /* generic interior loop (no else here!)*/
            const int u = nl + ns;
            energy = (u <= linearx::constants::energy::MAXLOOPSIZE)
                         ? ((*epm.internal_loop37)[u])
                         : ((*epm.internal_loop37)[30] + (int)(linearx::constants::energy::LXC37 * log((u) / 30.)));

            energy += std::min(epm.MAX_NINIO, (nl - ns) * epm.ninio37);

            energy += (*epm.mismatchI37)[type][nuci1][nucj_1] + (*epm.mismatchI37)[type_2][nucq1][nucp_1];
        }
        }

        // std::cout << i << " " << j << " " << p << " " << q << " " << energy << std::endl;
        return energy;
    }

    // multi_loop
    inline int E_MLstem(const int type, const int si1, const int sj1) const {
        int energy = 0;

        if (*epm.mismatchM37) {
            if (si1 >= 0 && sj1 >= 0)
                energy += (*epm.mismatchM37)[type][si1][sj1];
            else if (si1 >= 0)
                energy += (*epm.dangle5_37)[type][si1];
            else if (sj1 >= 0)
                energy += (*epm.dangle3_37)[type][sj1];
        } else {
            std::cout << "Warning: mismatchM37 is not defined, using dangle energies instead." << std::endl;
            if (si1 > 0) energy += (*epm.dangle5_37)[type][si1];
            if (sj1 > 0) energy += (*epm.dangle3_37)[type][sj1];
        }

        if (type > 2 || type == 0) {
            energy += epm.TerminalAU37;
        }

        energy += epm.ML_intern37;

        return energy;
    }

    inline int score_multi(const int i, const int j, const int nuci, const int nuci1, const int nucj_1, const int nucj,
                           const int len) const {
        const int tt = NUC_TO_PAIR(nucj, nuci);
        return E_MLstem(tt, nucj_1, nuci1) + epm.ML_closing37;
    }

    inline int score_M1(const int i, const int j, const int k, const int nuci_1, const int nuci, const int nuck,
                        const int nuck1, const int len) const {
        const int tt = NUC_TO_PAIR(nuci, nuck);
        return E_MLstem(tt, nuci_1, nuck1);
    }

    inline int score_multi_unpaired(const int i, const int j) const { return epm.ML_BASE37 * (j - i); }

    // exterior_loop
    inline int score_external_paired(const int i, const int j, const int nuci_1, const int nuci, const int nucj,
                                     const int nucj1, const int len) const {
        const int type = NUC_TO_PAIR(nuci, nucj);
        int energy = 0;

        if (*epm.mismatchExt37) {
            if (nuci_1 >= 0 && nucj1 >= 0)
                energy += (*epm.mismatchExt37)[type][nuci_1][nucj1];
            else if (nuci_1 >= 0)
                energy += (*epm.dangle5_37)[type][nuci_1];
            else if (nucj1 >= 0)
                energy += (*epm.dangle3_37)[type][nucj1];
        } else {
            if (nuci_1 > 0) energy += (*epm.dangle5_37)[type][nuci_1];
            if (nucj1 > 0) energy += (*epm.dangle3_37)[type][nucj1];
        }

        if (type > 2 || type == 0) {
            energy += epm.TerminalAU37;
        }

        return energy;
    }

    inline int score_external_unpaired(const int i, const int j) const { return 0; }
};
