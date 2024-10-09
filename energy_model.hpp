#ifndef ENERGY_MODEL_HPP
#define ENERGY_MODEL_HPP

#include "./energy_params/EnergyParams.hpp"
#include <vector>
#include <iostream>

#define NUC_TO_PAIR(x, y)                                                                                              \
    (x == 1 ? (y == 4 ? 5 : 0)                                                                                         \
            : (x == 2 ? (y == 3 ? 1 : 0)                                                                               \
                      : (x == 3 ? (y == 2 ? 2 : (y == 4 ? 3 : 0)) : (x == 4 ? (y == 3 ? 4 : (y == 1 ? 6 : 0)) : 0))))

class EnergyModel {

  private:
    inline const static char Triloops[241] = "CAACG "
                                             "GUUAC ";

    inline const static char Tetraloops[281] = "CAACGG "
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

    inline const static char Hexaloops[361] = "ACAGUACU "
                                              "ACAGUGAU "
                                              "ACAGUGCU "
                                              "ACAGUGUU ";

  public:
    const EnergyParams epm;
    const bool use_special_hairpin;

    EnergyModel(EnergyParamsType epm_type) : epm(epm_type), use_special_hairpin(true) {}

    void init_tetra_hex_tri(std::string &seq, int seq_length, std::vector<int> &if_tetraloops,
                            std::vector<int> &if_hexaloops, std::vector<int> &if_triloops);

    int score_hairpin(int i, int j, int nuci, int nuci1, int nucj_1, int nucj, int tetra_hex_tri_index = -1) {
        int size = j - i - 1;
        int type = NUC_TO_PAIR(nuci, nucj);

        int energy = size <= 30 ? (*epm.hairpin37)[size] : (*epm.hairpin37)[30] + (int)(LXC37 * log((size) / 30.));

        bool specialHP_exists = (*epm.Tetraloop37 && *epm.Hexaloop37 && *epm.Triloop37);
        if (specialHP_exists && size < 3)
            return energy; /* should only be the case when folding alignments */
        if (!specialHP_exists && size == 3)
            return energy + ((type > 2 || type == 0) ? epm.TerminalAU37 : 0);

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

    int score_single_loop(int i, int j, int p, int q, int nuci, int nuci1, int nucj_1, int nucj, int nucp_1, int nucp,
                          int nucq, int nucq1) {

        int type = NUC_TO_PAIR(nuci, nucj);
        int type_2 = NUC_TO_PAIR(nucq, nucp);
        int n1 = p - i - 1;
        int n2 = j - q - 1;
        int nl = (n1 > n2) ? n1 : n2;
        int ns = (n1 > n2) ? n2 : n1;
        int energy = 0;

        if (nl == 0)
            return (*epm.stack37)[type][type_2]; /* stack */

        if (ns == 0) { /* bulge */
            energy = (nl <= MAXLOOPSIZE) ? (*epm.bulge37)[nl] : ((*epm.bulge37)[30] + (int)(LXC37 * log(nl / 30.)));
            if (nl == 1)
                energy += (*epm.stack37)[type][type_2];
            else {
                if (type > 2 || type == 0)
                    energy += epm.TerminalAU37;
                if (type_2 > 2 || type == 0)
                    energy += epm.TerminalAU37;
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
                    energy = (nl + 1 <= MAXLOOPSIZE)
                                 ? ((*epm.internal_loop37)[nl + 1])
                                 : ((*epm.internal_loop37)[30] + (int)(LXC37 * log((nl + 1) / 30.)));
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
            int u = nl + ns;
            energy = (u <= MAXLOOPSIZE) ? ((*epm.internal_loop37)[u])
                                        : ((*epm.internal_loop37)[30] + (int)(LXC37 * log((u) / 30.)));

            energy += std::min(epm.MAX_NINIO, (nl - ns) * epm.ninio37);

            energy += (*epm.mismatchI37)[type][nuci1][nucj_1] + (*epm.mismatchI37)[type_2][nucq1][nucp_1];
        }
        }
        return energy;
    }

    // multi_loop
    int E_MLstem(int type, int si1, int sj1) {
        int energy = 0;

        if (*epm.mismatchM37) {
            if (si1 >= 0 && sj1 >= 0)
                energy += (*epm.mismatchM37)[type][si1][sj1];
            else if (si1 >= 0)
                energy += (*epm.dangle5_37)[type][si1];
            else if (sj1 >= 0)
                energy += (*epm.dangle3_37)[type][sj1];
        } else {
            if (si1 > 0)
                energy += (*epm.dangle5_37)[type][si1];
            if (sj1 > 0)
                energy += (*epm.dangle3_37)[type][sj1];
        }

        if (type > 2 || type == 0) {
            energy += epm.TerminalAU37;
        }

        energy += epm.ML_intern37;

        return energy;
    }

    int score_multi(int i, int j, int nuci, int nuci1, int nucj_1, int nucj, int len) {
        int tt = NUC_TO_PAIR(nucj, nuci);
        return E_MLstem(tt, nucj_1, nuci1) + epm.ML_closing37;
    }

    int score_M1(int i, int j, int k, int nuci_1, int nuci, int nuck, int nuck1, int len) {
        int tt = NUC_TO_PAIR(nuci, nuck);
        return E_MLstem(tt, nuci_1, nuck1);
    }

    int score_multi_unpaired(int i, int j) { return epm.ML_BASE37 * (j - i); }

    // exterior_loop
    int score_external_paired(int i, int j, int nuci_1, int nuci, int nucj, int nucj1, int len) {
        int type = NUC_TO_PAIR(nuci, nucj);
        int energy = 0;

        if (*epm.mismatchExt37) {
            if (nuci_1 >= 0 && nucj1 >= 0)
                energy += (*epm.mismatchExt37)[type][nuci_1][nucj1];
            else if (nuci_1 >= 0)
                energy += (*epm.dangle5_37)[type][nuci_1];
            else if (nucj1 >= 0)
                energy += (*epm.dangle3_37)[type][nucj1];
        } else {
            if (nuci_1 > 0)
                energy += (*epm.dangle5_37)[type][nuci_1];
            if (nucj1 > 0)
                energy += (*epm.dangle3_37)[type][nucj1];
        }

        if (type > 2 || type == 0)
            energy += epm.TerminalAU37;

        return energy;
    }

    int score_external_unpaired(int i, int j) { return 0; }
};

#endif // ENERGY_MODEL_HPP
