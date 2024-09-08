#ifndef ENERGY_PARAMS_HPP
#define ENERGY_PARAMS_HPP

#include "./BLParams.hpp"
#include "./ViennaParams.hpp"

enum EnergyParamsType { BL_STAR, VIENNA };

class EnergyParams {

public:
  int ML_intern37;
  int ML_closing37;
  int ML_BASE37;
  int MAX_NINIO;
  int ninio37;
  int TerminalAU37;

  const int (*hairpin37)[31] = nullptr;
  const int (*stack37)[8][8] = nullptr;
  const int (*bulge37)[31] = nullptr;
  const int (*internal_loop37)[31] = nullptr;

  const int (*mismatchH37)[8][5][5] = nullptr;
  const int (*mismatchM37)[8][5][5] = nullptr;
  const int (*mismatchExt37)[8][5][5] = nullptr;
  const int (*mismatchI37)[8][5][5] = nullptr;
  const int (*mismatch1nI37)[8][5][5] = nullptr;
  const int (*mismatch23I37)[8][5][5] = nullptr;

  const int (*dangle5_37)[8][5] = nullptr;
  const int (*dangle3_37)[8][5] = nullptr;

  const int (*int11_37)[8][8][5][5] = nullptr;
  const int (*int21_37)[8][8][5][5][5] = nullptr;
  const int (*int22_37)[8][8][5][5][5][5] = nullptr;

  const int (*Triloop37)[2] = nullptr;
  const int (*Tetraloop37)[16] = nullptr;
  const int (*Hexaloop37)[4] = nullptr;

  // constructor
  EnergyParams(EnergyParamsType model) {
    if (model == BL_STAR) {
      init_bl_params();
    } else if (model == VIENNA) {
      init_vienna_params();
    }
  }

private:
  void init_bl_params() {
    ML_intern37 = BLParams::ML_intern37;
    ML_closing37 = BLParams::ML_closing37;
    ML_BASE37 = BLParams::ML_BASE37;
    MAX_NINIO = BLParams::MAX_NINIO;
    ninio37 = BLParams::ninio37;
    TerminalAU37 = BLParams::TerminalAU37;

    hairpin37 = &BLParams::hairpin37;
    stack37 = &BLParams::stack37;
    bulge37 = &BLParams::bulge37;
    internal_loop37 = &BLParams::internal_loop37;

    mismatchH37 = &BLParams::mismatchH37;
    mismatchI37 = &BLParams::mismatchI37;

    dangle5_37 = &BLParams::dangle5_37;
    dangle3_37 = &BLParams::dangle3_37;

    int11_37 = &BLParams::int11_37;
    int21_37 = &BLParams::int21_37;
    int22_37 = &BLParams::int22_37;
  }

  void init_vienna_params() {
    ML_intern37 = ViennaParams::ML_intern37;
    ML_closing37 = ViennaParams::ML_closing37;
    ML_BASE37 = ViennaParams::ML_BASE37;
    MAX_NINIO = ViennaParams::MAX_NINIO;
    ninio37 = ViennaParams::ninio37;
    TerminalAU37 = ViennaParams::TerminalAU37;

    hairpin37 = &ViennaParams::hairpin37;
    stack37 = &ViennaParams::stack37;
    bulge37 = &ViennaParams::bulge37;
    internal_loop37 = &ViennaParams::internal_loop37;

    mismatchH37 = &ViennaParams::mismatchH37;
    mismatchM37 = &ViennaParams::mismatchM37;
    mismatchExt37 = &ViennaParams::mismatchExt37;
    mismatchI37 = &ViennaParams::mismatchI37;
    mismatch1nI37 = &ViennaParams::mismatch1nI37;
    mismatch23I37 = &ViennaParams::mismatch23I37;

    dangle5_37 = &ViennaParams::dangle5_37;
    dangle3_37 = &ViennaParams::dangle3_37;

    int11_37 = &ViennaParams::int11_37;
    int21_37 = &ViennaParams::int21_37;
    int22_37 = &ViennaParams::int22_37;

    Triloop37 = &ViennaParams::Triloop37;
    Tetraloop37 = &ViennaParams::Tetraloop37;
    Hexaloop37 = &ViennaParams::Hexaloop37;
  }
};

#endif // ENERGY_PARAMS_HPP
