// linearx/tests/test_multi_seq.cpp
#include <catch2/catch_all.hpp>
#include <linearx/turbofold/linear_turbofold.hpp>
#include <linearx/utility.hpp>

TEST_CASE("Simple 5 Sequence", "[tf][1]") {
    MultiSeq ms;
    ms.read_fasta("tests/data/s4.fasta", linearx::utils::VIENNA_NUC_ENCODING_SCHEME);
    EnergyModel energy_model(EnergyParamsType::VIENNA);
    LinearTurbofold ltf(ms, energy_model);
    ltf.run(3);
}
