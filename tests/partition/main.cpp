// linearx/tests/alignment/main.cpp
#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include <linearx/partition/linear_partition.hpp>

using namespace std;

TEST_CASE("[LinearPartition] MFE Basic Seq", "[mfe][1]") {
    Sequence seq("CCCAAAGGG");
    seq.set_encoding(linearx::utils::VIENNA_NUC_ENCODING_SCHEME);

    EnergyModel energy_model(EnergyParamsType::VIENNA);
    LinearPartition partition(seq, energy_model);
    partition.compute_inside(Mode::BEST, 100, true);
    Structure mfe_structure = partition.get_mfe_structure();
    std::cout << "MFE Structure: " << mfe_structure.getDotBracket() << std::endl;
}

TEST_CASE("[LinearPartition] Partition Basic Seq", "[part][1]") {
    Sequence seq("CCCAAAGGG");
    seq.set_encoding(linearx::utils::VIENNA_NUC_ENCODING_SCHEME);
    seq.print(std::cout, true);

    EnergyModel energy_model(EnergyParamsType::VIENNA);
    LinearPartition partition(seq, energy_model);
    partition.compute_inside(Mode::PARTITION, 100, true);
    partition.print_alpha_beta();
    Structure mfe_structure = partition.get_mfe_structure();
    partition.debug_states();
}

TEST_CASE("[LinearPartition] Big Seq", "[mfe][2]") {
    Sequence seq;
    seq.read_fasta("tests/data/s2.fasta");
    seq.set_encoding(linearx::utils::VIENNA_NUC_ENCODING_SCHEME);

    EnergyModel energy_model(EnergyParamsType::VIENNA);
    LinearPartition partition(seq, energy_model);
    partition.compute_inside(Mode::BEST, 100, true);
    Structure mfe_structure = partition.get_mfe_structure();
    std::cout << "MFE Structure: " << mfe_structure.getDotBracket() << std::endl;
}

TEST_CASE("[LinearPartition] Big Seq", "[part][2]") {
    Sequence seq;
    seq.read_fasta("tests/data/s2.fasta");
    seq.set_encoding(linearx::utils::VIENNA_NUC_ENCODING_SCHEME);

    EnergyModel energy_model(EnergyParamsType::VIENNA);
    LinearPartition partition(seq, energy_model);
    partition.compute_inside(Mode::PARTITION, 100, true);
    partition.compute_outside(linearx::constants::limits::DEVIATION_THRESHOLD, true);
    // partition.compute_outside(linearx::constants::math::VALUE_MAX * 100, true);
}

TEST_CASE("[LinearPartition] Covid Seq", "[mfe][3]") {
    Sequence seq;
    seq.read_fasta("tests/data/cvd_hbref.fasta");
    seq.set_encoding(linearx::utils::VIENNA_NUC_ENCODING_SCHEME);

    EnergyModel energy_model(EnergyParamsType::VIENNA);
    LinearPartition partition(seq, energy_model);
    partition.compute_inside(Mode::BEST, 100, true);
}

TEST_CASE("[LinearPartition] Covid Seq", "[part][3]") {
    Sequence seq;
    seq.read_fasta("tests/data/cvd_hbref.fasta");
    seq.set_encoding(linearx::utils::VIENNA_NUC_ENCODING_SCHEME);

    EnergyModel energy_model(EnergyParamsType::VIENNA);
    LinearPartition partition(seq, energy_model);
    partition.compute_inside(Mode::PARTITION, 100, true);
}

TEST_CASE("[LinearPartition] Covid Seq", "[part][out][3]") {
    Sequence seq;
    seq.read_fasta("tests/data/cvd_hbref.fasta");
    seq.set_encoding(linearx::utils::VIENNA_NUC_ENCODING_SCHEME);

    EnergyModel energy_model(EnergyParamsType::VIENNA);
    LinearPartition partition(seq, energy_model);
    partition.compute_inside(Mode::PARTITION, 100, true);
    // partition.compute_outside(linearx::constants::limits::DEVIATION_THRESHOLD, true);
    partition.compute_outside(linearx::constants::math::VALUE_MAX * 100, true);
}
