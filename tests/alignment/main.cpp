// linearx/tests/alignment/main.cpp
#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include <linearx/alignment/linear_align.hpp>

TEST_CASE("LinearAlignment basic functionality", "[align][1]") {
    MultiSeq msa;
    msa.read_fasta("tests/data/s2.fasta", linearx::utils::VIENNA_NUC_ENCODING_SCHEME);

    LinearAlignment la(msa[0], msa[1]);
    la.use_prob_set1();
    la.compute_inside(Mode::BEST);
    la.print_alpha_beta();

    auto aln = la.get_alignment();
    // aln.print(std::cout, false);
    const float similarity = aln.average_seq_identity();
    std::cout << "\nAverage sequence identity: " << similarity << std::endl << std::endl;

    la.reset_beams();
    la.use_prob_set2(similarity);
    la.compute_inside(Mode::PARTITION);
    la.compute_outside(true);
    // la.compute_outside(true, std::numeric_limits<value_type>::max());
    // la.compute_outside(false);
    la.compute_coincidence_probabilities();
}
