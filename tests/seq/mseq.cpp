// linearx/tests/test_multi_seq.cpp
#include <catch2/catch_all.hpp>
#include <linearx/sequence/multi_sequence.hpp>
#include <linearx/utility.hpp>

TEST_CASE("MultiSeq basic add/access", "[msq][1]") {
    MultiSeq ms;
    Sequence s1("ACGU"), s2("UUUU");

    ms.add_sequence(s1);
    ms.add_sequence(s2);

    REQUIRE(ms.size() == 2);
    REQUIRE(ms[1].at(0) == 'U');
}

TEST_CASE("MultiSeq find and replace", "[msq][2]") {
    MultiSeq ms;
    ms.add_sequence(Sequence("ACGU", "s1", 1));
    ms.add_sequence(Sequence("UUUU", "s2", 2));

    auto found = ms.find_sequence_by_name("s1");
    REQUIRE(found != nullptr);
    REQUIRE(found->seq == "ACGU");

    ms.replace_sequence(1, Sequence("CCCC", "s3", 3));
    REQUIRE(ms[1].seq == "CCCC");
}

TEST_CASE("MultiSeq FASTA read/write", "[msq][3]") {
    MultiSeq ms;
    REQUIRE(ms.read_fasta("tests/data/s1.fasta", linearx::utils::VIENNA_NUC_ENCODING_SCHEME) == true);
    REQUIRE(ms.size() > 0);

    bool write_ok = ms.write_fasta("tests/tmp/output.fa", 60);
    REQUIRE(write_ok == true);

    // print the sequences
    ms.print(std::cout, false);
    ms.print(std::cout, true);
}

TEST_CASE("MultiSeq 10 Sequences FASTA read, ID and Sequence Identity", "[msq][4]") {
    MultiSeq ms;
    REQUIRE(ms.read_fasta("tests/data/s3.fasta", linearx::utils::VIENNA_NUC_ENCODING_SCHEME) == true);
    REQUIRE(ms.size() == 20);

    // check if the id is correctly assigned
    for (size_t i = 0; i < ms.size(); i++) {
        REQUIRE(ms[i].id == static_cast<int>(i));
    }

    // get average sequence identity
    std::cout << "Average sequence identity: " << ms.average_seq_identity() << std::endl;
}

TEST_CASE("MultiSeq insert_sequence function", "[msq][5]") {
    MultiSeq ms;
    ms.add_sequence(Sequence("AAAA", "seq1", 1));
    ms.add_sequence(Sequence("CCCC", "seq3", 3));

    // Insert in middle
    ms.insert_sequence(1, Sequence("BBBB", "seq2", 2));

    REQUIRE(ms.size() == 3);
    REQUIRE(ms[0].seq == "AAAA");
    REQUIRE(ms[0].name == "seq1");
    REQUIRE(ms[1].seq == "BBBB");
    REQUIRE(ms[1].name == "seq2");
    REQUIRE(ms[2].seq == "CCCC");
    REQUIRE(ms[2].name == "seq3");

    // Insert at end (beyond current size)
    ms.insert_sequence(10, Sequence("DDDD", "seq4", 4));

    REQUIRE(ms.size() == 3);

    // Insert at beginning
    ms.insert_sequence(0, Sequence("ZZZZ", "seq0", 0));

    REQUIRE(ms.size() == 5);
    REQUIRE(ms[0].seq == "ZZZZ");
    REQUIRE(ms[0].name == "seq0");
    REQUIRE(ms[1].seq == "AAAA");
    REQUIRE(ms[1].name == "seq1");
}

TEST_CASE("MultiSeq delete_sequence function", "[msq][6]") {
    MultiSeq ms;
    ms.add_sequence(Sequence("AAAA", "seq1", 1));
    ms.add_sequence(Sequence("BBBB", "seq2", 2));
    ms.add_sequence(Sequence("CCCC", "seq3", 3));
    ms.add_sequence(Sequence("DDDD", "seq4", 4));

    REQUIRE(ms.size() == 4);

    // Delete from middle
    ms.delete_sequence(1);

    REQUIRE(ms.size() == 3);
    REQUIRE(ms[0].seq == "AAAA");
    REQUIRE(ms[0].name == "seq1");
    REQUIRE(ms[1].seq == "CCCC");
    REQUIRE(ms[1].name == "seq3");
    REQUIRE(ms[2].seq == "DDDD");
    REQUIRE(ms[2].name == "seq4");

    // Delete from beginning
    ms.delete_sequence(0);

    REQUIRE(ms.size() == 2);
    REQUIRE(ms[0].seq == "CCCC");
    REQUIRE(ms[0].name == "seq3");
    REQUIRE(ms[1].seq == "DDDD");
    REQUIRE(ms[1].name == "seq4");

    // Try to delete invalid index (should not change anything)
    ms.delete_sequence(10);  // out of bounds

    REQUIRE(ms.size() == 2);
}

TEST_CASE("MultiSeq Project function", "[msq][7]") {
    MultiSeq ms;

    // Create aligned sequences with some gaps
    // Position:  0123456789
    // seq1:      A-C-G-T-A-
    // seq2:      A-C-G---A-
    // seq3:      A---G-T-A-
    // seq4:      A-C-G-T---
    ms.add_sequence(Sequence("A-C-G-T-A-", "seq1", 1));
    ms.add_sequence(Sequence("A-C-G---A-", "seq2", 2));
    ms.add_sequence(Sequence("A---G-T-A-", "seq3", 3));
    ms.add_sequence(Sequence("A-C-G-T---", "seq4", 4));

    REQUIRE(ms.size() == 4);

    // Project sequences 0, 2, 3 (seq1, seq3, seq4)
    // Gap-only columns are positions 3,6,7,9 (0-indexed)
    // So positions 0,1,2,4,5,8 should remain
    // seq1: A-C-G-T-A- -> ACGTA (positions 0,2,4,6,8)
    // seq3: A---G-T-A- -> A-GTA  (positions 0,4,6,8)
    // seq4: A-C-G-T--- -> ACGT-  (positions 0,2,4,6)
    std::set<int> indices = {0, 2, 3};
    MultiSeq* projected = ms.Project(indices);

    REQUIRE(projected != nullptr);
    REQUIRE(projected->size() == 3);

    // Check the projected sequences
    REQUIRE((*projected)[0].seq == "ACGTA");
    REQUIRE((*projected)[0].name == "seq1");
    REQUIRE((*projected)[1].seq == "A-GTA");
    REQUIRE((*projected)[1].name == "seq3");
    REQUIRE((*projected)[2].seq == "ACGT-");
    REQUIRE((*projected)[2].name == "seq4");

    delete projected;

    // Test with all sequences - should remove only all-gap columns
    std::set<int> all_indices = {0, 1, 2, 3};
    MultiSeq* projected_all = ms.Project(all_indices);

    REQUIRE(projected_all->size() == 4);
    // Positions 3,6,7,9 are all gaps in all sequences, so they should be removed
    // Remaining positions: 0,1,2,4,5,8
    REQUIRE((*projected_all)[0].seq == "ACGTA");  // A-C-G-T-A-
    REQUIRE((*projected_all)[1].seq == "ACG-A");  // A-C-G---A-
    REQUIRE((*projected_all)[2].seq == "A-GTA");  // A---G-T-A-
    REQUIRE((*projected_all)[3].seq == "ACGT-");  // A-C-G-T---

    delete projected_all;
}
