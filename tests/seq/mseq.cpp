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
    REQUIRE(ms.read_fasta("tests/data/s1.fasta", linearx::utils::VIENNA_NUC_ENCODING_SCHEME, true) == true);
    REQUIRE(ms.size() > 0);

    bool write_ok = ms.write_fasta("tests/tmp/output.fa", 60);
    REQUIRE(write_ok == true);

    // print the sequences
    ms.print(std::cout, false);
    ms.print(std::cout, true);
}

TEST_CASE("MultiSeq 10 Sequences FASTA read, ID and Sequence Identity", "[msq][4]") {
    MultiSeq ms;
    REQUIRE(ms.read_fasta("tests/data/s3.fasta", linearx::utils::VIENNA_NUC_ENCODING_SCHEME, true) == true);
    REQUIRE(ms.size() == 20);

    // check if the id is correctly assigned
    for (size_t i = 0; i < ms.size(); i++) {
        REQUIRE(ms[i].id == static_cast<int>(i));
    }

    // get average sequence identity
    std::cout << "Average sequence identity: " << ms.average_seq_identity() << std::endl;
}
