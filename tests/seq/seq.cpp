// linearx/tests/test_seq.cpp
#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include <linearx/sequence/sequence.hpp>

using Catch::Approx;

TEST_CASE("Sequence basic construction and access", "[sq][1]"){
    Sequence s("ACGUACGU");

    REQUIRE(s.length() == 8);
    REQUIRE(s.at(0) == 'A');
    REQUIRE(s[3] == 'U');
    REQUIRE(s.at(5) == 'C');
    REQUIRE(s.id == 0);
}

TEST_CASE("Sequence add/insert/delete/replace nucleotides", "[sq][2]") {
    Sequence s("ACG");

    s.add_nuc('U');
    REQUIRE(s.length() == 4);
    REQUIRE(s[3] == 'U');

    s.insert_nuc(2, 'U');
    REQUIRE(s[2] == 'U');

    s.replace_nuc(0, 'G');
    REQUIRE(s[0] == 'G');

    s.delete_nuc(1);
    REQUIRE(s.length() == 4);

}

TEST_CASE("Sequence identity and gaps", "[sq][3]") {
    Sequence s1("ACGU");
    Sequence s2("ACGU");
    REQUIRE(s1.compute_seq_identity(s2) == Approx(1.0f));

    s1.add_nuc('-');
    s1.remove_gaps();
    REQUIRE(s1.length() == 4);
}

TEST_CASE("Sequence encoding and error handling", "[sq][4]") {
    std::unordered_map<char, int> scheme{{'A', 0}, {'C', 1}, {'G', 2}, {'U', 3}};
    Sequence s("ACGU");
    s.set_encoding(scheme);

    REQUIRE(s.enc.size() == 4);
    REQUIRE(s.enc[2] == 2);

    // print the sequence
    s.print(std::cout, false);
    s.print(std::cout, true);
}
