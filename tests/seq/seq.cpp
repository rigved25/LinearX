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

TEST_CASE("Sequence clone function", "[sq][5]") {
    Sequence original("ACGU", "test_seq", 42);
    std::unordered_map<char, int> scheme{{'A', 0}, {'C', 1}, {'G', 2}, {'U', 3}};
    original.set_encoding(scheme);

    Sequence* cloned = original.clone();

    REQUIRE(cloned != nullptr);
    REQUIRE(cloned != &original);  // Different memory location
    REQUIRE(cloned->seq == original.seq);
    REQUIRE(cloned->name == original.name);
    REQUIRE(cloned->id == original.id);
    REQUIRE(cloned->enc == original.enc);
    REQUIRE(cloned->length() == original.length());

    // Test that modifying clone doesn't affect original
    cloned->replace_nuc(0, 'U');
    REQUIRE(original.at(0) == 'A');
    REQUIRE(cloned->at(0) == 'U');

    delete cloned;  // Clean up
}

TEST_CASE("Sequence add_gaps function", "[sq][6]") {
    Sequence original("ATGCCGTCA", "test_seq", 1);
    std::string alignment = "XXXBBYYYBBYYXX";  // 14 characters
    char id = 'X';

    Sequence* gapped = original.add_gaps(&alignment, id);

    REQUIRE(gapped != nullptr);
    REQUIRE(gapped->name == original.name);
    REQUIRE(gapped->id == original.id);
    REQUIRE(gapped->length() == alignment.length());

    // According to the comment: "ATGCC---GT--CA"
    // Let's verify this matches the expected pattern
    REQUIRE(gapped->seq == "ATGCC---GT--CA");

    // Test with different alignment pattern
    std::string alignment2 = "ABCD";  // 4 characters
    Sequence original2("AT", "test2", 2);
    Sequence* gapped2 = original2.add_gaps(&alignment2, 'A');

    REQUIRE(gapped2->seq == "AT--");  // A from original, then gaps, then T

    delete gapped;
    delete gapped2;
}

TEST_CASE("Sequence get_mapping function", "[sq][7]") {
    // Test with gapped sequence
    Sequence gapped("ATGCC---GT--CA", "gapped_seq", 1);
    std::vector<int>* mapping = gapped.get_mapping();

    REQUIRE(mapping != nullptr);
    REQUIRE(mapping->size() == 10);  // 10 non-gap positions
    REQUIRE((*mapping)[0] == 0);    // starts with 0
    REQUIRE((*mapping)[1] == 0);    // position 0 again
    REQUIRE((*mapping)[2] == 1);    // position 1
    REQUIRE((*mapping)[3] == 2);    // position 2
    REQUIRE((*mapping)[4] == 3);    // position 3
    REQUIRE((*mapping)[5] == 4);    // position 4
    REQUIRE((*mapping)[6] == 8);    // position 8
    REQUIRE((*mapping)[7] == 9);    // position 9
    REQUIRE((*mapping)[8] == 12);   // position 12
    REQUIRE((*mapping)[9] == 13);   // position 13

    delete mapping;

    // Test with no gaps
    Sequence no_gaps("ATGC", "no_gaps", 2);
    std::vector<int>* mapping2 = no_gaps.get_mapping();

    REQUIRE(mapping2->size() == 5);  // {0,0,1,2,3}
    REQUIRE((*mapping2)[0] == 0);
    REQUIRE((*mapping2)[1] == 0);
    REQUIRE((*mapping2)[2] == 1);
    REQUIRE((*mapping2)[3] == 2);
    REQUIRE((*mapping2)[4] == 3);

    delete mapping2;
}

TEST_CASE("Sequence reverse function", "[sq][8]") {
    Sequence s("ATGC", "test", 1);
    s.reverse();

    REQUIRE(s.seq == "CGTA");
    REQUIRE(s.length() == 4);

    // Test with gapped sequence
    Sequence gapped("AT-GC", "gapped", 2);
    gapped.reverse();

    REQUIRE(gapped.seq == "CG-TA");
    REQUIRE(gapped.length() == 5);

    // Test with single character
    Sequence single("A", "single", 3);
    single.reverse();

    REQUIRE(single.seq == "A");
    REQUIRE(single.length() == 1);
}

TEST_CASE("Sequence read_fasta_stream function", "[sq][9]") {
    // Test successful parsing
    std::stringstream ss(">test_sequence\nATGCATGC\n");
    Sequence s;
    bool success = s.read_fasta_stream(ss);

    REQUIRE(success == true);
    REQUIRE(s.seq == "ATGCATGC");
    REQUIRE(s.name == "test_sequence");

    // Test with whitespace and case conversion
    std::stringstream ss2(">seq2\n a T g C \n");
    Sequence s2;
    success = s2.read_fasta_stream(ss2);

    REQUIRE(success == true);
    REQUIRE(s2.seq == "ATGC");
    REQUIRE(s2.name == "seq2");

    // Test with no sequence data
    std::stringstream ss3(">empty\n");
    Sequence s3;
    success = s3.read_fasta_stream(ss3);

    REQUIRE(success == false);
    REQUIRE(s3.seq.empty());
}

TEST_CASE("Sequence show_details function", "[sq][10]") {
    Sequence s("ATGC", "test_seq", 42);

    std::stringstream ss;
    s.show_details(ss);

    std::string output = ss.str();
    REQUIRE(output.find("test_seq") != std::string::npos);
    REQUIRE(output.find("42") != std::string::npos);
    REQUIRE(output.find("4") != std::string::npos);  // length
    REQUIRE(output.find("[Sequence Info]") != std::string::npos);
}
