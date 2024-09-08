#include <cassert>
#include <iostream>

#include "./../../Seq.hpp"
#include "./../../MultiSeq.hpp"

void test_sequence() {
    // Test default constructor and setters
    Seq seq1;
    seq1.set_id("seq1");
    seq1.set_sequence("AUGCUAGCUA");
    assert(seq1.get_id() == "seq1");
    assert(seq1.get_sequence() == "AUGCUAGCUA");

    // Test parameterized constructor
    Seq seq2("seq2", "CUGGAUCGUA");
    assert(seq2.get_id() == "seq2");
    assert(seq2.get_sequence() == "CUGGAUCGUA");

    // Test length
    assert(seq1.length() == 10);

    // Test character access with []
    assert(seq1[0] == 'A');
    assert(seq1[9] == 'A');

    // Test at() method
    assert(seq1.at(2) == 'G');

    // Test add_char
    seq1.add_nuc('U');
    assert(seq1.get_sequence() == "AUGCUAGCUAU");

    // Test insert_char
    seq1.insert_nuc(1, 'G');
    assert(seq1.get_sequence() == "AGUGCUAGCUAU");

    // Test delete_char
    seq1.delete_nuc(1);
    assert(seq1.get_sequence() == "AUGCUAGCUAU");

    // Test replace_char
    seq1.replace_nuc(0, 'U');
    assert(seq1.get_sequence() == "UUGCUAGCUAU");

    // Test read_from_fasta
    Seq seq3;
    bool read_success = seq3.read_fasta("./../test_data/test_seq.fasta");
    assert(read_success);
    assert(seq3.get_id() == "test_seq");
    assert(seq3.get_sequence() == "AUGCUAGCUA");

    // Test write_to_fasta
    bool write_success = seq3.write_fasta("output_seq.fasta");
    assert(write_success);

    std::cout << "All Seq tests passed!" << std::endl;
}

void test_seq_alignment() {
    // Test default constructor
    MultiSeq alignment;

    // Test add_sequence
    Seq seq1("seq1", "AUGCUAGCUA");
    Seq seq2("seq2", "CUGGAUCGUA");
    alignment.add_sequence(seq1);
    alignment.add_sequence(seq2);
    assert(alignment.size() == 2);

    // Test insert_sequence
    Seq seq3("seq3", "GGGGAUCGUA");
    alignment.insert_sequence(1, seq3);
    assert(alignment.size() == 3);
    assert(alignment[1].get_id() == "seq3");

    // Test delete_sequence
    alignment.delete_sequence(1);
    assert(alignment.size() == 2);
    assert(alignment[1].get_id() == "seq2");

    // Test replace_sequence
    alignment.replace_sequence(1, seq3);
    assert(alignment[1].get_id() == "seq3");

    // Test operator[]
    assert(alignment[0].get_sequence() == "AUGCUAGCUA");

    // Test alignment_length
    assert(alignment.alignment_length() == 10);

    // Test find_sequence_by_id
    Seq* found_seq = alignment.find_sequence_by_id("seq1");
    assert(found_seq != nullptr);
    assert(found_seq->get_id() == "seq1");

    // Test read_from_fasta
    MultiSeq alignment2;
    bool read_success = alignment2.read_fasta("./../test_data/test_aln.fasta");
    assert(read_success);
    assert(alignment2.size() == 2);
    assert(alignment2[0].get_id() == "test_seq1");

    // Test write_to_fasta
    bool write_success = alignment2.write_fasta("output_alignment.fasta");
    assert(write_success);

    std::cout << "All MultiSeq tests passed!" << std::endl;
}

int main() {
    test_sequence();
    test_seq_alignment();

    return 0;
}
