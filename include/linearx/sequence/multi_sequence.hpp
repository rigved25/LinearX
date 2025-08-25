// linearx/sequence/multi_seq.hpp
#pragma once
#include <linearx/sequence/sequence.hpp>
#include <string>
#include <unordered_map>
#include <vector>

class MultiSeq {
   public:
    // constructors
    MultiSeq();
    MultiSeq(const std::vector<Sequence>& sequences);

    // add, insert, delete, replace sequences
    void add_sequence(const Sequence& seq);
    void insert_sequence(size_t index, const Sequence& seq);
    void delete_sequence(size_t index);
    void replace_sequence(size_t index, const Sequence& seq);

    // access a sequence by index
    Sequence& operator[](size_t index);
    const Sequence& operator[](size_t index) const;

    // access a sequence with bounds checking
    Sequence& at(size_t index);
    const Sequence& at(size_t index) const;

    // get all sequences
    const std::vector<Sequence>& get_sequences() const;

    // get the number of sequences
    size_t size() const;

    // get the length of the alignment (assumes all sequences are same length)
    size_t alignment_length() const;

    // find a sequence by ID
    Sequence* find_sequence_by_id(const int id);

    // find a sequence by name
    Sequence* find_sequence_by_name(const std::string& name);

    // compute pairwise sequence identity between two indices
    float compute_seq_identity(size_t i, size_t j) const;

    // compute average pairwise sequence identity
    float average_seq_identity() const;

    // read sequences from multi-FASTA file
    bool read_fasta(const std::string& filepath, const std::unordered_map<char, int>& encoding_scheme = {},
                    bool randomize_N = false);

    // write sequences to multi-FASTA file
    bool write_fasta(const std::string& filepath, int max_line_length = 0) const;

    // print all sequences
    void print(std::ostream& out = std::cout, bool use_enc = false) const;

    // range-based for loop support
    inline std::vector<Sequence>::iterator begin() { return sequences_.begin(); }
    inline std::vector<Sequence>::iterator end() { return sequences_.end(); }

    inline std::vector<Sequence>::const_iterator begin() const { return sequences_.begin(); }
    inline std::vector<Sequence>::const_iterator end() const { return sequences_.end(); }

   private:
    std::vector<Sequence> sequences_;
};
