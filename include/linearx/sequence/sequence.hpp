// linearx/sequence/seq.hpp
#pragma once
#include <cstddef>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

class Sequence {
   public:
    std::string name;      // human-readable name
    int id;                // numeric ID
    std::string seq;       // sequence in uppercase
    std::vector<int> enc;  // encoded sequence (optional)

    // constructors
    Sequence();  // default
    explicit Sequence(const std::string& seq);
    Sequence(const std::string& seq, const std::string& name, int id);
    Sequence(const std::string& seq, const std::string& name, int id, const std::unordered_map<char, int>& encoding_scheme,
        bool randomize_N = false);

    // encoding
    void randomize_N();  // replace 'N' with random nucleotide
    void set_encoding(const std::unordered_map<char, int>& encoding_scheme);

    // length and access
    size_t length() const;
    char at(size_t index) const;
    char& operator[](size_t index);
    const char& operator[](size_t index) const;

    // comparison
    bool equals(const Sequence& other) const;

    // sequence utilities
    void remove_gaps();
    void swap_nuc(char nuc1, char nuc2);
    float compute_seq_identity(const Sequence& seq2) const;

    // edit operations
    void add_nuc(char nucleotide);
    void insert_nuc(size_t index, char nucleotide);
    void delete_nuc(size_t index);
    void replace_nuc(size_t index, char nucleotide);
    void reverse();

    // file operations
    bool read_fasta(const std::string& filepath);
    bool write_fasta(const std::string& filepath) const;

    // printing and details
    void show_details(std::ostream& os = std::cout) const;
    void print(std::ostream& os = std::cout, bool use_enc = false) const;

    // write a sequence to a FASTA file (optionally wrap lines)
    bool write_fasta(const std::string& filepath, size_t max_line_length = 0) const;
};
