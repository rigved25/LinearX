#ifndef SEQ_HPP
#define SEQ_HPP

#include <string>
#include <unordered_map>
#include <vector>

class Seq {
  public:
    std::string id;
    std::string sequence;
    int k_id;                 // optional integer id for the sequence
    std::vector<int> enc_seq; // optional encoded sequence

    // constructors
    Seq();
    Seq(const std::string &id, const std::string &sequence, const int k_id = 0,
        std::unordered_map<char, int> *encoding_scheme = nullptr);

    // encoding (encodes the sequence given a encoding scheme in dictionary)
    void set_encoding(std::unordered_map<char, int> &encoding_scheme);

    // additional functionalities
    size_t length() const;
    char at(size_t index) const;
    char &operator[](size_t index);
    const char &operator[](size_t index) const;

    void remove_gaps();
    void swap_nuc(const char nuc1, const char nuc2); // swap all occurences of nuc1 with nuc2
    float compute_seq_identity(const Seq &seq2);      // compute sequence identity with another sequence

    // add, insert, delete, replace characters
    void add_nuc(char nucleotide);
    void insert_nuc(size_t index, char nucleotide);
    void delete_nuc(size_t index);
    void replace_nuc(size_t index, char nucleotide);

    // read and Write functions
    bool read_fasta(const std::string &filepath);
    bool write_fasta(const std::string &filepath) const;
};

#endif // SEQ_HPP
