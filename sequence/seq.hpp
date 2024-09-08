#ifndef SEQ_HPP
#define SEQ_HPP

#include <string>
#include <unordered_map>
#include <vector>


class Seq {
public:
    // constructors
    Seq();
    Seq(const std::string& id, const std::string& sequence);

    // accessors
    const std::string& get_id() const;
    const std::string& get_sequence() const;

    // mutators
    void set_id(const std::string& id);
    void set_sequence(const std::string& sequence);

    // additional functionalities
    size_t length() const;
    char at(size_t index) const;
    char& operator[](size_t index);
    const char& operator[](size_t index) const;

    void remove_gaps();
    void swap_nuc(const char nuc1, const char nuc2);  // swap all occurences of nuc1 with nuc2

    // add, insert, delete, replace characters
    void add_nuc(char nucleotide);
    void insert_nuc(size_t index, char nucleotide);
    void delete_nuc(size_t index);
    void replace_nuc(size_t index, char nucleotide);

    // read and Write functions
    bool read_fasta(const std::string& filepath);
    bool write_fasta(const std::string& filepath) const;

    // static methods
    static std::vector<int> encode(const Seq &seq, std::unordered_map<char, int> scheme);  // encodes the sequence given a encoding scheme in dictionary

private:
    std::string id_;
    std::string sequence_;
};

#endif // SEQ_HPP
