#ifndef MULTI_SEQ_HPP
#define MULTI_SEQ_HPP

#include <vector>
#include "Seq.hpp"


class MultiSeq {
public:
    // constructors
    MultiSeq();

    // add, insert, delete, replace sequences
    void add_sequence(const Seq& seq);
    void insert_sequence(size_t index, const Seq& seq);
    void delete_sequence(size_t index);
    void replace_sequence(size_t index, const Seq& seq);

    // access a sequence by index
    Seq& operator[](size_t index);
    const Seq& operator[](size_t index) const;

    // get all sequences
    const std::vector<Seq>& get_sequences() const;

    // get the number of sequences in the alignment
    size_t size() const;

    // get the length of the alignment
    size_t alignment_length() const;

    // find a sequence by ID
    Seq* find_sequence_by_id(const std::string& id);

    // Read and Write functions
    bool read_fasta(const std::string& filepath);
    bool write_fasta(const std::string& filepath) const;

private:
    std::vector<Seq> sequences_;
};

#endif // MULTI_SEQ_HPP
