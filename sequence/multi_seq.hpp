#ifndef MULTI_SEQ_HPP
#define MULTI_SEQ_HPP

#include <vector>
#include <set>
#include "./seq.hpp"


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

    // access a sequence with bounds checking
    Seq& at(size_t index);
    const Seq& at(size_t index) const;

    //! Extracts all sequences from MultiSeq object whose index is given by a set. 
    //! Projects the multiple sequences to subset and returns as a new MultiSequence object.
    MultiSeq *Project(const std::set<int> &indices);

    // get all sequences
    const std::vector<Seq>& get_sequences() const;

    // get the number of sequences in the alignment
    size_t size() const;

    // get the length of the alignment
    size_t alignment_length() const;

    // find a sequence by ID
    Seq* find_sequence_by_id(const std::string& id);

    // compute pairwise sequence identity
    float get_seq_identity();

    // Read and Write functions
    bool read_fasta(const std::string& filepath);
    bool write_fasta(const std::string& filepath) const;
    bool write_fasta(std::ostream& out) const;

    // print all sequences
    void print_sequences();

    // Iterator support for range-based for loops
    inline std::vector<Seq>::iterator begin() { return sequences_.begin(); }
    inline std::vector<Seq>::iterator end() { return sequences_.end(); }

    inline std::vector<Seq>::const_iterator begin() const { return sequences_.begin(); }
    inline std::vector<Seq>::const_iterator end() const { return sequences_.end(); }

private:
    std::vector<Seq> sequences_;
};

#endif // MULTI_SEQ_HPP
