#include "./multi_seq.hpp"
#include <fstream>
#include <iostream>


// constructor
MultiSeq::MultiSeq() = default;

// add a sequence to the alignment
void MultiSeq::add_sequence(const Seq& seq) {
    sequences_.push_back(seq);
}

// insert a sequence at a specific index
void MultiSeq::insert_sequence(size_t index, const Seq& seq) {
    if (index <= sequences_.size()) {
        sequences_.insert(sequences_.begin() + index, seq);
    } else {
        std::cerr << "Warning: Invalid index " << index << " for insert_sequence. No action taken." << std::endl;
    }
}

// delete a sequence at a specific index
void MultiSeq::delete_sequence(size_t index) {
    if (index < sequences_.size()) {
        sequences_.erase(sequences_.begin() + index);
    } else {
        std::cerr << "Warning: Invalid index " << index << " for delete_sequence. No action taken." << std::endl;
    }
}

// replace a sequence at a specific index
void MultiSeq::replace_sequence(size_t index, const Seq& seq) {
    if (index < sequences_.size()) {
        sequences_[index] = seq;
    } else {
        std::cerr << "Warning: Invalid index " << index << " for replace_sequence. No action taken." << std::endl;
    }
}

// access a sequence by index
Seq& MultiSeq::operator[](size_t index) {
    if (index >= sequences_.size()) {
        std::cerr << "Warning: Invalid index " << index << " for operator[]. Returning the first sequence." << std::endl;
        return sequences_.at(0);
    }
    return sequences_.at(index);
}

const Seq& MultiSeq::operator[](size_t index) const {
    if (index >= sequences_.size()) {
        std::cerr << "Warning: Invalid index " << index << " for operator[]. Returning the first sequence." << std::endl;
        return sequences_.at(0);
    }
    return sequences_.at(index);
}

// get all sequences
const std::vector<Seq>& MultiSeq::get_sequences() const {
    return sequences_;
}

// get the number of sequences in the alignment
size_t MultiSeq::size() const {
    return sequences_.size();
}

// get the length of the alignment
size_t MultiSeq::alignment_length() const {
    if (sequences_.empty()) return 0;
    return sequences_.front().length();
}

// find a sequence by ID
Seq* MultiSeq::find_sequence_by_id(const std::string& id) {
    for (auto& seq : sequences_) {
        if (seq.get_id() == id) {
            return &seq;
        }
    }
    return nullptr;
}

// read sequences from a FASTA file and create an alignment
bool MultiSeq::read_fasta(const std::string& filepath) {
    std::ifstream infile(filepath);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << filepath << std::endl;
        return false;
    }

    std::string line;
    std::string id;
    std::string sequence;
    
    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!id.empty()) {
                // if we already have an ID, add the previous sequence to the alignment
                sequences_.emplace_back(id, sequence);
                sequence.clear();
            }
            id = line.substr(1); // remove the '>' character
        } else {
            sequence += line;
        }
    }

    // add the last sequence
    if (!id.empty() && !sequence.empty()) {
        sequences_.emplace_back(id, sequence);
    }

    infile.close();
    return !sequences_.empty();
}

// write all sequences in the alignment to a FASTA file
bool MultiSeq::write_fasta(const std::string& filepath) const {
    std::ofstream outfile(filepath);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filepath << " for writing" << std::endl;
        return false;
    }

    for (const auto& seq : sequences_) {
        outfile << ">" << seq.get_id() << "\n" << seq.get_sequence() << "\n";
    }

    outfile.close();
    return true;
}