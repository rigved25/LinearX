// linearx/sequence/multi_seq.cpp
#include <algorithm>
#include <fstream>
#include <iostream>
#include <linearx/sequence/multi_sequence.hpp>
#include <sstream>
#include <stdexcept>

// default constructor
MultiSeq::MultiSeq() = default;

// construct from vector of sequences
MultiSeq::MultiSeq(const std::vector<Sequence>& sequences) : sequences_(sequences) {}

// add a sequence
void MultiSeq::add_sequence(const Sequence& seq) { sequences_.push_back(seq); }

void MultiSeq::insert_sequence(size_t index, const Sequence& seq) {
    if (index <= sequences_.size()) {
        sequences_.insert(sequences_.begin() + index, seq);
    } else {
        std::cerr << "Warning: Invalid index " << index << " for insert_sequence.\n";
    }
}

void MultiSeq::delete_sequence(size_t index) {
    if (index < sequences_.size()) {
        sequences_.erase(sequences_.begin() + index);
    } else {
        std::cerr << "Warning: Invalid index " << index << " for delete_sequence.\n";
    }
}

void MultiSeq::replace_sequence(size_t index, const Sequence& seq) {
    if (index < sequences_.size()) {
        sequences_[index] = seq;
    } else {
        std::cerr << "Warning: Invalid index " << index << " for replace_sequence.\n";
    }
}

Sequence& MultiSeq::operator[](size_t index) {
    if (index >= sequences_.size()) {
        std::cerr << "Warning: Index " << index << " out of range for operator[]. Returning first.\n";
        return sequences_.at(0);
    }
    return sequences_.at(index);
}

const Sequence& MultiSeq::operator[](size_t index) const {
    if (index >= sequences_.size()) {
        std::cerr << "Warning: Index " << index << " out of range for operator[]. Returning first.\n";
        return sequences_.at(0);
    }
    return sequences_.at(index);
}

Sequence& MultiSeq::at(size_t index) {
    if (index >= sequences_.size()) throw std::out_of_range("Index out of range");
    return sequences_.at(index);
}

const Sequence& MultiSeq::at(size_t index) const {
    if (index >= sequences_.size()) throw std::out_of_range("Index out of range");
    return sequences_.at(index);
}

const std::vector<Sequence>& MultiSeq::get_sequences() const { return sequences_; }

size_t MultiSeq::size() const { return sequences_.size(); }

size_t MultiSeq::alignment_length() const { return sequences_.empty() ? 0 : sequences_[0].length(); }

Sequence* MultiSeq::find_sequence_by_id(const int id) {
    for (auto& s : sequences_) {
        if (s.id == id) return &s;
    }
    return nullptr;
}

Sequence* MultiSeq::find_sequence_by_name(const std::string& name) {
    for (auto& s : sequences_) {
        if (s.name == name) return &s;
    }
    return nullptr;
}

float MultiSeq::compute_seq_identity(size_t i, size_t j) const {
    if (i >= sequences_.size() || j >= sequences_.size()) {
        throw std::out_of_range("Index out of range in compute_seq_identity");
    }
    return sequences_[i].compute_seq_identity(sequences_[j]);
}

float MultiSeq::average_seq_identity() const {
    int n = sequences_.size();
    if (n < 2) return 1.0;
    float total = 0;
    int pairs = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            total += sequences_[i].compute_seq_identity(sequences_[j]);
            pairs++;
        }
    }
    return total / pairs;
}

bool MultiSeq::read_fasta(const std::string& filepath, const std::unordered_map<char, int>& encoding_scheme,
                          bool randomize_N) {
    std::ifstream infile(filepath);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open " << filepath << '\n';
        return false;
    }

    std::string line, name, seq;
    int k_id = 0;

    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!seq.empty()) {
                sequences_.emplace_back(seq, name, k_id++, encoding_scheme, randomize_N);
                seq.clear();
            }
            name = line.substr(1);
        } else {
            seq += line;
        }
    }

    if (!seq.empty()) {
        sequences_.emplace_back(seq, name, k_id++, encoding_scheme, randomize_N);
    }

    infile.close();
    return !sequences_.empty();
}

bool MultiSeq::write_fasta(const std::string& filepath, int max_line_length) const {
    std::ofstream outfile(filepath);
    if (!outfile.is_open()) {
        std::cerr << "Error: Cannot open file for writing: " << filepath << '\n';
        return false;
    }

    for (const auto& s : sequences_) {
        outfile << ">" << s.name << "\n";
        const std::string& str = s.seq;
        if (max_line_length <= 0) {
            outfile << str << "\n";
        } else {
            for (size_t i = 0; i < str.size(); i += max_line_length) {
                outfile << str.substr(i, max_line_length) << "\n";
            }
        }
    }

    outfile.close();
    return true;
}

void MultiSeq::print(std::ostream& out, bool use_enc) const {
    if (sequences_.empty()) {
        out << "No sequences available.\n";
        return;
    }
    for (const auto& s : sequences_) {
        s.print(out, use_enc);
    }
}
