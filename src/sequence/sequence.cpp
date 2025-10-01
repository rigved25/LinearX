// linearx/sequence/seq.cpp
#include <algorithm>
#include <fstream>
#include <linearx/sequence/sequence.hpp>
#include <stdexcept>

// default constructor
Sequence::Sequence() : name("seq"), id(0), seq("") {}

// constructor with only sequence
Sequence::Sequence(const std::string& seq_str) : name("seq"), id(0), seq(seq_str) {
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
}

// constructor with sequence, name, id
Sequence::Sequence(const std::string& seq_str, const std::string& name, int id) : name(name), id(id), seq(seq_str) {
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
}

// constructor with sequence, name, id, encoding
Sequence::Sequence(const std::string& seq_str, const std::string& name, int id,
                   const std::unordered_map<char, int>& encoding_scheme)
    : name(name), id(id), seq(seq_str) {
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    if (!encoding_scheme.empty()) {
        set_encoding(encoding_scheme);
    }
}

void Sequence::set_encoding(const std::unordered_map<char, int>& encoding_scheme) {
    enc.clear();
    for (char ch : seq) {
        if (encoding_scheme.count(ch)) {
            enc.push_back(encoding_scheme.at(ch));
        } else {
            throw std::invalid_argument(std::string("Nucleotide '") + ch + "' not in encoding scheme");
        }
    }
}

size_t Sequence::length() const { return seq.length(); }

char Sequence::at(size_t index) const {
    if (index >= seq.length()) throw std::out_of_range("Index out of range");
    return seq.at(index);
}

char& Sequence::operator[](size_t index) {
    if (index >= seq.length()) throw std::out_of_range("Index out of range");
    return seq[index];
}

const char& Sequence::operator[](size_t index) const {
    if (index >= seq.length()) throw std::out_of_range("Index out of range");
    return seq[index];
}

bool Sequence::equals(const Sequence& other) const { return seq == other.seq && id == other.id && name == other.name; }

void Sequence::remove_gaps() { seq.erase(std::remove(seq.begin(), seq.end(), '-'), seq.end()); }

void Sequence::swap_nuc(char nuc1, char nuc2) { std::replace(seq.begin(), seq.end(), nuc1, nuc2); }

float Sequence::compute_seq_identity(const Sequence& seq2) const {
    if (this->length() != seq2.length()) {
        throw std::invalid_argument("Sequences must be of equal length");
    }
    int valid_positions = 0, matches = 0;
    for (size_t i = 0; i < seq.size(); ++i) {
        if (seq[i] == '-' && seq2.seq[i] == '-') {
            continue;  // skip double-gap positions
        }
        ++valid_positions;
        if (seq[i] == seq2.seq[i]) {
            ++matches;
        }
    }
    if (valid_positions == 0) return 0.0f;  // avoid division by zero
    return static_cast<float>(matches) / valid_positions;
}

void Sequence::set_seq(const std::string& new_seq) {
    seq = new_seq;
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
}

void Sequence::add_nuc(char nucleotide) { seq += nucleotide; }

void Sequence::insert_nuc(size_t index, char nucleotide) {
    if (index > seq.length()) throw std::out_of_range("Insert index out of bounds");
    seq.insert(seq.begin() + index, nucleotide);
}

void Sequence::delete_nuc(size_t index) {
    if (index >= seq.length()) throw std::out_of_range("Delete index out of bounds");
    seq.erase(seq.begin() + index);
}

void Sequence::replace_nuc(size_t index, char nucleotide) {
    if (index >= seq.length()) throw std::out_of_range("Replace index out of bounds");
    seq[index] = nucleotide;
}

void Sequence::reverse() { std::reverse(seq.begin(), seq.end()); }

bool Sequence::read_fasta(const std::string& filepath) {
    std::ifstream infile(filepath);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << filepath << std::endl;
        return false;
    }

    std::string line;
    bool reading_sequence = false;

    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (reading_sequence) break;  // already read one sequence
            name = line.substr(1);        // remove '>'
            reading_sequence = true;
        } else if (reading_sequence) {
            for (char c : line) {
                seq += std::toupper(c);
            }
        }
    }

    infile.close();

    if (name.empty() || seq.empty()) {
        std::cerr << "Warning: No valid sequence found in file " << filepath << std::endl;
        return false;
    }

    return true;
}

bool Sequence::read_fasta_stream(std::istream& is) {
    std::string line, seq_data, seq_name;
    while (std::getline(is, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            seq_name = line.substr(1);  // everything after '>'
        } else {
            for (char c : line) {
                if (!std::isspace(static_cast<unsigned char>(c))) {
                    seq_data.push_back(std::toupper(static_cast<unsigned char>(c)));
                }
            }
        }
    }
    if (seq_data.empty()) {
        return false;  // no sequence parsed
    }
    this->seq = seq_data;

    return true;
}

bool Sequence::write_fasta(const std::string& filepath) const {
    std::ofstream outfile(filepath);
    if (!outfile.is_open()) return false;

    outfile << ">" << name << "\n" << seq << "\n";
    return true;
}

void Sequence::show_details(std::ostream& os) const {
    os << "[Sequence Info] Name: " << name << " | ID: " << id << " | Length: " << seq.length() << "\n";
}

void Sequence::print(std::ostream& os, bool use_enc) const {
    os << ">" << name << " (id: " << id << ")\n";
    if (use_enc && !enc.empty()) {
        for (int c : enc) {
            os << c;
        }
    } else {
        os << seq;
    }
    os << "\n";
}

bool Sequence::write_fasta(const std::string& filepath, size_t max_line_length) const {
    std::ofstream outfile(filepath);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filepath << " for writing" << std::endl;
        return false;
    }

    outfile << ">" << name << "\n";

    if (max_line_length <= 0) {
        outfile << seq << "\n";
    } else {
        for (size_t i = 0; i < seq.length(); i += max_line_length) {
            outfile << seq.substr(i, max_line_length) << "\n";
        }
    }

    outfile.close();
    return true;
}
