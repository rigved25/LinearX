#include "./seq.hpp"

#include <fstream>
#include <iostream>

// default constructor
Seq::Seq() = default;

// parameterized constructor
Seq::Seq(const std::string &id, const std::string &sequence, const int k_id,
         std::unordered_map<char, int> *encoding_scheme)
    : id(id), sequence(sequence), k_id(k_id) {
    if (encoding_scheme) {
        this->set_encoding(*encoding_scheme);
    }
}

void Seq::set_encoding(std::unordered_map<char, int> &encoding_scheme) {
    enc_seq.clear();  // clear any existing encoding
    for (int i = 0; i < sequence.length(); ++i) {
        char nuc = sequence[i];
        if (encoding_scheme.find(nuc) != encoding_scheme.end()) {
            enc_seq.push_back(encoding_scheme.at(nuc));
        } else {
            std::cerr << "Warning: Nucleotide " << nuc << " not found in encoding scheme. Skipping." << std::endl;
        }
    }
}

// additional functionalities
size_t Seq::length() const { return sequence.length(); }

char Seq::at(size_t index) const {
    if (index >= sequence.length()) {
        std::cerr << "Warning: Invalid index " << index << " for at(). Returning first character." << std::endl;
        return sequence.at(0);
    }
    return sequence.at(index);
}

// overloaded [] operator for accessing characters
char &Seq::operator[](size_t index) {
    if (index >= sequence.length()) {
        std::cerr << "Warning: Invalid index " << index << " for operator[]. Returning the first character."
                  << std::endl;
        return sequence.at(0);
    }
    return sequence[index];
}

const char &Seq::operator[](size_t index) const {
    if (index >= sequence.length()) {
        std::cerr << "Warning: Invalid index " << index << " for operator[]. Returning the first character."
                  << std::endl;
        return sequence.at(0);
    }
    return sequence[index];
}

// remove gaps ('-') from the sequence
void Seq::remove_gaps() { sequence.erase(std::remove(sequence.begin(), sequence.end(), '-'), sequence.end()); }

// swap all occurences of nuc1 with nuc2
void Seq::swap_nuc(const char nuc1, const char nuc2) { std::replace(sequence.begin(), sequence.end(), nuc1, nuc2); }

float Seq::compute_seq_identity(const Seq &seq2) {
    if (this->length() != seq2.sequence.length()) {
        std::invalid_argument("Sequences must be of the same length to compute sequence identity");
    }
    int match_pos = 0;
    for (int i = 0; i < this->length(); ++i) {
        if (this->sequence[i] == seq2.sequence[i]) {
            match_pos++;
        }
    }
    return float(match_pos) / this->length();
}

// add a character to the sequence at the end
void Seq::add_nuc(char nucleotide) { sequence += nucleotide; }

// insert a character at a specific index
void Seq::insert_nuc(size_t index, char nucleotide) {
    if (index <= sequence.length()) {
        sequence.insert(sequence.begin() + index, nucleotide);
    } else {
        std::cerr << "Warning: Invalid index " << index << " for insert_char. No action taken." << std::endl;
    }
}

// delete a character at a specific index
void Seq::delete_nuc(size_t index) {
    if (index < sequence.length()) {
        sequence.erase(sequence.begin() + index);
    } else {
        std::cerr << "Warning: Invalid index " << index << " for delete_char. No action taken." << std::endl;
    }
}

// replace a character at a specific index
void Seq::replace_nuc(size_t index, char nucleotide) {
    if (index < sequence.length()) {
        sequence[index] = nucleotide;
    } else {
        std::cerr << "Warning: Invalid index " << index << " for replace_char. No action taken." << std::endl;
    }
}

// read a sequence from a FASTA file
bool Seq::read_fasta(const std::string &filepath) {
    std::ifstream infile(filepath);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open file " << filepath << std::endl;
        return false;
    }

    std::string line;

    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!id.empty()) {
                // if we already have an ID, we've found the end of the previous sequence
                break;
            }
            id = line.substr(1);  // remove the '>' character
        } else {
            sequence += line;
        }
    }

    infile.close();

    if (id.empty() || sequence.empty()) {
        std::cerr << "Warning: No valid sequence found in file " << filepath << std::endl;
        return false;
    }

    return true;
}

// write a sequence to a FASTA file
bool Seq::write_fasta(const std::string &filepath) const {
    std::ofstream outfile(filepath);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filepath << " for writing" << std::endl;
        return false;
    }

    outfile << ">" << id << "\n" << sequence << "\n";
    outfile.close();

    return true;
}

// print the sequence
void Seq::print() const {
    printf(">%s (k_id: %d)\n", id.c_str(), k_id);
    std::cout << sequence << std::endl;
}