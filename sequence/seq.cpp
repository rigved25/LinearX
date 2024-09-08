#include "./seq.hpp"
#include <fstream>
#include <iostream>


// default constructor
Seq::Seq() = default;

// parameterized constructor
Seq::Seq(const std::string& id, const std::string& sequence)
    : id_(id), sequence_(sequence) {}

// accessors
const std::string& Seq::get_id() const { return id_; }

const std::string& Seq::get_sequence() const { return sequence_; }

// mutators
void Seq::set_id(const std::string& id) { id_ = id; }

void Seq::set_sequence(const std::string& sequence) { sequence_ = sequence; }

// additional functionalities
size_t Seq::length() const { return sequence_.length(); }

char Seq::at(size_t index) const {
    if (index >= sequence_.length()) {
        std::cerr << "Warning: Invalid index " << index << " for at(). Returning first character." << std::endl;
        return sequence_.at(0);
    }
    return sequence_.at(index);
}

// overloaded [] operator for accessing characters
char& Seq::operator[](size_t index) {
    if (index >= sequence_.length()) {
        std::cerr << "Warning: Invalid index " << index << " for operator[]. Returning the first character." << std::endl;
        return sequence_.at(0);
    }
    return sequence_[index];
}

const char& Seq::operator[](size_t index) const {
    if (index >= sequence_.length()) {
        std::cerr << "Warning: Invalid index " << index << " for operator[]. Returning the first character." << std::endl;
        return sequence_.at(0);
    }
    return sequence_[index];
}

// remove gaps ('-') from the sequence
void Seq::remove_gaps() {
    sequence_.erase(std::remove(sequence_.begin(), sequence_.end(), '-'), sequence_.end());
}

// swap all occurences of nuc1 with nuc2
void Seq::swap_nuc(const char nuc1, const char nuc2) {
    std::replace(sequence_.begin(), sequence_.end(), nuc1, nuc2);
}

// add a character to the sequence at the end
void Seq::add_nuc(char nucleotide) {
    sequence_ += nucleotide;
}

// insert a character at a specific index
void Seq::insert_nuc(size_t index, char nucleotide) {
    if (index <= sequence_.length()) {
        sequence_.insert(sequence_.begin() + index, nucleotide);
    } else {
        std::cerr << "Warning: Invalid index " << index << " for insert_char. No action taken." << std::endl;
    }
}

// delete a character at a specific index
void Seq::delete_nuc(size_t index) {
    if (index < sequence_.length()) {
        sequence_.erase(sequence_.begin() + index);
    } else {
        std::cerr << "Warning: Invalid index " << index << " for delete_char. No action taken." << std::endl;
    }
}

// replace a character at a specific index
void Seq::replace_nuc(size_t index, char nucleotide) {
    if (index < sequence_.length()) {
        sequence_[index] = nucleotide;
    } else {
        std::cerr << "Warning: Invalid index " << index << " for replace_char. No action taken." << std::endl;
    }
}

// read a sequence from a FASTA file
bool Seq::read_fasta(const std::string& filepath) {
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
                // if we already have an ID, we've found the end of the previous sequence
                break;
            }
            id = line.substr(1); // remove the '>' character
        } else {
            sequence += line;
        }
    }

    infile.close();

    if (id.empty() || sequence.empty()) {
        std::cerr << "Warning: No valid sequence found in file " << filepath << std::endl;
        return false;
    }

    set_id(id);
    set_sequence(sequence);
    return true;
}

// write a sequence to a FASTA file
bool Seq::write_fasta(const std::string& filepath) const {
    std::ofstream outfile(filepath);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file " << filepath << " for writing" << std::endl;
        return false;
    }

    outfile << ">" << id_ << "\n" << sequence_ << "\n";
    outfile.close();

    return true;
}

// static methods
// encode the sequence given an encoding scheme in a dictionary
std::vector<int> Seq::encode(const Seq &seq, std::unordered_map<char, int> scheme) {
    std::vector<int> encoded_seq;
    for(int i = 0; i < seq.length(); ++i) {
        char nuc = seq[i];
        if (scheme.find(nuc) != scheme.end()) {
            encoded_seq.push_back(scheme[nuc]);
        } else {
            std::cerr << "Warning: Nucleotide " << nuc << " not found in encoding scheme. Skipping." << std::endl;
        }
    }
    return encoded_seq;
}
