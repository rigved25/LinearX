#ifndef STRUCTURE_HPP
#define STRUCTURE_HPP

#include <iostream>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>

class Structure {
   private:
    int length;                  // Length of the sequence
    std::vector<int> basePairs;  // Base pair mappings

   public:
    // Constructor
    Structure(int len);

    // Get the length of the structure
    int getLength() const;

    // Get the paired index for a base
    int getPair(int index) const;

    // Set a pair (bi-directionally)
    void addPair(int i, int j);

    // Remove a pair (bi-directionally)
    void removePair(int i);

    // get the dot bracket representation of the structure
    std::string& getDotBracket() const;

    // Print all base pairs
    void printBasePairs() const;

    // Remove short helices
    void removeShortHelices(int minHelixLength);

    // Static methods below
    static std::string removeShortHelices(const std::string& dotBracket, int minHelixLength);
};

#endif  // STRUCTURE_HPP
