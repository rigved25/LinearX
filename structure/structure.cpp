#include "structure.hpp"

#include <stdexcept>

std::string Structure::removeShortHelices(const std::string& dotBracket, int minHelixSize) {
    // Map to match opening and closing brackets
    std::unordered_map<char, char> bracketPair = {{')', '('}, {']', '['}, {'}', '{'}, {'>', '<'}};

    // Reverse map for opening brackets
    std::unordered_map<char, char> reversePair = {{'(', ')'}, {'[', ']'}, {'{', '}'}, {'<', '>'}};

    // Vector to store base pairs as (open, close) indices
    std::vector<std::pair<int, int>> pairs;

    // Stacks for each bracket type
    std::unordered_map<char, std::stack<int>> stacks;
    for (const auto& entry : reversePair) {
        stacks[entry.first] = std::stack<int>();
    }

    // Step 1: Identify base pairs
    for (int i = 0; i < dotBracket.size(); i++) {
        char ch = dotBracket[i];
        if (reversePair.count(ch)) {  // Opening bracket
            stacks[ch].push(i);
        } else if (bracketPair.count(ch)) {  // Closing bracket
            char openBracket = bracketPair[ch];
            if (!stacks[openBracket].empty()) {
                int openIndex = stacks[openBracket].top();
                stacks[openBracket].pop();
                pairs.emplace_back(openIndex, i);
            }
        }
    }

    // Step 2: Group pairs into helices
    sort(pairs.begin(), pairs.end());  // Ensure pairs are sorted by open indices
    std::vector<std::vector<std::pair<int, int>>> helices;
    std::vector<std::pair<int, int>> currentHelix;

    for (const auto& pair : pairs) {
        if (!currentHelix.empty() && pair.first == currentHelix.back().first + 1 &&
            pair.second == currentHelix.back().second - 1) {
            currentHelix.push_back(pair);
        } else {
            if (!currentHelix.empty()) {
                helices.push_back(currentHelix);
            }
            currentHelix = {pair};
        }
    }
    if (!currentHelix.empty()) {
        helices.push_back(currentHelix);
    }

    // Step 3: Remove short helices
    std::string result = dotBracket;
    for (const auto& helix : helices) {
        if (helix.size() < minHelixSize) {
            for (const auto& [open, close] : helix) {
                result[open] = '.';
                result[close] = '.';
            }
        }
    }

    return result;
}

// Constructor
Structure::Structure(int len) : length(len) {
    if (len <= 0) {
        throw std::invalid_argument("Length must be positive.");
    }
    basePairs.resize(len, -1);  // Initialize with -1 (unpaired)
}

// Get the length of the structure
int Structure::getLength() const { return length; }

// Get the paired index for a base
int Structure::getPair(int index) const {
    if (index < 0 || index >= length) {
        throw std::out_of_range("Index out of range.");
    }
    return basePairs[index];
}

// Set a pair (bi-directionally)
void Structure::addPair(int i, int j) {
    if (i < 0 || i >= length || j < 0 || j >= length) {
        throw std::out_of_range("Indices out of range.");
    }
    basePairs[i] = j;
    basePairs[j] = i;
}

// Remove a pair (bi-directionally)
void Structure::removePair(int i) {
    if (i < 0 || i >= length) {
        throw std::out_of_range("Index out of range.");
    }
    int j = basePairs[i];
    if (j != -1) {  // If paired
        basePairs[i] = -1;
        basePairs[j] = -1;
    }
}

// Get the dot bracket representation of the structure
std::string& Structure::getDotBracket() const {
    // Vector to store base pairs as (open, close) indices
    std::vector<std::pair<int, int>> pairs;

    // Copy the base pairs
    for (int i = 0; i < length; ++i) {
        int j = basePairs[i];
        if (j > i) {
            pairs.emplace_back(i, j);
        }
    }

    std::vector<std::pair<int, int>> pseudo_pairs1;
    std::vector<std::pair<int, int>> pseudo_pairs2;
    std::vector<std::pair<int, int>> pseudo_pairs3;

    // check for pseudoknots 1
    for (int i = 0; i < pairs.size(); i++) {
        for (int j = i + 1; j < pairs.size(); j++) {
            if (pairs[i].first < pairs[j].first && pairs[j].first < pairs[i].second &&
                pairs[i].second < pairs[j].second) {
                pseudo_pairs1.push_back(std::make_pair(pairs[j].first, pairs[j].second));
            }
        }
    }

    // check for pseudoknots 2
    for (int i = 0; i < pseudo_pairs1.size(); i++) {
        for (int j = i + 1; j < pseudo_pairs1.size(); j++) {
            if (pseudo_pairs1[i].first < pseudo_pairs1[j].first && pseudo_pairs1[j].first < pseudo_pairs1[i].second &&
                pseudo_pairs1[i].second < pseudo_pairs1[j].second) {
                pseudo_pairs2.push_back(std::make_pair(pseudo_pairs1[j].first, pseudo_pairs1[j].second));
            }
        }
    }

    // check for pseudoknots 3
    for (int i = 0; i < pseudo_pairs2.size(); i++) {
        for (int j = i + 1; j < pseudo_pairs2.size(); j++) {
            if (pseudo_pairs2[i].first < pseudo_pairs2[j].first && pseudo_pairs2[j].first < pseudo_pairs2[i].second &&
                pseudo_pairs2[i].second < pseudo_pairs2[j].second) {
                pseudo_pairs3.push_back(std::make_pair(pseudo_pairs2[j].first, pseudo_pairs2[j].second));
            }
        }
    }

    std::string& structure = *new std::string(length, '.');

    // normal pairs
    for (int i = 0; i < pairs.size(); i++) {
        structure[pairs[i].first] = '(';
        structure[pairs[i].second] = ')';
    }

    // pseudoknots 1
    for (int i = 0; i < pseudo_pairs1.size(); i++) {
        structure[pseudo_pairs1[i].first] = '[';
        structure[pseudo_pairs1[i].second] = ']';
    }

    // pseudoknots 2
    for (int i = 0; i < pseudo_pairs2.size(); i++) {
        structure[pseudo_pairs2[i].first] = '{';
        structure[pseudo_pairs2[i].second] = '}';
    }

    // pseudoknots 3
    for (int i = 0; i < pseudo_pairs3.size(); i++) {
        structure[pseudo_pairs3[i].first] = '<';
        structure[pseudo_pairs3[i].second] = '>';
    }

    return structure;
}

// Print all base pairs
void Structure::printBasePairs() const {
    std::cout << "Base pairs:\n";
    for (int i = 0; i < length; ++i) {
        int j = basePairs[i];
        if (j > i) {  // Print each pair only once
            std::cout << i << " - " << j << "\n";
        }
    }
}

void Structure::removeShortHelices(int minHelixSize) {
    int pairs, i, j;

    for (i = 0; i < length; i++) {
        if (getPair(i) > i) {
            j = getPair(i);
            pairs = 1;

            // Count the helix length
            while ((i + 1 < length && getPair(i + 1) == j - 1) || (i + 2 < length && getPair(i + 2) == j - 1) ||
                   (i + 1 < length && getPair(i + 1) == j - 2)) {
                if (i + 1 < length && getPair(i + 1) == j - 1) {
                    i++;
                    j--;
                    pairs++;
                } else if (i + 2 < length && getPair(i + 2) == j - 1) {
                    if (getPair(i + 1) != -1) {
                        removePair(getPair(i + 1));
                        removePair(i + 1);
                    }
                    i += 2;
                    j--;
                    pairs++;
                } else {
                    i++;
                    j -= 2;
                    pairs++;
                }
            }

            // Remove helix if it's too short
            if (pairs < minHelixSize) {
                removePair(i);

                // Remove pairs backward
                while (i > 0 && ((getPair(i - 1) == j + 1) || (i - 2 >= 0 && getPair(i - 2) == j + 1) ||
                                 (getPair(i - 1) == j + 2))) {
                    if (getPair(i - 1) == j + 1) {
                        removePair(getPair(i - 1));
                        removePair(i - 1);
                        i--;
                        j++;
                    } else if (i - 2 >= 0 && getPair(i - 2) == j + 1) {
                        removePair(getPair(i - 2));
                        removePair(i - 2);
                        i -= 2;
                        j++;
                    } else {
                        removePair(getPair(i - 1));
                        removePair(i - 1);
                        i--;
                        j += 2;
                    }
                }
            }
        }
    }
}
