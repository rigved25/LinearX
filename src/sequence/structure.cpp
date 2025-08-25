// src/sequence/structure.cpp
#include <iostream>
#include <linearx/sequence/structure.hpp>
#include <linearx/utility/utils.hpp>
#include <set>
#include <stack>
#include <unordered_map>

using namespace std;

// Constructor
Structure::Structure(int len) : length(len) {
    if (len <= 0) {
        throw invalid_argument("Length must be positive.");
    }
    basePairs.resize(len, -1);  // Initialize with -1 (unpaired)
}

string Structure::removeShortHelices(const string& dotBracket, int minHelixSize) {
    // Map to match opening and closing brackets
    const unordered_map<char, char> bracketPair = {{')', '('}, {']', '['}, {'}', '{'}, {'>', '<'}};

    // Reverse map for opening brackets
    const unordered_map<char, char> reversePair = {{'(', ')'}, {'[', ']'}, {'{', '}'}, {'<', '>'}};

    // Vector to store base pairs as (open, close) indices
    vector<pair<int, int>> pairs;

    // Stacks for each bracket type
    unordered_map<char, stack<int>> stacks;
    for (const auto& entry : reversePair) {
        stacks[entry.first] = stack<int>();
    }

    // Step 1: Identify base pairs
    for (int i = 0; i < dotBracket.size(); i++) {
        char ch = dotBracket[i];
        if (reversePair.count(ch)) {  // Opening bracket
            stacks[ch].push(i);
        } else if (bracketPair.count(ch)) {  // Closing bracket
            char openBracket = bracketPair.at(ch);
            if (!stacks[openBracket].empty()) {
                int openIndex = stacks[openBracket].top();
                stacks[openBracket].pop();
                pairs.emplace_back(openIndex, i);
            }
        }
    }

    // Step 2: Group pairs into helices
    sort(pairs.begin(), pairs.end());  // Ensure pairs are sorted by open indices
    vector<vector<pair<int, int>>> helices;
    vector<pair<int, int>> currentHelix;

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
    string result = dotBracket;
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

// Get the length of the structure
int Structure::getLength() const { return length; }

// Get the paired index for a base
int Structure::getPair(int index) const {
    if (index < 0 || index >= length) {
        throw out_of_range("Index out of range.");
    }
    return basePairs[index];
}

// Set a pair (bi-directionally)
void Structure::addPair(int i, int j) {
    if (i < 0 || i >= length || j < 0 || j >= length) {
        throw out_of_range("Indices out of range.");
    }
    if (basePairs[i] != -1 || basePairs[j] != -1) {
        throw invalid_argument("One or both indices are already paired.");
    }
    basePairs[i] = j;
    basePairs[j] = i;
}

// Remove a pair (bi-directionally)
void Structure::removePair(int i) {
    if (i < 0 || i >= length) {
        throw out_of_range("Index out of range.");
    }
    int j = basePairs[i];
    if (j != -1) {  // If paired
        basePairs[i] = -1;
        basePairs[j] = -1;
    }
}

// Get the dot bracket representation of the structure
string Structure::getDotBracket() const {
    const vector<pair<int, int>> symbols = {{'(', ')'}, {'[', ']'}, {'{', '}'}, {'<', '>'}};
    set<pair<int, int>> remaining;
    for (int i = 0; i < length; ++i) {
        int j = basePairs[i];
        if (j > i) remaining.insert({i, j});
    }

    vector<set<pair<int, int>>> selections;
    for (int iter = 0; iter < 4 && !remaining.empty(); ++iter) {
        unordered_map<int, int> pairing;
        for (auto& [i, j] : remaining) {
            pairing[i] = j;
            pairing[j] = i;
        }

        vector<vector<int>> dp(length, vector<int>(length, 0));
        unordered_map<pair<int, int>, int, linearx::utils::PairHash> back;

        for (int span = 2; span <= length; ++span) {
            for (int i = 0; i + span - 1 < length; ++i) {
                int j = i + span - 1;
                dp[i][j] = dp[i][j - 1];
                back[{i, j}] = j;  // j unpaired
                auto it = pairing.find(j);
                if (it != pairing.end()) {
                    int t = it->second;
                    if (t >= i && t < j) {
                        int score = 1;
                        if (t > i) score += dp[i][t - 1];
                        if (t + 1 <= j - 1) score += dp[t + 1][j - 1];
                        if (score > dp[i][j]) {
                            dp[i][j] = score;
                            back[{i, j}] = t;
                        }
                    }
                }
            }
        }

        set<pair<int, int>> selected;
        function<void(int, int)> traceback = [&](int i, int j) {
            if (i > j || back.find({i, j}) == back.end()) return;
            int t = back[{i, j}];
            if (t == j) {
                traceback(i, j - 1);
            } else {
                selected.insert({t, j});
                traceback(i, t - 1);
                traceback(t + 1, j - 1);
            }
        };

        traceback(0, length - 1);
        selections.push_back(selected);
        for (auto& p : selected) remaining.erase(p);
    }

    if (!remaining.empty()) {
        cerr << "[WARNING] Could not incorporate all pairs into dot-bracket: " << remaining.size() << " left.\n";
    }

    string db(length, '.');
    for (int k = 0; k < selections.size(); ++k) {
        for (auto& [i, j] : selections[k]) {
            db[i] = symbols[k].first;
            db[j] = symbols[k].second;
        }
    }
    return db;
}

// Print all base pairs
void Structure::printBasePairs() const {
    cout << "Base pairs:\n";
    for (int i = 0; i < length; ++i) {
        int j = basePairs[i];
        if (j > i) {  // Print each pair only once
            cout << i << " - " << j << "\n";
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
