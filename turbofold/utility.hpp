#ifndef TF_UTILITIES_HPP
#define TF_UTILITIES_HPP

#include <unordered_map>
#include <tuple>

#include "./../utility/constants.hpp"


enum VerboseState {
    SILENT = 0,
    VERBOSE = 1,
    DEBUG = 2,
    DETAIL = 3
};

// used for hashing a triple of integers
struct TripleHash {
    std::size_t operator()(const std::tuple<int, int, int>& t) const {
        // Extract t1, t2, t3
        std::size_t t1 = static_cast<std::size_t>(std::get<0>(t));
        std::size_t t2 = static_cast<std::size_t>(std::get<1>(t));
        std::size_t t3 = static_cast<std::size_t>(std::get<2>(t));

        // Combine t2 and t3 using bit manipulation (as you already have)
        std::size_t hash = t2 ^ ((t3 << 16) | (t3 >> 16));

        // Mix in t1 with the current hash
        hash ^= t1 + GOLDEN_RATIO + (hash << 6) + (hash >> 2);

        // Apply the additional mix to further reduce collisions
        hash ^= (hash >> 13);
        hash *= GOLDEN_RATIO;
        hash ^= (hash >> 15);

        return hash;
    }
};


class Cache3D {
public:
    // Insert a value into the cache
    void insert(int k, int i, int j, double value) {
        cache[std::make_tuple(k, i, j)] = value;
    }

    // Check if a value exists in the cache
    bool exists(int k, int i, int j) {
        return cache.find(std::make_tuple(k, i, j)) != cache.end();
    }

    // Retrieve a value from the cache
    double get(int k, int i, int j) {
        return cache[std::make_tuple(k, i, j)];
    }

    // Clear the entire cache
    void clear() {
        cache.clear();
    }

private:
    // The cache data structure using a tuple (k, i, j) as the key
    std::unordered_map<std::tuple<int, int, int>, double, TripleHash> cache;
};

#endif // TF_UTILITIES_HPP
