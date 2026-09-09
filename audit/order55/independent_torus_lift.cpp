// Clean-room binary lift verifier.
//
// This intentionally uses a typed seven-coordinate key in an unordered_map
// with equality checks, rather than the production packed uint32 sorted
// table.  It also checks every shift, including shifts 11 and 22, after a
// join.  The input format is the committed lift_part*.txt format.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iterator>
#include <string>
#include <unordered_map>
#include <vector>

using u64 = std::uint64_t;

struct Signature {
    std::array<std::int16_t, 7> value{};
    bool operator==(const Signature& other) const { return value == other.value; }
};

struct SignatureHash {
    std::size_t operator()(const Signature& signature) const noexcept {
        std::size_t result = 1469598103934665603ULL;
        for (std::int16_t value : signature.value) {
            result ^= static_cast<std::uint16_t>(value);
            result *= 1099511628211ULL;
        }
        return result;
    }
};

using Bucket = std::unordered_map<Signature, std::vector<u64>, SignatureHash>;

int rows[5], columns[11], target[55], partial[7];
int split;
u64 nodes = 0, left_leaves = 0, right_leaves = 0, joined = 0, solutions = 0;
Bucket table;
std::string task;
std::ofstream word_output;

Signature signature(bool complement) {
    Signature result;
    for (int index = 0; index < 5; ++index) {
        int value = complement ? rows[index] - partial[index] : partial[index];
        result.value[static_cast<std::size_t>(index)] = static_cast<std::int16_t>(value);
    }
    for (int index = 5; index < 7; ++index) {
        int target_value = target[index == 5 ? 11 : 22];
        int value = complement ? target_value - partial[index] : partial[index];
        result.value[static_cast<std::size_t>(index)] = static_cast<std::int16_t>(value);
    }
    return result;
}

u64 place(int row, int column) {
    return static_cast<u64>((11 * row + 45 * column) % 55);
}

bool meets_target(u64 word) {
    const u64 mask = (u64{1} << 55) - 1;
    for (int shift = 1; shift <= 27; ++shift) {
        u64 rotated = ((word << shift) | (word >> (55 - shift))) & mask;
        if (std::popcount(word & rotated) != target[shift]) return false;
    }
    return true;
}

void visit_safe(int column, int end, u64 word, bool left) {
    ++nodes;
    if (column == end) {
        if (left) {
            ++left_leaves;
            table[signature(false)].push_back(word);
        } else {
            ++right_leaves;
            auto found = table.find(signature(true));
            if (found == table.end()) return;
            for (u64 left_word : found->second) {
                ++joined;
                u64 combined = word | left_word;
                if (!meets_target(combined)) continue;
                ++solutions;
                word_output << task << ' ';
                for (int position = 0; position < 55; ++position) {
                    word_output << (((combined >> position) & 1U) ? '1' : '0');
                }
                word_output << '\n';
            }
        }
        return;
    }
    for (unsigned mask = 0; mask < 32; ++mask) {
        if (std::popcount(mask) != columns[column]) continue;
        int near = std::popcount(mask & (((mask << 1) | (mask >> 4)) & 31U));
        int far = std::popcount(mask & (((mask << 2) | (mask >> 3)) & 31U));
        bool valid = partial[5] + near <= target[11]
                  && partial[6] + far <= target[22];
        for (int row = 0; row < 5; ++row) {
            int bit = static_cast<int>((mask >> row) & 1U);
            valid = valid && partial[row] + bit <= rows[row];
        }
        if (!valid) continue;
        for (int row = 0; row < 5; ++row) partial[row] += static_cast<int>((mask >> row) & 1U);
        partial[5] += near;
        partial[6] += far;
        u64 next = word;
        for (int row = 0; row < 5; ++row) {
            if ((mask >> row) & 1U) next |= u64{1} << place(row, column);
        }
        visit_safe(column + 1, end, next, left);
        for (int row = 0; row < 5; ++row) partial[row] -= static_cast<int>((mask >> row) & 1U);
        partial[5] -= near;
        partial[6] -= far;
    }
}

int main(int argc, char** argv) {
    if (argc != 5) return 2;
    std::ifstream input(argv[1]);
    std::ofstream audit(argv[2]);
    word_output.open(argv[3]);
    split = std::stoi(argv[4]);
    if (!input || !audit || !word_output || split < 1 || split > 10) return 3;
    while (input >> task) {
        for (int& value : rows) input >> value;
        for (int& value : columns) input >> value;
        for (int& value : target) input >> value;
        if (!input) return 4;
        table.clear();
        table.reserve(100000);
        std::fill(std::begin(partial), std::end(partial), 0);
        nodes = left_leaves = right_leaves = joined = solutions = 0;
        visit_safe(0, split, 0, true);
        // The map is intentionally unordered; equality is checked after
        // hashing and no sorted packed-key assumption is used.
        std::fill(std::begin(partial), std::end(partial), 0);
        visit_safe(split, 11, 0, false);
        audit << "{\"task\":\"" << task << "\",\"split\":" << split
              << ",\"nodes\":" << nodes << ",\"left_entries\":" << left_leaves
              << ",\"right_entries\":" << right_leaves << ",\"joined_words\":" << joined
              << ",\"solutions\":" << solutions << "}\n";
        audit.flush();
        word_output.flush();
    }
    return 0;
}
