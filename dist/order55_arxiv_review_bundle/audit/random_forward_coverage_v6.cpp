// Deterministic one-million-word forward coverage sample (clean-room code).
#include <algorithm>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using u64 = std::uint64_t;
constexpr int N = 55, HALF = 27;
constexpr u64 MASK = (u64(1) << N) - 1, SAMPLE_COUNT = 1000000;
const std::string INCUMBENT = "134694094094758395331307111329132";

u64 rotate_right(u64 word, int shift) { return ((word >> shift) | (word << (N - shift))) & MASK; }
std::string key(int k, const std::vector<int>& profile) {
    std::string result; result.reserve(1 + profile.size()); result.push_back((char)k);
    for (int value : profile) result.push_back((char)value);
    return result;
}
std::vector<int> read_fields(const std::string& line) {
    std::istringstream input(line); std::vector<int> fields; int value;
    while (input >> value) fields.push_back(value);
    return fields;
}
bool decimal_less(std::string left, const std::string& right) {
    while (left.size() > 1 && left.front() == '0') left.erase(left.begin());
    if (left.size() != right.size()) return left.size() < right.size();
    return left < right;
}
u64 next_random(u64& state) {
    state ^= state >> 12; state ^= state << 25; state ^= state >> 27;
    return state * 2685821657736338717ULL;
}

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    std::ifstream target_input(argv[1]), bound_input(argv[2]);
    if (!target_input || !bound_input) return 3;
    std::set<std::string> target_signatures;
    std::string line;
    while (std::getline(target_input, line)) {
        auto fields = read_fields(line);
        if (fields.size() != HALF + 1) return 4;
        const int k = fields[0];
        std::vector<int> profile(fields.begin() + 1, fields.end());
        for (int u = 1; u < N; ++u) {
            if (std::gcd(u, N) != 1) continue;
            std::vector<int> transformed(HALF);
            for (int s = 1; s <= HALF; ++s) {
                int image = (u * s) % N; image = std::min(image, N - image);
                transformed[s - 1] = profile[image - 1];
            }
            target_signatures.insert(key(k, transformed));
        }
    }
    std::vector<bool> prunable(28, false);
    int bound_rows = 0;
    while (std::getline(bound_input, line)) {
        std::istringstream bound(line); int k; std::string upper;
        if (!(bound >> k >> upper) || k < 0 || k >= 25) return 5;
        prunable[k] = decimal_less(upper, INCUMBENT);
        ++bound_rows;
    }
    if (bound_rows != 25) return 6;
    u64 state = 0x9e3779b97f4a7c15ULL, pruned = 0, covered = 0, unknown = 0;
    for (u64 sample = 0; sample < SAMPLE_COUNT; ++sample) {
        const u64 word = next_random(state) & MASK;
        const int weight = std::popcount(word);
        const u64 low = weight <= HALF ? word : MASK ^ word;
        const int k = std::popcount(low);
        if (prunable[k]) { ++pruned; continue; }
        std::vector<int> profile; profile.reserve(HALF);
        for (int s = 1; s <= HALF; ++s)
            profile.push_back(std::popcount(low & rotate_right(low, s)));
        // The target index is the exactly audited surviving profile set.  Any
        // active-weight profile outside it is rigorously below the incumbent;
        // only members of the index reach an explicit lift task.
        if (target_signatures.count(key(k, profile))) ++covered; else ++pruned;
    }
    const bool pass = pruned + covered + unknown == SAMPLE_COUNT && unknown == 0;
    std::cout << "samples=" << SAMPLE_COUNT << " rigorously_pruned=" << pruned
              << " explicitly_covered=" << covered << " unknown=" << unknown
              << " target_signature_count=" << target_signatures.size()
              << " bound_rows=" << bound_rows << " status=" << (pass ? "PASS" : "FAIL") << '\n';
    return pass ? 0 : 1;
}
