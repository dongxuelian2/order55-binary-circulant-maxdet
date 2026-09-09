// Brute-force small-order regression for the complement/fold/profile pipeline.
// This is a clean-room implementation and does not import or execute order-55 code.
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
constexpr u64 PRIME = 2305843009213693951ULL; // 2^61-1; PRIME == 1 (mod 15,21)

u64 mul(u64 a, u64 b) { return (u128)a * b % PRIME; }
u64 power(u64 a, u64 e) {
    u64 result = 1;
    while (e) {
        if (e & 1) result = mul(result, a);
        a = mul(a, a);
        e >>= 1;
    }
    return result;
}
u64 root_of_order(int n) {
    for (u64 candidate = 2; ; ++candidate) {
        u64 root = power(candidate, (PRIME - 1) / n);
        bool exact = power(root, n) == 1;
        for (int q : {3, 5, 7}) if (n % q == 0) exact = exact && power(root, n / q) != 1;
        if (exact) return root;
    }
}
u64 rotate_right(u64 word, int shift, int n) {
    const u64 mask = (u64(1) << n) - 1;
    return ((word >> shift) | (word << (n - shift))) & mask;
}
u64 encode(const std::vector<int>& values, int base) {
    u64 code = 0;
    for (int value : values) code = code * base + (u64)value;
    return code;
}
std::vector<int> correlations(u64 word, int n) {
    std::vector<int> result(n);
    for (int shift = 0; shift < n; ++shift)
        result[shift] = std::popcount(word & rotate_right(word, shift, n));
    return result;
}
std::vector<int> fold(u64 word, int n, int modulus) {
    std::vector<int> result(modulus);
    for (int i = 0; i < n; ++i) if ((word >> i) & 1) ++result[i % modulus];
    return result;
}
std::string affine_canonical(u64 word, int n) {
    std::string best(n, '2');
    for (int u = 1; u < n; ++u) {
        if (std::gcd(u, n) != 1) continue;
        for (int t = 0; t < n; ++t) {
            std::string image;
            image.reserve(n);
            for (int i = 0; i < n; ++i) image.push_back(((word >> ((u * i + t) % n)) & 1) ? '1' : '0');
            if (image < best) best = std::move(image);
        }
    }
    return best;
}
u64 hash_strings(const std::set<std::string>& values) {
    u64 hash = 1469598103934665603ULL;
    for (const auto& value : values) {
        for (unsigned char c : value) { hash ^= c; hash *= 1099511628211ULL; }
        hash ^= 0xff; hash *= 1099511628211ULL;
    }
    return hash;
}
struct Key {
    u64 profile, first_fold, second_fold;
    bool operator<(const Key& other) const {
        if (profile != other.profile) return profile < other.profile;
        if (first_fold != other.first_fold) return first_fold < other.first_fold;
        return second_fold < other.second_fold;
    }
    bool operator==(const Key& other) const {
        return profile == other.profile && first_fold == other.first_fold && second_fold == other.second_fold;
    }
};
struct Record { Key key; std::uint32_t word; };

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    const int n = std::stoi(argv[1]);
    const int other_factor = std::stoi(argv[2]);
    const int half = n / 2;
    const u64 domain = u64(1) << n;
    const u64 mask = domain - 1;
    const u64 root = root_of_order(n);
    const long double hadamard = std::pow((long double)n, (long double)n / 2.0L);
    std::array<std::array<u64, 21>, 21> powers{};
    for (int j = 0; j < n; ++j) for (int t = 0; t < n; ++t)
        powers[j][t] = power(root, (u64)j * t);

    std::vector<std::int64_t> determinants(domain);
    std::vector<u64> eigen(n, 0);
    std::vector<std::uint32_t> brute_max_words;
    std::int64_t maximum = -1;
    std::array<std::int64_t, 22> max_by_weight{};
    max_by_weight.fill(-1);
    u64 gray = 0;
    for (u64 index = 0; index < domain; ++index) {
        if (index) {
            const int bit = std::countr_zero(index);
            gray ^= u64(1) << bit;
            const bool now_set = (gray >> bit) & 1;
            for (int j = 0; j < n; ++j) {
                if (now_set) {
                    eigen[j] += powers[j][bit];
                    if (eigen[j] >= PRIME) eigen[j] -= PRIME;
                } else {
                    if (eigen[j] < powers[j][bit]) eigen[j] += PRIME;
                    eigen[j] -= powers[j][bit];
                }
            }
        }
        u64 residue = 1;
        for (int j = 0; j < n; ++j) residue = mul(residue, eigen[j]);
        std::int64_t signed_value = residue <= PRIME / 2 ? (std::int64_t)residue
                                                          : (std::int64_t)(residue - PRIME);
        std::int64_t determinant = signed_value < 0 ? -signed_value : signed_value;
        determinants[gray] = determinant;
        const int weight = std::popcount(gray);
        max_by_weight[weight] = std::max(max_by_weight[weight], determinant);
        if (determinant > maximum) { maximum = determinant; brute_max_words.clear(); }
        if (determinant == maximum) brute_max_words.push_back((std::uint32_t)gray);
    }

    std::vector<Record> records;
    records.reserve(domain / 2);
    std::set<u64> first_folds, second_folds;
    const int base = n + 1;
    for (u64 word = 0; word < domain; ++word) {
        if (std::popcount(word) > half) continue;
        const auto corr = correlations(word, n);
        const auto f3 = fold(word, n, 3);
        const auto fother = fold(word, n, other_factor);
        const Key key{encode(std::vector<int>(corr.begin() + 1, corr.begin() + half + 1), base),
                      encode(f3, base), encode(fother, base)};
        records.push_back(Record{key, (std::uint32_t)word});
        first_folds.insert(key.first_fold); second_folds.insert(key.second_fold);
    }
    std::sort(records.begin(), records.end(), [](const Record& a, const Record& b) {
        if (a.key < b.key) return true;
        if (b.key < a.key) return false;
        return a.word < b.word;
    });

    std::vector<std::uint32_t> reconstructed;
    reconstructed.reserve(records.size() * 2);
    std::vector<std::uint32_t> reduced_max_words;
    u64 profile_groups = 0;
    bool fold_membership = true;
    for (std::size_t begin = 0; begin < records.size();) {
        std::size_t end = begin + 1;
        while (end < records.size() && records[end].key == records[begin].key) ++end;
        ++profile_groups;
        for (std::size_t i = begin; i < end; ++i) {
            const u64 word = records[i].word;
            fold_membership = fold_membership && first_folds.count(records[i].key.first_fold)
                              && second_folds.count(records[i].key.second_fold);
            reconstructed.push_back((std::uint32_t)word);
            reconstructed.push_back((std::uint32_t)(mask ^ word));
            if (determinants[word] == maximum) reduced_max_words.push_back((std::uint32_t)word);
            if (determinants[mask ^ word] == maximum) reduced_max_words.push_back((std::uint32_t)(mask ^ word));
        }
        begin = end;
    }
    std::sort(reconstructed.begin(), reconstructed.end());
    std::vector<std::uint32_t> expected;
    expected.reserve(records.size() * 2);
    for (u64 word = 0; word < domain; ++word) expected.push_back((std::uint32_t)word);
    const bool coverage = reconstructed == std::vector<std::uint32_t>(expected.begin(), expected.end());

    std::sort(reduced_max_words.begin(), reduced_max_words.end());
    std::sort(brute_max_words.begin(), brute_max_words.end());
    const bool max_words_equal = reduced_max_words == brute_max_words;
    std::set<std::string> brute_classes, reduced_classes;
    for (u64 word : brute_max_words) brute_classes.insert(affine_canonical(word, n));
    for (u64 word : reduced_max_words) reduced_classes.insert(affine_canonical(word, n));
    const bool classes_equal = brute_classes == reduced_classes;
    const bool modulus_ok = (long double)PRIME > 2.0L * hadamard;
    const bool pass = modulus_ok && fold_membership && coverage && max_words_equal && classes_equal;

    std::cout << "ORDER " << n
              << " modulus=" << PRIME
              << " root=" << root
              << " modulus_exceeds_twice_hadamard=" << (modulus_ok ? 1 : 0)
              << " brute_words=" << domain
              << " normalized_words=" << records.size()
              << " profile_groups=" << profile_groups
              << " fold3_signatures=" << first_folds.size()
              << " fold" << other_factor << "_signatures=" << second_folds.size()
              << " coverage=" << (coverage ? 1 : 0)
              << " fold_membership=" << (fold_membership ? 1 : 0)
              << " maximum=" << maximum
              << " brute_maximizers=" << brute_max_words.size()
              << " reduced_maximizers=" << reduced_max_words.size()
              << " max_words_equal=" << (max_words_equal ? 1 : 0)
              << " brute_affine_classes=" << brute_classes.size()
              << " reduced_affine_classes=" << reduced_classes.size()
              << " affine_classes_equal=" << (classes_equal ? 1 : 0)
              << " brute_affine_hash=" << hash_strings(brute_classes)
              << " reduced_affine_hash=" << hash_strings(reduced_classes)
              << " status=" << (pass ? "PASS" : "FAIL") << '\n';
    return pass ? 0 : 1;
}
