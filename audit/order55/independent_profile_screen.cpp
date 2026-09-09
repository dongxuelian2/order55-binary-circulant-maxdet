// Clean-room exhaustive screen of the committed order-55 correlation
// partitions.  This deliberately uses next_permutation over each sorted
// multiset, rather than the production recursive counter, and recomputes the
// folded signatures, affine canonical test, and signed CRT product locally.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

constexpr u64 P = 2305843009213696591ULL;
constexpr u64 Q = 2305843009213697141ULL;
constexpr u128 MODULUS_PRODUCT = static_cast<u128>(P) * Q;
u128 maximum_product_threshold = 0;

u64 multiply(u64 left, u64 right, u64 modulus) {
    return static_cast<u64>(static_cast<u128>(left) * right % modulus);
}

u64 power(u64 base, u64 exponent, u64 modulus) {
    u64 result = 1;
    while (exponent != 0) {
        if ((exponent & 1U) != 0) result = multiply(result, base, modulus);
        base = multiply(base, base, modulus);
        exponent >>= 1U;
    }
    return result;
}

std::string decimal(u128 value) {
    std::string result;
    do {
        result.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    } while (value != 0);
    std::reverse(result.begin(), result.end());
    return result;
}

u128 parse_decimal(const std::string& text) {
    u128 value = 0;
    for (char digit : text) value = value * 10 + static_cast<unsigned>(digit - '0');
    return value;
}

template <std::size_t Size>
bool read_fold_file(const std::string& path, std::set<std::array<int, Size>>& output) {
    std::ifstream input(path);
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::istringstream parser(line);
        std::array<int, Size> signature{};
        for (int& value : signature) {
            if (!(parser >> value)) return false;
        }
        output.insert(signature);
    }
    return static_cast<bool>(input) || input.eof();
}

bool unit_canonical(const std::array<int, 27>& values) {
    for (int unit = 2; unit < 28; ++unit) {
        if (unit % 5 == 0 || unit % 11 == 0) continue;
        for (int shift = 1; shift <= 27; ++shift) {
            int transformed = (unit * shift) % 55;
            transformed = std::min(transformed, 55 - transformed);
            int candidate = values[static_cast<std::size_t>(transformed - 1)];
            int original = values[static_cast<std::size_t>(shift - 1)];
            if (candidate < original) return false;
            if (candidate > original) break;
        }
    }
    return true;
}

u128 absolute_profile_product(int k, const std::array<int, 27>& values,
                              const std::array<std::array<std::array<u64, 28>, 28>, 2>& coefficients) {
    const int base = k * (k - 1) / 54;
    std::array<u64, 2> residues{};
    const std::array<u64, 2> primes = {P, Q};
    for (std::size_t field = 0; field < primes.size(); ++field) {
        const u64 modulus = primes[field];
        u64 product = static_cast<u64>(55 - k);
        for (int frequency = 1; frequency <= 27; ++frequency) {
            u64 factor = static_cast<u64>(k - base) % modulus;
            for (int shift = 1; shift <= 27; ++shift) {
                const int deviation = values[static_cast<std::size_t>(shift - 1)] - base;
                const u64 coefficient = coefficients[field][static_cast<std::size_t>(frequency)]
                                                         [static_cast<std::size_t>(shift)];
                if (deviation >= 0) {
                    factor = (factor + multiply(static_cast<u64>(deviation), coefficient, modulus)) % modulus;
                } else {
                    factor = (factor + modulus
                              - multiply(static_cast<u64>(-deviation), coefficient, modulus)) % modulus;
                }
            }
            product = multiply(product, factor, modulus);
        }
        residues[field] = product;
    }
    const u64 inverse = power(P, Q - 2, Q);
    const u64 delta = (residues[1] + Q - residues[0]) % Q;
    const u128 positive = static_cast<u128>(residues[0])
                        + static_cast<u128>(P) * multiply(delta, inverse, Q);
    return positive > MODULUS_PRODUCT / 2 ? MODULUS_PRODUCT - positive : positive;
}

struct Counters {
    u64 ordered_profiles = 0;
    u64 folded_matching_profiles = 0;
    u64 canonical_determinants = 0;
    u64 above_screen = 0;
    u128 maximum_product = 0;
};

int main(int argc, char** argv) {
    if (argc != 7) return 2;
    const int k = std::stoi(argv[1]);
    maximum_product_threshold = parse_decimal(argv[2]);
    const std::string partition_path = argv[3];
    const std::string certificate_directory = argv[4];
    std::ofstream audit(argv[5]);
    std::ofstream targets(argv[6]);
    if (!audit || !targets) return 3;

    std::set<std::array<int, 3>> allowed5;
    std::set<std::array<int, 6>> allowed11;
    if (!read_fold_file(certificate_directory + "/foldcorr_m5_k" + std::to_string(k) + ".txt", allowed5)
        || !read_fold_file(certificate_directory + "/foldcorr_m11_k" + std::to_string(k) + ".txt", allowed11)
        || allowed5.empty() || allowed11.empty()) return 4;

    std::array<std::array<std::array<u64, 28>, 28>, 2> coefficients{};
    const std::array<u64, 2> primes = {P, Q};
    for (std::size_t field = 0; field < primes.size(); ++field) {
        const u64 modulus = primes[field];
        u64 root = 0;
        for (u64 generator = 2; generator < 1000 && root == 0; ++generator) {
            const u64 candidate = power(generator, (modulus - 1) / 55, modulus);
            if (power(candidate, 55, modulus) == 1
                && power(candidate, 5, modulus) != 1
                && power(candidate, 11, modulus) != 1) root = candidate;
        }
        if (root == 0) return 5;
        for (int frequency = 1; frequency <= 27; ++frequency) {
            for (int shift = 1; shift <= 27; ++shift) {
                const int exponent = (frequency * shift) % 55;
                const int opposite = (55 - exponent) % 55;
                coefficients[field][static_cast<std::size_t>(frequency)]
                                     [static_cast<std::size_t>(shift)] =
                    (power(root, static_cast<u64>(exponent), modulus)
                     + power(root, static_cast<u64>(opposite), modulus)) % modulus;
            }
        }
    }

    Counters counters;
    std::ifstream partitions(partition_path);
    std::string line;
    while (std::getline(partitions, line)) {
        if (line.empty()) continue;
        std::istringstream parser(line);
        std::array<int, 27> values{};
        for (int& value : values) {
            if (!(parser >> value)) return 6;
        }
        std::sort(values.begin(), values.end());
        do {
            ++counters.ordered_profiles;
            std::array<int, 3> fold5{};
            std::array<int, 6> fold11{};
            fold5[0] = k;
            fold11[0] = k;
            for (int shift = 1; shift <= 27; ++shift) {
                const int residue5 = shift % 5;
                const int residue11 = shift % 11;
                fold5[static_cast<std::size_t>(std::min(residue5, 5 - residue5))]
                    += values[static_cast<std::size_t>(shift - 1)] * (residue5 == 0 ? 2 : 1);
                fold11[static_cast<std::size_t>(std::min(residue11, 11 - residue11))]
                    += values[static_cast<std::size_t>(shift - 1)] * (residue11 == 0 ? 2 : 1);
            }
            if (!allowed5.count(fold5) || !allowed11.count(fold11)) continue;
            ++counters.folded_matching_profiles;
            if (!unit_canonical(values)) continue;
            ++counters.canonical_determinants;
            const u128 product = absolute_profile_product(k, values, coefficients);
            if (product > counters.maximum_product) counters.maximum_product = product;
            if (product < maximum_product_threshold) continue;
            ++counters.above_screen;
            targets << "{\"k\":" << k << ",\"absolute_profile_product\":\""
                    << decimal(product) << "\",\"correlations\":[";
            for (std::size_t index = 0; index < values.size(); ++index) {
                if (index != 0) targets << ',';
                targets << values[index];
            }
            targets << "]}\n";
        } while (std::next_permutation(values.begin(), values.end()));
    }
    if (!partitions && !partitions.eof()) return 7;
    audit << "{\"k\":" << k
          << ",\"ordered_profiles\":" << counters.ordered_profiles
          << ",\"folded_matching_profiles\":" << counters.folded_matching_profiles
          << ",\"canonical_determinants\":" << counters.canonical_determinants
          << ",\"above_screen\":" << counters.above_screen
          << ",\"max_absolute_profile_product\":\""
          << decimal(counters.maximum_product) << "\"}\n";
    return 0;
}
