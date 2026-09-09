// Clean-room exact brute-force ground truth for small odd composite orders.
// It enumerates every binary word in Gray-code order.  Fourier eigenvalues
// are updated one bit at a time in two independent prime fields, and the
// signed CRT product is exact under the checked Hadamard bound.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;
using i128 = __int128_t;

constexpr std::array<u64, 2> PRIMES = {1000000891ULL, 1000001311ULL};

u64 multiply(u64 left, u64 right, u64 modulus) {
    return static_cast<u64>(static_cast<u128>(left) * right % modulus);
}

u64 power(u64 base, u64 exponent, u64 modulus) {
    u64 result = 1;
    while (exponent != 0) {
        if (exponent & 1U) result = multiply(result, base, modulus);
        base = multiply(base, base, modulus);
        exponent >>= 1U;
    }
    return result;
}

u64 primitive_root(int order, u64 modulus) {
    std::array<int, 2> factors = {3, order == 15 ? 5 : 7};
    for (u64 generator = 2; generator < 1000; ++generator) {
        u64 root = power(generator, (modulus - 1) / static_cast<u64>(order), modulus);
        bool primitive = power(root, static_cast<u64>(order), modulus) == 1;
        for (int factor : factors) {
            primitive = primitive && power(root, static_cast<u64>(order / factor), modulus) != 1;
        }
        if (primitive) return root;
    }
    throw std::runtime_error("primitive root search failed");
}

std::string word_string(u64 word, int order) {
    std::string result(static_cast<std::size_t>(order), '0');
    for (int position = 0; position < order; ++position) {
        if ((word >> position) & 1U) result[static_cast<std::size_t>(position)] = '1';
    }
    return result;
}

u64 parse_u128_to_u64(i128 value) {
    if (value < 0 || static_cast<u128>(value) > static_cast<u128>(UINT64_MAX)) {
        throw std::runtime_error("small-order determinant outside u64");
    }
    return static_cast<u64>(value);
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const int order = std::stoi(argv[1]);
    if (order != 15 && order != 21) return 3;
    const u64 total = u64{1} << order;
    const u128 modulus_product = static_cast<u128>(PRIMES[0]) * PRIMES[1];
    const std::array<int, 2> prime_factors = {3, order == 15 ? 5 : 7};

    std::array<std::vector<std::vector<u64>>, 2> powers;
    std::array<std::vector<u64>, 2> eigenvalues;
    for (std::size_t field = 0; field < PRIMES.size(); ++field) {
        const u64 modulus = PRIMES[field];
        const u64 root = primitive_root(order, modulus);
        powers[field].assign(static_cast<std::size_t>(order), std::vector<u64>(static_cast<std::size_t>(order)));
        eigenvalues[field].assign(static_cast<std::size_t>(order), 0);
        for (int frequency = 0; frequency < order; ++frequency) {
            for (int position = 0; position < order; ++position) {
                powers[field][static_cast<std::size_t>(frequency)][static_cast<std::size_t>(position)] =
                    power(root, static_cast<u64>(frequency) * static_cast<u64>(position), modulus);
            }
        }
        if (power(root, static_cast<u64>(order), modulus) != 1
            || power(root, static_cast<u64>(order / prime_factors[0]), modulus) == 1
            || power(root, static_cast<u64>(order / prime_factors[1]), modulus) == 1) return 4;
    }

    u64 inverse = power(PRIMES[0], PRIMES[1] - 2, PRIMES[1]);
    auto reconstruct = [&](u64 first, u64 second) -> i128 {
        u64 delta = (second + PRIMES[1] - first) % PRIMES[1];
        u64 coefficient = multiply(delta, inverse, PRIMES[1]);
        u128 positive = static_cast<u128>(first) + static_cast<u128>(PRIMES[0]) * coefficient;
        u128 half = modulus_product / 2;
        if (positive > half) return -static_cast<i128>(modulus_product - positive);
        return static_cast<i128>(positive);
    };

    u128 hadamard_square = 1;
    for (int i = 0; i < order; ++i) hadamard_square *= static_cast<u64>(order);
    if (hadamard_square >= modulus_product * static_cast<u128>(modulus_product) / 4) return 5;

    i128 best = -1;
    std::vector<u64> maximizing_words;
    std::vector<i128> weight_best(static_cast<std::size_t>(order + 1), -1);
    std::vector<u64> weight_count(static_cast<std::size_t>(order + 1), 0);
    u64 previous_gray = 0;
    for (u64 index = 0; index < total; ++index) {
        u64 gray = index ^ (index >> 1U);
        if (index != 0) {
            u64 changed = previous_gray ^ gray;
            int position = std::countr_zero(changed);
            bool turned_on = (gray >> position) & 1U;
            for (std::size_t field = 0; field < PRIMES.size(); ++field) {
                const u64 modulus = PRIMES[field];
                for (int frequency = 0; frequency < order; ++frequency) {
                    u64 term = powers[field][static_cast<std::size_t>(frequency)][static_cast<std::size_t>(position)];
                    if (turned_on) {
                        eigenvalues[field][static_cast<std::size_t>(frequency)] =
                            (eigenvalues[field][static_cast<std::size_t>(frequency)] + term) % modulus;
                    } else {
                        u64 current = eigenvalues[field][static_cast<std::size_t>(frequency)];
                        eigenvalues[field][static_cast<std::size_t>(frequency)] =
                            current >= term ? current - term : modulus - (term - current);
                    }
                }
            }
        }
        previous_gray = gray;

        std::array<u64, 2> residues = {1, 1};
        for (std::size_t field = 0; field < PRIMES.size(); ++field) {
            const u64 modulus = PRIMES[field];
            u64 value = 1;
            for (u64 eigenvalue : eigenvalues[field]) value = multiply(value, eigenvalue, modulus);
            residues[field] = value;
        }
        i128 determinant = reconstruct(residues[0], residues[1]);
        if (determinant < 0) return 6;
        int weight = std::popcount(gray);
        if (determinant > best) {
            best = determinant;
            maximizing_words.clear();
        }
        if (determinant == best) maximizing_words.push_back(gray);
        if (determinant > weight_best[static_cast<std::size_t>(weight)]) {
            weight_best[static_cast<std::size_t>(weight)] = determinant;
            weight_count[static_cast<std::size_t>(weight)] = 1;
        } else if (determinant == weight_best[static_cast<std::size_t>(weight)]) {
            ++weight_count[static_cast<std::size_t>(weight)];
        }
    }

    std::cout << "{\"order\":" << order << ",\"words\":" << total
              << ",\"maximum\":" << parse_u128_to_u64(best)
              << ",\"maximizers\":[";
    for (std::size_t i = 0; i < maximizing_words.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << '"' << word_string(maximizing_words[i], order) << '"';
    }
    std::cout << "],\"weight_best\":[";
    for (int weight = 0; weight <= order; ++weight) {
        if (weight) std::cout << ',';
        std::cout << "{\"weight\":" << weight
                  << ",\"maximum\":" << parse_u128_to_u64(weight_best[static_cast<std::size_t>(weight)])
                  << ",\"count\":" << weight_count[static_cast<std::size_t>(weight)] << '}';
    }
    std::cout << "]}\n";
    return 0;
}
