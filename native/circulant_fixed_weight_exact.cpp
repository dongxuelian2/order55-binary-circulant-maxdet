#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {
using u64 = std::uint64_t;
using u128 = unsigned __int128;
using Clock = std::chrono::steady_clock;

constexpr int order = 55;
constexpr u64 prime1 = 2305843009213696591ULL;
constexpr u64 prime2 = 2305843009213697141ULL;

u64 multiply_mod(u64 a, u64 b, u64 modulus) {
  return static_cast<u64>((static_cast<u128>(a) * b) % modulus);
}

u64 power_mod(u64 base, u64 exponent, u64 modulus) {
  u64 result = 1;
  while (exponent) {
    if (exponent & 1) result = multiply_mod(result, base, modulus);
    base = multiply_mod(base, base, modulus);
    exponent >>= 1;
  }
  return result;
}

u64 primitive_55th_root(u64 prime) {
  for (u64 base = 2;; ++base) {
    const u64 root = power_mod(base, (prime - 1) / order, prime);
    if (root != 1 && power_mod(root, 5, prime) != 1 &&
        power_mod(root, 11, prime) != 1 && power_mod(root, 55, prime) == 1)
      return root;
  }
}

struct ModularSpectrum {
  u64 prime;
  u64 roots[order][order]{};

  explicit ModularSpectrum(u64 modulus) : prime(modulus) {
    const u64 omega = primitive_55th_root(prime);
    for (int frequency = 0; frequency < order; ++frequency)
      for (int position = 0; position < order; ++position)
        roots[frequency][position] = power_mod(omega, static_cast<u64>(frequency * position), prime);
  }

  u64 determinant(const std::vector<int>& positions) const {
    u64 result = 1;
    for (int frequency = 0; frequency < order; ++frequency) {
      u64 eigenvalue = 0;
      for (int position : positions) {
        eigenvalue += roots[frequency][position];
        if (eigenvalue >= prime) eigenvalue -= prime;
      }
      result = multiply_mod(result, eigenvalue, prime);
    }
    return result;
  }
};

std::string decimal(u128 value) {
  if (value == 0) return "0";
  std::string result;
  while (value) {
    result.push_back(static_cast<char>('0' + value % 10));
    value /= 10;
  }
  std::reverse(result.begin(), result.end());
  return result;
}

u128 absolute_crt(u64 residue1, u64 residue2) {
  const u64 inverse = power_mod(prime1 % prime2, prime2 - 2, prime2);
  const u64 difference = residue2 >= residue1 ? residue2 - residue1 : prime2 - (residue1 - residue2);
  const u64 multiplier = multiply_mod(difference, inverse, prime2);
  const u128 modulus = static_cast<u128>(prime1) * prime2;
  const u128 value = residue1 + static_cast<u128>(prime1) * multiplier;
  return value <= modulus / 2 ? value : modulus - value;
}

struct Best {
  u128 determinant = 0;
  u64 residue1 = 0, residue2 = 0;
  std::vector<int> positions;
};

struct Config { int weight = 0, threads = 1; };

Config parse_args(int argc, char** argv) {
  Config config;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    if (++index >= argc) throw std::invalid_argument("missing value for " + option);
    if (option == "--weight") config.weight = std::stoi(argv[index]);
    else if (option == "--threads") config.threads = std::stoi(argv[index]);
    else throw std::invalid_argument("unknown option: " + option);
  }
  if (config.weight < 1 || config.weight > 7) throw std::invalid_argument("weight must be 1..7");
  if (config.threads < 1 || config.threads > 256) throw std::invalid_argument("threads must be 1..256");
  return config;
}

void enumerate_tail(int next, int depth, int weight, std::vector<int>& positions,
                    const ModularSpectrum& spectrum1, const ModularSpectrum& spectrum2,
                    Best& local_best, u64& local_count) {
  if (depth == weight) {
    const u64 residue1 = spectrum1.determinant(positions);
    const u64 residue2 = spectrum2.determinant(positions);
    const u128 determinant = absolute_crt(residue1, residue2);
    ++local_count;
    if (determinant > local_best.determinant)
      local_best = Best{determinant, residue1, residue2, positions};
    return;
  }
  const int remaining = weight - depth;
  for (int position = next; position <= order - remaining; ++position) {
    positions[depth] = position;
    enumerate_tail(position + 1, depth + 1, weight, positions, spectrum1, spectrum2,
                   local_best, local_count);
  }
}

std::string word_from_positions(const std::vector<int>& positions) {
  std::string word(order, '0');
  for (int position : positions) word[position] = '1';
  return word;
}
}  // namespace

int main(int argc, char** argv) {
  try {
    const Config config = parse_args(argc, argv);
    const auto started = Clock::now();
    const ModularSpectrum spectrum1(prime1), spectrum2(prime2);
    std::atomic<int> next_first{1};
    std::mutex best_mutex;
    Best global_best;
    std::atomic<u64> candidates{0};
    std::vector<std::thread> workers;
    for (int worker = 0; worker < config.threads; ++worker) {
      workers.emplace_back([&, worker] {
        Best local_best;
        u64 local_count = 0;
        std::vector<int> positions(config.weight);
        positions[0] = 0;
        if (config.weight == 1 && worker == 0) {
          const u64 residue1 = spectrum1.determinant(positions);
          const u64 residue2 = spectrum2.determinant(positions);
          local_best = Best{absolute_crt(residue1, residue2), residue1, residue2, positions};
          local_count = 1;
        } else if (config.weight > 1) {
          for (;;) {
            const int first = next_first.fetch_add(1);
            if (first > order - (config.weight - 1)) break;
            positions[1] = first;
            enumerate_tail(first + 1, 2, config.weight, positions, spectrum1, spectrum2,
                           local_best, local_count);
          }
        }
        candidates.fetch_add(local_count);
        std::lock_guard<std::mutex> lock(best_mutex);
        if (local_best.determinant > global_best.determinant) global_best = local_best;
      });
    }
    for (auto& worker : workers) worker.join();
    const double elapsed = std::chrono::duration<double>(Clock::now() - started).count();
    const long double bound_log2 = 0.5L * order * std::log2(static_cast<long double>(config.weight));
    std::cout << "{\n  \"problem\": \"binary_circulant_fixed_weight\",\n"
              << "  \"order\": " << order << ",\n  \"weight\": " << config.weight
              << ",\n  \"absolute_determinant\": \"" << decimal(global_best.determinant)
              << "\",\n  \"word\": \"" << word_from_positions(global_best.positions)
              << "\",\n  \"residue_prime_1\": " << global_best.residue1
              << ",\n  \"residue_prime_2\": " << global_best.residue2
              << ",\n  \"prime_1\": " << prime1 << ",\n  \"prime_2\": " << prime2
              << ",\n  \"candidate_subsets_with_zero\": " << candidates.load()
              << ",\n  \"coverage\": \"every rotation orbit has a representative containing position 0\",\n"
              << "  \"hadamard_bound_log2\": " << std::setprecision(10)
              << static_cast<double>(bound_log2) << ",\n  \"crt_capacity_bits_lower_bound\": "
              << 121 << ",\n  \"threads\": " << config.threads
              << ",\n  \"elapsed_seconds\": " << elapsed << ",\n  \"candidates_per_second\": "
              << candidates.load() / elapsed << "\n}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << "\n";
    return 2;
  }
}
