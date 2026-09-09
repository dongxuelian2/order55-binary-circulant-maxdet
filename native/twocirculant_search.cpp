#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {
using Clock = std::chrono::steady_clock;
using Complex = std::complex<long double>;

struct Config {
  int n = 0, type = 1, threads = 1;
  std::uint64_t seed = 1, iterations = 0, cycle = 20000;
  double seconds = 10.0;
  std::string a, b;
};

struct Candidate {
  long double score = -std::numeric_limits<long double>::infinity();
  std::vector<signed char> a, b;
  std::uint64_t iteration = 0;
  int worker = -1;
};

std::uint64_t parse_u64(const char* value, const char* option) {
  try {
    std::size_t used = 0;
    const auto parsed = std::stoull(value, &used);
    if (used != std::string(value).size()) throw std::invalid_argument("tail");
    return parsed;
  } catch (...) {
    throw std::invalid_argument(std::string("invalid value for ") + option);
  }
}

int parse_int(const char* value, const char* option) {
  const auto parsed = parse_u64(value, option);
  if (parsed > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
    throw std::invalid_argument(std::string("value too large for ") + option);
  return static_cast<int>(parsed);
}

Config parse_args(int argc, char** argv) {
  Config config;
  for (int index = 1; index < argc; ++index) {
    const std::string option = argv[index];
    auto next = [&]() -> const char* {
      if (++index >= argc) throw std::invalid_argument("missing value for " + option);
      return argv[index];
    };
    if (option == "--n") config.n = parse_int(next(), "--n");
    else if (option == "--type") config.type = parse_int(next(), "--type");
    else if (option == "--threads") config.threads = parse_int(next(), "--threads");
    else if (option == "--seed") config.seed = parse_u64(next(), "--seed");
    else if (option == "--iterations") { config.iterations = parse_u64(next(), "--iterations"); config.seconds = 0; }
    else if (option == "--seconds") { config.seconds = std::stod(next()); config.iterations = 0; }
    else if (option == "--cycle") config.cycle = parse_u64(next(), "--cycle");
    else if (option == "--a") config.a = next();
    else if (option == "--b") config.b = next();
    else if (option == "--help") {
      std::cout << "twocirculant_search --n N --type {1|3} --a WORD --b WORD "
                   "[--threads T] [--seconds S|--iterations I] [--seed S] [--cycle I]\n";
      std::exit(0);
    } else throw std::invalid_argument("unknown option: " + option);
  }
  if (config.n < 3 || config.n % 2 == 0) throw std::invalid_argument("n must be odd and at least 3");
  const int m = (config.n - 1) / 2;
  if (config.type != 1 && config.type != 3) throw std::invalid_argument("type must be 1 or 3");
  if (static_cast<int>(config.a.size()) != m || static_cast<int>(config.b.size()) != m)
    throw std::invalid_argument("a and b must both have length (n-1)/2");
  for (char value : config.a + config.b)
    if (value != '+' && value != '-') throw std::invalid_argument("words use only + and -");
  if (config.threads < 1 || config.threads > 256) throw std::invalid_argument("threads must be 1..256");
  if (config.cycle < 10) throw std::invalid_argument("cycle must be at least 10");
  if (config.iterations == 0 && config.seconds <= 0) throw std::invalid_argument("seconds must be positive");
  return config;
}

std::vector<signed char> parse_word(const std::string& word) {
  std::vector<signed char> result;
  result.reserve(word.size());
  for (char value : word) result.push_back(value == '+' ? 1 : -1);
  return result;
}

std::string format_word(const std::vector<signed char>& word) {
  std::string result;
  result.reserve(word.size());
  for (auto value : word) result.push_back(value == 1 ? '+' : '-');
  return result;
}

class Objective {
 public:
  struct State {
    std::vector<signed char> a, b;
    std::vector<Complex> fa, fb;
    long double score;
  };

  Objective(int m, int type) : m_(m), type_(type), roots_(m * ((m + 1) / 2)) {
    const long double tau = 2 * std::acos(-1.0L);
    for (int frequency = 0; frequency < (m_ + 1) / 2; ++frequency)
      for (int position = 0; position < m_; ++position) {
        const long double angle = -tau * frequency * position / m_;
        roots_[frequency * m_ + position] = {std::cos(angle), std::sin(angle)};
      }
  }

  State make_state(std::vector<signed char> a, std::vector<signed char> b) const {
    State state{std::move(a), std::move(b), std::vector<Complex>((m_ + 1) / 2),
                std::vector<Complex>((m_ + 1) / 2), 0};
    for (int frequency = 0; frequency < (m_ + 1) / 2; ++frequency)
      for (int position = 0; position < m_; ++position) {
        state.fa[frequency] += static_cast<long double>(state.a[position]) * root(frequency, position);
        state.fb[frequency] += static_cast<long double>(state.b[position]) * root(frequency, position);
      }
    state.score = score(state.fa, state.fb);
    return state;
  }

  long double delta(const State& state, int which, int position) const {
    const signed char old = which == 0 ? state.a[position] : state.b[position];
    return changed_score(state, which, position, -2.0L * old) - state.score;
  }

  void flip(State& state, int which, int position) const {
    auto& word = which == 0 ? state.a : state.b;
    auto& spectrum = which == 0 ? state.fa : state.fb;
    const long double change = -2.0L * word[position];
    word[position] = -word[position];
    for (int frequency = 0; frequency < (m_ + 1) / 2; ++frequency)
      spectrum[frequency] += change * root(frequency, position);
    state.score = score(state.fa, state.fb);
  }

 private:
  const Complex& root(int frequency, int position) const { return roots_[frequency * m_ + position]; }

  long double zero_log(long double alpha, long double beta) const {
    const long double determinant = type_ == 1
        ? -(alpha * alpha + beta * beta + 2.0L * m_ * beta)
        : alpha * alpha + beta * beta + 2.0L * m_ * alpha;
    return std::abs(determinant) < 1e-30L
        ? -std::numeric_limits<long double>::infinity() : std::log(std::abs(determinant));
  }

  long double score(const std::vector<Complex>& fa, const std::vector<Complex>& fb) const {
    long double result = zero_log(fa[0].real(), fb[0].real());
    if (!std::isfinite(result)) return result;
    for (int frequency = 1; frequency < (m_ + 1) / 2; ++frequency) {
      const long double q = std::norm(fa[frequency]) + std::norm(fb[frequency]);
      if (q < 1e-30L) return -std::numeric_limits<long double>::infinity();
      result += 2 * std::log(q);
    }
    return result;
  }

  long double changed_score(const State& state, int which, int position, long double change) const {
    long double result = zero_log(state.fa[0].real() + (which == 0 ? change : 0),
                                  state.fb[0].real() + (which == 1 ? change : 0));
    if (!std::isfinite(result)) return result;
    for (int frequency = 1; frequency < (m_ + 1) / 2; ++frequency) {
      const Complex changed_a = state.fa[frequency] +
          (which == 0 ? change * root(frequency, position) : Complex(0, 0));
      const Complex changed_b = state.fb[frequency] +
          (which == 1 ? change * root(frequency, position) : Complex(0, 0));
      const long double q = std::norm(changed_a) + std::norm(changed_b);
      if (q < 1e-30L) return -std::numeric_limits<long double>::infinity();
      result += 2 * std::log(q);
    }
    return result;
  }

  int m_, type_;
  std::vector<Complex> roots_;
};

class SharedBest {
 public:
  void consider(const Objective::State& state, std::uint64_t iteration, int worker) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (state.score > best_.score + 1e-13L)
      best_ = Candidate{state.score, state.a, state.b, iteration, worker};
  }
  Candidate get() const { std::lock_guard<std::mutex> lock(mutex_); return best_; }
 private:
  mutable std::mutex mutex_;
  Candidate best_;
};

void quench(const Objective& objective, Objective::State& state,
            std::atomic<std::uint64_t>& evaluations) {
  for (;;) {
    long double best_delta = 1e-14L;
    int best_which = -1, best_position = -1;
    for (int which = 0; which < 2; ++which) {
      const int size = static_cast<int>(which == 0 ? state.a.size() : state.b.size());
      for (int position = 0; position < size; ++position) {
        const long double change = objective.delta(state, which, position);
        evaluations.fetch_add(1, std::memory_order_relaxed);
        if (change > best_delta) {
          best_delta = change;
          best_which = which;
          best_position = position;
        }
      }
    }
    if (best_which < 0) return;
    objective.flip(state, best_which, best_position);
  }
}

void run_worker(const Config& config, const Objective& objective, SharedBest& shared,
                std::atomic<std::uint64_t>& evaluations, Clock::time_point deadline,
                int worker_id) {
  std::mt19937_64 rng(config.seed + 0x9e3779b97f4a7c15ULL * (worker_id + 1));
  auto state = objective.make_state(parse_word(config.a), parse_word(config.b));
  std::uniform_int_distribution<int> which_distribution(0, 1);
  std::uniform_int_distribution<int> position_distribution(0, (config.n - 3) / 2);
  std::uniform_real_distribution<long double> unit(0, 1);
  for (int kick = 0; kick < (worker_id == 0 ? 0 : 2 + worker_id % 11); ++kick)
    objective.flip(state, which_distribution(rng), position_distribution(rng));
  quench(objective, state, evaluations);
  shared.consider(state, 0, worker_id);

  std::uint64_t iteration = 0;
  while (true) {
    if (config.iterations ? iteration >= config.iterations
                          : ((iteration & 1023ULL) == 0 && Clock::now() >= deadline)) break;
    const int which = which_distribution(rng), position = position_distribution(rng);
    const long double change = objective.delta(state, which, position);
    evaluations.fetch_add(1, std::memory_order_relaxed);
    const long double phase = static_cast<long double>(iteration % config.cycle) / config.cycle;
    const long double temperature = 0.35L * std::pow(0.003L / 0.35L, phase);
    if (change >= 0 || unit(rng) < std::exp(change / temperature))
      objective.flip(state, which, position);
    ++iteration;
    if (iteration % config.cycle == 0) {
      quench(objective, state, evaluations);
      shared.consider(state, iteration, worker_id);
      const int kicks = 2 + static_cast<int>((iteration / config.cycle + worker_id) % 13);
      for (int kick = 0; kick < kicks; ++kick)
        objective.flip(state, which_distribution(rng), position_distribution(rng));
    }
  }
  quench(objective, state, evaluations);
  shared.consider(state, iteration, worker_id);
}
}  // namespace

int main(int argc, char** argv) {
  try {
    const Config config = parse_args(argc, argv);
    Objective objective((config.n - 1) / 2, config.type);
    SharedBest shared;
    std::atomic<std::uint64_t> evaluations{0};
    const auto started = Clock::now();
    const auto deadline = started +
        std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(config.seconds));
    std::vector<std::thread> threads;
    for (int worker = 0; worker < config.threads; ++worker)
      threads.emplace_back(run_worker, std::cref(config), std::cref(objective),
                           std::ref(shared), std::ref(evaluations), deadline, worker);
    for (auto& thread : threads) thread.join();
    const double elapsed = std::chrono::duration<double>(Clock::now() - started).count();
    const Candidate best = shared.get();
    std::cout << "{\n  \"n\": " << config.n << ",\n  \"type\": " << config.type
              << ",\n  \"a\": \"" << format_word(best.a) << "\",\n  \"b\": \""
              << format_word(best.b) << "\",\n  \"logabsdet\": " << std::setprecision(20)
              << best.score << ",\n  \"seed\": " << config.seed << ",\n  \"threads\": "
              << config.threads << ",\n  \"evaluations\": " << evaluations.load()
              << ",\n  \"elapsed_seconds\": " << std::setprecision(10) << elapsed
              << ",\n  \"evaluations_per_second\": " << evaluations.load() / elapsed
              << ",\n  \"worker\": " << best.worker << ",\n  \"iteration\": "
              << best.iteration << "\n}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << "\n";
    return 2;
  }
}
