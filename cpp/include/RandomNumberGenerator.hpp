#pragma once

#include "Vector3D.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <vector>

/**
 * @brief Thread-safe random number generator for Monte Carlo radiation
 * transport simulation
 *
 * Modern C++17 implementation providing comprehensive random sampling
 * capabilities for Monte Carlo particle transport with high-quality random
 * number generation and thread safety for parallel simulations.
 */
class RandomNumberGenerator
{
public:
  using SeedType = std::uint64_t;
  using GeneratorType = std::mt19937_64;

private:
  mutable GeneratorType generator_;
  mutable std::mutex mutex_;
  SeedType current_seed_;
  static std::atomic<SeedType> global_seed_counter_;

  // Distribution caching for performance
  mutable std::uniform_real_distribution<double> uniform_dist_;
  mutable std::normal_distribution<double> normal_dist_;
  mutable std::exponential_distribution<double> exponential_dist_;

public:
  // Constructors
  RandomNumberGenerator();
  explicit RandomNumberGenerator(SeedType seed);

  // Copy constructor and assignment (creates independent generator with
  // different seed)
  RandomNumberGenerator(const RandomNumberGenerator &other);
  RandomNumberGenerator &operator=(const RandomNumberGenerator &other);

  // Move constructor and assignment
  RandomNumberGenerator(RandomNumberGenerator &&other) noexcept;
  RandomNumberGenerator &operator=(RandomNumberGenerator &&other) noexcept;

  // Seed management
  void setSeed(SeedType seed);
  SeedType getSeed() const { return current_seed_; }
  void randomSeed(); // Set random seed based on current time and thread
  static SeedType generateUniqueSeed();

  // Basic uniform sampling [0, 1)
  double uniform() const;
  double uniform(double min, double max) const;

  // Integer uniform sampling [min, max]
  int uniformInt(int min, int max) const;
  std::size_t uniformSize(std::size_t min, std::size_t max) const;

  // C++17 safe sampling with std::optional
  std::optional<double> tryUniform(double min, double max) const;
  std::optional<int> tryUniformInt(int min, int max) const;

  // Statistical distributions
  double exponential(double lambda = 1.0) const;
  std::optional<double> tryExponential(double lambda) const;

  double normal(double mean = 0.0, double stddev = 1.0) const;
  std::optional<double> tryNormal(double mean, double stddev) const;

  double gamma(double alpha, double beta = 1.0) const;
  std::optional<double> tryGamma(double alpha, double beta) const;

  double chi_squared(double degrees_of_freedom) const;
  double poisson(double mean) const;

  // Monte Carlo specific sampling
  Vector3D isotropicDirection() const;
  Vector3D cosineWeightedDirection(const Vector3D &normal) const;
  Vector3D uniformHemisphere(const Vector3D &normal) const;

  // Energy sampling
  double maxwellBoltzmann(double temperature) const;
  double powerLaw(double alpha, double x_min, double x_max) const;
  std::optional<double> tryPowerLaw(double alpha, double x_min,
                                    double x_max) const;

  // Discrete sampling
  template <typename T>
  std::size_t discreteSample(const std::vector<T> &weights) const
  {
    if(weights.empty())
      return 0;

    std::lock_guard<std::mutex> lock(mutex_);
    std::discrete_distribution<std::size_t> dist(weights.begin(),
                                                 weights.end());
    return dist(generator_);
  }

  template <typename T>
  std::optional<std::size_t>
  tryDiscreteSample(const std::vector<T> &weights) const
  {
    if(weights.empty())
      return std::nullopt;

    // Check for negative weights
    for(const auto &w : weights)
    {
      if(w < 0)
        return std::nullopt;
    }

    return discreteSample(weights);
  }

  // Array sampling utilities
  template <std::size_t N>
  std::size_t discreteSampleArray(const std::array<double, N> &weights) const
  {
    std::lock_guard<std::mutex> lock(mutex_);
    std::discrete_distribution<std::size_t> dist(weights.begin(),
                                                 weights.end());
    return dist(generator_);
  }

  // Selection and shuffling
  template <typename Iterator>
  Iterator selectRandom(Iterator first, Iterator last) const
  {
    auto dist = std::distance(first, last);
    if(dist <= 0)
      return last;

    auto advance_by = uniformSize(0, static_cast<std::size_t>(dist - 1));
    std::advance(first, advance_by);
    return first;
  }

  template <typename Container> void shuffle(Container &container) const
  {
    std::lock_guard<std::mutex> lock(mutex_);
    std::shuffle(container.begin(), container.end(), generator_);
  }

  // Physics-specific sampling
  double comptonScatteringAngle(double incident_energy) const;
  std::optional<double> tryComptonScatteringAngle(double incident_energy) const;

  double photoelectricEnergy(double binding_energy, double photon_energy) const;
  std::optional<double> tryPhotoelectricEnergy(double binding_energy,
                                               double photon_energy) const;

  // Advanced Monte Carlo techniques
  bool russianRoulette(double survival_probability) const;
  std::optional<bool> tryRussianRoulette(double survival_probability) const;

  std::size_t particleSplitting(double splitting_factor) const;
  std::optional<std::size_t>
  tryParticleSplitting(double splitting_factor) const;

  // Batch sampling for performance
  std::vector<double> uniformBatch(std::size_t count) const;
  std::vector<double> uniformBatch(std::size_t count, double min,
                                   double max) const;
  std::vector<double> exponentialBatch(std::size_t count, double lambda) const;
  std::vector<Vector3D> isotropicDirectionBatch(std::size_t count) const;

  // State management for reproducibility
  struct State
  {
    GeneratorType::result_type generator_state;
    SeedType seed;
  };

  State getState() const;
  void setState(const State &state);

  // Thread-local generators for high performance parallel Monte Carlo
  static thread_local std::unique_ptr<RandomNumberGenerator> thread_local_rng_;
  static RandomNumberGenerator &getThreadLocal();
  static void initThreadLocal(SeedType base_seed = 0);

  // Statistics and diagnostics
  struct Statistics
  {
    std::size_t uniform_calls = 0;
    std::size_t exponential_calls = 0;
    std::size_t normal_calls = 0;
    std::size_t direction_calls = 0;
    std::size_t total_calls = 0;
    double generation_time = 0.0; // In seconds
  };

  Statistics getStatistics() const;
  void resetStatistics();
  void enableStatistics(bool enable = true);

  // Validation and testing
  bool validateUniformity(std::size_t samples = 100000,
                          double tolerance = 0.01) const;
  bool validateExponential(double lambda, std::size_t samples = 100000,
                           double tolerance = 0.05) const;

  // Utility methods
  std::string toString() const;
  bool isSeeded() const { return current_seed_ != 0; }

  // Comparison operators
  bool operator==(const RandomNumberGenerator &other) const;
  bool operator!=(const RandomNumberGenerator &other) const;

  // Static utility methods
  static constexpr SeedType MAX_SEED = std::numeric_limits<SeedType>::max();
  static constexpr SeedType MIN_SEED = 1;

  static bool isValidSeed(SeedType seed)
  {
    return seed >= MIN_SEED && seed <= MAX_SEED;
  }
  static SeedType getTimeSeed();
  static SeedType getCryptoSeed(); // If available on platform

  // Global random number generator for convenience
  static RandomNumberGenerator &global();
  static void setGlobalSeed(SeedType seed);

private:
  mutable Statistics statistics_;
  mutable bool enable_statistics_;

  void initializeGenerator();
  void updateStatistics(const std::string &method) const;
  static SeedType getNextGlobalSeed();

  // Helper methods for complex sampling
  double sampleRejection(std::function<double()> proposal,
                         std::function<double(double)> target,
                         double max_target_value) const;
};