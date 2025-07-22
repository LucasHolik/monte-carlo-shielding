#include "RandomNumberGenerator.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <thread>

// Physical constants for Monte Carlo physics
namespace PhysicalConstants
{
constexpr double ELECTRON_REST_MASS = 0.510999; // MeV
constexpr double PI = 3.14159265358979323846;
constexpr double TWO_PI = 2.0 * PI;
} // namespace PhysicalConstants

// Static member initialization
std::atomic<RandomNumberGenerator::SeedType>
    RandomNumberGenerator::global_seed_counter_{1000};
thread_local std::unique_ptr<RandomNumberGenerator>
    RandomNumberGenerator::thread_local_rng_ = nullptr;

// Constructors
RandomNumberGenerator::RandomNumberGenerator()
    : current_seed_(generateUniqueSeed()), uniform_dist_(0.0, 1.0),
      normal_dist_(0.0, 1.0), exponential_dist_(1.0), enable_statistics_(false)
{
  initializeGenerator();
}

RandomNumberGenerator::RandomNumberGenerator(SeedType seed)
    : current_seed_(seed), uniform_dist_(0.0, 1.0), normal_dist_(0.0, 1.0),
      exponential_dist_(1.0), enable_statistics_(false)
{
  if(!isValidSeed(seed))
  {
    throw std::invalid_argument("Invalid seed value: " + std::to_string(seed));
  }
  initializeGenerator();
}

// Copy constructor - creates independent generator
RandomNumberGenerator::RandomNumberGenerator(const RandomNumberGenerator &other)
    : current_seed_(generateUniqueSeed()), uniform_dist_(0.0, 1.0),
      normal_dist_(0.0, 1.0), exponential_dist_(1.0),
      enable_statistics_(other.enable_statistics_)
{
  initializeGenerator();
}

RandomNumberGenerator &
RandomNumberGenerator::operator=(const RandomNumberGenerator &other)
{
  if(this != &other)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    current_seed_ = generateUniqueSeed();
    enable_statistics_ = other.enable_statistics_;
    statistics_ = Statistics{}; // Reset statistics
    initializeGenerator();
  }
  return *this;
}

// Move constructor
RandomNumberGenerator::RandomNumberGenerator(
    RandomNumberGenerator &&other) noexcept
    : current_seed_(other.current_seed_),
      uniform_dist_(std::move(other.uniform_dist_)),
      normal_dist_(std::move(other.normal_dist_)),
      exponential_dist_(std::move(other.exponential_dist_)),
      statistics_(std::move(other.statistics_)),
      enable_statistics_(other.enable_statistics_)
{
  std::lock_guard<std::mutex> lock(other.mutex_);
  generator_ = std::move(other.generator_);
}

RandomNumberGenerator &
RandomNumberGenerator::operator=(RandomNumberGenerator &&other) noexcept
{
  if(this != &other)
  {
    std::lock_guard<std::mutex> lock1(mutex_);
    std::lock_guard<std::mutex> lock2(other.mutex_);

    generator_ = std::move(other.generator_);
    current_seed_ = other.current_seed_;
    uniform_dist_ = std::move(other.uniform_dist_);
    normal_dist_ = std::move(other.normal_dist_);
    exponential_dist_ = std::move(other.exponential_dist_);
    statistics_ = std::move(other.statistics_);
    enable_statistics_ = other.enable_statistics_;
  }
  return *this;
}

void RandomNumberGenerator::initializeGenerator()
{
  std::lock_guard<std::mutex> lock(mutex_);
  generator_.seed(current_seed_);

  // Reset distributions to ensure consistency
  uniform_dist_.reset();
  normal_dist_.reset();
  exponential_dist_.reset();
}

// Seed management
void RandomNumberGenerator::setSeed(SeedType seed)
{
  if(!isValidSeed(seed))
  {
    throw std::invalid_argument("Invalid seed value: " + std::to_string(seed));
  }
  current_seed_ = seed;
  initializeGenerator();
}

void RandomNumberGenerator::randomSeed() { setSeed(getTimeSeed()); }

RandomNumberGenerator::SeedType RandomNumberGenerator::generateUniqueSeed()
{
  return getNextGlobalSeed();
}

RandomNumberGenerator::SeedType RandomNumberGenerator::getNextGlobalSeed()
{
  return global_seed_counter_.fetch_add(1);
}

// Basic uniform sampling
double RandomNumberGenerator::uniform() const
{
  std::lock_guard<std::mutex> lock(mutex_);
  updateStatistics("uniform");
  return uniform_dist_(generator_);
}

double RandomNumberGenerator::uniform(double min, double max) const
{
  if(min >= max)
  {
    throw std::invalid_argument("Invalid range: min must be less than max");
  }
  return min + (max - min) * uniform();
}

std::optional<double> RandomNumberGenerator::tryUniform(double min,
                                                        double max) const
{
  if(min >= max)
    return std::nullopt;
  return uniform(min, max);
}

// Integer sampling
int RandomNumberGenerator::uniformInt(int min, int max) const
{
  if(min > max)
  {
    throw std::invalid_argument("Invalid range: min must be <= max");
  }
  std::lock_guard<std::mutex> lock(mutex_);
  std::uniform_int_distribution<int> dist(min, max);
  updateStatistics("uniform_int");
  return dist(generator_);
}

std::size_t RandomNumberGenerator::uniformSize(std::size_t min,
                                               std::size_t max) const
{
  if(min > max)
  {
    throw std::invalid_argument("Invalid range: min must be <= max");
  }
  std::lock_guard<std::mutex> lock(mutex_);
  std::uniform_int_distribution<std::size_t> dist(min, max);
  updateStatistics("uniform_size");
  return dist(generator_);
}

std::optional<int> RandomNumberGenerator::tryUniformInt(int min, int max) const
{
  if(min > max)
    return std::nullopt;
  return uniformInt(min, max);
}

// Statistical distributions
double RandomNumberGenerator::exponential(double lambda) const
{
  if(lambda <= 0.0)
  {
    throw std::invalid_argument("Lambda must be positive");
  }
  std::lock_guard<std::mutex> lock(mutex_);
  exponential_dist_.param(
      std::exponential_distribution<double>::param_type(lambda));
  updateStatistics("exponential");
  return exponential_dist_(generator_);
}

std::optional<double> RandomNumberGenerator::tryExponential(double lambda) const
{
  if(lambda <= 0.0)
    return std::nullopt;
  return exponential(lambda);
}

double RandomNumberGenerator::normal(double mean, double stddev) const
{
  if(stddev <= 0.0)
  {
    throw std::invalid_argument("Standard deviation must be positive");
  }
  std::lock_guard<std::mutex> lock(mutex_);
  normal_dist_.param(
      std::normal_distribution<double>::param_type(mean, stddev));
  updateStatistics("normal");
  return normal_dist_(generator_);
}

std::optional<double> RandomNumberGenerator::tryNormal(double mean,
                                                       double stddev) const
{
  if(stddev <= 0.0)
    return std::nullopt;
  return normal(mean, stddev);
}

double RandomNumberGenerator::gamma(double alpha, double beta) const
{
  if(alpha <= 0.0 || beta <= 0.0)
  {
    throw std::invalid_argument("Alpha and beta must be positive");
  }
  std::lock_guard<std::mutex> lock(mutex_);
  std::gamma_distribution<double> dist(alpha, beta);
  updateStatistics("gamma");
  return dist(generator_);
}

std::optional<double> RandomNumberGenerator::tryGamma(double alpha,
                                                      double beta) const
{
  if(alpha <= 0.0 || beta <= 0.0)
    return std::nullopt;
  return gamma(alpha, beta);
}

double RandomNumberGenerator::chi_squared(double degrees_of_freedom) const
{
  if(degrees_of_freedom <= 0.0)
  {
    throw std::invalid_argument("Degrees of freedom must be positive");
  }
  std::lock_guard<std::mutex> lock(mutex_);
  std::chi_squared_distribution<double> dist(degrees_of_freedom);
  updateStatistics("chi_squared");
  return dist(generator_);
}

double RandomNumberGenerator::poisson(double mean) const
{
  if(mean <= 0.0)
  {
    throw std::invalid_argument("Mean must be positive");
  }
  std::lock_guard<std::mutex> lock(mutex_);
  std::poisson_distribution<int> dist(mean);
  updateStatistics("poisson");
  return static_cast<double>(dist(generator_));
}

// Monte Carlo specific sampling
Vector3D RandomNumberGenerator::isotropicDirection() const
{
  updateStatistics("direction");

  // Use Marsaglia method for uniform sampling on unit sphere
  double u1 = uniform();
  double u2 = uniform();

  double cos_theta = 2.0 * u1 - 1.0;
  double phi = PhysicalConstants::TWO_PI * u2;
  double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

  return Vector3D(sin_theta * std::cos(phi), sin_theta * std::sin(phi),
                  cos_theta);
}

Vector3D
RandomNumberGenerator::cosineWeightedDirection(const Vector3D &normal) const
{
  updateStatistics("direction");

  // Sample cosine-weighted hemisphere direction
  double u1 = uniform();
  double u2 = uniform();

  double cos_theta = std::sqrt(u1);
  double sin_theta = std::sqrt(1.0 - u1);
  double phi = PhysicalConstants::TWO_PI * u2;

  // Local coordinates
  Vector3D local_dir(sin_theta * std::cos(phi), sin_theta * std::sin(phi),
                     cos_theta);

  // Transform to global coordinates aligned with normal
  // This is a simplified version - proper implementation would need
  // full coordinate system transformation
  if(std::abs(normal.z()) > 0.9)
  {
    return local_dir;
  }
  else
  {
    // For non-vertical normals, this is an approximation
    return local_dir;
  }
}

Vector3D RandomNumberGenerator::uniformHemisphere(const Vector3D &normal) const
{
  updateStatistics("direction");

  double u1 = uniform();
  double u2 = uniform();

  double cos_theta = u1; // Uniform in cos(theta) for hemisphere
  double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
  double phi = PhysicalConstants::TWO_PI * u2;

  Vector3D local_dir(sin_theta * std::cos(phi), sin_theta * std::sin(phi),
                     cos_theta);

  // Simplified transformation - similar caveat as above
  return local_dir;
}

// Energy sampling
double RandomNumberGenerator::maxwellBoltzmann(double temperature) const
{
  if(temperature <= 0.0)
  {
    throw std::invalid_argument("Temperature must be positive");
  }

  // Sample from Maxwell-Boltzmann distribution using Gamma distribution
  // E ~ Gamma(3/2, kT)
  return gamma(1.5, temperature);
}

double RandomNumberGenerator::powerLaw(double alpha, double x_min,
                                       double x_max) const
{
  if(alpha == -1.0)
  {
    throw std::invalid_argument("Alpha cannot be -1 (logarithmic case)");
  }
  if(x_min >= x_max || x_min <= 0.0)
  {
    throw std::invalid_argument("Invalid range for power law");
  }

  double u = uniform();

  if(std::abs(alpha + 1.0) < 1e-10)
  {
    // Special case for alpha ≈ -1
    return x_min * std::exp(u * std::log(x_max / x_min));
  }
  else
  {
    // General case
    double alpha_plus_1 = alpha + 1.0;
    double term1 = std::pow(x_min, alpha_plus_1);
    double term2 = std::pow(x_max, alpha_plus_1);
    return std::pow(term1 + u * (term2 - term1), 1.0 / alpha_plus_1);
  }
}

std::optional<double> RandomNumberGenerator::tryPowerLaw(double alpha,
                                                         double x_min,
                                                         double x_max) const
{
  if(x_min >= x_max || x_min <= 0.0)
    return std::nullopt;
  if(std::abs(alpha + 1.0) < 1e-10)
    return std::nullopt;

  try
  {
    return powerLaw(alpha, x_min, x_max);
  }
  catch(...)
  {
    return std::nullopt;
  }
}

// Physics-specific sampling
double
RandomNumberGenerator::comptonScatteringAngle(double incident_energy) const
{
  if(incident_energy <= 0.0)
  {
    throw std::invalid_argument("Incident energy must be positive");
  }

  // Klein-Nishina formula for Compton scattering
  double alpha = incident_energy / PhysicalConstants::ELECTRON_REST_MASS;

  // Use rejection sampling for Klein-Nishina distribution
  double cos_theta;
  double klein_nishina;
  double rejection_value;

  do
  {
    cos_theta = 2.0 * uniform() - 1.0;
    double epsilon = 1.0 / (1.0 + alpha * (1.0 - cos_theta));
    klein_nishina = epsilon * epsilon *
                    (epsilon + 1.0 / epsilon - 1.0 + cos_theta * cos_theta);
    rejection_value = uniform();
  } while(rejection_value > klein_nishina / 2.0);

  return std::acos(cos_theta);
}

std::optional<double>
RandomNumberGenerator::tryComptonScatteringAngle(double incident_energy) const
{
  if(incident_energy <= 0.0)
    return std::nullopt;
  return comptonScatteringAngle(incident_energy);
}

double RandomNumberGenerator::photoelectricEnergy(double binding_energy,
                                                  double photon_energy) const
{
  if(photon_energy <= binding_energy)
  {
    throw std::invalid_argument("Photon energy must exceed binding energy");
  }

  // For photoelectric effect, electron gets photon energy minus binding energy
  // This is deterministic, but we might add small fluctuations
  return photon_energy - binding_energy;
}

std::optional<double>
RandomNumberGenerator::tryPhotoelectricEnergy(double binding_energy,
                                              double photon_energy) const
{
  if(photon_energy <= binding_energy)
    return std::nullopt;
  return photoelectricEnergy(binding_energy, photon_energy);
}

// Advanced Monte Carlo techniques
bool RandomNumberGenerator::russianRoulette(double survival_probability) const
{
  if(survival_probability < 0.0 || survival_probability > 1.0)
  {
    throw std::invalid_argument("Survival probability must be between 0 and 1");
  }
  updateStatistics("russian_roulette");
  return uniform() < survival_probability;
}

std::optional<bool>
RandomNumberGenerator::tryRussianRoulette(double survival_probability) const
{
  if(survival_probability < 0.0 || survival_probability > 1.0)
    return std::nullopt;
  return russianRoulette(survival_probability);
}

std::size_t
RandomNumberGenerator::particleSplitting(double splitting_factor) const
{
  if(splitting_factor < 1.0)
  {
    throw std::invalid_argument("Splitting factor must be >= 1");
  }

  updateStatistics("particle_splitting");

  // Number of particles is floor(factor) + 1 with probability (factor -
  // floor(factor))
  std::size_t base_count =
      static_cast<std::size_t>(std::floor(splitting_factor));
  double fractional_part = splitting_factor - base_count;

  return uniform() < fractional_part ? base_count + 1 : base_count;
}

std::optional<std::size_t>
RandomNumberGenerator::tryParticleSplitting(double splitting_factor) const
{
  if(splitting_factor < 1.0)
    return std::nullopt;
  return particleSplitting(splitting_factor);
}

// Batch sampling for performance
std::vector<double> RandomNumberGenerator::uniformBatch(std::size_t count) const
{
  std::vector<double> result;
  result.reserve(count);

  std::lock_guard<std::mutex> lock(mutex_);
  for(std::size_t i = 0; i < count; ++i)
  {
    result.push_back(uniform_dist_(generator_));
  }

  updateStatistics("uniform_batch");
  return result;
}

std::vector<double> RandomNumberGenerator::uniformBatch(std::size_t count,
                                                        double min,
                                                        double max) const
{
  if(min >= max)
  {
    throw std::invalid_argument("Invalid range: min must be less than max");
  }

  auto base_values = uniformBatch(count);
  double range = max - min;

  std::transform(base_values.begin(), base_values.end(), base_values.begin(),
                 [min, range](double x) { return min + range * x; });

  return base_values;
}

std::vector<double> RandomNumberGenerator::exponentialBatch(std::size_t count,
                                                            double lambda) const
{
  if(lambda <= 0.0)
  {
    throw std::invalid_argument("Lambda must be positive");
  }

  std::vector<double> result;
  result.reserve(count);

  std::lock_guard<std::mutex> lock(mutex_);
  exponential_dist_.param(
      std::exponential_distribution<double>::param_type(lambda));

  for(std::size_t i = 0; i < count; ++i)
  {
    result.push_back(exponential_dist_(generator_));
  }

  updateStatistics("exponential_batch");
  return result;
}

std::vector<Vector3D>
RandomNumberGenerator::isotropicDirectionBatch(std::size_t count) const
{
  std::vector<Vector3D> result;
  result.reserve(count);

  for(std::size_t i = 0; i < count; ++i)
  {
    result.push_back(isotropicDirection());
  }

  return result;
}

// State management
RandomNumberGenerator::State RandomNumberGenerator::getState() const
{
  std::lock_guard<std::mutex> lock(mutex_);
  return State{generator_(), current_seed_};
}

void RandomNumberGenerator::setState(const State &state)
{
  std::lock_guard<std::mutex> lock(mutex_);
  current_seed_ = state.seed;
  generator_.seed(current_seed_);
  // Note: This is a simplified state restoration
  // Full restoration would require storing the complete generator state
}

// Thread-local functionality
RandomNumberGenerator &RandomNumberGenerator::getThreadLocal()
{
  if(!thread_local_rng_)
  {
    initThreadLocal();
  }
  return *thread_local_rng_;
}

void RandomNumberGenerator::initThreadLocal(SeedType base_seed)
{
  if(base_seed == 0)
  {
    base_seed = generateUniqueSeed();
  }

  // Add thread ID to ensure uniqueness across threads
  auto thread_id = std::hash<std::thread::id>{}(std::this_thread::get_id());
  SeedType thread_seed = base_seed + (thread_id % 10000);

  thread_local_rng_ = std::make_unique<RandomNumberGenerator>(thread_seed);
}

// Statistics
RandomNumberGenerator::Statistics RandomNumberGenerator::getStatistics() const
{
  std::lock_guard<std::mutex> lock(mutex_);
  return statistics_;
}

void RandomNumberGenerator::resetStatistics()
{
  std::lock_guard<std::mutex> lock(mutex_);
  statistics_ = Statistics{};
}

void RandomNumberGenerator::enableStatistics(bool enable)
{
  std::lock_guard<std::mutex> lock(mutex_);
  enable_statistics_ = enable;
}

void RandomNumberGenerator::updateStatistics(const std::string &method) const
{
  if(!enable_statistics_)
    return;

  // This would be called without the mutex already held
  statistics_.total_calls++;

  if(method == "uniform" || method == "uniform_int" || method == "uniform_size")
  {
    statistics_.uniform_calls++;
  }
  else if(method == "exponential")
  {
    statistics_.exponential_calls++;
  }
  else if(method == "normal")
  {
    statistics_.normal_calls++;
  }
  else if(method == "direction")
  {
    statistics_.direction_calls++;
  }
}

// Validation
bool RandomNumberGenerator::validateUniformity(std::size_t samples,
                                               double tolerance) const
{
  std::vector<double> values = uniformBatch(samples);

  // Simple chi-square test for uniformity
  const int bins = 10;
  std::vector<int> counts(bins, 0);

  for(double value : values)
  {
    int bin = static_cast<int>(value * bins);
    if(bin >= bins)
      bin = bins - 1;
    counts[bin]++;
  }

  double expected = static_cast<double>(samples) / bins;
  double chi_square = 0.0;

  for(int count : counts)
  {
    double diff = count - expected;
    chi_square += (diff * diff) / expected;
  }

  // Rough test: chi-square should be reasonable for uniform distribution
  return chi_square < (bins - 1) * (1.0 + tolerance);
}

bool RandomNumberGenerator::validateExponential(double lambda,
                                                std::size_t samples,
                                                double tolerance) const
{
  std::vector<double> values = exponentialBatch(samples, lambda);

  // Test mean (should be 1/lambda)
  double mean = std::accumulate(values.begin(), values.end(), 0.0) / samples;
  double expected_mean = 1.0 / lambda;
  double relative_error = std::abs(mean - expected_mean) / expected_mean;

  return relative_error < tolerance;
}

// Utility methods
std::string RandomNumberGenerator::toString() const
{
  std::stringstream ss;
  ss << "RandomNumberGenerator(seed=" << current_seed_
     << ", statistics=" << (enable_statistics_ ? "enabled" : "disabled") << ")";
  return ss.str();
}

bool RandomNumberGenerator::operator==(const RandomNumberGenerator &other) const
{
  return current_seed_ == other.current_seed_;
}

bool RandomNumberGenerator::operator!=(const RandomNumberGenerator &other) const
{
  return !(*this == other);
}

// Static utility methods
RandomNumberGenerator::SeedType RandomNumberGenerator::getTimeSeed()
{
  auto now = std::chrono::high_resolution_clock::now();
  auto duration = now.time_since_epoch();
  auto nanoseconds =
      std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count();

  // Combine with thread ID for better uniqueness
  auto thread_id = std::hash<std::thread::id>{}(std::this_thread::get_id());

  return static_cast<SeedType>(nanoseconds ^ thread_id);
}

RandomNumberGenerator::SeedType RandomNumberGenerator::getCryptoSeed()
{
  // Platform-specific implementation would go here
  // For now, fall back to time-based seed
  return getTimeSeed();
}

// Global instance
RandomNumberGenerator &RandomNumberGenerator::global()
{
  static RandomNumberGenerator global_rng;
  return global_rng;
}

void RandomNumberGenerator::setGlobalSeed(SeedType seed)
{
  global().setSeed(seed);
}