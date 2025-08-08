#include "RandomNumberQuality.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>

namespace
{
// Statistical constants
constexpr double PI = 3.14159265358979323846;
constexpr double TWO_PI = 2.0 * PI;
const double SQRT_2PI = std::sqrt(TWO_PI);

// Compute gamma function approximation (Stirling's approximation for large n)
double gammaFunction(double x)
{
  if(x < 1.0)
    return gammaFunction(x + 1.0) / x;
  return std::sqrt(TWO_PI / x) * std::pow(x / std::exp(1.0), x);
}

// Incomplete gamma function (simplified implementation)
double incompleteGamma(double a, double x)
{
  if(x == 0.0)
    return 0.0;
  if(a == 0.0)
    return 1.0;

  // Use series expansion for small x
  double sum = 1.0;
  double term = 1.0;
  for(int n = 1; n < 100; ++n)
  {
    term *= x / (a + n - 1);
    sum += term;
    if(std::abs(term) < 1e-10)
      break;
  }

  return std::pow(x, a) * std::exp(-x) * sum / gammaFunction(a);
}

} // anonymous namespace

RandomNumberQuality::RandomNumberQuality(RandomNumberGenerator &rng,
                                         bool verbose)
    : rng_(rng), verbose_(false)  // Always silent in production
{}

// ============================================================================
// MEMBER FUNCTION IMPLEMENTATIONS FOR CRITICAL VALUES
// ============================================================================

double RandomNumberQuality::getChiSquareCriticalValue(std::size_t df,
                                                      double alpha)
{
  // Simplified critical values for chi-square distribution
  // In production, would use proper statistical library
  if(alpha > 0.1)
    return 2.0 * df;
  if(alpha > 0.05)
    return 2.5 * df;
  if(alpha > 0.01)
    return 3.0 * df;
  return 3.5 * df;
}

double RandomNumberQuality::getKSCriticalValue(std::size_t n, double alpha)
{
  // Kolmogorov-Smirnov critical values
  double sqrt_n = std::sqrt(static_cast<double>(n));
  if(alpha > 0.1)
    return 1.22 / sqrt_n;
  if(alpha > 0.05)
    return 1.36 / sqrt_n;
  if(alpha > 0.01)
    return 1.63 / sqrt_n;
  return 1.95 / sqrt_n;
}

double RandomNumberQuality::getNormalCriticalValue(double alpha)
{
  // Standard normal critical values
  if(alpha > 0.1)
    return 1.645;
  if(alpha > 0.05)
    return 1.96;
  if(alpha > 0.01)
    return 2.576;
  return 3.291;
}

// ============================================================================
// UNIFORMITY TESTS
// ============================================================================

TestResult RandomNumberQuality::chiSquareTest(std::size_t n_samples,
                                              std::size_t n_bins,
                                              double confidence_level)
{
  // Running chi-square test silently

  // Generate uniform samples
  std::vector<double> samples;
  samples.reserve(n_samples);
  for(std::size_t i = 0; i < n_samples; ++i)
  {
    samples.push_back(rng_.uniform());
  }

  // Count frequencies in bins
  std::vector<int> observed(n_bins, 0);
  for(double sample : samples)
  {
    std::size_t bin = static_cast<std::size_t>(sample * n_bins);
    if(bin >= n_bins)
      bin = n_bins - 1;
    observed[bin]++;
  }

  // Expected frequency per bin for uniform distribution
  double expected = static_cast<double>(n_samples) / n_bins;
  std::vector<double> expected_vec(n_bins, expected);

  // Calculate chi-square statistic
  double chi_square = calculateChiSquareStatistic(observed, expected_vec);

  // Degrees of freedom
  std::size_t df = n_bins - 1;

  // Critical value
  double alpha = 1.0 - confidence_level;
  double critical_value = getChiSquareCriticalValue(df, alpha);

  bool passed = chi_square <= critical_value;

  // Approximate p-value using incomplete gamma function
  double p_value = 1.0 - incompleteGamma(df / 2.0, chi_square / 2.0);

  TestResult result(passed, chi_square, p_value,
                    "Chi-square test for uniformity");

  result.details["critical_value"] = critical_value;
  result.details["degrees_of_freedom"] = static_cast<double>(df);
  result.details["expected_frequency"] = expected;

  // Test completed silently

  return result;
}

TestResult RandomNumberQuality::kolmogorovSmirnovTest(std::size_t n_samples,
                                                      double confidence_level)
{
  // Running KS test silently

  // Generate and sort samples
  std::vector<double> samples;
  samples.reserve(n_samples);
  for(std::size_t i = 0; i < n_samples; ++i)
  {
    samples.push_back(rng_.uniform());
  }
  std::sort(samples.begin(), samples.end());

  // Calculate KS statistic
  double ks_statistic = calculateKSStatistic(samples);

  // Critical value
  double alpha = 1.0 - confidence_level;
  double critical_value = getKSCriticalValue(n_samples, alpha);

  bool passed = ks_statistic <= critical_value;

  // Approximate p-value (simplified)
  double p_value =
      2.0 * std::exp(-2.0 * n_samples * ks_statistic * ks_statistic);
  p_value = std::max(0.0, std::min(1.0, p_value));

  TestResult result(passed, ks_statistic, p_value,
                    "Kolmogorov-Smirnov test for uniformity");

  result.details["critical_value"] = critical_value;
  result.details["sample_size"] = static_cast<double>(n_samples);

  // KS test completed silently

  return result;
}

TestResult RandomNumberQuality::andersonDarlingTest(std::size_t n_samples)
{
  if(verbose_)
  {
    std::cout << "Running Anderson-Darling test with " << n_samples
              << " samples...\n";
  }

  // Generate and sort samples
  std::vector<double> samples;
  samples.reserve(n_samples);
  for(std::size_t i = 0; i < n_samples; ++i)
  {
    samples.push_back(rng_.uniform());
  }
  std::sort(samples.begin(), samples.end());

  // Calculate Anderson-Darling statistic
  double ad_statistic = 0.0;
  double n = static_cast<double>(n_samples);

  for(std::size_t i = 0; i < n_samples; ++i)
  {
    double u_i = samples[i];
    double u_n_minus_i = samples[n_samples - 1 - i];

    // Avoid log(0) by clamping values
    u_i = std::max(1e-10, std::min(1.0 - 1e-10, u_i));
    u_n_minus_i = std::max(1e-10, std::min(1.0 - 1e-10, u_n_minus_i));

    double term =
        (2.0 * (i + 1) - 1.0) * (std::log(u_i) + std::log(1.0 - u_n_minus_i));
    ad_statistic += term;
  }

  ad_statistic = -n - ad_statistic / n;

  // Critical value for 5% significance level (approximate)
  double critical_value = 2.492;
  bool passed = ad_statistic <= critical_value;

  // Approximate p-value (very simplified)
  double p_value = passed ? 0.1 : 0.01;

  TestResult result(passed, ad_statistic, p_value,
                    "Anderson-Darling test for uniformity");

  result.details["critical_value"] = critical_value;
  result.details["sample_size"] = n;

  if(verbose_)
  {
    std::cout << "Anderson-Darling statistic: " << ad_statistic
              << ", Critical value: " << critical_value << std::endl;
  }

  return result;
}

// ============================================================================
// INDEPENDENCE TESTS
// ============================================================================

TestResult RandomNumberQuality::serialCorrelationTest(std::size_t n_samples,
                                                      std::size_t lag)
{
  if(verbose_)
  {
    std::cout << "Running serial correlation test with lag " << lag << " and "
              << n_samples << " samples...\n";
  }

  // Generate samples
  std::vector<double> samples;
  samples.reserve(n_samples);
  for(std::size_t i = 0; i < n_samples; ++i)
  {
    samples.push_back(rng_.uniform());
  }

  // Calculate autocorrelation
  double correlation = calculateAutoCorrelation(samples, lag);

  // For large samples, correlation should be approximately normal
  // with mean 0 and variance 1/n
  double std_error = 1.0 / std::sqrt(static_cast<double>(n_samples));
  double z_score = std::abs(correlation) / std_error;

  // Critical value for 5% significance level
  double critical_value = 1.96;
  bool passed = z_score <= critical_value;

  // P-value for two-tailed test
  double p_value =
      2.0 * (1.0 - 0.5 * (1.0 + std::erf(z_score / std::sqrt(2.0))));

  TestResult result(passed, correlation, p_value, "Serial correlation test");

  result.details["z_score"] = z_score;
  result.details["standard_error"] = std_error;
  result.details["lag"] = static_cast<double>(lag);

  if(verbose_)
  {
    std::cout << "Correlation: " << correlation << ", Z-score: " << z_score
              << ", P-value: " << p_value << std::endl;
  }

  return result;
}

TestResult RandomNumberQuality::runTest(std::size_t n_samples)
{
  if(verbose_)
  {
    std::cout << "Running run test with " << n_samples << " samples...\n";
  }

  // Generate samples and convert to binary sequence
  std::vector<double> samples;
  samples.reserve(n_samples);
  for(std::size_t i = 0; i < n_samples; ++i)
  {
    samples.push_back(rng_.uniform());
  }

  // Count runs (sequences above/below median)
  std::size_t n_runs = countRuns(samples, 0.5);

  // For large n, number of runs is approximately normal
  double n = static_cast<double>(n_samples);
  double expected_runs = (n + 1.0) / 2.0;
  double variance_runs = (n - 1.0) / 4.0;
  double std_runs = std::sqrt(variance_runs);

  double z_score =
      std::abs(static_cast<double>(n_runs) - expected_runs) / std_runs;

  // Critical value for 5% significance level
  double critical_value = 1.96;
  bool passed = z_score <= critical_value;

  // P-value for two-tailed test
  double p_value =
      2.0 * (1.0 - 0.5 * (1.0 + std::erf(z_score / std::sqrt(2.0))));

  TestResult result(passed, static_cast<double>(n_runs), p_value,
                    "Run test for independence");

  result.details["expected_runs"] = expected_runs;
  result.details["z_score"] = z_score;
  result.details["variance"] = variance_runs;

  if(verbose_)
  {
    std::cout << "Observed runs: " << n_runs << ", Expected: " << expected_runs
              << ", Z-score: " << z_score << std::endl;
  }

  return result;
}

TestResult RandomNumberQuality::gapTest(std::size_t n_samples, double alpha,
                                        double beta)
{
  if(verbose_)
  {
    std::cout << "Running gap test with interval [" << alpha << ", " << beta
              << "]...\n";
  }

  // Generate samples
  std::vector<double> samples;
  samples.reserve(n_samples);
  for(std::size_t i = 0; i < n_samples; ++i)
  {
    samples.push_back(rng_.uniform());
  }

  // Count gaps between occurrences in interval [alpha, beta]
  std::vector<std::size_t> gaps;
  std::size_t current_gap = 0;
  bool in_gap = true;

  for(double sample : samples)
  {
    if(sample >= alpha && sample <= beta)
    {
      if(in_gap)
      {
        gaps.push_back(current_gap);
        current_gap = 0;
        in_gap = false;
      }
    }
    else
    {
      if(!in_gap)
      {
        in_gap = true;
        current_gap = 0;
      }
      current_gap++;
    }
  }

  if(gaps.empty())
  {
    return TestResult(false, 0.0, 0.0, "Gap test failed - no gaps found");
  }

  // Expected gap length for geometric distribution
  double p = beta - alpha; // Probability of success
  double expected_gap = (1.0 - p) / p;

  // Calculate sample mean
  double sample_mean =
      std::accumulate(gaps.begin(), gaps.end(), 0.0) / gaps.size();

  // Test statistic (simplified)
  double test_statistic =
      std::abs(sample_mean - expected_gap) / std::sqrt(expected_gap);

  bool passed = test_statistic <= 2.0; // Simplified threshold
  double p_value = passed ? 0.1 : 0.01;

  TestResult result(passed, test_statistic, p_value, "Gap test");

  result.details["expected_gap"] = expected_gap;
  result.details["sample_mean"] = sample_mean;
  result.details["num_gaps"] = static_cast<double>(gaps.size());
  result.details["interval_width"] = beta - alpha;

  if(verbose_)
  {
    std::cout << "Sample mean gap: " << sample_mean
              << ", Expected: " << expected_gap
              << ", Test statistic: " << test_statistic << std::endl;
  }

  return result;
}

// ============================================================================
// DISTRIBUTION-SPECIFIC TESTS
// ============================================================================

TestResult
RandomNumberQuality::exponentialDistributionTest(double lambda,
                                                 std::size_t n_samples)
{
  if(verbose_)
  {
    std::cout << "Testing exponential distribution with λ=" << lambda << " and "
              << n_samples << " samples...\n";
  }

  // Generate exponential samples
  std::vector<double> samples;
  samples.reserve(n_samples);
  for(std::size_t i = 0; i < n_samples; ++i)
  {
    samples.push_back(rng_.exponential(lambda));
  }

  // Test sample mean (should be 1/lambda)
  double sample_mean =
      std::accumulate(samples.begin(), samples.end(), 0.0) / n_samples;
  double expected_mean = 1.0 / lambda;
  double mean_error = std::abs(sample_mean - expected_mean) / expected_mean;

  // Test sample variance (should be 1/lambda²)
  double variance = 0.0;
  for(double sample : samples)
  {
    double diff = sample - sample_mean;
    variance += diff * diff;
  }
  variance /= (n_samples - 1);

  double expected_variance = 1.0 / (lambda * lambda);
  double variance_error =
      std::abs(variance - expected_variance) / expected_variance;

  // Combined test statistic
  double test_statistic = std::max(mean_error, variance_error);

  bool passed = test_statistic < 0.1; // 10% tolerance
  double p_value = passed ? 0.1 : 0.01;

  TestResult result(passed, test_statistic, p_value,
                    "Exponential distribution validation");

  result.details["sample_mean"] = sample_mean;
  result.details["expected_mean"] = expected_mean;
  result.details["sample_variance"] = variance;
  result.details["expected_variance"] = expected_variance;
  result.details["mean_error"] = mean_error;
  result.details["variance_error"] = variance_error;

  if(verbose_)
  {
    std::cout << "Sample mean: " << sample_mean
              << " (expected: " << expected_mean << ")\n";
    std::cout << "Sample variance: " << variance
              << " (expected: " << expected_variance << ")\n";
    std::cout << "Test statistic: " << test_statistic << std::endl;
  }

  return result;
}

TestResult RandomNumberQuality::normalDistributionTest(double mean,
                                                       double stddev,
                                                       std::size_t n_samples)
{
  if(verbose_)
  {
    std::cout << "Testing normal distribution with μ=" << mean
              << ", σ=" << stddev << " and " << n_samples << " samples...\n";
  }

  // Generate normal samples
  std::vector<double> samples;
  samples.reserve(n_samples);
  for(std::size_t i = 0; i < n_samples; ++i)
  {
    samples.push_back(rng_.normal(mean, stddev));
  }

  // Test sample mean
  double sample_mean =
      std::accumulate(samples.begin(), samples.end(), 0.0) / n_samples;
  double mean_error =
      std::abs(sample_mean - mean) / (stddev / std::sqrt(n_samples));

  // Test sample standard deviation
  double variance = 0.0;
  for(double sample : samples)
  {
    double diff = sample - sample_mean;
    variance += diff * diff;
  }
  variance /= (n_samples - 1);
  double sample_stddev = std::sqrt(variance);
  double stddev_error = std::abs(sample_stddev - stddev) / stddev;

  // Combined test statistic
  double test_statistic = std::max(mean_error, stddev_error * 3.0);

  bool passed = (mean_error < 3.0) &&
                (stddev_error < 0.2); // z-score < 3, 20% tolerance on stddev
  double p_value = passed ? 0.1 : 0.01;

  TestResult result(passed, test_statistic, p_value,
                    "Normal distribution validation");

  result.details["sample_mean"] = sample_mean;
  result.details["expected_mean"] = mean;
  result.details["sample_stddev"] = sample_stddev;
  result.details["expected_stddev"] = stddev;
  result.details["mean_z_score"] = mean_error;
  result.details["stddev_error"] = stddev_error;

  if(verbose_)
  {
    std::cout << "Sample mean: " << sample_mean << " (expected: " << mean
              << ")\n";
    std::cout << "Sample stddev: " << sample_stddev << " (expected: " << stddev
              << ")\n";
    std::cout << "Mean z-score: " << mean_error << std::endl;
  }

  return result;
}

TestResult RandomNumberQuality::isotropicDirectionTest(std::size_t n_samples)
{
  if(verbose_)
  {
    std::cout << "Testing isotropic direction generation with " << n_samples
              << " samples...\n";
  }

  // Generate direction vectors
  std::vector<Vector3D> directions;
  directions.reserve(n_samples);
  for(std::size_t i = 0; i < n_samples; ++i)
  {
    directions.push_back(rng_.isotropicDirection());
  }

  // Test 1: All vectors should be unit vectors
  double max_length_error = 0.0;
  for(const auto &dir : directions)
  {
    double length = dir.magnitude();
    double error = std::abs(length - 1.0);
    max_length_error = std::max(max_length_error, error);
  }

  // Test 2: Components should have zero mean
  double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
  for(const auto &dir : directions)
  {
    sum_x += dir.x();
    sum_y += dir.y();
    sum_z += dir.z();
  }

  double mean_x = sum_x / n_samples;
  double mean_y = sum_y / n_samples;
  double mean_z = sum_z / n_samples;

  double mean_magnitude =
      std::sqrt(mean_x * mean_x + mean_y * mean_y + mean_z * mean_z);

  // Test 3: Check spherical uniformity using phi distribution
  std::vector<double> phi_values;
  phi_values.reserve(n_samples);
  for(const auto &dir : directions)
  {
    double phi = std::atan2(dir.y(), dir.x());
    if(phi < 0.0)
      phi += TWO_PI;
    phi_values.push_back(phi / TWO_PI); // Normalize to [0,1]
  }

  // Use KS test on phi values (should be uniform)
  std::sort(phi_values.begin(), phi_values.end());
  double ks_phi = calculateKSStatistic(phi_values);

  // Combined test statistic
  double test_statistic =
      std::max({max_length_error * 1000.0,
                mean_magnitude * std::sqrt(n_samples), ks_phi * 10.0});

  bool passed = (max_length_error < 1e-10) &&
                (mean_magnitude < 3.0 / std::sqrt(n_samples)) && (ks_phi < 0.1);

  double p_value = passed ? 0.1 : 0.01;

  TestResult result(passed, test_statistic, p_value,
                    "Isotropic direction validation");

  result.details["max_length_error"] = max_length_error;
  result.details["mean_magnitude"] = mean_magnitude;
  result.details["ks_phi_statistic"] = ks_phi;
  result.details["mean_x"] = mean_x;
  result.details["mean_y"] = mean_y;
  result.details["mean_z"] = mean_z;

  if(verbose_)
  {
    std::cout << "Max length error: " << max_length_error << "\n";
    std::cout << "Mean vector magnitude: " << mean_magnitude << "\n";
    std::cout << "KS statistic for phi: " << ks_phi << std::endl;
  }

  return result;
}

// ============================================================================
// MONTE CARLO SPECIFIC TESTS
// ============================================================================

TestResult RandomNumberQuality::randomWalkTest(std::size_t n_steps,
                                               std::size_t n_walks)
{
  if(verbose_)
  {
    std::cout << "Testing random walk with " << n_steps << " steps and "
              << n_walks << " walks...\n";
  }

  std::vector<double> final_distances;
  final_distances.reserve(n_walks);

  for(std::size_t walk = 0; walk < n_walks; ++walk)
  {
    Vector3D position(0.0, 0.0, 0.0);

    for(std::size_t step = 0; step < n_steps; ++step)
    {
      Vector3D direction = rng_.isotropicDirection();
      position = position + direction; // Unit step
    }

    final_distances.push_back(position.magnitude());
  }

  // For 3D isotropic random walk with unit steps, mean distance ≈ √(2n/3)
  // More precisely: E[R] = √(8n/(3π)) * Γ(3/2) = √(2n/3)
  double expected_mean_distance = std::sqrt(2.0 * n_steps / 3.0);
  double sample_mean =
      std::accumulate(final_distances.begin(), final_distances.end(), 0.0) /
      n_walks;

  double relative_error =
      std::abs(sample_mean - expected_mean_distance) / expected_mean_distance;

  bool passed =
      relative_error < 0.15; // 15% tolerance for statistical fluctuations
  double p_value = passed ? 0.1 : 0.01;

  TestResult result(passed, relative_error, p_value,
                    "Random walk diffusion test");

  result.details["sample_mean_distance"] = sample_mean;
  result.details["expected_mean_distance"] = expected_mean_distance;
  result.details["n_steps"] = static_cast<double>(n_steps);
  result.details["n_walks"] = static_cast<double>(n_walks);

  if(verbose_)
  {
    std::cout << "Sample mean distance: " << sample_mean
              << " (expected: " << expected_mean_distance << ")\n";
    std::cout << "Relative error: " << relative_error << std::endl;
  }

  return result;
}

TestResult RandomNumberQuality::particleSplittingTest(double splitting_factor,
                                                      std::size_t n_trials)
{
  if(verbose_)
  {
    std::cout << "Testing particle splitting with factor " << splitting_factor
              << " over " << n_trials << " trials...\n";
  }

  std::vector<std::size_t> split_counts;
  split_counts.reserve(n_trials);

  for(std::size_t i = 0; i < n_trials; ++i)
  {
    split_counts.push_back(rng_.particleSplitting(splitting_factor));
  }

  // Calculate sample mean
  double sample_mean =
      std::accumulate(split_counts.begin(), split_counts.end(), 0.0) / n_trials;
  double expected_mean = splitting_factor;
  double relative_error = std::abs(sample_mean - expected_mean) / expected_mean;

  bool passed = relative_error < 0.05; // 5% tolerance
  double p_value = passed ? 0.1 : 0.01;

  TestResult result(passed, relative_error, p_value,
                    "Particle splitting consistency test");

  result.details["sample_mean"] = sample_mean;
  result.details["expected_mean"] = expected_mean;
  result.details["splitting_factor"] = splitting_factor;
  result.details["n_trials"] = static_cast<double>(n_trials);

  if(verbose_)
  {
    std::cout << "Sample mean: " << sample_mean
              << " (expected: " << expected_mean << ")\n";
    std::cout << "Relative error: " << relative_error << std::endl;
  }

  return result;
}

TestResult RandomNumberQuality::russianRouletteTest(double survival_prob,
                                                    std::size_t n_trials)
{
  if(verbose_)
  {
    std::cout << "Testing Russian roulette with survival probability "
              << survival_prob << " over " << n_trials << " trials...\n";
  }

  std::size_t survival_count = 0;
  for(std::size_t i = 0; i < n_trials; ++i)
  {
    if(rng_.russianRoulette(survival_prob))
    {
      survival_count++;
    }
  }

  double sample_rate = static_cast<double>(survival_count) / n_trials;
  double expected_rate = survival_prob;
  double relative_error = std::abs(sample_rate - expected_rate) / expected_rate;

  bool passed = relative_error < 0.05; // 5% tolerance
  double p_value = passed ? 0.1 : 0.01;

  TestResult result(passed, relative_error, p_value,
                    "Russian roulette fairness test");

  result.details["sample_rate"] = sample_rate;
  result.details["expected_rate"] = expected_rate;
  result.details["survival_count"] = static_cast<double>(survival_count);
  result.details["n_trials"] = static_cast<double>(n_trials);

  if(verbose_)
  {
    std::cout << "Sample survival rate: " << sample_rate
              << " (expected: " << expected_rate << ")\n";
    std::cout << "Relative error: " << relative_error << std::endl;
  }

  return result;
}

// ============================================================================
// COMPREHENSIVE TEST SUITES
// ============================================================================

std::vector<TestResult>
RandomNumberQuality::runUniformityTests(std::size_t n_samples)
{
  if(verbose_)
  {
    std::cout << "\n=== Running Uniformity Tests ===\n";
  }

  std::vector<TestResult> results;
  results.push_back(chiSquareTest(n_samples));
  results.push_back(
      kolmogorovSmirnovTest(std::min(n_samples, std::size_t(1000))));
  results.push_back(
      andersonDarlingTest(std::min(n_samples, std::size_t(1000))));

  return results;
}

std::vector<TestResult>
RandomNumberQuality::runIndependenceTests(std::size_t n_samples)
{
  if(verbose_)
  {
    std::cout << "\n=== Running Independence Tests ===\n";
  }

  std::vector<TestResult> results;
  results.push_back(serialCorrelationTest(n_samples, 1));
  results.push_back(serialCorrelationTest(n_samples, 5));
  results.push_back(runTest(n_samples));
  results.push_back(gapTest(n_samples));

  return results;
}

std::vector<TestResult> RandomNumberQuality::runMonteCarloTests()
{
  if(verbose_)
  {
    std::cout << "\n=== Running Monte Carlo Specific Tests ===\n";
  }

  std::vector<TestResult> results;
  results.push_back(exponentialDistributionTest(1.0, 10000));
  results.push_back(exponentialDistributionTest(2.5, 10000));
  results.push_back(normalDistributionTest(0.0, 1.0, 10000));
  results.push_back(isotropicDirectionTest(5000));
  results.push_back(randomWalkTest(100, 1000));
  results.push_back(particleSplittingTest(2.5, 10000));
  results.push_back(russianRouletteTest(0.8, 10000));

  return results;
}

std::map<std::string, std::vector<TestResult>>
RandomNumberQuality::runCompleteTestSuite()
{
  if(verbose_)
  {
    std::cout << "\n========== COMPLETE RANDOM NUMBER QUALITY TEST SUITE "
                 "==========\n";
  }

  std::map<std::string, std::vector<TestResult>> all_results;

  all_results["Uniformity"] = runUniformityTests(10000);
  all_results["Independence"] = runIndependenceTests(10000);
  all_results["Monte Carlo"] = runMonteCarloTests();

  if(verbose_)
  {
    printTestReport(all_results);
  }

  return all_results;
}

// ============================================================================
// UTILITY METHODS
// ============================================================================

double RandomNumberQuality::calculateChiSquareStatistic(
    const std::vector<int> &observed, const std::vector<double> &expected)
{
  double chi_square = 0.0;
  for(std::size_t i = 0; i < observed.size(); ++i)
  {
    if(expected[i] > 0.0)
    {
      double diff = observed[i] - expected[i];
      chi_square += (diff * diff) / expected[i];
    }
  }
  return chi_square;
}

double
RandomNumberQuality::calculateKSStatistic(const std::vector<double> &samples)
{
  double max_diff = 0.0;
  std::size_t n = samples.size();

  for(std::size_t i = 0; i < n; ++i)
  {
    double empirical_cdf = static_cast<double>(i + 1) / n;
    double theoretical_cdf = samples[i]; // For uniform [0,1]
    double diff = std::abs(empirical_cdf - theoretical_cdf);
    max_diff = std::max(max_diff, diff);
  }

  return max_diff;
}

double
RandomNumberQuality::calculateAutoCorrelation(const std::vector<double> &data,
                                              std::size_t lag)
{
  if(lag >= data.size())
    return 0.0;

  std::size_t n = data.size() - lag;

  // Calculate means
  double mean1 = 0.0, mean2 = 0.0;
  for(std::size_t i = 0; i < n; ++i)
  {
    mean1 += data[i];
    mean2 += data[i + lag];
  }
  mean1 /= n;
  mean2 /= n;

  // Calculate correlation
  double numerator = 0.0;
  double denom1 = 0.0, denom2 = 0.0;

  for(std::size_t i = 0; i < n; ++i)
  {
    double diff1 = data[i] - mean1;
    double diff2 = data[i + lag] - mean2;
    numerator += diff1 * diff2;
    denom1 += diff1 * diff1;
    denom2 += diff2 * diff2;
  }

  double denominator = std::sqrt(denom1 * denom2);
  return (denominator > 0.0) ? numerator / denominator : 0.0;
}

std::size_t RandomNumberQuality::countRuns(const std::vector<double> &data,
                                           double threshold)
{
  if(data.empty())
    return 0;

  std::size_t runs = 1;
  bool above = data[0] > threshold;

  for(std::size_t i = 1; i < data.size(); ++i)
  {
    bool current_above = data[i] > threshold;
    if(current_above != above)
    {
      runs++;
      above = current_above;
    }
  }

  return runs;
}

std::string RandomNumberQuality::summarizeResult(const TestResult &result)
{
  std::stringstream ss;
  ss << result.description << ": " << (result.passed ? "PASS" : "FAIL")
     << " (statistic=" << std::fixed << std::setprecision(4)
     << result.test_statistic << ", p=" << result.p_value << ")";
  return ss.str();
}

void RandomNumberQuality::printTestReport(
    const std::map<std::string, std::vector<TestResult>> &results)
{
  std::size_t total_tests = 0;
  std::size_t passed_tests = 0;

  for(const auto &category : results)
  {
    std::cout << "\n--- " << category.first << " Tests ---\n";

    for(const auto &result : category.second)
    {
      total_tests++;
      if(result.passed)
        passed_tests++;

      std::cout << summarizeResult(result) << "\n";
    }
  }

  double pass_rate = static_cast<double>(passed_tests) / total_tests * 100.0;

  std::cout << "\n========== SUMMARY ==========\n";
  std::cout << "Total tests: " << total_tests << "\n";
  std::cout << "Passed: " << passed_tests << "\n";
  std::cout << "Failed: " << (total_tests - passed_tests) << "\n";
  std::cout << "Pass rate: " << std::fixed << std::setprecision(1) << pass_rate
            << "%\n";

  if(pass_rate >= 80.0)
  {
    std::cout << "Overall assessment: GOOD random number quality\n";
  }
  else if(pass_rate >= 60.0)
  {
    std::cout << "Overall assessment: ACCEPTABLE random number quality\n";
  }
  else
  {
    std::cout
        << "Overall assessment: POOR random number quality - investigate!\n";
  }

  std::cout << "===============================\n";
}