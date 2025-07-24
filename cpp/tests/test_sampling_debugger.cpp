#include "MonteCarloSampling.hpp"
#include "RandomNumberGenerator.hpp"
#include "SamplingDebugger.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

// Expected PDF for exponential distribution
double exponentialPDF(double x, double lambda = 1.0)
{
  return lambda * std::exp(-lambda * x);
}

// Expected CDF for exponential distribution
double exponentialCDF(double x, double lambda = 1.0)
{
  return 1.0 - std::exp(-lambda * x);
}

// Expected PDF for normal distribution (for testing isotropic direction
// components)
double normalPDF(double x, double mean = 0.0, double stddev = 1.0)
{
  double exponent = -0.5 * std::pow((x - mean) / stddev, 2);
  return (1.0 / (stddev * std::sqrt(2 * M_PI))) * std::exp(exponent);
}

void testExponentialSampling()
{
  std::cout << "\n=== Testing Exponential Sampling ===\n";

  RandomNumberGenerator rng(12345);
  MonteCarloSampling sampler(rng);
  SamplingDebugger debugger("debug_output");

  // Generate samples
  const size_t num_samples = 10000;
  const double lambda = 2.0;

  for(size_t i = 0; i < num_samples; ++i)
  {
    auto sample = sampler.tryExponentialSampling(lambda);
    if(sample)
    {
      debugger.recordSample("exponential_lambda2", sample.value());
    }
  }

  // Generate report
  debugger.generateReport("exponential_lambda2");

  // Perform chi-square test
  auto chi_result = debugger.performChiSquareTest(
      "exponential_lambda2",
      [lambda](double x) { return exponentialPDF(x, lambda); }, 20);

  std::cout << "\nChi-square test results:\n";
  std::cout << "  Chi-square statistic: " << chi_result.chi_square_statistic
            << "\n";
  std::cout << "  Degrees of freedom: " << chi_result.degrees_of_freedom
            << "\n";
  std::cout << "  p-value: " << chi_result.p_value << "\n";
  std::cout << "  Test passed: " << (chi_result.passes_test ? "YES" : "NO")
            << "\n";

  // Perform KS test
  auto [ks_statistic, ks_passed] =
      debugger.performKSTest("exponential_lambda2", [lambda](double x) {
        return exponentialCDF(x, lambda);
      });

  std::cout << "\nKolmogorov-Smirnov test results:\n";
  std::cout << "  KS statistic: " << ks_statistic << "\n";
  std::cout << "  Test passed: " << (ks_passed ? "YES" : "NO") << "\n";

  // Export data for external plotting
  debugger.exportHistogramToCSV("exponential_lambda2", 30,
                                "exponential_histogram.csv");
  debugger.exportSamplesToFile("exponential_lambda2",
                               "exponential_samples.csv");
}

void testIsotropicDirectionSampling()
{
  std::cout << "\n=== Testing Isotropic Direction Sampling ===\n";

  RandomNumberGenerator rng(54321);
  MonteCarloSampling sampler(rng);
  SamplingDebugger debugger("debug_output");

  // Generate samples and record x, y, z components separately
  const size_t num_samples = 10000;

  for(size_t i = 0; i < num_samples; ++i)
  {
    Vector3D direction = sampler.sampleIsotropicDirection();

    // Use getters to access components
    debugger.recordSample("isotropic_x", direction.x());
    debugger.recordSample("isotropic_y", direction.y());
    debugger.recordSample("isotropic_z", direction.z());

    // Also record the magnitude (should be ~1.0)
    double magnitude = direction.magnitude();
    debugger.recordSample("isotropic_magnitude", magnitude);
  }

  // Check statistics for each component
  std::cout << "\nComponent statistics (should have mean~0, variance~0.33):\n";
  for(const std::string &component :
      {"isotropic_x", "isotropic_y", "isotropic_z"})
  {
    auto stats = debugger.calculateStatistics(component);
    std::cout << component << ": mean=" << stats.mean
              << ", variance=" << stats.variance << "\n";
  }

  // Check magnitude statistics (should be very close to 1.0)
  auto mag_stats = debugger.calculateStatistics("isotropic_magnitude");
  std::cout << "\nMagnitude statistics (should be ~1.0):\n";
  std::cout << "  Mean: " << mag_stats.mean << "\n";
  std::cout << "  Min: " << mag_stats.min_value << "\n";
  std::cout << "  Max: " << mag_stats.max_value << "\n";

  // Export histograms
  debugger.exportHistogramToCSV("isotropic_x", 30, "isotropic_x_histogram.csv");
  debugger.exportHistogramToCSV("isotropic_z", 30, "isotropic_z_histogram.csv");
}

void testInteractionDistanceSampling()
{
  std::cout << "\n=== Testing Interaction Distance Sampling ===\n";

  RandomNumberGenerator rng(98765);
  MonteCarloSampling sampler(rng);
  SamplingDebugger debugger("debug_output");

  // Test with different attenuation coefficients
  const std::vector<double> mu_values = {0.5, 1.0, 2.0, 5.0}; // cm^-1

  for(double mu : mu_values)
  {
    std::string dist_name = "interaction_dist_mu" + std::to_string(mu);

    // Generate samples
    for(size_t i = 0; i < 5000; ++i)
    {
      double distance = sampler.sampleInteractionDistance(mu);
      debugger.recordSample(dist_name, distance);
    }

    // Check statistics (mean should be 1/mu)
    auto stats = debugger.calculateStatistics(dist_name);
    double expected_mean = 1.0 / mu;
    std::cout << "\nμ = " << mu << " cm^-1:\n";
    std::cout << "  Expected mean: " << expected_mean << " cm\n";
    std::cout << "  Actual mean: " << stats.mean << " cm\n";
    std::cout << "  Relative error: "
              << std::abs(stats.mean - expected_mean) / expected_mean * 100.0
              << "%\n";

    // Export histogram
    debugger.exportHistogramToCSV(
        dist_name, 25, "interaction_dist_mu" + std::to_string(mu) + ".csv");
  }
}

int main()
{
  std::cout << "=== Monte Carlo Sampling Debugger Tests ===\n";

  testExponentialSampling();
  testIsotropicDirectionSampling();
  testInteractionDistanceSampling();

  std::cout << "\nAll debugging tests completed!\n";
  std::cout << "Check the 'debug_output' directory for detailed reports and "
               "CSV files.\n";

  return 0;
}