#pragma once

#include "RandomNumberGenerator.hpp"
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <vector>

/**
 * @brief Result of a statistical test
 */
struct TestResult
{
  bool passed;
  double test_statistic;
  double p_value;
  std::string description;
  std::map<std::string, double> details;

  TestResult(bool p, double stat, double pval, const std::string &desc)
      : passed(p), test_statistic(stat), p_value(pval), description(desc)
  {}
};

/**
 * @brief Comprehensive random number quality testing
 *
 * Implements various statistical tests to verify the quality of random
 * number generators for Monte Carlo simulations.
 */
class RandomNumberQuality
{
private:
  RandomNumberGenerator &rng_;
  bool verbose_;

public:
  explicit RandomNumberQuality(RandomNumberGenerator &rng,
                               bool verbose = false);

  // ============================================================================
  // UNIFORMITY TESTS
  // ============================================================================

  /**
   * @brief Chi-square test for uniformity
   * @param n_samples Number of samples to generate
   * @param n_bins Number of bins to use
   * @param confidence_level Confidence level (e.g., 0.95 for 95%)
   * @return Test result with chi-square statistic and p-value
   */
  TestResult chiSquareTest(std::size_t n_samples = 10000,
                           std::size_t n_bins = 100,
                           double confidence_level = 0.95);

  /**
   * @brief Kolmogorov-Smirnov test for uniform distribution
   * @param n_samples Number of samples to test
   * @param confidence_level Confidence level
   * @return Test result with KS statistic and p-value
   */
  TestResult kolmogorovSmirnovTest(std::size_t n_samples = 1000,
                                   double confidence_level = 0.95);

  /**
   * @brief Anderson-Darling test for uniformity
   * @param n_samples Number of samples
   * @return Test result
   */
  TestResult andersonDarlingTest(std::size_t n_samples = 1000);

  // ============================================================================
  // INDEPENDENCE TESTS
  // ============================================================================

  /**
   * @brief Serial correlation test (autocorrelation)
   * @param n_samples Number of samples
   * @param lag Lag for correlation (default 1)
   * @return Test result with correlation coefficient
   */
  TestResult serialCorrelationTest(std::size_t n_samples = 10000,
                                   std::size_t lag = 1);

  /**
   * @brief Run test for independence
   * @param n_samples Number of samples
   * @return Test result with number of runs
   */
  TestResult runTest(std::size_t n_samples = 10000);

  /**
   * @brief Gap test for detecting patterns
   * @param n_samples Number of samples
   * @param alpha Lower bound of gap interval
   * @param beta Upper bound of gap interval
   * @return Test result
   */
  TestResult gapTest(std::size_t n_samples = 10000, double alpha = 0.0,
                     double beta = 0.5);

  // ============================================================================
  // DISTRIBUTION-SPECIFIC TESTS
  // ============================================================================

  /**
   * @brief Test exponential distribution quality
   * @param lambda Rate parameter
   * @param n_samples Number of samples
   * @return Test result
   */
  TestResult exponentialDistributionTest(double lambda,
                                         std::size_t n_samples = 10000);

  /**
   * @brief Test normal distribution quality
   * @param mean Expected mean
   * @param stddev Expected standard deviation
   * @param n_samples Number of samples
   * @return Test result
   */
  TestResult normalDistributionTest(double mean, double stddev,
                                    std::size_t n_samples = 10000);

  /**
   * @brief Test isotropic direction generation
   * @param n_samples Number of direction vectors
   * @return Test result
   */
  TestResult isotropicDirectionTest(std::size_t n_samples = 10000);

  // ============================================================================
  // MONTE CARLO SPECIFIC TESTS
  // ============================================================================

  /**
   * @brief Test random walk behaviour
   * @param n_steps Number of steps per walk
   * @param n_walks Number of random walks
   * @return Test result comparing to theoretical diffusion
   */
  TestResult randomWalkTest(std::size_t n_steps = 1000,
                            std::size_t n_walks = 1000);

  /**
   * @brief Test particle splitting consistency
   * @param splitting_factor Average splitting factor
   * @param n_trials Number of trials
   * @return Test result
   */
  TestResult particleSplittingTest(double splitting_factor = 2.5,
                                   std::size_t n_trials = 10000);

  /**
   * @brief Test Russian roulette fairness
   * @param survival_prob Survival probability
   * @param n_trials Number of trials
   * @return Test result
   */
  TestResult russianRouletteTest(double survival_prob = 0.8,
                                 std::size_t n_trials = 10000);

  // ============================================================================
  // COMPREHENSIVE TEST SUITES
  // ============================================================================

  /**
   * @brief Run all uniformity tests
   * @return Vector of test results
   */
  std::vector<TestResult> runUniformityTests(std::size_t n_samples = 10000);

  /**
   * @brief Run all independence tests
   * @return Vector of test results
   */
  std::vector<TestResult> runIndependenceTests(std::size_t n_samples = 10000);

  /**
   * @brief Run all Monte Carlo specific tests
   * @return Vector of test results
   */
  std::vector<TestResult> runMonteCarloTests();

  /**
   * @brief Run complete test suite
   * @return Map of test category to results
   */
  std::map<std::string, std::vector<TestResult>> runCompleteTestSuite();

  // ============================================================================
  // VISUAL AND DEBUG TOOLS
  // ============================================================================

  /**
   * @brief Generate data for uniformity histogram
   * @param n_samples Number of samples
   * @param n_bins Number of bins
   * @return Histogram data (bin edges and counts)
   */
  std::pair<std::vector<double>, std::vector<int>>
  generateUniformityHistogram(std::size_t n_samples = 10000,
                              std::size_t n_bins = 50);

  /**
   * @brief Generate scatter plot data for 2D uniformity
   * @param n_points Number of points
   * @return Pairs of (x, y) coordinates
   */
  std::vector<std::pair<double, double>>
  generate2DScatterPlot(std::size_t n_points = 5000);

  /**
   * @brief Generate serial correlation plot data
   * @param max_lag Maximum lag to test
   * @param n_samples Number of samples
   * @return Map of lag to correlation coefficient
   */
  std::map<std::size_t, double>
  generateCorrelationPlot(std::size_t max_lag = 50,
                          std::size_t n_samples = 10000);

  /**
   * @brief Export test results to file
   * @param results Test results to export
   * @param filename Output filename
   */
  void
  exportResults(const std::map<std::string, std::vector<TestResult>> &results,
                const std::string &filename);

  // ============================================================================
  // UTILITY METHODS
  // ============================================================================

  /**
   * @brief Set verbosity for detailed output
   */
  void setVerbose(bool verbose) { verbose_ = verbose; }

  /**
   * @brief Get string summary of test result
   */
  static std::string summarizeResult(const TestResult &result);

  /**
   * @brief Print detailed test report
   */
  static void printTestReport(
      const std::map<std::string, std::vector<TestResult>> &results);

private:
  // Statistical utility functions
  double calculateChiSquareStatistic(const std::vector<int> &observed,
                                     const std::vector<double> &expected);

  double calculateKSStatistic(const std::vector<double> &samples);

  double calculateAutoCorrelation(const std::vector<double> &data,
                                  std::size_t lag);

  std::size_t countRuns(const std::vector<double> &data,
                        double threshold = 0.5);

  // Critical value lookup (simplified - would use proper statistical library in
  // production)
  double getChiSquareCriticalValue(std::size_t degrees_of_freedom,
                                   double alpha);
  double getKSCriticalValue(std::size_t n_samples, double alpha);
  double getNormalCriticalValue(double alpha);
};

// ============================================================================
// STANDALONE QUALITY TESTING FUNCTIONS
// ============================================================================

namespace RandomQualityTests
{
/**
 * @brief Quick quality check for a random number generator
 * @param rng Random number generator to test
 * @return True if all basic tests pass
 */
bool quickQualityCheck(RandomNumberGenerator &rng);

/**
 * @brief Detailed quality assessment with report
 * @param rng Random number generator to test
 * @param output_file Optional file for detailed report
 * @return Overall pass/fail status
 */
bool detailedQualityAssessment(RandomNumberGenerator &rng,
                               const std::string &output_file = "");

/**
 * @brief Compare two random number generators
 * @param rng1 First generator
 * @param rng2 Second generator
 * @param n_samples Number of samples for comparison
 * @return Statistical comparison results
 */
std::map<std::string, double> compareGenerators(RandomNumberGenerator &rng1,
                                                RandomNumberGenerator &rng2,
                                                std::size_t n_samples = 10000);

/**
 * @brief Test generator performance
 * @param rng Generator to test
 * @param n_samples Number of samples to generate
 * @return Generation rate (samples per second)
 */
double testGeneratorPerformance(RandomNumberGenerator &rng,
                                std::size_t n_samples = 1000000);
} // namespace RandomQualityTests