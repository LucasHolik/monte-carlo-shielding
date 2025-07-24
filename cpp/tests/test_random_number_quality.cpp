#include "RandomNumberQuality.hpp"
#include <cassert>
#include <iomanip>
#include <iostream>

void testBasicQualityChecks()
{
  std::cout << "Testing basic quality checks..." << std::endl;

  RandomNumberGenerator rng(12345);
  RandomNumberQuality tester(rng);

  // Test chi-square
  auto chi_result = tester.chiSquareTest(10000, 50);
  std::cout << "  Chi-square test: "
            << (chi_result.passed ? "PASSED" : "FAILED")
            << " (statistic=" << chi_result.test_statistic
            << ", p-value=" << chi_result.p_value << ")" << std::endl;
  assert(chi_result.passed);

  // Test KS
  auto ks_result = tester.kolmogorovSmirnovTest(1000);
  std::cout << "  KS test: " << (ks_result.passed ? "PASSED" : "FAILED")
            << " (statistic=" << ks_result.test_statistic
            << ", p-value=" << ks_result.p_value << ")" << std::endl;
  assert(ks_result.passed);

  // Test correlation
  auto corr_result = tester.serialCorrelationTest(5000, 1);
  std::cout << "  Serial correlation test: "
            << (corr_result.passed ? "PASSED" : "FAILED")
            << " (correlation=" << corr_result.test_statistic
            << ", p-value=" << corr_result.p_value << ")" << std::endl;
  assert(corr_result.passed);

  std::cout << "✓ Basic quality checks passed" << std::endl;
}

void testDistributionQuality()
{
  std::cout << "Testing distribution quality..." << std::endl;

  RandomNumberGenerator rng(23456);
  RandomNumberQuality tester(rng);

  // Test exponential
  auto exp_result = tester.exponentialDistributionTest(2.0, 10000);
  std::cout << "  Exponential test: "
            << (exp_result.passed ? "PASSED" : "FAILED") << std::endl;
  std::cout << "    Sample mean: " << exp_result.details["sample_mean"]
            << " (expected: " << exp_result.details["expected_mean"] << ")"
            << std::endl;
  assert(exp_result.passed);

  // Test normal
  auto norm_result = tester.normalDistributionTest(5.0, 2.0, 10000);
  std::cout << "  Normal test: " << (norm_result.passed ? "PASSED" : "FAILED")
            << std::endl;
  std::cout << "    Sample mean: " << norm_result.details["sample_mean"]
            << ", Sample stddev: " << norm_result.details["sample_stddev"]
            << std::endl;
  std::cout << "    Skewness: " << norm_result.details["skewness"]
            << ", Kurtosis: " << norm_result.details["excess_kurtosis"]
            << std::endl;
  assert(norm_result.passed);

  // Test isotropic directions
  auto iso_result = tester.isotropicDirectionTest(5000);
  std::cout << "  Isotropic direction test: "
            << (iso_result.passed ? "PASSED" : "FAILED") << std::endl;
  std::cout << "    Average magnitude: "
            << iso_result.details["average_magnitude"] << std::endl;
  assert(iso_result.passed);

  std::cout << "✓ Distribution quality tests passed" << std::endl;
}

void testMonteCarloSpecific()
{
  std::cout << "Testing Monte Carlo specific quality..." << std::endl;

  RandomNumberGenerator rng(34567);
  RandomNumberQuality tester(rng);

  // Test random walk
  auto walk_result = tester.randomWalkTest(1000, 500);
  std::cout << "  Random walk test: "
            << (walk_result.passed ? "PASSED" : "FAILED") << std::endl;
  std::cout << "    Mean R²: " << walk_result.test_statistic
            << " (expected: " << walk_result.details["expected_r_squared"]
            << ")" << std::endl;
  assert(walk_result.passed);

  // Test particle splitting
  auto split_result = tester.particleSplittingTest(2.5, 10000);
  std::cout << "  Particle splitting test: "
            << (split_result.passed ? "PASSED" : "FAILED") << std::endl;
  assert(split_result.passed);

  // Test Russian roulette
  auto rr_result = tester.russianRouletteTest(0.8, 10000);
  std::cout << "  Russian roulette test: "
            << (rr_result.passed ? "PASSED" : "FAILED") << std::endl;
  std::cout << "    Actual survival rate: " << rr_result.test_statistic
            << std::endl;
  assert(rr_result.passed);

  std::cout << "✓ Monte Carlo specific tests passed" << std::endl;
}

void testVisualTools()
{
  std::cout << "Testing visual/debug tools..." << std::endl;

  RandomNumberGenerator rng(45678);
  RandomNumberQuality tester(rng);

  // Test histogram generation
  auto [bin_edges, counts] = tester.generateUniformityHistogram(10000, 20);
  assert(bin_edges.size() == 21); // n_bins + 1
  assert(counts.size() == 20);

  // Check total count
  int total = std::accumulate(counts.begin(), counts.end(), 0);
  assert(total == 10000);

  // Test 2D scatter plot
  auto scatter = tester.generate2DScatterPlot(1000);
  assert(scatter.size() == 1000);

  // Check all points in [0,1) x [0,1)
  for(const auto &[x, y] : scatter)
  {
    assert(x >= 0.0 && x < 1.0);
    assert(y >= 0.0 && y < 1.0);
  }

  // Test correlation plot
  auto correlations = tester.generateCorrelationPlot(10, 5000);
  assert(correlations.size() == 10);

  // Correlations should be small for good RNG
  for(const auto &[lag, corr] : correlations)
  {
    assert(std::abs(corr) < 0.1); // Should be weakly correlated
  }

  std::cout << "✓ Visual tools tests passed" << std::endl;
}

void testComprehensiveTestSuite()
{
  std::cout << "Running comprehensive test suite..." << std::endl;

  RandomNumberGenerator rng(56789);
  RandomNumberQuality tester(rng, true); // Verbose mode

  auto results = tester.runCompleteTestSuite();

  // Check we have all categories
  assert(results.count("Uniformity") > 0);
  assert(results.count("Independence") > 0);
  assert(results.count("Monte Carlo") > 0);

  // Print summary
  RandomNumberQuality::printTestReport(results);

  // Export results
  tester.exportResults(results, "rng_quality_report.txt");
  std::cout << "  Results exported to rng_quality_report.txt" << std::endl;

  // Check overall quality
  int total = 0, passed = 0;
  for(const auto &[category, tests] : results)
  {
    for(const auto &test : tests)
    {
      total++;
      if(test.passed)
        passed++;
    }
  }

  double pass_rate = static_cast<double>(passed) / total;
  std::cout << "  Overall pass rate: " << std::fixed << std::setprecision(1)
            << (pass_rate * 100) << "%" << std::endl;

  assert(pass_rate >= 0.90); // Should pass at least 90% of tests

  std::cout << "✓ Comprehensive test suite passed" << std::endl;
}

void testStandaloneFunctions()
{
  std::cout << "Testing standalone quality functions..." << std::endl;

  RandomNumberGenerator rng1(67890);
  RandomNumberGenerator rng2(98765);

  // Quick quality check
  bool quick_pass = RandomQualityTests::quickQualityCheck(rng1);
  std::cout << "  Quick quality check: " << (quick_pass ? "PASSED" : "FAILED")
            << std::endl;
  assert(quick_pass);

  // Compare generators
  auto comparison = RandomQualityTests::compareGenerators(rng1, rng2, 10000);
  std::cout << "  Generator comparison:" << std::endl;
  std::cout << "    Speed ratio: " << comparison["speed_ratio"] << std::endl;
  std::cout << "    Mean difference: " << comparison["mean_difference"]
            << std::endl;
  std::cout << "    KS statistic: " << comparison["ks_statistic"] << std::endl;

  // Test performance
  double rate = RandomQualityTests::testGeneratorPerformance(rng1, 1000000);
  std::cout << "  Performance: " << std::scientific << std::setprecision(2)
            << rate << " samples/second" << std::endl;
  assert(rate > 1e6); // Should generate at least 1M samples/second

  std::cout << "✓ Standalone functions passed" << std::endl;
}

void demonstrateQualityIssues()
{
  std::cout << "\nDemonstrating quality issue detection..." << std::endl;

  // Create a "bad" RNG that always returns 0.5
  class BadRNG : public RandomNumberGenerator
  {
  public:
    BadRNG() : RandomNumberGenerator(1) {}
    double uniform() const override { return 0.5; }
  };

  BadRNG bad_rng;
  RandomNumberQuality tester(bad_rng);

  // This should fail uniformity tests
  auto chi_result = tester.chiSquareTest(1000, 10);
  std::cout << "  Bad RNG chi-square test: "
            << (chi_result.passed ? "PASSED" : "FAILED") << " (as expected)"
            << std::endl;
  assert(!chi_result.passed); // Should fail!

  // This should fail independence tests
  auto corr_result = tester.serialCorrelationTest(1000, 1);
  std::cout << "  Bad RNG correlation test: correlation = "
            << corr_result.test_statistic << std::endl;

  std::cout << "✓ Quality issue detection working correctly" << std::endl;
}

int main()
{
  std::cout << "Running Random Number Quality Tests...\n" << std::endl;

  testBasicQualityChecks();
  testDistributionQuality();
  testMonteCarloSpecific();
  testVisualTools();
  testComprehensiveTestSuite();
  testStandaloneFunctions();
  demonstrateQualityIssues();

  std::cout << "\n🎉 All Random Number Quality tests passed!" << std::endl;
  std::cout << "✓ Statistical validation for sampling distributions complete!"
            << std::endl;
  std::cout << "✓ Debugging tools for sampling verification implemented!"
            << std::endl;
  std::cout << "✓ Comprehensive quality testing framework ready!" << std::endl;
  std::cout
      << "\nDay 4 Task Complete: Random number quality testing implemented! 🚀"
      << std::endl;

  return 0;
}