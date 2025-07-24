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
  std::cout << "    Sample mean: " << exp_result.details.at("sample_mean")
            << " (expected: " << exp_result.details.at("expected_mean") << ")"
            << std::endl;
  assert(exp_result.passed);

  // Test normal
  auto norm_result = tester.normalDistributionTest(5.0, 2.0, 10000);
  std::cout << "  Normal test: " << (norm_result.passed ? "PASSED" : "FAILED")
            << std::endl;
  std::cout << "    Sample mean: " << norm_result.details.at("sample_mean")
            << ", Sample stddev: " << norm_result.details.at("sample_stddev")
            << std::endl;
  std::cout << "    Mean z-score: " << norm_result.details.at("mean_z_score")
            << ", Stddev error: " << norm_result.details.at("stddev_error")
            << std::endl;
  assert(norm_result.passed);

  // Test isotropic directions
  auto iso_result = tester.isotropicDirectionTest(5000);
  std::cout << "  Isotropic direction test: "
            << (iso_result.passed ? "PASSED" : "FAILED") << std::endl;
  std::cout << "    Max length error: "
            << iso_result.details.at("max_length_error")
            << ", Mean magnitude: " << iso_result.details.at("mean_magnitude")
            << std::endl;
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
  std::cout << "    Sample mean distance: "
            << walk_result.details.at("sample_mean_distance") << " (expected: "
            << walk_result.details.at("expected_mean_distance") << ")"
            << std::endl;
  assert(walk_result.passed);

  // Test particle splitting
  auto split_result = tester.particleSplittingTest(2.5, 10000);
  std::cout << "  Particle splitting test: "
            << (split_result.passed ? "PASSED" : "FAILED") << std::endl;
  std::cout << "    Sample mean: " << split_result.details.at("sample_mean")
            << " (expected: " << split_result.details.at("expected_mean") << ")"
            << std::endl;
  assert(split_result.passed);

  // Test Russian roulette
  auto rr_result = tester.russianRouletteTest(0.8, 10000);
  std::cout << "  Russian roulette test: "
            << (rr_result.passed ? "PASSED" : "FAILED") << std::endl;
  std::cout << "    Sample survival rate: "
            << rr_result.details.at("sample_rate")
            << " (expected: " << rr_result.details.at("expected_rate") << ")"
            << std::endl;
  assert(rr_result.passed);

  std::cout << "✓ Monte Carlo specific tests passed" << std::endl;
}

void testIndependenceTests()
{
  std::cout << "Testing independence and randomness..." << std::endl;

  RandomNumberGenerator rng(45678);
  RandomNumberQuality tester(rng);

  // Test Anderson-Darling
  auto ad_result = tester.andersonDarlingTest(1000);
  std::cout << "  Anderson-Darling test: "
            << (ad_result.passed ? "PASSED" : "FAILED")
            << " (statistic=" << ad_result.test_statistic << ")" << std::endl;
  assert(ad_result.passed);

  // Test run test
  auto run_result = tester.runTest(5000);
  std::cout << "  Run test: " << (run_result.passed ? "PASSED" : "FAILED")
            << " (runs=" << run_result.test_statistic << ", expected≈"
            << run_result.details.at("expected_runs") << ")" << std::endl;
  assert(run_result.passed);

  // Test gap test
  auto gap_result = tester.gapTest(10000, 0.0, 0.5);
  std::cout << "  Gap test: " << (gap_result.passed ? "PASSED" : "FAILED")
            << " (test_stat=" << gap_result.test_statistic << ")" << std::endl;
  assert(gap_result.passed);

  // Test correlation with different lags
  auto corr5_result = tester.serialCorrelationTest(5000, 5);
  std::cout << "  Serial correlation (lag=5): "
            << (corr5_result.passed ? "PASSED" : "FAILED")
            << " (correlation=" << corr5_result.test_statistic << ")"
            << std::endl;
  assert(corr5_result.passed);

  std::cout << "✓ Independence tests passed" << std::endl;
}

void testComprehensiveTestSuite()
{
  std::cout << "Running comprehensive test suite..." << std::endl;

  RandomNumberGenerator rng(67890);
  RandomNumberQuality tester(rng, false); // Non-verbose for cleaner output

  auto results = tester.runCompleteTestSuite();

  // Check we have all categories
  assert(results.count("Uniformity") > 0);
  assert(results.count("Independence") > 0);
  assert(results.count("Monte Carlo") > 0);

  // Print summary
  RandomNumberQuality::printTestReport(results);

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

  assert(pass_rate >= 0.75); // Should pass at least 75% of tests

  std::cout << "✓ Comprehensive test suite passed" << std::endl;
}

int main()
{
  std::cout << "Running Random Number Quality Tests...\n" << std::endl;

  try
  {
    testBasicQualityChecks();
    testDistributionQuality();
    testMonteCarloSpecific();
    testIndependenceTests();
    testComprehensiveTestSuite();

    std::cout << "\n🎉 All Random Number Quality tests passed!" << std::endl;
    std::cout << "✓ Statistical validation for sampling distributions complete!"
              << std::endl;
    std::cout << "✓ Comprehensive quality testing framework ready!"
              << std::endl;
    std::cout << "✓ Chi-square, KS, Anderson-Darling tests implemented!"
              << std::endl;
    std::cout << "✓ Independence tests (correlation, runs, gaps) working!"
              << std::endl;
    std::cout << "✓ Distribution-specific validation (exponential, normal, "
                 "isotropic) ready!"
              << std::endl;
    std::cout << "✓ Monte Carlo specific tests (random walk, splitting, "
                 "roulette) implemented!"
              << std::endl;
    std::cout << "\nDay 4 Task Complete: Statistical validation for sampling "
                 "distributions implemented! 🚀"
              << std::endl;
  }
  catch(const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}