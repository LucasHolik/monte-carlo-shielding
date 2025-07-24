#pragma once

#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <optional>
#include <string>
#include <vector>

class SamplingDebugger
{
public:
  struct HistogramBin
  {
    double lower_edge;
    double upper_edge;
    size_t count;
    double expected_frequency;
  };

  struct DistributionStats
  {
    double mean;
    double variance;
    double skewness;
    double kurtosis;
    double min_value;
    double max_value;
    size_t sample_count;
  };

  struct ChiSquareResult
  {
    double chi_square_statistic;
    double p_value;
    size_t degrees_of_freedom;
    bool passes_test; // At 95% confidence level
  };

  // Constructor
  SamplingDebugger(const std::string &output_directory = "debug_output");

  // Record samples for analysis
  void recordSample(const std::string &distribution_name, double value);
  void recordSamples(const std::string &distribution_name,
                     const std::vector<double> &values);

  // Generate histogram for recorded samples
  std::vector<HistogramBin>
  generateHistogram(const std::string &distribution_name, size_t num_bins,
                    std::optional<double> min_value = std::nullopt,
                    std::optional<double> max_value = std::nullopt) const;

  // Calculate distribution statistics
  DistributionStats
  calculateStatistics(const std::string &distribution_name) const;

  // Perform chi-square goodness-of-fit test
  ChiSquareResult
  performChiSquareTest(const std::string &distribution_name,
                       std::function<double(double)> expected_pdf,
                       size_t num_bins = 20) const;

  // Kolmogorov-Smirnov test against expected CDF
  std::pair<double, bool>
  performKSTest(const std::string &distribution_name,
                std::function<double(double)> expected_cdf) const;

  // Export histogram to CSV for external plotting
  void exportHistogramToCSV(const std::string &distribution_name,
                            size_t num_bins, const std::string &filename) const;

  // Export all samples for a distribution
  void exportSamplesToFile(const std::string &distribution_name,
                           const std::string &filename) const;

  // Generate comprehensive report
  void generateReport(const std::string &distribution_name) const;

  // Clear recorded samples for a distribution
  void clearSamples(const std::string &distribution_name);

  // Get sample count
  size_t getSampleCount(const std::string &distribution_name) const;

private:
  std::map<std::string, std::vector<double>> samples_;
  std::string output_directory_;

  // Helper functions
  double calculateChiSquarePValue(double chi_square,
                                  size_t degrees_of_freedom) const;
  double calculateKSCriticalValue(size_t n, double alpha = 0.05) const;
};