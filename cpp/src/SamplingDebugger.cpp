#include "SamplingDebugger.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>

namespace fs = std::filesystem;

SamplingDebugger::SamplingDebugger(const std::string &output_directory)
    : output_directory_(output_directory)
{
  // Directory creation moved to ensureOutputDirectory() for lazy initialization
}

void SamplingDebugger::recordSample(const std::string &distribution_name,
                                    double value)
{
  samples_[distribution_name].push_back(value);
}

void SamplingDebugger::recordSamples(const std::string &distribution_name,
                                     const std::vector<double> &values)
{
  auto &dist_samples = samples_[distribution_name];
  dist_samples.insert(dist_samples.end(), values.begin(), values.end());
}

std::vector<SamplingDebugger::HistogramBin> SamplingDebugger::generateHistogram(
    const std::string &distribution_name, size_t num_bins,
    std::optional<double> min_value, std::optional<double> max_value) const
{

  auto it = samples_.find(distribution_name);
  if(it == samples_.end() || it->second.empty())
  {
    return {};
  }

  const auto &data = it->second;

  // Determine range
  double min_val =
      min_value.value_or(*std::min_element(data.begin(), data.end()));
  double max_val =
      max_value.value_or(*std::max_element(data.begin(), data.end()));

  // Add small margin to ensure all values fit
  double range = max_val - min_val;
  min_val -= range * 0.001;
  max_val += range * 0.001;

  double bin_width = (max_val - min_val) / num_bins;

  std::vector<HistogramBin> histogram(num_bins);

  // Initialize bins
  for(size_t i = 0; i < num_bins; ++i)
  {
    histogram[i].lower_edge = min_val + i * bin_width;
    histogram[i].upper_edge = min_val + (i + 1) * bin_width;
    histogram[i].count = 0;
    histogram[i].expected_frequency = 0.0;
  }

  // Count samples in each bin
  for(double value : data)
  {
    size_t bin_index = static_cast<size_t>((value - min_val) / bin_width);
    if(bin_index >= num_bins)
      bin_index = num_bins - 1;
    histogram[bin_index].count++;
  }

  return histogram;
}

SamplingDebugger::DistributionStats SamplingDebugger::calculateStatistics(
    const std::string &distribution_name) const
{

  DistributionStats stats{};

  auto it = samples_.find(distribution_name);
  if(it == samples_.end() || it->second.empty())
  {
    return stats;
  }

  const auto &data = it->second;
  stats.sample_count = data.size();

  // Basic statistics
  stats.min_value = *std::min_element(data.begin(), data.end());
  stats.max_value = *std::max_element(data.begin(), data.end());

  // Mean
  stats.mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();

  // Variance, skewness, kurtosis
  double sum_sq_diff = 0.0;
  double sum_cubed_diff = 0.0;
  double sum_fourth_diff = 0.0;

  for(double value : data)
  {
    double diff = value - stats.mean;
    double diff_sq = diff * diff;
    sum_sq_diff += diff_sq;
    sum_cubed_diff += diff_sq * diff;
    sum_fourth_diff += diff_sq * diff_sq;
  }

  stats.variance = sum_sq_diff / (data.size() - 1);
  double std_dev = std::sqrt(stats.variance);

  if(std_dev > 0)
  {
    stats.skewness = (sum_cubed_diff / data.size()) / std::pow(std_dev, 3);
    stats.kurtosis =
        (sum_fourth_diff / data.size()) / std::pow(std_dev, 4) - 3.0;
  }

  return stats;
}

SamplingDebugger::ChiSquareResult SamplingDebugger::performChiSquareTest(
    const std::string &distribution_name,
    std::function<double(double)> expected_pdf, size_t num_bins) const
{

  ChiSquareResult result{};

  auto histogram = generateHistogram(distribution_name, num_bins);
  if(histogram.empty())
    return result;

  auto it = samples_.find(distribution_name);
  if(it == samples_.end())
    return result;

  size_t total_samples = it->second.size();
  double chi_square = 0.0;
  size_t valid_bins = 0;

  for(auto &bin : histogram)
  {
    double bin_center = (bin.lower_edge + bin.upper_edge) / 2.0;
    double bin_width = bin.upper_edge - bin.lower_edge;

    // Expected frequency = pdf * bin_width * total_samples
    double expected_freq = expected_pdf(bin_center) * bin_width * total_samples;

    // Chi-square test requires expected frequency >= 5
    if(expected_freq >= 5.0)
    {
      double observed = static_cast<double>(bin.count);
      chi_square += std::pow(observed - expected_freq, 2) / expected_freq;
      valid_bins++;
    }
  }

  result.chi_square_statistic = chi_square;
  result.degrees_of_freedom = valid_bins - 1;
  result.p_value =
      calculateChiSquarePValue(chi_square, result.degrees_of_freedom);
  result.passes_test = result.p_value > 0.05;

  return result;
}

std::pair<double, bool> SamplingDebugger::performKSTest(
    const std::string &distribution_name,
    std::function<double(double)> expected_cdf) const
{

  auto it = samples_.find(distribution_name);
  if(it == samples_.end() || it->second.empty())
  {
    return {0.0, false};
  }

  std::vector<double> sorted_data = it->second;
  std::sort(sorted_data.begin(), sorted_data.end());

  double max_distance = 0.0;
  size_t n = sorted_data.size();

  for(size_t i = 0; i < n; ++i)
  {
    double empirical_cdf = static_cast<double>(i + 1) / n;
    double expected = expected_cdf(sorted_data[i]);
    double distance = std::abs(empirical_cdf - expected);
    max_distance = std::max(max_distance, distance);

    // Also check the distance just before this point
    double empirical_cdf_prev = static_cast<double>(i) / n;
    distance = std::abs(empirical_cdf_prev - expected);
    max_distance = std::max(max_distance, distance);
  }

  double critical_value = calculateKSCriticalValue(n);
  bool passes_test = max_distance < critical_value;

  return {max_distance, passes_test};
}

void SamplingDebugger::exportHistogramToCSV(
    const std::string &distribution_name, size_t num_bins,
    const std::string &filename) const
{

  auto histogram = generateHistogram(distribution_name, num_bins);
  if(histogram.empty())
    return;

  ensureOutputDirectory();
  fs::path filepath = fs::path(output_directory_) / filename;
  std::ofstream file(filepath);

  if(!file.is_open())
  {
    std::cerr << "Failed to open file: " << filepath << std::endl;
    return;
  }

  file << "bin_center,bin_lower,bin_upper,count,frequency\n";

  auto it = samples_.find(distribution_name);
  size_t total_samples = it->second.size();

  for(const auto &bin : histogram)
  {
    double bin_center = (bin.lower_edge + bin.upper_edge) / 2.0;
    double frequency = static_cast<double>(bin.count) / total_samples;

    file << std::fixed << std::setprecision(6) << bin_center << ","
         << bin.lower_edge << "," << bin.upper_edge << "," << bin.count << ","
         << frequency << "\n";
  }

  file.close();
}

void SamplingDebugger::exportSamplesToFile(const std::string &distribution_name,
                                           const std::string &filename) const
{

  auto it = samples_.find(distribution_name);
  if(it == samples_.end())
    return;

  ensureOutputDirectory();
  fs::path filepath = fs::path(output_directory_) / filename;
  std::ofstream file(filepath);

  if(!file.is_open())
  {
    std::cerr << "Failed to open file: " << filepath << std::endl;
    return;
  }

  file << "sample_index,value\n";

  size_t index = 0;
  for(double value : it->second)
  {
    file << index++ << "," << std::fixed << std::setprecision(10) << value
         << "\n";
  }

  file.close();
}

void SamplingDebugger::generateReport(
    const std::string &distribution_name) const
{
  auto stats = calculateStatistics(distribution_name);

  std::stringstream report;
  report << "=== Sampling Debug Report for: " << distribution_name
         << " ===\n\n";

  report << "Sample Statistics:\n";
  report << "  Total Samples: " << stats.sample_count << "\n";
  report << "  Mean: " << std::fixed << std::setprecision(6) << stats.mean
         << "\n";
  report << "  Variance: " << stats.variance << "\n";
  report << "  Std Dev: " << std::sqrt(stats.variance) << "\n";
  report << "  Skewness: " << stats.skewness << "\n";
  report << "  Kurtosis: " << stats.kurtosis << "\n";
  report << "  Min: " << stats.min_value << "\n";
  report << "  Max: " << stats.max_value << "\n\n";

  // Generate histogram
  auto histogram = generateHistogram(distribution_name, 20);
  report << "Histogram (20 bins):\n";

  // Find max count for scaling
  size_t max_count = 0;
  for(const auto &bin : histogram)
  {
    max_count = std::max(max_count, bin.count);
  }

  // ASCII histogram
  for(const auto &bin : histogram)
  {
    double bin_center = (bin.lower_edge + bin.upper_edge) / 2.0;
    report << std::setw(10) << std::fixed << std::setprecision(3) << bin_center
           << " | ";

    size_t bar_length = (bin.count * 50) / max_count;
    for(size_t i = 0; i < bar_length; ++i)
    {
      report << "#";
    }
    report << " (" << bin.count << ")\n";
  }

  // Save report
  ensureOutputDirectory();
  fs::path filepath =
      fs::path(output_directory_) / (distribution_name + "_report.txt");
  std::ofstream file(filepath);
  if(file.is_open())
  {
    file << report.str();
    file.close();
  }

  // Also print to console
  std::cout << report.str();
}

void SamplingDebugger::clearSamples(const std::string &distribution_name)
{
  samples_.erase(distribution_name);
}

size_t
SamplingDebugger::getSampleCount(const std::string &distribution_name) const
{
  auto it = samples_.find(distribution_name);
  return (it != samples_.end()) ? it->second.size() : 0;
}

double
SamplingDebugger::calculateChiSquarePValue(double chi_square,
                                           size_t degrees_of_freedom) const
{
  // Simplified p-value calculation using chi-square table critical values
  // For more accuracy, use a proper statistical library

  // Critical values at 95% confidence level
  std::vector<double> critical_values_95 = {3.841,  5.991,  7.815,  9.488,
                                            11.070, 12.592, 14.067, 15.507,
                                            16.919, 18.307};

  if(degrees_of_freedom > 0 && degrees_of_freedom <= 10)
  {
    double critical_value = critical_values_95[degrees_of_freedom - 1];
    return (chi_square < critical_value) ? 0.1 : 0.01; // Simplified
  }

  // For higher degrees of freedom, use approximation
  double z = std::sqrt(2 * chi_square) - std::sqrt(2 * degrees_of_freedom - 1);
  return 0.5 * std::erfc(z / std::sqrt(2));
}

double SamplingDebugger::calculateKSCriticalValue(size_t n, double alpha) const
{
  // Kolmogorov-Smirnov critical value approximation
  if(alpha == 0.05)
  {
    return 1.36 / std::sqrt(static_cast<double>(n));
  }
  else if(alpha == 0.01)
  {
    return 1.63 / std::sqrt(static_cast<double>(n));
  }
  return 1.36 / std::sqrt(static_cast<double>(n)); // Default to 0.05
}

void SamplingDebugger::ensureOutputDirectory() const
{
  try
  {
    // Check if directory already exists to avoid unnecessary filesystem calls
    if(!fs::exists(output_directory_))
    {
      fs::create_directories(output_directory_);
    }
  }
  catch(const std::filesystem::filesystem_error &e)
  {
    // Log error but continue execution - files will fail gracefully
    std::cerr << "Warning: Could not create output directory '" << output_directory_
              << "': " << e.what() << std::endl;
  }
}