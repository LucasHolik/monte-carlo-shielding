#include "../include/ProgressBar.hpp"

#include <iomanip>
#include <sstream>
#include <algorithm>
#include <cmath>

ProgressBar::ProgressBar(size_t total_items, size_t update_frequency, size_t bar_width, bool show_bar)
    : total_items_(total_items)
    , current_item_(0)
    , update_frequency_(update_frequency)
    , is_started_(false)
    , show_bar_(show_bar)
    , bar_width_(bar_width)
{
  // Ensure sensible defaults
  if (total_items_ == 0) {
    total_items_ = 1;
  }
  if (update_frequency_ == 0) {
    update_frequency_ = 1;
  }
  if (bar_width_ < 10) {
    bar_width_ = 10;
  }
}

void ProgressBar::start()
{
  start_time_ = std::chrono::steady_clock::now();
  last_update_ = start_time_;
  current_item_ = 0;
  is_started_ = true;
  
  // Initial display
  display();
}

void ProgressBar::update(bool force_display)
{
  if (!is_started_) {
    start();
  }
  
  current_item_++;
  
  if (force_display || shouldUpdate()) {
    display();
    last_update_ = std::chrono::steady_clock::now();
  }
}

void ProgressBar::updateTo(size_t current_item, bool force_display)
{
  if (!is_started_) {
    start();
  }
  
  current_item_ = std::min(current_item, total_items_);
  
  if (force_display || shouldUpdate()) {
    display();
    last_update_ = std::chrono::steady_clock::now();
  }
}

void ProgressBar::finish()
{
  if (!is_started_) {
    start();
  }
  
  current_item_ = total_items_;
  display();
  std::cout << std::endl;
  
  // Show completion summary
  double total_time = getElapsedSeconds();
  std::cout << "✅ Completed " << total_items_ << " items in " 
            << formatTime(total_time) << std::endl;
}

double ProgressBar::getProgressPercentage() const
{
  if (total_items_ == 0) return 1.0;
  return static_cast<double>(current_item_) / static_cast<double>(total_items_);
}

double ProgressBar::getElapsedSeconds() const
{
  if (!is_started_) return 0.0;
  
  auto now = std::chrono::steady_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time_);
  return elapsed.count() / 1000.0;
}

double ProgressBar::getETASeconds() const
{
  if (!is_started_ || current_item_ == 0) return 0.0;
  
  double elapsed = getElapsedSeconds();
  double rate = static_cast<double>(current_item_) / elapsed;
  double remaining_items = static_cast<double>(total_items_ - current_item_);
  
  if (rate <= 0.0) return 0.0;
  return remaining_items / rate;
}

void ProgressBar::display()
{
  std::ostringstream oss;
  
  // Progress percentage
  double percentage = getProgressPercentage();
  int percent_int = static_cast<int>(percentage * 100);
  
  // Visual progress bar
  if (show_bar_) {
    oss << createProgressBar() << " ";
  }
  
  // Percentage and count
  oss << std::setw(3) << percent_int << "% (" 
      << current_item_ << "/" << total_items_ << ")";
  
  // Timing information
  double elapsed = getElapsedSeconds();
  oss << " | " << formatTime(elapsed) << " elapsed";
  
  if (current_item_ > 0 && current_item_ < total_items_) {
    double eta = getETASeconds();
    oss << " | ETA: " << formatTime(eta);
  }
  
  // Clear line and display
  std::cout << "\r" << std::string(100, ' ') << "\r" << oss.str() << std::flush;
}

std::string ProgressBar::formatTime(double seconds) const
{
  std::ostringstream oss;
  
  if (seconds < 60.0) {
    oss << std::fixed << std::setprecision(1) << seconds << "s";
  } else if (seconds < 3600.0) {
    int minutes = static_cast<int>(seconds / 60);
    int secs = static_cast<int>(seconds) % 60;
    oss << minutes << "m " << secs << "s";
  } else {
    int hours = static_cast<int>(seconds / 3600);
    int minutes = static_cast<int>((seconds - hours * 3600) / 60);
    oss << hours << "h " << minutes << "m";
  }
  
  return oss.str();
}

std::string ProgressBar::createProgressBar() const
{
  std::string bar = "[";
  
  double percentage = getProgressPercentage();
  size_t filled_chars = static_cast<size_t>(percentage * bar_width_);
  
  // Filled portion
  for (size_t i = 0; i < filled_chars; ++i) {
    bar += "█";
  }
  
  // Empty portion
  for (size_t i = filled_chars; i < bar_width_; ++i) {
    bar += "░";
  }
  
  bar += "]";
  return bar;
}

bool ProgressBar::shouldUpdate() const
{
  return (current_item_ % update_frequency_ == 0) || 
         (current_item_ == total_items_);
}