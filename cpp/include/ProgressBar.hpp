#pragma once

#include <chrono>
#include <iostream>
#include <string>

/**
 * @brief A comprehensive progress bar for long-running Monte Carlo simulations
 * 
 * Features:
 * - Real-time percentage and count display
 * - Elapsed time tracking
 * - Estimated time remaining (ETA)
 * - Visual ASCII progress bar
 * - Configurable update frequency
 */
class ProgressBar
{
private:
  size_t total_items_;
  size_t current_item_;
  size_t update_frequency_;
  std::chrono::steady_clock::time_point start_time_;
  std::chrono::steady_clock::time_point last_update_;
  bool is_started_;
  bool show_bar_;
  size_t bar_width_;
  
public:
  /**
   * @brief Construct a new Progress Bar
   * @param total_items Total number of items to process
   * @param update_frequency Update display every N items (default: 10)
   * @param bar_width Width of visual progress bar in characters (default: 40)
   * @param show_bar Whether to show visual bar (default: true)
   */
  ProgressBar(size_t total_items, 
              size_t update_frequency = 10, 
              size_t bar_width = 40, 
              bool show_bar = true);
  
  /**
   * @brief Start the progress tracking
   */
  void start();
  
  /**
   * @brief Update progress by one item
   * @param force_display Force display even if update frequency not reached
   */
  void update(bool force_display = false);
  
  /**
   * @brief Update progress to specific item count
   * @param current_item Current item number (0-based)
   * @param force_display Force display even if update frequency not reached
   */
  void updateTo(size_t current_item, bool force_display = false);
  
  /**
   * @brief Mark as complete and show final summary
   */
  void finish();
  
  /**
   * @brief Get current progress as percentage (0.0 to 1.0)
   */
  double getProgressPercentage() const;
  
  /**
   * @brief Get elapsed time in seconds
   */
  double getElapsedSeconds() const;
  
  /**
   * @brief Get estimated time remaining in seconds
   */
  double getETASeconds() const;
  
  /**
   * @brief Get current item count
   */
  size_t getCurrentItem() const { return current_item_; }
  
  /**
   * @brief Get total item count
   */
  size_t getTotalItems() const { return total_items_; }
  
private:
  /**
   * @brief Display the current progress
   */
  void display();
  
  /**
   * @brief Format time duration as human readable string
   * @param seconds Time in seconds
   * @return Formatted string (e.g., "2m 30s", "45.2s")
   */
  std::string formatTime(double seconds) const;
  
  /**
   * @brief Create visual progress bar string
   * @return ASCII progress bar
   */
  std::string createProgressBar() const;
  
  /**
   * @brief Check if display should be updated
   */
  bool shouldUpdate() const;
};