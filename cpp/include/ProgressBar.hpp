#pragma once

#include <chrono>
#include <iostream>

class ProgressBar
{
private:
  size_t total_items_;
  size_t current_item_;
  std::chrono::steady_clock::time_point start_time_;
  bool is_started_;
  
public:
  ProgressBar(size_t total_items);
  
  void start();
  void update();
  void updateTo(size_t current_item);
  void finish();
  
  double getProgressPercentage() const;
  double getElapsedSeconds() const;
  double getETASeconds() const;
  
private:
  void display();
};