#include "ProgressBar.hpp"
#include <algorithm>
#include <iomanip>

ProgressBar::ProgressBar(size_t total_items)
    : total_items_(total_items), current_item_(0), is_started_(false)
{
  if(total_items_ == 0)
  {
    total_items_ = 1;
  }
}

void ProgressBar::start()
{
  start_time_ = std::chrono::steady_clock::now();
  current_item_ = 0;
  is_started_ = true;
  display();
}

void ProgressBar::update()
{
  if(!is_started_)
  {
    start();
  }

  current_item_++;
  display();
}

void ProgressBar::updateTo(size_t current_item)
{
  if(!is_started_)
  {
    start();
  }

  current_item_ = std::min(current_item, total_items_);
  display();
}

void ProgressBar::finish()
{
  if(!is_started_)
  {
    start();
  }

  current_item_ = total_items_;
  display();
  std::cout << std::endl;
}

double ProgressBar::getProgressPercentage() const
{
  if(total_items_ == 0)
    return 1.0;
  return static_cast<double>(current_item_) / static_cast<double>(total_items_);
}

double ProgressBar::getElapsedSeconds() const
{
  if(!is_started_)
    return 0.0;

  auto now = std::chrono::steady_clock::now();
  auto elapsed =
      std::chrono::duration_cast<std::chrono::milliseconds>(now - start_time_);
  return elapsed.count() / 1000.0;
}

double ProgressBar::getETASeconds() const
{
  if(!is_started_ || current_item_ == 0)
    return 0.0;

  double elapsed = getElapsedSeconds();
  double rate = static_cast<double>(current_item_) / elapsed;
  double remaining_items = static_cast<double>(total_items_ - current_item_);

  if(rate <= 0.0)
    return 0.0;
  return remaining_items / rate;
}

void ProgressBar::display()
{
  const int barWidth = 20;
  double progress = getProgressPercentage();

  std::cout << "\r[";
  int pos = static_cast<int>(barWidth * progress);
  for(int i = 0; i < barWidth; ++i)
  {
    if(i < pos)
      std::cout << "=";
    else if(i == pos)
      std::cout << ">";
    else
      std::cout << " ";
  }
  std::cout << "] " << std::setw(3) << static_cast<int>(progress * 100.0)
            << "%";

  if(is_started_)
  {
    double elapsed = getElapsedSeconds();
    std::cout << " | " << std::fixed << std::setprecision(1) << elapsed << "s";

    if(current_item_ > 0 && current_item_ < total_items_)
    {
      double eta = getETASeconds();
      std::cout << " | ETA: " << std::fixed << std::setprecision(1) << eta
                << "s";
    }
  }

  std::cout.flush();
}