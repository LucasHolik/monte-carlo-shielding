#include "cpp/include/ProgressBar.hpp"
#include <thread>
#include <chrono>

int main() {
    std::cout << "Testing Rewritten Progress Bar..." << std::endl;
    
    ProgressBar pb(100);
    pb.start();
    
    for (int i = 1; i <= 100; ++i) {
        pb.updateTo(i);
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
    
    pb.finish();
    
    std::cout << std::endl << "Test completed successfully!" << std::endl;
    return 0;
}