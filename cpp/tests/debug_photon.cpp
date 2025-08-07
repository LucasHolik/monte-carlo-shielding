#include "../include/PhotonPhysics.hpp"
#include "../include/RandomNumberGenerator.hpp"
#include "../include/MonteCarloSampling.hpp"
#include <iostream>
#include <cmath>

int main() {
    RandomNumberGenerator rng(12345);
    MonteCarloSampling sampling(rng);
    PhotonPhysics physics(sampling);
    
    // Test Klein-Nishina cross-section at low energy
    double low_energy = 0.001; // 1 keV
    double thomson_xs = physics.kleinNishinaCrossSection(low_energy);
    double expected_thomson = 8.0 / 3.0 * M_PI * std::pow(2.8179403262e-13, 2);
    
    std::cout << "Low energy: " << low_energy << " MeV" << std::endl;
    std::cout << "Calculated Klein-Nishina: " << thomson_xs << " cm²" << std::endl;
    std::cout << "Expected Thomson: " << expected_thomson << " cm²" << std::endl;
    std::cout << "Relative error: " << std::abs(thomson_xs - expected_thomson) / expected_thomson << std::endl;
    
    return 0;
}