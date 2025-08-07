#include "../include/PhotonPhysics.hpp"
#include "../include/Material.hpp"
#include "../include/RandomNumberGenerator.hpp"
#include "../include/MonteCarloSampling.hpp"
#include <iostream>
#include <map>

int main() {
    RandomNumberGenerator rng(12345);
    MonteCarloSampling sampling(rng);
    PhotonPhysics physics(sampling);
    
    // Create test material (aluminum - lighter element)
    Material aluminum("Aluminum", 2.70);
    aluminum.addElement(13, 26.98, 1.0, "Al");
    double energy = 1.0; // MeV
    
    // Check individual cross-sections
    std::cout << "Individual cross-sections for " << energy << " MeV photon in aluminum:" << std::endl;
    std::cout << "Photoelectric: " << physics.getPhotoelectricCrossSection(energy, 13) << std::endl;
    std::cout << "Compton: " << physics.getComptonCrossSection(energy, 13) << std::endl;
    std::cout << "Coherent: " << physics.getCoherentCrossSection(energy, 13) << std::endl;
    std::cout << "Pair production: " << physics.getPairProductionCrossSection(energy, 13) << std::endl;
    std::cout << "Total: " << physics.getTotalCrossSection(energy, aluminum) << std::endl;
    
    // Sample a few interactions to check variety
    std::map<PhotonInteractionType, int> interaction_counts;
    
    for(int i = 0; i < 100; i++) {
        auto interaction_type = physics.sampleInteractionType(energy, aluminum);
        interaction_counts[interaction_type]++;
    }
    
    std::cout << "\nSampled interactions (100 samples):" << std::endl;
    for(const auto& pair : interaction_counts) {
        std::cout << physics.getInteractionTypeName(pair.first) << ": " << pair.second << std::endl;
    }
    
    return 0;
}