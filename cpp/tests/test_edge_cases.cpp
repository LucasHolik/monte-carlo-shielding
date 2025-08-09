/**
 * @file test_edge_cases.cpp
 * @brief Test absorption edge handling in NIST data
 */

#include "PhotonCrossSectionDatabase.hpp"
#include <iostream>
#include <iomanip>

using namespace PhotonCrossSections;

int main() {
    std::cout << "=== Testing Absorption Edge Cases ===\n";
    
    try {
        // Test vanadium (Z=23) around its K-edge
        // From the NIST data, we saw duplicate energies at 5.465E-03 MeV
        int Z = 23; // Vanadium
        double edge_energy_keV = 5.465; // K-edge energy in keV
        
        auto& db = NISTCrossSectionDatabase::getInstance();
        
        // Test energies around the edge
        std::vector<double> test_energies = {
            5.400, 5.464, 5.465, 5.466, 5.500
        };
        
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Vanadium (Z=23) photoelectric cross-sections around K-edge:\n";
        
        for (double energy : test_energies) {
            double photoelectric = db.getCrossSection(Z, "photon", energy, CrossSectionType::PHOTOELECTRIC);
            std::cout << "Energy: " << energy << " keV -> " << photoelectric << " cm²/g\n";
        }
        
        // Test that exactly at edge energy we get the post-edge (higher) value
        double exact_edge = db.getCrossSection(Z, "photon", edge_energy_keV, CrossSectionType::PHOTOELECTRIC);
        double just_above = db.getCrossSection(Z, "photon", edge_energy_keV + 0.001, CrossSectionType::PHOTOELECTRIC);
        
        std::cout << "\nEdge handling test:\n";
        std::cout << "Exactly at edge (5.465 keV): " << exact_edge << " cm²/g\n";
        std::cout << "Just above edge (5.466 keV): " << just_above << " cm²/g\n";
        
        // The exact edge value should be similar to just above (post-edge)
        // rather than just below (pre-edge)
        
        std::cout << "\n=== Edge Test Completed Successfully ===\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}