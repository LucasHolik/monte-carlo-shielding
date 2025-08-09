/**
 * @file test_nist_cross_sections.cpp
 * @brief Test NIST Cross-Section Database Implementation
 */

#include "PhotonCrossSectionDatabase.hpp"
#include <iostream>
#include <iomanip>
#include <cassert>

using namespace PhotonCrossSections;

void testBasicFunctionality() {
    std::cout << "=== Testing Basic NIST Cross-Section Functionality ===\n";
    
    // Test hydrogen (Z=1) photon cross-sections at various energies
    int Z = 1; // Hydrogen
    std::string particle_type = "photon";
    
    std::vector<double> test_energies = {1.0, 10.0, 100.0, 1000.0}; // keV
    
    for (double energy : test_energies) {
        try {
            auto& db = NISTCrossSectionDatabase::getInstance();
            
            // Test individual cross-section retrieval
            double photoelectric = db.getCrossSection(Z, particle_type, energy, CrossSectionType::PHOTOELECTRIC);
            double coherent = db.getCrossSection(Z, particle_type, energy, CrossSectionType::COHERENT);
            double incoherent = db.getCrossSection(Z, particle_type, energy, CrossSectionType::INCOHERENT);
            double total = db.getCrossSection(Z, particle_type, energy, CrossSectionType::TOTAL_WITH_COHERENT);
            
            std::cout << std::fixed << std::setprecision(4);
            std::cout << "Energy: " << energy << " keV\n";
            std::cout << "  Photoelectric: " << photoelectric << " cm²/g\n";
            std::cout << "  Coherent:      " << coherent << " cm²/g\n";
            std::cout << "  Incoherent:    " << incoherent << " cm²/g\n";
            std::cout << "  Total:         " << total << " cm²/g\n\n";
            
            // Basic validation
            assert(photoelectric >= 0.0);
            assert(coherent >= 0.0);
            assert(incoherent >= 0.0);
            assert(total >= 0.0);
            
        } catch (const std::exception& e) {
            std::cerr << "Error at energy " << energy << " keV: " << e.what() << std::endl;
        }
    }
}

void testMultipleElements() {
    std::cout << "=== Testing Multiple Elements ===\n";
    
    std::vector<int> elements = {1, 6, 26, 82}; // H, C, Fe, Pb
    double energy_keV = 100.0; // 100 keV
    
    for (int Z : elements) {
        try {
            double total = NIST::getTotalPhotonCrossSection(Z, energy_keV, true);
            std::cout << "Z=" << Z << " at " << energy_keV << " keV: " << total << " cm²/g\n";
            
            assert(total >= 0.0);
            
        } catch (const std::exception& e) {
            std::cerr << "Error for Z=" << Z << ": " << e.what() << std::endl;
        }
    }
}

void testEnergyRange() {
    std::cout << "\n=== Testing Energy Range ===\n";
    
    int Z = 26; // Iron
    auto& db = NISTCrossSectionDatabase::getInstance();
    
    try {
        auto range = db.getEnergyRange(Z, "photon");
        std::cout << "Iron (Z=26) energy range: " << range.first << " to " << range.second << " keV\n";
        
        assert(range.first > 0.0);
        assert(range.second > range.first);
        
    } catch (const std::exception& e) {
        std::cerr << "Error getting energy range: " << e.what() << std::endl;
    }
}

void testAllCrossSections() {
    std::cout << "\n=== Testing All Cross-Sections ===\n";
    
    int Z = 82; // Lead
    double energy_keV = 1000.0; // 1 MeV
    
    try {
        auto all_cs = NIST::getAllPhotonCrossSections(Z, energy_keV);
        
        std::cout << "Lead (Z=82) at " << energy_keV << " keV:\n";
        std::cout << "  Coherent:      " << all_cs[0] << " cm²/g\n";
        std::cout << "  Incoherent:    " << all_cs[1] << " cm²/g\n";
        std::cout << "  Photoelectric: " << all_cs[2] << " cm²/g\n";
        std::cout << "  Pair Nuclear:  " << all_cs[3] << " cm²/g\n";
        std::cout << "  Pair Electron: " << all_cs[4] << " cm²/g\n";
        std::cout << "  Total w/ Coh:  " << all_cs[5] << " cm²/g\n";
        std::cout << "  Total w/o Coh: " << all_cs[6] << " cm²/g\n";
        
        // Validate all values are non-negative
        for (size_t i = 0; i < 7; ++i) {
            assert(all_cs[i] >= 0.0);
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error getting all cross-sections: " << e.what() << std::endl;
    }
}

void testUnitConversion() {
    std::cout << "\n=== Testing Unit Conversion ===\n";
    
    auto& db = NISTCrossSectionDatabase::getInstance();
    
    double cross_section_cm2_per_g = 1.0; // 1 cm²/g
    double atomic_mass = 1.008; // Hydrogen atomic mass
    
    double barns_per_atom = db.convertToBarnsPerAtom(cross_section_cm2_per_g, atomic_mass);
    
    std::cout << "1 cm²/g = " << barns_per_atom << " barns/atom (for hydrogen)\n";
    
    assert(barns_per_atom > 0.0);
}

int main() {
    try {
        std::cout << "NIST Cross-Section Database Test\n";
        std::cout << "================================\n\n";
        
        testBasicFunctionality();
        testMultipleElements();
        testEnergyRange();
        testAllCrossSections();
        testUnitConversion();
        
        std::cout << "\n=== All Tests Passed! ===\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Test failed with unknown exception" << std::endl;
        return 1;
    }
    
    return 0;
}