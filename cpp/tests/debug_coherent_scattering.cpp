#include "../include/PhotonCrossSectionDatabase.hpp"
#include <iostream>
#include <iomanip>

int main()
{
    std::cout << std::fixed << std::setprecision(6);
    
    // Debug coherent scattering for different elements
    std::vector<int> elements = {1, 6, 26, 82}; // H, C, Fe, Pb
    double energy_keV = 100.0;
    
    std::cout << "=== Coherent Scattering Debug @ 100 keV ===\n";
    
    for(int Z : elements) {
        double coherent_xs = PhotonCrossSections::CoherentScatteringData::getCoherentCrossSection(Z, energy_keV);
        double thomson_xs = PhotonCrossSections::CoherentScatteringData::getThomsonCrossSection();
        
        std::cout << "Z=" << Z << ": Coherent=" << coherent_xs 
                  << " barns, Thomson=" << thomson_xs << " barns\n";
        
        // Test form factor calculation
        double q = 0.1; // Small momentum transfer
        double form_factor = PhotonCrossSections::CoherentScatteringData::getAtomicFormFactor(Z, q);
        
        std::cout << "  Form factor @ q=0.1: " << form_factor << "\n";
        
        if(coherent_xs == 0.0 && Z > 10) {
            std::cout << "  ⚠ Zero coherent scattering for heavy element!\n";
        }
    }
    
    // Test form factor data availability
    std::cout << "\n=== Form Factor Data Test ===\n";
    for(int Z = 1; Z <= 10; Z++) {
        double ff = PhotonCrossSections::CoherentScatteringData::getAtomicFormFactor(Z, 0.1);
        std::cout << "Z=" << Z << " form factor: " << ff << "\n";
    }
    
    return 0;
}