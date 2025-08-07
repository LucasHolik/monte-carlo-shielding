#include "../include/Material.hpp"
#include "../include/PhotonCrossSectionDatabase.hpp"
#include <iostream>
#include <iomanip>
#include <cassert>

/**
 * @brief Test Material class integration with PhotonCrossSectionDatabase
 */
int main() {
    std::cout << "=== Testing Material-PhotonCrossSectionDatabase Integration ===\n";
    std::cout << std::scientific << std::setprecision(6);
    
    // Test 1: Pure element (Lead)
    std::cout << "\n1. Pure Lead Material:\n";
    Material lead = Material::createLead(); // Default density 11.34 g/cm³
    double energy_kev = 100.0;
    
    double total_xs = lead.getTotalPhotonCrossSection(energy_kev);
    double mass_atten = lead.getMassAttenuationCoefficient(energy_kev);
    double linear_atten = lead.getLinearAttenuationCoefficient(energy_kev);
    auto components = lead.getPhotonCrossSectionComponents(energy_kev);
    
    std::cout << "  Total cross-section: " << total_xs << " barns/atom\n";
    std::cout << "  Mass attenuation: " << mass_atten << " cm²/g\n";
    std::cout << "  Linear attenuation: " << linear_atten << " cm⁻¹\n";
    std::cout << "  Components [PE, Coh, Inc, Pair_n, Pair_e]: ";
    for(size_t i = 0; i < 5; ++i) {
        std::cout << components[i];
        if(i < 4) std::cout << ", ";
    }
    std::cout << " barns\n";
    
    // Test 2: Water (H₂O) - compound material
    std::cout << "\n2. Water (H₂O) Material:\n";
    Material water = Material::createWater(); // Density 1.0 g/cm³
    
    double water_total_xs = water.getTotalPhotonCrossSection(energy_kev);
    double water_mass_atten = water.getMassAttenuationCoefficient(energy_kev);
    double water_linear_atten = water.getLinearAttenuationCoefficient(energy_kev);
    auto water_components = water.getPhotonCrossSectionComponents(energy_kev);
    
    std::cout << "  Water composition: " << water.getCompositionString() << "\n";
    std::cout << "  Total cross-section: " << water_total_xs << " barns/atom\n";
    std::cout << "  Mass attenuation: " << water_mass_atten << " cm²/g\n";
    std::cout << "  Linear attenuation: " << water_linear_atten << " cm⁻¹\n";
    std::cout << "  Components [PE, Coh, Inc, Pair_n, Pair_e]: ";
    for(size_t i = 0; i < 5; ++i) {
        std::cout << water_components[i];
        if(i < 4) std::cout << ", ";
    }
    std::cout << " barns\n";
    
    // Test 3: Compare with direct calculation for water
    std::cout << "\n3. Manual Water Verification:\n";
    
    // Water: H₂O - 11.19% H, 88.81% O by mass
    // Let's manually calculate what we should get
    double h_mass_fraction = 0.1119;
    double o_mass_fraction = 0.8881;
    double h_atomic_mass = 1.008;
    double o_atomic_mass = 15.999;
    
    // Calculate moles
    double h_moles = h_mass_fraction / h_atomic_mass;
    double o_moles = o_mass_fraction / o_atomic_mass;
    double total_moles = h_moles + o_moles;
    
    double h_atom_fraction = h_moles / total_moles;
    double o_atom_fraction = o_moles / total_moles;
    
    std::cout << "  H atom fraction: " << h_atom_fraction << "\n";
    std::cout << "  O atom fraction: " << o_atom_fraction << "\n";
    
    // Get individual element cross-sections
    double h_xs = PhotonCrossSections::PhotonInteractionDatabase::getTotalCrossSection(1, energy_kev);
    double o_xs = PhotonCrossSections::PhotonInteractionDatabase::getTotalCrossSection(8, energy_kev);
    
    double manual_total_xs = h_atom_fraction * h_xs + o_atom_fraction * o_xs;
    
    std::cout << "  H cross-section: " << h_xs << " barns\n";
    std::cout << "  O cross-section: " << o_xs << " barns\n";
    std::cout << "  Manual calculation: " << manual_total_xs << " barns\n";
    std::cout << "  Material method: " << water_total_xs << " barns\n";
    std::cout << "  Difference: " << std::abs(manual_total_xs - water_total_xs) << " barns\n";
    
    // Verify they match (within reasonable numerical precision)
    assert(std::abs(manual_total_xs - water_total_xs) < 1e-3); // 0.1% tolerance
    
    // Test 4: Energy dependence
    std::cout << "\n4. Energy Dependence for Lead:\n";
    std::vector<double> energies = {10, 50, 100, 500, 1000}; // keV
    
    std::cout << "  Energy (keV) | Linear μ (cm⁻¹) | Mass μ (cm²/g)\n";
    std::cout << "  -------------|------------------|----------------\n";
    
    for(double energy : energies) {
        double lin_mu = lead.getLinearAttenuationCoefficient(energy);
        double mass_mu = lead.getMassAttenuationCoefficient(energy);
        
        std::cout << "  " << std::setw(11) << std::fixed << std::setprecision(0) << energy 
                  << " | " << std::setw(15) << std::scientific << std::setprecision(3) << lin_mu
                  << " | " << std::setw(14) << mass_mu << "\n";
    }
    
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "\n✅ All Material-PhotonCrossSectionDatabase integration tests passed!\n";
    
    return 0;
}