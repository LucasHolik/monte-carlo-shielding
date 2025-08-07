#include "../include/PhotonCrossSectionDatabase.hpp"
#include <iostream>
#include <iomanip>

/**
 * @brief Debug program to identify and fix coherent scattering and photoelectric issues
 */
int main() {
    std::cout << "=== Debugging Cross-Section Issues ===\n";
    std::cout << std::fixed << std::setprecision(6);
    
    // Test Iron (Z=26) photoelectric at 100 keV
    int Z_iron = 26;
    double energy_kev = 100.0;
    
    std::cout << "\n1. Iron (Z=26) binding energies:\n";
    double K_edge = PhotonCrossSections::ScofieldPhotoelectricData::getBindingEnergy(26, 'K');
    double L_edge = PhotonCrossSections::ScofieldPhotoelectricData::getBindingEnergy(26, 'L');
    double M_edge = PhotonCrossSections::ScofieldPhotoelectricData::getBindingEnergy(26, 'M');
    
    std::cout << "  K-edge: " << K_edge << " keV\n";
    std::cout << "  L-edge: " << L_edge << " keV\n";  
    std::cout << "  M-edge: " << M_edge << " keV\n";
    
    // Check which shells are accessible at 100 keV
    std::cout << "\n2. Shell accessibility at 100 keV:\n";
    std::cout << "  Above K-edge? " << (energy_kev >= K_edge) << " (need " << K_edge << " keV)\n";
    std::cout << "  Above L-edge? " << (energy_kev >= L_edge) << " (need " << L_edge << " keV)\n";
    std::cout << "  Above M-edge? " << (energy_kev >= M_edge) << " (need " << M_edge << " keV)\n";
    
    // Get individual shell cross-sections if we had them
    double pe_total = PhotonCrossSections::ScofieldPhotoelectricData::getPhotoelectricCrossSection(Z_iron, energy_kev);
    std::cout << "\n3. Total photoelectric cross-section: " << pe_total << " barns\n";
    
    // Test coherent scattering for elements with and without form factor data
    std::cout << "\n4. Coherent scattering test:\n";
    
    // Carbon (Z=6) - should have form factor data
    double carbon_coherent = PhotonCrossSections::CoherentScatteringData::getCoherentCrossSection(6, energy_kev);
    std::cout << "  Carbon (Z=6): " << carbon_coherent << " barns\n";
    
    // Iron (Z=26) - likely missing form factor data  
    double iron_coherent = PhotonCrossSections::CoherentScatteringData::getCoherentCrossSection(26, energy_kev);
    std::cout << "  Iron (Z=26): " << iron_coherent << " barns\n";
    
    // Test form factor directly
    double iron_form_factor = PhotonCrossSections::CoherentScatteringData::getAtomicFormFactor(26, 0.1);
    std::cout << "  Iron form factor at q=0.1: " << iron_form_factor << "\n";
    
    // Lead (Z=82) - definitely missing form factor data
    double lead_coherent = PhotonCrossSections::CoherentScatteringData::getCoherentCrossSection(82, energy_kev);
    std::cout << "  Lead (Z=82): " << lead_coherent << " barns\n";
    
    std::cout << "\n=== Analysis Complete ===\n";
    
    return 0;
}