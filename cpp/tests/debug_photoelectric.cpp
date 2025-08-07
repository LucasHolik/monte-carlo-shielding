#include "../include/PhotonCrossSectionDatabase.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

/**
 * @brief Debug photoelectric cross-section calculation step by step
 */
int main() {
    std::cout << "=== Debugging Photoelectric Cross-Section ===\n";
    std::cout << std::scientific << std::setprecision(6);
    
    int Z = 26; // Iron
    double energy_keV = 100.0;
    
    // Check constants
    std::cout << "Constants:\n";
    std::cout << "  Fine structure constant: " << PhotonCrossSections::Constants::FINE_STRUCTURE_CONSTANT << "\n";
    std::cout << "  Classical electron radius: " << PhotonCrossSections::Constants::CLASSICAL_ELECTRON_RADIUS_CM << " cm\n";
    std::cout << "  Pi: " << M_PI << "\n";
    
    // Get binding energies
    std::cout << "\nIron binding energies:\n";
    double K_edge = PhotonCrossSections::ScofieldPhotoelectricData::getBindingEnergy(26, 'K');
    double L_edge = PhotonCrossSections::ScofieldPhotoelectricData::getBindingEnergy(26, 'L'); 
    double M_edge = PhotonCrossSections::ScofieldPhotoelectricData::getBindingEnergy(26, 'M');
    
    std::cout << "  K-edge: " << K_edge << " keV\n";
    std::cout << "  L-edge: " << L_edge << " keV\n";
    std::cout << "  M-edge: " << M_edge << " keV\n";
    
    // Manual calculation for K-shell
    std::cout << "\nManual K-shell calculation:\n";
    double eta = energy_keV / K_edge;
    double alpha = PhotonCrossSections::Constants::FINE_STRUCTURE_CONSTANT;
    double r_e_squared = PhotonCrossSections::Constants::CLASSICAL_ELECTRON_RADIUS_CM * 
                        PhotonCrossSections::Constants::CLASSICAL_ELECTRON_RADIUS_CM;
    
    std::cout << "  eta (E/E_K): " << eta << "\n";
    std::cout << "  alpha: " << alpha << "\n";
    std::cout << "  r_e^2: " << r_e_squared << " cm^2\n";
    std::cout << "  (alpha*Z)^4: " << std::pow(alpha * Z, 4) << "\n";
    std::cout << "  eta^(-3.5): " << std::pow(eta, -3.5) << "\n";
    
    double sigma_K_cm2 = (32.0 * M_PI / 3.0) * r_e_squared * std::pow(alpha * Z, 4) * std::pow(eta, -3.5);
    double sigma_K_barns = sigma_K_cm2 * 1e24;
    
    std::cout << "  K-shell cross-section: " << sigma_K_cm2 << " cm^2\n";
    std::cout << "  K-shell cross-section: " << sigma_K_barns << " barns\n";
    
    // Test the actual function call  
    double total_pe = PhotonCrossSections::ScofieldPhotoelectricData::getPhotoelectricCrossSection(Z, energy_keV);
    std::cout << "\nTotal photoelectric from function: " << total_pe << " barns\n";
    
    return 0;
}