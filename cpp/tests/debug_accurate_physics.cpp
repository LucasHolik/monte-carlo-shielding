#include "../include/PhotonCrossSectionDatabase.hpp"
#include <iostream>

int main() {
    // Debug photoelectric cross-section calculation
    int Z_lead = 82;
    double energy_keV = 10.0;
    
    // Check binding energies first
    double K_edge = PhotonCrossSections::ScofieldPhotoelectricData::getBindingEnergy(Z_lead, 'K');
    double L_edge = PhotonCrossSections::ScofieldPhotoelectricData::getBindingEnergy(Z_lead, 'L');
    double M_edge = PhotonCrossSections::ScofieldPhotoelectricData::getBindingEnergy(Z_lead, 'M');
    
    std::cout << "Lead K-edge: " << K_edge << " keV" << std::endl;
    std::cout << "Lead L-edge: " << L_edge << " keV" << std::endl;
    std::cout << "Lead M-edge: " << M_edge << " keV" << std::endl;
    std::cout << "Test energy: " << energy_keV << " keV" << std::endl;
    std::cout << "Above K-edge? " << (energy_keV >= K_edge) << std::endl;
    std::cout << "Above L-edge? " << (energy_keV >= L_edge) << std::endl;
    std::cout << "Above M-edge? " << (energy_keV >= M_edge) << std::endl;
    
    // Try photoelectric calculation
    double xs = PhotonCrossSections::ScofieldPhotoelectricData::getPhotoelectricCrossSection(Z_lead, energy_keV);
    std::cout << "Photoelectric cross-section: " << xs << " barns" << std::endl;
    
    // Check if it's above absorption edge
    bool above_edge = PhotonCrossSections::ScofieldPhotoelectricData::isAboveAbsorptionEdge(Z_lead, energy_keV, 'K');
    std::cout << "Above absorption edge: " << above_edge << std::endl;
    
    return 0;
}