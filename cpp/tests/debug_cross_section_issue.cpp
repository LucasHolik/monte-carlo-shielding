#include "../include/PhotonCrossSectionDatabase.hpp"
#include <iostream>
#include <iomanip>

int main()
{
    std::cout << std::fixed << std::setprecision(6);
    
    // Test lead binding energy
    int Z_lead = 82;
    double K_edge = PhotonCrossSections::ScofieldPhotoelectricData::getBindingEnergy(Z_lead, 'K');
    double L_edge = PhotonCrossSections::ScofieldPhotoelectricData::getBindingEnergy(Z_lead, 'L');
    double M_edge = PhotonCrossSections::ScofieldPhotoelectricData::getBindingEnergy(Z_lead, 'M');
    
    std::cout << "Lead binding energies:" << std::endl;
    std::cout << "  K-edge: " << K_edge << " keV" << std::endl;
    std::cout << "  L-edge: " << L_edge << " keV" << std::endl;
    std::cout << "  M-edge: " << M_edge << " keV" << std::endl;
    
    // Test if 10 keV is above various edges
    std::cout << "\n10 keV above edges:" << std::endl;
    std::cout << "  K: " << (PhotonCrossSections::ScofieldPhotoelectricData::isAboveAbsorptionEdge(Z_lead, 10.0, 'K') ? "Yes" : "No") << std::endl;
    std::cout << "  L: " << (PhotonCrossSections::ScofieldPhotoelectricData::isAboveAbsorptionEdge(Z_lead, 10.0, 'L') ? "Yes" : "No") << std::endl;
    std::cout << "  M: " << (PhotonCrossSections::ScofieldPhotoelectricData::isAboveAbsorptionEdge(Z_lead, 10.0, 'M') ? "Yes" : "No") << std::endl;
    
    // Test photoelectric cross-section at various energies
    std::vector<double> energies = {5.0, 10.0, 20.0, 50.0, 100.0};
    std::cout << "\nLead photoelectric cross-sections:" << std::endl;
    for(double energy : energies) {
        double xs = PhotonCrossSections::ScofieldPhotoelectricData::getPhotoelectricCrossSection(Z_lead, energy);
        std::cout << "  " << energy << " keV: " << xs << " barns" << std::endl;
    }
    
    // Test a lighter element (carbon) for comparison
    int Z_carbon = 6;
    std::cout << "\nCarbon photoelectric @ 10 keV: " 
              << PhotonCrossSections::ScofieldPhotoelectricData::getPhotoelectricCrossSection(Z_carbon, 10.0) 
              << " barns" << std::endl;
    
    return 0;
}