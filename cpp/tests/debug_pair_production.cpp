#include "../include/PhotonCrossSectionDatabase.hpp"
#include <iostream>
#include <iomanip>

int main()
{
    std::cout << std::fixed << std::setprecision(6);
    
    // Debug pair production for iron at different energies
    int Z = 26; // Iron
    std::vector<double> energies_MeV = {0.5, 1.0, 1.5, 2.0, 5.0, 10.0, 20.0};
    
    std::cout << "=== Pair Production Debug for Iron (Z=26) ===\n";
    
    for(double energy : energies_MeV) {
        bool possible = PhotonCrossSections::PairProductionData::isPairProductionPossible(energy);
        double nuclear_xs = PhotonCrossSections::PairProductionData::getPairProductionNuclearCrossSection(Z, energy);
        double electronic_xs = PhotonCrossSections::PairProductionData::getPairProductionElectronCrossSection(Z, energy);
        double total_xs = PhotonCrossSections::PairProductionData::getTotalPairProductionCrossSection(Z, energy);
        
        std::cout << energy << " MeV: ";
        std::cout << "Possible=" << (possible ? "Yes" : "No") << " ";
        std::cout << "Nuclear=" << nuclear_xs << " ";
        std::cout << "Electronic=" << electronic_xs << " ";
        std::cout << "Total=" << total_xs << "\n";
        
        if(total_xs < 0.0 && possible) {
            std::cout << "  ⚠ NEGATIVE CROSS-SECTION!\n";
        }
    }
    
    // Test the threshold
    double threshold = PhotonCrossSections::Constants::PAIR_PRODUCTION_THRESHOLD_MEV;
    std::cout << "\nPair production threshold: " << threshold << " MeV\n";
    
    return 0;
}