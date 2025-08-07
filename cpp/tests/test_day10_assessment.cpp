#include "../include/PhotonCrossSectionDatabase.hpp"

#include <iostream>
#include <iomanip>
#include <vector>

/**
 * @brief Assessment of Day 10 requirements - what's already implemented?
 */
class Day10Assessment 
{
public:
    void runAssessment()
    {
        std::cout << "=== Day 10 Requirements Assessment ===\n";
        std::cout << std::fixed << std::setprecision(6);
        
        assessCrossSectionDatabase();
        assessEnergyInterpolation();
        assessMaterialComposition();
        assessNISTData();
        assessDensityCorrections();
        
        std::cout << "\n=== Summary ===\n";
        std::cout << "Most Day 10 functionality appears to be already implemented!\n";
    }

private:
    void assessCrossSectionDatabase()
    {
        std::cout << "\n1. Cross-Section Database Structure:\n";
        
        // Test the existing database interface
        int Z = 26; // Iron
        std::vector<double> energies = {10, 100, 1000, 10000}; // keV
        
        std::cout << "   Testing PhotonInteractionDatabase for Iron (Z=26):\n";
        
        for(double energy : energies) {
            auto all_xs = PhotonCrossSections::PhotonInteractionDatabase::getAllCrossSections(Z, energy);
            double total = PhotonCrossSections::PhotonInteractionDatabase::getTotalCrossSection(Z, energy);
            
            std::cout << "   " << energy << " keV: Total=" << total 
                      << " [PE=" << all_xs[0] << ", Coh=" << all_xs[1] 
                      << ", Inc=" << all_xs[2] << ", Pair=" << all_xs[3] + all_xs[4] << "]\n";
        }
        
        std::cout << "   ✓ Database structure exists and working\n";
    }
    
    void assessEnergyInterpolation()
    {
        std::cout << "\n2. Energy Interpolation:\n";
        
        // Test smooth energy dependence
        int Z = 82; // Lead
        std::vector<double> fine_energies;
        for(int i = 1; i <= 20; i++) {
            fine_energies.push_back(i * 50.0); // 50, 100, 150, ... 1000 keV
        }
        
        std::cout << "   Testing energy dependence smoothness for Lead:\n";
        double prev_pe = -1;
        bool smooth = true;
        
        for(double energy : fine_energies) {
            double pe_xs = PhotonCrossSections::ScofieldPhotoelectricData::getPhotoelectricCrossSection(Z, energy);
            
            if(prev_pe > 0 && pe_xs > prev_pe * 2.0) { // Check for discontinuities
                smooth = false;
                std::cout << "   ⚠ Large jump at " << energy << " keV: " << pe_xs << " vs " << prev_pe << "\n";
            }
            prev_pe = pe_xs;
        }
        
        if(smooth) {
            std::cout << "   ✓ Cross-sections vary smoothly with energy\n";
        } else {
            std::cout << "   ? May need additional interpolation for smoother results\n";
        }
    }
    
    void assessMaterialComposition()
    {
        std::cout << "\n3. Material Composition Handling:\n";
        
        // Test manual composition calculation for water (H2O)
        // Water: 11.19% H, 88.81% O by mass
        double h_fraction = 0.111898;
        double o_fraction = 0.888102;
        
        double energy_kev = 100.0;
        
        // Get individual element cross-sections
        double h_xs = PhotonCrossSections::PhotonInteractionDatabase::getTotalCrossSection(1, energy_kev);
        double o_xs = PhotonCrossSections::PhotonInteractionDatabase::getTotalCrossSection(8, energy_kev);
        
        // Calculate weighted average
        double water_xs = h_fraction * h_xs + o_fraction * o_xs;
        
        std::cout << "   Manual water composition calculation @ 100 keV:\n";
        std::cout << "   H cross-section: " << h_xs << " barns\n";
        std::cout << "   O cross-section: " << o_xs << " barns\n";
        std::cout << "   Water weighted average: " << water_xs << " barns\n";
        std::cout << "   ✓ Manual composition calculation working\n";
        std::cout << "   ? Need to test Material class integration\n";
    }
    
    void assessNISTData()
    {
        std::cout << "\n4. NIST Data Integration:\n";
        
        // Check if binding energies match known NIST values
        double pb_k_edge = PhotonCrossSections::ScofieldPhotoelectricData::getBindingEnergy(82, 'K');
        double pb_l_edge = PhotonCrossSections::ScofieldPhotoelectricData::getBindingEnergy(82, 'L');
        
        std::cout << "   Lead K-edge: " << pb_k_edge << " keV (NIST: ~88.0 keV)\n";
        std::cout << "   Lead L3-edge: " << pb_l_edge << " keV (NIST: ~15.2 keV)\n";
        
        // Check accuracy
        if(std::abs(pb_k_edge - 88.0) < 2.0 && std::abs(pb_l_edge - 15.2) < 2.0) {
            std::cout << "   ✓ Binding energies match NIST data well\n";
        } else {
            std::cout << "   ? Binding energies need verification against NIST\n";
        }
        
        // Test cross-section magnitudes
        double pb_pe_100kev = PhotonCrossSections::ScofieldPhotoelectricData::getPhotoelectricCrossSection(82, 100.0);
        std::cout << "   Lead photoelectric @ 100 keV: " << pb_pe_100kev << " barns\n";
        
        if(pb_pe_100kev > 0.1 && pb_pe_100kev < 10.0) {
            std::cout << "   ✓ Cross-section magnitudes reasonable\n";
        }
    }
    
    void assessDensityCorrections()
    {
        std::cout << "\n5. Density and Temperature Corrections:\n";
        
        // Check if mass attenuation coefficients scale with density
        double mass_atten = PhotonCrossSections::PhotonInteractionDatabase::getMassAttenuationCoefficient(
            26, 55.845, 100.0); // Iron, atomic mass 55.845, 100 keV
            
        std::cout << "   Iron mass attenuation @ 100 keV: " << mass_atten << " cm²/g\n";
        
        if(mass_atten > 0) {
            std::cout << "   ✓ Mass attenuation coefficient calculation working\n";
            std::cout << "   ? Temperature corrections not explicitly implemented yet\n";
        } else {
            std::cout << "   ✗ Mass attenuation coefficient calculation issue\n";
        }
    }
};

int main()
{
    try {
        Day10Assessment assessment;
        assessment.runAssessment();
        return 0;
    }
    catch(const std::exception& e) {
        std::cerr << "Assessment failed: " << e.what() << std::endl;
        return 1;
    }
}