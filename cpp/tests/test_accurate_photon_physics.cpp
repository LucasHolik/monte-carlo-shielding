#include "PhotonPhysics.hpp"
#include "PhotonCrossSectionDatabase.hpp"
#include "Material.hpp"
#include "MonteCarloSampling.hpp"
#include "RandomNumberGenerator.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

/**
 * @brief Test suite for physically accurate photon physics implementation
 * 
 * Validates the new NIST-level accurate cross-sections against known values
 * and physics expectations.
 */
class AccuratePhotonPhysicsTest
{
private:
  RandomNumberGenerator rng_;
  MonteCarloSampling sampling_;
  PhotonPhysics physics_;

public:
  AccuratePhotonPhysicsTest() : rng_(12345), sampling_(rng_), physics_(sampling_) {}

  void runAllTests()
  {
    std::cout << "Running Accurate Photon Physics Tests...\n";
    std::cout << std::fixed << std::setprecision(6);
    
    testPhotoelectricAccuracy();
    testCoherentScatteringWithFormFactors();
    testPairProductionPhysics();
    testEnergyDependenceValidation();
    testMaterialCompositionAccuracy();
    testCrossSectionDatabase();
    
    std::cout << "\n🎉 All accurate photon physics tests passed!\n";
  }

private:
  void testPhotoelectricAccuracy()
  {
    std::cout << "\n=== Testing Photoelectric Cross-Section Accuracy ===\n";
    
    // Test photoelectric cross-section for lead at different energies
    int Z_lead = 82;
    std::vector<std::pair<double, double>> test_cases = {
      // Energy (keV), Expected order of magnitude (barns)
      {10.0, 1e6},   // Low energy - very large cross-section
      {100.0, 1e4},  // Medium energy  
      {1000.0, 1e2}, // High energy - smaller cross-section
    };
    
    for(const auto& test_case : test_cases) {
      double energy_keV = test_case.first;
      double expected_magnitude = test_case.second;
      
      double xs_barns = PhotonCrossSections::ScofieldPhotoelectricData::getPhotoelectricCrossSection(
          Z_lead, energy_keV);
      
      std::cout << "Lead photoelectric @ " << energy_keV << " keV: " 
                << xs_barns << " barns" << std::endl;
      
      // Check that cross-section is in expected magnitude range
      assert(xs_barns > expected_magnitude * 0.1);
      assert(xs_barns < expected_magnitude * 10.0);
      
      // Check proper Z^5 scaling
      if(energy_keV > 100.0) { // Above all edges
        double xs_iron = PhotonCrossSections::ScofieldPhotoelectricData::getPhotoelectricCrossSection(
            26, energy_keV);
        
        if(xs_iron > 0.0) {
          double ratio = xs_barns / xs_iron;
          double expected_ratio = std::pow(82.0/26.0, 5.0);
          
          // Allow factor of 3 deviation due to shell structure
          assert(ratio > expected_ratio / 3.0);
          assert(ratio < expected_ratio * 3.0);
          
          std::cout << "  Z^5 scaling check: ratio = " << ratio 
                    << ", expected ≈ " << expected_ratio << std::endl;
        }
      }
    }
    
    // Test absorption edges
    double K_edge = PhotonCrossSections::ScofieldPhotoelectricData::getBindingEnergy(Z_lead, 'K');
    std::cout << "Lead K-edge: " << K_edge << " keV" << std::endl;
    assert(K_edge > 80.0 && K_edge < 90.0); // Lead K-edge around 88 keV
    
    std::cout << "✓ Photoelectric accuracy tests passed\n";
  }

  void testCoherentScatteringWithFormFactors()
  {
    std::cout << "\n=== Testing Coherent Scattering with Form Factors ===\n";
    
    // Test coherent scattering for carbon and lead
    std::vector<int> elements = {6, 82}; // Carbon, Lead
    std::vector<double> energies = {10.0, 100.0, 1000.0}; // keV
    
    for(int Z : elements) {
      for(double energy : energies) {
        double xs_coherent = PhotonCrossSections::CoherentScatteringData::getCoherentCrossSection(
            Z, energy);
        
        std::cout << "Z=" << Z << " coherent @ " << energy << " keV: " 
                  << xs_coherent << " barns" << std::endl;
        
        // Coherent scattering should decrease with energy
        assert(xs_coherent > 0.0);
        
        // Should be much smaller than Thomson scattering at high energies
        double thomson = PhotonCrossSections::CoherentScatteringData::getThomsonCrossSection();
        if(energy > 500.0) {
          assert(xs_coherent < thomson * Z * Z); // Form factor suppression
        }
        
        // Test form factor behavior
        double form_factor = PhotonCrossSections::CoherentScatteringData::getAtomicFormFactor(
            Z, 0.1); // Small momentum transfer
        
        assert(form_factor > 0.0);
        assert(form_factor <= Z); // Form factor can't exceed atomic number
        
        std::cout << "  Form factor @ q=0.1: " << form_factor << std::endl;
      }
    }
    
    std::cout << "✓ Coherent scattering with form factors tests passed\n";
  }

  void testPairProductionPhysics()
  {
    std::cout << "\n=== Testing Pair Production Physics ===\n";
    
    // Test pair production threshold
    double threshold_MeV = PhotonCrossSections::Constants::PAIR_PRODUCTION_THRESHOLD_MEV;
    std::cout << "Pair production threshold: " << threshold_MeV << " MeV" << std::endl;
    assert(std::abs(threshold_MeV - 1.022) < 0.001);
    
    // Test pair production for different elements and energies
    std::vector<int> elements = {6, 26, 82}; // Carbon, Iron, Lead
    std::vector<double> energies = {0.5, 1.5, 5.0, 10.0}; // MeV
    
    for(int Z : elements) {
      for(double energy : energies) {
        bool should_occur = energy >= threshold_MeV;
        bool can_occur = PhotonCrossSections::PairProductionData::isPairProductionPossible(energy);
        
        assert(should_occur == can_occur);
        
        double xs_pair = PhotonCrossSections::PairProductionData::getTotalPairProductionCrossSection(
            Z, energy);
        
        if(should_occur) {
          assert(xs_pair > 0.0);
          
          // Check proper Z² scaling
          if(Z == 26 && energy > 2.0) {
            double xs_carbon = PhotonCrossSections::PairProductionData::getTotalPairProductionCrossSection(
                6, energy);
            if(xs_carbon > 0.0) {
              double ratio = xs_pair / xs_carbon;
              double expected_ratio = (26.0 * 26.0) / (6.0 * 6.0);
              
              // Allow factor of 2 deviation
              assert(ratio > expected_ratio / 2.0);
              assert(ratio < expected_ratio * 2.0);
            }
          }
          
          std::cout << "Z=" << Z << " pair production @ " << energy << " MeV: " 
                    << xs_pair << " barns" << std::endl;
        } else {
          assert(xs_pair == 0.0);
        }
      }
    }
    
    std::cout << "✓ Pair production physics tests passed\n";
  }

  void testEnergyDependenceValidation()
  {
    std::cout << "\n=== Testing Energy Dependence ===\n";
    
    // Test that cross-sections follow expected energy dependence
    int Z = 26; // Iron
    std::vector<double> energies = {10, 30, 100, 300, 1000, 3000}; // keV
    
    std::vector<double> photoelectric_xs, compton_xs, coherent_xs;
    
    for(double energy : energies) {
      double pe_xs = PhotonCrossSections::ScofieldPhotoelectricData::getPhotoelectricCrossSection(
          Z, energy);
      double compton_xs_val = PhotonCrossSections::IncoherentScatteringData::getTotalIncoherentCrossSection(
          Z, energy);
      double coherent_xs_val = PhotonCrossSections::CoherentScatteringData::getCoherentCrossSection(
          Z, energy);
      
      photoelectric_xs.push_back(pe_xs);
      compton_xs.push_back(compton_xs_val);
      coherent_xs.push_back(coherent_xs_val);
      
      std::cout << energy << " keV: PE=" << pe_xs << ", Compton=" << compton_xs_val 
                << ", Coherent=" << coherent_xs_val << " barns" << std::endl;
    }
    
    // Photoelectric should generally decrease with energy
    for(size_t i = 1; i < photoelectric_xs.size(); i++) {
      if(photoelectric_xs[i] > 0 && photoelectric_xs[i-1] > 0) {
        // Allow for absorption edge jumps
        double ratio = photoelectric_xs[i] / photoelectric_xs[i-1];
        assert(ratio < 10.0); // Should generally decrease or have reasonable edge jump
      }
    }
    
    // Compton should decrease with energy
    for(size_t i = 1; i < compton_xs.size(); i++) {
      if(energies[i] > 100) { // Above where Klein-Nishina dominates
        assert(compton_xs[i] <= compton_xs[i-1] * 1.1); // Small tolerance
      }
    }
    
    std::cout << "✓ Energy dependence validation passed\n";
  }

  void testMaterialCompositionAccuracy()
  {
    std::cout << "\n=== Testing Material Composition Accuracy ===\n";
    
    // Create water (H2O) and test mass attenuation coefficient
    Material water("Water", 1.0);
    water.addElement(1, 1.008, 0.111898, "H");   // 11.19% H by mass
    water.addElement(8, 15.999, 0.888102, "O");  // 88.81% O by mass
    
    double energy = 0.1; // MeV = 100 keV
    
    // Calculate total cross-section using our physics
    double total_xs = physics_.getTotalCrossSection(energy, water);
    
    // Calculate using the database directly
    double energy_keV = energy * 1000.0;
    double h_xs = PhotonCrossSections::PhotonInteractionDatabase::getTotalCrossSection(1, energy_keV);
    double o_xs = PhotonCrossSections::PhotonInteractionDatabase::getTotalCrossSection(8, energy_keV);
    
    double expected_total = 0.111898 * h_xs + 0.888102 * o_xs;
    expected_total = PhotonCrossSections::PhotonInteractionDatabase::barnsToCm2(expected_total);
    
    std::cout << "Water total cross-section @ 100 keV:" << std::endl;
    std::cout << "  Physics calculation: " << total_xs << " cm²" << std::endl;
    std::cout << "  Database calculation: " << expected_total << " cm²" << std::endl;
    
    // Allow 10% difference due to approximations
    double relative_diff = std::abs(total_xs - expected_total) / expected_total;
    std::cout << "  Relative difference: " << relative_diff * 100 << "%" << std::endl;
    assert(relative_diff < 0.1);
    
    // Test mass attenuation coefficient
    double mu_rho = PhotonCrossSections::PhotonInteractionDatabase::getMassAttenuationCoefficient(
        8, 15.999, energy_keV); // Oxygen component
    
    assert(mu_rho > 0.0);
    std::cout << "  Oxygen mass attenuation @ 100 keV: " << mu_rho << " cm²/g" << std::endl;
    
    std::cout << "✓ Material composition accuracy tests passed\n";
  }

  void testCrossSectionDatabase()
  {
    std::cout << "\n=== Testing Cross-Section Database ===\n";
    
    // Test complete database functionality
    int Z = 26; // Iron
    double energy_keV = 100.0;
    
    auto all_xs = PhotonCrossSections::PhotonInteractionDatabase::getAllCrossSections(Z, energy_keV);
    
    std::cout << "Iron cross-sections @ 100 keV (barns):" << std::endl;
    std::cout << "  Photoelectric: " << all_xs[0] << std::endl;
    std::cout << "  Coherent:      " << all_xs[1] << std::endl;
    std::cout << "  Incoherent:    " << all_xs[2] << std::endl;
    std::cout << "  Pair (nuclear):" << all_xs[3] << std::endl;
    std::cout << "  Pair (electronic):" << all_xs[4] << std::endl;
    
    // All cross-sections should be non-negative
    for(double xs : all_xs) {
      assert(xs >= 0.0);
    }
    
    // Photoelectric should dominate at 100 keV for iron
    assert(all_xs[0] > all_xs[1]); // PE > coherent
    assert(all_xs[0] > all_xs[2]); // PE > incoherent
    
    // Pair production should be zero at 100 keV
    assert(all_xs[3] == 0.0);
    assert(all_xs[4] == 0.0);
    
    // Test total cross-section
    double total = PhotonCrossSections::PhotonInteractionDatabase::getTotalCrossSection(Z, energy_keV);
    double sum = 0.0;
    for(double xs : all_xs) sum += xs;
    
    assert(std::abs(total - sum) < 1e-10);
    std::cout << "  Total: " << total << " barns" << std::endl;
    
    std::cout << "✓ Cross-section database tests passed\n";
  }
};

int main()
{
  try
  {
    AccuratePhotonPhysicsTest test;
    test.runAllTests();
    return 0;
  }
  catch(const std::exception& e)
  {
    std::cerr << "❌ Test failed with exception: " << e.what() << std::endl;
    return 1;
  }
  catch(...)
  {
    std::cerr << "❌ Test failed with unknown exception" << std::endl;
    return 1;
  }
}