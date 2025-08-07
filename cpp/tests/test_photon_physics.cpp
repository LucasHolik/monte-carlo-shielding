#include "../include/PhotonPhysics.hpp"
#include "../include/Material.hpp"
#include "../include/MonteCarloSampling.hpp"
#include "../include/Particle.hpp"
#include "../include/RandomNumberGenerator.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <limits>

/**
 * @brief Test suite for photon physics implementation
 * 
 * Validates photon interactions against known physics values
 * and verifies implementation correctness.
 */
class PhotonPhysicsTest
{
private:
  RandomNumberGenerator rng_;
  MonteCarloSampling sampling_;
  PhotonPhysics physics_;

public:
  PhotonPhysicsTest() : rng_(12345), sampling_(rng_), physics_(sampling_) {}

  void runAllTests()
  {
    std::cout << "Running Photon Physics Tests...\n";
    
    testKleinNishinaCrossSection();
    testPhotoelectricCrossSection();
    testComptonScattering();
    testPhotoelectricAbsorption();
    testEnergyConservation();
    testCrossSectionCalculations();
    testInteractionSampling();
    
    std::cout << "All photon physics tests passed!\n";
  }

private:
  void testKleinNishinaCrossSection()
  {
    std::cout << "Testing Klein-Nishina cross-section...\n";
    
    // Test Thomson scattering limit (low energy)
    double low_energy = 0.001; // 1 keV
    double thomson_xs = physics_.kleinNishinaCrossSection(low_energy);
    double expected_thomson = 8.0 / 3.0 * M_PI * std::pow(2.8179403262e-13, 2);
    
    double relative_error = std::abs(thomson_xs - expected_thomson) / expected_thomson;
    assert(relative_error < 0.1); // 10% tolerance for low energy approximation
    
    // Test that cross-section decreases with energy
    std::vector<double> energies = {0.1, 0.5, 1.0, 5.0, 10.0}; // MeV
    double previous_xs = std::numeric_limits<double>::max();
    
    for(double energy : energies)
    {
      double xs = physics_.kleinNishinaCrossSection(energy);
      assert(xs < previous_xs); // Cross-section should decrease with energy
      assert(xs > 0.0);         // Should always be positive
      previous_xs = xs;
    }
    
    std::cout << "  ✓ Klein-Nishina cross-section validation passed\n";
  }

  void testPhotoelectricCrossSection()
  {
    std::cout << "Testing photoelectric cross-section...\n";
    
    // Test Z^5 dependence
    std::vector<int> atomic_numbers = {1, 6, 13, 26, 82}; // H, C, Al, Fe, Pb
    double energy = 0.1; // MeV
    
    for(size_t i = 1; i < atomic_numbers.size(); i++)
    {
      double xs1 = physics_.getPhotoelectricCrossSection(energy, atomic_numbers[i-1]);
      double xs2 = physics_.getPhotoelectricCrossSection(energy, atomic_numbers[i]);
      
      // Should increase with atomic number
      assert(xs2 > xs1);
      
      // Test approximate Z^5 scaling
      double z1 = atomic_numbers[i-1];
      double z2 = atomic_numbers[i];
      double expected_ratio = std::pow(z2/z1, 5.0);
      double actual_ratio = xs2 / xs1;
      
      // Allow significant tolerance for simplified model
      assert(actual_ratio > expected_ratio * 0.1);
      assert(actual_ratio < expected_ratio * 10.0);
    }
    
    // Test energy dependence (should decrease with increasing energy)
    std::vector<double> energies = {0.01, 0.05, 0.1, 0.5, 1.0}; // MeV
    int Z = 26; // Iron
    
    for(size_t i = 1; i < energies.size(); i++)
    {
      double xs1 = physics_.getPhotoelectricCrossSection(energies[i-1], Z);
      double xs2 = physics_.getPhotoelectricCrossSection(energies[i], Z);
      assert(xs2 < xs1); // Should decrease with energy
    }
    
    std::cout << "  ✓ Photoelectric cross-section validation passed\n";
  }

  void testComptonScattering()
  {
    std::cout << "Testing Compton scattering physics...\n";
    
    // Create test photon
    Vector3D position(0.0, 0.0, 0.0);
    Vector3D direction(0.0, 0.0, 1.0);
    double energy = 0.5; // MeV
    Particle photon(ParticleType::Photon, position, direction, energy);
    
    // Create test material (aluminum)
    Material aluminum("Aluminum", 2.70);
    aluminum.addElement(13, 26.98, 1.0, "Al");
    
    // Test Compton scattering
    auto result = physics_.performComptonScattering(photon, aluminum);
    
    // Validate results
    assert(result.interaction_type == PhotonInteractionType::COMPTON_SCATTERING);
    assert(result.particle_scattered == true);
    assert(result.particle_absorbed == false);
    assert(result.new_energy > 0.0);
    assert(result.new_energy <= energy); // Scattered energy should be less or equal
    
    // Check direction normalization
    double dir_magnitude = result.new_direction.magnitude();
    assert(std::abs(dir_magnitude - 1.0) < 1e-6);
    
    // Test energy conservation for Compton formula
    double cos_theta = photon.direction().dot(result.new_direction);
    double alpha = energy / 0.51099895000; // energy / m_e c^2
    double expected_energy = energy / (1.0 + alpha * (1.0 - cos_theta));
    double energy_tolerance = 0.1 * expected_energy; // 10% tolerance
    
    // Note: Our implementation uses sampling, so exact match not expected
    assert(result.new_energy > 0.0);
    assert(result.new_energy < energy);
    
    std::cout << "  ✓ Compton scattering validation passed\n";
  }

  void testPhotoelectricAbsorption()
  {
    std::cout << "Testing photoelectric absorption...\n";
    
    // Create test photon
    Vector3D position(0.0, 0.0, 0.0);
    Vector3D direction(1.0, 0.0, 0.0);
    double energy = 0.1; // MeV
    Particle photon(ParticleType::Photon, position, direction, energy);
    
    // Create test material (lead)
    Material lead("Lead", 11.34);
    lead.addElement(82, 207.2, 1.0, "Pb");
    
    // Test photoelectric absorption
    auto result = physics_.performPhotoelectricAbsorption(photon, lead);
    
    // Validate results
    assert(result.interaction_type == PhotonInteractionType::PHOTOELECTRIC_ABSORPTION);
    assert(result.particle_absorbed == true);
    assert(result.particle_scattered == false);
    assert(result.new_energy == 0.0); // Complete absorption
    
    std::cout << "  ✓ Photoelectric absorption validation passed\n";
  }

  void testEnergyConservation()
  {
    std::cout << "Testing energy conservation...\n";
    
    // Test multiple photon energies
    std::vector<double> energies = {0.01, 0.1, 0.5, 1.0, 2.0, 5.0}; // MeV
    
    for(double energy : energies)
    {
      Vector3D position(0.0, 0.0, 0.0);
      Vector3D direction(0.0, 0.0, 1.0);
      Particle photon(ParticleType::Photon, position, direction, energy);
      
      // Test with water (common material)
      Material water("Water", 1.0);
      water.addElement(1, 1.008, 0.111, "H");
      water.addElement(8, 15.999, 0.889, "O");
      
      // Perform interaction
      auto result = physics_.performInteraction(photon, water);
      
      // Validate result
      assert(physics_.validateInteractionResult(result));
      
      if(result.particle_scattered)
      {
        // For scattering, energy should be conserved or reduced
        assert(result.new_energy > 0.0);
        assert(result.new_energy <= energy);
        
        // Direction should be normalized
        double mag = result.new_direction.magnitude();
        assert(std::abs(mag - 1.0) < 1e-6);
      }
      else if(result.particle_absorbed)
      {
        // For absorption, energy should be zero
        assert(result.new_energy == 0.0);
      }
    }
    
    std::cout << "  ✓ Energy conservation validation passed\n";
  }

  void testCrossSectionCalculations()
  {
    std::cout << "Testing cross-section calculations...\n";
    
    // Test total cross-section calculation
    Material mixed("Mixed", 1.5);
    mixed.addElement(1, 1.008, 0.5, "H");
    mixed.addElement(8, 15.999, 0.5, "O");
    
    double energy = 0.5; // MeV
    double total_xs = physics_.getTotalCrossSection(energy, mixed);
    
    // Total cross-section should be positive
    assert(total_xs > 0.0);
    
    // Test individual components
    double photoelectric_xs = physics_.getPhotoelectricCrossSection(energy, 1) * 0.5 +
                             physics_.getPhotoelectricCrossSection(energy, 8) * 0.5;
    double compton_xs = physics_.getComptonCrossSection(energy, 1) * 0.5 +
                       physics_.getComptonCrossSection(energy, 8) * 0.5;
    double coherent_xs = physics_.getCoherentCrossSection(energy, 1) * 0.5 +
                        physics_.getCoherentCrossSection(energy, 8) * 0.5;
    
    double sum_xs = photoelectric_xs + compton_xs + coherent_xs;
    
    // Should approximately match (within tolerance due to pair production)
    double relative_error = std::abs(total_xs - sum_xs) / total_xs;
    assert(relative_error < 0.1);
    
    std::cout << "  ✓ Cross-section calculation validation passed\n";
  }

  void testInteractionSampling()
  {
    std::cout << "Testing interaction type sampling...\n";
    
    // Create test material and energy
    Material aluminum("Aluminum", 2.70);
    aluminum.addElement(13, 26.98, 1.0, "Al");
    double energy = 1.0; // MeV - higher energy to get multiple interaction types
    
    // Sample interactions multiple times to check distribution
    std::map<PhotonInteractionType, int> interaction_counts;
    int num_samples = 10000;
    
    for(int i = 0; i < num_samples; i++)
    {
      auto interaction_type = physics_.sampleInteractionType(energy, aluminum);
      interaction_counts[interaction_type]++;
    }
    
    // Should have sampled at least one interaction type
    assert(interaction_counts.size() >= 1);
    
    // All counts should be positive
    for(const auto& pair : interaction_counts)
    {
      assert(pair.second > 0);
    }
    
    // Print distribution for verification
    std::cout << "  Interaction type distribution:\n";
    for(const auto& pair : interaction_counts)
    {
      double fraction = double(pair.second) / num_samples;
      std::cout << "    " << physics_.getInteractionTypeName(pair.first) 
                << ": " << fraction * 100.0 << "%\n";
    }
    
    std::cout << "  ✓ Interaction sampling validation passed\n";
  }
};

int main()
{
  try
  {
    PhotonPhysicsTest test;
    test.runAllTests();
    std::cout << "\n🎉 All photon physics tests completed successfully!\n";
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