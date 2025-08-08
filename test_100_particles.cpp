#include "cpp/include/Box.hpp"
#include "cpp/include/Material.hpp"
#include "cpp/include/MonteCarloSampling.hpp"
#include "cpp/include/RandomNumberGenerator.hpp"
#include "cpp/include/Transport.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <chrono>

int main()
{
  std::cout << "Monte Carlo Shielding Analysis - Proper Air-Shield-Air Geometry\n";
  std::cout << "================================================================\n\n";

  // Setup
  const size_t num_runs = 10;
  const size_t particles_per_run = 10000;
  const double energy_MeV = 0.5;   // 500 keV
  const double thickness_cm = 1.0; // 1 cm lead

  std::cout << "Configuration:\n";
  std::cout << "  Energy: " << energy_MeV << " MeV\n";
  std::cout << "  Shield: " << thickness_cm << " cm lead\n";
  std::cout << "  Geometry: Air(-5 to 0) | Shield(0 to " << thickness_cm << ") | Air(" << thickness_cm << " to 5)\n";
  std::cout << "  Start position: z=-2.0 cm (in air before shield)\n";
  std::cout << "  Particles per run: " << particles_per_run << "\n";
  std::cout << "  Number of runs: " << num_runs << "\n\n";

  // Geometry setup: Proper air-shield-air configuration
  Material lead = Material::createLead();
  Geometry geometry = Geometry::createShieldingExperiment(
    thickness_cm,  // Shield thickness
    lead,         // Shield material  
    10.0,         // Total experiment length (z=-5 to z=5)
    20.0          // Cross-section (20x20 cm)
  );

  // Results storage for statistical analysis
  std::vector<double> transmission_rates;
  std::vector<size_t> transmitted_counts;
  std::vector<double> energy_deposited_per_run;

  auto start_time = std::chrono::steady_clock::now();

  std::cout << "Results:\n";

  // Run 10 independent simulations
  for(size_t run = 0; run < num_runs; ++run)
  {
    // Create new RNG with different seed for each run
    RandomNumberGenerator rng(12345 + run * 1000);
    MonteCarloSampling sampling(rng);
    Transport transport(geometry, sampling);

    // Statistics for this run
    size_t transmitted = 0;
    double total_energy_deposited = 0.0;

    // Transport particles for this run
    for(size_t i = 0; i < particles_per_run; ++i)
    {      
      // Start all particles in air region before shield
      Vector3D start_pos(0, 0, -2.0);  // In air, 2 cm before shield
      Vector3D direction(0, 0, 1);      // Moving toward shield
      
      Particle photon(ParticleType::Photon,
                      start_pos,
                      direction,
                      energy_MeV);

      double initial_energy = photon.energy();
      
      transport.trackParticle(photon);

      Vector3D final_pos = photon.position();
      double final_energy = photon.energy();
      
      // Energy deposited
      double energy_deposited = initial_energy - final_energy;
      total_energy_deposited += energy_deposited;

      // Classify particle fate
      // Particles transmitted if they reach the air region after the shield (z > thickness_cm)
      // In our geometry: air(-5 to 0), shield(0 to 1), air(1 to 5)
      if(final_pos.z() > thickness_cm) {
        transmitted++;  // Successfully transmitted through shield
      }
    }

    // Calculate results for this run
    double transmission_rate = static_cast<double>(transmitted) / particles_per_run * 100.0;
    
    // Store results
    transmission_rates.push_back(transmission_rate);
    transmitted_counts.push_back(transmitted);
    energy_deposited_per_run.push_back(total_energy_deposited);

    // Display run results
    std::cout << "Run " << std::setw(2) << (run + 1) << ": " 
              << std::fixed << std::setprecision(2) << transmission_rate << "% transmission ("
              << transmitted << "/" << particles_per_run << "), Energy deposited: "
              << std::setprecision(1) << total_energy_deposited << " MeV\n";
  }

  auto end_time = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

  // Calculate summary statistics across all runs
  double mean_transmission = 0.0;
  double mean_energy_deposited = 0.0;
  size_t total_transmitted = 0;
  
  for(size_t i = 0; i < num_runs; ++i) {
    mean_transmission += transmission_rates[i];
    mean_energy_deposited += energy_deposited_per_run[i];
    total_transmitted += transmitted_counts[i];
  }
  mean_transmission /= num_runs;
  mean_energy_deposited /= num_runs;
  
  // Calculate standard deviation
  double std_dev_transmission = 0.0;
  if(num_runs > 1) {
    for(double rate : transmission_rates) {
      double diff = rate - mean_transmission;
      std_dev_transmission += diff * diff;
    }
    std_dev_transmission = std::sqrt(std_dev_transmission / (num_runs - 1));
  }
  
  // Statistical uncertainty (standard error of mean)
  double statistical_uncertainty = (num_runs > 1) ? std_dev_transmission / std::sqrt(num_runs) : 0.0;

  // Theoretical value
  double theoretical = std::exp(-lead.getLinearAttenuationCoefficient(energy_MeV * 1000.0) * thickness_cm) * 100.0;

  // Display summary statistics
  std::cout << "\n========== SUMMARY STATISTICS ==========\n";
  std::cout << "Total particles: " << (num_runs * particles_per_run) << " (" << num_runs << " runs of " << particles_per_run << " each)\n";
  std::cout << "Simulation time: " << duration.count() << " seconds\n\n";
  
  std::cout << "Transmission Analysis:\n";
  std::cout << "  Mean transmission:     " << std::fixed << std::setprecision(3) << mean_transmission << " ± " << statistical_uncertainty << "%\n";
  std::cout << "  Standard deviation:    " << std::setprecision(3) << std_dev_transmission << "%\n";
  std::cout << "  Total transmitted:     " << total_transmitted << "/" << (num_runs * particles_per_run) << "\n";
  std::cout << "  Theoretical (Beer-L):  " << std::setprecision(3) << theoretical << "%\n";
  std::cout << "  Difference from theory:" << std::setprecision(3) << std::abs(mean_transmission - theoretical) << " percentage points\n\n";
  
  std::cout << "Energy Analysis:\n";
  std::cout << "  Mean energy deposited per run: " << std::setprecision(1) << mean_energy_deposited << " MeV\n";
  std::cout << "  Total energy deposited: " << std::setprecision(1) << (mean_energy_deposited * num_runs) << " MeV\n\n";
  
  std::cout << "Statistical Quality:\n";
  std::cout << "  Relative uncertainty:  " << std::setprecision(2) << (statistical_uncertainty / mean_transmission * 100.0) << "%\n";
  std::cout << "  Performance:           " << std::setprecision(0) << ((num_runs * particles_per_run) / duration.count()) << " particles/second\n";
  
  std::cout << "========================================\n";

  return 0;
}