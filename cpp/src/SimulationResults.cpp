#include "SimulationResults.hpp"
#include "PhotonPhysics.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

// =============================================================================
// SimulationResults Implementation
// =============================================================================

double SimulationResults::getCompletionPercentage() const
{
  if(total_histories_requested_ == 0)
    return 0.0;
  return 100.0 * total_histories_completed_ / total_histories_requested_;
}

void SimulationResults::recordEnergyDeposition(
    const std::string &material_name, double energy_deposited,
    PhotonInteractionType interaction_type)
{
  if(energy_deposited <= 0.0)
    return;

  // Record in material-specific tally
  auto &material_data = energy_deposition_by_material_[material_name];
  material_data.total_energy_deposited.add(energy_deposited);
  material_data.by_interaction_type[interaction_type].add(energy_deposited);

  // Record in overall tally
  total_energy_deposited_.add(energy_deposited);
}

void SimulationResults::recordFlux(const std::string &region_name,
                                   const Particle &particle,
                                   double detector_area)
{
  auto &flux_data = flux_by_region_[region_name];

  // Record particle count
  flux_data.particle_count.add(1.0);

  // Record fluence (particles per unit area)
  if(detector_area > 0.0)
  {
    flux_data.fluence.add(1.0 / detector_area);
    flux_data.detector_area = detector_area;
  }

  // Record energy flux
  double energy_MeV = particle.energy();
  flux_data.energy_flux.add(energy_MeV / detector_area);

  // Record energy spectrum (binned by energy)
  double energy_bin = std::floor(energy_MeV * 10.0) / 10.0; // 0.1 MeV bins
  flux_data.energy_spectrum[energy_bin].add(1.0);
}

void SimulationResults::recordInteraction(const InteractionRecord &record)
{
  // Record energy deposition from this interaction
  if(record.energy_deposited > 0.0)
  {
    recordEnergyDeposition(std::string(record.material.name()),
                           record.energy_deposited, record.interaction_type);
  }

  // Update performance counters
  performance_.total_interactions++;
  performance_.interactions_by_type[record.interaction_type]++;
}

void SimulationResults::recordParticleHistory(const ParticleHistory &history)
{
  total_histories_completed_++;

  // Record final particle state
  switch(history.final_state)
  {
  case ParticleState::Absorbed:
    performance_.particles_absorbed++;
    break;
  case ParticleState::Escaped:
    performance_.particles_escaped++;
    break;
  case ParticleState::Terminated:
    performance_.particles_terminated++;
    break;
  default:
    break;
  }

  // Store complete history if enabled
  if(store_complete_histories_)
  {
    complete_histories_.push_back(history);
  }
}

void SimulationResults::recordPerformanceData(const PerformanceData &perf_data)
{
  performance_ = perf_data;

  // Calculate derived performance metrics
  if(performance_.simulation_time_seconds > 0.0)
  {
    performance_.particles_per_second =
        total_histories_completed_ / performance_.simulation_time_seconds;
  }
}

EnergyDepositionData SimulationResults::getEnergyDepositionForMaterial(
    const std::string &material_name) const
{
  auto it = energy_deposition_by_material_.find(material_name);
  return (it != energy_deposition_by_material_.end()) ? it->second
                                                      : EnergyDepositionData{};
}

FluxData
SimulationResults::getFluxForRegion(const std::string &region_name) const
{
  auto it = flux_by_region_.find(region_name);
  return (it != flux_by_region_.end()) ? it->second : FluxData{};
}

bool SimulationResults::hasConverged(double target_relative_error) const
{
  // Check convergence of key quantities
  bool total_energy_converged =
      total_energy_deposited_.relativeError() < target_relative_error;
  bool transmission_converged = true;
  bool reflection_converged = true;

  if(transmission_coefficient_.count > 0)
  {
    transmission_converged =
        transmission_coefficient_.relativeError() < target_relative_error;
  }

  if(reflection_coefficient_.count > 0)
  {
    reflection_converged =
        reflection_coefficient_.relativeError() < target_relative_error;
  }

  return total_energy_converged && transmission_converged &&
         reflection_converged;
}

std::vector<std::string>
SimulationResults::getConvergenceSummary(double target_relative_error) const
{
  std::vector<std::string> summary;

  double total_energy_error = total_energy_deposited_.relativeError();
  summary.push_back(
      "Total Energy Deposition: " + std::to_string(total_energy_error * 100.0) +
      "% " +
      (total_energy_error < target_relative_error ? "(CONVERGED)"
                                                  : "(NOT CONVERGED)"));

  if(transmission_coefficient_.count > 0)
  {
    double trans_error = transmission_coefficient_.relativeError();
    summary.push_back(
        "Transmission Coefficient: " + std::to_string(trans_error * 100.0) +
        "% " +
        (trans_error < target_relative_error ? "(CONVERGED)"
                                             : "(NOT CONVERGED)"));
  }

  if(reflection_coefficient_.count > 0)
  {
    double refl_error = reflection_coefficient_.relativeError();
    summary.push_back(
        "Reflection Coefficient: " + std::to_string(refl_error * 100.0) + "% " +
        (refl_error < target_relative_error ? "(CONVERGED)"
                                            : "(NOT CONVERGED)"));
  }

  return summary;
}

std::string SimulationResults::getSummaryReport() const
{
  std::stringstream ss;
  ss << std::fixed << std::setprecision(4);

  ss << "========================================\n";
  ss << "     MONTE CARLO SIMULATION RESULTS    \n";
  ss << "========================================\n\n";

  ss << "Simulation: " << simulation_name_ << "\n";
  ss << "Histories Completed: " << total_histories_completed_ << " / "
     << total_histories_requested_ << " (" << getCompletionPercentage()
     << "%)\n\n";

  // Energy Deposition Summary
  ss << "ENERGY DEPOSITION SUMMARY:\n";
  ss << "  Total Energy Deposited: " << total_energy_deposited_.mean() << " ± "
     << total_energy_deposited_.standardError() << " MeV\n";
  ss << "  Relative Error: " << total_energy_deposited_.relativeError() * 100.0
     << "%\n\n";

  // Material-specific energy deposition
  if(!energy_deposition_by_material_.empty())
  {
    ss << "ENERGY DEPOSITION BY MATERIAL:\n";
    for(const auto &[material_name, data] : energy_deposition_by_material_)
    {
      ss << "  " << material_name << ": " << data.total_energy_deposited.mean()
         << " ± " << data.total_energy_deposited.standardError() << " MeV\n";
    }
    ss << "\n";
  }

  // Flux Summary
  if(!flux_by_region_.empty())
  {
    ss << "FLUX SUMMARY:\n";
    for(const auto &[region_name, data] : flux_by_region_)
    {
      ss << "  " << region_name << " - Particles: " << data.particle_count.sum
         << ", Fluence: " << data.fluence.mean() << " ± "
         << data.fluence.standardError() << " particles/cm²\n";
    }
    ss << "\n";
  }

  // Shielding Analysis
  if(transmission_coefficient_.count > 0)
  {
    ss << "SHIELDING ANALYSIS:\n";
    ss << "  Transmission Coefficient: " << transmission_coefficient_.mean()
       << " ± " << transmission_coefficient_.standardError() << "\n";
    ss << "  Relative Error: "
       << transmission_coefficient_.relativeError() * 100.0 << "%\n";

    if(reflection_coefficient_.count > 0)
    {
      ss << "  Reflection Coefficient: " << reflection_coefficient_.mean()
         << " ± " << reflection_coefficient_.standardError() << "\n";
    }
    ss << "\n";
  }

  // Performance Summary
  ss << "PERFORMANCE SUMMARY:\n";
  ss << "  Simulation Time: " << performance_.simulation_time_seconds
     << " seconds\n";
  ss << "  Particles/Second: " << performance_.particles_per_second << "\n";
  ss << "  Total Interactions: " << performance_.total_interactions << "\n";
  ss << "  Particles Absorbed: " << performance_.particles_absorbed << "\n";
  ss << "  Particles Escaped: " << performance_.particles_escaped << "\n";
  ss << "  Particles Terminated: " << performance_.particles_terminated
     << "\n\n";

  // Convergence Status
  ss << "CONVERGENCE STATUS (5% target):\n";
  auto convergence_summary = getConvergenceSummary(0.05);
  for(const auto &line : convergence_summary)
  {
    ss << "  " << line << "\n";
  }

  return ss.str();
}

void SimulationResults::clear()
{
  simulation_name_.clear();
  total_histories_requested_ = 0;
  total_histories_completed_ = 0;
  energy_deposition_by_material_.clear();
  flux_by_region_.clear();
  complete_histories_.clear();
  performance_.clear();
  total_energy_deposited_.clear();
  transmission_coefficient_.clear();
  reflection_coefficient_.clear();
}

void SimulationResults::merge(const SimulationResults &other)
{
  // Merge basic counters
  total_histories_completed_ += other.total_histories_completed_;

  // Merge energy deposition data
  for(const auto &[material_name, other_data] :
      other.energy_deposition_by_material_)
  {
    auto &my_data = energy_deposition_by_material_[material_name];
    my_data.total_energy_deposited.add(other_data.total_energy_deposited);
    for(const auto &[interaction_type, other_tally] :
        other_data.by_interaction_type)
    {
      my_data.by_interaction_type[interaction_type].add(other_tally);
    }
  }

  // Merge flux data
  for(const auto &[region_name, other_flux] : other.flux_by_region_)
  {
    auto &my_flux = flux_by_region_[region_name];
    my_flux.particle_count.add(other_flux.particle_count);
    my_flux.fluence.add(other_flux.fluence);
    my_flux.energy_flux.add(other_flux.energy_flux);
  }

  // Merge overall tallies
  total_energy_deposited_.add(other.total_energy_deposited_);
  transmission_coefficient_.add(other.transmission_coefficient_);
  reflection_coefficient_.add(other.reflection_coefficient_);

  // Merge performance data
  performance_.total_interactions += other.performance_.total_interactions;
  performance_.particles_absorbed += other.performance_.particles_absorbed;
  performance_.particles_escaped += other.performance_.particles_escaped;
  performance_.particles_terminated += other.performance_.particles_terminated;

  for(const auto &[interaction_type, count] :
      other.performance_.interactions_by_type)
  {
    performance_.interactions_by_type[interaction_type] += count;
  }
}

bool SimulationResults::isValid() const
{
  // Check basic consistency
  if(total_histories_completed_ > total_histories_requested_)
  {
    return false;
  }

  // Check that particle counts add up
  size_t total_particles = performance_.particles_absorbed +
                           performance_.particles_escaped +
                           performance_.particles_terminated;

  if(total_particles > total_histories_completed_)
  {
    return false;
  }

  // Check for reasonable values
  if(total_energy_deposited_.sum < 0.0)
  {
    return false;
  }

  return true;
}

void SimulationResults::exportToCSV(const std::string &filename) const
{
  std::ofstream file(filename);
  file << std::fixed << std::setprecision(6);

  file << "# Monte Carlo Simulation Results - " << simulation_name_ << "\n";
  file << "# Histories Completed: " << total_histories_completed_ << "\n\n";

  // Energy deposition data
  file << "Material,Total_Energy_Deposited_MeV,Energy_Uncertainty_MeV,Relative_"
          "Error_Percent\n";
  for(const auto &[material_name, data] : energy_deposition_by_material_)
  {
    file << material_name << "," << data.total_energy_deposited.mean() << ","
         << data.total_energy_deposited.standardError() << ","
         << data.total_energy_deposited.relativeError() * 100.0 << "\n";
  }

  file << "\n# Flux Data\n";
  file << "Region,Particle_Count,Fluence_per_cm2,Fluence_Uncertainty,Energy_"
          "Flux_MeV_per_cm2\n";
  for(const auto &[region_name, data] : flux_by_region_)
  {
    file << region_name << "," << data.particle_count.sum << ","
         << data.fluence.mean() << "," << data.fluence.standardError() << ","
         << data.energy_flux.mean() << "\n";
  }

  file.close();
}

// =============================================================================
// SimulationAnalysis Helper Functions
// =============================================================================

namespace SimulationAnalysis
{

double calculateDoseRate(double energy_deposited_MeV, double mass_g,
                         double time_s)
{
  if(mass_g <= 0.0 || time_s <= 0.0)
    return 0.0;

  // Convert MeV to Joules: 1 MeV = 1.602176634e-13 J
  const double MeV_to_Joules = 1.602176634e-13;
  double energy_J = energy_deposited_MeV * MeV_to_Joules;

  // Convert mass to kg
  double mass_kg = mass_g / 1000.0;

  // Dose rate in Gy/s (1 Gy = 1 J/kg)
  return energy_J / (mass_kg * time_s);
}

double calculateAttenuationCoefficient(double transmission_coefficient,
                                       double thickness_cm)
{
  if(transmission_coefficient <= 0.0 || transmission_coefficient > 1.0 ||
     thickness_cm <= 0.0)
  {
    return 0.0;
  }

  return -std::log(transmission_coefficient) / thickness_cm;
}

double calculateBuildupFactor(double total_flux, double uncollided_flux)
{
  if(uncollided_flux <= 0.0)
    return 1.0;
  return total_flux / uncollided_flux;
}

bool validateEnergyConservation(const SimulationResults &results,
                                double tolerance)
{
  // This would need access to initial particle energies to fully validate
  // For now, check that energy deposition is reasonable
  double total_deposited = results.getTotalEnergyDeposited();

  // Basic sanity checks
  if(total_deposited < 0.0)
    return false;
  if(std::isnan(total_deposited) || std::isinf(total_deposited))
    return false;

  return true;
}

double compareToBeerLambert(const SimulationResults &results,
                            const Material &material, double thickness_cm,
                            double energy_keV)
{
  double simulated_transmission = results.getTransmissionCoefficient();

  // Calculate theoretical transmission using Beer-Lambert law
  double linear_mu = material.getLinearAttenuationCoefficient(energy_keV);
  double theoretical_transmission = std::exp(-linear_mu * thickness_cm);

  if(theoretical_transmission > 0.0)
  {
    return std::abs(simulated_transmission - theoretical_transmission) /
           theoretical_transmission;
  }

  return 0.0;
}

} // namespace SimulationAnalysis