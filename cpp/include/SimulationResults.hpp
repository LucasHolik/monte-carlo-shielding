#pragma once

#include "Material.hpp"
#include "Particle.hpp"
#include "Vector3D.hpp"
#include "PhotonPhysics.hpp"

#include <array>
#include <map>
#include <string>
#include <vector>

/**
 * @brief Statistical tally for accumulating simulation results with uncertainty
 */
struct StatisticalTally 
{
  double sum = 0.0;           // Sum of values
  double sum_squared = 0.0;   // Sum of squared values  
  size_t count = 0;           // Number of entries
  double min_value = 1e30;    // Minimum observed value
  double max_value = -1e30;   // Maximum observed value
  
  // Add a value to the tally
  void add(double value) {
    sum += value;
    sum_squared += value * value;
    count++;
    if (value < min_value) min_value = value;
    if (value > max_value) max_value = value;
  }
  
  // Add another tally (for combining results)
  void add(const StatisticalTally& other) {
    sum += other.sum;
    sum_squared += other.sum_squared;
    count += other.count;
    if (other.min_value < min_value) min_value = other.min_value;
    if (other.max_value > max_value) max_value = other.max_value;
  }
  
  // Calculate mean
  double mean() const {
    return count > 0 ? sum / count : 0.0;
  }
  
  // Calculate standard deviation
  double standardDeviation() const {
    if (count <= 1) return 0.0;
    double mean_val = mean();
    return std::sqrt((sum_squared / count - mean_val * mean_val) * count / (count - 1));
  }
  
  // Calculate relative error (standard error / mean)
  double relativeError() const {
    double mean_val = mean();
    if (mean_val == 0.0 || count <= 1) return 0.0;
    return standardDeviation() / (std::sqrt(count) * std::abs(mean_val));
  }
  
  // Calculate standard error of mean
  double standardError() const {
    return count > 1 ? standardDeviation() / std::sqrt(count) : 0.0;
  }
  
  // Reset tally
  void clear() {
    sum = 0.0;
    sum_squared = 0.0;
    count = 0;
    min_value = 1e30;
    max_value = -1e30;
  }
};

/**
 * @brief Energy deposition data for a specific material/region
 */
struct EnergyDepositionData
{
  StatisticalTally total_energy_deposited;     // Total energy absorbed (MeV)
  StatisticalTally dose_rate;                  // Dose rate (Gy/s)
  StatisticalTally kerma_rate;                 // Kinetic energy released per unit mass (Gy/s)
  std::map<PhotonInteractionType, StatisticalTally> by_interaction_type; // Energy dep by interaction type
  std::vector<double> energy_spectrum;         // Energy deposition spectrum
  double material_mass = 0.0;                  // Mass of material region (g)
  
  void clear() {
    total_energy_deposited.clear();
    dose_rate.clear();
    kerma_rate.clear();
    by_interaction_type.clear();
    energy_spectrum.clear();
    material_mass = 0.0;
  }
};

/**
 * @brief Particle flux data for a specific region/detector
 */
struct FluxData
{
  StatisticalTally particle_count;             // Number of particles
  StatisticalTally fluence;                    // Particles per unit area (particles/cm²)
  StatisticalTally energy_flux;               // Energy flux (MeV/cm²)
  std::map<double, StatisticalTally> energy_spectrum; // Flux vs energy
  std::vector<double> angular_distribution;   // Angular flux distribution
  double detector_area = 0.0;                 // Detector area (cm²)
  
  void clear() {
    particle_count.clear();
    fluence.clear();
    energy_flux.clear();
    energy_spectrum.clear();
    angular_distribution.clear();
    detector_area = 0.0;
  }
};

/**
 * @brief Individual particle interaction record
 */
struct InteractionRecord
{
  Vector3D position;                           // Interaction position (cm)
  double energy_before = 0.0;                 // Energy before interaction (MeV)
  double energy_after = 0.0;                  // Energy after interaction (MeV) 
  double energy_deposited = 0.0;              // Energy deposited (MeV)
  PhotonInteractionType interaction_type;     // Type of interaction
  Material material;                           // Material where interaction occurred
  Vector3D direction_before;                   // Direction before interaction
  Vector3D direction_after;                    // Direction after interaction
  int generation = 0;                          // Particle generation (0=primary)
  bool particle_absorbed = false;              // True if particle was absorbed
  std::vector<Particle> secondary_particles;  // Any secondary particles created
};

/**
 * @brief Complete particle history from birth to termination
 */
struct ParticleHistory
{
  int history_id = 0;                          // Unique history identifier
  Particle initial_particle;                  // Initial particle state
  std::vector<InteractionRecord> interactions; // All interactions
  Vector3D final_position;                     // Final position
  double final_energy = 0.0;                  // Final energy (MeV)
  ParticleState final_state;                   // Final state (absorbed, escaped, etc.)
  double total_path_length = 0.0;             // Total distance traveled (cm)
  double total_energy_deposited = 0.0;        // Total energy deposited by this history (MeV)
  
  void clear() {
    history_id = 0;
    interactions.clear();
    final_position = Vector3D(0, 0, 0);
    final_energy = 0.0;
    final_state = ParticleState::Alive;
    total_path_length = 0.0;
    total_energy_deposited = 0.0;
  }
};

/**
 * @brief Simulation timing and performance data
 */
struct PerformanceData
{
  double simulation_time_seconds = 0.0;       // Total simulation time
  double particles_per_second = 0.0;          // Particle transport rate
  size_t total_interactions = 0;              // Total number of interactions
  size_t particles_absorbed = 0;              // Number absorbed
  size_t particles_escaped = 0;               // Number escaped
  size_t particles_terminated = 0;            // Number terminated (low energy, etc.)
  std::map<PhotonInteractionType, size_t> interactions_by_type; // Interaction counts
  
  void clear() {
    simulation_time_seconds = 0.0;
    particles_per_second = 0.0;
    total_interactions = 0;
    particles_absorbed = 0;
    particles_escaped = 0;
    particles_terminated = 0;
    interactions_by_type.clear();
  }
};

/**
 * @brief Complete simulation results with comprehensive data collection
 * 
 * Modern C++17 implementation providing professional-grade Monte Carlo
 * simulation results with statistical analysis and uncertainty quantification.
 */
class SimulationResults
{
private:
  std::string simulation_name_;
  size_t total_histories_requested_ = 0;
  size_t total_histories_completed_ = 0;
  
  // Energy deposition tallies by material
  std::map<std::string, EnergyDepositionData> energy_deposition_by_material_;
  
  // Flux data by detector/region
  std::map<std::string, FluxData> flux_by_region_;
  
  // Complete particle histories (optional - can be memory intensive)
  std::vector<ParticleHistory> complete_histories_;
  bool store_complete_histories_ = false;
  
  // Performance and timing data
  PerformanceData performance_;
  
  // Overall simulation statistics
  StatisticalTally total_energy_deposited_;    // Sum over all materials
  StatisticalTally transmission_coefficient_;   // For shielding calculations
  StatisticalTally reflection_coefficient_;     // For backscatter calculations
  
public:
  // Constructors
  SimulationResults() = default;
  explicit SimulationResults(const std::string& name) : simulation_name_(name) {}
  
  // Configuration
  void setSimulationName(const std::string& name) { simulation_name_ = name; }
  void setTotalHistoriesRequested(size_t n) { total_histories_requested_ = n; }
  void enableCompleteHistoryStorage(bool enable) { store_complete_histories_ = enable; }
  
  // Data recording methods
  void recordEnergyDeposition(const std::string& material_name, 
                            double energy_deposited, 
                            PhotonInteractionType interaction_type = PhotonInteractionType::PHOTOELECTRIC_ABSORPTION);
  
  void recordFlux(const std::string& region_name, 
                 const Particle& particle,
                 double detector_area = 1.0);
  
  void recordInteraction(const InteractionRecord& record);
  void recordParticleHistory(const ParticleHistory& history);
  void recordPerformanceData(const PerformanceData& perf_data);
  
  // Batch statistics methods (for proper uncertainty estimation)
  void startNewBatch();
  void completeBatch();
  
  // Accessors
  const std::string& getSimulationName() const { return simulation_name_; }
  size_t getTotalHistoriesRequested() const { return total_histories_requested_; }
  size_t getTotalHistoriesCompleted() const { return total_histories_completed_; }
  double getCompletionPercentage() const;
  
  // Energy deposition results
  const std::map<std::string, EnergyDepositionData>& getEnergyDepositionData() const {
    return energy_deposition_by_material_;
  }
  
  EnergyDepositionData getEnergyDepositionForMaterial(const std::string& material_name) const;
  double getTotalEnergyDeposited() const { return total_energy_deposited_.sum; }
  double getTotalEnergyDepositedUncertainty() const { return total_energy_deposited_.standardError(); }
  
  // Flux results  
  const std::map<std::string, FluxData>& getFluxData() const { return flux_by_region_; }
  FluxData getFluxForRegion(const std::string& region_name) const;
  
  // Shielding analysis results
  double getTransmissionCoefficient() const { return transmission_coefficient_.mean(); }
  double getTransmissionUncertainty() const { return transmission_coefficient_.standardError(); }
  double getReflectionCoefficient() const { return reflection_coefficient_.mean(); }
  double getReflectionUncertainty() const { return reflection_coefficient_.standardError(); }
  
  // Performance data
  const PerformanceData& getPerformanceData() const { return performance_; }
  double getSimulationTime() const { return performance_.simulation_time_seconds; }
  double getParticlesPerSecond() const { return performance_.particles_per_second; }
  
  // Complete histories (if enabled)
  const std::vector<ParticleHistory>& getCompleteHistories() const { return complete_histories_; }
  bool isCompleteHistoryStorageEnabled() const { return store_complete_histories_; }
  
  // Statistical analysis
  bool hasConverged(double target_relative_error = 0.05) const;
  std::vector<std::string> getConvergenceSummary(double target_relative_error = 0.05) const;
  
  // Data export methods
  void exportToCSV(const std::string& filename) const;
  void exportSummaryToJSON(const std::string& filename) const;
  std::string getSummaryReport() const;
  
  // Utility methods
  void clear();
  void merge(const SimulationResults& other);  // Combine results from multiple runs
  bool isValid() const;
  
  // Debugging and diagnostics
  std::string getDetailedSummary() const;
  void printPerformanceReport() const;
};

/**
 * @brief Helper functions for simulation analysis
 */
namespace SimulationAnalysis 
{
  /**
   * @brief Calculate dose rate from energy deposition
   * @param energy_deposited_MeV Energy deposited (MeV)
   * @param mass_g Mass of material (g) 
   * @param time_s Simulation time (s)
   * @return Dose rate (Gy/s)
   */
  double calculateDoseRate(double energy_deposited_MeV, double mass_g, double time_s);
  
  /**
   * @brief Calculate attenuation coefficient from transmission data
   * @param transmission_coefficient Measured transmission coefficient
   * @param thickness_cm Material thickness (cm)
   * @return Linear attenuation coefficient (cm⁻¹)
   */
  double calculateAttenuationCoefficient(double transmission_coefficient, double thickness_cm);
  
  /**
   * @brief Calculate buildup factor from simulation results
   * @param total_flux Total flux (uncollided + scattered)
   * @param uncollided_flux Uncollided flux only  
   * @return Buildup factor (dimensionless)
   */
  double calculateBuildupFactor(double total_flux, double uncollided_flux);
  
  /**
   * @brief Validate energy conservation in simulation
   * @param results Simulation results to validate
   * @return True if energy is conserved within tolerance
   */
  bool validateEnergyConservation(const SimulationResults& results, double tolerance = 1e-6);
  
  /**
   * @brief Compare simulation results with Beer-Lambert law
   * @param results Simulation results
   * @param material Material properties
   * @param thickness_cm Material thickness (cm)
   * @param energy_keV Photon energy (keV)
   * @return Relative difference from theoretical value
   */
  double compareToBeerLambert(const SimulationResults& results, 
                            const Material& material,
                            double thickness_cm, 
                            double energy_keV);
}