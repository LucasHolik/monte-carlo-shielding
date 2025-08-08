#pragma once

#include "Geometry.hpp"
#include "MonteCarloSampling.hpp"
#include "ParticleSource.hpp"
#include "RandomNumberGenerator.hpp"
#include "SimulationResults.hpp"
#include "Transport.hpp"

#include <chrono>
#include <functional>
#include <memory>
#include <string>

/**
 * @brief Simulation configuration parameters
 */
struct SimulationConfiguration
{
  std::string name = "Photon Monte Carlo Simulation";
  size_t number_of_histories = 10000;           // Total particle histories to simulate
  bool enable_detailed_logging = false;         // Enable detailed interaction logging
  bool enable_variance_reduction = false;       // Enable variance reduction techniques
  double energy_cutoff_keV = 1.0;              // Energy below which particles are killed
  double time_cutoff_seconds = 3600.0;         // Maximum simulation time
  size_t batch_size = 1000;                    // Number of histories per batch for statistics
  double target_statistical_uncertainty = 0.05; // Target relative error (5%)
  
  // Progress reporting
  std::function<void(size_t, size_t, double)> progress_callback; // callback(completed, total, time_elapsed)
  size_t progress_reporting_frequency = 1000;   // Report progress every N histories
  
  // Validation and debugging
  bool validate_energy_conservation = true;     // Check energy conservation
  bool enable_particle_tracking = false;        // Store complete particle tracks
  
  // Parallel processing (for future extension)
  size_t number_of_threads = 1;                // Number of parallel threads
  
  bool isValid() const;
  std::string toString() const;
};

/**
 * @brief Complete photon Monte Carlo simulation engine
 * 
 * Professional implementation integrating all components:
 * - ParticleSource for initial particle generation
 * - Transport for particle movement and interactions  
 * - PhotonPhysics for realistic interaction sampling
 * - SimulationResults for comprehensive data collection
 * - Statistical analysis with uncertainty quantification
 */
class PhotonSimulation
{
private:
  std::shared_ptr<RandomNumberGenerator> rng_;
  std::shared_ptr<MonteCarloSampling> sampling_;
  std::shared_ptr<ParticleSource> source_;
  std::shared_ptr<Transport> transport_;
  std::shared_ptr<SimulationResults> results_;
  
  SimulationConfiguration config_;
  
  // Internal state
  bool is_initialized_ = false;
  bool is_running_ = false;
  bool is_paused_ = false;
  std::chrono::steady_clock::time_point start_time_;
  
  // Batch processing for proper statistics
  size_t current_batch_ = 0;
  size_t histories_completed_ = 0;
  std::vector<double> batch_results_;  // For batch-wise statistics
  
public:
  // Constructors
  PhotonSimulation();
  explicit PhotonSimulation(const SimulationConfiguration& config);
  
  // Setup and configuration
  void setConfiguration(const SimulationConfiguration& config);
  void setSource(std::shared_ptr<ParticleSource> source);
  void setGeometry(std::shared_ptr<Geometry> geometry);
  void setRandomSeed(uint64_t seed);
  
  // Main simulation methods
  std::shared_ptr<SimulationResults> runSimulation();
  std::shared_ptr<SimulationResults> runSimulation(const SourceConfiguration& source_config,
                                                   std::shared_ptr<Geometry> geometry);
  
  // Batch processing methods
  std::shared_ptr<SimulationResults> runBatchedSimulation();
  void runSingleBatch(size_t batch_size);
  bool hasConverged(double target_uncertainty = 0.05) const;
  
  // Control methods
  void pauseSimulation();
  void resumeSimulation();
  void stopSimulation();
  bool isRunning() const { return is_running_; }
  bool isPaused() const { return is_paused_; }
  
  // Progress and status
  double getCompletionPercentage() const;
  double getElapsedTime() const;
  double getEstimatedRemainingTime() const;
  std::string getStatusSummary() const;
  
  // Results access
  std::shared_ptr<SimulationResults> getResults() const { return results_; }
  
  // Validation methods
  bool validateConfiguration() const;
  bool validateEnergyConservation() const;
  std::vector<std::string> runDiagnostics() const;
  
  // Advanced features
  void enableVarianceReduction(bool enable);
  void setProgressCallback(std::function<void(size_t, size_t, double)> callback);
  
  // Factory methods for common simulation types
  static PhotonSimulation createShieldingAnalysis(
    const Vector3D& source_position, double source_energy_MeV,
    std::shared_ptr<Geometry> shielding_geometry,
    size_t num_histories = 100000);
    
  static PhotonSimulation createDoseCalculation(
    const SourceConfiguration& source_config,
    std::shared_ptr<Geometry> phantom_geometry,
    size_t num_histories = 1000000);
    
  static PhotonSimulation createBuildupFactorAnalysis(
    const Vector3D& source_position, double energy_MeV,
    const Material& shield_material, double thickness_cm,
    size_t num_histories = 100000);
    
  static PhotonSimulation createSpectrumAnalysis(
    const SourceConfiguration& source_config,
    std::shared_ptr<Geometry> detector_geometry,
    size_t num_histories = 500000);

private:
  // Internal simulation methods
  void initialize();
  void validateSetup();
  ParticleHistory transportSingleParticle(const SourceParticle& source_particle);
  void processParticleHistory(const ParticleHistory& history);
  void updateProgress(size_t histories_completed);
  void finalizeBatch();
  
  // Statistics and analysis helpers
  void calculateBatchStatistics();
  void updateConvergenceMetrics();
  bool checkTerminationCriteria();
  
  // Validation helpers
  bool validateParticleHistory(const ParticleHistory& history) const;
  double calculateTotalEnergyBalance() const;
  
  // Performance monitoring
  void recordPerformanceMetrics();
  double calculateParticlesPerSecond() const;
};

/**
 * @brief Convenience function for simple photon transport simulations
 * 
 * This function provides a simple interface for common simulation scenarios
 * without requiring detailed setup of all components.
 */
std::shared_ptr<SimulationResults> runPhotonSimulation(
  const SourceConfiguration& source_config,
  std::shared_ptr<Geometry> geometry,
  size_t num_histories = 10000,
  std::function<void(size_t, size_t)> progress_callback = nullptr);

/**
 * @brief Simplified shielding analysis function
 * 
 * Calculates transmission coefficient through a material shield
 */
struct ShieldingAnalysisResult
{
  double transmission_coefficient = 0.0;
  double transmission_uncertainty = 0.0;
  double buildup_factor = 0.0;
  double dose_reduction_factor = 0.0;
  std::shared_ptr<SimulationResults> detailed_results;
};

ShieldingAnalysisResult analyzeShielding(
  const Vector3D& source_position,
  double source_energy_MeV,
  const Material& shield_material, 
  double shield_thickness_cm,
  size_t num_histories = 100000);

/**
 * @brief Simplified dose rate calculation function
 */
struct DoseRateResult
{
  double dose_rate_Gy_per_s = 0.0;
  double dose_rate_uncertainty = 0.0;
  double total_energy_deposited_MeV = 0.0;
  std::map<std::string, double> dose_by_material;
  std::shared_ptr<SimulationResults> detailed_results;
};

DoseRateResult calculateDoseRate(
  const SourceConfiguration& source_config,
  std::shared_ptr<Geometry> phantom_geometry,
  size_t num_histories = 100000);

/**
 * @brief Helper functions for simulation analysis
 */
namespace SimulationUtils
{
  /**
   * @brief Compare simulation results with analytical solution
   */
  double compareWithBeerLambert(
    const SimulationResults& results,
    double initial_intensity,
    const Material& material,
    double thickness_cm,
    double energy_keV);
  
  /**
   * @brief Estimate required number of histories for target uncertainty
   */
  size_t estimateRequiredHistories(
    double target_relative_error,
    const SimulationResults& pilot_results);
  
  /**
   * @brief Generate comprehensive simulation report
   */
  std::string generateSimulationReport(
    const SimulationResults& results,
    const SimulationConfiguration& config);
  
  /**
   * @brief Validate simulation results for physical consistency
   */
  bool validateSimulationPhysics(
    const SimulationResults& results,
    double tolerance = 1e-6);
  
  /**
   * @brief Calculate figure of merit for simulation efficiency
   */
  double calculateFigureOfMerit(
    const SimulationResults& results,
    double simulation_time_seconds);
}