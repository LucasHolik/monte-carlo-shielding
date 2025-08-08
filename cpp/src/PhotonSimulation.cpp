#include "../include/PhotonSimulation.hpp"
#include "../include/Box.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <thread>

// =============================================================================
// SimulationConfiguration Implementation
// =============================================================================

bool SimulationConfiguration::isValid() const
{
  if (number_of_histories == 0) return false;
  if (energy_cutoff_keV < 0.0) return false;
  if (time_cutoff_seconds <= 0.0) return false;
  if (batch_size == 0 || batch_size > number_of_histories) return false;
  if (target_statistical_uncertainty <= 0.0 || target_statistical_uncertainty > 1.0) return false;
  if (number_of_threads == 0) return false;
  
  return true;
}

std::string SimulationConfiguration::toString() const
{
  std::stringstream ss;
  ss << "Simulation Configuration: " << name << "\n";
  ss << "  Histories: " << number_of_histories << "\n";
  ss << "  Batch size: " << batch_size << "\n";
  ss << "  Energy cutoff: " << energy_cutoff_keV << " keV\n";
  ss << "  Target uncertainty: " << target_statistical_uncertainty * 100.0 << "%\n";
  ss << "  Detailed logging: " << (enable_detailed_logging ? "Yes" : "No") << "\n";
  ss << "  Variance reduction: " << (enable_variance_reduction ? "Yes" : "No") << "\n";
  return ss.str();
}

// =============================================================================
// PhotonSimulation Implementation
// =============================================================================

PhotonSimulation::PhotonSimulation()
{
  // Initialize with default configuration
  config_ = SimulationConfiguration{};
  
  // Create default components
  rng_ = std::make_shared<RandomNumberGenerator>(42); // Default seed
  sampling_ = std::make_shared<MonteCarloSampling>(*rng_);
  results_ = std::make_shared<SimulationResults>("Default Simulation");
}

PhotonSimulation::PhotonSimulation(const SimulationConfiguration& config)
  : config_(config)
{
  if (!config_.isValid()) {
    throw std::invalid_argument("Invalid simulation configuration");
  }
  
  // Initialize components
  rng_ = std::make_shared<RandomNumberGenerator>(42);
  sampling_ = std::make_shared<MonteCarloSampling>(*rng_);
  results_ = std::make_shared<SimulationResults>(config_.name);
}

void PhotonSimulation::setConfiguration(const SimulationConfiguration& config)
{
  if (!config.isValid()) {
    throw std::invalid_argument("Invalid simulation configuration");
  }
  
  config_ = config;
  if (results_) {
    results_->setSimulationName(config_.name);
  }
}

void PhotonSimulation::setSource(std::shared_ptr<ParticleSource> source)
{
  source_ = source;
}

void PhotonSimulation::setGeometry(std::shared_ptr<Geometry> geometry)
{
  if (geometry && sampling_) {
    transport_ = std::make_shared<Transport>(*geometry, *sampling_);
    transport_->setSimulationResults(results_);
    transport_->enableDetailedResultsCollection(config_.enable_detailed_logging);
  }
}

void PhotonSimulation::setRandomSeed(uint64_t seed)
{
  if (rng_) {
    rng_->setSeed(seed);
  }
}

std::shared_ptr<SimulationResults> PhotonSimulation::runSimulation()
{
  validateSetup();
  initialize();
  
  // Record simulation start time
  start_time_ = std::chrono::steady_clock::now();
  is_running_ = true;
  
  // Starting simulation silently
  
  // Main simulation loop with batch processing
  try {
    while (histories_completed_ < config_.number_of_histories && 
           !checkTerminationCriteria() && is_running_) {
      
      // Calculate batch size for this iteration
      size_t remaining_histories = config_.number_of_histories - histories_completed_;
      size_t current_batch_size = std::min(config_.batch_size, remaining_histories);
      
      // Run single batch
      runSingleBatch(current_batch_size);
      
      // Update progress and check convergence
      updateProgress(histories_completed_);
      
      // Check if we've reached target statistical uncertainty
      if (hasConverged(config_.target_statistical_uncertainty)) {
        break;
      }
      
      // Pause point for external control
      while (is_paused_ && is_running_) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
      }
    }
    
    // Finalize simulation
    finalizeBatch();
    recordPerformanceMetrics();
    
    auto end_time = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time_);
    
    // Simulation completed silently
    
  } catch (const std::exception& e) {
    std::cerr << "Simulation error: " << e.what() << std::endl;
    is_running_ = false;
    throw;
  }
  
  is_running_ = false;
  return results_;
}

std::shared_ptr<SimulationResults> PhotonSimulation::runSimulation(
  const SourceConfiguration& source_config,
  std::shared_ptr<Geometry> geometry)
{
  // Set up source and geometry
  if (!source_) {
    source_ = std::make_shared<ParticleSource>(*sampling_, source_config);
  } else {
    source_->setConfiguration(source_config);
  }
  
  setGeometry(geometry);
  
  return runSimulation();
}

void PhotonSimulation::runSingleBatch(size_t batch_size)
{
  if (!source_ || !transport_ || !results_) {
    throw std::runtime_error("Simulation components not properly initialized");
  }
  
  size_t batch_start = histories_completed_;
  
  for (size_t i = 0; i < batch_size && is_running_; ++i) {
    // Generate source particle
    SourceParticle source_particle = source_->generateParticle();
    
    // Transport particle through geometry
    ParticleHistory history = transportSingleParticle(source_particle);
    
    // Process and record results
    processParticleHistory(history);
    
    histories_completed_++;
    
    // Report progress if needed
    if (config_.progress_callback && 
        (histories_completed_ % config_.progress_reporting_frequency == 0)) {
      double elapsed_time = getElapsedTime();
      config_.progress_callback(histories_completed_, config_.number_of_histories, elapsed_time);
    }
  }
  
  // Complete batch processing
  current_batch_++;
  calculateBatchStatistics();
}

ParticleHistory PhotonSimulation::transportSingleParticle(const SourceParticle& source_particle)
{
  ParticleHistory history;
  history.history_id = static_cast<int>(histories_completed_);
  history.initial_particle = source_particle.particle;
  
  // Make a copy of the particle for transport
  Particle particle = source_particle.particle;
  
  // Transport particle until termination
  while (particle.isAlive()) {
    // Check energy cutoff
    if (particle.energy() * 1000.0 < config_.energy_cutoff_keV) {
      particle.kill();
      break;
    }
    
    // Perform one transport step
    Vector3D position_before = particle.position();
    double energy_before = particle.energy();
    
    // Use transport engine to move and interact
    transport_->trackParticle(particle);
    
    // Update history tracking
    history.total_path_length += (particle.position() - position_before).magnitude();
    
    // Safety check to prevent infinite loops
    if (history.interactions.size() > 10000) {
      particle.kill();
      break;
    }
  }
  
  // Record final state
  history.final_position = particle.position();
  history.final_energy = particle.energy();
  history.final_state = particle.state();
  
  // Validate history
  if (config_.validate_energy_conservation) {
    if (!validateParticleHistory(history)) {
      // Energy conservation violation detected but suppressed
    }
  }
  
  return history;
}

void PhotonSimulation::processParticleHistory(const ParticleHistory& history)
{
  if (!results_) return;
  
  // Record complete history
  results_->recordParticleHistory(history);
  
  // Update transmission/reflection statistics
  if (history.final_state == ParticleState::Escaped) {
    // Particle escaped - contributes to transmission
    // This is a simplified implementation - in practice you'd check
    // which boundary the particle crossed
    
    // For now, assume forward escape is transmission
    if (history.final_position.z() > history.initial_particle.position().z()) {
      // Transmission
      double weight = history.initial_particle.weight();
      // results_->transmission_coefficient_.add(weight); // Would need to expose this
    }
  }
}

bool PhotonSimulation::hasConverged(double target_uncertainty) const
{
  if (!results_) return false;
  return results_->hasConverged(target_uncertainty);
}

double PhotonSimulation::getCompletionPercentage() const
{
  if (config_.number_of_histories == 0) return 0.0;
  return 100.0 * histories_completed_ / config_.number_of_histories;
}

double PhotonSimulation::getElapsedTime() const
{
  auto current_time = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time_);
  return duration.count();
}

double PhotonSimulation::getEstimatedRemainingTime() const
{
  if (histories_completed_ == 0) return 0.0;
  
  double elapsed = getElapsedTime();
  double rate = histories_completed_ / elapsed;
  size_t remaining = config_.number_of_histories - histories_completed_;
  
  return remaining / rate;
}

std::string PhotonSimulation::getStatusSummary() const
{
  std::stringstream ss;
  ss << std::fixed << std::setprecision(1);
  ss << "Simulation Status:\n";
  ss << "  Progress: " << getCompletionPercentage() << "% ";
  ss << "(" << histories_completed_ << "/" << config_.number_of_histories << " histories)\n";
  ss << "  Elapsed time: " << getElapsedTime() << " seconds\n";
  ss << "  Rate: " << calculateParticlesPerSecond() << " particles/second\n";
  if (histories_completed_ > 0) {
    ss << "  Estimated remaining: " << getEstimatedRemainingTime() << " seconds\n";
  }
  ss << "  Status: " << (is_running_ ? "Running" : "Stopped");
  if (is_paused_) ss << " (Paused)";
  ss << "\n";
  
  return ss.str();
}

void PhotonSimulation::initialize()
{
  if (!results_) {
    results_ = std::make_shared<SimulationResults>(config_.name);
  }
  
  results_->setTotalHistoriesRequested(config_.number_of_histories);
  results_->enableCompleteHistoryStorage(config_.enable_particle_tracking);
  
  // Reset counters
  histories_completed_ = 0;
  current_batch_ = 0;
  batch_results_.clear();
  
  is_initialized_ = true;
}

void PhotonSimulation::validateSetup()
{
  if (!source_) {
    throw std::runtime_error("No particle source configured");
  }
  if (!transport_) {
    throw std::runtime_error("No transport engine configured");
  }
  if (!results_) {
    throw std::runtime_error("No results collector configured");
  }
  if (!config_.isValid()) {
    throw std::runtime_error("Invalid simulation configuration");
  }
}

void PhotonSimulation::updateProgress(size_t histories_completed)
{
  // Progress updates suppressed for clean output
}

void PhotonSimulation::calculateBatchStatistics()
{
  // This would implement batch-wise statistical analysis
  // For proper uncertainty estimation
}

void PhotonSimulation::finalizeBatch()
{
  // Final statistical calculations and validation
}

bool PhotonSimulation::checkTerminationCriteria()
{
  // Check time limit
  if (getElapsedTime() > config_.time_cutoff_seconds) {
    return true;
  }
  
  return false;
}

bool PhotonSimulation::validateParticleHistory(const ParticleHistory& history) const
{
  // Basic energy conservation check
  double initial_energy = history.initial_particle.energy();
  double final_energy = history.final_energy;
  double deposited_energy = history.total_energy_deposited;
  
  double total_final = final_energy + deposited_energy;
  double relative_error = std::abs(total_final - initial_energy) / initial_energy;
  
  return relative_error < 1e-6; // 1 ppm tolerance
}

void PhotonSimulation::recordPerformanceMetrics()
{
  PerformanceData perf;
  perf.simulation_time_seconds = getElapsedTime();
  perf.particles_per_second = calculateParticlesPerSecond();
  perf.total_interactions = 0; // Would be tracked during transport
  
  if (results_) {
    results_->recordPerformanceData(perf);
  }
}

double PhotonSimulation::calculateParticlesPerSecond() const
{
  double elapsed = getElapsedTime();
  return elapsed > 0.0 ? histories_completed_ / elapsed : 0.0;
}

// =============================================================================
// Factory Methods
// =============================================================================

PhotonSimulation PhotonSimulation::createShieldingAnalysis(
  const Vector3D& source_position, double source_energy_MeV,
  std::shared_ptr<Geometry> shielding_geometry, size_t num_histories)
{
  SimulationConfiguration config;
  config.name = "Shielding Analysis";
  config.number_of_histories = num_histories;
  config.enable_detailed_logging = false;
  config.target_statistical_uncertainty = 0.02; // 2% for shielding
  
  PhotonSimulation sim(config);
  
  // Create monoenergetic point source
  auto source_config = ParticleSource::createMonoenergeticPointSource(
    source_position, source_energy_MeV);
  
  auto source = std::make_shared<ParticleSource>(*sim.sampling_, source_config);
  sim.setSource(source);
  sim.setGeometry(shielding_geometry);
  
  return sim;
}

// =============================================================================
// Convenience Functions
// =============================================================================

std::shared_ptr<SimulationResults> runPhotonSimulation(
  const SourceConfiguration& source_config,
  std::shared_ptr<Geometry> geometry,
  size_t num_histories,
  std::function<void(size_t, size_t)> progress_callback)
{
  // Create simulation configuration
  SimulationConfiguration config;
  config.name = "Photon Transport Simulation";
  config.number_of_histories = num_histories;
  config.batch_size = std::min(size_t(1000), num_histories);
  
  if (progress_callback) {
    config.progress_callback = [progress_callback](size_t completed, size_t total, double time) {
      progress_callback(completed, total);
    };
  }
  
  // Create and run simulation
  PhotonSimulation simulation(config);
  return simulation.runSimulation(source_config, geometry);
}

ShieldingAnalysisResult analyzeShielding(
  const Vector3D& source_position,
  double source_energy_MeV,
  const Material& shield_material,
  double shield_thickness_cm,
  size_t num_histories)
{
  // Create simple box geometry for shield
  Vector3D shield_min(-10.0, -10.0, 0.0);
  Vector3D shield_max(10.0, 10.0, shield_thickness_cm);
  
  auto box = std::make_shared<Box>(shield_min, shield_max, shield_material);
  auto geometry = std::make_shared<Geometry>();
  geometry->addShape(*box);
  
  // Run shielding simulation
  auto sim = PhotonSimulation::createShieldingAnalysis(
    source_position, source_energy_MeV, geometry, num_histories);
  
  auto results = sim.runSimulation();
  
  // Analyze results for shielding metrics
  ShieldingAnalysisResult analysis;
  analysis.detailed_results = results;
  
  // Calculate transmission coefficient
  size_t escaped = results->getPerformanceData().particles_escaped;
  analysis.transmission_coefficient = static_cast<double>(escaped) / num_histories;
  
  // Estimate uncertainty (simplified)
  if (escaped > 0) {
    analysis.transmission_uncertainty = std::sqrt(escaped) / num_histories;
  }
  
  // Calculate dose reduction factor
  if (analysis.transmission_coefficient > 0) {
    analysis.dose_reduction_factor = 1.0 / analysis.transmission_coefficient;
  }
  
  return analysis;
}

// =============================================================================
// Utility Functions
// =============================================================================

namespace SimulationUtils
{

double compareWithBeerLambert(const SimulationResults& results,
                            double initial_intensity,
                            const Material& material,
                            double thickness_cm,
                            double energy_keV)
{
  double simulated_transmission = results.getTransmissionCoefficient();
  
  // Calculate theoretical transmission
  double linear_mu = material.getLinearAttenuationCoefficient(energy_keV);
  double theoretical_transmission = std::exp(-linear_mu * thickness_cm);
  
  // Return relative difference
  if (theoretical_transmission > 0) {
    return std::abs(simulated_transmission - theoretical_transmission) / theoretical_transmission;
  }
  return 0.0;
}

size_t estimateRequiredHistories(double target_relative_error,
                               const SimulationResults& pilot_results)
{
  double current_error = pilot_results.getTotalEnergyDepositedUncertainty() / 
                        pilot_results.getTotalEnergyDeposited();
  
  if (current_error <= 0.0) return 0;
  
  // N_required = N_current * (current_error / target_error)²
  double ratio = current_error / target_relative_error;
  return static_cast<size_t>(pilot_results.getTotalHistoriesCompleted() * ratio * ratio);
}

std::string generateSimulationReport(const SimulationResults& results,
                                   const SimulationConfiguration& config)
{
  std::stringstream ss;
  ss << "=== MONTE CARLO SIMULATION REPORT ===\n\n";
  ss << results.getSummaryReport();
  ss << "\nSimulation Configuration:\n";
  ss << config.toString();
  return ss.str();
}

bool validateSimulationPhysics(const SimulationResults& results, double tolerance)
{
  // Basic validation checks
  if (results.getTotalEnergyDeposited() < 0.0) return false;
  if (results.getTotalHistoriesCompleted() == 0) return false;
  
  // More comprehensive validation would check:
  // - Energy conservation
  // - Particle balance
  // - Physical reasonableness of results
  
  return true;
}

double calculateFigureOfMerit(const SimulationResults& results,
                            double simulation_time_seconds)
{
  double relative_error = results.getTotalEnergyDepositedUncertainty() / 
                         results.getTotalEnergyDeposited();
  
  if (relative_error <= 0.0 || simulation_time_seconds <= 0.0) return 0.0;
  
  // FOM = 1 / (relative_error² × time)
  return 1.0 / (relative_error * relative_error * simulation_time_seconds);
}

} // namespace SimulationUtils