#pragma once

#include "MonteCarloSampling.hpp"
#include "Particle.hpp"
#include "Vector3D.hpp"

#include <functional>
#include <optional>
#include <string>
#include <variant>
#include <vector>

/**
 * @brief Source geometry types for particle emission
 */
enum class SourceGeometryType
{
  Point,           // Point source at specific location
  Plane,           // Planar source (infinite or finite)
  Sphere,          // Spherical source (surface or volume)
  Cylinder,        // Cylindrical source
  Cone,            // Conical emission
  Rectangle,       // Rectangular source
  Custom           // User-defined geometry function
};

/**
 * @brief Source intensity variation types
 */
enum class IntensityDistributionType
{
  Uniform,         // Uniform intensity
  Gaussian,        // Gaussian (normal) distribution
  Exponential,     // Exponential falloff
  Custom           // User-defined intensity function
};

/**
 * @brief Geometric source definition
 */
struct SourceGeometry
{
  SourceGeometryType type = SourceGeometryType::Point;
  Vector3D center = Vector3D(0, 0, 0);          // Source center position
  Vector3D direction = Vector3D(0, 0, 1);       // Reference direction
  double parameter1 = 0.0;                       // Radius, width, etc.
  double parameter2 = 0.0;                       // Height, length, etc.
  double parameter3 = 0.0;                       // Additional parameter
  std::function<Vector3D()> custom_position_sampler; // For custom geometries
  
  // Factory methods
  static SourceGeometry createPoint(const Vector3D& position);
  static SourceGeometry createPlane(const Vector3D& center, const Vector3D& normal, 
                                   double width = 1e30, double height = 1e30);
  static SourceGeometry createSphere(const Vector3D& center, double radius, bool surface_only = true);
  static SourceGeometry createCylinder(const Vector3D& center, const Vector3D& axis, 
                                      double radius, double height, bool surface_only = true);
  static SourceGeometry createRectangle(const Vector3D& corner, const Vector3D& side1, 
                                       const Vector3D& side2);
  static SourceGeometry createCone(const Vector3D& apex, const Vector3D& axis, 
                                  double base_radius, double height);
};

/**
 * @brief Intensity distribution definition
 */
struct IntensityDistribution
{
  IntensityDistributionType type = IntensityDistributionType::Uniform;
  Vector3D center = Vector3D(0, 0, 0);          // Center of distribution
  double parameter1 = 1.0;                       // Standard deviation, decay constant, etc.
  double parameter2 = 1.0;                       // Additional parameter
  std::function<double(const Vector3D&)> custom_intensity_function; // For custom distributions
  
  // Factory methods
  static IntensityDistribution createUniform();
  static IntensityDistribution createGaussian(const Vector3D& center, double sigma);
  static IntensityDistribution createExponential(const Vector3D& center, double lambda);
};

/**
 * @brief Complete source configuration combining all aspects
 */
struct SourceConfiguration
{
  std::string name = "Default Source";
  ParticleType particle_type = ParticleType::Photon;
  
  // Energy distribution (using existing EnergySpectrum from MonteCarloSampling.hpp)
  EnergySpectrum energy_spectrum;
  
  // Angular distribution (using existing AngularDistribution from MonteCarloSampling.hpp)  
  AngularDistribution angular_distribution;
  
  // Spatial distribution
  SourceGeometry geometry;
  IntensityDistribution intensity;
  
  // Source strength and timing
  double activity = 1.0;                         // Particles per second
  double total_particles = 1000;                 // Total particles to emit
  double simulation_time = 1.0;                  // Simulation time (seconds)
  
  // Statistical weights
  double base_statistical_weight = 1.0;          // Base weight for variance reduction
  
  // Validation
  bool isValid() const;
  std::string toString() const;
};

/**
 * @brief Sampled particle from source with complete initial state
 */
struct SourceParticle
{
  Particle particle;                             // Complete particle object
  double emission_time = 0.0;                   // Time of emission (seconds)
  double statistical_weight = 1.0;              // Statistical weight
  Vector3D emission_position;                   // Position where emitted
  Vector3D emission_direction;                  // Initial direction
  double emission_energy = 0.0;                 // Initial energy (MeV)
  
  SourceParticle() = default;
  SourceParticle(const Particle& p, double time = 0.0) 
    : particle(p), emission_time(time) {}
};

/**
 * @brief Professional particle source for Monte Carlo simulations
 * 
 * Modern C++17 implementation supporting various source geometries,
 * energy spectra, and angular distributions with proper statistical
 * sampling and validation.
 */
class ParticleSource
{
private:
  MonteCarloSampling& sampling_;
  SourceConfiguration config_;
  size_t particles_generated_ = 0;
  double current_simulation_time_ = 0.0;
  
  // Internal sampling helpers
  Vector3D sampleSourcePosition() const;
  Vector3D sampleSourceDirection(const Vector3D& position) const;
  double sampleSourceEnergy() const;
  double sampleIntensityWeight(const Vector3D& position) const;
  double sampleEmissionTime() const;
  
public:
  // Constructors
  explicit ParticleSource(MonteCarloSampling& sampling);
  ParticleSource(MonteCarloSampling& sampling, const SourceConfiguration& config);
  
  // Configuration
  void setConfiguration(const SourceConfiguration& config);
  const SourceConfiguration& getConfiguration() const { return config_; }
  
  // Individual source component setters
  void setParticleType(ParticleType type) { config_.particle_type = type; }
  void setEnergySpectrum(const EnergySpectrum& spectrum) { config_.energy_spectrum = spectrum; }
  void setAngularDistribution(const AngularDistribution& distribution) { config_.angular_distribution = distribution; }
  void setSourceGeometry(const SourceGeometry& geometry) { config_.geometry = geometry; }
  void setIntensityDistribution(const IntensityDistribution& intensity) { config_.intensity = intensity; }
  void setActivity(double activity) { config_.activity = activity; }
  void setTotalParticles(size_t total) { config_.total_particles = total; }
  
  // Particle generation
  SourceParticle generateParticle();
  std::vector<SourceParticle> generateParticles(size_t count);
  std::vector<SourceParticle> generateParticlesForTime(double time_duration);
  
  // Batch generation with progress callback
  std::vector<SourceParticle> generateParticlesBatch(size_t count, 
    std::function<void(size_t, size_t)> progress_callback = nullptr);
  
  // Source statistics and information
  size_t getParticlesGenerated() const { return particles_generated_; }
  double getCurrentSimulationTime() const { return current_simulation_time_; }
  double getExpectedEmissionRate() const { return config_.activity; }
  
  // Validation and diagnostics
  bool isConfigurationValid() const { return config_.isValid(); }
  std::vector<std::string> validateConfiguration() const;
  std::string getSourceSummary() const;
  
  // Reset and cleanup
  void reset();
  void setSimulationTime(double time) { current_simulation_time_ = time; }
  
  // Common source factory methods
  static SourceConfiguration createMonoenergeticPointSource(
    const Vector3D& position, double energy_MeV, ParticleType type = ParticleType::Photon);
  
  static SourceConfiguration createIsotropicPointSource(
    const Vector3D& position, const EnergySpectrum& spectrum, ParticleType type = ParticleType::Photon);
  
  static SourceConfiguration createBeamSource(
    const Vector3D& position, const Vector3D& direction, double energy_MeV,
    double beam_divergence = 0.0, ParticleType type = ParticleType::Photon);
  
  static SourceConfiguration createPlanarSource(
    const Vector3D& center, const Vector3D& normal, double width, double height,
    const EnergySpectrum& spectrum, ParticleType type = ParticleType::Photon);
  
  static SourceConfiguration createSphericalSource(
    const Vector3D& center, double radius, const EnergySpectrum& spectrum,
    bool surface_only = true, ParticleType type = ParticleType::Photon);
  
  // Advanced source configurations
  static SourceConfiguration createGammaCalibratedSource(
    const Vector3D& position, const std::vector<double>& energies_keV, 
    const std::vector<double>& intensities);
  
  static SourceConfiguration createMedicalLinacSource(
    const Vector3D& position, const Vector3D& beam_direction, double nominal_energy_MeV);
  
  static SourceConfiguration createCosmicRaySource(
    const Vector3D& detector_center, double detector_radius);
  
private:
  // Internal helper methods
  void validateAndThrow() const;
  Vector3D samplePointSourcePosition() const;
  Vector3D samplePlaneSourcePosition() const;
  Vector3D sampleSphereSourcePosition() const;
  Vector3D sampleCylinderSourcePosition() const;
  Vector3D sampleRectangleSourcePosition() const;
  Vector3D sampleConeSourcePosition() const;
  
  bool isPositionWithinGeometry(const Vector3D& position) const;
  double calculateGeometryArea() const;
  double calculateGeometryVolume() const;
};

/**
 * @brief Helper functions for source analysis and validation
 */
namespace SourceAnalysis
{
  /**
   * @brief Calculate total source strength for given configuration
   * @param config Source configuration
   * @return Total particles per second
   */
  double calculateTotalSourceStrength(const SourceConfiguration& config);
  
  /**
   * @brief Estimate memory requirements for storing source particles
   * @param config Source configuration
   * @return Estimated memory usage in MB
   */
  double estimateMemoryRequirement(const SourceConfiguration& config);
  
  /**
   * @brief Validate source geometry for physical reasonableness
   * @param geometry Source geometry to validate
   * @return True if geometry is physically reasonable
   */
  bool validateSourceGeometry(const SourceGeometry& geometry);
  
  /**
   * @brief Calculate effective source area for flux calculations
   * @param geometry Source geometry
   * @return Effective area in cm²
   */
  double calculateEffectiveSourceArea(const SourceGeometry& geometry);
  
  /**
   * @brief Generate quality assurance report for source configuration
   * @param config Source configuration to analyze
   * @return QA report string
   */
  std::string generateQualityAssuranceReport(const SourceConfiguration& config);
}