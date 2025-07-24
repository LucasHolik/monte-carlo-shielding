#pragma once

#include "Material.hpp"
#include "RandomNumberGenerator.hpp"
#include "Vector3D.hpp"

#include <functional>
#include <optional>
#include <string>
#include <variant>
#include <vector>

/**
 * @brief Energy spectrum types for particle source sampling
 */
enum class SpectrumType
{
  Monoenergetic, // Single energy
  Uniform,       // Uniform distribution
  Maxwell,       // Maxwell-Boltzmann
  Exponential,   // Exponential decay
  PowerLaw,      // Power law spectrum
  Discrete,      // Discrete energy lines
  Custom         // User-defined function
};

/**
 * @brief Energy spectrum definition for source sampling
 */
struct EnergySpectrum
{
  SpectrumType type;
  double min_energy = 0.0;            // Minimum energy (MeV)
  double max_energy = 10.0;           // Maximum energy (MeV)
  double characteristic_energy = 1.0; // Characteristic energy parameter
  double parameter = 1.0; // Additional parameter (alpha, temperature, etc.)
  std::vector<double> discrete_energies; // For discrete spectrum
  std::vector<double> discrete_weights;  // Weights for discrete energies
  std::function<double(double)>
      custom_pdf; // Custom probability density function

  // Factory methods
  static EnergySpectrum monoenergetic(double energy);
  static EnergySpectrum uniform(double min_energy, double max_energy);
  static EnergySpectrum maxwell(double temperature);
  static EnergySpectrum exponential(double characteristic_energy);
  static EnergySpectrum powerLaw(double alpha, double min_energy,
                                 double max_energy);
  static EnergySpectrum discrete(const std::vector<double> &energies,
                                 const std::vector<double> &weights);
};

/**
 * @brief Direction sampling types for particle sources
 */
enum class DirectionType
{
  Isotropic,      // Uniform over 4π steradians
  Beam,           // Pencil beam in specific direction
  Cone,           // Uniform within cone
  CosineWeighted, // Cosine-weighted (Lambert's law)
  Custom          // User-defined angular distribution
};

/**
 * @brief Angular distribution definition
 */
struct AngularDistribution
{
  DirectionType type;
  Vector3D reference_direction = Vector3D::UNIT_Z; // Reference direction
  double cone_angle = 0.0;                         // Cone half-angle (radians)
  std::function<Vector3D()> custom_sampler;        // Custom direction sampler

  // Factory methods
  static AngularDistribution isotropic();
  static AngularDistribution beam(const Vector3D &direction);
  static AngularDistribution cone(const Vector3D &axis, double half_angle);
  static AngularDistribution cosineWeighted(const Vector3D &normal);
};

/**
 * @brief Interaction sampling result
 */
struct InteractionSample
{
  double distance;            // Distance to interaction (cm)
  Vector3D interaction_point; // Position of interaction
  Material material;          // Material where interaction occurs
  bool escaped = false;       // True if particle escaped without interaction

  InteractionSample(double dist, const Vector3D &point, const Material &mat)
      : distance(dist), interaction_point(point), material(mat)
  {}

  InteractionSample() : distance(0.0), escaped(true) {} // Escape constructor
};

/**
 * @brief Monte Carlo sampling utilities for radiation transport
 *
 * Modern C++17 implementation providing comprehensive sampling techniques
 * specifically designed for Monte Carlo radiation transport simulations.
 */
class MonteCarloSampling
{
private:
  RandomNumberGenerator &rng_;

public:
  explicit MonteCarloSampling(RandomNumberGenerator &rng);

  // ============================================================================
  // EXPONENTIAL SAMPLING FOR INTERACTION DISTANCES
  // ============================================================================

  /**
   * @brief Sample distance to next interaction using exponential distribution
   * @param linear_attenuation_coefficient μ (cm⁻¹)
   * @return Distance to interaction (cm)
   */
  double sampleInteractionDistance(double linear_attenuation_coefficient) const;

  /**
   * @brief Safe exponential sampling with validation
   */
  std::optional<double> tryExponentialSampling(double lambda) const;

  /**
   * @brief Sample interaction distance in material with given cross-section
   * @param material Material properties
   * @param particle_energy Energy of particle (MeV)
   * @param cross_section_function Function to calculate cross-section
   * @return Distance to interaction
   */
  double
  sampleInteractionDistance(const Material &material, double particle_energy,
                            std::function<double(const Material &, double)>
                                cross_section_function) const;

  /**
   * @brief Sample interaction distance with multiple interaction types
   * @param total_cross_section Total macroscopic cross-section (cm⁻¹)
   * @return Distance to interaction
   */
  double sampleTotalInteractionDistance(double total_cross_section) const;

  // ============================================================================
  // UNIFORM SPHERICAL DIRECTION SAMPLING
  // ============================================================================

  /**
   * @brief Sample isotropic direction uniformly over 4π steradians
   * @return Unit vector with isotropic direction
   */
  Vector3D sampleIsotropicDirection() const;

  /**
   * @brief Sample direction from angular distribution
   * @param distribution Angular distribution specification
   * @return Unit vector with sampled direction
   */
  Vector3D sampleDirection(const AngularDistribution &distribution) const;

  /**
   * @brief Sample direction within cone
   * @param axis Cone axis direction
   * @param half_angle Cone half-angle (radians)
   * @return Unit vector within cone
   */
  Vector3D sampleConeDirection(const Vector3D &axis, double half_angle) const;

  /**
   * @brief Sample cosine-weighted direction (Lambert's law)
   * @param normal Surface normal
   * @return Unit vector with cosine-weighted distribution
   */
  Vector3D sampleCosineDirection(const Vector3D &normal) const;

  /**
   * @brief Sample direction in hemisphere uniformly
   * @param normal Hemisphere normal
   * @return Unit vector in hemisphere
   */
  Vector3D sampleHemisphereDirection(const Vector3D &normal) const;

  // ============================================================================
  // ENERGY DISTRIBUTION SAMPLING
  // ============================================================================

  /**
   * @brief Sample energy from spectrum
   * @param spectrum Energy spectrum definition
   * @return Sampled energy (MeV)
   */
  double sampleEnergy(const EnergySpectrum &spectrum) const;

  /**
   * @brief Sample monoenergetic particle
   * @param energy Fixed energy (MeV)
   * @return Energy value
   */
  double sampleMonoenergeticEnergy(double energy) const;

  /**
   * @brief Sample energy from uniform distribution
   * @param min_energy Minimum energy (MeV)
   * @param max_energy Maximum energy (MeV)
   * @return Sampled energy
   */
  double sampleUniformEnergy(double min_energy, double max_energy) const;

  /**
   * @brief Sample energy from Maxwell-Boltzmann distribution
   * @param temperature Temperature parameter (MeV)
   * @return Sampled energy
   */
  double sampleMaxwellEnergy(double temperature) const;

  /**
   * @brief Sample energy from exponential distribution
   * @param characteristic_energy Characteristic energy (MeV)
   * @return Sampled energy
   */
  double sampleExponentialEnergy(double characteristic_energy) const;

  /**
   * @brief Sample energy from power law distribution
   * @param alpha Power law index
   * @param min_energy Minimum energy (MeV)
   * @param max_energy Maximum energy (MeV)
   * @return Sampled energy
   */
  double samplePowerLawEnergy(double alpha, double min_energy,
                              double max_energy) const;

  /**
   * @brief Sample energy from discrete spectrum
   * @param energies Available energy values (MeV)
   * @param weights Relative weights for each energy
   * @return Sampled energy
   */
  double sampleDiscreteEnergy(const std::vector<double> &energies,
                              const std::vector<double> &weights) const;

  // ============================================================================
  // ADVANCED MONTE CARLO TECHNIQUES
  // ============================================================================

  /**
   * @brief Sample scattering angle using rejection method
   * @param pdf Probability density function for angle
   * @param max_pdf Maximum value of PDF (for rejection sampling)
   * @return Scattering angle (radians)
   */
  double sampleScatteringAngle(std::function<double(double)> pdf,
                               double max_pdf) const;

  /**
   * @brief Sample azimuthal angle uniformly
   * @return Azimuthal angle [0, 2π) radians
   */
  double sampleAzimuthalAngle() const;

  /**
   * @brief Transform direction by scattering angles
   * @param incident_direction Original direction
   * @param polar_angle Polar scattering angle (radians)
   * @param azimuthal_angle Azimuthal scattering angle (radians)
   * @return New direction after scattering
   */
  Vector3D applyScattering(const Vector3D &incident_direction,
                           double polar_angle, double azimuthal_angle) const;

  /**
   * @brief Sample interaction type from cross-sections
   * @param cross_sections Vector of cross-section values
   * @return Index of selected interaction type
   */
  std::size_t
  sampleInteractionType(const std::vector<double> &cross_sections) const;

  /**
   * @brief Russian roulette for variance reduction
   * @param survival_probability Probability of survival [0,1]
   * @return True if particle survives
   */
  bool russianRoulette(double survival_probability) const;

  /**
   * @brief Particle splitting for variance reduction
   * @param splitting_factor Average number of split particles
   * @return Number of particles to create
   */
  std::size_t particleSplitting(double splitting_factor) const;

  // ============================================================================
  // PHYSICS-SPECIFIC SAMPLING
  // ============================================================================

  /**
   * @brief Sample Compton scattering angle using Klein-Nishina formula
   * @param incident_energy Incident photon energy (MeV)
   * @return Scattering angle (radians)
   */
  double sampleComptonAngle(double incident_energy) const;

  /**
   * @brief Sample scattered photon energy from Compton scattering
   * @param incident_energy Incident photon energy (MeV)
   * @param scattering_angle Scattering angle (radians)
   * @return Scattered photon energy (MeV)
   */
  double sampleComptonEnergy(double incident_energy,
                             double scattering_angle) const;

  /**
   * @brief Sample photoelectron direction (isotropic approximation)
   * @return Photoelectron direction
   */
  Vector3D samplePhotoelectronDirection() const;

  /**
   * @brief Sample pair production angles (simplified)
   * @param photon_energy Incident photon energy (MeV)
   * @return Pair of directions for electron and positron
   */
  std::pair<Vector3D, Vector3D>
  samplePairProductionDirections(double photon_energy) const;

  // ============================================================================
  // UTILITY METHODS
  // ============================================================================

  /**
   * @brief Validate sampling parameters
   */
  bool validateEnergySpectrum(const EnergySpectrum &spectrum) const;
  bool
  validateAngularDistribution(const AngularDistribution &distribution) const;

  /**
   * @brief Get random number generator reference
   */
  RandomNumberGenerator &getRNG() { return rng_; }
  const RandomNumberGenerator &getRNG() const { return rng_; }

  /**
   * @brief Sample multiple values efficiently
   */
  std::vector<double>
  sampleInteractionDistances(double linear_attenuation_coefficient,
                             std::size_t count) const;
  std::vector<Vector3D> sampleIsotropicDirections(std::size_t count) const;
  std::vector<double> sampleEnergies(const EnergySpectrum &spectrum,
                                     std::size_t count) const;
};

// ============================================================================
// CONVENIENCE FUNCTIONS
// ============================================================================

namespace MonteCarloUtils
{
/**
 * @brief Calculate Klein-Nishina differential cross-section
 * @param incident_energy Incident photon energy (MeV)
 * @param scattering_angle Scattering angle (radians)
 * @return Differential cross-section (relative units)
 */
double kleinNishinaCrossSection(double incident_energy,
                                double scattering_angle);

/**
 * @brief Convert between polar and Cartesian coordinates
 */
Vector3D polarToCartesian(double theta, double phi);
std::tuple<double, double> cartesianToPolar(const Vector3D &direction);

/**
 * @brief Transform vector from local to global coordinates
 * @param local_vector Vector in local coordinate system
 * @param local_z Local Z-axis in global coordinates
 * @return Vector in global coordinates
 */
Vector3D localToGlobal(const Vector3D &local_vector, const Vector3D &local_z);

/**
 * @brief Generate orthonormal basis from single vector
 * @param z Z-axis of the basis
 * @return Pair of orthonormal X and Y axes
 */
std::pair<Vector3D, Vector3D> generateOrthonormalBasis(const Vector3D &z);
} // namespace MonteCarloUtils