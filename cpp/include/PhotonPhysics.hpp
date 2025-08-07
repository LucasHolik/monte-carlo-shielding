#pragma once

#include "Material.hpp"
#include "MonteCarloSampling.hpp"
#include "Particle.hpp"
#include "Vector3D.hpp"
#include "PhotonCrossSectionDatabase.hpp"

#include <optional>
#include <vector>

/**
 * @brief Photon interaction types based on NIST XCOM database
 */
enum class PhotonInteractionType
{
  PHOTOELECTRIC_ABSORPTION,
  COMPTON_SCATTERING,
  COHERENT_SCATTERING, // Rayleigh scattering
  PAIR_PRODUCTION_NUCLEAR,
  PAIR_PRODUCTION_ELECTRON
};

/**
 * @brief Result of a photon interaction
 */
struct PhotonInteractionResult
{
  PhotonInteractionType interaction_type;
  bool particle_absorbed = false;
  bool particle_scattered = false;
  Vector3D new_direction;
  double new_energy = 0.0;
  std::vector<Particle> secondary_particles;

  PhotonInteractionResult() = default;
  PhotonInteractionResult(PhotonInteractionType type) : interaction_type(type) {}
};

/**
 * @brief Comprehensive photon physics implementation based on NIST XCOM data
 *
 * Implements the major photon interaction processes:
 * - Photoelectric absorption (Scofield calculations)
 * - Compton scattering (Klein-Nishina + corrections)
 * - Coherent (Rayleigh) scattering
 * - Pair production (nuclear and electronic fields)
 *
 * Cross-section data and formulas based on NIST XCOM database
 * Reference: https://physics.nist.gov/PhysRefData/Xcom/Text/intro.html
 */
class PhotonPhysics
{
private:
  MonteCarloSampling &sampling_;

  // Physical constants
  static constexpr double ELECTRON_REST_MASS_MEV = 0.51099895000; // MeV
  static constexpr double CLASSICAL_ELECTRON_RADIUS = 2.8179403262e-13; // cm
  static constexpr double FINE_STRUCTURE_CONSTANT = 7.2973525693e-3;
  static constexpr double PAIR_PRODUCTION_THRESHOLD = 2.0 * ELECTRON_REST_MASS_MEV; // MeV

public:
  // Constructor
  explicit PhotonPhysics(MonteCarloSampling &sampling);

  // Main interaction methods
  PhotonInteractionResult performInteraction(const Particle &photon,
                                             const Material &material);
  PhotonInteractionType sampleInteractionType(double energy,
                                              const Material &material);

  // Individual interaction implementations
  PhotonInteractionResult performPhotoelectricAbsorption(const Particle &photon,
                                                         const Material &material);
  PhotonInteractionResult performComptonScattering(const Particle &photon,
                                                   const Material &material);
  PhotonInteractionResult performCoherentScattering(const Particle &photon,
                                                    const Material &material);
  PhotonInteractionResult performPairProduction(const Particle &photon,
                                                const Material &material);

  // Cross-section calculations (based on NIST XCOM)
  double getPhotoelectricCrossSection(double energy, int atomic_number) const;
  double getComptonCrossSection(double energy, int atomic_number) const;
  double getCoherentCrossSection(double energy, int atomic_number) const;
  double getPairProductionCrossSection(double energy, int atomic_number) const;
  double getTotalCrossSection(double energy, const Material &material) const;

  // Klein-Nishina formula implementation
  double kleinNishinaCrossSection(double energy) const;
  double kleinNishinaDifferential(double energy, double cos_theta) const;

  // Energy and angle sampling for interactions
  double sampleComptonScatteredEnergy(double incident_energy);
  double sampleComptonScatteringAngle(double incident_energy);
  Vector3D sampleScatteredDirection(const Vector3D &incident_direction,
                                    double cos_theta, double phi);

  // Utility methods
  bool isPhotoelectricPossible(double energy, const Material &material) const;
  bool isPairProductionPossible(double energy) const;
  double calculateMassAttenuationCoefficient(double energy,
                                             const Material &material) const;

  // Validation methods
  bool validateInteractionResult(const PhotonInteractionResult &result) const;
  std::string getInteractionTypeName(PhotonInteractionType type) const;

private:
  // Internal calculation helpers
  double calculatePhotoelectricCoefficient(double energy, int Z) const;
  double calculateComptonCoefficient(double energy, int Z) const;
  double calculateCoherentCoefficient(double energy, int Z) const;
  double calculatePairProductionCoefficient(double energy, int Z) const;

  // Simplified approximations for cross-sections (to be replaced with tabulated data)
  double approximatePhotoelectricCrossSection(double energy, int Z) const;
  double approximateCoherentCrossSection(double energy, int Z) const;
  double approximatePairProductionCrossSection(double energy, int Z) const;

  // Energy and angle sampling helpers
  double samplePhiAngle();
  Vector3D rotateDirection(const Vector3D &original, double cos_theta, double phi);
};

/**
 * @brief Constants and utility functions for photon physics
 */
namespace PhotonConstants
{
// Energy thresholds (MeV)
constexpr double PHOTOELECTRIC_THRESHOLD = 0.001; // 1 keV
constexpr double COMPTON_THRESHOLD = 0.001;       // 1 keV  
constexpr double PAIR_PRODUCTION_THRESHOLD = 1.022; // 2 * m_e c^2

// Cross-section scaling factors
constexpr double BARN_TO_CM2 = 1e-24; // barn to cm²
constexpr double AVOGADRO = 6.02214076e23; // mol⁻¹

/**
 * @brief Convert mass attenuation coefficient to linear attenuation coefficient
 * @param mass_attenuation Mass attenuation coefficient (cm²/g)
 * @param density Material density (g/cm³)
 * @return Linear attenuation coefficient (cm⁻¹)
 */
inline double massToLinearAttenuation(double mass_attenuation, double density)
{
  return mass_attenuation * density;
}

/**
 * @brief Convert photon energy between units
 */
inline double kevToMeV(double energy_kev) { return energy_kev / 1000.0; }
inline double mevToKeV(double energy_mev) { return energy_mev * 1000.0; }

} // namespace PhotonConstants