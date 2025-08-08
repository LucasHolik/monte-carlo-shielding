#include "PhotonPhysics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

PhotonPhysics::PhotonPhysics(MonteCarloSampling &sampling) : sampling_(sampling) {}

PhotonInteractionResult PhotonPhysics::performInteraction(const Particle &photon,
                                                          const Material &material)
{
  if(photon.energy() <= 0.0)
  {
    throw std::invalid_argument("Photon energy must be positive");
  }

  // Sample interaction type based on relative cross-sections
  auto interaction_type = sampleInteractionType(photon.energy(), material);

  // Perform the specific interaction
  switch(interaction_type)
  {
  case PhotonInteractionType::PHOTOELECTRIC_ABSORPTION:
    return performPhotoelectricAbsorption(photon, material);

  case PhotonInteractionType::COMPTON_SCATTERING:
    return performComptonScattering(photon, material);

  case PhotonInteractionType::COHERENT_SCATTERING:
    return performCoherentScattering(photon, material);

  case PhotonInteractionType::PAIR_PRODUCTION_NUCLEAR:
  case PhotonInteractionType::PAIR_PRODUCTION_ELECTRON:
    return performPairProduction(photon, material);

  default:
    throw std::runtime_error("Unknown photon interaction type");
  }
}

PhotonInteractionType PhotonPhysics::sampleInteractionType(double energy,
                                                           const Material &material)
{
  // Get cross-sections for all interaction types
  double photoelectric_xs = 0.0;
  double compton_xs = 0.0;
  double coherent_xs = 0.0;
  double pair_xs = 0.0;

  // Calculate cross-sections for material composition
  for(const auto &element : material.composition())
  {
    int Z = element.atomic_number;
    double fraction = element.weight_fraction;

    photoelectric_xs += fraction * getPhotoelectricCrossSection(energy, Z);
    compton_xs += fraction * getComptonCrossSection(energy, Z);
    coherent_xs += fraction * getCoherentCrossSection(energy, Z);
    
    if(isPairProductionPossible(energy))
    {
      pair_xs += fraction * getPairProductionCrossSection(energy, Z);
    }
  }

  double total_xs = photoelectric_xs + compton_xs + coherent_xs + pair_xs;

  if(total_xs <= 0.0)
  {
    // Fallback to Compton scattering
    return PhotonInteractionType::COMPTON_SCATTERING;
  }

  // Sample interaction type based on relative probabilities
  double random = sampling_.getRNG().uniform() * total_xs;

  if(random < photoelectric_xs)
  {
    return PhotonInteractionType::PHOTOELECTRIC_ABSORPTION;
  }
  else if(random < photoelectric_xs + compton_xs)
  {
    return PhotonInteractionType::COMPTON_SCATTERING;
  }
  else if(random < photoelectric_xs + compton_xs + coherent_xs)
  {
    return PhotonInteractionType::COHERENT_SCATTERING;
  }
  else
  {
    return PhotonInteractionType::PAIR_PRODUCTION_NUCLEAR;
  }
}

PhotonInteractionResult
PhotonPhysics::performPhotoelectricAbsorption(const Particle &photon,
                                              const Material &material)
{
  PhotonInteractionResult result(PhotonInteractionType::PHOTOELECTRIC_ABSORPTION);

  // Photoelectric absorption completely absorbs the photon
  result.particle_absorbed = true;
  result.new_energy = 0.0;

  // In a full implementation, would produce:
  // - Photoelectron with energy = E_photon - E_binding
  // - Characteristic X-rays and Auger electrons
  // For now, we assume all energy is deposited locally

  return result;
}

PhotonInteractionResult PhotonPhysics::performComptonScattering(const Particle &photon,
                                                                const Material &material)
{
  PhotonInteractionResult result(PhotonInteractionType::COMPTON_SCATTERING);

  double incident_energy = photon.energy();

  // Sample scattering angle using Klein-Nishina distribution
  double cos_theta = sampleComptonScatteringAngle(incident_energy);
  double phi = samplePhiAngle();

  // Calculate scattered photon energy using Compton formula
  double scattered_energy = sampleComptonScatteredEnergy(incident_energy);

  // Calculate new direction
  Vector3D new_direction = sampleScatteredDirection(photon.direction(), cos_theta, phi);

  // Set result
  result.particle_scattered = true;
  result.new_direction = new_direction;
  result.new_energy = scattered_energy;

  // In a full implementation, would also create recoil electron
  // with energy = incident_energy - scattered_energy

  return result;
}

PhotonInteractionResult PhotonPhysics::performCoherentScattering(const Particle &photon,
                                                                 const Material &material)
{
  PhotonInteractionResult result(PhotonInteractionType::COHERENT_SCATTERING);

  // Coherent (Rayleigh) scattering: photon changes direction but not energy
  result.particle_scattered = true;
  result.new_energy = photon.energy(); // No energy change

  // Sample scattering angle (simpler distribution for coherent scattering)
  // For coherent scattering, use simple isotropic approximation
  double cos_theta = 2.0 * sampling_.getRNG().uniform() - 1.0;
  double phi = samplePhiAngle();

  result.new_direction = sampleScatteredDirection(photon.direction(), cos_theta, phi);

  return result;
}

PhotonInteractionResult PhotonPhysics::performPairProduction(const Particle &photon,
                                                             const Material &material)
{
  PhotonInteractionResult result(PhotonInteractionType::PAIR_PRODUCTION_NUCLEAR);

  if(!isPairProductionPossible(photon.energy()))
  {
    throw std::runtime_error("Pair production requires energy > 1.022 MeV");
  }

  // Pair production completely absorbs the photon
  result.particle_absorbed = true;
  result.new_energy = 0.0;

  // Available kinetic energy for the pair
  double kinetic_energy = photon.energy() - PAIR_PRODUCTION_THRESHOLD;

  // In a full implementation, would create electron-positron pair
  // For now, assume energy is deposited locally

  return result;
}

double PhotonPhysics::getPhotoelectricCrossSection(double energy, int atomic_number) const
{
  // Convert energy from MeV to keV for the database
  double energy_keV = energy * 1000.0;
  
  // Use NIST-accurate Scofield photoelectric data
  double cross_section_barns = PhotonCrossSections::ScofieldPhotoelectricData::getPhotoelectricCrossSection(
      atomic_number, energy_keV);
  
  // Convert from barns to cm²
  return PhotonCrossSections::PhotonInteractionDatabase::barnsToCm2(cross_section_barns);
}

double PhotonPhysics::getComptonCrossSection(double energy, int atomic_number) const
{
  // Convert energy from MeV to keV for the database
  double energy_keV = energy * 1000.0;
  
  // Use accurate incoherent scattering with binding corrections
  double cross_section_barns = PhotonCrossSections::IncoherentScatteringData::getTotalIncoherentCrossSection(
      atomic_number, energy_keV);
  
  // Convert from barns to cm²
  return PhotonCrossSections::PhotonInteractionDatabase::barnsToCm2(cross_section_barns);
}

double PhotonPhysics::getCoherentCrossSection(double energy, int atomic_number) const
{
  // Convert energy from MeV to keV for the database
  double energy_keV = energy * 1000.0;
  
  // Use accurate coherent scattering with atomic form factors
  double cross_section_barns = PhotonCrossSections::CoherentScatteringData::getCoherentCrossSection(
      atomic_number, energy_keV);
  
  // Convert from barns to cm²
  return PhotonCrossSections::PhotonInteractionDatabase::barnsToCm2(cross_section_barns);
}

double PhotonPhysics::getPairProductionCrossSection(double energy, int atomic_number) const
{
  if(!isPairProductionPossible(energy))
  {
    return 0.0;
  }
  
  // Use accurate Bethe-Heitler pair production with screening corrections
  double cross_section_barns = PhotonCrossSections::PairProductionData::getTotalPairProductionCrossSection(
      atomic_number, energy);
  
  // Convert from barns to cm²
  return PhotonCrossSections::PhotonInteractionDatabase::barnsToCm2(cross_section_barns);
}

double PhotonPhysics::getTotalCrossSection(double energy, const Material &material) const
{
  double total = 0.0;

  for(const auto &element : material.composition())
  {
    int Z = element.atomic_number;
    double fraction = element.weight_fraction;

    total += fraction * (getPhotoelectricCrossSection(energy, Z) +
                        getComptonCrossSection(energy, Z) +
                        getCoherentCrossSection(energy, Z) +
                        getPairProductionCrossSection(energy, Z));
  }

  return total;
}

double PhotonPhysics::kleinNishinaCrossSection(double energy) const
{
  double alpha = energy / ELECTRON_REST_MASS_MEV;
  
  if(alpha < 1e-6)
  {
    // Thomson scattering limit for low energies
    return (8.0 / 3.0) * M_PI * std::pow(CLASSICAL_ELECTRON_RADIUS, 2);
  }

  // Klein-Nishina formula (correct total cross-section)
  // σ = 2πr_e² [(1+α)/α³] [(2α(1+α))/(1+2α) - ln(1+2α)] + ln(1+2α)/(2α) - (1+3α)/(1+2α)²
  
  double r_e_squared = std::pow(CLASSICAL_ELECTRON_RADIUS, 2);
  
  double bracket1 = (1.0 + alpha) / (alpha * alpha * alpha);
  double bracket2_term1 = (2.0 * alpha * (1.0 + alpha)) / (1.0 + 2.0 * alpha);
  double bracket2_term2 = std::log(1.0 + 2.0 * alpha);
  double bracket2 = bracket2_term1 - bracket2_term2;
  
  double term3 = std::log(1.0 + 2.0 * alpha) / (2.0 * alpha);
  double term4 = (1.0 + 3.0 * alpha) / std::pow(1.0 + 2.0 * alpha, 2);
  
  return 2.0 * M_PI * r_e_squared * (bracket1 * bracket2 + term3 - term4);
}

double PhotonPhysics::sampleComptonScatteredEnergy(double incident_energy)
{
  double alpha = incident_energy / ELECTRON_REST_MASS_MEV;
  
  // Sample cosine of scattering angle using Klein-Nishina distribution
  double cos_theta = sampleComptonScatteringAngle(incident_energy);
  
  // Calculate scattered energy using Compton formula
  double scattered_energy = incident_energy / (1.0 + alpha * (1.0 - cos_theta));
  
  return scattered_energy;
}

double PhotonPhysics::sampleComptonScatteringAngle(double incident_energy)
{
  double alpha = incident_energy / ELECTRON_REST_MASS_MEV;
  
  // Simplified sampling - use rejection method for Klein-Nishina distribution
  // For demonstration purposes, using approximate sampling
  
  double cos_theta;
  double rejection_value;
  double klein_nishina_value;
  
  do
  {
    cos_theta = 2.0 * sampling_.getRNG().uniform() - 1.0; // Uniform in [-1, 1]
    rejection_value = sampling_.getRNG().uniform();
    klein_nishina_value = kleinNishinaDifferential(incident_energy, cos_theta);
  } while(rejection_value > klein_nishina_value);
  
  return cos_theta;
}

double PhotonPhysics::kleinNishinaDifferential(double energy, double cos_theta) const
{
  double alpha = energy / ELECTRON_REST_MASS_MEV;
  double ratio = 1.0 / (1.0 + alpha * (1.0 - cos_theta));
  
  return ratio * ratio * (ratio + 1.0 / ratio - 1.0 + cos_theta * cos_theta);
}

Vector3D PhotonPhysics::sampleScatteredDirection(const Vector3D &incident_direction,
                                                 double cos_theta, double phi)
{
  return rotateDirection(incident_direction, cos_theta, phi);
}

bool PhotonPhysics::isPhotoelectricPossible(double energy, const Material &material) const
{
  return energy >= PhotonConstants::PHOTOELECTRIC_THRESHOLD;
}

bool PhotonPhysics::isPairProductionPossible(double energy) const
{
  return energy >= PAIR_PRODUCTION_THRESHOLD;
}

double PhotonPhysics::calculateMassAttenuationCoefficient(double energy,
                                                         const Material &material) const
{
  double total_cross_section = getTotalCrossSection(energy, material);
  
  // Convert to mass attenuation coefficient
  // Note: This is a simplified conversion - full implementation would account for
  // proper atomic mass and Avogadro's number scaling
  return total_cross_section / material.density();
}

bool PhotonPhysics::validateInteractionResult(const PhotonInteractionResult &result) const
{
  // Check energy conservation (simplified)
  if(result.new_energy < 0.0)
  {
    return false;
  }
  
  // Check direction normalization if scattered
  if(result.particle_scattered)
  {
    double magnitude = result.new_direction.magnitude();
    return std::abs(magnitude - 1.0) < 1e-6;
  }
  
  return true;
}

std::string PhotonPhysics::getInteractionTypeName(PhotonInteractionType type) const
{
  switch(type)
  {
  case PhotonInteractionType::PHOTOELECTRIC_ABSORPTION:
    return "Photoelectric Absorption";
  case PhotonInteractionType::COMPTON_SCATTERING:
    return "Compton Scattering";
  case PhotonInteractionType::COHERENT_SCATTERING:
    return "Coherent (Rayleigh) Scattering";
  case PhotonInteractionType::PAIR_PRODUCTION_NUCLEAR:
    return "Pair Production (Nuclear)";
  case PhotonInteractionType::PAIR_PRODUCTION_ELECTRON:
    return "Pair Production (Electronic)";
  default:
    return "Unknown Interaction";
  }
}

// Private helper methods

double PhotonPhysics::approximatePhotoelectricCrossSection(double energy, int Z) const
{
  if(energy <= 0.0) return 0.0;
  
  // Simplified approximation: σ ∝ Z^5 / E^3.5
  // Real implementation would use Scofield calculations or tabulated data
  double z_factor = std::pow(Z, 5.0);
  double energy_factor = std::pow(energy, -3.5);
  
  return 1e-24 * z_factor * energy_factor; // Approximate in cm²
}

double PhotonPhysics::approximateCoherentCrossSection(double energy, int Z) const
{
  if(energy <= 0.0) return 0.0;
  
  // Simplified approximation: σ ∝ Z² / E²
  double z_factor = Z * Z;
  double energy_factor = 1.0 / (energy * energy);
  
  return 1e-25 * z_factor * energy_factor; // Approximate in cm²
}

double PhotonPhysics::approximatePairProductionCrossSection(double energy, int Z) const
{
  if(energy < PAIR_PRODUCTION_THRESHOLD) return 0.0;
  
  // Simplified approximation: σ ∝ Z² * ln(E/m_e c²)
  double z_factor = Z * Z;
  double energy_ratio = energy / ELECTRON_REST_MASS_MEV;
  double log_factor = std::log(energy_ratio);
  
  return 1e-26 * z_factor * log_factor; // Approximate in cm²
}

double PhotonPhysics::samplePhiAngle()
{
  return 2.0 * M_PI * sampling_.getRNG().uniform();
}

Vector3D PhotonPhysics::rotateDirection(const Vector3D &original, double cos_theta, double phi)
{
  double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
  double cos_phi = std::cos(phi);
  double sin_phi = std::sin(phi);
  
  // Create local coordinate system
  Vector3D u = original;
  Vector3D v, w;
  
  // Find perpendicular vectors
  if(std::abs(u.x()) > 0.1)
  {
    v = Vector3D(0.0, 1.0, 0.0).cross(u).normalise();
  }
  else
  {
    v = Vector3D(1.0, 0.0, 0.0).cross(u).normalise();
  }
  w = u.cross(v);
  
  // Rotate in local coordinates
  return u * cos_theta + v * (sin_theta * cos_phi) + w * (sin_theta * sin_phi);
}

// Additional missing private methods
double PhotonPhysics::calculatePhotoelectricCoefficient(double energy, int Z) const
{
  return approximatePhotoelectricCrossSection(energy, Z);
}

double PhotonPhysics::calculateComptonCoefficient(double energy, int Z) const
{
  return kleinNishinaCrossSection(energy) * Z;
}

double PhotonPhysics::calculateCoherentCoefficient(double energy, int Z) const
{
  return approximateCoherentCrossSection(energy, Z);
}

double PhotonPhysics::calculatePairProductionCoefficient(double energy, int Z) const
{
  return approximatePairProductionCrossSection(energy, Z);
}