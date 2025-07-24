#include "MonteCarloSampling.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

// Physical constants
namespace PhysicalConstants
{
constexpr double PI = 3.14159265358979323846;
constexpr double TWO_PI = 2.0 * PI;
constexpr double ELECTRON_REST_MASS = 0.510999; // MeV
} // namespace PhysicalConstants

// ============================================================================
// EnergySpectrum Factory Methods
// ============================================================================

EnergySpectrum EnergySpectrum::monoenergetic(double energy)
{
  EnergySpectrum spectrum;
  spectrum.type = SpectrumType::Monoenergetic;
  spectrum.characteristic_energy = energy;
  spectrum.min_energy = energy;
  spectrum.max_energy = energy;
  return spectrum;
}

EnergySpectrum EnergySpectrum::uniform(double min_energy, double max_energy)
{
  if(min_energy >= max_energy)
  {
    throw std::invalid_argument(
        "Minimum energy must be less than maximum energy");
  }

  EnergySpectrum spectrum;
  spectrum.type = SpectrumType::Uniform;
  spectrum.min_energy = min_energy;
  spectrum.max_energy = max_energy;
  return spectrum;
}

EnergySpectrum EnergySpectrum::maxwell(double temperature)
{
  if(temperature <= 0.0)
  {
    throw std::invalid_argument("Temperature must be positive");
  }

  EnergySpectrum spectrum;
  spectrum.type = SpectrumType::Maxwell;
  spectrum.characteristic_energy = temperature;
  spectrum.parameter = temperature;
  spectrum.max_energy = 10.0 * temperature; // Reasonable cutoff
  return spectrum;
}

EnergySpectrum EnergySpectrum::exponential(double characteristic_energy)
{
  if(characteristic_energy <= 0.0)
  {
    throw std::invalid_argument("Characteristic energy must be positive");
  }

  EnergySpectrum spectrum;
  spectrum.type = SpectrumType::Exponential;
  spectrum.characteristic_energy = characteristic_energy;
  spectrum.max_energy = 10.0 * characteristic_energy; // Reasonable cutoff
  return spectrum;
}

EnergySpectrum EnergySpectrum::powerLaw(double alpha, double min_energy,
                                        double max_energy)
{
  if(min_energy >= max_energy)
  {
    throw std::invalid_argument(
        "Minimum energy must be less than maximum energy");
  }
  if(std::abs(alpha + 1.0) < 1e-10)
  {
    throw std::invalid_argument("Alpha cannot be -1");
  }

  EnergySpectrum spectrum;
  spectrum.type = SpectrumType::PowerLaw;
  spectrum.parameter = alpha;
  spectrum.min_energy = min_energy;
  spectrum.max_energy = max_energy;
  return spectrum;
}

EnergySpectrum EnergySpectrum::discrete(const std::vector<double> &energies,
                                        const std::vector<double> &weights)
{
  if(energies.empty() || weights.empty())
  {
    throw std::invalid_argument("Energies and weights cannot be empty");
  }
  if(energies.size() != weights.size())
  {
    throw std::invalid_argument("Energies and weights must have same size");
  }

  EnergySpectrum spectrum;
  spectrum.type = SpectrumType::Discrete;
  spectrum.discrete_energies = energies;
  spectrum.discrete_weights = weights;
  spectrum.min_energy = *std::min_element(energies.begin(), energies.end());
  spectrum.max_energy = *std::max_element(energies.begin(), energies.end());
  return spectrum;
}

// ============================================================================
// AngularDistribution Factory Methods
// ============================================================================

AngularDistribution AngularDistribution::isotropic()
{
  AngularDistribution distribution;
  distribution.type = DirectionType::Isotropic;
  return distribution;
}

AngularDistribution AngularDistribution::beam(const Vector3D &direction)
{
  auto normalised = direction.tryNormalise();
  if(!normalised)
  {
    throw std::invalid_argument("Beam direction cannot be zero vector");
  }

  AngularDistribution distribution;
  distribution.type = DirectionType::Beam;
  distribution.reference_direction = *normalised;
  return distribution;
}

AngularDistribution AngularDistribution::cone(const Vector3D &axis,
                                              double half_angle)
{
  if(half_angle < 0.0 || half_angle > PhysicalConstants::PI)
  {
    throw std::invalid_argument("Cone half-angle must be between 0 and π");
  }

  auto normalised = axis.tryNormalise();
  if(!normalised)
  {
    throw std::invalid_argument("Cone axis cannot be zero vector");
  }

  AngularDistribution distribution;
  distribution.type = DirectionType::Cone;
  distribution.reference_direction = *normalised;
  distribution.cone_angle = half_angle;
  return distribution;
}

AngularDistribution AngularDistribution::cosineWeighted(const Vector3D &normal)
{
  auto normalised = normal.tryNormalise();
  if(!normalised)
  {
    throw std::invalid_argument("Normal vector cannot be zero vector");
  }

  AngularDistribution distribution;
  distribution.type = DirectionType::CosineWeighted;
  distribution.reference_direction = *normalised;
  return distribution;
}

// ============================================================================
// MonteCarloSampling Implementation
// ============================================================================

MonteCarloSampling::MonteCarloSampling(RandomNumberGenerator &rng) : rng_(rng)
{}

// ============================================================================
// EXPONENTIAL SAMPLING FOR INTERACTION DISTANCES
// ============================================================================

double MonteCarloSampling::sampleInteractionDistance(
    double linear_attenuation_coefficient) const
{
  if(linear_attenuation_coefficient <= 0.0)
  {
    throw std::invalid_argument(
        "Linear attenuation coefficient must be positive");
  }

  return rng_.exponential(linear_attenuation_coefficient);
}

std::optional<double>
MonteCarloSampling::tryExponentialSampling(double lambda) const
{
  return rng_.tryExponential(lambda);
}

double MonteCarloSampling::sampleInteractionDistance(
    const Material &material, double particle_energy,
    std::function<double(const Material &, double)> cross_section_function)
    const
{
  double cross_section = cross_section_function(material, particle_energy);
  if(cross_section <= 0.0)
  {
    throw std::invalid_argument("Cross-section must be positive");
  }

  // Convert microscopic cross-section to macroscopic
  double atom_density = material.getAtomDensity();
  double macroscopic_cross_section = cross_section * atom_density;

  return sampleInteractionDistance(macroscopic_cross_section);
}

double MonteCarloSampling::sampleTotalInteractionDistance(
    double total_cross_section) const
{
  return sampleInteractionDistance(total_cross_section);
}

// ============================================================================
// UNIFORM SPHERICAL DIRECTION SAMPLING
// ============================================================================

Vector3D MonteCarloSampling::sampleIsotropicDirection() const
{
  return rng_.isotropicDirection();
}

Vector3D MonteCarloSampling::sampleDirection(
    const AngularDistribution &distribution) const
{
  switch(distribution.type)
  {
  case DirectionType::Isotropic:
    return sampleIsotropicDirection();

  case DirectionType::Beam:
    return distribution.reference_direction;

  case DirectionType::Cone:
    return sampleConeDirection(distribution.reference_direction,
                               distribution.cone_angle);

  case DirectionType::CosineWeighted:
    return sampleCosineDirection(distribution.reference_direction);

  case DirectionType::Custom:
    if(distribution.custom_sampler)
    {
      return distribution.custom_sampler();
    }
    throw std::invalid_argument("Custom sampler not provided");

  default:
    throw std::invalid_argument("Unknown direction type");
  }
}

Vector3D MonteCarloSampling::sampleConeDirection(const Vector3D &axis,
                                                 double half_angle) const
{
  if(half_angle < 0.0)
  {
    throw std::invalid_argument("Cone half-angle cannot be negative");
  }

  if(half_angle == 0.0)
  {
    return axis; // Degenerate case - pencil beam
  }

  // Sample uniformly within cone
  double cos_half_angle = std::cos(half_angle);
  double cos_theta = cos_half_angle + rng_.uniform() * (1.0 - cos_half_angle);
  double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
  double phi = PhysicalConstants::TWO_PI * rng_.uniform();

  // Local coordinates (cone axis is Z)
  Vector3D local_direction(sin_theta * std::cos(phi), sin_theta * std::sin(phi),
                           cos_theta);

  // Transform to global coordinates
  return MonteCarloUtils::localToGlobal(local_direction, axis);
}

Vector3D MonteCarloSampling::sampleCosineDirection(const Vector3D &normal) const
{
  // Sample cosine-weighted hemisphere
  double cos_theta = std::sqrt(rng_.uniform()); // √ξ for cosine weighting
  double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
  double phi = PhysicalConstants::TWO_PI * rng_.uniform();

  Vector3D local_direction(sin_theta * std::cos(phi), sin_theta * std::sin(phi),
                           cos_theta);

  return MonteCarloUtils::localToGlobal(local_direction, normal);
}

Vector3D
MonteCarloSampling::sampleHemisphereDirection(const Vector3D &normal) const
{
  // Sample uniformly in hemisphere
  double cos_theta = rng_.uniform(); // Uniform in cos(θ)
  double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
  double phi = PhysicalConstants::TWO_PI * rng_.uniform();

  Vector3D local_direction(sin_theta * std::cos(phi), sin_theta * std::sin(phi),
                           cos_theta);

  return MonteCarloUtils::localToGlobal(local_direction, normal);
}

// ============================================================================
// ENERGY DISTRIBUTION SAMPLING
// ============================================================================

double MonteCarloSampling::sampleEnergy(const EnergySpectrum &spectrum) const
{
  switch(spectrum.type)
  {
  case SpectrumType::Monoenergetic:
    return sampleMonoenergeticEnergy(spectrum.characteristic_energy);

  case SpectrumType::Uniform:
    return sampleUniformEnergy(spectrum.min_energy, spectrum.max_energy);

  case SpectrumType::Maxwell:
    return sampleMaxwellEnergy(spectrum.characteristic_energy);

  case SpectrumType::Exponential:
    return sampleExponentialEnergy(spectrum.characteristic_energy);

  case SpectrumType::PowerLaw:
    return samplePowerLawEnergy(spectrum.parameter, spectrum.min_energy,
                                spectrum.max_energy);

  case SpectrumType::Discrete:
    return sampleDiscreteEnergy(spectrum.discrete_energies,
                                spectrum.discrete_weights);

  case SpectrumType::Custom:
    if(spectrum.custom_pdf)
    {
      // Use rejection sampling for custom PDF
      double max_energy = spectrum.max_energy;
      double min_energy = spectrum.min_energy;
      double max_pdf = 1.0; // Assume normalised PDF

      for(int attempt = 0; attempt < 10000; ++attempt)
      {
        double energy = rng_.uniform(min_energy, max_energy);
        double pdf_value = spectrum.custom_pdf(energy);
        if(rng_.uniform() * max_pdf <= pdf_value)
        {
          return energy;
        }
      }
      throw std::runtime_error(
          "Failed to sample from custom energy distribution");
    }
    throw std::invalid_argument("Custom PDF not provided");

  default:
    throw std::invalid_argument("Unknown spectrum type");
  }
}

double MonteCarloSampling::sampleMonoenergeticEnergy(double energy) const
{
  if(energy <= 0.0)
  {
    throw std::invalid_argument("Energy must be positive");
  }
  return energy;
}

double MonteCarloSampling::sampleUniformEnergy(double min_energy,
                                               double max_energy) const
{
  if(min_energy >= max_energy)
  {
    throw std::invalid_argument(
        "Minimum energy must be less than maximum energy");
  }
  return rng_.uniform(min_energy, max_energy);
}

double MonteCarloSampling::sampleMaxwellEnergy(double temperature) const
{
  if(temperature <= 0.0)
  {
    throw std::invalid_argument("Temperature must be positive");
  }
  return rng_.maxwellBoltzmann(temperature);
}

double
MonteCarloSampling::sampleExponentialEnergy(double characteristic_energy) const
{
  if(characteristic_energy <= 0.0)
  {
    throw std::invalid_argument("Characteristic energy must be positive");
  }
  return rng_.exponential(1.0 / characteristic_energy);
}

double MonteCarloSampling::samplePowerLawEnergy(double alpha, double min_energy,
                                                double max_energy) const
{
  if(min_energy >= max_energy)
  {
    throw std::invalid_argument(
        "Minimum energy must be less than maximum energy");
  }
  return rng_.powerLaw(alpha, min_energy, max_energy);
}

double MonteCarloSampling::sampleDiscreteEnergy(
    const std::vector<double> &energies,
    const std::vector<double> &weights) const
{
  if(energies.empty() || weights.empty())
  {
    throw std::invalid_argument("Energies and weights cannot be empty");
  }
  if(energies.size() != weights.size())
  {
    throw std::invalid_argument("Energies and weights must have same size");
  }

  std::size_t index = rng_.discreteSample(weights);
  return energies[index];
}

// ============================================================================
// ADVANCED MONTE CARLO TECHNIQUES
// ============================================================================

double
MonteCarloSampling::sampleScatteringAngle(std::function<double(double)> pdf,
                                          double max_pdf) const
{
  // Rejection sampling for arbitrary PDF
  for(int attempt = 0; attempt < 10000; ++attempt)
  {
    double angle = rng_.uniform(0.0, PhysicalConstants::PI);
    double pdf_value = pdf(angle);
    if(rng_.uniform() * max_pdf <= pdf_value)
    {
      return angle;
    }
  }
  throw std::runtime_error("Failed to sample scattering angle");
}

double MonteCarloSampling::sampleAzimuthalAngle() const
{
  return rng_.uniform(0.0, PhysicalConstants::TWO_PI);
}

Vector3D MonteCarloSampling::applyScattering(const Vector3D &incident_direction,
                                             double polar_angle,
                                             double azimuthal_angle) const
{
  // Create coordinate system with incident direction as Z-axis
  auto [x_axis, y_axis] =
      MonteCarloUtils::generateOrthonormalBasis(incident_direction);

  // Calculate scattered direction in local coordinates
  double sin_theta = std::sin(polar_angle);
  double cos_theta = std::cos(polar_angle);
  double sin_phi = std::sin(azimuthal_angle);
  double cos_phi = std::cos(azimuthal_angle);

  Vector3D scattered_local(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);

  // Transform back to global coordinates
  return scattered_local.x() * x_axis + scattered_local.y() * y_axis +
         scattered_local.z() * incident_direction;
}

std::size_t MonteCarloSampling::sampleInteractionType(
    const std::vector<double> &cross_sections) const
{
  return rng_.discreteSample(cross_sections);
}

bool MonteCarloSampling::russianRoulette(double survival_probability) const
{
  return rng_.russianRoulette(survival_probability);
}

std::size_t MonteCarloSampling::particleSplitting(double splitting_factor) const
{
  return rng_.particleSplitting(splitting_factor);
}

// ============================================================================
// PHYSICS-SPECIFIC SAMPLING
// ============================================================================

double MonteCarloSampling::sampleComptonAngle(double incident_energy) const
{
  return rng_.comptonScatteringAngle(incident_energy);
}

double MonteCarloSampling::sampleComptonEnergy(double incident_energy,
                                               double scattering_angle) const
{
  // Compton formula: E' = E / (1 + (E/m_e c²)(1 - cos θ))
  double alpha = incident_energy / PhysicalConstants::ELECTRON_REST_MASS;
  double energy_ratio =
      1.0 / (1.0 + alpha * (1.0 - std::cos(scattering_angle)));
  return incident_energy * energy_ratio;
}

Vector3D MonteCarloSampling::samplePhotoelectronDirection() const
{
  // Simple isotropic approximation for photoelectron direction
  return sampleIsotropicDirection();
}

std::pair<Vector3D, Vector3D>
MonteCarloSampling::samplePairProductionDirections(double photon_energy) const
{
  // Simplified pair production - electrons share energy roughly equally
  // and are emitted in roughly forward direction

  double electron_angle = rng_.normal(0.0, 0.1); // Small forward scatter
  double positron_angle = rng_.normal(0.0, 0.1);

  // Sample azimuthal angles
  double electron_phi = sampleAzimuthalAngle();
  double positron_phi = sampleAzimuthalAngle();

  Vector3D electron_dir =
      MonteCarloUtils::polarToCartesian(electron_angle, electron_phi);
  Vector3D positron_dir =
      MonteCarloUtils::polarToCartesian(positron_angle, positron_phi);

  return std::make_pair(electron_dir, positron_dir);
}

// ============================================================================
// UTILITY METHODS
// ============================================================================

bool MonteCarloSampling::validateEnergySpectrum(
    const EnergySpectrum &spectrum) const
{
  if(spectrum.min_energy < 0.0)
    return false;
  if(spectrum.max_energy <= spectrum.min_energy)
    return false;
  if(spectrum.characteristic_energy <= 0.0)
    return false;

  if(spectrum.type == SpectrumType::Discrete)
  {
    if(spectrum.discrete_energies.empty() || spectrum.discrete_weights.empty())
      return false;
    if(spectrum.discrete_energies.size() != spectrum.discrete_weights.size())
      return false;

    for(double energy : spectrum.discrete_energies)
    {
      if(energy <= 0.0)
        return false;
    }
    for(double weight : spectrum.discrete_weights)
    {
      if(weight < 0.0)
        return false;
    }
  }

  return true;
}

bool MonteCarloSampling::validateAngularDistribution(
    const AngularDistribution &distribution) const
{
  if(!distribution.reference_direction.isUnit(1e-10))
    return false;

  if(distribution.type == DirectionType::Cone)
  {
    if(distribution.cone_angle < 0.0 ||
       distribution.cone_angle > PhysicalConstants::PI)
    {
      return false;
    }
  }

  return true;
}

std::vector<double> MonteCarloSampling::sampleInteractionDistances(
    double linear_attenuation_coefficient, std::size_t count) const
{
  return rng_.exponentialBatch(count, linear_attenuation_coefficient);
}

std::vector<Vector3D>
MonteCarloSampling::sampleIsotropicDirections(std::size_t count) const
{
  return rng_.isotropicDirectionBatch(count);
}

std::vector<double>
MonteCarloSampling::sampleEnergies(const EnergySpectrum &spectrum,
                                   std::size_t count) const
{
  std::vector<double> energies;
  energies.reserve(count);

  for(std::size_t i = 0; i < count; ++i)
  {
    energies.push_back(sampleEnergy(spectrum));
  }

  return energies;
}

// ============================================================================
// CONVENIENCE FUNCTIONS
// ============================================================================

namespace MonteCarloUtils
{
double kleinNishinaCrossSection(double incident_energy, double scattering_angle)
{
  double alpha = incident_energy / PhysicalConstants::ELECTRON_REST_MASS;
  double cos_theta = std::cos(scattering_angle);
  double energy_ratio = 1.0 / (1.0 + alpha * (1.0 - cos_theta));

  return energy_ratio * energy_ratio *
         (energy_ratio + 1.0 / energy_ratio - 1.0 + cos_theta * cos_theta);
}

Vector3D polarToCartesian(double theta, double phi)
{
  double sin_theta = std::sin(theta);
  return Vector3D(sin_theta * std::cos(phi), sin_theta * std::sin(phi),
                  std::cos(theta));
}

std::tuple<double, double> cartesianToPolar(const Vector3D &direction)
{
  double theta = std::acos(std::max(-1.0, std::min(1.0, direction.z())));
  double phi = std::atan2(direction.y(), direction.x());
  if(phi < 0.0)
    phi += PhysicalConstants::TWO_PI;
  return std::make_tuple(theta, phi);
}

Vector3D localToGlobal(const Vector3D &local_vector, const Vector3D &local_z)
{
  auto [x_axis, y_axis] = generateOrthonormalBasis(local_z);

  return local_vector.x() * x_axis + local_vector.y() * y_axis +
         local_vector.z() * local_z;
}

std::pair<Vector3D, Vector3D> generateOrthonormalBasis(const Vector3D &z)
{
  Vector3D x_axis;

  // Choose initial vector that's not parallel to z
  if(std::abs(z.x()) < 0.9)
  {
    x_axis = Vector3D(1.0, 0.0, 0.0);
  }
  else
  {
    x_axis = Vector3D(0.0, 1.0, 0.0);
  }

  // Gram-Schmidt orthogonalisation
  x_axis = x_axis - z * z.dot(x_axis);
  x_axis = x_axis.normalise();

  Vector3D y_axis = z.cross(x_axis);

  return std::make_pair(x_axis, y_axis);
}
} // namespace MonteCarloUtils