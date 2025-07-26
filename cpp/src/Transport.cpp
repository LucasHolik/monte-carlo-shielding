#include "Transport.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

Transport::Transport(const Geometry &geometry, MonteCarloSampling &sampling)
    : geometry_(geometry), sampling_(sampling), maximum_step_size_(1e30),
      enable_boundary_crossing_(true), enable_interaction_sampling_(true)
{
  // Set default cross-section function for photons
  cross_section_function_ = DefaultCrossSections::photonTotalCrossSection;
}

void Transport::trackParticle(Particle &particle)
{
  while(particle.isAlive())
  {
    // Sample distance to next interaction
    auto current_material = getCurrentMaterial(particle);
    double distance = sampleDistanceToInteraction(particle, current_material);

    // Check if particle hits boundary before interaction
    if(hitBoundary(particle, distance))
    {
      handleBoundary(particle);
    }
    else
    {
      // Move particle and perform interaction
      moveParticle(particle, distance);
      performInteraction(particle, current_material);
    }

    // Safety check to prevent infinite loops
    if(!isValidTransportState(particle))
    {
      particle.kill();
      break;
    }
  }
}

TransportStep Transport::transportParticleStep(Particle &particle)
{
  if(!particle.isAlive())
  {
    return TransportStep(0.0, particle.position(), Material::createVacuum());
  }

  auto current_material = getCurrentMaterial(particle);
  double distance = sampleDistanceToInteraction(particle, current_material);

  TransportStep step(distance, particle.position(), current_material);

  if(hitBoundary(particle, distance))
  {
    step.hit_boundary = true;
    step.boundary_info = findNextBoundary(particle);
    handleBoundary(particle);
  }
  else
  {
    moveParticle(particle, distance);
    step.interaction_occurred = true;
    performInteraction(particle, current_material);
  }

  step.final_position = particle.position();
  return step;
}

double Transport::sampleDistanceToInteraction(const Particle &particle,
                                              const Material &material)
{
  if(!enable_interaction_sampling_)
  {
    return maximum_step_size_;
  }

  // Calculate linear attenuation coefficient
  double mu =
      calculateLinearAttenuationCoefficient(material, particle.energy());

  if(mu <= 0.0)
  {
    return maximum_step_size_; // No interaction possible
  }

  // Sample exponential distribution
  double distance = sampling_.sampleInteractionDistance(mu);

  // Apply maximum step size limit
  return std::min(distance, maximum_step_size_);
}

void Transport::moveParticle(Particle &particle, double distance)
{
  if(distance <= 0.0)
  {
    return;
  }

  // Calculate new position
  Vector3D displacement = particle.direction() * distance;
  Vector3D new_position = particle.position() + displacement;

  // Update particle position
  particle.setPosition(new_position);
}

bool Transport::hitBoundary(const Particle &particle, double proposed_distance)
{
  if(!enable_boundary_crossing_)
  {
    return false;
  }

  double distance_to_boundary = getDistanceToBoundary(particle);
  return distance_to_boundary < proposed_distance;
}

void Transport::handleBoundary(Particle &particle)
{
  auto boundary_intersection = findNextBoundary(particle);

  if(!boundary_intersection)
  {
    // No boundary found - particle escapes
    particle.escape();
    return;
  }

  // Move particle to boundary
  double distance_to_boundary = boundary_intersection->distance;
  moveParticle(particle, distance_to_boundary);

  // Handle geometry transition
  handleGeometryTransition(particle, *boundary_intersection);

  // Check if particle is now outside all geometry
  if(!geometry_.isInsideAnyShape(particle.position()))
  {
    particle.escape();
  }
}

void Transport::performInteraction(Particle &particle, const Material &material)
{
  // Record interaction point
  particle.recordInteraction(particle.position());

  // For now, implement simple absorption
  // In a full implementation, this would sample interaction type
  // and handle scattering, secondary particle production, etc.
  particle.absorb();
}

std::vector<TransportStep> Transport::getParticleTrack(Particle particle,
                                                       double max_distance)
{
  std::vector<TransportStep> track;
  double total_distance = 0.0;

  while(particle.isAlive() && total_distance < max_distance)
  {
    auto step = transportParticleStep(particle);
    track.push_back(step);

    total_distance += step.distance_travelled;

    if(!isValidTransportState(particle))
    {
      break;
    }
  }

  return track;
}

bool Transport::willParticleEscape(const Particle &particle)
{
  return geometry_.willParticleEscape(particle.position(),
                                      particle.direction());
}

std::optional<Vector3D> Transport::findEscapePoint(const Particle &particle)
{
  auto intersection = geometry_.findClosestIntersection(particle.position(),
                                                        particle.direction());

  if(!intersection)
  {
    return std::nullopt; // Already outside or parallel to boundaries
  }

  return intersection->intersection_point;
}

Material Transport::getCurrentMaterial(const Particle &particle) const
{
  return geometry_.getMaterialAtPoint(particle.position());
}

std::optional<GeometryIntersection>
Transport::findNextBoundary(const Particle &particle) const
{
  return geometry_.findClosestIntersection(particle.position(),
                                           particle.direction());
}

double Transport::getDistanceToBoundary(const Particle &particle) const
{
  auto boundary = findNextBoundary(particle);
  if(boundary)
  {
    return boundary->distance;
  }
  return std::numeric_limits<double>::infinity();
}

bool Transport::isValidTransportState(const Particle &particle) const
{
  // Check for valid position
  auto pos = particle.position();
  if(!pos.isFinite())
  {
    return false;
  }

  // Check for valid direction
  auto dir = particle.direction();
  if(!dir.isFinite() || std::abs(dir.magnitude() - 1.0) > 1e-10)
  {
    return false;
  }

  // Check for valid energy
  if(particle.energy() < 0.0 || !std::isfinite(particle.energy()))
  {
    return false;
  }

  return true;
}

std::string Transport::getTransportDiagnostics(const Particle &particle) const
{
  std::ostringstream oss;
  oss << "Transport Diagnostics:\n";
  oss << "  Position: " << particle.position().toString() << "\n";
  oss << "  Direction: " << particle.direction().toString() << "\n";
  oss << "  Energy: " << particle.energy() << " MeV\n";
  oss << "  State: " << particle.getStateName() << "\n";
  oss << "  Distance to boundary: " << getDistanceToBoundary(particle)
      << " cm\n";

  auto material = getCurrentMaterial(particle);
  oss << "  Current material: " << material.name() << "\n";

  return oss.str();
}

double
Transport::calculateLinearAttenuationCoefficient(const Material &material,
                                                 double energy) const
{
  if(!cross_section_function_)
  {
    return 0.0; // No interactions
  }

  // Get macroscopic cross-section (cm⁻¹)
  double cross_section = cross_section_function_(material, energy);

  // Convert to linear attenuation coefficient
  // For now, assume cross_section is already in cm⁻¹
  return cross_section;
}

void Transport::validateParticleState(const Particle &particle) const
{
  if(!isValidTransportState(particle))
  {
    throw std::runtime_error(
        "Invalid particle state detected during transport");
  }
}

void Transport::handleGeometryTransition(
    Particle &particle, const GeometryIntersection &intersection)
{
  // Set particle position exactly on boundary
  particle.setPosition(intersection.intersection_point);

  // In a full implementation, this would handle:
  // - Reflection/transmission at interfaces
  // - Energy threshold checks
  // - Material property changes
  // For now, particle continues in same direction
}

// Default cross-section implementations
namespace DefaultCrossSections
{

double photonTotalCrossSection(const Material &material, double energy)
{
  // Simplified photon cross-section (placeholder)
  // In reality, would use tabulated data or analytical approximations

  if(energy <= 0.0)
  {
    return 0.0;
  }

  // Simple approximation: μ = ρ * (μ/ρ)
  // Where (μ/ρ) is mass attenuation coefficient
  double density = material.density(); // g/cm³

  // Very rough approximation for demonstration
  // Real implementation would use NIST data
  double mass_attenuation = 0.1 / energy; // cm²/g (simplified)

  return density * mass_attenuation; // cm⁻¹
}

double neutronTotalCrossSection(const Material &material, double energy)
{
  // Placeholder neutron cross-section
  if(energy <= 0.0)
  {
    return 0.0;
  }

  // Simple 1/v dependence for thermal neutrons
  double reference_energy = 0.025;      // eV
  double reference_cross_section = 0.1; // cm⁻¹

  return reference_cross_section * std::sqrt(reference_energy / energy);
}

double electronCrossSection(const Material &material, double energy)
{
  // Placeholder electron cross-section
  if(energy <= 0.0)
  {
    return 0.0;
  }

  // Simple energy dependence
  return material.density() * 0.01 / energy; // cm⁻¹
}

} // namespace DefaultCrossSections