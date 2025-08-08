#include "Transport.hpp"
#include "PhotonCrossSectionDatabase.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

Transport::Transport(const Geometry &geometry, MonteCarloSampling &sampling)
    : geometry_(geometry), sampling_(sampling), 
      photon_physics_(std::make_unique<PhotonPhysics>(sampling)),  // Enable photon physics
      maximum_step_size_(1e30),
      enable_boundary_crossing_(true), enable_interaction_sampling_(true),
      results_(nullptr), collect_detailed_results_(false)
{
  // Set default cross-section function for photons
  cross_section_function_ = DefaultCrossSections::photonTotalCrossSection;
}

void Transport::trackParticle(Particle &particle)
{
  int iteration_count = 0;
  const int MAX_ITERATIONS = 1000; // Safety limit
  
  while(particle.isAlive())
  {
    iteration_count++;
    
    if (iteration_count > MAX_ITERATIONS) {
      particle.kill();
      break;
    }
    
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
      auto interaction_result = performInteraction(particle, current_material);
      
      // Handle any secondary particles created from the interaction
      if (!interaction_result.secondary_particles.empty()) {
        handleSecondaryParticles(interaction_result.secondary_particles);
      }
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
    auto interaction_result = performInteraction(particle, current_material);
    
    // Handle any secondary particles created from the interaction
    if (!interaction_result.secondary_particles.empty()) {
      handleSecondaryParticles(interaction_result.secondary_particles);
    }
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
  double final_distance = std::min(distance, maximum_step_size_);
  return final_distance;
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

PhotonInteractionResult Transport::performInteraction(Particle &particle, const Material &material)
{
  // Record interaction point
  particle.recordInteraction(particle.position());
  
  // Initialize default result
  PhotonInteractionResult result;
  result.particle_absorbed = false;
  result.particle_scattered = false;
  
  // Handle physics-specific interactions based on particle type
  if (particle.type() == ParticleType::Photon && photon_physics_) {
    // Use realistic photon physics
    result = photon_physics_->performInteraction(particle, material);
    
    // Update particle state based on interaction result
    if (result.particle_absorbed) {
      particle.absorb();
    } else if (result.particle_scattered) {
      // Update particle direction and energy
      particle.setDirection(result.new_direction);
      particle.setEnergy(result.new_energy);
    }
    
    // Record interaction for statistics if results collection is enabled
    if (results_) {
      recordInteraction(particle, material, result);
    }
    
  } else {
    // Fallback for other particle types or when photon physics is not available
    // Simple absorption for now
    particle.absorb();
    result.particle_absorbed = true;
    result.interaction_type = PhotonInteractionType::PHOTOELECTRIC_ABSORPTION;
    
    // Record energy deposition
    if (results_) {
      recordEnergyDeposition(material, particle.energy(), PhotonInteractionType::PHOTOELECTRIC_ABSORPTION);
    }
  }
  
  return result;
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

// =============================================================================
// New helper method implementations for Day 11 integration
// =============================================================================

void Transport::handleSecondaryParticles(const std::vector<Particle>& secondaries)
{
  // For now, we'll just track secondary particles in the results
  // In a full implementation, you might want to queue them for transport
  // or handle them recursively based on their importance
  
  if (results_) {
    for (const auto& secondary : secondaries) {
      // Record secondary particle creation
      // This would typically involve adding them to a particle queue
      // For statistical purposes, we can record their energy contributions
      
      if (collect_detailed_results_) {
        // Log secondary particle creation
        // In practice, you might transport these particles as well
        // depending on their type and energy threshold
      }
    }
  }
}

void Transport::recordInteraction(const Particle& particle, const Material& material,
                                const PhotonInteractionResult& result)
{
  if (!results_) return;
  
  // Calculate energy deposited in this interaction
  double energy_deposited = 0.0;
  
  if (result.particle_absorbed) {
    // All particle energy is deposited
    energy_deposited = particle.energy();
  } else if (result.particle_scattered) {
    // Energy difference is deposited (recoil electron energy, etc.)
    energy_deposited = particle.energy() - result.new_energy;
  }
  
  // Record energy deposition
  recordEnergyDeposition(material, energy_deposited, result.interaction_type);
  
  // Create detailed interaction record if requested
  if (collect_detailed_results_) {
    InteractionRecord record = createInteractionRecord(particle, material, result);
    results_->recordInteraction(record);
  }
}

void Transport::recordEnergyDeposition(const Material& material, double energy_deposited,
                                     PhotonInteractionType interaction_type)
{
  if (!results_ || energy_deposited <= 0.0) return;
  
  // Record in simulation results
  results_->recordEnergyDeposition(std::string(material.name()), energy_deposited, interaction_type);
}

InteractionRecord Transport::createInteractionRecord(const Particle& particle, 
                                                   const Material& material,
                                                   const PhotonInteractionResult& result) const
{
  InteractionRecord record;
  record.position = particle.position();
  record.energy_before = particle.energy();
  record.energy_after = result.new_energy;
  record.energy_deposited = particle.energy() - result.new_energy;
  record.interaction_type = result.interaction_type;
  record.material = material;
  record.direction_before = particle.direction();
  record.direction_after = result.new_direction;
  record.generation = particle.generation();
  record.particle_absorbed = result.particle_absorbed;
  record.secondary_particles = result.secondary_particles;
  
  return record;
}

// Default cross-section implementations
namespace DefaultCrossSections
{

double photonTotalCrossSection(const Material &material, double energy)
{
  // Accurate photon cross-section using NIST-based PhotonCrossSectionDatabase
  if(energy <= 0.0)
  {
    return 0.0;
  }

  // Convert energy to keV (assuming input is in MeV)
  double energy_keV = energy * 1000.0;
  
  // Use Material's integrated photon cross-section methods
  return material.getLinearAttenuationCoefficient(energy_keV); // cm⁻¹
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