#include "../include/ParticleSource.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// =============================================================================
// SourceGeometry Factory Methods
// =============================================================================

SourceGeometry SourceGeometry::createPoint(const Vector3D& position)
{
  SourceGeometry geom;
  geom.type = SourceGeometryType::Point;
  geom.center = position;
  return geom;
}

SourceGeometry SourceGeometry::createPlane(const Vector3D& center, const Vector3D& normal,
                                         double width, double height)
{
  SourceGeometry geom;
  geom.type = SourceGeometryType::Plane;
  geom.center = center;
  geom.direction = normal.normalise();
  geom.parameter1 = width;
  geom.parameter2 = height;
  return geom;
}

SourceGeometry SourceGeometry::createSphere(const Vector3D& center, double radius, bool surface_only)
{
  SourceGeometry geom;
  geom.type = SourceGeometryType::Sphere;
  geom.center = center;
  geom.parameter1 = radius;
  geom.parameter2 = surface_only ? 1.0 : 0.0;  // Use parameter2 as flag
  return geom;
}

SourceGeometry SourceGeometry::createCylinder(const Vector3D& center, const Vector3D& axis,
                                            double radius, double height, bool surface_only)
{
  SourceGeometry geom;
  geom.type = SourceGeometryType::Cylinder;
  geom.center = center;
  geom.direction = axis.normalise();
  geom.parameter1 = radius;
  geom.parameter2 = height;
  geom.parameter3 = surface_only ? 1.0 : 0.0;
  return geom;
}

SourceGeometry SourceGeometry::createRectangle(const Vector3D& corner, const Vector3D& side1,
                                              const Vector3D& side2)
{
  SourceGeometry geom;
  geom.type = SourceGeometryType::Rectangle;
  geom.center = corner;
  geom.direction = side1.normalise();  // Store first side direction
  geom.parameter1 = side1.magnitude(); // First side length
  geom.parameter2 = side2.magnitude(); // Second side length
  geom.parameter3 = side1.dot(side2) / (side1.magnitude() * side2.magnitude()); // cos(angle)
  return geom;
}

SourceGeometry SourceGeometry::createCone(const Vector3D& apex, const Vector3D& axis,
                                         double base_radius, double height)
{
  SourceGeometry geom;
  geom.type = SourceGeometryType::Cone;
  geom.center = apex;
  geom.direction = axis.normalise();
  geom.parameter1 = base_radius;
  geom.parameter2 = height;
  return geom;
}

// =============================================================================
// IntensityDistribution Factory Methods  
// =============================================================================

IntensityDistribution IntensityDistribution::createUniform()
{
  IntensityDistribution dist;
  dist.type = IntensityDistributionType::Uniform;
  return dist;
}

IntensityDistribution IntensityDistribution::createGaussian(const Vector3D& center, double sigma)
{
  IntensityDistribution dist;
  dist.type = IntensityDistributionType::Gaussian;
  dist.center = center;
  dist.parameter1 = sigma;
  return dist;
}

IntensityDistribution IntensityDistribution::createExponential(const Vector3D& center, double lambda)
{
  IntensityDistribution dist;
  dist.type = IntensityDistributionType::Exponential;
  dist.center = center;
  dist.parameter1 = lambda;
  return dist;
}

// =============================================================================
// SourceConfiguration Methods
// =============================================================================

bool SourceConfiguration::isValid() const
{
  // Check basic parameters
  if (activity <= 0.0 || total_particles <= 0 || simulation_time <= 0.0) {
    return false;
  }
  
  // Check energy spectrum validity (basic checks)
  if (energy_spectrum.min_energy < 0.0 || energy_spectrum.max_energy < energy_spectrum.min_energy) {
    return false;
  }
  
  // Check geometry parameters
  switch (geometry.type) {
    case SourceGeometryType::Sphere:
      if (geometry.parameter1 <= 0.0) return false;
      break;
    case SourceGeometryType::Cylinder:
      if (geometry.parameter1 <= 0.0 || geometry.parameter2 <= 0.0) return false;
      break;
    case SourceGeometryType::Plane:
      if (geometry.parameter1 <= 0.0 || geometry.parameter2 <= 0.0) return false;
      break;
    case SourceGeometryType::Rectangle:
      if (geometry.parameter1 <= 0.0 || geometry.parameter2 <= 0.0) return false;
      break;
    case SourceGeometryType::Cone:
      if (geometry.parameter1 <= 0.0 || geometry.parameter2 <= 0.0) return false;
      break;
    default:
      break;
  }
  
  return true;
}

std::string SourceConfiguration::toString() const
{
  std::stringstream ss;
  ss << "Source: " << name << "\n";
  ss << "  Particle Type: " << static_cast<int>(particle_type) << "\n";
  ss << "  Activity: " << activity << " particles/s\n";
  ss << "  Total Particles: " << total_particles << "\n";
  ss << "  Energy Range: " << energy_spectrum.min_energy << " - " << energy_spectrum.max_energy << " MeV\n";
  ss << "  Position: (" << geometry.center.x() << ", " << geometry.center.y() << ", " << geometry.center.z() << ")\n";
  return ss.str();
}

// =============================================================================
// ParticleSource Implementation
// =============================================================================

ParticleSource::ParticleSource(MonteCarloSampling& sampling)
  : sampling_(sampling)
{
  // Set up default configuration
  config_.name = "Default Point Source";
  config_.particle_type = ParticleType::Photon;
  config_.energy_spectrum = EnergySpectrum::monoenergetic(1.0); // 1 MeV
  config_.angular_distribution = AngularDistribution::isotropic();
  config_.geometry = SourceGeometry::createPoint(Vector3D(0, 0, 0));
  config_.intensity = IntensityDistribution::createUniform();
  config_.activity = 1.0;
  config_.total_particles = 1000;
  config_.simulation_time = 1.0;
}

ParticleSource::ParticleSource(MonteCarloSampling& sampling, const SourceConfiguration& config)
  : sampling_(sampling), config_(config)
{
  validateAndThrow();
}

void ParticleSource::setConfiguration(const SourceConfiguration& config)
{
  config_ = config;
  validateAndThrow();
}

SourceParticle ParticleSource::generateParticle()
{
  // Sample source position
  Vector3D position = sampleSourcePosition();
  
  // Sample direction based on angular distribution
  Vector3D direction = sampleSourceDirection(position);
  
  // Sample energy from spectrum
  double energy = sampleSourceEnergy();
  
  // Sample intensity weight
  double intensity_weight = sampleIntensityWeight(position);
  
  // Sample emission time
  double emission_time = sampleEmissionTime();
  
  // Create particle with sampled properties
  Particle particle(config_.particle_type, position, direction, energy, 
                   config_.base_statistical_weight * intensity_weight);
  
  // Create source particle record
  SourceParticle source_particle(particle, emission_time);
  source_particle.emission_position = position;
  source_particle.emission_direction = direction;
  source_particle.emission_energy = energy;
  source_particle.statistical_weight = config_.base_statistical_weight * intensity_weight;
  
  // Update counters
  particles_generated_++;
  current_simulation_time_ = emission_time;
  
  return source_particle;
}

std::vector<SourceParticle> ParticleSource::generateParticles(size_t count)
{
  std::vector<SourceParticle> particles;
  particles.reserve(count);
  
  for (size_t i = 0; i < count; ++i) {
    particles.push_back(generateParticle());
  }
  
  return particles;
}

std::vector<SourceParticle> ParticleSource::generateParticlesForTime(double time_duration)
{
  std::vector<SourceParticle> particles;
  double start_time = current_simulation_time_;
  double end_time = start_time + time_duration;
  
  // Estimate number of particles for this time duration
  size_t estimated_count = static_cast<size_t>(config_.activity * time_duration * 1.1); // 10% buffer
  particles.reserve(estimated_count);
  
  while (current_simulation_time_ < end_time) {
    SourceParticle particle = generateParticle();
    if (particle.emission_time <= end_time) {
      particles.push_back(particle);
    } else {
      break;
    }
  }
  
  return particles;
}

std::vector<SourceParticle> ParticleSource::generateParticlesBatch(size_t count,
  std::function<void(size_t, size_t)> progress_callback)
{
  std::vector<SourceParticle> particles;
  particles.reserve(count);
  
  for (size_t i = 0; i < count; ++i) {
    particles.push_back(generateParticle());
    
    // Report progress every 1000 particles
    if (progress_callback && (i % 1000 == 0 || i == count - 1)) {
      progress_callback(i + 1, count);
    }
  }
  
  return particles;
}

Vector3D ParticleSource::sampleSourcePosition() const
{
  switch (config_.geometry.type) {
    case SourceGeometryType::Point:
      return samplePointSourcePosition();
    case SourceGeometryType::Plane:
      return samplePlaneSourcePosition();
    case SourceGeometryType::Sphere:
      return sampleSphereSourcePosition();
    case SourceGeometryType::Cylinder:
      return sampleCylinderSourcePosition();
    case SourceGeometryType::Rectangle:
      return sampleRectangleSourcePosition();
    case SourceGeometryType::Cone:
      return sampleConeSourcePosition();
    case SourceGeometryType::Custom:
      if (config_.geometry.custom_position_sampler) {
        return config_.geometry.custom_position_sampler();
      }
      break;
  }
  
  // Fallback to point source
  return config_.geometry.center;
}

Vector3D ParticleSource::sampleSourceDirection(const Vector3D& position) const
{
  return sampling_.sampleDirection(config_.angular_distribution);
}

double ParticleSource::sampleSourceEnergy() const
{
  return sampling_.sampleEnergy(config_.energy_spectrum);
}

double ParticleSource::sampleIntensityWeight(const Vector3D& position) const
{
  switch (config_.intensity.type) {
    case IntensityDistributionType::Uniform:
      return 1.0;
      
    case IntensityDistributionType::Gaussian: {
      Vector3D diff = position - config_.intensity.center;
      double distance_squared = diff.magnitudeSquared();
      double sigma = config_.intensity.parameter1;
      return std::exp(-distance_squared / (2.0 * sigma * sigma));
    }
    
    case IntensityDistributionType::Exponential: {
      Vector3D diff = position - config_.intensity.center;
      double distance = diff.magnitude();
      double lambda = config_.intensity.parameter1;
      return std::exp(-lambda * distance);
    }
    
    case IntensityDistributionType::Custom:
      if (config_.intensity.custom_intensity_function) {
        return config_.intensity.custom_intensity_function(position);
      }
      break;
  }
  
  return 1.0;
}

double ParticleSource::sampleEmissionTime() const
{
  // For steady-state sources, use uniform time sampling
  // For time-dependent sources, this could be enhanced
  return current_simulation_time_ + sampling_.getRNG().uniform(0.0, 1.0 / config_.activity);
}

Vector3D ParticleSource::samplePointSourcePosition() const
{
  return config_.geometry.center;
}

Vector3D ParticleSource::samplePlaneSourcePosition() const
{
  double width = config_.geometry.parameter1;
  double height = config_.geometry.parameter2;
  
  // Sample uniform position on rectangle
  double u = sampling_.getRNG().uniform(-width / 2.0, width / 2.0);
  double v = sampling_.getRNG().uniform(-height / 2.0, height / 2.0);
  
  // Create orthonormal basis from normal
  Vector3D normal = config_.geometry.direction;
  Vector3D u_axis, v_axis;
  
  // Find perpendicular vectors
  if (std::abs(normal.z()) < 0.9) {
    u_axis = Vector3D(0, 0, 1).cross(normal).normalise();
  } else {
    u_axis = Vector3D(1, 0, 0).cross(normal).normalise();
  }
  v_axis = normal.cross(u_axis);
  
  return config_.geometry.center + u_axis * u + v_axis * v;
}

Vector3D ParticleSource::sampleSphereSourcePosition() const
{
  double radius = config_.geometry.parameter1;
  bool surface_only = (config_.geometry.parameter2 > 0.5);
  
  if (surface_only) {
    // Sample on sphere surface
    Vector3D direction = sampling_.sampleIsotropicDirection();
    return config_.geometry.center + direction * radius;
  } else {
    // Sample in sphere volume
    double r = radius * std::pow(sampling_.getRNG().uniform(0.0, 1.0), 1.0/3.0);
    Vector3D direction = sampling_.sampleIsotropicDirection();
    return config_.geometry.center + direction * r;
  }
}

Vector3D ParticleSource::sampleCylinderSourcePosition() const
{
  double radius = config_.geometry.parameter1;
  double height = config_.geometry.parameter2;
  bool surface_only = (config_.geometry.parameter3 > 0.5);
  
  Vector3D axis = config_.geometry.direction;
  
  if (surface_only) {
    // Sample on cylinder surface (sides + ends)
    double total_area = 2.0 * M_PI * radius * height + 2.0 * M_PI * radius * radius;
    double side_area = 2.0 * M_PI * radius * height;
    
    if (sampling_.getRNG().uniform(0.0, total_area) < side_area) {
      // Sample on cylindrical side
      double phi = sampling_.getRNG().uniform(0.0, 2.0 * M_PI);
      double z = sampling_.getRNG().uniform(-height/2.0, height/2.0);
      
      Vector3D u_axis, v_axis;
      if (std::abs(axis.z()) < 0.9) {
        u_axis = Vector3D(0, 0, 1).cross(axis).normalise();
      } else {
        u_axis = Vector3D(1, 0, 0).cross(axis).normalise();
      }
      v_axis = axis.cross(u_axis);
      
      Vector3D radial = u_axis * std::cos(phi) + v_axis * std::sin(phi);
      return config_.geometry.center + radial * radius + axis * z;
    } else {
      // Sample on end caps
      double r = radius * std::sqrt(sampling_.getRNG().uniform(0.0, 1.0));
      double phi = sampling_.getRNG().uniform(0.0, 2.0 * M_PI);
      double z = sampling_.getRNG().uniform(0.0, 1.0) > 0.5 ? height/2.0 : -height/2.0;
      
      Vector3D u_axis, v_axis;
      if (std::abs(axis.z()) < 0.9) {
        u_axis = Vector3D(0, 0, 1).cross(axis).normalise();
      } else {
        u_axis = Vector3D(1, 0, 0).cross(axis).normalise();
      }
      v_axis = axis.cross(u_axis);
      
      Vector3D radial = u_axis * std::cos(phi) + v_axis * std::sin(phi);
      return config_.geometry.center + radial * r + axis * z;
    }
  } else {
    // Sample in cylinder volume
    double r = radius * std::sqrt(sampling_.getRNG().uniform(0.0, 1.0));
    double phi = sampling_.getRNG().uniform(0.0, 2.0 * M_PI);
    double z = sampling_.getRNG().uniform(-height/2.0, height/2.0);
    
    Vector3D u_axis, v_axis;
    if (std::abs(axis.z()) < 0.9) {
      u_axis = Vector3D(0, 0, 1).cross(axis).normalise();
    } else {
      u_axis = Vector3D(1, 0, 0).cross(axis).normalise();
    }
    v_axis = axis.cross(u_axis);
    
    Vector3D radial = u_axis * std::cos(phi) + v_axis * std::sin(phi);
    return config_.geometry.center + radial * r + axis * z;
  }
}

Vector3D ParticleSource::sampleRectangleSourcePosition() const
{
  double length1 = config_.geometry.parameter1;
  double length2 = config_.geometry.parameter2;
  
  double u = sampling_.getRNG().uniform(0.0, 1.0);
  double v = sampling_.getRNG().uniform(0.0, 1.0);
  
  // For simplicity, assume rectangle is aligned with coordinate axes
  // In a full implementation, you'd construct the proper orthogonal basis
  Vector3D side1 = config_.geometry.direction * length1;
  Vector3D side2 = Vector3D(0, 0, 1).cross(side1).normalise() * length2;
  
  return config_.geometry.center + side1 * u + side2 * v;
}

Vector3D ParticleSource::sampleConeSourcePosition() const
{
  // Sample on cone surface (simplified implementation)
  double base_radius = config_.geometry.parameter1;
  double height = config_.geometry.parameter2;
  
  double z = sampling_.getRNG().uniform(0.0, height);
  double r = base_radius * (1.0 - z / height); // Linear taper
  double phi = sampling_.getRNG().uniform(0.0, 2.0 * M_PI);
  
  Vector3D axis = config_.geometry.direction;
  Vector3D u_axis, v_axis;
  
  if (std::abs(axis.z()) < 0.9) {
    u_axis = Vector3D(0, 0, 1).cross(axis).normalise();
  } else {
    u_axis = Vector3D(1, 0, 0).cross(axis).normalise();
  }
  v_axis = axis.cross(u_axis);
  
  Vector3D radial = u_axis * std::cos(phi) + v_axis * std::sin(phi);
  return config_.geometry.center + radial * r + axis * z;
}

void ParticleSource::validateAndThrow() const
{
  if (!config_.isValid()) {
    throw std::invalid_argument("Invalid source configuration");
  }
}

void ParticleSource::reset()
{
  particles_generated_ = 0;
  current_simulation_time_ = 0.0;
}

std::string ParticleSource::getSourceSummary() const
{
  std::stringstream ss;
  ss << "Source Summary:\n";
  ss << config_.toString();
  ss << "  Particles Generated: " << particles_generated_ << "\n";
  ss << "  Current Time: " << current_simulation_time_ << " s\n";
  return ss.str();
}

// =============================================================================
// Factory Methods for Common Source Configurations
// =============================================================================

SourceConfiguration ParticleSource::createMonoenergeticPointSource(const Vector3D& position,
                                                                  double energy_MeV,
                                                                  ParticleType type)
{
  SourceConfiguration config;
  config.name = "Monoenergetic Point Source";
  config.particle_type = type;
  config.energy_spectrum = EnergySpectrum::monoenergetic(energy_MeV);
  config.angular_distribution = AngularDistribution::isotropic();
  config.geometry = SourceGeometry::createPoint(position);
  config.intensity = IntensityDistribution::createUniform();
  return config;
}

SourceConfiguration ParticleSource::createIsotropicPointSource(const Vector3D& position,
                                                             const EnergySpectrum& spectrum,
                                                             ParticleType type)
{
  SourceConfiguration config;
  config.name = "Isotropic Point Source";
  config.particle_type = type;
  config.energy_spectrum = spectrum;
  config.angular_distribution = AngularDistribution::isotropic();
  config.geometry = SourceGeometry::createPoint(position);
  config.intensity = IntensityDistribution::createUniform();
  return config;
}

SourceConfiguration ParticleSource::createBeamSource(const Vector3D& position,
                                                   const Vector3D& direction,
                                                   double energy_MeV,
                                                   double beam_divergence,
                                                   ParticleType type)
{
  SourceConfiguration config;
  config.name = "Beam Source";
  config.particle_type = type;
  config.energy_spectrum = EnergySpectrum::monoenergetic(energy_MeV);
  
  if (beam_divergence > 0.0) {
    config.angular_distribution = AngularDistribution::cone(direction, beam_divergence);
  } else {
    config.angular_distribution = AngularDistribution::beam(direction);
  }
  
  config.geometry = SourceGeometry::createPoint(position);
  config.intensity = IntensityDistribution::createUniform();
  return config;
}

// =============================================================================
// SourceAnalysis Helper Functions
// =============================================================================

namespace SourceAnalysis
{

double calculateTotalSourceStrength(const SourceConfiguration& config)
{
  return config.activity; // Particles per second
}

double estimateMemoryRequirement(const SourceConfiguration& config)
{
  // Rough estimate: each SourceParticle is about 200 bytes
  size_t particles = static_cast<size_t>(config.total_particles);
  return particles * 200.0 / (1024.0 * 1024.0); // MB
}

bool validateSourceGeometry(const SourceGeometry& geometry)
{
  switch (geometry.type) {
    case SourceGeometryType::Sphere:
      return geometry.parameter1 > 0.0; // Radius must be positive
    case SourceGeometryType::Cylinder:
      return geometry.parameter1 > 0.0 && geometry.parameter2 > 0.0; // Radius and height
    case SourceGeometryType::Plane:
      return geometry.parameter1 > 0.0 && geometry.parameter2 > 0.0; // Width and height
    default:
      return true;
  }
}

double calculateEffectiveSourceArea(const SourceGeometry& geometry)
{
  switch (geometry.type) {
    case SourceGeometryType::Point:
      return 0.0;
    case SourceGeometryType::Plane:
    case SourceGeometryType::Rectangle:
      return geometry.parameter1 * geometry.parameter2;
    case SourceGeometryType::Sphere:
      return 4.0 * M_PI * geometry.parameter1 * geometry.parameter1;
    case SourceGeometryType::Cylinder: {
      double radius = geometry.parameter1;
      double height = geometry.parameter2;
      return 2.0 * M_PI * radius * height + 2.0 * M_PI * radius * radius;
    }
    default:
      return 1.0; // Default area
  }
}

std::string generateQualityAssuranceReport(const SourceConfiguration& config)
{
  std::stringstream ss;
  ss << "=== SOURCE QUALITY ASSURANCE REPORT ===\n\n";
  ss << config.toString();
  ss << "\nValidation Results:\n";
  ss << "  Configuration Valid: " << (config.isValid() ? "PASS" : "FAIL") << "\n";
  ss << "  Geometry Valid: " << (validateSourceGeometry(config.geometry) ? "PASS" : "FAIL") << "\n";
  ss << "  Estimated Memory: " << estimateMemoryRequirement(config) << " MB\n";
  ss << "  Effective Area: " << calculateEffectiveSourceArea(config.geometry) << " cm²\n";
  return ss.str();
}

} // namespace SourceAnalysis