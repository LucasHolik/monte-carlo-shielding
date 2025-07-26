#pragma once

#include "Geometry.hpp"
#include "Material.hpp"
#include "MonteCarloSampling.hpp"
#include "Particle.hpp"
#include "Vector3D.hpp"

#include <functional>
#include <optional>
#include <string>
#include <vector>

/**
 * @brief Result of a particle transport step
 */
struct TransportStep
{
  double distance_travelled = 0.0;
  Vector3D final_position;
  Material material_encountered;
  bool hit_boundary = false;
  bool interaction_occurred = false;
  std::optional<GeometryIntersection> boundary_info;

  TransportStep() = default;
  TransportStep(double dist, const Vector3D &pos, const Material &mat)
      : distance_travelled(dist), final_position(pos), material_encountered(mat)
  {}
};

/**
 * @brief Cross-section function type for calculating interaction probabilities
 */
using CrossSectionFunction = std::function<double(const Material &, double)>;

/**
 * @brief Monte Carlo particle transport engine
 *
 * Modern C++17 implementation providing comprehensive particle transport
 * through complex geometries with realistic physics interactions.
 */
class Transport
{
private:
  const Geometry &geometry_;
  MonteCarloSampling &sampling_;
  CrossSectionFunction cross_section_function_;
  double maximum_step_size_;
  bool enable_boundary_crossing_;
  bool enable_interaction_sampling_;

public:
  // Constructor
  Transport(const Geometry &geometry, MonteCarloSampling &sampling);

  // Configuration
  void setCrossSectionFunction(const CrossSectionFunction &func)
  {
    cross_section_function_ = func;
  }
  void setMaximumStepSize(double max_step) { maximum_step_size_ = max_step; }
  void enableBoundaryCrossing(bool enable)
  {
    enable_boundary_crossing_ = enable;
  }
  void enableInteractionSampling(bool enable)
  {
    enable_interaction_sampling_ = enable;
  }

  // Main transport methods
  void trackParticle(Particle &particle);
  TransportStep transportParticleStep(Particle &particle);

  // Core transport components
  double sampleDistanceToInteraction(const Particle &particle,
                                     const Material &material);
  void moveParticle(Particle &particle, double distance);
  bool hitBoundary(const Particle &particle, double proposed_distance);
  void handleBoundary(Particle &particle);
  void performInteraction(Particle &particle, const Material &material);

  // Advanced transport features
  std::vector<TransportStep> getParticleTrack(Particle particle,
                                              double max_distance = 1e30);
  bool willParticleEscape(const Particle &particle);
  std::optional<Vector3D> findEscapePoint(const Particle &particle);

  // Utility methods
  Material getCurrentMaterial(const Particle &particle) const;
  std::optional<GeometryIntersection>
  findNextBoundary(const Particle &particle) const;
  double getDistanceToBoundary(const Particle &particle) const;

  // Validation and diagnostics
  bool isValidTransportState(const Particle &particle) const;
  std::string getTransportDiagnostics(const Particle &particle) const;

private:
  // Internal helper methods
  double calculateLinearAttenuationCoefficient(const Material &material,
                                               double energy) const;
  void validateParticleState(const Particle &particle) const;
  void handleGeometryTransition(Particle &particle,
                                const GeometryIntersection &intersection);
};

// Default cross-section functions for common particle types
namespace DefaultCrossSections
{
/**
 * @brief Simple photon total cross-section (placeholder)
 */
double photonTotalCrossSection(const Material &material, double energy);

/**
 * @brief Simple neutron total cross-section (placeholder)
 */
double neutronTotalCrossSection(const Material &material, double energy);

/**
 * @brief Simple electron cross-section (placeholder)
 */
double electronCrossSection(const Material &material, double energy);
} // namespace DefaultCrossSections