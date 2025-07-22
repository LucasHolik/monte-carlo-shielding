#pragma once

#include "Vector3D.hpp"

#include <optional>
#include <variant>
#include <string>

/**
 * @brief Particle types supported by the Monte Carlo simulator
 */
enum class ParticleType
{
  Photon,
  Neutron,
  Electron,
  Proton
};

/**
 * @brief Particle state during transport simulation
 */
enum class ParticleState
{
  Alive,
  Absorbed,
  Escaped,
  Scattered,
  Split,
  Terminated
};

/**
 * @brief Particle class for Monte Carlo radiation transport simulation
 *
 * Modern C++17 implementation supporting multi-particle transport with
 * type-safe operations and comprehensive state management.
 */
class Particle
{
private:
  Vector3D position_;
  Vector3D direction_;
  double energy_; // Energy in MeV
  double weight_; // Statistical weight for variance reduction
  ParticleType type_;
  ParticleState state_;
  std::optional<Vector3D> last_interaction_point_;
  std::optional<double> birth_energy_; // Initial energy for analysis
  int generation_;                     // Track particle generation (0 = primary)
  int history_id_;                     // Unique identifier for this particle history

public:
  // Constructors
  Particle();
  Particle(ParticleType type, const Vector3D &position, const Vector3D &direction,
           double energy, double weight = 1.0);
  Particle(const Particle &other) = default;
  Particle &operator=(const Particle &other) = default;
  Particle(Particle &&other) noexcept = default;
  Particle &operator=(Particle &&other) noexcept = default;

  // Accessors
  Vector3D position() const { return position_; }
  Vector3D direction() const { return direction_; }
  double energy() const { return energy_; }
  double weight() const { return weight_; }
  ParticleType type() const { return type_; }
  ParticleState state() const { return state_; }
  int generation() const { return generation_; }
  int historyId() const { return history_id_; }

  // C++17 structured binding support
  auto getPositionAndDirection() const
  {
    return std::make_pair(position_, direction_);
  }

  auto getEnergyAndWeight() const
  {
    return std::make_tuple(energy_, weight_);
  }

  auto getTypeAndState() const
  {
    return std::make_pair(type_, state_);
  }

  // Safe mutators with C++17 std::optional
  std::optional<double> trySetEnergy(double new_energy);
  std::optional<double> trySetWeight(double new_weight);
  std::optional<Vector3D> trySetDirection(const Vector3D &new_direction);

  // Direct mutators (for performance-critical sections)
  void setPosition(const Vector3D &position) { position_ = position; }
  void setDirection(const Vector3D &direction) { direction_ = direction; }
  void setEnergy(double energy) { energy_ = energy; }
  void setWeight(double weight) { weight_ = weight; }
  void setState(ParticleState state) { state_ = state; }
  void setGeneration(int generation) { generation_ = generation; }
  void setHistoryId(int history_id) { history_id_ = history_id; }

  // Interaction tracking
  void recordInteraction(const Vector3D &interaction_point);
  std::optional<Vector3D> getLastInteractionPoint() const { return last_interaction_point_; }

  std::optional<double> getBirthEnergy() const { return birth_energy_; }
  void setBirthEnergy(double energy) { birth_energy_ = energy; }

  // Particle movement
  void move(double distance);
  void moveToPosition(const Vector3D &new_position);

  // State management
  bool isAlive() const { return state_ == ParticleState::Alive; }
  void kill() { state_ = ParticleState::Terminated; }
  void absorb() { state_ = ParticleState::Absorbed; }
  void escape() { state_ = ParticleState::Escaped; }
  void scatter() { state_ = ParticleState::Scattered; }

  // Energy operations
  std::optional<double> loseEnergy(double energy_loss);
  double getKineticEnergy() const;
  double getRestMassEnergy() const;
  double getTotalEnergy() const;

  // Particle physics utilities
  double getSpeed() const;    // Fraction of speed of light
  double getBeta() const;     // v/c
  double getGamma() const;    // Lorentz factor
  double getMomentum() const; // MeV/c

  // Utility methods
  bool isPhoton() const { return type_ == ParticleType::Photon; }
  bool isNeutron() const { return type_ == ParticleType::Neutron; }
  bool isElectron() const { return type_ == ParticleType::Electron; }
  bool isProton() const { return type_ == ParticleType::Proton; }
  bool isCharged() const;

  // Particle creation for secondaries
  Particle createSecondary(ParticleType type, const Vector3D &position,
                           const Vector3D &direction, double energy) const;

  // String representation
  std::string toString() const;
  std::string getTypeName() const;
  std::string getStateName() const;

  // Comparison operators
  bool operator==(const Particle &other) const;
  bool operator!=(const Particle &other) const;

  // Static factory methods
  static Particle createPhoton(const Vector3D &position, const Vector3D &direction, double energy);
  static Particle createNeutron(const Vector3D &position, const Vector3D &direction, double energy);
  static Particle createElectron(const Vector3D &position, const Vector3D &direction, double energy);
  static Particle createProton(const Vector3D &position, const Vector3D &direction, double energy);

  // Physical constants access
  static double getRestMass(ParticleType type); // MeV/c²
  static int getCharge(ParticleType type);      // Elementary charges
  static std::string getParticleSymbol(ParticleType type);

private:
  void validateConstruction();
  static int getNextHistoryId();
};

// Utility functions for particle type conversions
std::string particleTypeToString(ParticleType type);
std::string particleStateToString(ParticleState state);
std::optional<ParticleType> stringToParticleType(const std::string &type_str);