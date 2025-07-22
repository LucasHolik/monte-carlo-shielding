#include "Particle.hpp"

#include <stdexcept>
#include <sstream>
#include <atomic>
#include <cmath>
#include <algorithm>

// Physical constants (in appropriate units)
namespace PhysicalConstants
{
  constexpr double SPEED_OF_LIGHT = 299792458.0;  // m/s
  constexpr double ELECTRON_REST_MASS = 0.510999; // MeV/c²
  constexpr double PROTON_REST_MASS = 938.272;    // MeV/c²
  constexpr double NEUTRON_REST_MASS = 939.565;   // MeV/c²
}

// Static history ID counter
static std::atomic<int> global_history_id{1};

// Constructors
Particle::Particle()
    : position_(Vector3D::ZERO),
      direction_(Vector3D::UNIT_Z),
      energy_(0.0),
      weight_(1.0),
      type_(ParticleType::Photon),
      state_(ParticleState::Alive),
      generation_(0),
      history_id_(getNextHistoryId())
{
  validateConstruction();
}

Particle::Particle(ParticleType type, const Vector3D &position, const Vector3D &direction,
                   double energy, double weight)
    : position_(position),
      direction_(direction),
      energy_(energy),
      weight_(weight),
      type_(type),
      state_(ParticleState::Alive),
      birth_energy_(energy),
      generation_(0),
      history_id_(getNextHistoryId())
{
  validateConstruction();
}

void Particle::validateConstruction()
{
  if (energy_ < 0.0)
  {
    throw std::invalid_argument("Particle energy cannot be negative");
  }
  if (weight_ <= 0.0)
  {
    throw std::invalid_argument("Particle weight must be positive");
  }
  if (direction_.isZero())
  {
    throw std::invalid_argument("Particle direction cannot be zero vector");
  }

  // Ensure direction is normalised
  auto normalised_direction = direction_.tryNormalise();
  if (!normalised_direction)
  {
    throw std::invalid_argument("Cannot normalise particle direction vector");
  }
  direction_ = *normalised_direction;
}

// Safe mutators with std::optional
std::optional<double> Particle::trySetEnergy(double new_energy)
{
  if (new_energy < 0.0)
  {
    return std::nullopt;
  }

  double old_energy = energy_;
  energy_ = new_energy;
  return old_energy;
}

std::optional<double> Particle::trySetWeight(double new_weight)
{
  if (new_weight <= 0.0)
  {
    return std::nullopt;
  }

  double old_weight = weight_;
  weight_ = new_weight;
  return old_weight;
}

std::optional<Vector3D> Particle::trySetDirection(const Vector3D &new_direction)
{
  auto normalised = new_direction.tryNormalise();
  if (!normalised)
  {
    return std::nullopt;
  }

  Vector3D old_direction = direction_;
  direction_ = *normalised;
  return old_direction;
}

// Interaction tracking
void Particle::recordInteraction(const Vector3D &interaction_point)
{
  last_interaction_point_ = interaction_point;
}

// Particle movement
void Particle::move(double distance)
{
  if (distance < 0.0)
  {
    throw std::invalid_argument("Cannot move particle negative distance");
  }

  if (isAlive())
  {
    position_ += direction_ * distance;
  }
}

void Particle::moveToPosition(const Vector3D &new_position)
{
  if (isAlive())
  {
    position_ = new_position;
  }
}

// Energy operations
std::optional<double> Particle::loseEnergy(double energy_loss)
{
  if (energy_loss < 0.0 || energy_loss > energy_)
  {
    return std::nullopt;
  }

  double old_energy = energy_;
  energy_ -= energy_loss;

  // Check if particle should be terminated due to low energy
  if (energy_ < 1e-6)
  { // 1 eV threshold
    kill();
  }

  return old_energy;
}

double Particle::getKineticEnergy() const
{
  if (type_ == ParticleType::Photon)
  {
    return energy_; // Photons: E = pc, no rest mass
  }
  else
  {
    // For massive particles: KE = Total Energy - Rest Mass Energy
    return energy_ - getRestMassEnergy();
  }
}

double Particle::getRestMassEnergy() const
{
  return getRestMass(type_);
}

double Particle::getTotalEnergy() const
{
  if (type_ == ParticleType::Photon)
  {
    return energy_;
  }
  else
  {
    return energy_; // We store total energy for massive particles
  }
}

// Relativistic calculations
double Particle::getSpeed() const
{
  if (type_ == ParticleType::Photon)
  {
    return 1.0; // Speed of light in natural units (c = 1)
  }

  double gamma = getGamma();
  return std::sqrt(1.0 - 1.0 / (gamma * gamma));
}

double Particle::getBeta() const
{
  return getSpeed(); // β = v/c
}

double Particle::getGamma() const
{
  if (type_ == ParticleType::Photon)
  {
    return std::numeric_limits<double>::infinity();
  }

  double rest_mass = getRestMassEnergy();
  return getTotalEnergy() / rest_mass;
}

double Particle::getMomentum() const
{
  if (type_ == ParticleType::Photon)
  {
    return energy_; // p = E/c for photons
  }
  else
  {
    double gamma = getGamma();
    double rest_mass = getRestMassEnergy();
    return gamma * getBeta() * rest_mass;
  }
}

// Utility methods
bool Particle::isCharged() const
{
  return getCharge(type_) != 0;
}

// Secondary particle creation
Particle Particle::createSecondary(ParticleType type, const Vector3D &position,
                                   const Vector3D &direction, double energy) const
{
  Particle secondary(type, position, direction, energy);
  secondary.setGeneration(generation_ + 1);
  secondary.setBirthEnergy(energy);
  return secondary;
}

// String representations
std::string Particle::toString() const
{
  std::stringstream ss;
  ss << getTypeName() << " [" << getStateName() << "] "
     << "at (" << position_.x() << ", " << position_.y() << ", " << position_.z() << ") "
     << "with energy " << energy_ << " MeV, "
     << "weight " << weight_ << ", "
     << "generation " << generation_;
  return ss.str();
}

std::string Particle::getTypeName() const
{
  return particleTypeToString(type_);
}

std::string Particle::getStateName() const
{
  return particleStateToString(state_);
}

// Comparison operators
bool Particle::operator==(const Particle &other) const
{
  constexpr double tolerance = 1e-10;
  return type_ == other.type_ &&
         state_ == other.state_ &&
         position_ == other.position_ &&
         direction_ == other.direction_ &&
         std::abs(energy_ - other.energy_) < tolerance &&
         std::abs(weight_ - other.weight_) < tolerance &&
         generation_ == other.generation_;
}

bool Particle::operator!=(const Particle &other) const
{
  return !(*this == other);
}

// Static factory methods
Particle Particle::createPhoton(const Vector3D &position, const Vector3D &direction, double energy)
{
  return Particle(ParticleType::Photon, position, direction, energy);
}

Particle Particle::createNeutron(const Vector3D &position, const Vector3D &direction, double energy)
{
  return Particle(ParticleType::Neutron, position, direction, energy);
}

Particle Particle::createElectron(const Vector3D &position, const Vector3D &direction, double energy)
{
  return Particle(ParticleType::Electron, position, direction, energy);
}

Particle Particle::createProton(const Vector3D &position, const Vector3D &direction, double energy)
{
  return Particle(ParticleType::Proton, position, direction, energy);
}

// Physical constants
double Particle::getRestMass(ParticleType type)
{
  switch (type)
  {
  case ParticleType::Photon:
    return 0.0; // Massless
  case ParticleType::Electron:
    return PhysicalConstants::ELECTRON_REST_MASS;
  case ParticleType::Proton:
    return PhysicalConstants::PROTON_REST_MASS;
  case ParticleType::Neutron:
    return PhysicalConstants::NEUTRON_REST_MASS;
  default:
    throw std::invalid_argument("Unknown particle type in getRestMass");
  }
}

int Particle::getCharge(ParticleType type)
{
  switch (type)
  {
  case ParticleType::Photon:
    return 0;
  case ParticleType::Electron:
    return -1;
  case ParticleType::Proton:
    return +1;
  case ParticleType::Neutron:
    return 0;
  default:
    throw std::invalid_argument("Unknown particle type in getCharge");
  }
}

std::string Particle::getParticleSymbol(ParticleType type)
{
  switch (type)
  {
  case ParticleType::Photon:
    return "γ";
  case ParticleType::Electron:
    return "e⁻";
  case ParticleType::Proton:
    return "p";
  case ParticleType::Neutron:
    return "n";
  default:
    return "?";
  }
}

// Private methods
int Particle::getNextHistoryId()
{
  return global_history_id.fetch_add(1);
}

// Utility functions
std::string particleTypeToString(ParticleType type)
{
  switch (type)
  {
  case ParticleType::Photon:
    return "Photon";
  case ParticleType::Neutron:
    return "Neutron";
  case ParticleType::Electron:
    return "Electron";
  case ParticleType::Proton:
    return "Proton";
  default:
    return "Unknown";
  }
}

std::string particleStateToString(ParticleState state)
{
  switch (state)
  {
  case ParticleState::Alive:
    return "Alive";
  case ParticleState::Absorbed:
    return "Absorbed";
  case ParticleState::Escaped:
    return "Escaped";
  case ParticleState::Scattered:
    return "Scattered";
  case ParticleState::Split:
    return "Split";
  case ParticleState::Terminated:
    return "Terminated";
  default:
    return "Unknown";
  }
}

std::optional<ParticleType> stringToParticleType(const std::string &type_str)
{
  std::string lower_str = type_str;
  std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);

  if (lower_str == "photon" || lower_str == "gamma")
  {
    return ParticleType::Photon;
  }
  else if (lower_str == "neutron")
  {
    return ParticleType::Neutron;
  }
  else if (lower_str == "electron")
  {
    return ParticleType::Electron;
  }
  else if (lower_str == "proton")
  {
    return ParticleType::Proton;
  }
  else
  {
    return std::nullopt;
  }
}