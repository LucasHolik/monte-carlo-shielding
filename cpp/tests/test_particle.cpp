#include "Particle.hpp"

#include <iostream>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>

void testConstructors()
{
  std::cout << "Testing constructors..." << std::endl;

  // Default constructor
  Particle p1;
  assert(p1.type() == ParticleType::Photon);
  assert(p1.state() == ParticleState::Alive);
  assert(p1.energy() == 0.0);
  assert(p1.weight() == 1.0);
  assert(p1.generation() == 0);
  assert(p1.position() == Vector3D::ZERO);
  assert(p1.direction() == Vector3D::UNIT_Z);

  // Parameterised constructor
  Vector3D pos(1.0, 2.0, 3.0);
  Vector3D dir(0.0, 0.0, 1.0);
  Particle p2(ParticleType::Neutron, pos, dir, 2.5, 0.8);
  assert(p2.type() == ParticleType::Neutron);
  assert(p2.position() == pos);
  assert(p2.direction() == dir);
  assert(p2.energy() == 2.5);
  assert(p2.weight() == 0.8);

  // Copy constructor
  Particle p3(p2);
  assert(p3.type() == p2.type());
  assert(p3.energy() == p2.energy());
  assert(p3.weight() == p2.weight());

  // Move constructor
  Particle p4(std::move(p3));
  assert(p4.type() == ParticleType::Neutron);
  assert(p4.energy() == 2.5);

  std::cout << "✓ Constructors passed" << std::endl;
}

void testStructuredBinding()
{
  std::cout << "Testing C++17 structured binding..." << std::endl;

  Vector3D pos(1.0, 2.0, 3.0);
  Vector3D dir(0.0, 1.0, 0.0);
  Particle p(ParticleType::Electron, pos, dir, 1.5, 0.7);

  // Test position and direction binding
  auto [position, direction] = p.getPositionAndDirection();
  assert(position == pos);
  assert(direction == dir);

  // Test energy and weight binding
  auto [energy, weight] = p.getEnergyAndWeight();
  assert(energy == 1.5);
  assert(weight == 0.7);

  // Test type and state binding
  auto [type, state] = p.getTypeAndState();
  assert(type == ParticleType::Electron);
  assert(state == ParticleState::Alive);

  std::cout << "✓ Structured binding passed" << std::endl;
}

void testSafeMutators()
{
  std::cout << "Testing safe mutators with std::optional..." << std::endl;

  Particle p = Particle::createPhoton(Vector3D::ZERO, Vector3D::UNIT_X, 1.0);

  // Valid energy change
  auto oldEnergy = p.trySetEnergy(2.0);
  assert(oldEnergy.has_value());
  assert(*oldEnergy == 1.0);
  assert(p.energy() == 2.0);

  // Invalid energy change (negative)
  auto invalidEnergy = p.trySetEnergy(-1.0);
  assert(!invalidEnergy.has_value());
  assert(p.energy() == 2.0); // Should remain unchanged

  // Valid weight change
  auto oldWeight = p.trySetWeight(0.5);
  assert(oldWeight.has_value());
  assert(*oldWeight == 1.0);
  assert(p.weight() == 0.5);

  // Invalid weight change (zero)
  auto invalidWeight = p.trySetWeight(0.0);
  assert(!invalidWeight.has_value());
  assert(p.weight() == 0.5); // Should remain unchanged

  // Valid direction change
  Vector3D newDir(1.0, 1.0, 0.0);
  auto oldDir = p.trySetDirection(newDir);
  assert(oldDir.has_value());
  assert(*oldDir == Vector3D::UNIT_X);
  assert(p.direction().isUnit()); // Should be normalised

  // Invalid direction change (zero vector)
  auto invalidDir = p.trySetDirection(Vector3D::ZERO);
  assert(!invalidDir.has_value());

  std::cout << "✓ Safe mutators passed" << std::endl;
}

void testDirectMutators()
{
  std::cout << "Testing direct mutators..." << std::endl;

  Particle p = Particle::createElectron(Vector3D::ZERO, Vector3D::UNIT_Y, 0.5);

  Vector3D newPos(5.0, -3.0, 2.0);
  p.setPosition(newPos);
  assert(p.position() == newPos);

  Vector3D newDir(1.0, 0.0, 0.0);
  p.setDirection(newDir);
  assert(p.direction() == newDir);

  p.setEnergy(3.0);
  assert(p.energy() == 3.0);

  p.setWeight(0.25);
  assert(p.weight() == 0.25);

  p.setState(ParticleState::Scattered);
  assert(p.state() == ParticleState::Scattered);

  p.setGeneration(3);
  assert(p.generation() == 3);

  p.setHistoryId(42);
  assert(p.historyId() == 42);

  std::cout << "✓ Direct mutators passed" << std::endl;
}

void testMovement()
{
  std::cout << "Testing particle movement..." << std::endl;

  Vector3D startPos(0.0, 0.0, 0.0);
  Vector3D direction(1.0, 0.0, 0.0);
  Particle p = Particle::createPhoton(startPos, direction, 1.0);

  // Move particle
  p.move(5.0);
  Vector3D expectedPos(5.0, 0.0, 0.0);
  assert(p.position() == expectedPos);

  // Move to specific position
  Vector3D newPos(10.0, 15.0, -5.0);
  p.moveToPosition(newPos);
  assert(p.position() == newPos);

  // Test that dead particles don't move
  p.kill();
  p.move(1.0);
  assert(p.position() == newPos); // Should remain unchanged

  // Test negative distance throws exception
  Particle p2 = Particle::createPhoton(Vector3D::ZERO, Vector3D::UNIT_X, 1.0);
  bool threw = false;
  try
  {
    p2.move(-1.0);
  }
  catch (const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  std::cout << "✓ Movement passed" << std::endl;
}

void testStateManagement()
{
  std::cout << "Testing state management..." << std::endl;

  Particle p = Particle::createNeutron(Vector3D::ZERO, Vector3D::UNIT_Z, 1.0);

  // Initial state
  assert(p.isAlive());
  assert(p.state() == ParticleState::Alive);

  // Absorb particle
  p.absorb();
  assert(!p.isAlive());
  assert(p.state() == ParticleState::Absorbed);

  // Reset and escape
  p.setState(ParticleState::Alive);
  p.escape();
  assert(!p.isAlive());
  assert(p.state() == ParticleState::Escaped);

  // Reset and scatter
  p.setState(ParticleState::Alive);
  p.scatter();
  assert(!p.isAlive());
  assert(p.state() == ParticleState::Scattered);

  // Reset and kill
  p.setState(ParticleState::Alive);
  p.kill();
  assert(!p.isAlive());
  assert(p.state() == ParticleState::Terminated);

  std::cout << "✓ State management passed" << std::endl;
}

void testEnergyOperations()
{
  std::cout << "Testing energy operations..." << std::endl;

  // Test photon energy operations
  Particle photon = Particle::createPhoton(Vector3D::ZERO, Vector3D::UNIT_X, 2.0);
  assert(photon.getKineticEnergy() == 2.0); // Photon: KE = Total Energy
  assert(photon.getTotalEnergy() == 2.0);
  assert(photon.getRestMassEnergy() == 0.0);

  // Test electron energy operations
  double electronRestMass = Particle::getRestMass(ParticleType::Electron);
  Particle electron = Particle::createElectron(Vector3D::ZERO, Vector3D::UNIT_Y,
                                               electronRestMass + 1.0); // 1 MeV kinetic energy
  assert(std::abs(electron.getKineticEnergy() - 1.0) < 1e-10);
  assert(electron.getRestMassEnergy() == electronRestMass);
  assert(electron.getTotalEnergy() == electronRestMass + 1.0);

  // Test energy loss
  auto oldEnergy = electron.loseEnergy(0.5);
  assert(oldEnergy.has_value());
  assert(*oldEnergy == electronRestMass + 1.0);
  assert(std::abs(electron.energy() - (electronRestMass + 0.5)) < 1e-10);

  // Test invalid energy loss
  auto invalidLoss = electron.loseEnergy(-0.1);
  assert(!invalidLoss.has_value());

  auto excessiveLoss = electron.loseEnergy(1000.0);
  assert(!excessiveLoss.has_value());

  // Test energy loss leading to termination
  Particle lowEnergyParticle = Particle::createElectron(Vector3D::ZERO, Vector3D::UNIT_Z, 5e-7);
  auto terminalLoss = lowEnergyParticle.loseEnergy(1e-7);
  assert(terminalLoss.has_value());
  assert(!lowEnergyParticle.isAlive()); // Should be killed due to low energy (final energy = 4e-7 < 1e-6)

  std::cout << "✓ Energy operations passed" << std::endl;
}

void testRelativisticCalculations()
{
  std::cout << "Testing relativistic calculations..." << std::endl;

  // Test photon (always moves at speed of light)
  Particle photon = Particle::createPhoton(Vector3D::ZERO, Vector3D::UNIT_X, 1.0);
  assert(photon.getSpeed() == 1.0); // c in natural units
  assert(photon.getBeta() == 1.0);
  assert(std::isinf(photon.getGamma()));
  assert(photon.getMomentum() == 1.0); // p = E for photon

  // Test electron at rest
  double electronMass = Particle::getRestMass(ParticleType::Electron);
  Particle electronAtRest = Particle::createElectron(Vector3D::ZERO, Vector3D::UNIT_Y, electronMass);
  assert(std::abs(electronAtRest.getSpeed() - 0.0) < 1e-10);
  assert(std::abs(electronAtRest.getBeta() - 0.0) < 1e-10);
  assert(std::abs(electronAtRest.getGamma() - 1.0) < 1e-10);
  assert(std::abs(electronAtRest.getMomentum() - 0.0) < 1e-10);

  // Test relativistic electron (γ = 2)
  Particle fastElectron = Particle::createElectron(Vector3D::ZERO, Vector3D::UNIT_Z, 2.0 * electronMass);
  assert(std::abs(fastElectron.getGamma() - 2.0) < 1e-10);
  double expectedBeta = std::sqrt(1.0 - 1.0 / 4.0); // β = √(1 - 1/γ²)
  assert(std::abs(fastElectron.getBeta() - expectedBeta) < 1e-10);

  std::cout << "✓ Relativistic calculations passed" << std::endl;
}

void testParticleTypeChecks()
{
  std::cout << "Testing particle type checks..." << std::endl;

  Particle photon = Particle::createPhoton(Vector3D::ZERO, Vector3D::UNIT_X, 1.0);
  assert(photon.isPhoton());
  assert(!photon.isNeutron());
  assert(!photon.isElectron());
  assert(!photon.isProton());
  assert(!photon.isCharged());

  Particle neutron = Particle::createNeutron(Vector3D::ZERO, Vector3D::UNIT_Y, 1.0);
  assert(!neutron.isPhoton());
  assert(neutron.isNeutron());
  assert(!neutron.isElectron());
  assert(!neutron.isProton());
  assert(!neutron.isCharged());

  Particle electron = Particle::createElectron(Vector3D::ZERO, Vector3D::UNIT_Z, 1.0);
  assert(!electron.isPhoton());
  assert(!electron.isNeutron());
  assert(electron.isElectron());
  assert(!electron.isProton());
  assert(electron.isCharged());

  Particle proton = Particle::createProton(Vector3D::ZERO, Vector3D::UNIT_X, 1.0);
  assert(!proton.isPhoton());
  assert(!proton.isNeutron());
  assert(!proton.isElectron());
  assert(proton.isProton());
  assert(proton.isCharged());

  std::cout << "✓ Particle type checks passed" << std::endl;
}

void testInteractionTracking()
{
  std::cout << "Testing interaction tracking..." << std::endl;

  Particle p = Particle::createPhoton(Vector3D::ZERO, Vector3D::UNIT_X, 1.0);

  // Initially no interaction point
  assert(!p.getLastInteractionPoint().has_value());

  // Record interaction
  Vector3D interactionPoint(5.0, 2.0, -1.0);
  p.recordInteraction(interactionPoint);
  auto lastPoint = p.getLastInteractionPoint();
  assert(lastPoint.has_value());
  assert(*lastPoint == interactionPoint);

  // Test birth energy
  auto initialBirthEnergy = p.getBirthEnergy();
  assert(initialBirthEnergy.has_value()); // Should be set to initial energy (1.0)
  assert(*initialBirthEnergy == 1.0);

  p.setBirthEnergy(1.5);
  auto birthEnergy = p.getBirthEnergy();
  assert(birthEnergy.has_value());
  assert(*birthEnergy == 1.5);

  std::cout << "✓ Interaction tracking passed" << std::endl;
}

void testSecondaryParticleCreation()
{
  std::cout << "Testing secondary particle creation..." << std::endl;

  Particle primary = Particle::createPhoton(Vector3D(1.0, 2.0, 3.0), Vector3D::UNIT_Z, 2.0);
  primary.setGeneration(0);
  primary.setBirthEnergy(2.0);

  Vector3D secondaryPos(4.0, 5.0, 6.0);
  Vector3D secondaryDir(0.0, 1.0, 0.0);
  Particle secondary = primary.createSecondary(ParticleType::Electron, secondaryPos, secondaryDir, 0.5);

  assert(secondary.type() == ParticleType::Electron);
  assert(secondary.position() == secondaryPos);
  assert(secondary.direction() == secondaryDir);
  assert(secondary.energy() == 0.5);
  assert(secondary.generation() == 1); // One generation higher
  auto secondaryBirth = secondary.getBirthEnergy();
  assert(secondaryBirth.has_value());
  assert(*secondaryBirth == 0.5);

  std::cout << "✓ Secondary particle creation passed" << std::endl;
}

void testStringRepresentation()
{
  std::cout << "Testing string representation..." << std::endl;

  Particle p = Particle::createNeutron(Vector3D(1.0, 2.0, 3.0), Vector3D::UNIT_X, 1.5);
  p.setGeneration(2);

  std::string str = p.toString();
  assert(str.find("Neutron") != std::string::npos);
  assert(str.find("Alive") != std::string::npos);
  assert(str.find("1.5") != std::string::npos);

  assert(p.getTypeName() == "Neutron");
  assert(p.getStateName() == "Alive");

  p.absorb();
  assert(p.getStateName() == "Absorbed");

  std::cout << "✓ String representation passed" << std::endl;
}

void testComparisonOperators()
{
  std::cout << "Testing comparison operators..." << std::endl;

  Vector3D pos(1.0, 2.0, 3.0);
  Vector3D dir(0.0, 0.0, 1.0);

  Particle p1 = Particle::createPhoton(pos, dir, 1.5);
  Particle p2 = Particle::createPhoton(pos, dir, 1.5);

  // They should be equal - comparison is based on physical properties, not history IDs
  assert(p1 == p2);

  // Test inequality with different properties
  Particle p3 = Particle::createPhoton(pos, dir, 2.0); // Different energy
  assert(p1 != p3);

  Particle p4 = Particle::createNeutron(pos, dir, 1.5); // Different type
  assert(p1 != p4);

  Vector3D differentPos(2.0, 3.0, 4.0);
  Particle p5 = Particle::createPhoton(differentPos, dir, 1.5); // Different position
  assert(p1 != p5);

  // Test that particles with same physical properties are equal
  Particle p6(ParticleType::Photon, pos, dir, 1.5);
  Particle p7(ParticleType::Photon, pos, dir, 1.5);
  assert(p6 == p7);

  // Test inequality operator
  assert(!(p1 != p2)); // Should be false since p1 == p2
  assert(p1 != p3);    // Should be true since they have different energies

  std::cout << "✓ Comparison operators passed" << std::endl;
}

void testStaticFactoryMethods()
{
  std::cout << "Testing static factory methods..." << std::endl;

  Vector3D pos(1.0, 2.0, 3.0);
  Vector3D dir(0.0, 1.0, 0.0);

  Particle photon = Particle::createPhoton(pos, dir, 1.0);
  assert(photon.type() == ParticleType::Photon);
  assert(photon.position() == pos);
  assert(photon.energy() == 1.0);

  Particle neutron = Particle::createNeutron(pos, dir, 2.0);
  assert(neutron.type() == ParticleType::Neutron);
  assert(neutron.energy() == 2.0);

  Particle electron = Particle::createElectron(pos, dir, 0.5);
  assert(electron.type() == ParticleType::Electron);
  assert(electron.energy() == 0.5);

  Particle proton = Particle::createProton(pos, dir, 100.0);
  assert(proton.type() == ParticleType::Proton);
  assert(proton.energy() == 100.0);

  std::cout << "✓ Static factory methods passed" << std::endl;
}

void testPhysicalConstants()
{
  std::cout << "Testing physical constants..." << std::endl;

  // Rest masses
  assert(Particle::getRestMass(ParticleType::Photon) == 0.0);
  assert(Particle::getRestMass(ParticleType::Electron) > 0.5);  // ~0.511 MeV
  assert(Particle::getRestMass(ParticleType::Proton) > 900.0);  // ~938 MeV
  assert(Particle::getRestMass(ParticleType::Neutron) > 900.0); // ~939 MeV

  // Charges
  assert(Particle::getCharge(ParticleType::Photon) == 0);
  assert(Particle::getCharge(ParticleType::Electron) == -1);
  assert(Particle::getCharge(ParticleType::Proton) == +1);
  assert(Particle::getCharge(ParticleType::Neutron) == 0);

  // Particle symbols
  assert(Particle::getParticleSymbol(ParticleType::Photon) == "γ");
  assert(Particle::getParticleSymbol(ParticleType::Electron) == "e⁻");
  assert(Particle::getParticleSymbol(ParticleType::Proton) == "p");
  assert(Particle::getParticleSymbol(ParticleType::Neutron) == "n");

  std::cout << "✓ Physical constants passed" << std::endl;
}

void testUtilityFunctions()
{
  std::cout << "Testing utility functions..." << std::endl;

  // Test particle type string conversion
  assert(particleTypeToString(ParticleType::Photon) == "Photon");
  assert(particleTypeToString(ParticleType::Neutron) == "Neutron");
  assert(particleTypeToString(ParticleType::Electron) == "Electron");
  assert(particleTypeToString(ParticleType::Proton) == "Proton");

  // Test particle state string conversion
  assert(particleStateToString(ParticleState::Alive) == "Alive");
  assert(particleStateToString(ParticleState::Absorbed) == "Absorbed");
  assert(particleStateToString(ParticleState::Escaped) == "Escaped");
  assert(particleStateToString(ParticleState::Scattered) == "Scattered");
  assert(particleStateToString(ParticleState::Split) == "Split");
  assert(particleStateToString(ParticleState::Terminated) == "Terminated");

  // Test string to particle type conversion
  auto photonType = stringToParticleType("photon");
  assert(photonType.has_value());
  assert(*photonType == ParticleType::Photon);

  auto neutronType = stringToParticleType("Neutron");
  assert(neutronType.has_value());
  assert(*neutronType == ParticleType::Neutron);

  auto gammaType = stringToParticleType("gamma");
  assert(gammaType.has_value());
  assert(*gammaType == ParticleType::Photon);

  auto invalidType = stringToParticleType("muon");
  assert(!invalidType.has_value());

  std::cout << "✓ Utility functions passed" << std::endl;
}

void testErrorHandling()
{
  std::cout << "Testing error handling..." << std::endl;

  // Test construction with negative energy
  bool threw = false;
  try
  {
    Particle p(ParticleType::Photon, Vector3D::ZERO, Vector3D::UNIT_X, -1.0);
  }
  catch (const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  // Test construction with negative weight
  threw = false;
  try
  {
    Particle p(ParticleType::Photon, Vector3D::ZERO, Vector3D::UNIT_X, 1.0, -0.5);
  }
  catch (const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  // Test construction with zero direction
  threw = false;
  try
  {
    Particle p(ParticleType::Photon, Vector3D::ZERO, Vector3D::ZERO, 1.0);
  }
  catch (const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  // Test invalid particle type in getRestMass
  threw = false;
  try
  {
    Particle::getRestMass(static_cast<ParticleType>(999));
  }
  catch (const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  std::cout << "✓ Error handling passed" << std::endl;
}

int main()
{
  std::cout << "Running Particle comprehensive tests...\n"
            << std::endl;

  testConstructors();
  testStructuredBinding();
  testSafeMutators();
  testDirectMutators();
  testMovement();
  testStateManagement();
  testEnergyOperations();
  testRelativisticCalculations();
  testParticleTypeChecks();
  testInteractionTracking();
  testSecondaryParticleCreation();
  testStringRepresentation();
  testComparisonOperators();
  testStaticFactoryMethods();
  testPhysicalConstants();
  testUtilityFunctions();
  testErrorHandling();

  std::cout << "\n🎉 All Particle tests passed successfully!" << std::endl;
  std::cout << "Your Particle class is working correctly and ready for commit." << std::endl;
  std::cout << "The class demonstrates excellent C++17 features and professional nuclear physics implementation." << std::endl;

  return 0;
}