#define _USE_MATH_DEFINES // M_PI is used

#include "MonteCarloSampling.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

constexpr double TOLERANCE = 1e-10;
constexpr double STATISTICAL_TOLERANCE =
    0.1; // 10% tolerance for statistical tests

void testEnergySpectrumFactories()
{
  std::cout << "Testing EnergySpectrum factory methods..." << std::endl;

  // Test monoenergetic
  auto mono = EnergySpectrum::monoenergetic(1.5);
  assert(mono.type == SpectrumType::Monoenergetic);
  assert(std::abs(mono.characteristic_energy - 1.5) < TOLERANCE);

  // Test uniform
  auto uniform = EnergySpectrum::uniform(0.5, 2.0);
  assert(uniform.type == SpectrumType::Uniform);
  assert(std::abs(uniform.min_energy - 0.5) < TOLERANCE);
  assert(std::abs(uniform.max_energy - 2.0) < TOLERANCE);

  // Test Maxwell
  auto maxwell = EnergySpectrum::maxwell(1.0);
  assert(maxwell.type == SpectrumType::Maxwell);
  assert(std::abs(maxwell.characteristic_energy - 1.0) < TOLERANCE);

  // Test exponential
  auto exponential = EnergySpectrum::exponential(2.0);
  assert(exponential.type == SpectrumType::Exponential);
  assert(std::abs(exponential.characteristic_energy - 2.0) < TOLERANCE);

  // Test power law
  auto powerlaw = EnergySpectrum::powerLaw(-2.0, 1.0, 10.0);
  assert(powerlaw.type == SpectrumType::PowerLaw);
  assert(std::abs(powerlaw.parameter + 2.0) < TOLERANCE);

  // Test discrete
  std::vector<double> energies = {1.0, 2.0, 3.0};
  std::vector<double> weights = {0.5, 0.3, 0.2};
  auto discrete = EnergySpectrum::discrete(energies, weights);
  assert(discrete.type == SpectrumType::Discrete);
  assert(discrete.discrete_energies.size() == 3);
  assert(discrete.discrete_weights.size() == 3);

  // Test error handling
  bool threw = false;
  try
  {
    EnergySpectrum::uniform(2.0, 1.0); // Invalid range
  }
  catch(const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  std::cout << "✓ EnergySpectrum factory methods passed" << std::endl;
}

void testAngularDistributionFactories()
{
  std::cout << "Testing AngularDistribution factory methods..." << std::endl;

  // Test isotropic
  auto isotropic = AngularDistribution::isotropic();
  assert(isotropic.type == DirectionType::Isotropic);

  // Test beam
  Vector3D beam_dir(1.0, 0.0, 0.0);
  auto beam = AngularDistribution::beam(beam_dir);
  assert(beam.type == DirectionType::Beam);
  assert(beam.reference_direction == beam_dir);

  // Test cone
  Vector3D cone_axis(0.0, 0.0, 1.0);
  double half_angle = 0.5;
  auto cone = AngularDistribution::cone(cone_axis, half_angle);
  assert(cone.type == DirectionType::Cone);
  assert(cone.reference_direction == cone_axis);
  assert(std::abs(cone.cone_angle - half_angle) < TOLERANCE);

  // Test cosine weighted
  Vector3D normal(0.0, 1.0, 0.0);
  auto cosine = AngularDistribution::cosineWeighted(normal);
  assert(cosine.type == DirectionType::CosineWeighted);
  assert(cosine.reference_direction == normal);

  // Test error handling
  bool threw = false;
  try
  {
    AngularDistribution::beam(Vector3D::ZERO); // Zero vector
  }
  catch(const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  std::cout << "✓ AngularDistribution factory methods passed" << std::endl;
}

void testExponentialSampling()
{
  std::cout << "Testing exponential sampling for interaction distances..."
            << std::endl;

  RandomNumberGenerator rng(12345);
  MonteCarloSampling mc(rng);

  // Test basic exponential sampling
  double lambda = 2.0;
  std::vector<double> samples;

  for(int i = 0; i < 10000; ++i)
  {
    double distance = mc.sampleInteractionDistance(lambda);
    assert(distance >= 0.0);
    samples.push_back(distance);
  }

  // Check mean is approximately 1/lambda
  double mean =
      std::accumulate(samples.begin(), samples.end(), 0.0) / samples.size();
  double expected_mean = 1.0 / lambda;
  assert(std::abs(mean - expected_mean) < STATISTICAL_TOLERANCE);

  // Test safe exponential sampling
  auto safe_sample = mc.tryExponentialSampling(1.0);
  assert(safe_sample.has_value());
  assert(*safe_sample >= 0.0);

  auto invalid_sample = mc.tryExponentialSampling(-1.0);
  assert(!invalid_sample.has_value());

  // Test total interaction distance
  double total_cross_section = 1.5;
  double total_distance =
      mc.sampleTotalInteractionDistance(total_cross_section);
  assert(total_distance >= 0.0);

  // Test error handling
  bool threw = false;
  try
  {
    mc.sampleInteractionDistance(-1.0);
  }
  catch(const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  std::cout << "✓ Exponential sampling passed" << std::endl;
}

void testDirectionSampling()
{
  std::cout << "Testing uniform spherical direction sampling..." << std::endl;

  RandomNumberGenerator rng(23456);
  MonteCarloSampling mc(rng);

  // Test isotropic direction sampling
  std::vector<Vector3D> directions;
  for(int i = 0; i < 1000; ++i)
  {
    Vector3D dir = mc.sampleIsotropicDirection();
    assert(dir.isUnit(1e-10));
    directions.push_back(dir);
  }

  // Check isotropy - average should be close to zero
  Vector3D average(0, 0, 0);
  for(const auto &dir : directions)
  {
    average += dir;
  }
  average = average * (1.0 / directions.size());
  assert(average.magnitude() < 0.1); // Should be small for isotropic

  // Test angular distribution sampling
  auto isotropic_dist = AngularDistribution::isotropic();
  Vector3D iso_dir = mc.sampleDirection(isotropic_dist);
  assert(iso_dir.isUnit(1e-10));

  auto beam_dist = AngularDistribution::beam(Vector3D::UNIT_Z);
  Vector3D beam_dir = mc.sampleDirection(beam_dist);
  assert(beam_dir == Vector3D::UNIT_Z);

  // Test cone direction sampling
  Vector3D cone_axis(0.0, 0.0, 1.0);
  double half_angle = 0.5;
  Vector3D cone_dir = mc.sampleConeDirection(cone_axis, half_angle);
  assert(cone_dir.isUnit(1e-10));

  // Check angle constraint
  double angle = std::acos(cone_dir.dot(cone_axis));
  assert(angle <= half_angle + 1e-10);

  // Test cosine weighted direction
  Vector3D normal(0.0, 0.0, 1.0);
  Vector3D cosine_dir = mc.sampleCosineDirection(normal);
  assert(cosine_dir.isUnit(1e-10));
  assert(cosine_dir.dot(normal) >= 0.0); // Should be in hemisphere

  // Test hemisphere direction
  Vector3D hemisphere_dir = mc.sampleHemisphereDirection(normal);
  assert(hemisphere_dir.isUnit(1e-10));
  assert(hemisphere_dir.dot(normal) >= 0.0);

  std::cout << "✓ Direction sampling passed" << std::endl;
}

void testEnergySampling()
{
  std::cout << "Testing energy distribution sampling..." << std::endl;

  RandomNumberGenerator rng(34567);
  MonteCarloSampling mc(rng);

  // Test monoenergetic sampling
  double fixed_energy = 1.5;
  auto mono_spectrum = EnergySpectrum::monoenergetic(fixed_energy);

  for(int i = 0; i < 100; ++i)
  {
    double energy = mc.sampleEnergy(mono_spectrum);
    assert(std::abs(energy - fixed_energy) < TOLERANCE);
  }

  // Test uniform energy sampling
  double min_energy = 0.5;
  double max_energy = 2.0;
  auto uniform_spectrum = EnergySpectrum::uniform(min_energy, max_energy);

  std::vector<double> uniform_energies;
  for(int i = 0; i < 1000; ++i)
  {
    double energy = mc.sampleEnergy(uniform_spectrum);
    assert(energy >= min_energy);
    assert(energy <= max_energy);
    uniform_energies.push_back(energy);
  }

  // Check uniform distribution mean
  double uniform_mean =
      std::accumulate(uniform_energies.begin(), uniform_energies.end(), 0.0) /
      uniform_energies.size();
  double expected_uniform_mean = (min_energy + max_energy) / 2.0;
  assert(std::abs(uniform_mean - expected_uniform_mean) <
         STATISTICAL_TOLERANCE);

  // Test Maxwell energy sampling
  double temperature = 1.0;
  auto maxwell_spectrum = EnergySpectrum::maxwell(temperature);

  for(int i = 0; i < 100; ++i)
  {
    double energy = mc.sampleEnergy(maxwell_spectrum);
    assert(energy >= 0.0);
  }

  // Test exponential energy sampling
  double char_energy = 2.0;
  auto exp_spectrum = EnergySpectrum::exponential(char_energy);

  for(int i = 0; i < 100; ++i)
  {
    double energy = mc.sampleEnergy(exp_spectrum);
    assert(energy >= 0.0);
  }

  // Test power law energy sampling
  double alpha = -2.0;
  auto power_spectrum = EnergySpectrum::powerLaw(alpha, 1.0, 10.0);

  for(int i = 0; i < 100; ++i)
  {
    double energy = mc.sampleEnergy(power_spectrum);
    assert(energy >= 1.0);
    assert(energy <= 10.0);
  }

  // Test discrete energy sampling
  std::vector<double> energies = {1.0, 2.0, 3.0};
  std::vector<double> weights = {0.5, 0.3, 0.2};
  auto discrete_spectrum = EnergySpectrum::discrete(energies, weights);

  for(int i = 0; i < 100; ++i)
  {
    double energy = mc.sampleEnergy(discrete_spectrum);
    assert(std::find(energies.begin(), energies.end(), energy) !=
           energies.end());
  }

  // Test direct sampling methods
  assert(mc.sampleMonoenergeticEnergy(1.0) == 1.0);

  double uniform_energy = mc.sampleUniformEnergy(1.0, 2.0);
  assert(uniform_energy >= 1.0 && uniform_energy <= 2.0);

  double maxwell_energy = mc.sampleMaxwellEnergy(1.0);
  assert(maxwell_energy >= 0.0);

  std::cout << "✓ Energy sampling passed" << std::endl;
}

void testAdvancedTechniques()
{
  std::cout << "Testing advanced Monte Carlo techniques..." << std::endl;

  RandomNumberGenerator rng(45678);
  MonteCarloSampling mc(rng);

  // Test azimuthal angle sampling
  for(int i = 0; i < 100; ++i)
  {
    double phi = mc.sampleAzimuthalAngle();
    assert(phi >= 0.0);
    assert(phi < 2.0 * M_PI);
  }

  // Test scattering application
  Vector3D incident(1.0, 0.0, 0.0);
  double polar = M_PI / 4.0;
  double azimuthal = M_PI / 2.0;

  Vector3D scattered = mc.applyScattering(incident, polar, azimuthal);
  assert(scattered.isUnit(1e-10));

  // Test interaction type sampling
  std::vector<double> cross_sections = {1.0, 2.0, 3.0};
  for(int i = 0; i < 100; ++i)
  {
    std::size_t type = mc.sampleInteractionType(cross_sections);
    assert(type < cross_sections.size());
  }

  // Test Russian roulette
  double survival_prob = 0.8;
  int survivors = 0;
  int trials = 1000;

  for(int i = 0; i < trials; ++i)
  {
    if(mc.russianRoulette(survival_prob))
    {
      survivors++;
    }
  }

  double actual_survival = static_cast<double>(survivors) / trials;
  assert(std::abs(actual_survival - survival_prob) < 0.1); // Within 10%

  // Test particle splitting
  double splitting_factor = 2.5;
  for(int i = 0; i < 100; ++i)
  {
    std::size_t count = mc.particleSplitting(splitting_factor);
    assert(count >= 2); // At least floor(2.5) = 2
    assert(count <= 3); // At most floor(2.5) + 1 = 3
  }

  std::cout << "✓ Advanced techniques passed" << std::endl;
}

void testPhysicsSpecificSampling()
{
  std::cout << "Testing physics-specific sampling..." << std::endl;

  RandomNumberGenerator rng(56789);
  MonteCarloSampling mc(rng);

  // Test Compton scattering
  double photon_energy = 1.0; // MeV

  for(int i = 0; i < 100; ++i)
  {
    double compton_angle = mc.sampleComptonAngle(photon_energy);
    assert(compton_angle >= 0.0);
    assert(compton_angle <= M_PI);

    double scattered_energy =
        mc.sampleComptonEnergy(photon_energy, compton_angle);
    assert(scattered_energy > 0.0);
    assert(scattered_energy <= photon_energy); // Energy can't increase
  }

  // Test photoelectron direction
  for(int i = 0; i < 100; ++i)
  {
    Vector3D pe_dir = mc.samplePhotoelectronDirection();
    assert(pe_dir.isUnit(1e-10));
  }

  // Test pair production
  double high_energy = 5.0; // MeV
  auto [electron_dir, positron_dir] =
      mc.samplePairProductionDirections(high_energy);
  assert(electron_dir.isUnit(1e-10));
  assert(positron_dir.isUnit(1e-10));

  std::cout << "✓ Physics-specific sampling passed" << std::endl;
}

void testValidation()
{
  std::cout << "Testing validation methods..." << std::endl;

  RandomNumberGenerator rng(67890);
  MonteCarloSampling mc(rng);

  // Test valid energy spectrum
  auto valid_spectrum = EnergySpectrum::uniform(1.0, 2.0);
  assert(mc.validateEnergySpectrum(valid_spectrum));

  // Test invalid energy spectrum
  EnergySpectrum invalid_spectrum;
  invalid_spectrum.min_energy = -1.0; // Invalid
  assert(!mc.validateEnergySpectrum(invalid_spectrum));

  // Test valid angular distribution
  auto valid_dist = AngularDistribution::isotropic();
  assert(mc.validateAngularDistribution(valid_dist));

  // Test invalid angular distribution
  AngularDistribution invalid_dist;
  invalid_dist.reference_direction = Vector3D(2.0, 0.0, 0.0); // Not unit vector
  assert(!mc.validateAngularDistribution(invalid_dist));

  std::cout << "✓ Validation passed" << std::endl;
}

void testBatchSampling()
{
  std::cout << "Testing batch sampling methods..." << std::endl;

  RandomNumberGenerator rng(78901);
  MonteCarloSampling mc(rng);

  std::size_t count = 1000;

  // Test batch interaction distances
  double lambda = 1.5;
  auto distances = mc.sampleInteractionDistances(lambda, count);
  assert(distances.size() == count);

  for(double distance : distances)
  {
    assert(distance >= 0.0);
  }

  // Test batch isotropic directions
  auto directions = mc.sampleIsotropicDirections(count);
  assert(directions.size() == count);

  for(const Vector3D &dir : directions)
  {
    assert(dir.isUnit(1e-10));
  }

  // Test batch energies
  auto spectrum = EnergySpectrum::uniform(1.0, 2.0);
  auto energies = mc.sampleEnergies(spectrum, count);
  assert(energies.size() == count);

  for(double energy : energies)
  {
    assert(energy >= 1.0);
    assert(energy <= 2.0);
  }

  std::cout << "✓ Batch sampling passed" << std::endl;
}

void testUtilityFunctions()
{
  std::cout << "Testing utility functions..." << std::endl;

  // Test Klein-Nishina cross-section
  double incident_energy = 1.0;
  double scattering_angle = M_PI / 2.0;
  double cross_section = MonteCarloUtils::kleinNishinaCrossSection(
      incident_energy, scattering_angle);
  assert(cross_section > 0.0);

  // Test coordinate transformations
  double theta = M_PI / 4.0;
  double phi = M_PI / 3.0;
  Vector3D cartesian = MonteCarloUtils::polarToCartesian(theta, phi);
  assert(cartesian.isUnit(1e-10));

  auto [theta_back, phi_back] = MonteCarloUtils::cartesianToPolar(cartesian);
  assert(std::abs(theta_back - theta) < TOLERANCE);
  assert(std::abs(phi_back - phi) < TOLERANCE);

  // Test local to global transformation
  Vector3D local(1.0, 0.0, 0.0);
  Vector3D local_z(0.0, 0.0, 1.0);
  Vector3D global = MonteCarloUtils::localToGlobal(local, local_z);
  assert(global.isUnit(1e-10));

  // Test orthonormal basis generation
  Vector3D z_axis(1.0, 1.0, 0.0);
  z_axis = z_axis.normalise();
  auto [x_axis, y_axis] = MonteCarloUtils::generateOrthonormalBasis(z_axis);

  assert(x_axis.isUnit(1e-10));
  assert(y_axis.isUnit(1e-10));
  assert(std::abs(x_axis.dot(y_axis)) < TOLERANCE); // Orthogonal
  assert(std::abs(x_axis.dot(z_axis)) < TOLERANCE); // Orthogonal
  assert(std::abs(y_axis.dot(z_axis)) < TOLERANCE); // Orthogonal

  std::cout << "✓ Utility functions passed" << std::endl;
}

void testErrorHandling()
{
  std::cout << "Testing error handling..." << std::endl;

  RandomNumberGenerator rng(89012);
  MonteCarloSampling mc(rng);

  bool threw = false;

  // Test invalid energy parameters
  try
  {
    mc.sampleMonoenergeticEnergy(-1.0);
  }
  catch(const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  threw = false;
  try
  {
    mc.sampleUniformEnergy(2.0, 1.0); // Invalid range
  }
  catch(const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  threw = false;
  try
  {
    mc.sampleMaxwellEnergy(-1.0); // Invalid temperature
  }
  catch(const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  // Test invalid cone parameters
  threw = false;
  try
  {
    mc.sampleConeDirection(Vector3D::UNIT_Z, -1.0); // Invalid angle
  }
  catch(const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  // Test empty discrete energy vectors
  threw = false;
  try
  {
    std::vector<double> empty_energies;
    std::vector<double> empty_weights;
    mc.sampleDiscreteEnergy(empty_energies, empty_weights);
  }
  catch(const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  std::cout << "✓ Error handling passed" << std::endl;
}

void testIntegrationWithExistingClasses()
{
  std::cout << "Testing integration with existing classes..." << std::endl;

  RandomNumberGenerator rng(90123);
  MonteCarloSampling mc(rng);

  // Test with Material class
  Material water = Material::createWater();
  double particle_energy = 1.0;

  // Create a simple cross-section function
  auto cross_section_func = [](const Material &mat, double energy) -> double {
    return 0.1 * mat.density(); // Simple linear relationship
  };

  double distance =
      mc.sampleInteractionDistance(water, particle_energy, cross_section_func);
  assert(distance > 0.0);

  // Test with Vector3D class
  Vector3D test_dir(1.0, 1.0, 1.0);
  test_dir = test_dir.normalise();

  Vector3D cone_sample = mc.sampleConeDirection(test_dir, 0.1);
  assert(cone_sample.isUnit(1e-10));

  double angle_between = std::acos(cone_sample.dot(test_dir));
  assert(angle_between <= 0.1 + 1e-10);

  std::cout << "✓ Integration with existing classes passed" << std::endl;
}

int main()
{
  std::cout << "Running MonteCarloSampling comprehensive tests...\n"
            << std::endl;

  testEnergySpectrumFactories();
  testAngularDistributionFactories();
  testExponentialSampling();
  testDirectionSampling();
  testEnergySampling();
  testAdvancedTechniques();
  testPhysicsSpecificSampling();
  testValidation();
  testBatchSampling();
  testUtilityFunctions();
  testErrorHandling();
  testIntegrationWithExistingClasses();

  std::cout << "\n🎉 All MonteCarloSampling tests passed successfully!"
            << std::endl;
  std::cout
      << "✓ Exponential sampling for interaction distances working perfectly!"
      << std::endl;
  std::cout << "✓ Uniform spherical direction sampling implemented and tested!"
            << std::endl;
  std::cout << "✓ Energy distribution sampling covering all major spectra!"
            << std::endl;
  std::cout << "✓ Physics-specific sampling (Compton, photoelectric, pair "
               "production) ready!"
            << std::endl;
  std::cout << "✓ Advanced Monte Carlo techniques (Russian roulette, "
               "splitting) implemented!"
            << std::endl;
  std::cout << "\nDay 4 Monte Carlo Foundation COMPLETE! 🚀" << std::endl;
  std::cout << "Ready for photon physics implementation next!" << std::endl;

  return 0;
}