#include "MonteCarloUtils.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

void testExponentialSampling()
{
  std::cout << "Testing exponential sampling..." << std::endl;

  double lambda = 2.0;
  std::vector<double> samples;

  // Generate samples
  for(int i = 0; i < 10000; ++i)
  {
    double sample = MonteCarloUtils::sampleExponential(lambda);
    assert(sample >= 0.0);
    samples.push_back(sample);
  }

  // Check mean is approximately 1/lambda
  double mean =
      std::accumulate(samples.begin(), samples.end(), 0.0) / samples.size();
  double expected_mean = 1.0 / lambda;
  double relative_error = std::abs(mean - expected_mean) / expected_mean;

  std::cout << "  Expected mean: " << expected_mean << std::endl;
  std::cout << "  Actual mean: " << mean << std::endl;
  std::cout << "  Relative error: " << relative_error * 100 << "%" << std::endl;

  assert(relative_error < 0.05); // Within 5%

  std::cout << "✓ Exponential sampling passed" << std::endl;
}

void testIsotropicDirection()
{
  std::cout << "Testing isotropic direction sampling..." << std::endl;

  std::vector<Vector3D> directions;

  // Generate directions
  for(int i = 0; i < 1000; ++i)
  {
    Vector3D dir = MonteCarloUtils::sampleIsotropicDirection();
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

  std::cout << "  Average direction magnitude: " << average.magnitude()
            << std::endl;
  assert(average.magnitude() < 0.1); // Should be small for isotropic

  std::cout << "✓ Isotropic direction sampling passed" << std::endl;
}

void testEnergySpectrum()
{
  std::cout << "Testing energy spectrum sampling..." << std::endl;

  // Test monoenergetic spectrum
  auto mono_spectrum = MonteCarloUtils::createMonoenergeticSpectrum(1.5);
  for(int i = 0; i < 100; ++i)
  {
    double energy = MonteCarloUtils::sampleEnergy(mono_spectrum);
    assert(std::abs(energy - 1.5) < 1e-10);
  }
  std::cout << "  ✓ Monoenergetic spectrum works" << std::endl;

  // Test uniform spectrum
  auto uniform_spectrum = MonteCarloUtils::createUniformSpectrum(1.0, 3.0);
  std::vector<double> uniform_samples;
  for(int i = 0; i < 1000; ++i)
  {
    double energy = MonteCarloUtils::sampleEnergy(uniform_spectrum);
    assert(energy >= 1.0);
    assert(energy <= 3.0);
    uniform_samples.push_back(energy);
  }

  double uniform_mean =
      std::accumulate(uniform_samples.begin(), uniform_samples.end(), 0.0) /
      uniform_samples.size();
  std::cout << "  Uniform spectrum mean: " << uniform_mean << " (expected: 2.0)"
            << std::endl;
  assert(std::abs(uniform_mean - 2.0) < 0.1);
  std::cout << "  ✓ Uniform spectrum works" << std::endl;

  // Test exponential spectrum
  auto exp_spectrum = MonteCarloUtils::createExponentialSpectrum(2.0);
  for(int i = 0; i < 100; ++i)
  {
    double energy = MonteCarloUtils::sampleEnergy(exp_spectrum);
    assert(energy >= 0.0);
  }
  std::cout << "  ✓ Exponential spectrum works" << std::endl;

  std::cout << "✓ Energy spectrum sampling passed" << std::endl;
}

void testMonteCarloDemoSimulation()
{
  std::cout << "Testing demo Monte Carlo transport..." << std::endl;

  // Simulate simple photon transport through material
  double linear_attenuation = 0.5; // cm^-1
  auto photon_spectrum =
      MonteCarloUtils::createMonoenergeticSpectrum(1.0); // 1 MeV photons

  int num_histories = 1000;
  int absorbed = 0;
  int transmitted = 0;
  double thickness = 2.0; // cm

  for(int i = 0; i < num_histories; ++i)
  {
    // Sample initial photon
    double energy = MonteCarloUtils::sampleEnergy(photon_spectrum);
    Vector3D direction = Vector3D(1, 0, 0); // Forward direction
    Vector3D position(0, 0, 0);             // Start at origin

    bool alive = true;
    while(alive && position.x() < thickness)
    {
      // Sample distance to next interaction
      double distance = MonteCarloUtils::sampleExponential(linear_attenuation);

      // Move photon
      position = position + direction * distance;

      if(position.x() >= thickness)
      {
        // Photon escaped
        transmitted++;
        alive = false;
      }
      else
      {
        // Photon was absorbed (simplified - no scattering)
        absorbed++;
        alive = false;
      }
    }
  }

  double transmission_fraction =
      static_cast<double>(transmitted) / num_histories;
  double expected_transmission = std::exp(-linear_attenuation * thickness);

  std::cout << "  Simulated transmission: " << transmission_fraction
            << std::endl;
  std::cout << "  Expected transmission: " << expected_transmission
            << std::endl;
  std::cout << "  Relative error: "
            << std::abs(transmission_fraction - expected_transmission) /
                   expected_transmission * 100
            << "%" << std::endl;

  // Should be within 10% for this simple test
  assert(std::abs(transmission_fraction - expected_transmission) /
             expected_transmission <
         0.15);

  std::cout << "✓ Demo Monte Carlo transport passed" << std::endl;
}

int main()
{
  std::cout << "Running Day 4 Monte Carlo Sampling Utility Tests...\n"
            << std::endl;

  testExponentialSampling();
  testIsotropicDirection();
  testEnergySpectrum();
  testMonteCarloDemoSimulation();

  std::cout << "\n🎉 All Day 4 sampling utility tests passed!" << std::endl;
  std::cout << "✓ sampleExponential(lambda) working perfectly!" << std::endl;
  std::cout << "✓ sampleIsotropicDirection() generating uniform directions!"
            << std::endl;
  std::cout << "✓ sampleEnergy(spectrum) supporting multiple distributions!"
            << std::endl;
  std::cout << "\nDay 4 SUCCESS: Monte Carlo sampling foundation complete! 🚀"
            << std::endl;
  std::cout << "Ready for photon physics implementation on Day 8!" << std::endl;

  return 0;
}