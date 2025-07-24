#pragma once

#include "RandomNumberGenerator.hpp"
#include "Vector3D.hpp"

#include <cmath>
#include <functional>
#include <optional>
#include <vector>

/**
 * @brief Simple energy spectrum for utility sampling
 */
struct EnergySpectrum
{
  enum class Type
  {
    Monoenergetic,
    Uniform,
    Exponential,
    Maxwell,
    Discrete
  };

  Type type = Type::Monoenergetic;
  double min_energy = 0.0;
  double max_energy = 1.0;
  double characteristic_energy = 1.0;
  std::vector<double> discrete_energies;
  std::vector<double> discrete_weights;
};

/**
 * @brief Monte Carlo sampling utility functions
 *
 * These are the core sampling functions needed for Day 4 of the
 * Monte Carlo radiation transport simulator development.
 */
namespace MonteCarloUtils
{
/**
 * @brief Sample from exponential distribution for interaction distances
 * @param lambda Rate parameter (1/mean)
 * @return Sampled distance
 */
inline double sampleExponential(double lambda)
{
  if(lambda <= 0.0)
  {
    throw std::invalid_argument("Lambda must be positive");
  }

  auto &rng = RandomNumberGenerator::getThreadLocal();
  return rng.exponential(lambda);
}

/**
 * @brief Sample isotropic direction uniformly over 4π steradians
 * @return Unit vector with isotropic direction
 */
inline Vector3D sampleIsotropicDirection()
{
  auto &rng = RandomNumberGenerator::getThreadLocal();
  return rng.isotropicDirection();
}

/**
 * @brief Sample energy from spectrum
 * @param spectrum Energy spectrum definition
 * @return Sampled energy (MeV)
 */
inline double sampleEnergy(const EnergySpectrum &spectrum)
{
  auto &rng = RandomNumberGenerator::getThreadLocal();

  switch(spectrum.type)
  {
  case EnergySpectrum::Type::Monoenergetic:
    return spectrum.characteristic_energy;

  case EnergySpectrum::Type::Uniform:
    return rng.uniform(spectrum.min_energy, spectrum.max_energy);

  case EnergySpectrum::Type::Exponential:
    return rng.exponential(1.0 / spectrum.characteristic_energy);

  case EnergySpectrum::Type::Maxwell:
    return rng.maxwellBoltzmann(spectrum.characteristic_energy);

  case EnergySpectrum::Type::Discrete:
    if(spectrum.discrete_energies.empty())
    {
      return spectrum.characteristic_energy;
    }
    {
      std::size_t index = rng.discreteSample(spectrum.discrete_weights);
      return spectrum.discrete_energies[index];
    }

  default:
    return spectrum.characteristic_energy;
  }
}

// Helper factory functions for energy spectra
inline EnergySpectrum createMonoenergeticSpectrum(double energy)
{
  EnergySpectrum spectrum;
  spectrum.type = EnergySpectrum::Type::Monoenergetic;
  spectrum.characteristic_energy = energy;
  return spectrum;
}

inline EnergySpectrum createUniformSpectrum(double min_energy,
                                            double max_energy)
{
  EnergySpectrum spectrum;
  spectrum.type = EnergySpectrum::Type::Uniform;
  spectrum.min_energy = min_energy;
  spectrum.max_energy = max_energy;
  return spectrum;
}

inline EnergySpectrum createExponentialSpectrum(double characteristic_energy)
{
  EnergySpectrum spectrum;
  spectrum.type = EnergySpectrum::Type::Exponential;
  spectrum.characteristic_energy = characteristic_energy;
  return spectrum;
}
} // namespace MonteCarloUtils