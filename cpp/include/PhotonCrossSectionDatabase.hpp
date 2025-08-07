#pragma once

#include <array>
#include <map>
#include <vector>

/**
 * @brief Comprehensive photon cross-section database with NIST-level accuracy
 * 
 * Implementation based on:
 * - NIST XCOM database
 * - Scofield photoelectric calculations
 * - Hubbell-Øverøø atomic form factors
 * - Bethe-Heitler pair production with screening
 */

namespace PhotonCrossSections
{

/**
 * @brief Atomic binding energies for K, L, M shells (keV)
 * Data from NIST X-ray transition energies database
 */
struct AtomicBindingEnergies
{
  double K_edge = 0.0;     // K absorption edge (keV)
  double L1_edge = 0.0;    // L1 absorption edge (keV)
  double L2_edge = 0.0;    // L2 absorption edge (keV) 
  double L3_edge = 0.0;    // L3 absorption edge (keV)
  double M1_edge = 0.0;    // M1 absorption edge (keV)
  // Additional edges can be added as needed
};

/**
 * @brief Atomic form factor data for coherent scattering
 * Based on Hubbell and Øverøø calculations
 */
struct AtomicFormFactor
{
  double a[4];  // Form factor coefficients a1, a2, a3, a4
  double b[4];  // Form factor coefficients b1, b2, b3, b4
  double c;     // Additional parameter
};

/**
 * @brief Photoelectric cross-section data based on Scofield calculations
 * Provides accurate cross-sections with proper shell structure
 */
class ScofieldPhotoelectricData
{
private:
  // Pre-computed photoelectric data for elements Z=1-100
  static const std::array<std::vector<std::pair<double, double>>, 100> PHOTOELECTRIC_DATA;
  static const std::array<AtomicBindingEnergies, 100> BINDING_ENERGIES;

public:
  /**
   * @brief Get photoelectric cross-section using Scofield calculations
   * @param Z Atomic number
   * @param energy Photon energy (keV)
   * @return Cross-section per atom (barns)
   */
  static double getPhotoelectricCrossSection(int Z, double energy_keV);
  
  /**
   * @brief Get photoelectric cross-section above K-edge
   * @param Z Atomic number
   * @param energy Photon energy (keV)
   * @return Cross-section per atom (barns)
   */
  static double getPhotoelectricCrossSectionAboveKEdge(int Z, double energy_keV);
  
  /**
   * @brief Check if energy is above absorption edge
   * @param Z Atomic number
   * @param energy Photon energy (keV)
   * @param shell Shell designation ('K', 'L1', 'L2', 'L3', 'M1')
   * @return True if above edge
   */
  static bool isAboveAbsorptionEdge(int Z, double energy_keV, char shell);
  
  /**
   * @brief Get binding energy for specific shell
   * @param Z Atomic number
   * @param shell Shell designation
   * @return Binding energy (keV)
   */
  static double getBindingEnergy(int Z, char shell);

private:
  /**
   * @brief Calculate K-shell photoelectric cross-section
   */
  static double calculateKShellCrossSection(int Z, double energy_keV, double edge_keV);
  
  /**
   * @brief Calculate L-shell photoelectric cross-section
   */
  static double calculateLShellCrossSection(int Z, double energy_keV, double edge_keV, int subshell);
  
  /**
   * @brief Calculate M-shell photoelectric cross-section
   */
  static double calculateMShellCrossSection(int Z, double energy_keV, double edge_keV);
  /**
   * @brief Interpolate photoelectric data
   * @param data Energy-cross-section pairs
   * @param energy Target energy
   * @return Interpolated cross-section
   */
  static double interpolateLogLog(const std::vector<std::pair<double, double>>& data, 
                                  double energy);
};

/**
 * @brief Coherent (Rayleigh) scattering with atomic form factors
 * Based on Hubbell and Øverøø tabulations
 */
class CoherentScatteringData  
{
private:
  // Atomic form factor data for Z=1-100
  static const std::array<AtomicFormFactor, 100> FORM_FACTOR_DATA;

public:
  /**
   * @brief Get coherent scattering cross-section with form factors
   * @param Z Atomic number
   * @param energy Photon energy (keV)
   * @return Cross-section per atom (barns)
   */
  static double getCoherentCrossSection(int Z, double energy_keV);
  
  /**
   * @brief Get atomic form factor at momentum transfer
   * @param Z Atomic number
   * @param momentum_transfer q = (sin(θ/2))/λ (Å⁻¹)
   * @return Form factor F(q)
   */
  static double getAtomicFormFactor(int Z, double momentum_transfer);
  
  /**
   * @brief Calculate Thomson scattering cross-section
   * @return Thomson cross-section (barns)
   */
  static double getThomsonCrossSection();

private:
  /**
   * @brief Calculate form factor using analytical approximation
   * @param ff Form factor parameters
   * @param q Momentum transfer
   * @return Form factor value
   */
  static double calculateFormFactor(const AtomicFormFactor& ff, double q);
};

/**
 * @brief Pair production cross-sections with screening corrections
 * Based on Bethe-Heitler theory with Coulomb corrections
 */
class PairProductionData
{
public:
  /**
   * @brief Get pair production cross-section in nuclear field
   * @param Z Atomic number
   * @param energy Photon energy (MeV)
   * @return Cross-section per atom (barns)
   */
  static double getPairProductionNuclearCrossSection(int Z, double energy_MeV);
  
  /**
   * @brief Get pair production cross-section in electron field
   * @param Z Atomic number  
   * @param energy Photon energy (MeV)
   * @return Cross-section per atom (barns)
   */
  static double getPairProductionElectronCrossSection(int Z, double energy_MeV);
  
  /**
   * @brief Get total pair production cross-section
   * @param Z Atomic number
   * @param energy Photon energy (MeV)
   * @return Cross-section per atom (barns)
   */
  static double getTotalPairProductionCrossSection(int Z, double energy_MeV);
  
  /**
   * @brief Check if pair production is energetically possible
   * @param energy Photon energy (MeV)
   * @return True if energy > 1.022 MeV
   */
  static bool isPairProductionPossible(double energy_MeV);

private:
  /**
   * @brief Calculate screening functions for pair production
   * @param Z Atomic number
   * @param energy Photon energy (MeV)
   * @return Screening correction factors
   */
  static std::pair<double, double> calculateScreeningFunctions(int Z, double energy_MeV);
  
  /**
   * @brief Calculate Coulomb correction factor
   * @param Z Atomic number
   * @return Coulomb correction
   */
  static double calculateCoulombCorrection(int Z);
};

/**
 * @brief Incoherent (Compton) scattering functions
 * Hubbell-Øverøø tabulated incoherent scattering functions
 */
class IncoherentScatteringData
{
private:
  // Incoherent scattering function data
  static const std::array<std::vector<std::pair<double, double>>, 100> SCATTERING_FUNCTIONS;

public:
  /**
   * @brief Get incoherent scattering function S(q,Z)
   * @param Z Atomic number
   * @param momentum_transfer q = (sin(θ/2))/λ (Å⁻¹)
   * @return Scattering function value
   */
  static double getIncoherentScatteringFunction(int Z, double momentum_transfer);
  
  /**
   * @brief Get total incoherent (Compton) cross-section
   * Includes binding effects through scattering function
   * @param Z Atomic number
   * @param energy Photon energy (keV)
   * @return Cross-section per atom (barns)
   */
  static double getTotalIncoherentCrossSection(int Z, double energy_keV);
  
  /**
   * @brief Get Klein-Nishina cross-section per electron
   * @param energy Photon energy (keV)
   * @return Cross-section per electron (barns)
   */
  static double getKleinNishinaCrossSection(double energy_keV);
};

/**
 * @brief Complete photon interaction database
 * Combines all interaction types with proper energy dependence
 */
class PhotonInteractionDatabase
{
public:
  /**
   * @brief Get total photon interaction cross-section
   * @param Z Atomic number
   * @param energy Photon energy (keV)
   * @return Total cross-section per atom (barns)
   */
  static double getTotalCrossSection(int Z, double energy_keV);
  
  /**
   * @brief Get individual cross-section components
   * @param Z Atomic number
   * @param energy Photon energy (keV)
   * @return Array of [photoelectric, coherent, incoherent, pair_nuclear, pair_electron]
   */
  static std::array<double, 5> getAllCrossSections(int Z, double energy_keV);
  
  /**
   * @brief Get mass attenuation coefficient (cm²/g)
   * @param Z Atomic number
   * @param atomic_mass Atomic mass (u)
   * @param energy Photon energy (keV)
   * @return Mass attenuation coefficient
   */
  static double getMassAttenuationCoefficient(int Z, double atomic_mass, double energy_keV);
  
  /**
   * @brief Convert between energy units
   */
  static double keVToMeV(double energy_keV) { return energy_keV / 1000.0; }
  static double MeVTokeV(double energy_MeV) { return energy_MeV * 1000.0; }
  
  /**
   * @brief Convert cross-section units
   */
  static double barnsToCm2(double barns) { return barns * 1e-24; }
  static double cm2ToBarns(double cm2) { return cm2 * 1e24; }
};

/**
 * @brief Physical constants for photon interactions
 */
namespace Constants
{
  constexpr double ELECTRON_REST_MASS_KEV = 510.99895000;    // keV
  constexpr double ELECTRON_REST_MASS_MEV = 0.51099895000;   // MeV
  constexpr double CLASSICAL_ELECTRON_RADIUS_CM = 2.8179403262e-13; // cm
  constexpr double FINE_STRUCTURE_CONSTANT = 7.2973525693e-3;
  constexpr double AVOGADRO_NUMBER = 6.02214076e23;          // mol⁻¹
  constexpr double HBAR_C_KEV_FEMTOMETER = 197.3269804;     // keV⋅fm
  constexpr double PAIR_PRODUCTION_THRESHOLD_KEV = 1021.9979; // keV (2⋅mₑc²)
  constexpr double PAIR_PRODUCTION_THRESHOLD_MEV = 1.0219979; // MeV
  constexpr double THOMSON_CROSS_SECTION_BARNS = 0.6652458734; // barns
}

} // namespace PhotonCrossSections