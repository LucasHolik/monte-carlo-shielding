#pragma once

#include <array>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

/**
 * @brief NIST XCOM Cross-Section Database for Nuclear Industry Applications
 * 
 * Direct implementation using authoritative NIST XCOM photon cross-section data.
 * Provides nuclear industry accuracy standards with proper edge handling.
 * 
 * Data format: Energy (MeV), Coherent, Incoherent, Photoelectric, Pair Nuclear, Pair Electron, Total w/ Coherent, Total w/o Coherent (all in cm²/g)
 */

namespace PhotonCrossSections
{

/**
 * @brief Cross-section interaction types
 */
enum class CrossSectionType {
    COHERENT = 0,           // Coherent (Rayleigh) scattering
    INCOHERENT = 1,         // Incoherent (Compton) scattering
    PHOTOELECTRIC = 2,      // Photoelectric absorption
    PAIR_NUCLEAR = 3,       // Pair production in nuclear field
    PAIR_ELECTRON = 4,      // Pair production in electron field
    TOTAL_WITH_COHERENT = 5,    // Total cross-section including coherent
    TOTAL_WITHOUT_COHERENT = 6  // Total cross-section excluding coherent
};

/**
 * @brief Physical constants for cross-section calculations
 */
namespace Constants {
    constexpr double AVOGADRO_NUMBER = 6.02214076e23;        // mol⁻¹
    constexpr double BARNS_TO_CM2 = 1e-24;                   // cm²/barn
    constexpr double CM2_TO_BARNS = 1e24;                    // barn/cm²
    constexpr double MEV_TO_KEV = 1000.0;                    // keV/MeV
    constexpr double KEV_TO_MEV = 1e-3;                      // MeV/keV
}

/**
 * @brief NIST Cross-Section Database Manager
 * 
 * Thread-safe singleton that loads and manages NIST XCOM data for all elements.
 * Uses lazy loading and caching for optimal performance.
 */
class NISTCrossSectionDatabase {
private:
    static std::unique_ptr<NISTCrossSectionDatabase> instance_;
    static std::mutex mutex_;
    
    // Data structure: particle_type -> atomic_number -> [energy, cross_section_data...]
    // Each data point is array<double,8>: [energy_MeV, coherent, incoherent, photoelectric, pair_nuclear, pair_electron, total_w_coherent, total_wo_coherent]
    std::map<std::string, std::map<int, std::vector<std::array<double, 8>>>> cross_section_data_;
    
    NISTCrossSectionDatabase() = default;
    
public:
    /**
     * @brief Get singleton instance (thread-safe)
     */
    static NISTCrossSectionDatabase& getInstance() {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!instance_) {
            instance_ = std::unique_ptr<NISTCrossSectionDatabase>(new NISTCrossSectionDatabase());
        }
        return *instance_;
    }
    
    /**
     * @brief Get cross-section for specified interaction type
     * @param Z Atomic number (1-100)
     * @param particle_type Particle type ("photon", extensible to "neutron", "electron")
     * @param energy_keV Photon energy in keV
     * @param type Cross-section interaction type
     * @return Cross-section in cm²/g
     */
    double getCrossSection(int Z, const std::string& particle_type, double energy_keV, CrossSectionType type);
    
    /**
     * @brief Get all cross-section components
     * @param Z Atomic number (1-100)
     * @param particle_type Particle type ("photon")
     * @param energy_keV Photon energy in keV
     * @return Array of all cross-sections: [coherent, incoherent, photoelectric, pair_nuclear, pair_electron, total_w_coherent, total_wo_coherent]
     */
    std::array<double, 7> getAllCrossSections(int Z, const std::string& particle_type, double energy_keV);
    
    /**
     * @brief Convert cross-section to mass attenuation coefficient
     * @param cross_section_cm2_per_g Cross-section in cm²/g
     * @return Mass attenuation coefficient in cm²/g (same units, for consistency)
     */
    double getMassAttenuationCoefficient(double cross_section_cm2_per_g) const {
        return cross_section_cm2_per_g;  // Already in correct units
    }
    
    /**
     * @brief Convert cross-section units cm²/g to barns/atom
     * @param cross_section_cm2_per_g Cross-section in cm²/g
     * @param atomic_mass_u Atomic mass in atomic mass units
     * @return Cross-section in barns/atom
     */
    double convertToBarnsPerAtom(double cross_section_cm2_per_g, double atomic_mass_u) const;
    
    /**
     * @brief Check if element data is loaded
     * @param Z Atomic number
     * @param particle_type Particle type
     * @return True if data is available
     */
    bool isDataLoaded(int Z, const std::string& particle_type) const;
    
    /**
     * @brief Get energy range for loaded data
     * @param Z Atomic number
     * @param particle_type Particle type
     * @return Pair of [min_energy_keV, max_energy_keV]
     */
    std::pair<double, double> getEnergyRange(int Z, const std::string& particle_type) const;

private:
    /**
     * @brief Load element data from NIST XCOM file
     * @param Z Atomic number
     * @param particle_type Particle type
     * @return True if successful
     */
    bool loadElementData(int Z, const std::string& particle_type);
    
    /**
     * @brief Parse NIST XCOM data file
     * @param filename Path to NIST data file
     * @return Vector of data points [energy_MeV, cross_sections...]
     */
    std::vector<std::array<double, 8>> parseNISTFile(const std::string& filename);
    
    /**
     * @brief Interpolate cross-section using log-log interpolation
     * @param data Vector of data points
     * @param energy_MeV Target energy in MeV
     * @param cross_section_index Index of cross-section type (0-6)
     * @return Interpolated cross-section in cm²/g
     */
    double interpolateCrossSection(const std::vector<std::array<double, 8>>& data, 
                                 double energy_MeV, 
                                 size_t cross_section_index) const;
    
    /**
     * @brief Find bracketing indices for interpolation
     * @param data Vector of data points
     * @param energy_MeV Target energy in MeV
     * @return Pair of [lower_index, upper_index]
     */
    std::pair<size_t, size_t> findBracketingIndices(const std::vector<std::array<double, 8>>& data, 
                                                   double energy_MeV) const;
    
    /**
     * @brief Handle absorption edge cases (duplicate energies with different cross-sections)
     * @param data Vector of data points
     * @param energy_MeV Exact energy match
     * @param cross_section_index Cross-section type index
     * @return Cross-section value (uses post-edge value for exact matches)
     */
    double handleEdgeCase(const std::vector<std::array<double, 8>>& data,
                         double energy_MeV,
                         size_t cross_section_index) const;
    
    /**
     * @brief Get data file path for element
     * @param Z Atomic number
     * @param particle_type Particle type
     * @return Full path to data file
     */
    std::string getDataFilePath(int Z, const std::string& particle_type) const;
    
    /**
     * @brief Validate atomic number range
     * @param Z Atomic number
     * @return True if valid (1-100)
     */
    bool isValidAtomicNumber(int Z) const {
        return Z >= 1 && Z <= 100;
    }
    
    /**
     * @brief Validate energy range
     * @param energy_keV Energy in keV
     * @return True if positive
     */
    bool isValidEnergy(double energy_keV) const {
        return energy_keV > 0.0;
    }
};

/**
 * @brief Convenient wrapper functions for common operations
 */
namespace NIST {
    
    /**
     * @brief Get photon cross-section (convenience function)
     * @param Z Atomic number
     * @param energy_keV Photon energy in keV
     * @param type Cross-section type
     * @return Cross-section in cm²/g
     */
    inline double getPhotonCrossSection(int Z, double energy_keV, CrossSectionType type) {
        return NISTCrossSectionDatabase::getInstance().getCrossSection(Z, "photon", energy_keV, type);
    }
    
    /**
     * @brief Get total photon attenuation coefficient
     * @param Z Atomic number
     * @param energy_keV Photon energy in keV
     * @param include_coherent Whether to include coherent scattering
     * @return Total cross-section in cm²/g
     */
    inline double getTotalPhotonCrossSection(int Z, double energy_keV, bool include_coherent = true) {
        CrossSectionType type = include_coherent ? CrossSectionType::TOTAL_WITH_COHERENT : CrossSectionType::TOTAL_WITHOUT_COHERENT;
        return getPhotonCrossSection(Z, energy_keV, type);
    }
    
    /**
     * @brief Get all photon cross-sections
     * @param Z Atomic number
     * @param energy_keV Photon energy in keV
     * @return Array of [coherent, incoherent, photoelectric, pair_nuclear, pair_electron, total_w_coherent, total_wo_coherent]
     */
    inline std::array<double, 7> getAllPhotonCrossSections(int Z, double energy_keV) {
        return NISTCrossSectionDatabase::getInstance().getAllCrossSections(Z, "photon", energy_keV);
    }
}

} // namespace PhotonCrossSections