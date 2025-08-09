/**
 * @file PhotonCrossSectionDatabase.cpp
 * @brief NIST XCOM Cross-Section Database Implementation
 * 
 * Direct implementation using authoritative NIST XCOM photon cross-section data.
 * Replaces all manual physics calculations with NIST data interpolation.
 */

#include "PhotonCrossSectionDatabase.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <filesystem>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace PhotonCrossSections
{

// Static member definitions
std::unique_ptr<NISTCrossSectionDatabase> NISTCrossSectionDatabase::instance_ = nullptr;
std::mutex NISTCrossSectionDatabase::mutex_;

// Element symbol lookup for file names
const std::array<std::string, 101> ELEMENT_SYMBOLS = {
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm"
};

double NISTCrossSectionDatabase::getCrossSection(int Z, const std::string& particle_type, 
                                               double energy_keV, CrossSectionType type) {
    // Input validation
    if (!isValidAtomicNumber(Z)) {
        throw std::invalid_argument("Invalid atomic number: " + std::to_string(Z) + " (must be 1-100)");
    }
    if (!isValidEnergy(energy_keV)) {
        throw std::invalid_argument("Invalid energy: " + std::to_string(energy_keV) + " (must be positive)");
    }
    
    // Load data if not already loaded
    if (!isDataLoaded(Z, particle_type)) {
        if (!loadElementData(Z, particle_type)) {
            throw std::runtime_error("Failed to load data for element Z=" + std::to_string(Z) + 
                                   ", particle=" + particle_type);
        }
    }
    
    // Get data reference
    const auto& data = cross_section_data_[particle_type][Z];
    if (data.empty()) {
        throw std::runtime_error("No data available for Z=" + std::to_string(Z) + 
                               ", particle=" + particle_type);
    }
    
    // Convert energy to MeV for consistency with data
    double energy_MeV = energy_keV * Constants::KEV_TO_MEV;
    
    // Get cross-section index
    size_t index = static_cast<size_t>(type) + 1; // +1 because energy is at index 0
    if (index > 7) {
        throw std::invalid_argument("Invalid cross-section type");
    }
    
    // Interpolate and return
    return interpolateCrossSection(data, energy_MeV, index);
}

std::array<double, 7> NISTCrossSectionDatabase::getAllCrossSections(int Z, const std::string& particle_type, 
                                                                   double energy_keV) {
    // Input validation
    if (!isValidAtomicNumber(Z)) {
        throw std::invalid_argument("Invalid atomic number: " + std::to_string(Z));
    }
    if (!isValidEnergy(energy_keV)) {
        throw std::invalid_argument("Invalid energy: " + std::to_string(energy_keV));
    }
    
    // Load data if needed
    if (!isDataLoaded(Z, particle_type)) {
        if (!loadElementData(Z, particle_type)) {
            throw std::runtime_error("Failed to load data for element Z=" + std::to_string(Z));
        }
    }
    
    const auto& data = cross_section_data_[particle_type][Z];
    double energy_MeV = energy_keV * Constants::KEV_TO_MEV;
    
    std::array<double, 7> result;
    for (size_t i = 0; i < 7; ++i) {
        result[i] = interpolateCrossSection(data, energy_MeV, i + 1); // +1 for energy offset
    }
    
    return result;
}

double NISTCrossSectionDatabase::convertToBarnsPerAtom(double cross_section_cm2_per_g, 
                                                     double atomic_mass_u) const {
    return cross_section_cm2_per_g * Constants::CM2_TO_BARNS * atomic_mass_u / Constants::AVOGADRO_NUMBER;
}

bool NISTCrossSectionDatabase::isDataLoaded(int Z, const std::string& particle_type) const {
    auto particle_it = cross_section_data_.find(particle_type);
    if (particle_it == cross_section_data_.end()) {
        return false;
    }
    
    auto element_it = particle_it->second.find(Z);
    return element_it != particle_it->second.end() && !element_it->second.empty();
}

std::pair<double, double> NISTCrossSectionDatabase::getEnergyRange(int Z, const std::string& particle_type) const {
    if (!isDataLoaded(Z, particle_type)) {
        return {0.0, 0.0};
    }
    
    const auto& data = cross_section_data_.at(particle_type).at(Z);
    if (data.empty()) {
        return {0.0, 0.0};
    }
    
    double min_energy_keV = data.front()[0] * Constants::MEV_TO_KEV;
    double max_energy_keV = data.back()[0] * Constants::MEV_TO_KEV;
    
    return {min_energy_keV, max_energy_keV};
}

bool NISTCrossSectionDatabase::loadElementData(int Z, const std::string& particle_type) {
    std::string filepath = getDataFilePath(Z, particle_type);
    
    try {
        auto data = parseNISTFile(filepath);
        if (data.empty()) {
            std::cerr << "Warning: No data loaded from " << filepath << std::endl;
            return false;
        }
        
        // Store the data
        cross_section_data_[particle_type][Z] = std::move(data);
        
        std::cout << "Loaded " << cross_section_data_[particle_type][Z].size() 
                  << " data points for Z=" << Z << " (" << ELEMENT_SYMBOLS[Z] << "), " 
                  << particle_type << std::endl;
        
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error loading " << filepath << ": " << e.what() << std::endl;
        return false;
    }
}

std::vector<std::array<double, 8>> NISTCrossSectionDatabase::parseNISTFile(const std::string& filename) {
    std::vector<std::array<double, 8>> data;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::string line;
    size_t line_number = 0;
    
    while (std::getline(file, line)) {
        ++line_number;
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue;
        }
        
        // Skip header lines (first few lines typically contain column headers)
        if (line_number <= 3 || line.find("Photon") != std::string::npos || 
            line.find("Energy") != std::string::npos || line.find("Scatter") != std::string::npos) {
            continue;
        }
        
        // Parse data line
        std::istringstream iss(line);
        std::array<double, 8> point;
        
        // Read all 8 columns: energy + 7 cross-sections
        bool parse_success = true;
        for (size_t i = 0; i < 8; ++i) {
            if (!(iss >> point[i])) {
                parse_success = false;
                break;
            }
        }
        
        if (parse_success) {
            // Validate data
            if (point[0] > 0.0) { // Energy must be positive
                data.push_back(point);
            }
        }
    }
    
    file.close();
    
    // Sort by energy (should already be sorted, but ensure it)
    std::sort(data.begin(), data.end(), 
              [](const std::array<double, 8>& a, const std::array<double, 8>& b) {
                  return a[0] < b[0]; // Compare energy (index 0)
              });
    
    return data;
}

double NISTCrossSectionDatabase::interpolateCrossSection(const std::vector<std::array<double, 8>>& data,
                                                       double energy_MeV,
                                                       size_t cross_section_index) const {
    if (data.empty() || cross_section_index >= 8) {
        return 0.0;
    }
    
    // Handle boundary cases
    if (energy_MeV <= data.front()[0]) {
        return data.front()[cross_section_index];
    }
    if (energy_MeV >= data.back()[0]) {
        return data.back()[cross_section_index];
    }
    
    // Find bracketing points
    auto indices = findBracketingIndices(data, energy_MeV);
    size_t lower_idx = indices.first;
    size_t upper_idx = indices.second;
    
    // Check for exact match (absorption edge case)
    if (std::abs(data[lower_idx][0] - energy_MeV) < 1e-12 || 
        std::abs(data[upper_idx][0] - energy_MeV) < 1e-12) {
        return handleEdgeCase(data, energy_MeV, cross_section_index);
    }
    
    // Extract values for interpolation
    double E1 = data[lower_idx][0];
    double E2 = data[upper_idx][0];
    double sigma1 = data[lower_idx][cross_section_index];
    double sigma2 = data[upper_idx][cross_section_index];
    
    // Handle zero or negative cross-sections (use linear interpolation)
    if (sigma1 <= 0.0 || sigma2 <= 0.0) {
        double f = (energy_MeV - E1) / (E2 - E1);
        return sigma1 + f * (sigma2 - sigma1);
    }
    
    // Log-log interpolation for positive values
    double log_E1 = std::log(E1);
    double log_E2 = std::log(E2);
    double log_sigma1 = std::log(sigma1);
    double log_sigma2 = std::log(sigma2);
    double log_E_target = std::log(energy_MeV);
    
    double slope = (log_sigma2 - log_sigma1) / (log_E2 - log_E1);
    double log_sigma_target = log_sigma1 + slope * (log_E_target - log_E1);
    
    return std::exp(log_sigma_target);
}

std::pair<size_t, size_t> NISTCrossSectionDatabase::findBracketingIndices(
    const std::vector<std::array<double, 8>>& data, double energy_MeV) const {
    
    // Binary search for upper bound
    auto upper_it = std::lower_bound(data.begin(), data.end(), energy_MeV,
                                   [](const std::array<double, 8>& point, double energy) {
                                       return point[0] < energy;
                                   });
    
    if (upper_it == data.begin()) {
        return {0, 0};
    }
    if (upper_it == data.end()) {
        size_t last_idx = data.size() - 1;
        return {last_idx, last_idx};
    }
    
    size_t upper_idx = std::distance(data.begin(), upper_it);
    size_t lower_idx = upper_idx - 1;
    
    return {lower_idx, upper_idx};
}

double NISTCrossSectionDatabase::handleEdgeCase(const std::vector<std::array<double, 8>>& data,
                                              double energy_MeV,
                                              size_t cross_section_index) const {
    // Find all points with matching energy (within tolerance)
    constexpr double tolerance = 1e-12;
    
    std::vector<size_t> matching_indices;
    for (size_t i = 0; i < data.size(); ++i) {
        if (std::abs(data[i][0] - energy_MeV) < tolerance) {
            matching_indices.push_back(i);
        }
    }
    
    if (matching_indices.empty()) {
        // No exact match, this shouldn't happen
        return 0.0;
    }
    
    if (matching_indices.size() == 1) {
        // Single point, just return it
        return data[matching_indices[0]][cross_section_index];
    }
    
    // Multiple points at same energy (absorption edge)
    // Use the higher cross-section value (post-edge)
    double max_cross_section = 0.0;
    for (size_t idx : matching_indices) {
        max_cross_section = std::max(max_cross_section, data[idx][cross_section_index]);
    }
    
    return max_cross_section;
}

std::string NISTCrossSectionDatabase::getDataFilePath(int Z, const std::string& particle_type) const {
    if (particle_type == "photon") {
        // NIST XCOM file naming convention: element_XXX_Symbol.txt
        char filename[50];
        std::snprintf(filename, sizeof(filename), "element_%03d_%s.txt", Z, ELEMENT_SYMBOLS[Z].c_str());
        
        // Build relative path from source location
        std::string base_path = "cpp/data/nist_xcom/photon/";
        return base_path + filename;
    }
    
    throw std::invalid_argument("Unsupported particle type: " + particle_type);
}

} // namespace PhotonCrossSections