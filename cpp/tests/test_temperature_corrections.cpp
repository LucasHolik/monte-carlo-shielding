#include "Material.hpp"

#include <iomanip>
#include <iostream>
#include <vector>

/**
 * @brief Test temperature-dependent density corrections
 */
int main()
{
  std::cout << "=== Testing Temperature-Dependent Material Properties ===\n";
  std::cout << std::fixed << std::setprecision(6);

  // Test with Lead material
  Material lead = Material::createLead(); // Default density 11.34 g/cm³
  double energy_kev = 100.0;

  std::cout << "\nLead Material Properties:\n";
  std::cout << "  Name: " << lead.name() << "\n";
  std::cout << "  Reference density: " << lead.density() << " g/cm³\n";
  std::cout << "  Reference temperature: " << lead.getReferenceTemperature()
            << " K\n";
  std::cout << "  Thermal expansion coeff: "
            << lead.getThermalExpansionCoefficient() * 1e6 << " × 10⁻⁶ K⁻¹\n";

  // Reference attenuation at room temperature
  double ref_linear_atten = lead.getLinearAttenuationCoefficient(energy_kev);
  std::cout << "\nReference linear attenuation @ " << energy_kev << " keV:\n";
  std::cout << "  " << ref_linear_atten << " cm⁻¹\n";

  // Test temperature range
  std::vector<double> temperatures = {223.15, 273.15, 293.15,
                                      373.15, 473.15, 573.15}; // K
  std::vector<std::string> temp_labels = {
      "Cold (-50°C)", "Freezing (0°C)",   "Room (20°C)",
      "Hot (100°C)",  "Very Hot (200°C)", "Extreme (300°C)"};

  std::cout << "\nTemperature Effects on Lead:\n";
  std::cout
      << "Temperature | Density (g/cm³) | Linear μ (cm⁻¹) | Relative Change\n";
  std::cout << "------------|------------------|------------------|------------"
               "-----\n";

  for(size_t i = 0; i < temperatures.size(); ++i)
  {
    double temp_K = temperatures[i];
    double corrected_density = lead.getTemperatureCorrectedDensity(temp_K);
    double corrected_atten =
        lead.getLinearAttenuationCoefficientAtTemperature(energy_kev, temp_K);
    double relative_change =
        (corrected_atten - ref_linear_atten) / ref_linear_atten * 100.0;

    std::cout << std::left << std::setw(12) << temp_labels[i] << "| "
              << std::setw(16) << corrected_density << " | " << std::setw(16)
              << corrected_atten << " | " << std::setw(15);

    if(std::abs(relative_change) < 0.01)
    {
      std::cout << "0.00%";
    }
    else
    {
      std::cout << std::showpos << relative_change << "%";
    }
    std::cout << std::noshowpos << "\n";
  }

  // Test custom material with high thermal expansion
  std::cout << "\n\nCustom High-Expansion Material Test:\n";
  Material aluminum = Material::createAluminium(); // Default density 2.70 g/cm³
  aluminum.setThermalExpansionCoefficient(69.0e-6); // Aluminum: 69.0 × 10⁻⁶ K⁻¹

  double room_temp = 293.15; // K
  double hot_temp = 573.15;  // 300°C

  double room_density = aluminum.getTemperatureCorrectedDensity(room_temp);
  double hot_density = aluminum.getTemperatureCorrectedDensity(hot_temp);

  double room_atten = aluminum.getLinearAttenuationCoefficientAtTemperature(
      energy_kev, room_temp);
  double hot_atten = aluminum.getLinearAttenuationCoefficientAtTemperature(
      energy_kev, hot_temp);

  std::cout << "  Aluminum @ 20°C:  ρ = " << room_density
            << " g/cm³, μ = " << room_atten << " cm⁻¹\n";
  std::cout << "  Aluminum @ 300°C: ρ = " << hot_density
            << " g/cm³, μ = " << hot_atten << " cm⁻¹\n";
  std::cout << "  Density change: "
            << (hot_density - room_density) / room_density * 100.0 << "%\n";
  std::cout << "  Attenuation change: "
            << (hot_atten - room_atten) / room_atten * 100.0 << "%\n";

  std::cout
      << "\n✅ Temperature-dependent material properties test completed!\n";

  return 0;
}