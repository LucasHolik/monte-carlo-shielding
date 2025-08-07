#include "Material.hpp"
#include "PhotonCrossSectionDatabase.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <array>

// Physical constants
namespace PhysicalConstants
{
constexpr double AVOGADRO_NUMBER = 6.02214076e23;              // mol⁻¹
constexpr double CLASSICAL_ELECTRON_RADIUS = 2.8179403262e-13; // cm
constexpr double ELECTRON_REST_MASS = 0.510999;                // MeV
constexpr double FINE_STRUCTURE_CONSTANT = 1.0 / 137.036;      // dimensionless
} // namespace PhysicalConstants

// Standard atomic masses (u) - first 30 elements
static const double STANDARD_ATOMIC_MASSES[] = {
    0.0,    // Z=0 (placeholder)
    1.008,  // H
    4.003,  // He
    6.939,  // Li
    9.012,  // Be
    10.81,  // B
    12.01,  // C
    14.01,  // N
    16.00,  // O
    19.00,  // F
    20.18,  // Ne
    22.99,  // Na
    24.31,  // Mg
    26.98,  // Al
    28.09,  // Si
    30.97,  // P
    32.07,  // S
    35.45,  // Cl
    39.95,  // Ar
    39.10,  // K
    40.08,  // Ca
    44.96,  // Sc
    47.87,  // Ti
    50.94,  // V
    52.00,  // Cr
    54.94,  // Mn
    55.85,  // Fe
    58.93,  // Co
    58.69,  // Ni
    63.55,  // Cu
    65.38,  // Zn
    69.72,  // Ga
    72.64,  // Ge
    74.92,  // As
    78.96,  // Se
    79.90,  // Br
    83.80,  // Kr
    85.47,  // Rb
    87.62,  // Sr
    88.91,  // Y
    91.22,  // Zr
    92.91,  // Nb
    95.96,  // Mo
    98.0,   // Tc
    101.07, // Ru
    102.91, // Rh
    106.42, // Pd
    107.87, // Ag
    112.41, // Cd
    114.82, // In
    118.71, // Sn
    121.76, // Sb
    127.60, // Te
    126.90, // I
    131.29, // Xe
    132.91, // Cs
    137.33, // Ba
    138.91, // La
    140.12, // Ce
    140.91, // Pr
    144.24, // Nd
    145.0,  // Pm
    150.36, // Sm
    151.96, // Eu
    157.25, // Gd
    158.93, // Tb
    162.50, // Dy
    164.93, // Ho
    167.26, // Er
    168.93, // Tm
    173.05, // Yb
    174.97, // Lu
    178.49, // Hf
    180.95, // Ta
    183.84, // W
    186.21, // Re
    190.23, // Os
    192.22, // Ir
    195.08, // Pt
    196.97, // Au
    200.59, // Hg
    204.38, // Tl
    207.2,  // Pb
    208.98, // Bi
    209.0,  // Po
    210.0,  // At
    222.0,  // Rn
    223.0,  // Fr
    226.0,  // Ra
    227.0,  // Ac
    232.04, // Th
    231.04, // Pa
    238.03  // U
};

// Element symbols
static const std::string ELEMENT_SYMBOLS[] = {
    "", // Z=0
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg",
    "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr",
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U"};

// Element names
static const std::string ELEMENT_NAMES[] = {
    "", // Z=0
    "Hydrogen",   "Helium",     "Lithium",      "Beryllium",
    "Boron",      "Carbon",     "Nitrogen",     "Oxygen",
    "Fluorine",   "Neon",       "Sodium",       "Magnesium",
    "Aluminium",  "Silicon",    "Phosphorus",   "Sulphur",
    "Chlorine",   "Argon",      "Potassium",    "Calcium",
    "Scandium",   "Titanium",   "Vanadium",     "Chromium",
    "Manganese",  "Iron",       "Cobalt",       "Nickel",
    "Copper",     "Zinc",       "Gallium",      "Germanium",
    "Arsenic",    "Selenium",   "Bromine",      "Krypton",
    "Rubidium",   "Strontium",  "Yttrium",      "Zirconium",
    "Niobium",    "Molybdenum", "Technetium",   "Ruthenium",
    "Rhodium",    "Palladium",  "Silver",       "Cadmium",
    "Indium",     "Tin",        "Antimony",     "Tellurium",
    "Iodine",     "Xenon",      "Caesium",      "Barium",
    "Lanthanum",  "Cerium",     "Praseodymium", "Neodymium",
    "Promethium", "Samarium",   "Europium",     "Gadolinium",
    "Terbium",    "Dysprosium", "Holmium",      "Erbium",
    "Thulium",    "Ytterbium",  "Lutetium",     "Hafnium",
    "Tantalum",   "Tungsten",   "Rhenium",      "Osmium",
    "Iridium",    "Platinum",   "Gold",         "Mercury",
    "Thallium",   "Lead",       "Bismuth",      "Polonium",
    "Astatine",   "Radon",      "Francium",     "Radium",
    "Actinium",   "Thorium",    "Protactinium", "Uranium"};

constexpr size_t MAX_ATOMIC_NUMBER = 92;

// ElementComposition implementation
ElementComposition::ElementComposition(int z, double a, double fraction,
                                       std::string_view sym)
    : atomic_number(z), atomic_mass(a), weight_fraction(fraction), symbol(sym)
{
  if(z < 0 || z > static_cast<int>(MAX_ATOMIC_NUMBER))
  {
    throw std::invalid_argument("Invalid atomic number: " + std::to_string(z));
  }
  if(a <= 0.0)
  {
    throw std::invalid_argument("Atomic mass must be positive");
  }
  if(fraction < 0.0 || fraction > 1.0)
  {
    throw std::invalid_argument("Weight fraction must be between 0.0 and 1.0");
  }
}

// Material implementation  
Material::Material() : name_("Unknown"), density_(0.0), reference_temperature_(293.15), thermal_expansion_coeff_(0.0) {}

Material::Material(std::string_view name, double density)
    : name_(name), density_(density), reference_temperature_(293.15), thermal_expansion_coeff_(0.0)
{
  if(density < 0.0)
  {
    throw std::invalid_argument("Material density cannot be negative");
  }
}

void Material::setDensity(double density)
{
  if(density < 0.0)
  {
    throw std::invalid_argument("Material density cannot be negative");
  }
  density_ = density;
  invalidateCache();
}

void Material::addElement(int atomic_number, double atomic_mass,
                          double weight_fraction, std::string_view symbol)
{
  composition_.emplace_back(atomic_number, atomic_mass, weight_fraction,
                            symbol);
  invalidateCache();
}

void Material::addElement(const ElementComposition &element)
{
  composition_.push_back(element);
  invalidateCache();
}

void Material::clearComposition()
{
  composition_.clear();
  invalidateCache();
}

void Material::normaliseWeightFractions()
{
  if(composition_.empty())
    return;

  double total_weight =
      std::accumulate(composition_.begin(), composition_.end(), 0.0,
                      [](double sum, const ElementComposition &elem) {
                        return sum + elem.weight_fraction;
                      });

  if(total_weight > 0.0)
  {
    for(auto &element : composition_)
    {
      element.weight_fraction /= total_weight;
    }
    invalidateCache();
  }
}

bool Material::hasProperty(std::string_view name) const
{
  return properties_.find(std::string(name)) != properties_.end();
}

void Material::removeProperty(std::string_view name)
{
  properties_.erase(std::string(name));
}

std::vector<std::string> Material::getPropertyNames() const
{
  std::vector<std::string> names;
  names.reserve(properties_.size());
  for(const auto &[name, value] : properties_)
  {
    names.push_back(name);
  }
  return names;
}

double Material::getEffectiveAtomicNumber() const
{
  if(!effective_z_.has_value())
  {
    effective_z_ = calculateEffectiveZ();
  }
  return *effective_z_;
}

double Material::getEffectiveAtomicMass() const
{
  if(!effective_a_.has_value())
  {
    effective_a_ = calculateEffectiveA();
  }
  return *effective_a_;
}

double Material::getElectronDensity() const
{
  if(composition_.empty())
    return 0.0;

  double electron_density = 0.0;
  for(const auto &element : composition_)
  {
    double number_density = density_ * PhysicalConstants::AVOGADRO_NUMBER *
                            element.weight_fraction / element.atomic_mass;
    electron_density += number_density * element.atomic_number;
  }
  return electron_density;
}

double Material::getAtomDensity() const
{
  if(composition_.empty())
    return 0.0;

  double atom_density = 0.0;
  for(const auto &element : composition_)
  {
    atom_density += density_ * PhysicalConstants::AVOGADRO_NUMBER *
                    element.weight_fraction / element.atomic_mass;
  }
  return atom_density;
}

double Material::getNumberDensity(int atomic_number) const
{
  for(const auto &element : composition_)
  {
    if(element.atomic_number == atomic_number)
    {
      return density_ * PhysicalConstants::AVOGADRO_NUMBER *
             element.weight_fraction / element.atomic_mass;
    }
  }
  return 0.0;
}

double Material::getRadiationLength() const
{
  if(composition_.empty())
    return std::numeric_limits<double>::infinity();

  // Approximate formula for radiation length (PDG approximation)
  double z_eff = getEffectiveAtomicNumber();
  double a_eff = getEffectiveAtomicMass();

  if(z_eff < 1.0 || a_eff < 1.0)
    return std::numeric_limits<double>::infinity();

  // X₀ ≈ 716.4 A / (Z(Z+1)ln(287/√Z)) g/cm²
  double log_term = std::log(287.0 / std::sqrt(z_eff));
  return 716.4 * a_eff / (z_eff * (z_eff + 1.0) * log_term);
}

double Material::getNuclearInteractionLength() const
{
  if(composition_.empty())
    return std::numeric_limits<double>::infinity();

  // Approximate nuclear interaction length
  double a_eff = getEffectiveAtomicMass();
  return a_eff * 1.4; // Very rough approximation
}

std::optional<double> Material::getMeanExcitationEnergy() const
{
  if(composition_.empty())
    return std::nullopt;

  // Rough approximation: I ≈ 10 * Z eV for Z > 1, 19.2 eV for hydrogen
  double z_eff = getEffectiveAtomicNumber();
  if(z_eff < 1.0)
    return std::nullopt;

  if(z_eff < 1.5)
  {
    return 19.2e-6; // 19.2 eV in MeV for hydrogen
  }
  else
  {
    return 10.0 * z_eff * 1e-6; // 10*Z eV in MeV
  }
}

bool Material::isValid() const
{
  if(density_ < 0.0)
    return false;
  if(composition_.empty())
    return true; // Vacuum is valid

  for(const auto &element : composition_)
  {
    if(element.atomic_number < 0 ||
       element.atomic_number > static_cast<int>(MAX_ATOMIC_NUMBER))
      return false;
    if(element.atomic_mass <= 0.0)
      return false;
    if(element.weight_fraction < 0.0 || element.weight_fraction > 1.0)
      return false;
  }

  return isNormalised();
}

bool Material::isNormalised(double tolerance) const
{
  if(composition_.empty())
    return true;

  double total_weight =
      std::accumulate(composition_.begin(), composition_.end(), 0.0,
                      [](double sum, const ElementComposition &elem) {
                        return sum + elem.weight_fraction;
                      });

  return std::abs(total_weight - 1.0) <= tolerance;
}

std::vector<std::string> Material::validateComposition() const
{
  std::vector<std::string> errors;

  if(density_ < 0.0)
  {
    errors.push_back("Negative density");
  }

  if(composition_.empty() && density_ > 0.0)
  {
    errors.push_back("Non-zero density but empty composition");
  }

  for(size_t i = 0; i < composition_.size(); ++i)
  {
    const auto &element = composition_[i];
    std::string prefix = "Element " + std::to_string(i) + ": ";

    if(element.atomic_number < 1 ||
       element.atomic_number > static_cast<int>(MAX_ATOMIC_NUMBER))
    {
      errors.push_back(prefix + "Invalid atomic number " +
                       std::to_string(element.atomic_number));
    }
    if(element.atomic_mass <= 0.0)
    {
      errors.push_back(prefix + "Non-positive atomic mass");
    }
    if(element.weight_fraction < 0.0)
    {
      errors.push_back(prefix + "Negative weight fraction");
    }
    if(element.weight_fraction > 1.0)
    {
      errors.push_back(prefix + "Weight fraction > 1.0");
    }
  }

  if(!isNormalised(1e-6))
  {
    errors.push_back("Weight fractions do not sum to 1.0");
  }

  return errors;
}

std::optional<ElementComposition> Material::getElement(int atomic_number) const
{
  auto it = std::find_if(composition_.begin(), composition_.end(),
                         [atomic_number](const ElementComposition &elem) {
                           return elem.atomic_number == atomic_number;
                         });

  return (it != composition_.end()) ? std::optional<ElementComposition>(*it)
                                    : std::nullopt;
}

std::optional<ElementComposition>
Material::getElement(std::string_view symbol) const
{
  auto it = std::find_if(composition_.begin(), composition_.end(),
                         [symbol](const ElementComposition &elem) {
                           return elem.symbol == symbol;
                         });

  return (it != composition_.end()) ? std::optional<ElementComposition>(*it)
                                    : std::nullopt;
}

double Material::getWeightFraction(int atomic_number) const
{
  auto element = getElement(atomic_number);
  return element ? element->weight_fraction : 0.0;
}

double Material::getWeightFraction(std::string_view symbol) const
{
  auto element = getElement(symbol);
  return element ? element->weight_fraction : 0.0;
}

std::string Material::toString() const
{
  std::stringstream ss;
  ss << "Material '" << name_ << "' (density: " << density_ << " g/cm³)\n";
  ss << getCompositionString();

  if(!properties_.empty())
  {
    ss << "\nProperties:\n" << getPropertiesString();
  }

  return ss.str();
}

std::string Material::getCompositionString() const
{
  if(composition_.empty())
  {
    return "Empty composition (vacuum)";
  }

  std::stringstream ss;
  ss << "Composition:\n";

  for(size_t i = 0; i < composition_.size(); ++i)
  {
    const auto &element = composition_[i];
    ss << "  " << element.symbol << " (Z=" << element.atomic_number
       << ", A=" << std::fixed << std::setprecision(3) << element.atomic_mass
       << "): " << std::setprecision(1) << element.weight_fraction * 100.0
       << "%\n";
  }

  return ss.str();
}

std::string Material::getPropertiesString() const
{
  if(properties_.empty())
  {
    return "No additional properties";
  }

  std::stringstream ss;
  for(const auto &[name, value] : properties_)
  {
    ss << "  " << name << ": ";

    std::visit([&ss](const auto &v) { ss << v; }, value);

    ss << "\n";
  }

  return ss.str();
}

bool Material::operator==(const Material &other) const
{
  constexpr double tolerance = 1e-10;

  if(name_ != other.name_ || std::abs(density_ - other.density_) > tolerance)
  {
    return false;
  }

  if(composition_.size() != other.composition_.size())
  {
    return false;
  }

  for(size_t i = 0; i < composition_.size(); ++i)
  {
    const auto &elem1 = composition_[i];
    const auto &elem2 = other.composition_[i];

    if(elem1.atomic_number != elem2.atomic_number ||
       std::abs(elem1.atomic_mass - elem2.atomic_mass) > tolerance ||
       std::abs(elem1.weight_fraction - elem2.weight_fraction) > tolerance ||
       elem1.symbol != elem2.symbol)
    {
      return false;
    }
  }

  return properties_ == other.properties_;
}

bool Material::operator!=(const Material &other) const
{
  return !(*this == other);
}

// Static factory methods
Material Material::createVacuum()
{
  Material vacuum("Vacuum", 0.0);
  return vacuum;
}

Material Material::createAir(double density)
{
  Material air("Air", density);
  air.addElement(7, 14.01, 0.755, "N");   // Nitrogen
  air.addElement(8, 16.00, 0.232, "O");   // Oxygen
  air.addElement(18, 39.95, 0.013, "Ar"); // Argon
  return air;
}

Material Material::createWater(double density)
{
  Material water("Water", density);
  water.addElement(1, 1.008, 0.112, "H"); // Hydrogen
  water.addElement(8, 16.00, 0.888, "O"); // Oxygen
  return water;
}

Material Material::createLead(double density)
{
  Material lead = createElement(82, density, "Lead", "Pb");
  lead.setThermalExpansionCoefficient(87.4e-6); // Lead: 87.4 × 10⁻⁶ K⁻¹
  return lead;
}

Material Material::createConcrete(double density)
{
  Material concrete("Concrete", density);
  // Typical concrete composition
  concrete.addElement(1, 1.008, 0.022, "H");   // Hydrogen
  concrete.addElement(6, 12.01, 0.002, "C");   // Carbon
  concrete.addElement(8, 16.00, 0.575, "O");   // Oxygen
  concrete.addElement(11, 22.99, 0.015, "Na"); // Sodium
  concrete.addElement(12, 24.31, 0.001, "Mg"); // Magnesium
  concrete.addElement(13, 26.98, 0.019, "Al"); // Aluminium
  concrete.addElement(14, 28.09, 0.304, "Si"); // Silicon
  concrete.addElement(19, 39.10, 0.010, "K");  // Potassium
  concrete.addElement(20, 40.08, 0.050, "Ca"); // Calcium
  concrete.addElement(26, 55.85, 0.002, "Fe"); // Iron
  return concrete;
}

Material Material::createSteel(double density)
{
  Material steel("Steel", density);
  steel.addElement(26, 55.85, 0.98, "Fe"); // Iron
  steel.addElement(6, 12.01, 0.02, "C");   // Carbon
  return steel;
}

Material Material::createAluminium(double density)
{
  return createElement(13, density, "Aluminium", "Al");
}

Material Material::createPolyethylene(double density)
{
  Material polyethylene("Polyethylene", density);
  polyethylene.addElement(1, 1.008, 0.143, "H"); // Hydrogen
  polyethylene.addElement(6, 12.01, 0.857, "C"); // Carbon (C₂H₄)ₙ
  return polyethylene;
}

Material Material::createElement(int atomic_number, double density,
                                 std::string_view name, std::string_view symbol)
{
  if(atomic_number < 1 || atomic_number > static_cast<int>(MAX_ATOMIC_NUMBER))
  {
    throw std::invalid_argument("Invalid atomic number: " +
                                std::to_string(atomic_number));
  }

  std::string element_name =
      name.empty() ? getElementName(atomic_number) : std::string(name);
  std::string element_symbol =
      symbol.empty() ? getElementSymbol(atomic_number) : std::string(symbol);
  double atomic_mass = getStandardAtomicMass(atomic_number);

  Material material(element_name, density);
  material.addElement(atomic_number, atomic_mass, 1.0, element_symbol);

  return material;
}

Material
Material::createCompound(std::string_view name, double density,
                         const std::vector<ElementComposition> &elements)
{
  Material material(name, density);

  for(const auto &element : elements)
  {
    material.addElement(element);
  }

  return material;
}

// Static utility methods
double Material::getStandardAtomicMass(int atomic_number)
{
  if(atomic_number < 0 || atomic_number > static_cast<int>(MAX_ATOMIC_NUMBER))
  {
    throw std::invalid_argument("Invalid atomic number: " +
                                std::to_string(atomic_number));
  }
  return STANDARD_ATOMIC_MASSES[atomic_number];
}

std::string Material::getElementSymbol(int atomic_number)
{
  if(atomic_number < 0 || atomic_number > static_cast<int>(MAX_ATOMIC_NUMBER))
  {
    return "?";
  }
  return ELEMENT_SYMBOLS[atomic_number];
}

std::string Material::getElementName(int atomic_number)
{
  if(atomic_number < 0 || atomic_number > static_cast<int>(MAX_ATOMIC_NUMBER))
  {
    return "Unknown";
  }
  return ELEMENT_NAMES[atomic_number];
}

std::optional<int> Material::getAtomicNumber(std::string_view symbol)
{
  for(size_t i = 1; i <= MAX_ATOMIC_NUMBER; ++i)
  {
    if(ELEMENT_SYMBOLS[i] == symbol)
    {
      return static_cast<int>(i);
    }
  }
  return std::nullopt;
}

// Private methods
void Material::invalidateCache() const
{
  effective_z_.reset();
  effective_a_.reset();
}

double Material::calculateEffectiveZ() const
{
  if(composition_.empty())
    return 0.0;

  // Calculate effective Z using electron fraction weighting
  double total_electrons = 0.0;
  double weighted_z = 0.0;

  for(const auto &element : composition_)
  {
    double electrons_per_gram =
        element.weight_fraction * element.atomic_number / element.atomic_mass;
    total_electrons += electrons_per_gram;
    weighted_z += electrons_per_gram * element.atomic_number;
  }

  return (total_electrons > 0.0) ? weighted_z / total_electrons : 0.0;
}

double Material::calculateEffectiveA() const
{
  if(composition_.empty())
    return 0.0;

  // Calculate effective A using weight fraction
  return std::accumulate(composition_.begin(), composition_.end(), 0.0,
                         [](double sum, const ElementComposition &elem) {
                           return sum + elem.weight_fraction * elem.atomic_mass;
                         });
}

void Material::validateAndThrow() const
{
  auto errors = validateComposition();
  if(!errors.empty())
  {
    std::stringstream ss;
    ss << "Material validation failed:\n";
    for(const auto &error : errors)
    {
      ss << "  - " << error << "\n";
    }
    throw std::invalid_argument(ss.str());
  }
}

// Utility functions
double
calculateMeanAtomicMass(const std::vector<ElementComposition> &composition)
{
  return std::accumulate(composition.begin(), composition.end(), 0.0,
                         [](double sum, const ElementComposition &elem) {
                           return sum + elem.weight_fraction * elem.atomic_mass;
                         });
}

double
calculateMeanAtomicNumber(const std::vector<ElementComposition> &composition)
{
  double total_atoms = 0.0;
  double weighted_z = 0.0;

  for(const auto &element : composition)
  {
    double atom_fraction = element.weight_fraction / element.atomic_mass;
    total_atoms += atom_fraction;
    weighted_z += atom_fraction * element.atomic_number;
  }

  return (total_atoms > 0.0) ? weighted_z / total_atoms : 0.0;
}

std::string
formatChemicalFormula(const std::vector<ElementComposition> &composition)
{
  if(composition.empty())
  {
    return "Vacuum";
  }

  if(composition.size() == 1)
  {
    return composition[0].symbol;
  }

  std::stringstream ss;
  for(size_t i = 0; i < composition.size(); ++i)
  {
    if(i > 0)
      ss << " + ";
    ss << composition[i].symbol << "(" << std::fixed << std::setprecision(1)
       << composition[i].weight_fraction * 100.0 << "%)";
  }

  return ss.str();
}

// =============================================================================
// PHOTON INTERACTION CROSS-SECTIONS
// =============================================================================

double Material::getTotalPhotonCrossSection(double energy_keV) const
{
  if(isEmpty() || energy_keV <= 0.0)
  {
    return 0.0;
  }

  double weighted_cross_section = 0.0;
  
  // Sum cross-sections weighted by atom fractions
  for(const auto& element : composition_)
  {
    double atom_fraction = getAtomFraction(element.atomic_number);
    double element_xs = PhotonCrossSections::PhotonInteractionDatabase::getTotalCrossSection(
        element.atomic_number, energy_keV);
    weighted_cross_section += atom_fraction * element_xs;
  }

  return weighted_cross_section;
}

double Material::getMassAttenuationCoefficient(double energy_keV) const
{
  if(isEmpty() || energy_keV <= 0.0)
  {
    return 0.0;
  }

  double mass_attenuation = 0.0;

  // Sum mass attenuation coefficients weighted by mass fractions
  for(const auto& element : composition_)
  {
    double mass_atten_element = PhotonCrossSections::PhotonInteractionDatabase::getMassAttenuationCoefficient(
        element.atomic_number, element.atomic_mass, energy_keV);
    mass_attenuation += element.weight_fraction * mass_atten_element;
  }

  return mass_attenuation;
}

double Material::getLinearAttenuationCoefficient(double energy_keV) const
{
  double mass_atten = getMassAttenuationCoefficient(energy_keV);
  return mass_atten * density_; // μ = (μ/ρ) × ρ
}

std::array<double, 5> Material::getPhotonCrossSectionComponents(double energy_keV) const
{
  std::array<double, 5> total_components = {0.0, 0.0, 0.0, 0.0, 0.0};
  
  if(isEmpty() || energy_keV <= 0.0)
  {
    return total_components;
  }

  // Sum components weighted by atom fractions
  for(const auto& element : composition_)
  {
    double atom_fraction = getAtomFraction(element.atomic_number);
    auto element_components = PhotonCrossSections::PhotonInteractionDatabase::getAllCrossSections(
        element.atomic_number, energy_keV);
    
    for(size_t i = 0; i < 5; ++i)
    {
      total_components[i] += atom_fraction * element_components[i];
    }
  }

  return total_components;
}

// Helper method to calculate atom fraction from weight fraction
double Material::getAtomFraction(int atomic_number) const
{
  auto element = std::find_if(composition_.begin(), composition_.end(),
      [atomic_number](const ElementComposition& elem) {
        return elem.atomic_number == atomic_number;
      });
  
  if(element == composition_.end())
  {
    return 0.0;
  }

  // Calculate total moles
  double total_moles = 0.0;
  for(const auto& elem : composition_)
  {
    total_moles += elem.weight_fraction / elem.atomic_mass;
  }

  if(total_moles == 0.0)
  {
    return 0.0;
  }

  // Calculate atom fraction for this element
  double element_moles = element->weight_fraction / element->atomic_mass;
  return element_moles / total_moles;
}

// =============================================================================
// TEMPERATURE-DEPENDENT DENSITY CORRECTIONS
// =============================================================================

double Material::getTemperatureCorrectedDensity(double temperature_K) const
{
  if(temperature_K <= 0.0)
  {
    return density_;
  }

  // ρ(T) = ρ₀ / [1 + α(T - T₀)]
  // where α is volumetric thermal expansion coefficient
  double temperature_diff = temperature_K - reference_temperature_;
  double volume_expansion = 1.0 + thermal_expansion_coeff_ * temperature_diff;
  
  if(volume_expansion <= 0.0)
  {
    // Prevent negative/zero densities
    return density_;
  }

  return density_ / volume_expansion;
}

double Material::getLinearAttenuationCoefficientAtTemperature(double energy_keV, double temperature_K) const
{
  double mass_atten = getMassAttenuationCoefficient(energy_keV);
  double temperature_corrected_density = getTemperatureCorrectedDensity(temperature_K);
  return mass_atten * temperature_corrected_density;
}

void Material::setThermalExpansionCoefficient(double alpha_per_K)
{
  thermal_expansion_coeff_ = alpha_per_K;
}

double Material::getThermalExpansionCoefficient() const
{
  return thermal_expansion_coeff_;
}

void Material::setReferenceTemperature(double temperature_K)
{
  reference_temperature_ = temperature_K;
}

double Material::getReferenceTemperature() const
{
  return reference_temperature_;
}