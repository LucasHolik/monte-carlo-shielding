#pragma once

#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <variant>
#include <vector>

/**
 * @brief Type-safe material property storage using C++17 std::variant
 */
using MaterialProperty = std::variant<double, int, std::string>;

/**
 * @brief Element composition in a material
 */
struct ElementComposition
{
  int atomic_number;      // Z (number of protons)
  double atomic_mass;     // A (atomic mass in u)
  double weight_fraction; // Weight fraction in material (0.0 to 1.0)
  std::string symbol;     // Chemical symbol (e.g., "H", "O", "Pb")

  ElementComposition(int z, double a, double fraction, std::string_view sym);
  ElementComposition() = default;
};

/**
 * @brief Material class for Monte Carlo radiation transport simulation
 *
 * Modern C++17 implementation supporting both pure elements and compound
 * materials with type-safe property storage and comprehensive nuclear physics
 * calculations.
 */
class Material
{
private:
  std::string name_;
  double density_;                              // g/cm³
  std::vector<ElementComposition> composition_; // Elements in material
  std::unordered_map<std::string, MaterialProperty>
      properties_;                            // Additional properties
  mutable std::optional<double> effective_z_; // Cached effective Z
  mutable std::optional<double> effective_a_; // Cached effective A

public:
  // Constructors
  Material();
  Material(std::string_view name, double density);
  Material(const Material &other) = default;
  Material &operator=(const Material &other) = default;
  Material(Material &&other) noexcept = default;
  Material &operator=(Material &&other) noexcept = default;

  // Basic accessors
  std::string_view name() const { return name_; }
  double density() const { return density_; }
  const std::vector<ElementComposition> &composition() const
  {
    return composition_;
  }
  size_t getNumberOfElements() const { return composition_.size(); }
  bool isEmpty() const { return composition_.empty(); }

  // Basic mutators
  void setName(std::string_view name) { name_ = name; }
  void setDensity(double density);

  // Element composition management
  void addElement(int atomic_number, double atomic_mass, double weight_fraction,
                  std::string_view symbol);
  void addElement(const ElementComposition &element);
  void clearComposition();
  void normaliseWeightFractions();

  // Property management with C++17 template features
  template <typename T> void setProperty(std::string_view name, const T &value)
  {
    if constexpr(std::is_same_v<T, double> || std::is_same_v<T, int>)
    {
      properties_[std::string(name)] = value;
    }
    else if constexpr(std::is_same_v<T, std::string>)
    {
      properties_[std::string(name)] = value;
    }
    else if constexpr(std::is_convertible_v<T, std::string>)
    {
      properties_[std::string(name)] = std::string(value);
    }
    // Constexpr if ensures only valid types compile
  }

  template <typename T>
  std::optional<T> getProperty(std::string_view name) const
  {
    auto it = properties_.find(std::string(name));
    if(it != properties_.end())
    {
      if constexpr(std::is_same_v<T, double>)
      {
        if(std::holds_alternative<double>(it->second))
        {
          return std::get<double>(it->second);
        }
      }
      else if constexpr(std::is_same_v<T, int>)
      {
        if(std::holds_alternative<int>(it->second))
        {
          return std::get<int>(it->second);
        }
      }
      else if constexpr(std::is_same_v<T, std::string>)
      {
        if(std::holds_alternative<std::string>(it->second))
        {
          return std::get<std::string>(it->second);
        }
      }
    }
    return std::nullopt;
  }

  bool hasProperty(std::string_view name) const;
  void removeProperty(std::string_view name);
  std::vector<std::string> getPropertyNames() const;

  // Nuclear physics calculations
  double getEffectiveAtomicNumber() const;
  double getEffectiveAtomicMass() const;
  double getElectronDensity() const; // electrons/cm³
  double getAtomDensity() const;     // atoms/cm³
  double
  getNumberDensity(int atomic_number) const; // atoms/cm³ for specific element

  // Advanced material properties
  double getRadiationLength() const;          // g/cm² (approximate)
  double getNuclearInteractionLength() const; // g/cm² (approximate)
  std::optional<double>
  getMeanExcitationEnergy() const; // I-value for energy loss calculations

  // Validation
  bool isValid() const;
  bool isNormalised(double tolerance = 1e-6) const;
  std::vector<std::string> validateComposition() const;

  // Utility methods
  bool isPureElement() const { return composition_.size() == 1; }
  bool isCompound() const { return composition_.size() > 1; }

  std::optional<ElementComposition> getElement(int atomic_number) const;
  std::optional<ElementComposition> getElement(std::string_view symbol) const;
  double getWeightFraction(int atomic_number) const;
  double getWeightFraction(std::string_view symbol) const;

  // String representation
  std::string toString() const;
  std::string getCompositionString() const;
  std::string getPropertiesString() const;

  // Comparison operators
  bool operator==(const Material &other) const;
  bool operator!=(const Material &other) const;

  // Static factory methods for common materials
  static Material createVacuum();
  static Material createAir(double density = 0.001225); // g/cm³ at STP
  static Material createWater(double density = 1.0);
  static Material createLead(double density = 11.34);
  static Material createConcrete(double density = 2.3);
  static Material createSteel(double density = 7.87);
  static Material createAluminium(double density = 2.70);
  static Material createPolyethylene(double density = 0.94);

  // Pure element factory
  static Material createElement(int atomic_number, double density,
                                std::string_view name = "",
                                std::string_view symbol = "");

  // Compound material factory
  static Material
  createCompound(std::string_view name, double density,
                 const std::vector<ElementComposition> &elements);

  // Physical constants access
  static double getStandardAtomicMass(int atomic_number);
  static std::string getElementSymbol(int atomic_number);
  static std::string getElementName(int atomic_number);
  static std::optional<int> getAtomicNumber(std::string_view symbol);

private:
  void invalidateCache() const;
  double calculateEffectiveZ() const;
  double calculateEffectiveA() const;
  void validateAndThrow() const;
};

// Utility functions for material calculations
double
calculateMeanAtomicMass(const std::vector<ElementComposition> &composition);
double
calculateMeanAtomicNumber(const std::vector<ElementComposition> &composition);
std::string
formatChemicalFormula(const std::vector<ElementComposition> &composition);