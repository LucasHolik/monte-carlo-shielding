/*
 * Monte Carlo Radiation Transport Simulator - Hello World Test
 *
 * Modern C++17 pybind11 test to verify C++/Python integration.
 * Demonstrates C++17 features and professional coding practices.
 *
 * Author: Lucas Holik
 * Date: January 2025
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <string>
#include <vector>
#include <cmath>
#include <optional>
#include <string_view>
#include <array>

namespace py = pybind11;

/**
 * Modern greeting function using string_view for efficiency
 */
std::string say_hello(std::string_view name = "Monte Carlo Simulator")
{
  return std::string("Hello, ") + std::string(name) + "! C++17 to Python binding is working perfectly.";
}

/**
 * Modern energy calculation with C++17 features
 */
std::optional<double> calculate_exponential_attenuation(double thickness, double linear_attenuation_coefficient)
{
  // Input validation using modern error handling
  if (thickness < 0.0 || linear_attenuation_coefficient < 0.0)
  {
    return std::nullopt; // C++17 std::optional for error handling
  }

  // Simple Beer-Lambert law: I = I₀ * exp(-μt)
  // Returns the attenuation factor (I/I₀)
  return std::exp(-linear_attenuation_coefficient * thickness);
}

/**
 * C++17 structured bindings demonstration
 */
std::pair<double, std::string> analyse_attenuation(double attenuation_factor)
{
  std::string category;

  // C++17 if constexpr could be used here for compile-time decisions
  if (attenuation_factor > 0.9)
  {
    category = "Low attenuation";
  }
  else if (attenuation_factor > 0.1)
  {
    category = "Moderate attenuation";
  }
  else
  {
    category = "High attenuation";
  }

  return {attenuation_factor, category};
}

/**
 * Modern vector generation using C++17 features
 */
std::vector<double> generate_test_energies(int count, double max_energy = 1.0)
{
  std::vector<double> energies;
  energies.reserve(count);

  // C++17 range-based loop with structured binding (if we had pairs)
  for (int i = 0; i < count; ++i)
  {
    // Modern initialization
    auto energy = max_energy * static_cast<double>(i) / static_cast<double>(count - 1);
    energies.push_back(energy);
  }

  return energies;
}

/**
 * Simple 3D vector class prototype (foundation for future Vector3D)
 */
class SimpleVector3D
{
public:
  double x, y, z;

  SimpleVector3D(double x = 0.0, double y = 0.0, double z = 0.0)
      : x(x), y(y), z(z) {}

  // Basic vector operations
  SimpleVector3D operator+(const SimpleVector3D &other) const
  {
    return SimpleVector3D(x + other.x, y + other.y, z + other.z);
  }

  SimpleVector3D operator*(double scalar) const
  {
    return SimpleVector3D(x * scalar, y * scalar, z * scalar);
  }

  double magnitude() const
  {
    return std::sqrt(x * x + y * y + z * z);
  }

  double dot(const SimpleVector3D &other) const
  {
    return x * other.x + y * other.y + z * other.z;
  }

  std::string to_string() const
  {
    return "Vector3D(" + std::to_string(x) + ", " +
           std::to_string(y) + ", " + std::to_string(z) + ")";
  }
};

/**
 * Test class for future particle implementation
 */
class TestParticle
{
private:
  SimpleVector3D position_;
  SimpleVector3D direction_;
  double energy_;
  bool alive_;

public:
  TestParticle(const SimpleVector3D &pos, const SimpleVector3D &dir, double energy)
      : position_(pos), direction_(dir), energy_(energy), alive_(true) {}

  // Getters
  SimpleVector3D get_position() const { return position_; }
  SimpleVector3D get_direction() const { return direction_; }
  double get_energy() const { return energy_; }
  bool is_alive() const { return alive_; }

  // Setters
  void set_position(const SimpleVector3D &pos) { position_ = pos; }
  void set_energy(double energy) { energy_ = energy; }
  void kill() { alive_ = false; }

  // Simple movement function
  void move(double distance)
  {
    if (alive_)
    {
      position_ = position_ + direction_ * distance;
    }
  }

  std::string info() const
  {
    return "Particle at " + position_.to_string() +
           " with energy " + std::to_string(energy_) + " MeV" +
           (alive_ ? " (alive)" : " (terminated)");
  }
};

/**
 * Pybind11 module definition
 * This exposes all our C++ functions and classes to Python
 */
PYBIND11_MODULE(mcshield, m)
{
  m.doc() = "Monte Carlo Radiation Transport Simulator - Hello World Test Module";

  // Module information
  m.attr("__version__") = "0.1.0";
  m.attr("__author__") = "Lucas Holik";

  // Basic functions
  m.def("say_hello", &say_hello,
        "A simple greeting function to test C++/Python integration",
        py::arg("name") = "Monte Carlo Simulator");

  m.def("calculate_exponential_attenuation", &calculate_exponential_attenuation,
        "Calculate radiation attenuation using Beer-Lambert law (returns None for invalid inputs)",
        py::arg("thickness"), py::arg("linear_attenuation_coefficient"));

  m.def("analyse_attenuation", &analyse_attenuation,
        "Analyse attenuation factor and categorise the result",
        py::arg("attenuation_factor"));

  m.def("generate_test_energies", &generate_test_energies,
        "Generate a vector of test energy values",
        py::arg("count"), py::arg("max_energy") = 1.0);

  // SimpleVector3D class binding
  py::class_<SimpleVector3D>(m, "SimpleVector3D")
      .def(py::init<double, double, double>(),
           py::arg("x") = 0.0, py::arg("y") = 0.0, py::arg("z") = 0.0)
      .def_readwrite("x", &SimpleVector3D::x)
      .def_readwrite("y", &SimpleVector3D::y)
      .def_readwrite("z", &SimpleVector3D::z)
      .def("__add__", &SimpleVector3D::operator+)
      .def("__mul__", &SimpleVector3D::operator*)
      .def("magnitude", &SimpleVector3D::magnitude)
      .def("dot", &SimpleVector3D::dot)
      .def("__str__", &SimpleVector3D::to_string)
      .def("__repr__", &SimpleVector3D::to_string);

  // TestParticle class binding
  py::class_<TestParticle>(m, "TestParticle")
      .def(py::init<const SimpleVector3D &, const SimpleVector3D &, double>())
      .def("get_position", &TestParticle::get_position)
      .def("get_direction", &TestParticle::get_direction)
      .def("get_energy", &TestParticle::get_energy)
      .def("is_alive", &TestParticle::is_alive)
      .def("set_position", &TestParticle::set_position)
      .def("set_energy", &TestParticle::set_energy)
      .def("kill", &TestParticle::kill)
      .def("move", &TestParticle::move)
      .def("info", &TestParticle::info)
      .def("__str__", &TestParticle::info);
}