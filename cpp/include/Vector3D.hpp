#pragma once

#include <optional>
#include <tuple>

/**
 * @brief 3D vector class for Monte Carlo radiation transport simulation
 *
 * Provides comprehensive 3D vector operations with C++17 features including
 * std::optional for safe operations and structured binding support.
 */
class Vector3D
{
private:
  double x_, y_, z_;

public:
  // Constructors
  Vector3D();
  Vector3D(double x, double y, double z);
  Vector3D(const Vector3D &other) = default;
  Vector3D &operator=(const Vector3D &other) = default;

  // Accessors
  double x() const { return x_; }
  double y() const { return y_; }
  double z() const { return z_; }

  // Mutators
  void setX(double x) { x_ = x; }
  void setY(double y) { y_ = y; }
  void setZ(double z) { z_ = z; }
  void set(double x, double y, double z);

  // C++17 structured binding support
  auto getComponents() const { return std::make_tuple(x_, y_, z_); }

  // Vector operations
  Vector3D operator+(const Vector3D &other) const;
  Vector3D operator-(const Vector3D &other) const;
  Vector3D operator*(double scalar) const;
  Vector3D operator/(double scalar) const;
  Vector3D &operator+=(const Vector3D &other);
  Vector3D &operator-=(const Vector3D &other);
  Vector3D &operator*=(double scalar);
  Vector3D &operator/=(double scalar);

  // Unary operators
  Vector3D operator-() const;

  // Comparison operators
  bool operator==(const Vector3D &other) const;
  bool operator!=(const Vector3D &other) const;

  // Vector calculations
  double dot(const Vector3D &other) const;
  Vector3D cross(const Vector3D &other) const;
  double magnitude() const;
  double magnitudeSquared() const;
  double distance(const Vector3D &other) const;
  double distanceSquared(const Vector3D &other) const;

  // Safe normalisation with C++17 std::optional
  std::optional<Vector3D> tryNormalise() const;
  Vector3D normalise() const; // Throws if magnitude is zero

  // Unit vector creation
  static std::optional<Vector3D> createUnitVector(double x, double y, double z);

  // Utility methods
  bool isZero(double tolerance = 1e-10) const;
  bool isUnit(double tolerance = 1e-10) const;

  // Static utility vectors
  static const Vector3D ZERO;
  static const Vector3D UNIT_X;
  static const Vector3D UNIT_Y;
  static const Vector3D UNIT_Z;
};

// Non-member operators
Vector3D operator*(double scalar, const Vector3D &vector);

// Utility functions
double angleBetween(const Vector3D &a, const Vector3D &b);
std::optional<double> safeAngleBetween(const Vector3D &a, const Vector3D &b);