#include "Vector3D.hpp"
#include <stdexcept>
#include <cmath>
#include <algorithm>

// Static constants
const Vector3D Vector3D::ZERO(0.0, 0.0, 0.0);
const Vector3D Vector3D::UNIT_X(1.0, 0.0, 0.0);
const Vector3D Vector3D::UNIT_Y(0.0, 1.0, 0.0);
const Vector3D Vector3D::UNIT_Z(0.0, 0.0, 1.0);

// Constructors
Vector3D::Vector3D() : x_(0.0), y_(0.0), z_(0.0) {}

Vector3D::Vector3D(double x, double y, double z) : x_(x), y_(y), z_(z) {}

void Vector3D::set(double x, double y, double z)
{
  x_ = x;
  y_ = y;
  z_ = z;
}

// Addition operators
Vector3D Vector3D::operator+(const Vector3D &other) const
{
  return Vector3D(x_ + other.x_, y_ + other.y_, z_ + other.z_);
}

Vector3D &Vector3D::operator+=(const Vector3D &other)
{
  x_ += other.x_;
  y_ += other.y_;
  z_ += other.z_;
  return *this;
}

// Subtraction operators
Vector3D Vector3D::operator-(const Vector3D &other) const
{
  return Vector3D(x_ - other.x_, y_ - other.y_, z_ - other.z_);
}

Vector3D &Vector3D::operator-=(const Vector3D &other)
{
  x_ -= other.x_;
  y_ -= other.y_;
  z_ -= other.z_;
  return *this;
}

// Scalar multiplication operators
Vector3D Vector3D::operator*(double scalar) const
{
  return Vector3D(x_ * scalar, y_ * scalar, z_ * scalar);
}

Vector3D &Vector3D::operator*=(double scalar)
{
  x_ *= scalar;
  y_ *= scalar;
  z_ *= scalar;
  return *this;
}

// Scalar division operators
Vector3D Vector3D::operator/(double scalar) const
{
  if (std::abs(scalar) < 1e-15)
  {
    throw std::invalid_argument("Division by zero in Vector3D");
  }
  return Vector3D(x_ / scalar, y_ / scalar, z_ / scalar);
}

Vector3D &Vector3D::operator/=(double scalar)
{
  if (std::abs(scalar) < 1e-15)
  {
    throw std::invalid_argument("Division by zero in Vector3D");
  }
  x_ /= scalar;
  y_ /= scalar;
  z_ /= scalar;
  return *this;
}

// Unary minus
Vector3D Vector3D::operator-() const
{
  return Vector3D(-x_, -y_, -z_);
}

// Comparison operators
bool Vector3D::operator==(const Vector3D &other) const
{
  constexpr double tolerance = 1e-10;
  return std::abs(x_ - other.x_) < tolerance &&
         std::abs(y_ - other.y_) < tolerance &&
         std::abs(z_ - other.z_) < tolerance;
}

bool Vector3D::operator!=(const Vector3D &other) const
{
  return !(*this == other);
}

// Dot product
double Vector3D::dot(const Vector3D &other) const
{
  return x_ * other.x_ + y_ * other.y_ + z_ * other.z_;
}

// Cross product
Vector3D Vector3D::cross(const Vector3D &other) const
{
  return Vector3D(
      y_ * other.z_ - z_ * other.y_,
      z_ * other.x_ - x_ * other.z_,
      x_ * other.y_ - y_ * other.x_);
}

// Magnitude calculations
double Vector3D::magnitude() const
{
  return std::sqrt(x_ * x_ + y_ * y_ + z_ * z_);
}

double Vector3D::magnitudeSquared() const
{
  return x_ * x_ + y_ * y_ + z_ * z_;
}

// Distance calculations
double Vector3D::distance(const Vector3D &other) const
{
  return (*this - other).magnitude();
}

double Vector3D::distanceSquared(const Vector3D &other) const
{
  return (*this - other).magnitudeSquared();
}

// Safe normalisation with std::optional (C++17 feature)
std::optional<Vector3D> Vector3D::tryNormalise() const
{
  double mag = magnitude();
  if (mag < 1e-15)
  {
    return std::nullopt; // Cannot normalise zero vector
  }
  return Vector3D(x_ / mag, y_ / mag, z_ / mag);
}

// Normalisation that throws on zero vector
Vector3D Vector3D::normalise() const
{
  auto result = tryNormalise();
  if (!result)
  {
    throw std::invalid_argument("Cannot normalise zero vector");
  }
  return *result;
}

// Static unit vector creation
std::optional<Vector3D> Vector3D::createUnitVector(double x, double y, double z)
{
  Vector3D temp(x, y, z);
  return temp.tryNormalise();
}

// Utility methods
bool Vector3D::isZero(double tolerance) const
{
  return magnitude() < tolerance;
}

bool Vector3D::isUnit(double tolerance) const
{
  return std::abs(magnitude() - 1.0) < tolerance;
}

// Non-member operators
Vector3D operator*(double scalar, const Vector3D &vector)
{
  return vector * scalar;
}

// Utility functions
double angleBetween(const Vector3D &a, const Vector3D &b)
{
  double magA = a.magnitude();
  double magB = b.magnitude();

  if (magA < 1e-15 || magB < 1e-15)
  {
    throw std::invalid_argument("Cannot calculate angle with zero vector");
  }

  double cosTheta = a.dot(b) / (magA * magB);

  // Clamp to [-1, 1] to handle floating point errors
  cosTheta = std::max(-1.0, std::min(1.0, cosTheta));

  return std::acos(cosTheta);
}

std::optional<double> safeAngleBetween(const Vector3D &a, const Vector3D &b)
{
  double magA = a.magnitude();
  double magB = b.magnitude();

  if (magA < 1e-15 || magB < 1e-15)
  {
    return std::nullopt;
  }

  double cosTheta = a.dot(b) / (magA * magB);

  // Clamp to [-1, 1] to handle floating point errors
  cosTheta = std::max(-1.0, std::min(1.0, cosTheta));

  return std::acos(cosTheta);
}