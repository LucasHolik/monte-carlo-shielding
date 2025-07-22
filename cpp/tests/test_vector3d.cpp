#define _USE_MATH_DEFINES // M_PI is used

#include "Particle.hpp"
#include "Vector3D.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

void testConstructors()
{
  std::cout << "Testing constructors..." << std::endl;

  // Default constructor
  Vector3D v1;
  assert(v1.x() == 0.0);
  assert(v1.y() == 0.0);
  assert(v1.z() == 0.0);

  // Parameterised constructor
  Vector3D v2(1.0, 2.0, 3.0);
  assert(v2.x() == 1.0);
  assert(v2.y() == 2.0);
  assert(v2.z() == 3.0);

  // Copy constructor
  Vector3D v3(v2);
  assert(v3.x() == 1.0);
  assert(v3.y() == 2.0);
  assert(v3.z() == 3.0);

  std::cout << "✓ Constructors passed" << std::endl;
}

void testAccessorsAndMutators()
{
  std::cout << "Testing accessors and mutators..." << std::endl;

  Vector3D v;
  v.setX(5.0);
  v.setY(-2.0);
  v.setZ(7.5);

  assert(v.x() == 5.0);
  assert(v.y() == -2.0);
  assert(v.z() == 7.5);

  v.set(1.0, 2.0, 3.0);
  assert(v.x() == 1.0);
  assert(v.y() == 2.0);
  assert(v.z() == 3.0);

  std::cout << "✓ Accessors and mutators passed" << std::endl;
}

void testStructuredBinding()
{
  std::cout << "Testing C++17 structured binding..." << std::endl;

  Vector3D v(1.5, 2.5, 3.5);
  auto [x, y, z] = v.getComponents();

  assert(x == 1.5);
  assert(y == 2.5);
  assert(z == 3.5);

  std::cout << "✓ Structured binding passed" << std::endl;
}

void testArithmeticOperators()
{
  std::cout << "Testing arithmetic operators..." << std::endl;

  Vector3D v1(1.0, 2.0, 3.0);
  Vector3D v2(4.0, 5.0, 6.0);

  // Addition
  Vector3D v3 = v1 + v2;
  assert(v3.x() == 5.0);
  assert(v3.y() == 7.0);
  assert(v3.z() == 9.0);

  // Subtraction
  Vector3D v4 = v2 - v1;
  assert(v4.x() == 3.0);
  assert(v4.y() == 3.0);
  assert(v4.z() == 3.0);

  // Scalar multiplication
  Vector3D v5 = v1 * 2.0;
  assert(v5.x() == 2.0);
  assert(v5.y() == 4.0);
  assert(v5.z() == 6.0);

  // Scalar division
  Vector3D v6 = v5 / 2.0;
  assert(v6.x() == 1.0);
  assert(v6.y() == 2.0);
  assert(v6.z() == 3.0);

  // Non-member scalar multiplication
  Vector3D v7 = 3.0 * v1;
  assert(v7.x() == 3.0);
  assert(v7.y() == 6.0);
  assert(v7.z() == 9.0);

  std::cout << "✓ Arithmetic operators passed" << std::endl;
}

void testCompoundAssignmentOperators()
{
  std::cout << "Testing compound assignment operators..." << std::endl;

  Vector3D v1(1.0, 2.0, 3.0);
  Vector3D v2(2.0, 3.0, 4.0);

  // +=
  v1 += v2;
  assert(v1.x() == 3.0);
  assert(v1.y() == 5.0);
  assert(v1.z() == 7.0);

  // -=
  v1 -= v2;
  assert(v1.x() == 1.0);
  assert(v1.y() == 2.0);
  assert(v1.z() == 3.0);

  // *=
  v1 *= 2.0;
  assert(v1.x() == 2.0);
  assert(v1.y() == 4.0);
  assert(v1.z() == 6.0);

  // /=
  v1 /= 2.0;
  assert(v1.x() == 1.0);
  assert(v1.y() == 2.0);
  assert(v1.z() == 3.0);

  std::cout << "✓ Compound assignment operators passed" << std::endl;
}

void testUnaryOperators()
{
  std::cout << "Testing unary operators..." << std::endl;

  Vector3D v1(1.0, -2.0, 3.0);
  Vector3D v2 = -v1;

  assert(v2.x() == -1.0);
  assert(v2.y() == 2.0);
  assert(v2.z() == -3.0);

  std::cout << "✓ Unary operators passed" << std::endl;
}

void testComparisonOperators()
{
  std::cout << "Testing comparison operators..." << std::endl;

  Vector3D v1(1.0, 2.0, 3.0);
  Vector3D v2(1.0, 2.0, 3.0);
  Vector3D v3(1.0, 2.0, 3.1);

  assert(v1 == v2);
  assert(v1 != v3);
  assert(!(v1 == v3));
  assert(!(v1 != v2));

  std::cout << "✓ Comparison operators passed" << std::endl;
}

void testVectorCalculations()
{
  std::cout << "Testing vector calculations..." << std::endl;

  Vector3D v1(1.0, 0.0, 0.0);
  Vector3D v2(0.0, 1.0, 0.0);

  // Dot product
  double dot = v1.dot(v2);
  assert(std::abs(dot - 0.0) < 1e-10);

  Vector3D v3(1.0, 1.0, 0.0);
  double dot2 = v1.dot(v3);
  assert(std::abs(dot2 - 1.0) < 1e-10);

  // Cross product
  Vector3D cross = v1.cross(v2);
  assert(std::abs(cross.x() - 0.0) < 1e-10);
  assert(std::abs(cross.y() - 0.0) < 1e-10);
  assert(std::abs(cross.z() - 1.0) < 1e-10);

  // Magnitude
  Vector3D v4(3.0, 4.0, 0.0);
  assert(std::abs(v4.magnitude() - 5.0) < 1e-10);
  assert(std::abs(v4.magnitudeSquared() - 25.0) < 1e-10);

  // Distance
  Vector3D v5(0.0, 0.0, 0.0);
  Vector3D v6(3.0, 4.0, 0.0);
  assert(std::abs(v5.distance(v6) - 5.0) < 1e-10);
  assert(std::abs(v5.distanceSquared(v6) - 25.0) < 1e-10);

  std::cout << "✓ Vector calculations passed" << std::endl;
}

void testNormalisation()
{
  std::cout << "Testing normalisation..." << std::endl;

  // Normal case
  Vector3D v1(3.0, 4.0, 0.0);
  auto normalised = v1.tryNormalise();
  assert(normalised.has_value());
  assert(std::abs(normalised->x() - 0.6) < 1e-10);
  assert(std::abs(normalised->y() - 0.8) < 1e-10);
  assert(std::abs(normalised->z() - 0.0) < 1e-10);
  assert(normalised->isUnit());

  // Zero vector case
  Vector3D v2(0.0, 0.0, 0.0);
  auto zeroNormalised = v2.tryNormalise();
  assert(!zeroNormalised.has_value());

  // Test normalise() method
  Vector3D v3(1.0, 1.0, 1.0);
  Vector3D normalised2 = v3.normalise();
  assert(normalised2.isUnit());

  // Test that normalise() throws for zero vector
  bool threw = false;
  try
  {
    v2.normalise();
  }
  catch(const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  std::cout << "✓ Normalisation passed" << std::endl;
}

void testStaticMethods()
{
  std::cout << "Testing static methods..." << std::endl;

  auto unitVector = Vector3D::createUnitVector(3.0, 4.0, 0.0);
  assert(unitVector.has_value());
  assert(unitVector->isUnit());

  auto zeroUnit = Vector3D::createUnitVector(0.0, 0.0, 0.0);
  assert(!zeroUnit.has_value());

  std::cout << "✓ Static methods passed" << std::endl;
}

void testUtilityMethods()
{
  std::cout << "Testing utility methods..." << std::endl;

  Vector3D v1(0.0, 0.0, 0.0);
  assert(v1.isZero());

  Vector3D v2(1e-12, 1e-12, 1e-12);
  assert(v2.isZero());

  Vector3D v3(1.0, 0.0, 0.0);
  assert(v3.isUnit());

  Vector3D v4(1.0000001, 0.0, 0.0);
  assert(v4.isUnit(1e-6));

  std::cout << "✓ Utility methods passed" << std::endl;
}

void testStaticConstants()
{
  std::cout << "Testing static constants..." << std::endl;

  assert(Vector3D::ZERO.isZero());
  assert(Vector3D::UNIT_X.isUnit());
  assert(Vector3D::UNIT_Y.isUnit());
  assert(Vector3D::UNIT_Z.isUnit());

  assert(Vector3D::UNIT_X.x() == 1.0);
  assert(Vector3D::UNIT_Y.y() == 1.0);
  assert(Vector3D::UNIT_Z.z() == 1.0);

  std::cout << "✓ Static constants passed" << std::endl;
}

void testErrorHandling()
{
  std::cout << "Testing error handling..." << std::endl;

  Vector3D v1(1.0, 2.0, 3.0);

  // Division by zero
  bool threw = false;
  try
  {
    v1 / 0.0;
  }
  catch(const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  threw = false;
  try
  {
    v1 /= 0.0;
  }
  catch(const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  std::cout << "✓ Error handling passed" << std::endl;
}

void testUtilityFunctions()
{
  std::cout << "Testing utility functions..." << std::endl;

  Vector3D v1(1.0, 0.0, 0.0);
  Vector3D v2(0.0, 1.0, 0.0);

  double angle = angleBetween(v1, v2);
  assert(std::abs(angle - M_PI / 2.0) < 1e-10);

  auto safeAngle = safeAngleBetween(v1, v2);
  assert(safeAngle.has_value());
  assert(std::abs(*safeAngle - M_PI / 2.0) < 1e-10);

  Vector3D zero(0.0, 0.0, 0.0);
  auto zeroAngle = safeAngleBetween(v1, zero);
  assert(!zeroAngle.has_value());

  bool threw = false;
  try
  {
    angleBetween(v1, zero);
  }
  catch(const std::invalid_argument &)
  {
    threw = true;
  }
  assert(threw);

  std::cout << "✓ Utility functions passed" << std::endl;
}

int main()
{
  std::cout << "Running Vector3D comprehensive tests...\n" << std::endl;

  testConstructors();
  testAccessorsAndMutators();
  testStructuredBinding();
  testArithmeticOperators();
  testCompoundAssignmentOperators();
  testUnaryOperators();
  testComparisonOperators();
  testVectorCalculations();
  testNormalisation();
  testStaticMethods();
  testUtilityMethods();
  testStaticConstants();
  testErrorHandling();
  testUtilityFunctions();

  std::cout << "\n🎉 All Vector3D tests passed successfully!" << std::endl;
  std::cout << "Your Vector3D class is working correctly and ready for commit."
            << std::endl;

  return 0;
}