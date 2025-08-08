#pragma once

#include "Box.hpp"
#include "Material.hpp"
#include "Vector3D.hpp"

#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <variant>
#include <vector>

// Forward declarations for future geometry types
// class Cylinder;
// class Sphere;

/**
 * @brief Type-safe geometry variant using C++17 std::variant
 *
 * Currently supports Box geometry, easily extensible for future shapes.
 * This approach provides compile-time type safety and efficient dispatch.
 */
using GeometryShape =
    std::variant<Box>; // Future: std::variant<Box, Cylinder, Sphere>

/**
 * @brief Result of geometry intersection query
 */
struct GeometryIntersection
{
  bool intersects = false;
  double distance = 0.0;
  Vector3D intersection_point;
  Vector3D surface_normal;
  Material material = Material::createVacuum();
  std::size_t shape_index = 0;
  std::string shape_name = "";
  bool entering = true;

  GeometryIntersection() = default;

  GeometryIntersection(double dist, const Vector3D &point,
                       const Vector3D &normal, const Material &mat,
                       std::size_t index, std::string_view name, bool enter)
      : intersects(true), distance(dist), intersection_point(point),
        surface_normal(normal), material(mat), shape_index(index),
        shape_name(name), entering(enter)
  {}
};

/**
 * @brief Result of point location query
 */
struct LocationResult
{
  bool inside_geometry = false;
  Material material = Material::createVacuum();
  std::size_t shape_index = 0;
  std::string shape_name = "";
  Vector3D local_position;

  LocationResult() = default;

  LocationResult(const Material &mat, std::size_t index, std::string_view name,
                 const Vector3D &pos)
      : inside_geometry(true), material(mat), shape_index(index),
        shape_name(name), local_position(pos)
  {}
};

/**
 * @brief Comprehensive geometry container for Monte Carlo radiation transport
 *
 * Modern C++17 implementation providing multi-region geometry support with
 * type-safe shape storage, efficient ray intersection, and comprehensive
 * boundary detection for particle transport simulations.
 */
class Geometry
{
private:
  std::vector<GeometryShape> shapes_;
  std::vector<std::string> shape_names_;
  std::unordered_map<std::string, std::size_t> name_to_index_;
  Material default_material_;
  std::string name_;

  // Performance optimisation: spatial indexing (future enhancement)
  mutable bool spatial_index_valid_ = false;

public:
  // Constructors
  Geometry();
  explicit Geometry(std::string_view name);
  explicit Geometry(const Material &default_material);
  Geometry(std::string_view name, const Material &default_material);

  // Copy and move semantics
  Geometry(const Geometry &other) = default;
  Geometry &operator=(const Geometry &other) = default;
  Geometry(Geometry &&other) noexcept = default;
  Geometry &operator=(Geometry &&other) noexcept = default;

  // Basic accessors
  std::string_view name() const { return name_; }
  const Material &defaultMaterial() const { return default_material_; }
  std::size_t getNumberOfShapes() const { return shapes_.size(); }
  bool isEmpty() const { return shapes_.empty(); }

  // Basic mutators
  void setName(std::string_view name) { name_ = name; }
  void setDefaultMaterial(const Material &material);

  // Shape management
  std::size_t addShape(const GeometryShape &shape, std::string_view name = "");
  std::size_t addBox(const Box &box, std::string_view name = "");

  // Future shape additions:
  // std::size_t addCylinder(const Cylinder& cylinder, std::string_view name =
  // ""); std::size_t addSphere(const Sphere& sphere, std::string_view name =
  // "");

  bool removeShape(std::size_t index);
  bool removeShape(std::string_view name);
  void clearShapes();

  // Shape access
  std::optional<GeometryShape> getShape(std::size_t index) const;
  std::optional<GeometryShape> getShape(std::string_view name) const;
  std::optional<std::size_t> getShapeIndex(std::string_view name) const;
  std::string getShapeName(std::size_t index) const;

  // Template-based shape retrieval with C++17 type safety
  template <typename ShapeType>
  std::optional<ShapeType> getShapeAs(std::size_t index) const
  {
    if(index >= shapes_.size())
    {
      return std::nullopt;
    }

    if(std::holds_alternative<ShapeType>(shapes_[index]))
    {
      return std::get<ShapeType>(shapes_[index]);
    }

    return std::nullopt;
  }

  template <typename ShapeType>
  std::optional<ShapeType> getShapeAs(std::string_view name) const
  {
    auto index = getShapeIndex(name);
    if(!index)
    {
      return std::nullopt;
    }
    return getShapeAs<ShapeType>(*index);
  }

  // Ray intersection methods
  std::optional<GeometryIntersection>
  findClosestIntersection(const Vector3D &ray_origin,
                          const Vector3D &ray_direction) const;

  std::vector<GeometryIntersection>
  findAllIntersections(const Vector3D &ray_origin,
                       const Vector3D &ray_direction) const;

  std::optional<double> distanceToSurface(const Vector3D &position,
                                          const Vector3D &direction) const;

  // Point location methods
  LocationResult locatePoint(const Vector3D &point) const;
  Material getMaterialAtPoint(const Vector3D &point) const;
  bool isInsideAnyShape(const Vector3D &point) const;

  // Boundary detection
  bool isOnBoundary(const Vector3D &point, double tolerance = 1e-10) const;
  Vector3D getSurfaceNormal(const Vector3D &surface_point,
                            double tolerance = 1e-10) const;

  // Particle transport utilities
  std::optional<GeometryIntersection>
  transportParticle(const Vector3D &start_position, const Vector3D &direction,
                    double max_distance = 1e30) const;

  bool willParticleEscape(const Vector3D &position,
                          const Vector3D &direction) const;

  std::vector<GeometryIntersection>
  getParticlePath(const Vector3D &start_position, const Vector3D &direction,
                  double max_distance = 1e30) const;

  // Advanced intersection queries
  std::vector<std::size_t>
  findIntersectingShapes(const Vector3D &ray_origin,
                         const Vector3D &ray_direction) const;

  std::optional<GeometryIntersection>
  findIntersectionWithShape(std::size_t shape_index, const Vector3D &ray_origin,
                            const Vector3D &ray_direction) const;

  // Geometry analysis
  std::vector<Vector3D>
  getBoundingBox() const; // Returns {min_corner, max_corner}
  double getTotalVolume() const;
  std::vector<Material> getAllMaterials() const;

  // Validation
  bool isValid() const;
  std::vector<std::string> validate() const;
  bool hasOverlappingShapes(double tolerance = 1e-10) const;

  // Statistical and sampling methods
  Vector3D generateRandomPoint() const;
  Vector3D generateRandomPointInShape(std::size_t shape_index) const;
  Vector3D generateRandomPointOnSurface() const;

  // Geometric transformations
  void translateAll(const Vector3D &offset);
  void scaleAll(double factor);
  void rotateAll(const Vector3D &axis, double angle); // Future enhancement

  // I/O and string representation
  std::string toString() const;
  std::string getStatistics() const;

  // Comparison operators
  bool operator==(const Geometry &other) const;
  bool operator!=(const Geometry &other) const;

  // Static factory methods for common geometries
  static Geometry createSimpleShield(const Vector3D &shield_min,
                                     const Vector3D &shield_max,
                                     const Material &shield_material,
                                     const Material &surrounding_material);

  static Geometry createLayeredShield(const std::vector<Box> &layers,
                                      const Material &surrounding_material);

  static Geometry createDetectorGeometry(const Vector3D &detector_position,
                                         const Vector3D &detector_size,
                                         const Material &detector_material,
                                         const Material &surrounding_material);

  static Geometry createShieldingExperiment(double shield_thickness_cm,
                                           const Material &shield_material,
                                           double experiment_length = 10.0,
                                           double cross_section = 20.0);

  // Advanced features for future development
  void enableSpatialIndexing(bool enable = true); // Future optimisation
  void rebuildSpatialIndex() const;               // Future optimisation

private:
  void invalidateSpatialIndex() const { spatial_index_valid_ = false; }
  void updateNameIndex();
  std::string generateDefaultName(std::size_t index) const;

  // Helper methods for visitor pattern with std::variant
  GeometryIntersection
  computeShapeIntersection(const GeometryShape &shape, std::size_t index,
                           const Vector3D &ray_origin,
                           const Vector3D &ray_direction) const;

  bool isPointInsideShape(const GeometryShape &shape,
                          const Vector3D &point) const;
  Material getShapeMaterial(const GeometryShape &shape) const;
  Vector3D getShapeBoundingBoxMin(const GeometryShape &shape) const;
  Vector3D getShapeBoundingBoxMax(const GeometryShape &shape) const;
  double getShapeVolume(const GeometryShape &shape) const;
};

// Utility functions for geometry operations
namespace GeometryUtils
{

// Intersection utilities
bool geometriesIntersect(const Geometry &geom1, const Geometry &geom2);
Geometry combineGeometries(const Geometry &geom1, const Geometry &geom2);

// Analysis utilities
double calculateOverlapVolume(const Geometry &geom1, const Geometry &geom2);
std::vector<Vector3D> findIntersectionPoints(const Geometry &geom1,
                                             const Geometry &geom2);

// Validation utilities
bool validateGeometryForTransport(const Geometry &geometry);
std::vector<std::string> checkGeometryIntegrity(const Geometry &geometry);
} // namespace GeometryUtils

// Advanced geometry intersection result
struct DetailedIntersectionResult
{
  std::vector<GeometryIntersection> intersections;
  std::vector<double> distances;
  std::vector<Material> materials;
  double total_distance = 0.0;
  bool particle_escapes = false;

  bool hasIntersections() const { return !intersections.empty(); }
  std::size_t getNumberOfIntersections() const { return intersections.size(); }

  std::optional<GeometryIntersection> getClosestIntersection() const
  {
    if(intersections.empty())
    {
      return std::nullopt;
    }
    return intersections[0];
  }
};