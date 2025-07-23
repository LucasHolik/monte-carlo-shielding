#pragma once

#include "Material.hpp"
#include "Vector3D.hpp"

#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

/**
 * @brief Surface normal directions for box faces
 */
enum class BoxFace
{
  XMin, // -X face (left)
  XMax, // +X face (right)
  YMin, // -Y face (bottom)
  YMax, // +Y face (top)
  ZMin, // -Z face (back)
  ZMax  // +Z face (front)
};

/**
 * @brief Result of ray-box intersection calculation
 */
struct RayBoxIntersection
{
  bool intersects;             // Whether ray intersects box
  double distance;             // Distance to intersection point
  Vector3D intersection_point; // Point of intersection
  BoxFace face;                // Which face was intersected
  Vector3D surface_normal;     // Outward surface normal at intersection
  bool entering;               // True if entering box, false if exiting

  RayBoxIntersection()
      : intersects(false), distance(0.0), face(BoxFace::XMin), entering(true)
  {}

  RayBoxIntersection(double dist, const Vector3D &point, BoxFace f,
                     const Vector3D &normal, bool enter)
      : intersects(true), distance(dist), intersection_point(point), face(f),
        surface_normal(normal), entering(enter)
  {}
};

/**
 * @brief Axis-aligned box geometry class for Monte Carlo radiation transport
 *
 * Modern C++17 implementation supporting ray-box intersection, boundary
 * detection, and material assignment for rectangular shielding configurations.
 */
class Box
{
private:
  Vector3D min_corner_;
  Vector3D max_corner_;
  Material material_;
  std::string name_;
  bool is_valid_;

  // Cached values for performance
  mutable std::optional<Vector3D> centre_;
  mutable std::optional<Vector3D> dimensions_;
  mutable std::optional<double> volume_;

public:
  // Constructors
  Box();
  Box(const Vector3D &min_corner, const Vector3D &max_corner);
  Box(const Vector3D &min_corner, const Vector3D &max_corner,
      const Material &material);
  Box(const Vector3D &min_corner, const Vector3D &max_corner,
      const Material &material, std::string_view name);
  Box(const Box &other) = default;
  Box &operator=(const Box &other) = default;
  Box(Box &&other) noexcept = default;
  Box &operator=(Box &&other) noexcept = default;

  // Basic accessors
  const Vector3D &minCorner() const { return min_corner_; }
  const Vector3D &maxCorner() const { return max_corner_; }
  const Material &material() const { return material_; }
  std::string_view name() const { return name_; }
  bool isValid() const { return is_valid_; }

  // Computed properties
  Vector3D centre() const;
  Vector3D dimensions() const; // Width, height, depth
  double volume() const;
  double surfaceArea() const;

  // Individual dimension accessors
  double width() const { return max_corner_.x() - min_corner_.x(); }
  double height() const { return max_corner_.y() - min_corner_.y(); }
  double depth() const { return max_corner_.z() - min_corner_.z(); }

  // Basic mutators
  void setMinCorner(const Vector3D &min_corner);
  void setMaxCorner(const Vector3D &max_corner);
  void setCorners(const Vector3D &min_corner, const Vector3D &max_corner);
  void setMaterial(const Material &material) { material_ = material; }
  void setName(std::string_view name) { name_ = name; }

  // Geometry queries
  bool isInside(const Vector3D &point, double tolerance = 0.0) const;
  bool isOnSurface(const Vector3D &point, double tolerance = 1e-10) const;
  bool isOutside(const Vector3D &point, double tolerance = 0.0) const;

  // Distance calculations
  double distanceToPoint(const Vector3D &point) const;
  double distanceToPointSquared(const Vector3D &point) const;

  // C++17 safe ray intersection with std::optional
  std::optional<RayBoxIntersection>
  tryRayIntersection(const Vector3D &ray_origin,
                     const Vector3D &ray_direction) const;

  // Traditional ray intersection (throws on invalid input)
  RayBoxIntersection rayIntersection(const Vector3D &ray_origin,
                                     const Vector3D &ray_direction) const;

  // Distance to surface along ray direction
  std::optional<double> distanceToSurface(const Vector3D &position,
                                          const Vector3D &direction) const;

  // Multiple intersection points (entry and exit)
  std::vector<RayBoxIntersection>
  allRayIntersections(const Vector3D &ray_origin,
                      const Vector3D &ray_direction) const;

  // Surface normal calculation
  Vector3D getSurfaceNormal(const Vector3D &surface_point,
                            double tolerance = 1e-10) const;
  Vector3D getSurfaceNormal(BoxFace face) const;

  // Boundary analysis
  BoxFace getClosestFace(const Vector3D &point) const;
  std::optional<BoxFace> getFaceContainingPoint(const Vector3D &point,
                                                double tolerance = 1e-10) const;

  // Geometric transformations
  Box translate(const Vector3D &offset) const;
  Box scale(double factor) const;
  Box scale(const Vector3D &factors) const;
  Box expand(double amount) const;
  Box expand(const Vector3D &amounts) const;

  // Intersection and union with other boxes
  std::optional<Box> intersection(const Box &other) const;
  Box boundingBox(const Box &other) const;
  bool intersects(const Box &other) const;
  bool contains(const Box &other) const;

  // Point generation
  Vector3D randomPointInside() const; // Requires RNG in implementation
  Vector3D randomPointOnSurface() const;
  std::vector<Vector3D> gridPoints(int nx, int ny, int nz) const;

  // Validation
  bool hasValidDimensions() const;
  std::vector<std::string> validate() const;

  // String representation
  std::string toString() const;
  std::string getGeometryString() const;

  // Comparison operators
  bool operator==(const Box &other) const;
  bool operator!=(const Box &other) const;

  // Static factory methods
  static Box createCube(const Vector3D &centre, double side_length);
  static Box createCube(const Vector3D &centre, double side_length,
                        const Material &material);
  static Box createFromCentreAndDimensions(const Vector3D &centre,
                                           const Vector3D &dimensions);
  static Box createFromCentreAndDimensions(const Vector3D &centre,
                                           const Vector3D &dimensions,
                                           const Material &material);

  // Utility methods for common geometries
  static Box createLeadShield(const Vector3D &min_corner,
                              const Vector3D &max_corner,
                              double density = 11.34);
  static Box createConcreteShield(const Vector3D &min_corner,
                                  const Vector3D &max_corner,
                                  double density = 2.3);
  static Box createWaterPhantom(const Vector3D &min_corner,
                                const Vector3D &max_corner);
  static Box createAirRegion(const Vector3D &min_corner,
                             const Vector3D &max_corner);

private:
  void invalidateCache() const;
  void validateAndSetState();
  RayBoxIntersection computeIntersection(const Vector3D &ray_origin,
                                         const Vector3D &ray_direction,
                                         double t, BoxFace face,
                                         bool entering) const;
  std::tuple<double, double> intersectSlab(double ray_origin, double ray_dir,
                                           double slab_min,
                                           double slab_max) const;
};

// Utility functions
std::string boxFaceToString(BoxFace face);
Vector3D getBoxFaceNormal(BoxFace face);
std::optional<BoxFace> stringToBoxFace(const std::string &face_str);

// Box intersection utilities
bool boxesIntersect(const Box &box1, const Box &box2);
std::optional<Box> intersectBoxes(const Box &box1, const Box &box2);
Box boundingBoxOf(const std::vector<Box> &boxes);