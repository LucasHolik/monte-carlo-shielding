#include "Geometry.hpp"
#include "RandomNumberGenerator.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

// Constructors
Geometry::Geometry()
    : default_material_(Material::createVacuum()), name_("Geometry")
{}

Geometry::Geometry(std::string_view name)
    : default_material_(Material::createVacuum()), name_(name)
{}

Geometry::Geometry(const Material &default_material)
    : default_material_(default_material), name_("Geometry")
{}

Geometry::Geometry(std::string_view name, const Material &default_material)
    : default_material_(default_material), name_(name)
{}

// Basic mutators
void Geometry::setDefaultMaterial(const Material &material)
{
  default_material_ = material;
}

// Shape management
std::size_t Geometry::addShape(const GeometryShape &shape,
                               std::string_view name)
{
  std::size_t index = shapes_.size();
  shapes_.push_back(shape);

  std::string shape_name =
      name.empty() ? generateDefaultName(index) : std::string(name);
  shape_names_.push_back(shape_name);

  updateNameIndex();
  invalidateSpatialIndex();

  return index;
}

std::size_t Geometry::addBox(const Box &box, std::string_view name)
{
  return addShape(GeometryShape{box}, name);
}

bool Geometry::removeShape(std::size_t index)
{
  if(index >= shapes_.size())
  {
    return false;
  }

  shapes_.erase(shapes_.begin() + index);
  shape_names_.erase(shape_names_.begin() + index);

  updateNameIndex();
  invalidateSpatialIndex();

  return true;
}

bool Geometry::removeShape(std::string_view name)
{
  auto it = name_to_index_.find(std::string(name));
  if(it == name_to_index_.end())
  {
    return false;
  }

  return removeShape(it->second);
}

void Geometry::clearShapes()
{
  shapes_.clear();
  shape_names_.clear();
  name_to_index_.clear();
  invalidateSpatialIndex();
}

// Shape access
std::optional<GeometryShape> Geometry::getShape(std::size_t index) const
{
  if(index >= shapes_.size())
  {
    return std::nullopt;
  }
  return shapes_[index];
}

std::optional<GeometryShape> Geometry::getShape(std::string_view name) const
{
  auto index = getShapeIndex(name);
  if(!index)
  {
    return std::nullopt;
  }
  return getShape(*index);
}

std::optional<std::size_t> Geometry::getShapeIndex(std::string_view name) const
{
  auto it = name_to_index_.find(std::string(name));
  if(it == name_to_index_.end())
  {
    return std::nullopt;
  }
  return it->second;
}

std::string Geometry::getShapeName(std::size_t index) const
{
  if(index >= shape_names_.size())
  {
    return "";
  }
  return shape_names_[index];
}

// Ray intersection methods
std::optional<GeometryIntersection>
Geometry::findClosestIntersection(const Vector3D &ray_origin,
                                  const Vector3D &ray_direction) const
{
  if(shapes_.empty())
  {
    return std::nullopt;
  }

  std::optional<GeometryIntersection> closest_intersection;
  double min_distance = std::numeric_limits<double>::infinity();

  for(std::size_t i = 0; i < shapes_.size(); ++i)
  {
    auto intersection =
        computeShapeIntersection(shapes_[i], i, ray_origin, ray_direction);

    if(intersection.intersects && intersection.distance >= 0.0 &&
       intersection.distance < min_distance)
    {
      min_distance = intersection.distance;
      closest_intersection = intersection;
    }
  }

  return closest_intersection;
}

std::vector<GeometryIntersection>
Geometry::findAllIntersections(const Vector3D &ray_origin,
                               const Vector3D &ray_direction) const
{
  std::vector<GeometryIntersection> intersections;

  for(std::size_t i = 0; i < shapes_.size(); ++i)
  {
    auto intersection =
        computeShapeIntersection(shapes_[i], i, ray_origin, ray_direction);

    if(intersection.intersects && intersection.distance >= 0.0)
    {
      intersections.push_back(intersection);
    }
  }

  // Sort by distance
  std::sort(intersections.begin(), intersections.end(),
            [](const GeometryIntersection &a, const GeometryIntersection &b) {
              return a.distance < b.distance;
            });

  return intersections;
}

std::optional<double>
Geometry::distanceToSurface(const Vector3D &position,
                            const Vector3D &direction) const
{
  auto intersection = findClosestIntersection(position, direction);
  if(intersection && intersection->intersects)
  {
    return intersection->distance;
  }
  return std::nullopt;
}

// Point location methods
LocationResult Geometry::locatePoint(const Vector3D &point) const
{
  // Check shapes in order - first match wins
  for(std::size_t i = 0; i < shapes_.size(); ++i)
  {
    if(isPointInsideShape(shapes_[i], point))
    {
      Material material = getShapeMaterial(shapes_[i]);
      std::string name = getShapeName(i);
      return LocationResult(material, i, name, point);
    }
  }

  // Point not inside any shape - return default material
  return LocationResult();
}

Material Geometry::getMaterialAtPoint(const Vector3D &point) const
{
  auto location = locatePoint(point);
  return location.inside_geometry ? location.material : default_material_;
}

bool Geometry::isInsideAnyShape(const Vector3D &point) const
{
  return locatePoint(point).inside_geometry;
}

// Boundary detection
bool Geometry::isOnBoundary(const Vector3D &point, double tolerance) const
{
  for(const auto &shape : shapes_)
  {
    bool on_boundary = std::visit(
        [&](const auto &s) { return s.isOnSurface(point, tolerance); }, shape);

    if(on_boundary)
    {
      return true;
    }
  }
  return false;
}

Vector3D Geometry::getSurfaceNormal(const Vector3D &surface_point,
                                    double tolerance) const
{
  for(const auto &shape : shapes_)
  {
    bool on_surface = std::visit(
        [&](const auto &s) { return s.isOnSurface(surface_point, tolerance); },
        shape);

    if(on_surface)
    {
      return std::visit(
          [&](const auto &s) {
            return s.getSurfaceNormal(surface_point, tolerance);
          },
          shape);
    }
  }

  // Default to upward normal if not on any surface
  return Vector3D::UNIT_Z;
}

// Particle transport utilities
std::optional<GeometryIntersection>
Geometry::transportParticle(const Vector3D &start_position,
                            const Vector3D &direction,
                            double max_distance) const
{
  auto intersection = findClosestIntersection(start_position, direction);

  if(intersection && intersection->distance <= max_distance)
  {
    return intersection;
  }

  return std::nullopt;
}

bool Geometry::willParticleEscape(const Vector3D &position,
                                  const Vector3D &direction) const
{
  // If no intersection found, particle escapes
  auto intersection = findClosestIntersection(position, direction);
  return !intersection.has_value();
}

std::vector<GeometryIntersection>
Geometry::getParticlePath(const Vector3D &start_position,
                          const Vector3D &direction, double max_distance) const
{
  std::vector<GeometryIntersection> path;

  Vector3D current_position = start_position;
  double remaining_distance = max_distance;
  const double epsilon = 1e-10; // Small step to avoid numerical issues

  while(remaining_distance > 0.0)
  {
    auto intersection = findClosestIntersection(current_position, direction);

    if(!intersection || intersection->distance > remaining_distance)
    {
      // No more intersections within remaining distance
      break;
    }

    path.push_back(*intersection);

    // Move slightly past the intersection point
    current_position = intersection->intersection_point + direction * epsilon;
    remaining_distance -= (intersection->distance + epsilon);
  }

  return path;
}

// Advanced intersection queries
std::vector<std::size_t>
Geometry::findIntersectingShapes(const Vector3D &ray_origin,
                                 const Vector3D &ray_direction) const
{
  std::vector<std::size_t> intersecting_shapes;

  for(std::size_t i = 0; i < shapes_.size(); ++i)
  {
    auto intersection =
        computeShapeIntersection(shapes_[i], i, ray_origin, ray_direction);

    if(intersection.intersects && intersection.distance >= 0.0)
    {
      intersecting_shapes.push_back(i);
    }
  }

  return intersecting_shapes;
}

std::optional<GeometryIntersection>
Geometry::findIntersectionWithShape(std::size_t shape_index,
                                    const Vector3D &ray_origin,
                                    const Vector3D &ray_direction) const
{
  if(shape_index >= shapes_.size())
  {
    return std::nullopt;
  }

  auto intersection = computeShapeIntersection(
      shapes_[shape_index], shape_index, ray_origin, ray_direction);

  return intersection.intersects
             ? std::optional<GeometryIntersection>{intersection}
             : std::nullopt;
}

// Geometry analysis
std::vector<Vector3D> Geometry::getBoundingBox() const
{
  if(shapes_.empty())
  {
    return {Vector3D::ZERO, Vector3D::ZERO};
  }

  Vector3D min_corner = getShapeBoundingBoxMin(shapes_[0]);
  Vector3D max_corner = getShapeBoundingBoxMax(shapes_[0]);

  for(std::size_t i = 1; i < shapes_.size(); ++i)
  {
    Vector3D shape_min = getShapeBoundingBoxMin(shapes_[i]);
    Vector3D shape_max = getShapeBoundingBoxMax(shapes_[i]);

    min_corner = Vector3D(std::min(min_corner.x(), shape_min.x()),
                          std::min(min_corner.y(), shape_min.y()),
                          std::min(min_corner.z(), shape_min.z()));

    max_corner = Vector3D(std::max(max_corner.x(), shape_max.x()),
                          std::max(max_corner.y(), shape_max.y()),
                          std::max(max_corner.z(), shape_max.z()));
  }

  return {min_corner, max_corner};
}

double Geometry::getTotalVolume() const
{
  double total_volume = 0.0;
  for(const auto &shape : shapes_)
  {
    total_volume += getShapeVolume(shape);
  }
  return total_volume;
}

std::vector<Material> Geometry::getAllMaterials() const
{
  std::vector<Material> materials;
  materials.push_back(default_material_);

  for(const auto &shape : shapes_)
  {
    Material shape_material = getShapeMaterial(shape);

    // Add if not already present
    bool found = false;
    for(const auto &existing : materials)
    {
      if(existing == shape_material)
      {
        found = true;
        break;
      }
    }

    if(!found)
    {
      materials.push_back(shape_material);
    }
  }

  return materials;
}

// Validation
bool Geometry::isValid() const
{
  auto errors = validate();
  return errors.empty();
}

std::vector<std::string> Geometry::validate() const
{
  std::vector<std::string> errors;

  if(!default_material_.isValid())
  {
    errors.push_back("Invalid default material");
  }

  for(std::size_t i = 0; i < shapes_.size(); ++i)
  {
    bool shape_valid = std::visit(
        [](const auto &shape) { return shape.isValid(); }, shapes_[i]);

    if(!shape_valid)
    {
      errors.push_back("Invalid shape at index " + std::to_string(i));
    }
  }

  if(shapes_.size() != shape_names_.size())
  {
    errors.push_back("Inconsistent shape and name vectors");
  }

  return errors;
}

bool Geometry::hasOverlappingShapes(double tolerance) const
{
  for(std::size_t i = 0; i < shapes_.size(); ++i)
  {
    for(std::size_t j = i + 1; j < shapes_.size(); ++j)
    {
      // Check if shapes i and j overlap
      // For now, we'll do a simple bounding box check
      Vector3D min1 = getShapeBoundingBoxMin(shapes_[i]);
      Vector3D max1 = getShapeBoundingBoxMax(shapes_[i]);
      Vector3D min2 = getShapeBoundingBoxMin(shapes_[j]);
      Vector3D max2 = getShapeBoundingBoxMax(shapes_[j]);

      bool overlap = (min1.x() <= max2.x() + tolerance &&
                      max1.x() >= min2.x() - tolerance) &&
                     (min1.y() <= max2.y() + tolerance &&
                      max1.y() >= min2.y() - tolerance) &&
                     (min1.z() <= max2.z() + tolerance &&
                      max1.z() >= min2.z() - tolerance);

      if(overlap)
      {
        return true;
      }
    }
  }
  return false;
}

// Statistical and sampling methods
Vector3D Geometry::generateRandomPoint() const
{
  if(shapes_.empty())
  {
    return Vector3D::ZERO;
  }

  auto &rng = RandomNumberGenerator::getThreadLocal();

  // Simple approach: pick random shape and generate point inside it
  std::size_t shape_index = rng.uniformSize(0, shapes_.size() - 1);
  return generateRandomPointInShape(shape_index);
}

Vector3D Geometry::generateRandomPointInShape(std::size_t shape_index) const
{
  if(shape_index >= shapes_.size())
  {
    return Vector3D::ZERO;
  }

  return std::visit([](const auto &shape) { return shape.randomPointInside(); },
                    shapes_[shape_index]);
}

Vector3D Geometry::generateRandomPointOnSurface() const
{
  if(shapes_.empty())
  {
    return Vector3D::ZERO;
  }

  auto &rng = RandomNumberGenerator::getThreadLocal();

  // Weight by surface area and pick random shape
  std::vector<double> surface_areas;
  for(const auto &shape : shapes_)
  {
    double area =
        std::visit([](const auto &s) { return s.surfaceArea(); }, shape);
    surface_areas.push_back(area);
  }

  std::size_t shape_index = rng.discreteSample(surface_areas);

  return std::visit(
      [](const auto &shape) { return shape.randomPointOnSurface(); },
      shapes_[shape_index]);
}

// Geometric transformations
void Geometry::translateAll(const Vector3D &offset)
{
  for(auto &shape : shapes_)
  {
    std::visit([&offset](auto &s) { s = s.translate(offset); }, shape);
  }
  invalidateSpatialIndex();
}

void Geometry::scaleAll(double factor)
{
  for(auto &shape : shapes_)
  {
    std::visit([factor](auto &s) { s = s.scale(factor); }, shape);
  }
  invalidateSpatialIndex();
}

// I/O and string representation
std::string Geometry::toString() const
{
  std::stringstream ss;
  ss << "Geometry '" << name_ << "'\n";
  ss << "  Number of shapes: " << shapes_.size() << "\n";
  ss << "  Default material: " << default_material_.name() << "\n";

  if(!shapes_.empty())
  {
    ss << "  Total volume: " << getTotalVolume() << "\n";

    auto bounding_box = getBoundingBox();
    ss << "  Bounding box: (" << bounding_box[0].x() << ", "
       << bounding_box[0].y() << ", " << bounding_box[0].z() << ") to ("
       << bounding_box[1].x() << ", " << bounding_box[1].y() << ", "
       << bounding_box[1].z() << ")\n";

    ss << "  Shapes:\n";
    for(std::size_t i = 0; i < shapes_.size(); ++i)
    {
      ss << "    [" << i << "] " << shape_names_[i] << ": ";
      std::visit([&ss](const auto &shape) { ss << shape.getGeometryString(); },
                 shapes_[i]);
      ss << "\n";
    }
  }

  return ss.str();
}

std::string Geometry::getStatistics() const
{
  std::stringstream ss;
  ss << "Geometry Statistics:\n";
  ss << "  Total shapes: " << shapes_.size() << "\n";
  ss << "  Total volume: " << getTotalVolume() << "\n";

  auto materials = getAllMaterials();
  ss << "  Unique materials: " << materials.size() << "\n";

  if(hasOverlappingShapes())
  {
    ss << "  Warning: Overlapping shapes detected\n";
  }

  return ss.str();
}

// Comparison operators
bool Geometry::operator==(const Geometry &other) const
{
  return name_ == other.name_ && default_material_ == other.default_material_ &&
         shapes_ == other.shapes_ && shape_names_ == other.shape_names_;
}

bool Geometry::operator!=(const Geometry &other) const
{
  return !(*this == other);
}

// Static factory methods
Geometry Geometry::createSimpleShield(const Vector3D &shield_min,
                                      const Vector3D &shield_max,
                                      const Material &shield_material,
                                      const Material &surrounding_material)
{
  Geometry geometry("Simple Shield", surrounding_material);

  Box shield_box(shield_min, shield_max, shield_material, "Shield");
  geometry.addBox(shield_box, "Shield");

  return geometry;
}

Geometry Geometry::createLayeredShield(const std::vector<Box> &layers,
                                       const Material &surrounding_material)
{
  Geometry geometry("Layered Shield", surrounding_material);

  for(std::size_t i = 0; i < layers.size(); ++i)
  {
    std::string layer_name = "Layer_" + std::to_string(i);
    geometry.addBox(layers[i], layer_name);
  }

  return geometry;
}

Geometry Geometry::createDetectorGeometry(const Vector3D &detector_position,
                                          const Vector3D &detector_size,
                                          const Material &detector_material,
                                          const Material &surrounding_material)
{
  Geometry geometry("Detector Geometry", surrounding_material);

  Vector3D half_size = detector_size * 0.5;
  Vector3D detector_min = detector_position - half_size;
  Vector3D detector_max = detector_position + half_size;

  Box detector_box(detector_min, detector_max, detector_material, "Detector");
  geometry.addBox(detector_box, "Detector");

  return geometry;
}

Geometry Geometry::createShieldingExperiment(double shield_thickness_cm,
                                            const Material &shield_material,
                                            double experiment_length,
                                            double cross_section)
{
  // Create air-shield-air configuration for transmission experiments
  Material air = Material::createAir();
  Geometry geometry("Shielding Experiment", air);

  double half_length = experiment_length * 0.5;
  double half_cross = cross_section * 0.5;

  // Region 1: Air before shield (z=-half_length to z=0)
  Box air_before(Vector3D(-half_cross, -half_cross, -half_length),
                 Vector3D(half_cross, half_cross, 0.0),
                 air, "Air_Before");
  
  // Region 2: Shield material (z=0 to z=shield_thickness_cm)  
  Box shield(Vector3D(-half_cross, -half_cross, 0.0),
             Vector3D(half_cross, half_cross, shield_thickness_cm),
             shield_material, "Shield");
             
  // Region 3: Air after shield (z=shield_thickness_cm to z=half_length)
  Box air_after(Vector3D(-half_cross, -half_cross, shield_thickness_cm),
                Vector3D(half_cross, half_cross, half_length),
                air, "Air_After");

  geometry.addBox(air_before, "Air_Before");
  geometry.addBox(shield, "Shield");  
  geometry.addBox(air_after, "Air_After");

  return geometry;
}

// Advanced features
void Geometry::enableSpatialIndexing(bool enable)
{
  // Future enhancement - spatial indexing for performance
  // For now, just invalidate the index
  invalidateSpatialIndex();
}

void Geometry::rebuildSpatialIndex() const
{
  // Future enhancement - rebuild spatial acceleration structure
  spatial_index_valid_ = true;
}

// Private helper methods
void Geometry::updateNameIndex()
{
  name_to_index_.clear();
  for(std::size_t i = 0; i < shape_names_.size(); ++i)
  {
    name_to_index_[shape_names_[i]] = i;
  }
}

std::string Geometry::generateDefaultName(std::size_t index) const
{
  return "Shape_" + std::to_string(index);
}

GeometryIntersection Geometry::computeShapeIntersection(
    const GeometryShape &shape, std::size_t index, const Vector3D &ray_origin,
    const Vector3D &ray_direction) const
{
  return std::visit(
      [&](const auto &s) -> GeometryIntersection {
        auto intersection = s.tryRayIntersection(ray_origin, ray_direction);

        if(!intersection || !intersection->intersects)
        {
          return GeometryIntersection();
        }

        Material material = s.material();
        std::string name = getShapeName(index);

        return GeometryIntersection(intersection->distance,
                                    intersection->intersection_point,
                                    intersection->surface_normal, material,
                                    index, name, intersection->entering);
      },
      shape);
}

bool Geometry::isPointInsideShape(const GeometryShape &shape,
                                  const Vector3D &point) const
{
  return std::visit([&point](const auto &s) { return s.isInside(point); },
                    shape);
}

Material Geometry::getShapeMaterial(const GeometryShape &shape) const
{
  return std::visit([](const auto &s) { return s.material(); }, shape);
}

Vector3D Geometry::getShapeBoundingBoxMin(const GeometryShape &shape) const
{
  return std::visit([](const auto &s) { return s.minCorner(); }, shape);
}

Vector3D Geometry::getShapeBoundingBoxMax(const GeometryShape &shape) const
{
  return std::visit([](const auto &s) { return s.maxCorner(); }, shape);
}

double Geometry::getShapeVolume(const GeometryShape &shape) const
{
  return std::visit([](const auto &s) { return s.volume(); }, shape);
}

// Utility functions
namespace GeometryUtils
{

bool geometriesIntersect(const Geometry &geom1, const Geometry &geom2)
{
  auto bb1 = geom1.getBoundingBox();
  auto bb2 = geom2.getBoundingBox();

  // Simple bounding box intersection test
  return (bb1[0].x() <= bb2[1].x() && bb1[1].x() >= bb2[0].x()) &&
         (bb1[0].y() <= bb2[1].y() && bb1[1].y() >= bb2[0].y()) &&
         (bb1[0].z() <= bb2[1].z() && bb1[1].z() >= bb2[0].z());
}

Geometry combineGeometries(const Geometry &geom1, const Geometry &geom2)
{
  Geometry combined("Combined Geometry", geom1.defaultMaterial());

  // Add all shapes from first geometry
  for(std::size_t i = 0; i < geom1.getNumberOfShapes(); ++i)
  {
    auto shape = geom1.getShape(i);
    if(shape)
    {
      std::string name = geom1.getShapeName(i) + "_from_geom1";
      combined.addShape(*shape, name);
    }
  }

  // Add all shapes from second geometry
  for(std::size_t i = 0; i < geom2.getNumberOfShapes(); ++i)
  {
    auto shape = geom2.getShape(i);
    if(shape)
    {
      std::string name = geom2.getShapeName(i) + "_from_geom2";
      combined.addShape(*shape, name);
    }
  }

  return combined;
}

double calculateOverlapVolume(const Geometry &geom1, const Geometry &geom2)
{
  // Simplified implementation - would need more sophisticated algorithm
  // for accurate overlap calculation
  return 0.0;
}

std::vector<Vector3D> findIntersectionPoints(const Geometry &geom1,
                                             const Geometry &geom2)
{
  // Simplified implementation - would need sophisticated algorithm
  // to find actual intersection points between geometries
  return {};
}

bool validateGeometryForTransport(const Geometry &geometry)
{
  return geometry.isValid() && !geometry.hasOverlappingShapes();
}

std::vector<std::string> checkGeometryIntegrity(const Geometry &geometry)
{
  std::vector<std::string> issues;

  auto validation_errors = geometry.validate();
  issues.insert(issues.end(), validation_errors.begin(),
                validation_errors.end());

  if(geometry.hasOverlappingShapes())
  {
    issues.push_back("Geometry contains overlapping shapes");
  }

  if(geometry.isEmpty())
  {
    issues.push_back("Geometry contains no shapes");
  }

  return issues;
}

} // namespace GeometryUtils