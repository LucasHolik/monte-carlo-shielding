#include "Box.hpp"
#include "RandomNumberGenerator.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

// Constructors
Box::Box()
    : min_corner_(Vector3D::ZERO), max_corner_(Vector3D(1.0, 1.0, 1.0)),
      material_(Material::createAir()), name_("Default Box"), is_valid_(true)
{
  validateAndSetState();
}

Box::Box(const Vector3D &min_corner, const Vector3D &max_corner)
    : min_corner_(min_corner), max_corner_(max_corner),
      material_(Material::createAir()), name_("Box"), is_valid_(false)
{
  validateAndSetState();
}

Box::Box(const Vector3D &min_corner, const Vector3D &max_corner,
         const Material &material)
    : min_corner_(min_corner), max_corner_(max_corner), material_(material),
      name_("Box"), is_valid_(false)
{
  validateAndSetState();
}

Box::Box(const Vector3D &min_corner, const Vector3D &max_corner,
         const Material &material, std::string_view name)
    : min_corner_(min_corner), max_corner_(max_corner), material_(material),
      name_(name), is_valid_(false)
{
  validateAndSetState();
}

// Computed properties
Vector3D Box::centre() const
{
  if(!centre_.has_value())
  {
    centre_ = (min_corner_ + max_corner_) * 0.5;
  }
  return *centre_;
}

Vector3D Box::dimensions() const
{
  if(!dimensions_.has_value())
  {
    dimensions_ = max_corner_ - min_corner_;
  }
  return *dimensions_;
}

double Box::volume() const
{
  if(!volume_.has_value())
  {
    Vector3D dim = dimensions();
    volume_ = dim.x() * dim.y() * dim.z();
  }
  return *volume_;
}

double Box::surfaceArea() const
{
  Vector3D dim = dimensions();
  return 2.0 * (dim.x() * dim.y() + dim.y() * dim.z() + dim.z() * dim.x());
}

// Basic mutators
void Box::setMinCorner(const Vector3D &min_corner)
{
  min_corner_ = min_corner;
  invalidateCache();
  validateAndSetState();
}

void Box::setMaxCorner(const Vector3D &max_corner)
{
  max_corner_ = max_corner;
  invalidateCache();
  validateAndSetState();
}

void Box::setCorners(const Vector3D &min_corner, const Vector3D &max_corner)
{
  min_corner_ = min_corner;
  max_corner_ = max_corner;
  invalidateCache();
  validateAndSetState();
}

// Geometry queries
bool Box::isInside(const Vector3D &point, double tolerance) const
{
  return point.x() >= (min_corner_.x() - tolerance) &&
         point.x() <= (max_corner_.x() + tolerance) &&
         point.y() >= (min_corner_.y() - tolerance) &&
         point.y() <= (max_corner_.y() + tolerance) &&
         point.z() >= (min_corner_.z() - tolerance) &&
         point.z() <= (max_corner_.z() + tolerance);
}

bool Box::isOnSurface(const Vector3D &point, double tolerance) const
{
  if(!isInside(point, tolerance))
  {
    return false;
  }

  // Check if point is close to any face
  bool near_x_min = std::abs(point.x() - min_corner_.x()) <= tolerance;
  bool near_x_max = std::abs(point.x() - max_corner_.x()) <= tolerance;
  bool near_y_min = std::abs(point.y() - min_corner_.y()) <= tolerance;
  bool near_y_max = std::abs(point.y() - max_corner_.y()) <= tolerance;
  bool near_z_min = std::abs(point.z() - min_corner_.z()) <= tolerance;
  bool near_z_max = std::abs(point.z() - max_corner_.z()) <= tolerance;

  return near_x_min || near_x_max || near_y_min || near_y_max || near_z_min ||
         near_z_max;
}

bool Box::isOutside(const Vector3D &point, double tolerance) const
{
  return !isInside(point, tolerance);
}

// Distance calculations
double Box::distanceToPoint(const Vector3D &point) const
{
  return std::sqrt(distanceToPointSquared(point));
}

double Box::distanceToPointSquared(const Vector3D &point) const
{
  double dx =
      std::max({0.0, min_corner_.x() - point.x(), point.x() - max_corner_.x()});
  double dy =
      std::max({0.0, min_corner_.y() - point.y(), point.y() - max_corner_.y()});
  double dz =
      std::max({0.0, min_corner_.z() - point.z(), point.z() - max_corner_.z()});

  return dx * dx + dy * dy + dz * dz;
}

// Ray intersection implementation using slab method
std::optional<RayBoxIntersection>
Box::tryRayIntersection(const Vector3D &ray_origin,
                        const Vector3D &ray_direction) const
{
  if(!is_valid_)
  {
    return std::nullopt;
  }

  // Check for zero direction vector
  if(ray_direction.isZero())
  {
    return std::nullopt;
  }

  // Use slab method for ray-box intersection
  double t_min = -std::numeric_limits<double>::infinity();
  double t_max = std::numeric_limits<double>::infinity();

  BoxFace entry_face = BoxFace::XMin;
  BoxFace exit_face = BoxFace::XMax;

  // Check intersection with each pair of parallel planes (slabs)
  for(int axis = 0; axis < 3; ++axis)
  {
    double ray_orig, ray_dir, box_min, box_max;

    if(axis == 0)
    { // X axis
      ray_orig = ray_origin.x();
      ray_dir = ray_direction.x();
      box_min = min_corner_.x();
      box_max = max_corner_.x();
    }
    else if(axis == 1)
    { // Y axis
      ray_orig = ray_origin.y();
      ray_dir = ray_direction.y();
      box_min = min_corner_.y();
      box_max = max_corner_.y();
    }
    else
    { // Z axis
      ray_orig = ray_origin.z();
      ray_dir = ray_direction.z();
      box_min = min_corner_.z();
      box_max = max_corner_.z();
    }

    if(std::abs(ray_dir) < 1e-15)
    {
      // Ray is parallel to the slab planes
      if(ray_orig < box_min || ray_orig > box_max)
      {
        return std::nullopt; // Ray misses the box
      }
    }
    else
    {
      // Calculate intersection distances
      double t1 = (box_min - ray_orig) / ray_dir;
      double t2 = (box_max - ray_orig) / ray_dir;

      // Ensure t1 is the near intersection and t2 is the far intersection
      if(t1 > t2)
      {
        std::swap(t1, t2);
      }

      // Update the intersection interval
      if(t1 > t_min)
      {
        t_min = t1;
        // Determine which face we're hitting
        if(axis == 0)
        {
          entry_face = (ray_dir > 0) ? BoxFace::XMin : BoxFace::XMax;
        }
        else if(axis == 1)
        {
          entry_face = (ray_dir > 0) ? BoxFace::YMin : BoxFace::YMax;
        }
        else
        {
          entry_face = (ray_dir > 0) ? BoxFace::ZMin : BoxFace::ZMax;
        }
      }

      if(t2 < t_max)
      {
        t_max = t2;
        // Determine exit face
        if(axis == 0)
        {
          exit_face = (ray_dir > 0) ? BoxFace::XMax : BoxFace::XMin;
        }
        else if(axis == 1)
        {
          exit_face = (ray_dir > 0) ? BoxFace::YMax : BoxFace::YMin;
        }
        else
        {
          exit_face = (ray_dir > 0) ? BoxFace::ZMax : BoxFace::ZMin;
        }
      }

      // Check if intersection interval is empty
      if(t_min > t_max)
      {
        return std::nullopt;
      }
    }
  }

  // Determine the closest intersection (prefer positive t values)
  double intersection_distance;
  BoxFace intersection_face;
  bool entering;

  if(t_min >= 0.0)
  {
    // Ray starts outside the box, hits entry face
    intersection_distance = t_min;
    intersection_face = entry_face;
    entering = true;
  }
  else if(t_max >= 0.0)
  {
    // Ray starts inside the box, hits exit face
    intersection_distance = t_max;
    intersection_face = exit_face;
    entering = false;
  }
  else
  {
    // Box is behind the ray
    return std::nullopt;
  }

  return computeIntersection(ray_origin, ray_direction, intersection_distance,
                             intersection_face, entering);
}

RayBoxIntersection Box::rayIntersection(const Vector3D &ray_origin,
                                        const Vector3D &ray_direction) const
{
  auto result = tryRayIntersection(ray_origin, ray_direction);
  if(!result)
  {
    throw std::invalid_argument(
        "Ray does not intersect box or invalid parameters");
  }
  return *result;
}

std::optional<double> Box::distanceToSurface(const Vector3D &position,
                                             const Vector3D &direction) const
{
  auto intersection = tryRayIntersection(position, direction);
  if(intersection && intersection->intersects)
  {
    return intersection->distance;
  }
  return std::nullopt;
}

std::vector<RayBoxIntersection>
Box::allRayIntersections(const Vector3D &ray_origin,
                         const Vector3D &ray_direction) const
{
  std::vector<RayBoxIntersection> intersections;

  if(!is_valid_ || ray_direction.isZero())
  {
    return intersections;
  }

  // Calculate all potential intersection points
  std::vector<std::pair<double, BoxFace>> candidates;

  // Check each face
  for(int axis = 0; axis < 3; ++axis)
  {
    double ray_orig, ray_dir, box_min, box_max;

    if(axis == 0)
    {
      ray_orig = ray_origin.x();
      ray_dir = ray_direction.x();
      box_min = min_corner_.x();
      box_max = max_corner_.x();
    }
    else if(axis == 1)
    {
      ray_orig = ray_origin.y();
      ray_dir = ray_direction.y();
      box_min = min_corner_.y();
      box_max = max_corner_.y();
    }
    else
    {
      ray_orig = ray_origin.z();
      ray_dir = ray_direction.z();
      box_min = min_corner_.z();
      box_max = max_corner_.z();
    }

    if(std::abs(ray_dir) > 1e-15)
    {
      double t_min = (box_min - ray_orig) / ray_dir;
      double t_max = (box_max - ray_orig) / ray_dir;

      BoxFace min_face, max_face;
      if(axis == 0)
      {
        min_face = BoxFace::XMin;
        max_face = BoxFace::XMax;
      }
      else if(axis == 1)
      {
        min_face = BoxFace::YMin;
        max_face = BoxFace::YMax;
      }
      else
      {
        min_face = BoxFace::ZMin;
        max_face = BoxFace::ZMax;
      }

      candidates.emplace_back(t_min, min_face);
      candidates.emplace_back(t_max, max_face);
    }
  }

  // Check which candidates are valid intersections
  for(const auto &[t, face] : candidates)
  {
    if(t >= 0.0)
    { // Only forward intersections
      Vector3D intersection_point = ray_origin + ray_direction * t;

      // Check if intersection point is actually on the box surface
      constexpr double tolerance = 1e-10;
      if(isOnSurface(intersection_point, tolerance))
      {
        bool entering = !isInside(ray_origin);
        intersections.emplace_back(
            computeIntersection(ray_origin, ray_direction, t, face, entering));
      }
    }
  }

  // Sort by distance
  std::sort(intersections.begin(), intersections.end(),
            [](const RayBoxIntersection &a, const RayBoxIntersection &b) {
              return a.distance < b.distance;
            });

  return intersections;
}

// Surface normal calculation
Vector3D Box::getSurfaceNormal(const Vector3D &surface_point,
                               double tolerance) const
{
  // Find which face the point is closest to
  auto face = getFaceContainingPoint(surface_point, tolerance);
  if(face)
  {
    return getSurfaceNormal(*face);
  }

  // If not on surface, find closest face
  BoxFace closest_face = getClosestFace(surface_point);
  return getSurfaceNormal(closest_face);
}

Vector3D Box::getSurfaceNormal(BoxFace face) const
{
  return getBoxFaceNormal(face);
}

// Boundary analysis
BoxFace Box::getClosestFace(const Vector3D &point) const
{
  double min_distance = std::numeric_limits<double>::infinity();
  BoxFace closest_face = BoxFace::XMin;

  // Calculate distance to each face surface (not just the plane)
  std::array<std::pair<double, BoxFace>, 6> face_distances;

  // XMin face: project point onto face and calculate distance
  Vector3D xmin_projected(
      min_corner_.x(),
      std::max(min_corner_.y(), std::min(max_corner_.y(), point.y())),
      std::max(min_corner_.z(), std::min(max_corner_.z(), point.z())));
  face_distances[0] = {point.distance(xmin_projected), BoxFace::XMin};

  // XMax face
  Vector3D xmax_projected(
      max_corner_.x(),
      std::max(min_corner_.y(), std::min(max_corner_.y(), point.y())),
      std::max(min_corner_.z(), std::min(max_corner_.z(), point.z())));
  face_distances[1] = {point.distance(xmax_projected), BoxFace::XMax};

  // YMin face
  Vector3D ymin_projected(
      std::max(min_corner_.x(), std::min(max_corner_.x(), point.x())),
      min_corner_.y(),
      std::max(min_corner_.z(), std::min(max_corner_.z(), point.z())));
  face_distances[2] = {point.distance(ymin_projected), BoxFace::YMin};

  // YMax face
  Vector3D ymax_projected(
      std::max(min_corner_.x(), std::min(max_corner_.x(), point.x())),
      max_corner_.y(),
      std::max(min_corner_.z(), std::min(max_corner_.z(), point.z())));
  face_distances[3] = {point.distance(ymax_projected), BoxFace::YMax};

  // ZMin face
  Vector3D zmin_projected(
      std::max(min_corner_.x(), std::min(max_corner_.x(), point.x())),
      std::max(min_corner_.y(), std::min(max_corner_.y(), point.y())),
      min_corner_.z());
  face_distances[4] = {point.distance(zmin_projected), BoxFace::ZMin};

  // ZMax face
  Vector3D zmax_projected(
      std::max(min_corner_.x(), std::min(max_corner_.x(), point.x())),
      std::max(min_corner_.y(), std::min(max_corner_.y(), point.y())),
      max_corner_.z());
  face_distances[5] = {point.distance(zmax_projected), BoxFace::ZMax};

  // Find the face with minimum distance
  for(const auto &[distance, face] : face_distances)
  {
    if(distance < min_distance)
    {
      min_distance = distance;
      closest_face = face;
    }
  }

  return closest_face;
}

std::optional<BoxFace> Box::getFaceContainingPoint(const Vector3D &point,
                                                   double tolerance) const
{
  if(!isOnSurface(point, tolerance))
  {
    return std::nullopt;
  }

  // Check which face the point lies on
  if(std::abs(point.x() - min_corner_.x()) <= tolerance)
    return BoxFace::XMin;
  if(std::abs(point.x() - max_corner_.x()) <= tolerance)
    return BoxFace::XMax;
  if(std::abs(point.y() - min_corner_.y()) <= tolerance)
    return BoxFace::YMin;
  if(std::abs(point.y() - max_corner_.y()) <= tolerance)
    return BoxFace::YMax;
  if(std::abs(point.z() - min_corner_.z()) <= tolerance)
    return BoxFace::ZMin;
  if(std::abs(point.z() - max_corner_.z()) <= tolerance)
    return BoxFace::ZMax;

  return std::nullopt;
}

// Geometric transformations
Box Box::translate(const Vector3D &offset) const
{
  return Box(min_corner_ + offset, max_corner_ + offset, material_, name_);
}

Box Box::scale(double factor) const
{
  Vector3D centre_point = centre();
  Vector3D new_dimensions = dimensions() * factor;
  Vector3D half_dims = new_dimensions * 0.5;

  return Box(centre_point - half_dims, centre_point + half_dims, material_,
             name_);
}

Box Box::scale(const Vector3D &factors) const
{
  Vector3D centre_point = centre();
  Vector3D current_dims = dimensions();
  Vector3D new_dimensions(current_dims.x() * factors.x(),
                          current_dims.y() * factors.y(),
                          current_dims.z() * factors.z());
  Vector3D half_dims = new_dimensions * 0.5;

  return Box(centre_point - half_dims, centre_point + half_dims, material_,
             name_);
}

Box Box::expand(double amount) const
{
  Vector3D expansion(amount, amount, amount);
  return Box(min_corner_ - expansion, max_corner_ + expansion, material_,
             name_);
}

Box Box::expand(const Vector3D &amounts) const
{
  return Box(min_corner_ - amounts, max_corner_ + amounts, material_, name_);
}

// Box-box operations
std::optional<Box> Box::intersection(const Box &other) const
{
  Vector3D new_min(std::max(min_corner_.x(), other.min_corner_.x()),
                   std::max(min_corner_.y(), other.min_corner_.y()),
                   std::max(min_corner_.z(), other.min_corner_.z()));

  Vector3D new_max(std::min(max_corner_.x(), other.max_corner_.x()),
                   std::min(max_corner_.y(), other.max_corner_.y()),
                   std::min(max_corner_.z(), other.max_corner_.z()));

  // Check if intersection is valid
  if(new_min.x() >= new_max.x() || new_min.y() >= new_max.y() ||
     new_min.z() >= new_max.z())
  {
    return std::nullopt;
  }

  return Box(new_min, new_max, material_, "Intersection");
}

Box Box::boundingBox(const Box &other) const
{
  Vector3D new_min(std::min(min_corner_.x(), other.min_corner_.x()),
                   std::min(min_corner_.y(), other.min_corner_.y()),
                   std::min(min_corner_.z(), other.min_corner_.z()));

  Vector3D new_max(std::max(max_corner_.x(), other.max_corner_.x()),
                   std::max(max_corner_.y(), other.max_corner_.y()),
                   std::max(max_corner_.z(), other.max_corner_.z()));

  return Box(new_min, new_max, Material::createAir(), "Bounding Box");
}

bool Box::intersects(const Box &other) const
{
  return min_corner_.x() <= other.max_corner_.x() &&
         max_corner_.x() >= other.min_corner_.x() &&
         min_corner_.y() <= other.max_corner_.y() &&
         max_corner_.y() >= other.min_corner_.y() &&
         min_corner_.z() <= other.max_corner_.z() &&
         max_corner_.z() >= other.min_corner_.z();
}

bool Box::contains(const Box &other) const
{
  return min_corner_.x() <= other.min_corner_.x() &&
         min_corner_.y() <= other.min_corner_.y() &&
         min_corner_.z() <= other.min_corner_.z() &&
         max_corner_.x() >= other.max_corner_.x() &&
         max_corner_.y() >= other.max_corner_.y() &&
         max_corner_.z() >= other.max_corner_.z();
}

// Point generation
Vector3D Box::randomPointInside() const
{
  auto &rng = RandomNumberGenerator::getThreadLocal();
  return Vector3D(rng.uniform(min_corner_.x(), max_corner_.x()),
                  rng.uniform(min_corner_.y(), max_corner_.y()),
                  rng.uniform(min_corner_.z(), max_corner_.z()));
}

Vector3D Box::randomPointOnSurface() const
{
  auto &rng = RandomNumberGenerator::getThreadLocal();

  // Choose a random face weighted by area
  Vector3D dim = dimensions();
  double area_xy = dim.x() * dim.y();
  double area_xz = dim.x() * dim.z();
  double area_yz = dim.y() * dim.z();
  double total_area = 2.0 * (area_xy + area_xz + area_yz);

  std::vector<double> face_weights = {area_yz, area_yz, area_xz,
                                      area_xz, area_xy, area_xy};

  std::size_t face_index = rng.discreteSample(face_weights);
  BoxFace face = static_cast<BoxFace>(face_index);

  Vector3D point;
  switch(face)
  {
  case BoxFace::XMin:
    point =
        Vector3D(min_corner_.x(), rng.uniform(min_corner_.y(), max_corner_.y()),
                 rng.uniform(min_corner_.z(), max_corner_.z()));
    break;
  case BoxFace::XMax:
    point =
        Vector3D(max_corner_.x(), rng.uniform(min_corner_.y(), max_corner_.y()),
                 rng.uniform(min_corner_.z(), max_corner_.z()));
    break;
  case BoxFace::YMin:
    point =
        Vector3D(rng.uniform(min_corner_.x(), max_corner_.x()), min_corner_.y(),
                 rng.uniform(min_corner_.z(), max_corner_.z()));
    break;
  case BoxFace::YMax:
    point =
        Vector3D(rng.uniform(min_corner_.x(), max_corner_.x()), max_corner_.y(),
                 rng.uniform(min_corner_.z(), max_corner_.z()));
    break;
  case BoxFace::ZMin:
    point = Vector3D(rng.uniform(min_corner_.x(), max_corner_.x()),
                     rng.uniform(min_corner_.y(), max_corner_.y()),
                     min_corner_.z());
    break;
  case BoxFace::ZMax:
    point = Vector3D(rng.uniform(min_corner_.x(), max_corner_.x()),
                     rng.uniform(min_corner_.y(), max_corner_.y()),
                     max_corner_.z());
    break;
  }

  return point;
}

std::vector<Vector3D> Box::gridPoints(int nx, int ny, int nz) const
{
  std::vector<Vector3D> points;
  points.reserve(nx * ny * nz);

  Vector3D dim = dimensions();

  for(int i = 0; i < nx; ++i)
  {
    for(int j = 0; j < ny; ++j)
    {
      for(int k = 0; k < nz; ++k)
      {
        double x = min_corner_.x() + (dim.x() * i) / (nx - 1);
        double y = min_corner_.y() + (dim.y() * j) / (ny - 1);
        double z = min_corner_.z() + (dim.z() * k) / (nz - 1);
        points.emplace_back(x, y, z);
      }
    }
  }

  return points;
}

// Validation
bool Box::hasValidDimensions() const
{
  return min_corner_.x() < max_corner_.x() &&
         min_corner_.y() < max_corner_.y() && min_corner_.z() < max_corner_.z();
}

std::vector<std::string> Box::validate() const
{
  std::vector<std::string> errors;

  if(!hasValidDimensions())
  {
    errors.push_back("Invalid dimensions: min corner must be less than max "
                     "corner in all axes");
  }

  if(!material_.isValid())
  {
    errors.push_back("Invalid material");
  }

  Vector3D dim = dimensions();
  if(dim.x() <= 0.0 || dim.y() <= 0.0 || dim.z() <= 0.0)
  {
    errors.push_back("Non-positive dimensions");
  }

  return errors;
}

// String representation
std::string Box::toString() const
{
  std::stringstream ss;
  ss << "Box '" << name_ << "'\n";
  ss << "  Min corner: (" << min_corner_.x() << ", " << min_corner_.y() << ", "
     << min_corner_.z() << ")\n";
  ss << "  Max corner: (" << max_corner_.x() << ", " << max_corner_.y() << ", "
     << max_corner_.z() << ")\n";
  ss << "  Dimensions: " << width() << " × " << height() << " × " << depth()
     << "\n";
  ss << "  Volume: " << volume() << "\n";
  ss << "  Surface area: " << surfaceArea() << "\n";
  ss << "  Material: " << material_.name()
     << " (density: " << material_.density() << " g/cm³)\n";
  ss << "  Valid: " << (is_valid_ ? "Yes" : "No");

  return ss.str();
}

std::string Box::getGeometryString() const
{
  std::stringstream ss;
  ss << std::fixed << std::setprecision(3);
  ss << "Box[(" << min_corner_.x() << "," << min_corner_.y() << ","
     << min_corner_.z() << ") to (" << max_corner_.x() << "," << max_corner_.y()
     << "," << max_corner_.z() << ")]";
  return ss.str();
}

// Comparison operators
bool Box::operator==(const Box &other) const
{
  constexpr double tolerance = 1e-10;
  return min_corner_ == other.min_corner_ && max_corner_ == other.max_corner_ &&
         material_ == other.material_ && name_ == other.name_;
}

bool Box::operator!=(const Box &other) const { return !(*this == other); }

// Static factory methods
Box Box::createCube(const Vector3D &centre, double side_length)
{
  double half_side = side_length * 0.5;
  Vector3D offset(half_side, half_side, half_side);
  return Box(centre - offset, centre + offset);
}

Box Box::createCube(const Vector3D &centre, double side_length,
                    const Material &material)
{
  double half_side = side_length * 0.5;
  Vector3D offset(half_side, half_side, half_side);
  return Box(centre - offset, centre + offset, material);
}

Box Box::createFromCentreAndDimensions(const Vector3D &centre,
                                       const Vector3D &dimensions)
{
  Vector3D half_dims = dimensions * 0.5;
  return Box(centre - half_dims, centre + half_dims);
}

Box Box::createFromCentreAndDimensions(const Vector3D &centre,
                                       const Vector3D &dimensions,
                                       const Material &material)
{
  Vector3D half_dims = dimensions * 0.5;
  return Box(centre - half_dims, centre + half_dims, material);
}

Box Box::createLeadShield(const Vector3D &min_corner,
                          const Vector3D &max_corner, double density)
{
  return Box(min_corner, max_corner, Material::createLead(density),
             "Lead Shield");
}

Box Box::createConcreteShield(const Vector3D &min_corner,
                              const Vector3D &max_corner, double density)
{
  return Box(min_corner, max_corner, Material::createConcrete(density),
             "Concrete Shield");
}

Box Box::createWaterPhantom(const Vector3D &min_corner,
                            const Vector3D &max_corner)
{
  return Box(min_corner, max_corner, Material::createWater(), "Water Phantom");
}

Box Box::createAirRegion(const Vector3D &min_corner, const Vector3D &max_corner)
{
  return Box(min_corner, max_corner, Material::createAir(), "Air Region");
}

// Private methods
void Box::invalidateCache() const
{
  centre_.reset();
  dimensions_.reset();
  volume_.reset();
}

void Box::validateAndSetState()
{
  auto errors = validate();
  is_valid_ = errors.empty();
}

RayBoxIntersection Box::computeIntersection(const Vector3D &ray_origin,
                                            const Vector3D &ray_direction,
                                            double t, BoxFace face,
                                            bool entering) const
{
  Vector3D intersection_point = ray_origin + ray_direction * t;
  Vector3D surface_normal = getSurfaceNormal(face);

  // Ensure normal points outward
  if(entering && surface_normal.dot(ray_direction) > 0)
  {
    surface_normal = surface_normal * -1.0;
  }
  else if(!entering && surface_normal.dot(ray_direction) < 0)
  {
    surface_normal = surface_normal * -1.0;
  }

  return RayBoxIntersection(t, intersection_point, face, surface_normal,
                            entering);
}

std::tuple<double, double> Box::intersectSlab(double ray_origin, double ray_dir,
                                              double slab_min,
                                              double slab_max) const
{
  if(std::abs(ray_dir) < 1e-15)
  {
    // Ray is parallel to slab
    if(ray_origin < slab_min || ray_origin > slab_max)
    {
      return {1.0, -1.0}; // Invalid interval
    }
    else
    {
      return {-std::numeric_limits<double>::infinity(),
              std::numeric_limits<double>::infinity()};
    }
  }

  double t1 = (slab_min - ray_origin) / ray_dir;
  double t2 = (slab_max - ray_origin) / ray_dir;

  if(t1 > t2)
  {
    std::swap(t1, t2);
  }

  return {t1, t2};
}

// Utility functions
std::string boxFaceToString(BoxFace face)
{
  switch(face)
  {
  case BoxFace::XMin:
    return "X-Min";
  case BoxFace::XMax:
    return "X-Max";
  case BoxFace::YMin:
    return "Y-Min";
  case BoxFace::YMax:
    return "Y-Max";
  case BoxFace::ZMin:
    return "Z-Min";
  case BoxFace::ZMax:
    return "Z-Max";
  default:
    return "Unknown";
  }
}

Vector3D getBoxFaceNormal(BoxFace face)
{
  switch(face)
  {
  case BoxFace::XMin:
    return Vector3D(-1.0, 0.0, 0.0);
  case BoxFace::XMax:
    return Vector3D(1.0, 0.0, 0.0);
  case BoxFace::YMin:
    return Vector3D(0.0, -1.0, 0.0);
  case BoxFace::YMax:
    return Vector3D(0.0, 1.0, 0.0);
  case BoxFace::ZMin:
    return Vector3D(0.0, 0.0, -1.0);
  case BoxFace::ZMax:
    return Vector3D(0.0, 0.0, 1.0);
  default:
    return Vector3D::ZERO;
  }
}

std::optional<BoxFace> stringToBoxFace(const std::string &face_str)
{
  std::string lower_str = face_str;
  std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(),
                 ::tolower);

  if(lower_str == "xmin" || lower_str == "x-min")
    return BoxFace::XMin;
  if(lower_str == "xmax" || lower_str == "x-max")
    return BoxFace::XMax;
  if(lower_str == "ymin" || lower_str == "y-min")
    return BoxFace::YMin;
  if(lower_str == "ymax" || lower_str == "y-max")
    return BoxFace::YMax;
  if(lower_str == "zmin" || lower_str == "z-min")
    return BoxFace::ZMin;
  if(lower_str == "zmax" || lower_str == "z-max")
    return BoxFace::ZMax;

  return std::nullopt;
}

bool boxesIntersect(const Box &box1, const Box &box2)
{
  return box1.intersects(box2);
}

std::optional<Box> intersectBoxes(const Box &box1, const Box &box2)
{
  return box1.intersection(box2);
}

Box boundingBoxOf(const std::vector<Box> &boxes)
{
  if(boxes.empty())
  {
    throw std::invalid_argument(
        "Cannot compute bounding box of empty box list");
  }

  Vector3D min_corner = boxes[0].minCorner();
  Vector3D max_corner = boxes[0].maxCorner();

  for(std::size_t i = 1; i < boxes.size(); ++i)
  {
    const Vector3D &box_min = boxes[i].minCorner();
    const Vector3D &box_max = boxes[i].maxCorner();

    min_corner = Vector3D(std::min(min_corner.x(), box_min.x()),
                          std::min(min_corner.y(), box_min.y()),
                          std::min(min_corner.z(), box_min.z()));

    max_corner = Vector3D(std::max(max_corner.x(), box_max.x()),
                          std::max(max_corner.y(), box_max.y()),
                          std::max(max_corner.z(), box_max.z()));
  }

  return Box(min_corner, max_corner, Material::createAir(), "Bounding Box");
}