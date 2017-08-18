using CoordinateTransformations, Rotations, StaticArrays;
using IterTools;
include("raytrace-demo.jl")
# -------------------------------------------------
# Setup world coordinates, location of scene objects and camera

polar_t = SphericalFromCartesian()
cartesian_t = CartesianFromSpherical()

# detector config
detector_npixels = SVector(71, 71) # number of pixels in each dimension (x, y)
detector_pixel_size = 0.0001 # single pixel size along edge in mm
detector_size = detector_npixels*detector_pixel_size

detector_to_camera_t =
    CylindricalFromCartesian() ∘
    LinearMap(SMatrix{3, 3}(1, 0, 0, 0, -1, 0, 0, 0, 1)) ∘
    Translation(SVector(-0.5detector_size[1] - 0.5*detector_pixel_size,
                        -0.5detector_size[2] - 0.5*detector_pixel_size, 0.)) ∘
    LinearMap(detector_pixel_size*SMatrix{3, 3}(1, 0, 0, 0, 1, 0, 0, 0, 1))

# Camera to world transformation
const camera_origin = SVector(0., 0., 0.)
const camera_t = inv(
    Translation(camera_origin) ∘
    LinearMap(RotZXY(0., 0., 0.))
    )
# Transformation from detector cyclindrical coordinates to camera spherical coordinates
function camera_model(coords::Cylindrical, f=0.3, f0=1., a1=0., a2=0.)
  r = coords.r/f0
  ϕ = 2.*atan(0.5f0/f*(r + a1*r*r*r + a2*r*r*r*r*r)) + π/2
  return Spherical(1., -coords.θ, ϕ)
end

# -------------------------------------------------
# For each ray, find the intersection with the plane and obtain texture at intersection point

tri = Triangle(SVector(-0.5, -0.5*tan(π/6), -0.1),
         SVector(0.5, -0.5*tan(π/6), -0.1),
         SVector(0., sqrt(3/4) - 0.5*tan(π/6), -0.1))

#ray = Ray(SVector(0., 0., 1.), SVector(0., 0., 0.))

# Plane in scene transformation
grid_t = inv(Translation(SVector(0., 0., 1)) ∘ LinearMap(RotZXY(0., 0., 0.)))

# -------------------------------------------------
# Generate rays from the detector in pixel coordinates, ray trace them from
# the detector plane to the scene

make_ray(x, y) = (camera_t ∘ cartesian_t ∘ camera_model ∘ detector_to_camera_t)(SVector(x, y, 0))
pt = map(v->SVector(v...), product(1:detector_npixels[1], 1:detector_npixels[2], 0))
camera_pt = map(CartesianFromCylindrical() ∘ detector_to_camera_t, pt)

transformed_pt = map(CartesianFromSpherical() ∘ camera_model ∘ detector_to_camera_t, pt)
trays = map(d -> Ray(camera_origin, camera_origin + d), transformed_pt)
pts_at_1 = map(i -> [i.point[1:2], i.is_hit ? "g" : "r"], map(r -> trace_ray(r, 1., :z), trays))

z = map(vec -> Ray(camera_origin, camera_origin + vec), make_ray(v[1], v[2]) for v in pt)
rays = reshape(z, prod(size(z)))

#points = []
#foreach(rays) do ray
#  t = 1/ray.direction[3]
#  v = t*ray.direction + ray.origin
#  res = intersect_ray(ray, tri)
#  println("$(ray.direction) → $(res) → $v")
#  if res.is_hit
#      push!(points, SVector(res.point[1], res.point[2]))
#  end
#end

using PyPlot
using PyCall
@pyimport matplotlib.patches as patch

a, b, c = tri.a[1:2], (tri.a + tri.u)[1:2], (tri.a + tri.v)[1:2]
tpatch = patch.Polygon([a, b, c], true, fill=false)
ax = gca()
ax[:add_artist](tpatch)
x, y, c = map(p -> p[1][1], pts_at_1), map(p -> p[1][2], pts_at_1), map(p -> p[2], pts_at_1)
dx, dy = map(p->p[1], camera_pt), map(p -> p[2], camera_pt)
scatter(x, y, color=c, s=2)
scatter(dx, dy, s=1.)
#xlim(-0.6, 0.6)
#ylim(-0.3, 0.7)
ax[:set_aspect]("equal")
