using StaticArrays;

abstract type SceneObject end

immutable Ray
  origin::SVector{3, Float64}
  direction::SVector{3, Float64}
  Ray(a, b) = new(a, normalize(b - a))
end

immutable Triangle <: SceneObject
  a::SVector{3, Float64} # point
  u::SVector{3, Float64} # basis vector
  v::SVector{3, Float64} # basis vector
  n::SVector{3, Float64} # normal vector

  function Triangle(a::SVector{3, Float64}, b::SVector{3, Float64}, c::SVector{3, Float64})
    u = b - a
    v = c - a
    n = cross(u, v)
    new(a, u, v, normalize(n))
  end
end

immutable Intersection
    point::SVector{3, Float64}
    distance::Float64
    is_hit::Bool
end

const no_intersection = Intersection((0., 0., 0.), 0., false)

function intersect_ray(ray::Ray, obj::Triangle)
  M = det(@SMatrix [-obj.u[1] -obj.v[1] ray.direction[1];
                    -obj.u[2] -obj.v[2] ray.direction[2];
                    -obj.u[3] -obj.v[3] ray.direction[3]])

  M == 0.0 && return no_intersection

  β = det(@SMatrix [(obj.a[1] - ray.origin[1]) -obj.v[1] ray.direction[1];
                    (obj.a[2] - ray.origin[2]) -obj.v[2] ray.direction[2];
                    (obj.a[3] - ray.origin[3]) -obj.v[3] ray.direction[3]])/M

  (β < 0.) && return no_intersection

  γ = det(@SMatrix [-obj.u[1] (obj.a[1] - ray.origin[1]) ray.direction[1];
                    -obj.u[2] (obj.a[2] - ray.origin[2]) ray.direction[2];
                    -obj.u[3] (obj.a[3] - ray.origin[3]) ray.direction[3]])/M

  (γ < 0.) && return no_intersection

  t = det(@SMatrix [-obj.u[1] -obj.v[1] (obj.a[1] - ray.origin[1]);
                    -obj.u[2] -obj.v[2] (obj.a[2] - ray.origin[2]);
                    -obj.u[3] -obj.v[3] (obj.a[3] - ray.origin[3])])/M

  ((β + γ) > 1.) && return no_intersection
  t > 0. && return no_intersection
  return Intersection(tri.a + β*tri.u + γ*tri.v, t, true)
end

function trace_ray(ray::Ray, value::Float64, direction::Symbol=:z)
  const idx_map = Dict{Symbol, Int64}(:x => 1, :y => 2, :z => 3)
  idx = idx_map[direction]
  (ray.direction[idx] == 0.0) && return no_intersection

  t::Float64 = (value - ray.origin[idx])/ray.direction[idx]
  t < 0. && return no_intersection
  return Intersection(ray.origin + t*ray.direction, t, true)
end

#tri = Triangle(SVector(-1., 0., 0.), SVector(1., 0., 0.), SVector(0., 1., 0.))

#ray = Ray(SVector(0., 0., 0.), SVector(0., 0., -1))

#println(intersect_ray(ray, tri))
