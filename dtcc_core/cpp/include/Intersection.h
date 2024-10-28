//
// Created by Dag Wästberg on 2024-09-30.
//

#ifndef DTCC_BUILDER_INTERSECTION_H
#define DTCC_BUILDER_INTERSECTION_H

#include <limits>
#include <cmath>

#include "Geometry.h"
#include "Logging.h"
#include "MeshBuilder.h"

#include "model/Mesh.h"
#include "model/Vector.h"
#include "model/Surface.h"
namespace DTCC_BUILDER
{
    class Intersection
    {
    public:
        static Vector3D ray_triangle_intersection(const std::vector<Vector3D> &triangle,
                                                  const Vector3D &ray_origin,
                                                  const Vector3D &ray_vector)
        {
            constexpr float epsilon = std::numeric_limits<float>::epsilon();
            const Vector3D MISS = Vector3D(NAN, NAN, NAN);
            Vector3D edge1 = triangle[1] - triangle[0];
            Vector3D edge2 = triangle[2] - triangle[0];
            Vector3D ray_cross_e2 = Geometry::cross_3d(ray_vector, edge2);
            auto det = static_cast<float>(edge1.dot(ray_cross_e2));

            if (fabs(det) < epsilon)
                return MISS; // This ray is parallel to this triangle.

            auto inv_det = static_cast<float>(1.0 / det);
            auto s = ray_origin - triangle[0];
            float u = inv_det * static_cast<float>(s.dot(ray_cross_e2));

            if (u < 0 || u > 1)
                return MISS;

            auto s_cross_e1 = Geometry::cross_3d(s, edge1);
            float v = inv_det * static_cast<float>(ray_vector.dot(s_cross_e1));

            if (v < 0 || u + v > 1)
                return MISS;

            // At this stage we can compute t to find out where the intersection point is on the line.
            float t = inv_det * static_cast<float>(edge2.dot(s_cross_e1));

            if (t > epsilon) // ray intersection
            {
                return (ray_origin + ray_vector * t);
            }
            else // This means that there is a line intersection but not a ray intersection.
                return MISS;
        }

        static Vector3D ray_surface_intersection(const Surface &surface,
                                                 const Vector3D &ray_origin,
                                                 const Vector3D &ray_vector)
        {

            const Vector3D MISS = Vector3D(NAN, NAN, NAN);
            if (surface.vertices.size() < 3)
            {
                return MISS;
            }
            if (surface.vertices.size() == 3 && surface.holes.size() == 0)
            {

                return ray_triangle_intersection(surface.vertices, ray_origin, ray_vector);
            }

            Mesh surface_mesh = MeshBuilder::mesh_surface(surface);
            for (const auto &t : surface_mesh.triangles())
            {
                Vector3D intersection = ray_triangle_intersection(t, ray_origin, ray_vector);
                if (!std::isnan(intersection[0]))
                {
                    return intersection;
                }
            }
            return MISS;
        }

        static Vector3D ray_multisurface_intersection(const MultiSurface &multisurface,
                                                      const Vector3D &ray_origin,
                                                      const Vector3D &ray_vector)
        {
            Vector3D closest_hit = Vector3D(NAN, NAN, NAN);
            double closest_distance = std::numeric_limits<double>::max();
            for (const auto &s : multisurface.surfaces)
            {
                Vector3D intersection = ray_surface_intersection(s, ray_origin, ray_vector);
                if (std::isnan(intersection[0]))
                {
                    continue;
                }
                double distance = (intersection - ray_origin).magnitude();
                if (distance < closest_distance)
                {
                    closest_distance = distance;
                    closest_hit = intersection;
                }
            }
            return closest_hit;
        }
    };
}

#endif // DTCC_BUILDER_INTERSECTION_H
