#ifndef DTCC_VOLUME_MESH_IMPROVEMENT_H
#define DTCC_VOLUME_MESH_IMPROVEMENT_H

#include "Geometry.h"
#include "../model/VolumeMesh.h"
#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

namespace DTCC_BUILDER
{

class VolumeMeshImprovement
{
private:

  static std::vector<std::pair<size_t,size_t>> get_candidate_vertices_pairs(const VolumeMesh &volume_mesh,
                                                                          const std::vector<size_t> tetrahedra_to_remove){
    // Ideally, each vertex should merge into itself
    std::vector<std::pair<size_t, size_t>> vertex_pairs;

    for (auto c : tetrahedra_to_remove)
    {
      // Get the smallest edge of the tetrahedron
      std::array<size_t, 4> vertex_indices = {volume_mesh.cells[c].v0, volume_mesh.cells[c].v1,
                                              volume_mesh.cells[c].v2, volume_mesh.cells[c].v3};

      const std::array<Vector3D, 4> vertices = {volume_mesh.vertices[volume_mesh.cells[c].v0],
                                          volume_mesh.vertices[volume_mesh.cells[c].v1],
                                          volume_mesh.vertices[volume_mesh.cells[c].v2],
                                          volume_mesh.vertices[volume_mesh.cells[c].v3]};

      std::pair<size_t, size_t> smallest_edge = get_smallest_edge(vertices);

      // Merge the two vertices of the smallest edge

      auto rv1 = vertex_indices[smallest_edge.first];
      auto rv2 = vertex_indices[smallest_edge.second];

      // rv2 should always have marker greater or equal to rv1
      if (volume_mesh.markers[rv1] > volume_mesh.markers[rv2])
        std::swap(rv1, rv2);
      if (!(volume_mesh.markers[rv1] > -3))
      {
        vertex_pairs.push_back({rv1, rv2});
      }

    }

    for (auto &p : vertex_pairs)
    {
      if (p.second < p.first)
      {
        std::swap(p.first, p.second);
      }
    }

  
    // Sort the vector so that duplicate pairs are adjacent.
    std::sort(vertex_pairs.begin(), vertex_pairs.end());

    // Remove duplicates.
    auto last = std::unique(vertex_pairs.begin(), vertex_pairs.end());
    vertex_pairs.erase(last, vertex_pairs.end());

    return vertex_pairs;
  }

  // Helper function: computes the union-find mapping given the merging pairs.
  static std::vector<size_t> computeUnionFindMapping(
      const VolumeMesh &volume_mesh,
      const std::vector<std::pair<size_t, size_t>> &vertexPairs)
  {
    size_t n = volume_mesh.vertices.size();
    std::vector<size_t> parent(n);
    std::vector<size_t> rank(n, 0);

    // Initialize each vertex as its own parent.
    for (size_t i = 0; i < n; ++i)
    {
      parent[i] = i;
    }

    // Lambda for find with path compression.
    std::function<size_t(size_t)> find = [&](size_t x) -> size_t {
      if (parent[x] != x)
        parent[x] = find(parent[x]);
      return parent[x];
    };

    // Union function that uses both rank and marker criteria.
    // auto unionSets = [&](size_t x, size_t y) {
    //   size_t rootX = find(x);
    //   size_t rootY = find(y);
    //   if (rootX == rootY)
    //     return;

    //   // Use marker preference: assume higher marker means we want that vertex to be the representative.
    //   if (volume_mesh.markers[rootX] < volume_mesh.markers[rootY])
    //     std::swap(rootX, rootY);

    //   // Then do union by rank.
    //   if (rank[rootX] < rank[rootY])
    //   {
    //     parent[rootX] = rootY;
    //   }
    //   else if (rank[rootX] > rank[rootY])
    //   {
    //     parent[rootY] = rootX;
    //   }
    //   else
    //   {
    //     parent[rootY] = rootX;
    //     rank[rootX]++;
    //   }
    // };

    auto union_sets = [&](size_t i, size_t j)
    {
      size_t root_i = find(i);
      size_t root_j = find(j);
      if (root_i == root_j)
        return;
      // Prefer the vertex with the higher marker (assuming positive markers indicate domain
      // boundaries). Adjust the comparison as needed based on your marker convention.
      if (volume_mesh.markers[root_i] < volume_mesh.markers[root_j])
        std::swap(root_i, root_j);
      parent[root_j] = root_i;
    };

    // Process each merging pair.
    for (const auto &p : vertexPairs)
    {
      union_sets(p.first, p.second);
    }

    // Build and return the mapping: for each vertex, find its representative.
    std::vector<size_t> vertexMap(n);
    for (size_t i = 0; i < n; ++i)
    {
      vertexMap[i] = find(i);
    }
    return vertexMap;
  }

  static std::pair<size_t, size_t> get_smallest_edge(const std::array<Vector3D, 4> &vertices)
  {
    std::pair<size_t, size_t> best_pair = {0, 1};
    double min_distance = (vertices[1] - vertices[0]).magnitude(); // initial edge distance

    // Compare all unique pairs (i, j) with i < j
    for (size_t i = 0; i < vertices.size(); ++i)
    {
      for (size_t j = i + 1; j < vertices.size(); ++j)
      {
        double _distance = (vertices[j] - vertices[i]).magnitude();
        if (_distance < min_distance)
        {
          min_distance = _distance;
          best_pair = {i, j};
        }
      }
    }
    return best_pair;
  }

  static size_t count_unique(size_t a, size_t b, size_t c, size_t d) {
    std::array<size_t, 4> arr = {a, b, c, d};
    std::sort(arr.begin(), arr.end());
    auto unique_end = std::unique(arr.begin(), arr.end());
    return std::distance(arr.begin(), unique_end);
  }

  static void _update_mesh(VolumeMesh &volume_mesh){
    info("Updating mesh vertices");
    std::vector<bool> keep_vertex(volume_mesh.vertices.size(),false);
    for (const auto cell : volume_mesh.cells){
      keep_vertex[cell.v0] = true;
      keep_vertex[cell.v1] = true;
      keep_vertex[cell.v2] = true;
      keep_vertex[cell.v3] = true;
    }

    std::vector<size_t> new_index(volume_mesh.vertices.size(), 0);
    std::vector<Vector3D> new_vertices;
    std::vector<int> new_markers;
    size_t counter = 0;
    for (size_t i = 0; i < volume_mesh.vertices.size(); i++)
    {
      if (keep_vertex[i])
      {
        new_index[i] = counter;
        new_vertices.push_back(volume_mesh.vertices[i]);
        new_markers.push_back(volume_mesh.markers[i]);
        ++counter;
      }
    }

    for (auto &cell : volume_mesh.cells)
    {
      cell.v0 = new_index[cell.v0];
      cell.v1 = new_index[cell.v1];
      cell.v2 = new_index[cell.v2];
      cell.v3 = new_index[cell.v3];
    }

    volume_mesh.vertices = new_vertices;
    volume_mesh.markers = new_markers;
  }



public:
  static VolumeMesh remove_tetrahedra(const VolumeMesh &volume_mesh,
                                const double aspect_ratio_threshold = 10.0)
  {

    const auto _aspect_ratios = Geometry::aspect_ratios(volume_mesh);

    std::vector<size_t> _tetrahedra_to_remove;
    for (size_t i = 0; i < volume_mesh.cells.size(); i++)
    {
      if (_aspect_ratios[i] > aspect_ratio_threshold)
      {
        _tetrahedra_to_remove.push_back(i);
      }
    }
    // Sort the tetrahedra to remove based on their aspect ratios
    // This is optional, but it can help in prioritizing which tetrahedra to remove first.
    std::sort(_tetrahedra_to_remove.begin(), _tetrahedra_to_remove.end(),
          [&_aspect_ratios](size_t a, size_t b) {
              return _aspect_ratios[a] > _aspect_ratios[b]; // descending order
          });

    std::cout << "Number of tetrahedra to remove: " << _tetrahedra_to_remove.size() << std::endl;
  
    std::vector<std::pair<size_t,size_t>> vertex_pairs = get_candidate_vertices_pairs(volume_mesh,_tetrahedra_to_remove);
    
    std::vector<size_t> _vertex_map = computeUnionFindMapping(volume_mesh,vertex_pairs);

    // Now _vertex_map contains the mapping from old vertex indices to their new representative.
    // The next step would be to update the connectivity of your tetrahedral cells:
    // For each cell, replace its vertex indices with _vertex_map[vertex_index].
    // Also, you might need to remove any tetrahedra that become degenerate
    // (i.e., if after mapping two or more of its vertices are the same).

    // === End Union-Find ===

    // At this point, you can continue by rebuilding the mesh:
    //   - Create a new vertex list that contains only the unique representatives.
    //   - Update each cell’s connectivity using _vertex_map.
    //   - Remove cells that are now degenerate (e.g., have duplicate vertex indices).
    VolumeMesh _refined_volume_mesh; 
    _refined_volume_mesh.vertices = volume_mesh.vertices;
    _refined_volume_mesh.markers = volume_mesh.markers;
    
    for (size_t i = 0; i < volume_mesh.cells.size(); i++)
    {
       size_t v0 =_vertex_map[volume_mesh.cells[i].v0];
       size_t v1 =_vertex_map[volume_mesh.cells[i].v1];
       size_t v2 =_vertex_map[volume_mesh.cells[i].v2];
       size_t v3 =_vertex_map[volume_mesh.cells[i].v3];

      if (count_unique(v0,v1,v2,v3) == 4)
      {
        _refined_volume_mesh.cells.push_back(Simplex3D(v0,v1,v2,v3));
      }
    }
    
    // Erase non-connected vertices and renumber mesh
    _update_mesh(_refined_volume_mesh);

    return _refined_volume_mesh;
  }
};

} // namespace DTCC_BUILDER

#endif