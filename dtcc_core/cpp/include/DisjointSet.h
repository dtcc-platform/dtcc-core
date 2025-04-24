#include <set>
#include <unordered_map>
#include <vector>

#ifndef DTCC_DISJOINSET_H
#define DTCC_DISJOINSET_H

namespace DTCC_BUILDER
{

class DisjointSet
{
private:
  std::vector<int> parent;
  std::vector<int> rank;
  int n;

public:
  // Constructor
  DisjointSet(int size) : n(size)
  {
    parent.resize(size);
    rank.resize(size, 0);

    // Initialize each element as a separate set
    for (int i = 0; i < size; i++)
    {
      parent[i] = i;
    }
  }

  // Find the representative (root) of the set containing x
  int find(int x)
  {
    // Path compression: make every examined node point directly to the root
    if (parent[x] != x)
    {
      parent[x] = find(parent[x]);
    }
    return parent[x];
  }

  // Merge the sets containing x and y
  void unionSets(int x, int y)
  {
    int rootX = find(x);
    int rootY = find(y);

    // If x and y are already in the same set
    if (rootX == rootY)
    {
      return;
    }

    // Union by rank: attach smaller rank tree under root of higher rank tree
    if (rank[rootX] < rank[rootY])
    {
      parent[rootX] = rootY;
    }
    else if (rank[rootX] > rank[rootY])
    {
      parent[rootY] = rootX;
    }
    else
    {
      // If ranks are equal, make one as root and increment its rank
      parent[rootY] = rootX;
      rank[rootX]++;
    }
  }

  // Get all sets as a map where keys are representatives and values are vectors of elements
  std::unordered_map<int, std::vector<int>> getSets()
  {
    // Make sure all paths are compressed
    for (int i = 0; i < n; i++)
    {
      find(i);
    }

    // Group elements by their representatives
    std::unordered_map<int, std::vector<int>> sets;
    for (int i = 0; i < n; i++)
    {
      int root = parent[i];
      sets[root].push_back(i);
    }

    return sets;
  }
};
} // namespace DTCC_BUILDER
#endif // DTCC_DISJOINSET_H