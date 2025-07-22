#ifndef BD_H
#define BD_H

#include <vector>
#include "mfem.hpp"

/// Extract
/// - Locations (idx) of boundary points
/// - Points on the boundary (vector of true and false)
/// - No. points on the boundary
/// from the mesh
std::pair<std::vector<bool>,
          unsigned int>
extractBDPts(mfem::Mesh& mesh);

/// Extract
/// - Locations (idx) of boundary edges
/// - Edges on the boundary (also true false)
/// - No. edges on the boundary
/// from the mesh
std::pair<std::vector<bool>,
          unsigned int>
extractBDEdges(mfem::Mesh& mesh);

#endif
