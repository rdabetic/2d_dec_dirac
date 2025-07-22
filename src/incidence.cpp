#include "incidence.hpp"
#include "bd.hpp"
#include <vector>

Eigen::SparseMatrix<double> VertexEdgeIncidence(mfem::Mesh& mesh) {
  // Number of edges and vertices
  const int nv = mesh.GetNV(),
            ne = mesh.GetNEdges();

  std::vector<Eigen::Triplet<double>> triplets; 
  triplets.reserve(2 * nv);
  // Store the vertices 
  mfem::Array<int> vert(2);
  for(int ei = 0; ei < ne; ++ei){
    // Get vertices of edge `ei`
    mesh.GetEdgeVertices(ei, vert);
    const int sign = (vert[0] < vert[1] ? 1 : -1);
    // Add to matrix
    triplets.emplace_back(ei, vert[0], -sign);
    triplets.emplace_back(ei, vert[1], sign);
  }
  // Tell Eigen to assemble the matrix and return
  Eigen::SparseMatrix<double> D(ne, nv);
  D.setFromTriplets(triplets.begin(), triplets.end());
  D.makeCompressed();
  return D;
}

Eigen::SparseMatrix<double> EdgeFaceIncidence(mfem::Mesh& mesh) {
  // Number of edges and faces (triangles)
  const int ne = mesh.GetNEdges(),
            nf = mesh.GetNE();
  std::vector<Eigen::Triplet<double>> triplets; 
  triplets.reserve(3 * ne);

  mfem::Array<int> edges(3), cor(3);
  for(int fi = 0; fi < nf; ++fi){
    // Get edges and orientations of element fi
    mesh.GetElementEdges(fi, edges, cor);
    // Add to matrix
    for(unsigned int k = 0; k < 3; ++k)
      triplets.emplace_back(fi, edges[k], cor[k]);
  }
  // Tell Eigen to assemble the matrix and return
  Eigen::SparseMatrix<double> D(nf, ne);
  D.setFromTriplets(triplets.begin(), triplets.end());
  D.makeCompressed();
  return D;
}
