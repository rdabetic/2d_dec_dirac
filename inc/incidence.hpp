#ifndef INCIDENCE_H
#define INCIDENCE_H

#include "mfem.hpp"
#include <eigen3/Eigen/Eigen>

/// Extract the Vertex-Edge incidence matrix
/// Assume edges are oriented according to vertex numbers
Eigen::SparseMatrix<double> VertexEdgeIncidence(mfem::Mesh& mesh);

/// Extract the Vertex-Edge incidence matrix
/// Assume edges are oriented according to vertex numbers
Eigen::SparseMatrix<double> EdgeFaceIncidence(mfem::Mesh& mesh);

#endif
