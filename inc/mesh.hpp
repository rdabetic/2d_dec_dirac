#ifndef MESH_H
#define MESH_H

#include "mfem.hpp"
#include <eigen3/Eigen/Eigen>

/// Custom mesh class that inherits from mfem::Mesh
/// We need more functionality than they do (or I just couldn't find it)
struct Mesh2D: public mfem::Mesh {
  using mfem::Mesh::vertices;
  // Incidence matrices
  Eigen::SparseMatrix<double> vertex_edge, edge_face;
  Eigen::VectorXd triangle_areas;
  std::vector<bool> ess_v, ess_e;
  unsigned int n_ess_v, n_ess_e;

  /// Constructor, call with a move, takes over the mfem::Mesh
  Mesh2D(mfem::Mesh&& mesh);
  /// Update mesh attributes that come from the derived class
  void update();
  /// Refine using UniformRefinement
  void refine();

  /// Get element center
  Eigen::Vector2d center(int elem_id);
  Eigen::Vector2d barycenter(int elem_id);
  /// Get vertex coordinates as a vector with Eigen
  Eigen::Vector2d vert2vec(int vid) const;
  /// Check if all the triangles are acute
  void assertCircumcentric() const;
  /// Get the mesh-width
  double meshWidth() const;
};

#endif
