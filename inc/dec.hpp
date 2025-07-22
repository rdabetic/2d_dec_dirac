#ifndef DEC_H
#define DEC_H

#include "mfem.hpp"
#include "mesh.hpp"
#include "cochain.hpp"
#include <eigen3/Eigen/Eigen>

struct DEC {
  // We do not copy
  Mesh2D& mesh;
  // Exterior derivatives, basically refs to mesh incidence matrices
  const Eigen::SparseMatrix<double>& d0, d1;
  // Codifferentials
  Eigen::SparseMatrix<double>  delta1, delta2;
  Eigen::VectorXd star0, star1, star2,
                  star0_inv, star1_inv, star2_inv;

  DEC(Mesh2D& mesh_);

  Cochain createCochain() const;
  Cochain getRandomCochain() const;
  void assembleHodgeStars();
  void applyHodgeStar(Cochain& u) const;
  void applyInvHodgeStar(Cochain& u) const;
  /// Set the entries in the vectors on the boundary to zero
  void zeroCochainBoundary(Cochain& c) const;
  // Apply the exterior derivative
  void diff(const Cochain& src,
                  Cochain& dest) const;

  void codiff(const Cochain& src,
                    Cochain& dest) const;

  /// deRham map using midpoint rule
  /// Functors have to evaluate like
  /// - f0: Eigen::Vector2d -> double
  /// - f1: Eigen::Vector2d -> Eigen::Vector2d
  /// - f2: Eigen::Vector2d -> double
  template <typename FUNC0,
            typename FUNC1,
            typename FUNC2>
  Cochain deRham(FUNC0 f0, FUNC1 f1, FUNC2 f2) const {
    // Needs to be implemented here, as it is a template
    const unsigned int nv = mesh.GetNV(),
                       ne = mesh.GetNEdges(),
                       nt = mesh.GetNE();
    // Result
    Cochain c(nv, ne, nt);
    // Vertices, point evaluation, no quadrature necessary
    for(unsigned int k = 0; k < nv; ++k)
      c.u0(k) = f0(mesh.vert2vec(k));
    // Edges
    mfem::Array<int> vert;
    for(unsigned int k = 0; k < ne; ++k) {
      mesh.GetEdgeVertices(k, vert);
      const Eigen::Vector2d 
        v0 = mesh.vert2vec(vert[0]),
        v1 = mesh.vert2vec(vert[1]),
        midpoint = (v0 + v1) / 2.;
      const double sign = (vert[0] < vert[1] ? 1 : -1);
      c.u1(k) = f1(midpoint).dot(v1 - v0) * sign;
    }
    // Triangles
    for(unsigned int k = 0; k < nt; ++k)
      c.u2(k) = mesh.triangle_areas(k) * 
                       f2(mesh.barycenter(k));

    // Essential boundary conditions and zero mean
    zeroCochainBoundary(c);
    return c;  
  }

  inline double inner(const Cochain& u,
                      const Cochain& v) const {
    return u.u0.dot((star0.array() * v.u0.array()).matrix()) + 
           u.u1.dot((star1.array() * v.u1.array()).matrix()) + 
           u.u2.dot((star2.array() * v.u2.array()).matrix());
  }

  inline double norm(const Cochain& u) const {
    return std::sqrt(inner(u, u));
  }

  double sobolevInner(const Cochain& u,
                      const Cochain& v) const;
  double sobolevNorm(const Cochain& u) const;
};

#endif
