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
  template <bool MIDPOINT = false,
            typename FUNC0,
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
    mfem::Array<int> vert;
    if constexpr(MIDPOINT) {
      // Edges
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
    } else {
      /* Gauss Quadrature for the edges */
      const Eigen::Vector3d 
        gauss_pts({
          -std::sqrt(3. / 5),
          0.,
          std::sqrt(3. / 5),
        }),
        gauss_weights({
            5. / 9.,
            8. / 9.,
            5. / 9.
        });
      auto gauss_quad = 
        [&]<typename F> (F f, const double& a, const double& b) -> double {
          double s = 0;
          for(unsigned int k = 0; k < 3; ++k)
            s += f((b - a) / 2. * gauss_pts(k) + (a + b) / 2.) * 
                  gauss_weights(k);
          //const Eigen::Vector3d f_eval = 
          //  Eigen::Vector3d::NullaryExpr(
          //      [&](const unsigned int i) -> double {
          //        return f((b - a) / 2. * gauss_pts(i) + (a + b) / 2.);
          //      }
          //  );
          //return (b - a) / 2. * gauss_weights.dot(f_eval);
          return s * (b - a) / 2.;
        }
      ;
      // Edges
      for(unsigned int k = 0; k < ne; ++k) {
        mesh.GetEdgeVertices(k, vert);
        const Eigen::Vector2d 
          v0 = mesh.vert2vec(vert[0]),
          v1 = mesh.vert2vec(vert[1]);
        auto integrand = [&](const double& s) -> double {
          return f1(v0 + s * (v1 - v0)).dot(v1 - v0);
        };
        c.u1(k) = gauss_quad(integrand, 0, 1); 
      }
      /* Quadrature for the faces 
       * See https://github.com/pratyushpotu/DEC_convergence_tests/blob/main/dec/integrators.py
       * */
      Eigen::Vector<double, 7> tria_weights;
      tria_weights << 
        0.225,
        0.132394152788, 0.132394152788, 0.132394152788,
        0.125939180544, 0.125939180544, 0.125939180544;
      Eigen::Matrix<double, 3, 7> tria_pts_bary;
      {
        Eigen::Matrix<double, 7, 3> tria_pts_bary_;
        tria_pts_bary_ << 
          1./3., 1./3., 1./3.,
          0.0597158717898, 0.470142064105, 0.470142064105,
          0.470142064105, 0.0597158717898, 0.470142064105,
          0.470142064105, 0.470142064105, 0.0597158717898,
          0.797426985353, 0.101286507323, 0.101286507323,
          0.101286507323, 0.797426985353, 0.101286507323,
          0.101286507323, 0.101286507323, 0.797426985353;
        // Tranpose for in-stride access
        tria_pts_bary = tria_pts_bary_.transpose();
      }
      auto triaArea = [&] (const Eigen::Vector2d& A, 
                           const Eigen::Vector2d& B,
                           const Eigen::Vector2d& C) -> double {
        const auto AB = B - A,
                   AC = C - A;
        return 
          std::abs(
            AB(0) * AC(1) - AB(1) * AC(0)
          ) / 2.;
      };
      auto tria_quad = 
        [&]<typename F> (F f, const Eigen::Vector2d& A, 
                              const Eigen::Vector2d& B,
                              const Eigen::Vector2d& C) -> double {
          double s = 0;
          for(unsigned int k = 0; k < 7; ++k) {
            const Eigen::Vector2d p = 
              tria_pts_bary(0, k) * A + 
              tria_pts_bary(1, k) * B + 
              tria_pts_bary(2, k) * C;
            s += tria_weights(k) * f(p);
          }
          // Triangle area
          return s * triaArea(A, B, C);
        }
      ;
      // Triangles
      for(unsigned int k = 0; k < nt; ++k) {
        mesh.GetElementVertices(k, vert);
        const Eigen::Vector2d A = mesh.vert2vec(vert[0]),
                              B = mesh.vert2vec(vert[1]),
                              C = mesh.vert2vec(vert[2]);
        c.u2(k) = tria_quad(f2, A, B, C);
      }
    }

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
  double sobolevSemiNorm(const Cochain& u) const;
  double sobolevNorm(const Cochain& u) const;
};

#endif
