#ifndef DIRAC_H
#define DIRAC_H

#include "dec.hpp"
#include "cochain.hpp"
#include <eigen3/Eigen/Eigen>

struct DECDirac: public DEC {
  Eigen::SparseMatrix<double> L0, L1, L2;
#ifdef COARSE_SOLVE_PINV
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> direct_solver;
#else
  Eigen::SparseLU<Eigen::SparseMatrix<double>> direct_solver;
#endif

  DECDirac(Mesh2D& mesh): DEC(mesh) {
    buildLaplace();
  }

  void buildCoarseSolver();
  void applyDirac(const Cochain& src, Cochain& dest) const;
  void buildLaplace();
  Cochain residual(const Cochain& u, const Cochain& rhs) const;
  void DGS(Cochain& u, const Cochain& rhs) const;
  Cochain directSolve(const Cochain& rhs) const; 
};

#endif
