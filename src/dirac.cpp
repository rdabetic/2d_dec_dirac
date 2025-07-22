#include "dirac.hpp"
#include <eigen3/unsupported/Eigen/KroneckerProduct>

void DECDirac::applyDirac(const Cochain& src,
                                Cochain& dest) const {
  dest.u0 = delta1 * src.u1;
  dest.u1 = d0 * src.u0 + delta2 * src.u2;
  dest.u2 = d1 * src.u1;
}

void DECDirac::buildLaplace() {
  L0 = delta1 * d0;
  L1 = delta2 * d1 + d0 * delta1;
  L2 = d1 * delta2;
  // Prune entries with boundary
  L0.prune([&](const int r, const int c, const double&) -> bool {
        // Keep?
        return !(mesh.ess_v[r] || mesh.ess_v[c]);
      }
  );
  L1.prune([&](const int r, const int c, const double&) -> bool {
        // Keep?
        return !(mesh.ess_e[r] || mesh.ess_e[c]);
      }
  );
  // Set diagonal entries for boundary DOFS to one so 
  // we can do DGS with Eigen's triangularView
  //
  // Vertices
  Eigen::VectorXd tmp;
  tmp = Eigen::VectorXd::NullaryExpr(mesh.GetNV(),
      [&](int vid) -> double {
        return mesh.ess_v[vid];
      }
  );
  L0 += tmp.asDiagonal();
  // Edges
  tmp = Eigen::VectorXd::NullaryExpr(mesh.GetNEdges(),
      [&](int eid) -> double {
        return mesh.ess_e[eid];
      }
  );
  L1 += tmp.asDiagonal();

  L0.makeCompressed();
  L1.makeCompressed();
  L2.makeCompressed();
  delta1.makeCompressed();
  delta2.makeCompressed();
}


Cochain DECDirac::residual(const Cochain& u, 
                           const Cochain& rhs) const{
  Cochain residual(u.n0, u.n1, u.n2);
  applyDirac(u, residual);
  residual -= rhs;

  return residual;
}

void DECDirac::DGS(      Cochain& u,
                   const Cochain& rhs) const {
  Cochain corr(u.n0, u.n1, u.n2);
  Cochain r = residual(u, rhs);
  // Go to the Laplacian and do GS
  L0.triangularView<Eigen::Lower>().solveInPlace(r.u0);
  L1.triangularView<Eigen::Lower>().solveInPlace(r.u1);
  L2.triangularView<Eigen::Lower>().solveInPlace(r.u2);
  // Distribute
  applyDirac(r, corr);
  u -= corr;
  // Try to combat rounding errors
  u.setZeroMean();
}

void DECDirac::buildCoarseSolver() {
  const unsigned int nv = mesh.GetNV(),
                     ne = mesh.GetNEdges(),
                     nt = mesh.GetNE(),
                     n = nv + ne + nt;
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(4 * (nv + ne) + 2 * nt);

  // d0
  for(unsigned int k = 0; k < d0.outerSize(); ++k)
    for(Eigen::SparseMatrix<double>::InnerIterator it(d0, k); it; ++it)
      if(!mesh.ess_v[it.col()] && !mesh.ess_e[it.row()])
        triplets.emplace_back(it.row() + nv, it.col(), it.value());
  // d1
  for(unsigned int k = 0; k < d1.outerSize(); ++k)
    for(Eigen::SparseMatrix<double>::InnerIterator it(d1, k); it; ++it)
      if(!mesh.ess_e[it.col()])
        triplets.emplace_back(it.row() + nv + ne, it.col() + nv, it.value());
  // delta1
  for(unsigned int k = 0; k < delta1.outerSize(); ++k)
    for(Eigen::SparseMatrix<double>::InnerIterator it(delta1, k); it; ++it)
      if(!mesh.ess_v[it.row()] && !mesh.ess_e[it.col()])
        triplets.emplace_back(it.row(), it.col() + nv, it.value());
  // delta2
  for(unsigned int k = 0; k < delta2.outerSize(); ++k)
    for(Eigen::SparseMatrix<double>::InnerIterator it(delta2, k); it; ++it)
      if(!mesh.ess_e[it.row()])
        triplets.emplace_back(it.row() + nv, it.col() + nv + ne, it.value());

  // Fill the rows on the boundary with 
  // the ones from the identity matrix
  // (the system is singular otherwise)
  for(unsigned int k = 0; k < nv; ++k)
    if(mesh.ess_v[k])
      triplets.emplace_back(k, k, 1);
  for(unsigned int k = 0; k < ne; ++k)
    if(mesh.ess_e[k])
      triplets.emplace_back(k + nv, k + nv, 1);

#ifdef COARSE_SOLVE_PINV
  Eigen::SparseMatrix<double> full_system(n, n);
#else
  // Zero mean condition
  for(unsigned int k = 0; k < nt; ++k) {
    triplets.emplace_back(n, nv + ne + k, 1);
    triplets.emplace_back(nv + ne + k, n, 1);
  }
  Eigen::SparseMatrix<double> full_system(n + 1, n + 1);
#endif
  full_system.setFromTriplets(triplets.begin(), triplets.end());
  direct_solver.compute(full_system);
}

Cochain DECDirac::directSolve(const Cochain& rhs) const {
  const unsigned int nv = mesh.GetNV(),
                     ne = mesh.GetNEdges(),
                     nt = mesh.GetNE(),
                     n = nv + ne + nt;
#ifdef COARSE_SOLVE_PINV
  Eigen::VectorXd u(n);
  u << rhs.u0, rhs.u1, rhs.u2;
#else
  Eigen::VectorXd u(n + 1);
  u << rhs.u0, rhs.u1, rhs.u2, 0;
#endif

  u = direct_solver.solve(u);
  return Cochain(u.segment(0, nv),
                 u.segment(nv, ne),
                 u.segment(nv + ne, nt));
}
