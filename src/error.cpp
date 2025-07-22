#include "error.hpp"
#include "arnoldi.hpp"
#include <Spectra/GenEigsSolver.h>
#include <Spectra/Util/CompInfo.h>

unsigned int MGErrorMatrix::rows() const {
  return MG.meshes.back()->GetNV() + 
         MG.meshes.back()->GetNEdges() + 
         MG.meshes.back()->GetNE();
}

unsigned int MGErrorMatrix::cols() const {
  return MG.meshes.back()->GetNV() + 
         MG.meshes.back()->GetNEdges() + 
         MG.meshes.back()->GetNE();
}

void MGErrorMatrix::perform_op(const MGErrorMatrix::Scalar* x_in,
                                     MGErrorMatrix::Scalar* x_out) const {
  const unsigned int nv = MG.meshes.back()->GetNV(),
                     ne = MG.meshes.back()->GetNEdges(),
                     nt = MG.meshes.back()->GetNE();
                     
  Eigen::Map<const Eigen::VectorXd> v0(x_in, nv),
                                    v1(x_in + nv, ne), 
                                    v2(x_in + nv + ne, nt);

  Eigen::Map<Eigen::VectorXd> w0(x_out, nv),
                              w1(x_out + nv, ne), 
                              w2(x_out + nv + ne, nt);
  Cochain v(nv, ne, nt);
  Cochain Av(nv, ne, nt);

  v.u0 = v0;
  v.u1 = v1;
  v.u2 = v2;

  MG.diracs.back()->zeroCochainBoundary(v);

  MG.diracs.back()->applyDirac(v, Av);
  Cochain tmp = MG.createCochain();
  MG.MGCycle(Av, tmp, 1);

  v -= tmp;
  MG.diracs.back()->zeroCochainBoundary(v);

  w0 = v.u0;
  w1 = v.u1;
  w2 = v.u2;
}

MGErrorMatrix::Matrix
MGErrorMatrix::operator*(const Eigen::Ref<const MGErrorMatrix::Matrix>& mat_in) const {
  Eigen::MatrixXd mat_out(mat_in.rows(), mat_in.cols());
  // Multiply each column
  for(unsigned int c = 0; c < mat_in.cols(); ++c) {
    const Eigen::VectorXd col = mat_in.col(c);
    perform_op(&col(0), &mat_out(0, c));
  }

   return mat_out;
}

double MGErrorMatrix::operator()(unsigned int i,
                                 unsigned int j) {
  // Panic, undefined
  std::cerr << std::endl
            << "Indexing this operator is not allowed!" 
            << std::endl;
  std::abort();
  return 0;
}

Eigen::VectorXcd getEW(DiracMG& mg, unsigned int nev, double tol) {
  MGErrorMatrix E(mg);

  Spectra::GenEigsSolver<MGErrorMatrix> es(E, nev, 2 * nev + 1);
  es.init();
  es.compute(Spectra::SortRule::LargestMagn, 
      1000, tol, Spectra::SortRule::LargestMagn);

  assert(es.info() == Spectra::CompInfo::Successful);

  return es.eigenvalues();
}

Eigen::VectorXd deflatedPowerIt(DiracMG& MG, unsigned int nev,
                                double tol) {
  const DECDirac& D = *MG.diracs.back();

  Eigen::VectorXd ews(nev);
  std::vector<Cochain> evs;
  Cochain Av = MG.createCochain();

  for(unsigned int n = 0; n < nev; ++n) {
    auto deflated_op = 
      [&](Cochain& v) {
        // Apply the multigrid operator
        MG.diracs.back()->zeroCochainBoundary(v);

        MG.diracs.back()->applyDirac(v, Av);
        Cochain tmp = MG.createCochain();
        MG.MGCycle(Av, tmp, 1);

        v -= tmp;
        MG.diracs.back()->zeroCochainBoundary(v);
        // Deflate
        for(unsigned int k = 0; k < n; ++k)
          v -= evs[k] * D.inner(v, evs[k]);
      }
    ;
    auto [ew, ev] = powerIt(MG, deflated_op, tol);

    ews(n) = ew;
    const double ev_n = D.norm(ev);
    ev *= 1. / ev_n;

    evs.push_back(std::move(ev));
  }

  return ews;
}


Eigen::VectorXcd deflatedArnoldiMG(DiracMG& MG, unsigned int nev,
                                   double tol, int k) {

  const unsigned int nv = MG.meshes.back()->GetNV(),
                     ne = MG.meshes.back()->GetNEdges(),
                     nt = MG.meshes.back()->GetNE();

  Cochain Av = MG.createCochain();
  auto mvec_inplace = 
    [&] <typename VEC, typename VEC_>
    (const VEC& src, VEC_& target) {
      Cochain v(nv, ne, nt);
      v.u0 = src.head(nv);
      v.u1 = src.segment(nv, ne);
      v.u2 = src.segment(nv + ne, nt);

      MG.diracs.back()->zeroCochainBoundary(v);

      MG.diracs.back()->applyDirac(v, Av);
      Cochain tmp = MG.createCochain();
      MG.MGCycle(Av, tmp, 1);

      v -= tmp;
      MG.diracs.back()->zeroCochainBoundary(v);

      target.head(nv) = v.u0;
      target.segment(nv, ne) = v.u1;
      target.segment(nv + ne, nt) = v.u2;
  };

  return deflatedArnoldi(mvec_inplace, nv + ne + nt, nev, tol, k);
}
