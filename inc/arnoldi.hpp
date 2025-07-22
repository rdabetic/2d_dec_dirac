#include <Spectra/MatOp/DenseGenMatProd.h>
#include <Spectra/Util/CompInfo.h>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include <Spectra/GenEigsSolver.h>

template <typename MVEC>
std::pair<Eigen::MatrixXcd,
          Eigen::MatrixXcd>
arnoldi(MVEC mvec, unsigned int size, 
        Eigen::VectorXcd& v0, unsigned int k = 8) {

  Eigen::MatrixXcd H(k, k), 
                   V(size, k);

  V.setZero();
  H.setZero();

  V.col(0) = v0;

  Eigen::VectorXcd v(size);

  for(unsigned int l = 0; l + 1 < k; ++l) {
    v = mvec(V.col(l));
    for(unsigned int j = 0; j < l + 1; ++j) {
      H(j, l) = V.col(j).dot(v);
      v -= H(j, l) * V.col(j);
    }
    H(l + 1, l) = v.norm();
    if(std::abs(H(l + 1, l)) < 1e-15)
      break;
    V.col(l + 1) = v / H(l + 1, l);
  }
  
  return std::make_pair(std::move(V), std::move(H));
}

template <typename MVEC>
std::pair<Eigen::MatrixXd,
          Eigen::MatrixXd>
realArnoldi(MVEC mvec, unsigned int size, 
            unsigned int k) {

  Eigen::MatrixXd H(k, k), 
                  V(size, k);

  V.col(0).setRandom().normalize();
  H.setZero();

  Eigen::VectorXd v(size);

  for(unsigned int l = 0; l + 1 < k; ++l) {
    mvec(V.col(l), v);
    for(unsigned int j = 0; j < l + 1; ++j) {
      H(j, l) = V.col(j).dot(v);
      v -= H(j, l) * V.col(j);
    }
    H(l + 1, l) = v.norm();
    if(std::abs(H(l + 1, l)) < 1e-15)
      break;
    V.col(l + 1) = v / H(l + 1, l);
  }
  
  return std::make_pair(V, H);
}

template <typename MVEC>
Eigen::VectorXcd
krylovEW(MVEC mvec, unsigned int size,
         unsigned int nev = 1,
         unsigned int k = 12) {
  auto [V, H] = realArnoldi(mvec, size, k);

  Spectra::DenseGenMatProd<double> op(H);
  Spectra::GenEigsSolver es(op, nev, std::min(H.rows(), 4 * nev));

  es.init();
  es.compute();

  return es.eigenvalues();
}

template <typename MATRIX>
std::pair<Eigen::VectorXcd,
          std::complex<double>>
maxEWV(const MATRIX& M) {
  Spectra::DenseGenMatProd<std::complex<double>> op(M);
  Spectra::GenEigsSolver es(op, 1, std::min(8, (int) M.rows()));

  es.init();
  es.compute();
  assert(es.info() == Spectra::CompInfo::Successful);

  return std::make_pair(es.eigenvectors().col(0), es.eigenvalues()(0));
}


template <typename MVEC>
std::pair<Eigen::VectorXcd, std::complex<double>>
arnoldiLargest(MVEC mvec, unsigned int size,
               unsigned int k = 4, double tol = 1e-3) {
  Eigen::VectorXcd v(size);
  std::complex<double> lam = 1. + tol, lam_old = 0;

  v.setRandom().normalize();

  do {
    auto [V, H] = arnoldi(mvec, size, v, k);
    // Eigenvector
    auto [ev, ew] = maxEWV(H);

    lam_old = lam;
    lam = std::move(ew);

    v = V * ev;
    v.normalize();
  } while(std::abs(lam - lam_old) > tol);

  return std::make_pair(v, lam);
}


template <typename MVEC>
Eigen::VectorXcd
deflatedArnoldi(MVEC mvec, 
                unsigned int size,
                unsigned int nev,
                double tol = 1e-3,
                int k = -1) {
  if(k < 0)
    k = 2 * nev + 1;
  assert(k >= 2 * nev + 1);

  // Perform mvec with complex vectors
  auto mvec_ = [&](const Eigen::VectorXcd& v) -> Eigen::VectorXcd {
    Eigen::VectorXd Av_r(size), Av_i(size);
    mvec(v.real(), Av_r);
    mvec(v.imag(), Av_i);
    return std::move(Av_r + std::complex<double>(0, 1) * Av_i);
  };

  Eigen::VectorXcd ews(nev);
  Eigen::MatrixXcd evs(size, nev);

  for(unsigned int n = 0; n < nev; ++n) {
    // Operator, deflated with a projection
    auto mvec_deflated =
      [&] <typename VEC> (const VEC& v) {
        Eigen::VectorXcd w = mvec_(v);
        // Deflate
        for(unsigned int l = 0; l < n; ++l)
          w -= evs.col(l).dot(w) * evs.col(l);
        return w;
      }
    ;
    // Compute the eigenvalues of the deflated operator
    auto [ev, ew] = arnoldiLargest(mvec_deflated, size, k, tol);
    // Update vectors
    evs.col(n) = ev;
    ews(n) = ew;
  }

  return ews;
}
