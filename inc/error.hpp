#ifndef MGERROR_H
#define MGERROR_H

#include <eigen3/Eigen/Eigen>
#include "diracmg.hpp"

struct MGErrorMatrix {
  using Scalar = double;
  using Matrix = Eigen::MatrixXd;

  DiracMG& MG;

  MGErrorMatrix(DiracMG& mg): MG(mg) {};

  unsigned int rows() const;
  unsigned int cols() const;
  void perform_op(const Scalar* x_in,
                        Scalar* x_out) const;
  Matrix operator*(const Eigen::Ref<const Matrix>& mat_in) const;
  Scalar operator()(unsigned int i, 
                    unsigned int j);
};

Eigen::VectorXcd getEW(DiracMG& mg, unsigned int nev,
                       double tol = 1e-3);

template <typename OP_T>
std::pair<double, Cochain>
powerIt(DiracMG& MG, OP_T inplace_op, double tol = 1e-3) {
  Cochain v   = MG.getRandomCochain();
  double l_old = 0, 
         l_new = 1 + tol;

  do {
    l_old = l_new;
    const double vnorm = v.norm();
    // Normalize
    v.u0.array() /= vnorm;
    v.u1.array() /= vnorm;
    v.u2.array() /= vnorm;
    // Apply operator
    inplace_op(v);

    l_new = v.norm();
    v.setZeroMean();
  } while(std::abs(l_new - l_old) > tol * std::abs(l_new));

  const double vnorm = v.norm();
  v.u0.array() /= vnorm;
  v.u1.array() /= vnorm;
  v.u2.array() /= vnorm;

  return std::make_pair(l_new, v);
}

Eigen::VectorXd deflatedPowerIt(DiracMG& MG, unsigned int nev,
                                double tol = 1e-3);

Eigen::VectorXcd deflatedArnoldiMG(DiracMG& MG, unsigned int nev,
                                   double tol = 1e-3, int k = -1);

#endif
