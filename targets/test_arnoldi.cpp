#include "arnoldi.hpp"
#include <eigen3/Eigen/Eigenvalues>

int main() {
  unsigned int s = 32;
  Eigen::MatrixXd R(s, s);
  R = Eigen::VectorXd::LinSpaced(s, 0, s).asDiagonal();

  auto mvec = 
    [&] <typename VEC, typename VEC_>
    (const VEC& src, VEC_& target) {
    target = R * src;
  };

  auto mvec_ = 
    [&] <typename VEC>
    (const VEC& src) {
    return R * src;
  };

  Eigen::VectorXcd ew = deflatedArnoldi(mvec, s, 6, 1e-12);
  std::cout << std::endl
            << ew.cwiseAbs().transpose()
            //<< bruteForceArnoldi(mvec, s, 2)
            << std::endl;
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(R);
  std::cout << es.eigenvalues().cwiseAbs().maxCoeff() << std::endl;
  return 0;
}
