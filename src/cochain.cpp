#include "cochain.hpp"

Cochain::Cochain(unsigned int n0_, unsigned int n1_, unsigned int n2_):
    n0(n0_), n1(n1_), n2(n2_) {
  u0 = Eigen::VectorXd::Zero(n0);
  u1 = Eigen::VectorXd::Zero(n1);
  u2 = Eigen::VectorXd::Zero(n2);
  // Share the data
  u0_mfem = mfem::Vector(&u0[0], n0);
  u1_mfem = mfem::Vector(&u1[0], n1);
  u2_mfem = mfem::Vector(&u2[0], n2);
}

Cochain::Cochain(const Eigen::VectorXd& u0_,
                 const Eigen::VectorXd& u1_,
                 const Eigen::VectorXd& u2_) {
  n0 = u0_.size();
  n1 = u1_.size();
  n2 = u2_.size();

  u0 = u0_;
  u1 = u1_;
  u2 = u2_;
  // Share the data
  u0_mfem = mfem::Vector(&u0[0], n0);
  u1_mfem = mfem::Vector(&u1[0], n1);
  u2_mfem = mfem::Vector(&u2[0], n2);
}

Cochain& Cochain::operator=(const Cochain& other) {
  n0 = other.n0;
  n1 = other.n1;
  n2 = other.n2;

  u0 = other.u0;
  u1 = other.u1;
  u2 = other.u2;

  // We cannoy copy the vector, as
  // we want the mfem::Vector and Eigen's vectors
  // to share data
  u0_mfem = mfem::Vector(&u0[0], n0);
  u1_mfem = mfem::Vector(&u1[0], n1);
  u2_mfem = mfem::Vector(&u2[0], n2);

  return *this;
}

std::ostream& operator<<(std::ostream& out,
                         const Cochain& c) {
  out << c.u0.transpose() << std::endl 
      << c.u1.transpose() << std::endl
      << c.u2.transpose() << std::endl;
  return out;
}
