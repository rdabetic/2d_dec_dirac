#ifndef COCHAIN_H
#define COCHAIN_H

#include <eigen3/Eigen/Eigen>
#include "mfem.hpp"
#include <iostream>

/// Helper struct to store the load and solution vector etc
struct Cochain {
  unsigned int n0, n1, n2;
  // Solution components
  // They OWN the data
  Eigen::VectorXd u0, u1, u2;
  // Easy interoparability between 
  // mfem and Eigen, mfem DOES NOT own the data
  // Reason: We want to be able to prolong and restrict 
  //         using mfem, but mfem wants mfem::Vector
  mfem::Vector u0_mfem, u1_mfem, u2_mfem;
  Cochain() = default;
  /// Constructor, init with zeros given shapes
  Cochain(unsigned int n0_, unsigned int n1_, unsigned int n2_);
  /// Construct from vectors
  Cochain(const Eigen::VectorXd& u0_,
          const Eigen::VectorXd& u1_,
          const Eigen::VectorXd& u2_);

  /// Deepcopy `other` (except for the mfem vectors)
  Cochain& operator=(const Cochain& other);

  /* Convenience */

  inline Cochain& operator+=(const Cochain& other) {
    u0 += other.u0;
    u1 += other.u1;
    u2 += other.u2;
    return *this;
  }

  inline Cochain& operator-=(const Cochain& other) {
    u0 -= other.u0;
    u1 -= other.u1;
    u2 -= other.u2;
    return *this;
  }

  inline Cochain& operator*=(const double& d) {
    u0 *= d;
    u1 *= d;
    u2 *= d;
    return *this;
  }

  inline Cochain operator*(const double& d) {
    Cochain res(*this);
    res *= d;
    return res;
  }

  inline void setZeroMean() {
    u2.array() -= u2.mean();
  }

  inline double norm() const {
    return std::sqrt(u0.squaredNorm() + 
                     u1.squaredNorm() +
                     u2.squaredNorm());
  }

  friend std::ostream& operator<<(std::ostream& out,
                                  const Cochain& c);
};

#endif
