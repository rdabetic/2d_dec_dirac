#include "bary.hpp"

double l0(const Eigen::Vector2d& v) {
  return 1 - v(0) - v(1) / std::sqrt(3);
}

double l1(const Eigen::Vector2d& v) {
  return v(0) - v(1) / std::sqrt(3);
}

double l2(const Eigen::Vector2d& v) {
  return 2 * v(1) / std::sqrt(3);
}


Eigen::Vector2d gradL0(const Eigen::Vector2d v) {
  return Eigen::Vector2d({
      -1., 
      -1. / std::sqrt(3)
  });
}

Eigen::Vector2d gradL1(const Eigen::Vector2d v) {
  return Eigen::Vector2d({
      1., 
      -1. / std::sqrt(3)
  });
}

Eigen::Vector2d gradL2(const Eigen::Vector2d v) {
  return Eigen::Vector2d({
      0, 
      2. / std::sqrt(3)
  });
}


double u0(const Eigen::Vector2d& v) {
  return l0(v) * l1(v) * l2(v);
}

Eigen::Vector2d u1(const Eigen::Vector2d& v) {
  return Eigen::Vector2d({
      u0(v),
      u0(v),
  });
}

double u2(const Eigen::Vector2d& v) {
  // Subtract the mean, which is 1 / 60
  return u0(v) - 1. / 60.;
}


Eigen::Vector2d gradU0(const Eigen::Vector2d& v) {
  return l0(v) * l1(v) * gradL2(v) + 
         l0(v) * l2(v) * gradL1(v) + 
         l1(v) * l2(v) * gradL0(v);
}

Eigen::Vector2d curlU2(const Eigen::Vector2d& v) {
  Eigen::Vector2d d = gradU0(v);
  // Rotate by 90 degrees
  d(0) = -d(0);
  std::swap(d(0), d(1));

  return d;
}

double rotU1(const Eigen::Vector2d& v) {
  const Eigen::Vector2d d = gradU0(v);
  // \partial_x u_0 - \partial_y u_0
  return d(0) - d(1);
}

double divU1(const Eigen::Vector2d& v) {
  const Eigen::Vector2d d = gradU0(v);
  // \partial_x u_0 + \partial_y u_0
  return d(0) + d(1);
}


double f0(const Eigen::Vector2d& v) {
  return -divU1(v);
}

Eigen::Vector2d f1(const Eigen::Vector2d& v) {
  return gradU0(v) + curlU2(v);
}

double f2(const Eigen::Vector2d& v) {
  return rotU1(v);
}
