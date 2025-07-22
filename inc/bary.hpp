#ifndef BARY_H
#define BARY_H

/*
 *  Implements barycentric coodrinate functions for 
 *  an equilateral triangle. These are helper functions for 
 *  manufacturing a solution.
 */

#include <Eigen/Eigen>

// Barucentric coodrinate functions
double l0(const Eigen::Vector2d& v);
double l1(const Eigen::Vector2d& v);
double l2(const Eigen::Vector2d& v);

// Gradients
Eigen::Vector2d gradL0(const Eigen::Vector2d v);
Eigen::Vector2d gradL1(const Eigen::Vector2d v);
Eigen::Vector2d gradL2(const Eigen::Vector2d v);

// Manufactured solution ...
double u0(const Eigen::Vector2d& v);
Eigen::Vector2d u1(const Eigen::Vector2d& v);
// CAREFUL: This does not give something with 
// mean zero, instead, we project to mean zero when 
// assembling the load vector
double u2(const Eigen::Vector2d& v);

Eigen::Vector2d gradU0(const Eigen::Vector2d& v);
Eigen::Vector2d curlU2(const Eigen::Vector2d& v);
double rotU1(const Eigen::Vector2d& v);
double divU1(const Eigen::Vector2d& v);

// ... and RHS
double f0(const Eigen::Vector2d& v);
Eigen::Vector2d f1(const Eigen::Vector2d& v);
double f2(const Eigen::Vector2d& v);

#endif
