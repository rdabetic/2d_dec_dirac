#include "mesh.hpp"
#include "incidence.hpp"
#include "bd.hpp"

void Mesh2D::update() {
  auto triangleArea = [&](int elem_id) -> double {
    // Get vertex ids
    mfem::Array<int> vids;
    this->GetElementVertices(elem_id, vids);

    Eigen::Matrix3d D;
    // Iterate over vertices and add to matrix
    for(int i = 0; i < 3; ++i) {
      const int vid = vids[i];
      D(0, i) = vertices[vid](0);
      D(1, i) = vertices[vid](1);
      D(2, i) = 1;
    }
    return std::abs(D.determinant()) / 2.;
  };

  vertex_edge = VertexEdgeIncidence(*this);
  edge_face = EdgeFaceIncidence(*this);

  triangle_areas = Eigen::VectorXd::NullaryExpr(this->GetNE(), 
      [&](int i) {return triangleArea(i);}
  );
  // Get essential boundary info
  auto [is_essential_pt_, n_essential_pts_] = 
    extractBDPts(*this);

  ess_v = std::move(is_essential_pt_);
  n_ess_v = n_essential_pts_;

  auto [is_essential_edge_, n_essential_edges_] = 
    extractBDEdges(*this);

  ess_e = std::move(is_essential_edge_);
  n_ess_e = n_essential_edges_;
#ifndef BARYCENTRIC_DUAL
  assertCircumcentric();
#endif
}

Mesh2D::Mesh2D(mfem::Mesh&& mesh): mfem::Mesh(mesh) {
  update();
}

void Mesh2D::refine() {
  UniformRefinement();
  update();
}

/// Get barycenter of element `elem_id`
Eigen::Vector2d Mesh2D::barycenter(int elem_id) {
  Eigen::Vector2d c;
  c.setZero();
  // Get vertex ids
  mfem::Array<int> vids;
  this->GetElementVertices(elem_id, vids);

  // Iterate over vertices and add
  for(int i = 0; i < 3; ++i) {
    const int vid = vids[i];
    c(0) += vertices[vid](0) / 3.;
    c(1) += vertices[vid](1) / 3.;
  }
  return c;
}

/// Circumcenter (or barycenter if BARYCENTRIC_DUAL is defined)
/// of element `elem_id`
Eigen::Vector2d Mesh2D::center(int elem_id) {
#ifdef BARYCENTRIC_DUAL
  return barycenter(elem_id);
#else
  mfem::Array<int> vids;
  this->GetElementVertices(elem_id, vids);
  const double ax = vertices[vids[0]](0),
               ay = vertices[vids[0]](1),
               bx = vertices[vids[1]](0),
               by = vertices[vids[1]](1),
               cx = vertices[vids[2]](0),
               cy = vertices[vids[2]](1);
  const double d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
  const double ux = ((ax * ax + ay * ay) * (by - cy) +
                     (bx * bx + by * by) * (cy - ay) +
                     (cx * cx + cy * cy) * (ay - by)) / d;
  const double uy = ((ax * ax + ay * ay) * (cx - bx) +
                     (bx * bx + by * by) * (ax - cx) + 
                     (cx * cx + cy * cy) * (bx - ax)) / d;
  return Eigen::Vector2d({ux, uy});
#endif
}

Eigen::Vector2d Mesh2D::vert2vec(int vid) const {
  return Eigen::Vector2d({vertices[vid](0), 
                          vertices[vid](1)});
}


void Mesh2D::assertCircumcentric() const {
  double s0, s1, s2; // side lengths SQUARED
  mfem::Array<int> vert;
  for(unsigned int elem_id = 0; elem_id < NumOfElements; ++elem_id) {
    GetElementVertices(elem_id, vert);
    const Eigen::Vector2d
      v0 = vert2vec(vert[0]),
      v1 = vert2vec(vert[1]),
      v2 = vert2vec(vert[2]);
    s0 = (v1 - v0).squaredNorm();
    s1 = (v2 - v1).squaredNorm();
    s2 = (v0 - v2).squaredNorm();
    const double max_len = std::max({s0, s1, s2});
    if(2 * max_len >= s0 + s1 + s2) {
      std::cerr << "Triangle with id " << elem_id << " not acute!" << std::endl
                << "Vertices are: (rowwise)" << std::endl
                << v0.transpose() << std::endl
                << v1.transpose() << std::endl
                << v2.transpose() << std::endl;
      std::abort();
    }
  }
}


double Mesh2D::meshWidth() const {
  double h = 0;
  mfem::Array<int> vert(2);
  for(int ei = 0; ei < NumOfEdges; ++ei){
    // Get vertices of edge `ei`
    GetEdgeVertices(ei, vert);
    h = std::max(h, (vert2vec(vert[1]) - 
                     vert2vec(vert[0]))
                    .norm());
  }
  return h;
}
