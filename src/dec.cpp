#include "dec.hpp"

DEC::DEC(Mesh2D& mesh_): 
    mesh(mesh_), d0(mesh_.vertex_edge), d1(mesh_.edge_face) {
  this->assembleHodgeStars();
  // Codifferentials
  delta1 = star0_inv.asDiagonal() * d0.transpose() * star1.asDiagonal();
  delta2 = star1_inv.asDiagonal() * d1.transpose() * star2.asDiagonal();
}

Cochain DEC::createCochain() const {
  return Cochain(mesh.GetNV(), 
                 mesh.GetNEdges(),
                 mesh.GetNE());
}

Cochain DEC::getRandomCochain() const {
  Cochain c(Eigen::VectorXd::Random(mesh.GetNV()),
            Eigen::VectorXd::Random(mesh.GetNEdges()),
            Eigen::VectorXd::Random(mesh.GetNE())
           );
  zeroCochainBoundary(c);
  return c;
}

//void DEC::assembleHodgeStars() {
//  // The 0-star
//  mfem::Table* v2el = mesh.GetVertexToElementTable();
//  
//  delete v2el;
//  // The 1-star
//  star1.setZero();
//  mfem::Table* e2el = mesh.GetFaceToElementTable();
//  mfem::Array<int> vert(2);
//  for(int eid = 0; eid < mesh.GetNEdges(); ++eid) {
//    // SKIP the ones on the boundary
//    if(mesh.ess_e[eid])
//      continue;
//    // else
//    mesh.GetEdgeVertices(eid, vert);
//    const double primal_len = std::sqrt( 
//        std::pow(mesh.vertices[vert[1]](0) - mesh.vertices[vert[0]](0), 2) + 
//        std::pow(mesh.vertices[vert[1]](1) - mesh.vertices[vert[0]](1), 2) 
//        );
//    // Get adjacent triangles
//    const int* adj_tria = e2el->GetRow(eid);
//    const double dual_len = (mesh.center(adj_tria[1]) - 
//                             mesh.center(adj_tria[0]))
//                            .norm();
//
//    star1[eid] = primal_len / dual_len;
//  }
//  delete e2el;
//  
//  // The 2-star is just the inverse of the triangle areas
//  star2 = mesh.triangle_areas.cwiseInverse();
//}

void DEC::assembleHodgeStars() {
  // Area of triangle
  auto area = [&](Eigen::Vector2d a,
                  Eigen::Vector2d b,
                  Eigen::Vector2d c) {
      Eigen::Matrix3d tmp;
      tmp << a(0), b(0), c(0),
             a(1), b(1), c(1),
             1, 1, 1;
      return std::abs(tmp.determinant()) / 2.;
    }
  ;
  // The 0- and 1-star
  Eigen::VectorXd primal_len(mesh.GetNEdges()),
                  dual_len(mesh.GetNEdges()),
                  dual_vol(mesh.GetNV());
  dual_vol.setZero();
  primal_len.setZero();
  dual_len.setZero();

  mfem::Array<int> edges(3), cor(3),
                   vert(2);
  for(int elem_id = 0; elem_id < mesh.GetNE(); ++elem_id){
    // Center of the element
    const Eigen::Vector2d center = mesh.center(elem_id);
    mesh.GetElementEdges(elem_id, edges, cor);
    // Iterate over edges
    for(int k = 0; k < edges.Size(); ++k) {
      const int e_id = edges[k];
      // Get edge vertex ids
      mesh.GetEdgeVertices(e_id, vert);
      const Eigen::Vector2d v0 = mesh.vert2vec(vert[0]),
                            v1 = mesh.vert2vec(vert[1]);
      // Midpoint of the edge
      const Eigen::Vector2d midpoint = (v0 + v1) / 2.;
      // Length of the edge
      primal_len(e_id) = (v1 - v0).norm();
      // Add part of the dual edge in this element
      dual_len(e_id) += (center - midpoint).norm();
      // Add parts of dual vol of this element to vertices
      dual_vol(vert[0]) += area(v0, center, midpoint);
      dual_vol(vert[1]) += area(v1, center, midpoint);
    }
  }
  
  star0 = Eigen::VectorXd::NullaryExpr(mesh.GetNV(),
      [&](int vid) -> double {
        return mesh.ess_v[vid] ? 0 : dual_vol(vid);
      }
  );
  star1 = Eigen::VectorXd::NullaryExpr(mesh.GetNEdges(),
      [&](int eid) -> double {
        return mesh.ess_e[eid] ? 0 : dual_len(eid) / primal_len(eid);
      }
  );
  // The inverses
  star0_inv = Eigen::VectorXd::NullaryExpr(mesh.GetNV(),
      [&](int vid) -> double {
        return mesh.ess_v[vid] ? 0 : 1. / star0(vid);
      }
  );
  star1_inv = Eigen::VectorXd::NullaryExpr(mesh.GetNEdges(),
      [&](int eid) -> double {
        return mesh.ess_e[eid] ? 0 : 1. / star1(eid);
      }
  );
  // The 2-star is just the inverse of the triangle areas
  star2 = mesh.triangle_areas.cwiseInverse();
  star2_inv = mesh.triangle_areas;
}

void DEC::applyHodgeStar(Cochain& u) const {
  u.u0.array() *= star0.array();
  u.u1.array() *= star1.array();
  u.u2.array() *= star2.array();
}

void DEC::applyInvHodgeStar(Cochain& u) const {
  u.u0.array() *= star0_inv.array();
  u.u1.array() *= star1_inv.array();
  u.u2.array() *= star2_inv.array();
}

void DEC::zeroCochainBoundary(Cochain& c) const {
  // Zero vertexes
  for(unsigned int i = 0; i < c.u0.size(); ++i)
    if(mesh.ess_v[i])
      c.u0[i] = 0;
  // Zero edges
  for(unsigned int i = 0; i < c.u1.size(); ++i)
    if(mesh.ess_e[i])
      c.u1[i] = 0;
  // For good measure, mean zero
  c.setZeroMean();
}




void DEC::diff(const Cochain& src,
                     Cochain& dest) const {
  // 0-forms do not get mapped to
  dest.u0.setZero();
  dest.u1 = d0 * src.u0;
  dest.u2 = d1 * src.u1;
}

void DEC::codiff(const Cochain& src,
                       Cochain& dest) const {
  dest.u0 = delta1 * src.u1;
  dest.u1 = delta2 * src.u2;
  dest.u2.setZero();
}


double DEC::sobolevInner(const Cochain& u,
                         const Cochain& v) const {
  Cochain du(u.n0, u.n1, u.n2),
          dv(v.n0, v.n1, v.n2); 
  diff(u, du);
  diff(v, dv);
  return inner(u, v) + inner(du, dv);
}

double DEC::sobolevNorm(const Cochain& u) const {
  Cochain du(u.n0, u.n1, u.n2);
  diff(u, du);
  return std::sqrt(inner(u, u) + inner(du, du));
} 
