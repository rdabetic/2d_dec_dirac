#include "mfem.hpp"
#include <mfem/fem/fe_coll.hpp>
#include <mfem/fem/fespace.hpp>
#include <vector>

std::pair<std::vector<bool>,
          unsigned int>
extractBDPts(mfem::Mesh& mesh) {
  // Results
  std::vector<bool> points(mesh.GetNV(), false);
  unsigned int n_bd_pts = 0;

  mfem::FiniteElementCollection* H1_ = new mfem::H1_FECollection(1, 2);
  mfem::FiniteElementSpace* H1 = new mfem::FiniteElementSpace(&mesh, H1_);

  mfem::Array<int> bd_v;
  H1->GetBoundaryTrueDofs(bd_v);
  for(const auto& idx: bd_v) {
    points[idx] = true;
    ++n_bd_pts;
  }
  delete H1; delete H1_;

  return std::make_pair(points, n_bd_pts);
}

std::pair<std::vector<bool>,
          unsigned int>
extractBDEdges(mfem::Mesh& mesh) {
  // Results
  std::vector<bool> edges(mesh.GetNEdges(), false);
  unsigned int n_bd_edges = 0;

  mfem::FiniteElementCollection* HCurl_ = new mfem::ND_FECollection(1, 2);
  mfem::FiniteElementSpace* HCurl = new mfem::FiniteElementSpace(&mesh, HCurl_);

  mfem::Array<int> bd_e;
  HCurl->GetBoundaryTrueDofs(bd_e);
  for(const auto& idx: bd_e) {
    edges[idx] = true;
    ++n_bd_edges;
  }
  delete HCurl; delete HCurl_;

  return std::make_pair(edges, n_bd_edges);
}


