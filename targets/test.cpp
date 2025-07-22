#include "incidence.hpp"
#include "bd.hpp"
#include "mesh.hpp"
#include "dec.hpp"
#include "dirac.hpp"
#include <iostream>

int main (){
  const unsigned int nx = 5,
                     ny = nx;
  mfem::Mesh mesh("/home/me/Downloads/square.msh");

  mesh.Save("fname.msh");
  mesh.UniformRefinement();
  mesh.Save("fname_refined.msh");


  Mesh2D m(std::move(mesh));
  std::cout << "No. elements: "
            << m.GetNE() << std::endl
            << "No. essential vertices: "
            << m.n_ess_v << std::endl
            << "No. essential edges: " 
            << m.n_ess_e << std::endl;


  DECDirac D(m);

  const unsigned int nv = m.GetNV(),
                     ne = m.GetNEdges(),
                     nt = m.GetNE();
  Cochain u(nv, ne, nt);
  Cochain rhs(nv, ne, nt);
  rhs.u0.array() = 1;
  rhs.u1.array() = 1;

  D.zeroCochainBoundary(rhs);
  std::cout << "-- Testing DGSR ---" << std::endl;

  for(unsigned int i = 0; i < 1024; ++i) {
    Cochain res = D.residual(u, rhs);
    //std::cout << res.u0.transpose() << std::endl;
    //D.zeroCochainBoundary(res);
    //std::cout << res.u0.transpose() << std::endl;
    //std::cout << std::endl;
    std::cout << "it " << i << " residual.norm() \t" << res.norm() << std::endl;
    D.DGS(u, rhs);
  }

  return 0;
}
