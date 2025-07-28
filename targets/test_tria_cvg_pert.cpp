#include "mesh.hpp"
#include "bd.hpp"
#include "cochain.hpp"
#include "dirac.hpp"
#include <eigen3/Eigen/Eigen>
#include <ios>
#include "diracmg.hpp"
#include "bary.hpp"

int main(int argc, char* argv[]) {
  // Get mesh file and number of refinements from arguments
  if(argc != 3) {
    std::cerr << "Wrong number of arguments"
              << std::endl
              << "Arguments: [Mesh file as *.msh] [Number of refinements]"
              << std::endl;
    std::abort();
  }
  const std::string fname = argv[1];
  const unsigned int nref = atoi(argv[2]);

  std::cout << std::scientific;

  // Mesh
  mfem::Mesh mesh__(fname);
  /* Perturb mesh */
  const double pert_size = 0.025;
  Eigen::VectorXd perturbation_ = 
    Eigen::VectorXd::Random(2 * mesh__.GetNV()) * pert_size;
  mfem::Vector perturbation(&perturbation_[0], 2 * mesh__.GetNV());
  // Zero out the perturbations on the boundary (we want the same shape)
  auto [bd_v, n_bd_v] = extractBDPts(mesh__);
  for(unsigned int k = 0; k < bd_v.size(); ++k) {
    if(bd_v[k]) {
      perturbation(k) = 0; 
      perturbation(k + mesh__.GetNV()) = 0; 
    }
  }
  // Move the points
  mesh__.MoveVertices(perturbation);
  // Print the perturbed mesh
  std::ofstream file;
  file.open("perturbed_tria.vtk");
  mesh__.PrintVTK(file);
  file.close();
  /* Convergence analysis */
  Mesh2D mesh_(std::move(mesh__));
  DiracMG MG(mesh_, 3);
  // Find errors
  double errn_old, h_old;
  for(unsigned int lvl = 0; lvl < nref + 1; ++lvl) {
    const Mesh2D& mesh = *MG.meshes.back();
    std::cout << "Refinement level " << lvl << std::endl;
    const double h = mesh.meshWidth();
    std::cout << "\tMesh-width: " 
              << h << std::endl
              << "\tNo. triangles: "
              << mesh.GetNE() << std::endl;

    const DECDirac& D = *MG.diracs.back();
    // Assemble load
    Cochain load = D.deRham(f0, f1, f2),
            err  = D.deRham(u0, u1, u2);
    // Solve
    //D.buildCoarseSolver();
    //u = D.directSolve(load);
    Cochain u = MG.FMG(load, 1e-2 * std::pow(h, 2), 100'000);
    std::cout << "\tSolver residual norm: " 
              << D.residual(u, load).norm() << std::endl;
    err -= u;
    const double errn = D.sobolevNorm(err);
    std::cout << "\tError in DEC L2 norm: " << D.norm(err) << std::endl;
    std::cout << "\tError in DEC Sobolev norm: " << errn << std::endl;

    if(lvl > 0) {
      const double eoc = std::log(errn_old / errn) / std::log(h_old / h);
      std::cout << "\tEOC: " << eoc << std::endl;
    }

    errn_old = errn;
    h_old = h;
    // To avoid out-of-memory errors
    if(lvl < nref)
      MG.addMGLevel();
  }

  return 0;
}
