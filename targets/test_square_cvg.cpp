#include "mesh.hpp"
#include "cochain.hpp"
#include "dirac.hpp"
#include <eigen3/Eigen/Eigen>
#include <ios>
#include "diracmg.hpp"

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
  // Load `functions`
  auto f0 = [](const Eigen::Vector2d& v) -> double {
    return 0;
  };
  auto f1 = [](const Eigen::Vector2d& v) -> Eigen::Vector2d {
    const double& x = v(0),
                  y = v(1);
    return Eigen::Vector2d({
        0,
        4 * M_PI * std::sin(2 * M_PI * x) * std::cos(2 * M_PI * y)
    });
  };
  auto f2 = [](const Eigen::Vector2d v) -> double {
    const double& x = v(0),
                  y = v(1);
    return 2 * M_PI * (std::cos(2 * M_PI * x) - std::cos(2 * M_PI * y));
  };
  // The exact solution
  auto u0 = [](const Eigen::Vector2d& v) -> double {
    const double& x = v(0),
                  y = v(1);
    return std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y);
  };
  auto u1 = [](const Eigen::Vector2d& v) -> Eigen::Vector2d {
    const double& x = v(0),
                  y = v(1);
    return Eigen::Vector2d({
      std::sin(2 * M_PI * y),
      std::sin(2 * M_PI * x)
    });
  };
  auto u2 = [](const Eigen::Vector2d& v) -> double {
    const double& x = v(0),
                  y = v(1);
    return std::cos(2 * M_PI * x) * std::cos(2 * M_PI * y);
  };
  // Mesh
  mfem::Mesh mesh__(fname);
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
    Cochain u = MG.FMG(load, 1e-3 * h, 10'000);
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
