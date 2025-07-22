#include "mesh.hpp"
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
  Mesh2D mesh_(std::move(mesh__));
  DiracMG MG(mesh_, 4);
  // Find errors
  double errn_old, h_old;
  for(unsigned int lvl = 0; lvl < nref + 1; ++lvl) {
    const Mesh2D& mesh = *MG.meshes.back();
    std::cout << "Refinement level " << lvl << std::endl;
    const double h = mesh.meshWidth();
    std::cout << "\tMesh-width: " 
              << h << std::endl;

    const DECDirac& D = *MG.diracs.back();

    auto zero = [&] (const Eigen::Vector2d& v) {
      return 0.;
    };
    auto zero_vec = [&] (const Eigen::Vector2d& v) {
      Eigen::Vector2d w;
      w.setZero();
      return w;
    };

    Cochain pi_du = D.deRham(zero, gradU0, rotU1);
    Cochain piu = D.deRham(u0, u1, zero);

    //Cochain pi_du = D.deRham(zero, zero_vec, curlU1, zero);
    //Cochain piu = D.deRham(zero, u1, zero_vec, zero);

    //Cochain pi_du = D.deRham(zero, zero_vec, zero_vec, divU2);
    //Cochain piu = D.deRham(zero, zero_vec, u1, zero);

    Cochain err(piu.n0, piu.n1, piu.n2);
    D.diff(piu, err);
    // err = \dd\Pi u - \Pi\dd u 
    err -= pi_du;

    const double errn = D.sobolevNorm(err);
    std::cout << "\tError in infty norm / h: " << 
      err.u1.cwiseAbs().maxCoeff() / h << std::endl;
    std::cout << "\tError in DEC L2 norm: " << D.norm(err) << std::endl;
    std::cout << "\tError in Sobolev norm: " << errn << std::endl;

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
