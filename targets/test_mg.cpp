#include "cochain.hpp"
#include "diracmg.hpp"
#include <ios>
#include "error.hpp"


int main(int argc, char* argv[]) {
  std::cout << std::scientific;
  if(argc < 3) {
    std::cerr << "Wrong number of arguments"
              << std::endl
              << "Arguments: [Mesh file as *.msh] [Number of refinements] {Optional: Number of eigenvalues to esimate, defaults to 1} {Optional: 'V'- or 'W'-cycle, defaults to a V-cycle} {Optional: number of smoothening steps, defaults to 1}"
              << std::endl;
    std::abort();
  }
  const std::string fname = argv[1];
  const unsigned int nref = atoi(argv[2]);
  unsigned int rec_steps = 1,
               nsmooth = 1,
               nev = 1;

  if(argc >= 4)
    nev = atoi(argv[3]);
  if(argc >= 5)
    rec_steps = (*argv[4] == 'V' ? 1 : 2);
  if(argc >= 6)
    nsmooth = atoi(argv[5]);

  mfem::Mesh mesh_(fname);
  //mfem::Mesh mesh_ = mfem::Mesh::MakeCartesian2D(2, 2,
  //    mfem::Element::TRIANGLE);
  Mesh2D mesh(std::move(mesh_));

  std::cout << "--- Testing "
            << (rec_steps == 1 ? "V-Cycle" : "W-Cycle")
            << " with "
            << nsmooth
            << " smoothening step(s) ---"
            << std::endl;

  DiracMG MG(mesh, 0, nsmooth);
  for(unsigned int lvl = 1; lvl < nref + 1; ++lvl) {
    MG.addMGLevel();

    //Eigen::VectorXcd ew = getEW(MG, nev, 1e-2);
    //Eigen::VectorXcd ew = deflatedArnoldiMG(MG, nev);
    Eigen::VectorXd ew = deflatedPowerIt(MG, nev);

    std::cout << "Refinement level " << lvl << std::endl
              << "\tMesh-width: " 
              << MG.meshes.back()->meshWidth() << std::endl
              << "\tNo. triangles: " 
              << MG.meshes.back()->GetNE() << std::endl
              << "\tlambda(s) in absolute value: " << ew.transpose().cwiseAbs()
              << std::endl;
  }
  return 0;
}
