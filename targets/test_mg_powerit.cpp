#include "diracmg.hpp"
#include <ios>


int main(int argc, char* argv[]) {
  std::cout << std::scientific;
  if(argc < 3) {
    std::cerr << "Wrong number of arguments"
              << std::endl
              << "Arguments: [Mesh file as *.msh] [Number of refinements] {Optional: 'V'- or 'W'-cycle, defaults to a V-cycle} {Optional: number of smoothening steps, defaults to 1}"
              << std::endl;
    std::abort();
  }
  const std::string fname = argv[1];
  const unsigned int nref = atoi(argv[2]);
  unsigned int rec_steps = 1,
               nsmooth = 1;

  if(argc >= 4)
    rec_steps = (*argv[3] == 'V' ? 1 : 2);
  if(argc >= 5)
    nsmooth = atoi(argv[4]);

  mfem::Mesh mesh_(fname);
  //mfem::Mesh mesh_ = mfem::Mesh::MakeCartesian2D(2, 2,
  //    mfem::Element::TRIANGLE);
  Mesh2D mesh(std::move(mesh_));

  const double tol = 1e-3;

  std::cout << "--- Testing "
            << (rec_steps == 1 ? "V-Cycle" : "W-Cycle")
            << " with "
            << nsmooth
            << " smoothening step(s) ---"
            << std::endl;

  DiracMG MG(mesh, 0, nsmooth);
  for(unsigned int lvl = 1; lvl < nref + 1; ++lvl) {
    MG.addMGLevel();
    // Random initial guess for the eigenvector
    Cochain v   = MG.getRandomCochain(),
            Av  = MG.createCochain();
    double l_old = 0, 
           l_new = 1;
    // Power iteration
    do {
      l_old = l_new;
      const double vnorm = v.norm();
      // Normalize
      v.u0.array() /= vnorm;
      v.u1.array() /= vnorm;
      v.u2.array() /= vnorm;
      // Apply operator
      MG.diracs.back()->applyDirac(v, Av);
      // Zero by default
      Cochain tmp = MG.createCochain();
      // W-cycle
      MG.MGCycle(Av, tmp, rec_steps);

      v -= tmp;
      l_new = v.norm();

      v.setZeroMean();
    } while(std::abs(l_new - l_old) > tol * std::abs(l_new));
    std::cout << "Refinement level " << lvl << std::endl
              << "\tMesh-width: " 
              << MG.meshes.back()->meshWidth() << std::endl
              << "\tNo. triangles: " 
              << MG.meshes.back()->GetNE() << std::endl
              << "\tlambda: " << l_new << std::endl;
  }
  return 0;
}
