#include "diracmg.hpp"
#include "dirac.hpp"
#include "hierarchy.hpp"

DiracMG::DiracMG(Mesh2D& mesh_,
                 unsigned int nsmooth):
    DeRhamFESpaceHierarchy(mesh_),
    smoothening_steps(nsmooth) {
  // Create new Dirac op
  DECDirac* new_dirac = new DECDirac(*meshes.back());
  diracs.push_back(new_dirac);
  // We are at the coarse level, so we need to
  // assemble a coarse system for multigrid
  new_dirac->buildCoarseSolver();
}

DiracMG::DiracMG(Mesh2D& mesh_, 
                 unsigned int nref,
                 unsigned int nsmooth): 
    DiracMG(mesh_, nsmooth) {
  for(unsigned int k = 0; k < nref; ++k)
    addMGLevel();
}

void DiracMG::addMGLevel(){
  addLevel();
  DECDirac* new_dirac = new DECDirac(*meshes.back());
  diracs.push_back(new_dirac);
}

DiracMG::~DiracMG() {
  // Delete the diracs
  for(auto dptr: diracs)
    delete dptr;
}

// TODO: Optionally leave src modified, the
//       contents of it don't really matter in the
//       multigrid iteration

Cochain DiracMG::prolongCochain(Cochain& src,
                                unsigned int lvl) const {
  assert(lvl < max_lvl);
  Cochain dest = diracs[lvl + 1]->createCochain();
  // Re-scale to get to L2 (not whitney 2-forms)
  src.u2.array() /= meshes[lvl]->triangle_areas.array();
  // Prolong src
  prolong(src.u0_mfem, 
          src.u1_mfem,
          src.u2_mfem,
          dest.u0_mfem,
          dest.u1_mfem,
          dest.u2_mfem,
          lvl);
  // Scale back
  dest.u2.array() *= meshes[lvl + 1]->triangle_areas.array();
  // Scale back src
  src.u2.array() *= meshes[lvl]->triangle_areas.array();

  return dest;
}

Cochain DiracMG::restrictCochain(Cochain& src,
                                 unsigned int lvl) const {
  assert(lvl > 0);
  Cochain dest = diracs[lvl - 1]->createCochain();
  // Hodge
  diracs[lvl]->applyHodgeStar(src);
  // Transpose of the prolongation
  src.u2.array() *= meshes[lvl]->triangle_areas.array();
  restrict(src.u0_mfem, 
           src.u1_mfem,
           src.u2_mfem,
           dest.u0_mfem,
           dest.u1_mfem,
           dest.u2_mfem,
           lvl);
  dest.u2.array() /= meshes[lvl - 1]->triangle_areas.array();
  // Inverse Hodge
  diracs[lvl - 1]->applyInvHodgeStar(dest);
  // Scale back src
  src.u2.array() /= meshes[lvl]->triangle_areas.array();
  diracs[lvl]->applyInvHodgeStar(src);

  return dest;
}

void DiracMG::MGIt(const Cochain& rhs,
                         Cochain& u,
                   double tol,
                   unsigned int max_steps,
                   unsigned int rec_steps,
                   unsigned int lvl) const {
  if(lvl == 0) {
    u = diracs[0]->directSolve(rhs);
  }
  else {
    for(unsigned int n = 0; n < max_steps; ++n) {
      Cochain u_old(u);
      // Pre-smoothening
      for(unsigned int k = 0;
          k < smoothening_steps;
          ++k)
        diracs[lvl]->DGS(u, rhs);
      // Residual
      Cochain residual = diracs[lvl]->residual(u, rhs);
      residual = restrictCochain(residual, lvl);
      // Initial guess, zero
      Cochain corr(residual.n0, residual.n1, residual.n2);
      // Recursion
      MGIt(residual, corr, 0,
           rec_steps, rec_steps, lvl - 1);
      // Prolong and update
      corr = prolongCochain(corr, lvl - 1);
      u -= corr;
      // Post-smoothening
      for(unsigned int k = 0;
          k < smoothening_steps;
          ++k)
        diracs[lvl]->DGS(u, rhs);
      // Check for termination
      if(tol != 0) {
        u_old -= u;
        if(u_old.norm() < tol * u.norm())
          break;
      }
    }
  }
}

void DiracMG::MGCycle(const Cochain &rhs, 
                            Cochain& u, 
                            unsigned int rec_steps) const {
  MGIt(rhs, u, 0, 1, rec_steps, max_lvl);
}

Cochain DiracMG::FMG_(Cochain& rhs, 
                      double tol,
                      unsigned int max_steps,
                      unsigned int lvl) const {
  if(lvl == 0) {
    return diracs[0]->directSolve(rhs);
  } else {
    Cochain u;
    { // rhs_restricted can go out of scope, we don't need it after this
      Cochain rhs_restricted = restrictCochain(rhs, lvl);
      u = FMG_(rhs_restricted,
                       2 * tol,
                       max_steps,
                       lvl - 1);
    }
    u = prolongCochain(u, lvl - 1);
    MGIt(rhs, u, tol, max_steps, 1, lvl);
    return u;
  }
}

Cochain DiracMG::FMG(Cochain& rhs,
                     double tol,
                     unsigned int max_steps) const {
  return FMG_(rhs, tol, max_steps, max_lvl);
}


Cochain DiracMG::createCochain() const {
  return diracs.back()->createCochain();
}
Cochain DiracMG::getRandomCochain() const {
  return diracs.back()->getRandomCochain();
}
