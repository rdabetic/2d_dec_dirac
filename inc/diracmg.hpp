#ifndef DIRACMG_H
#define DIRACMG_H

#include "hierarchy.hpp"
#include "dirac.hpp"

struct DiracMG: DeRhamFESpaceHierarchy {
  unsigned int smoothening_steps = 1;
  std::vector<DECDirac*> diracs;
  /// Construct with a mesh; takes ownership
  DiracMG(Mesh2D& mesh_, unsigned int nsmooth = 1);
  /// Construct with a mesh; takes ownership and calls addMGLevel nref times
  DiracMG(Mesh2D& mesh_, 
          unsigned int nref,
          unsigned int nsmooth);
  ~DiracMG();
  /// Add level by uniform refinement and update Dirac operators
  /// Beware: addLevel() does NOT update the Dirac operators
  void addMGLevel();
  Cochain createCochain() const;
  Cochain getRandomCochain() const;

  Cochain prolongCochain(Cochain& src,
                         unsigned int lvl) const;
  Cochain restrictCochain(Cochain& src,
                          unsigned int lvl) const;

  void MGIt(const Cochain& rhs,
                  Cochain& u,
            double tol,
            unsigned int max_steps = 1,
            unsigned int rec_steps = 1,
            unsigned int lvl = 0) const;

  void MGCycle(const Cochain& rhs,
                     Cochain& u,
                     unsigned int rec_steps = 1) const;

  Cochain FMG(Cochain& rhs,
              double tol,
              unsigned int max_steps) const;
  Cochain FMG_(Cochain& rhs, 
               double tol,
               unsigned int max_steps,
               unsigned int lvl) const;

  template <typename FUNC0,
            typename FUNC1,
            typename FUNC2>
  Cochain getLoadVector(FUNC0 f0, FUNC1 f1, FUNC2 f2) const {
    return diracs.back()->deRham(f0, f1, f2);
  }
};

#endif
