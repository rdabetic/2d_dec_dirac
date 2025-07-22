#ifndef HIERARCHY_H
#define HIERARCHY_H

#include "mfem.hpp"
#include "mesh.hpp"

struct DeRhamFESpaceHierarchy {
  unsigned int max_lvl = 0;
  std::vector<Mesh2D*> meshes;
  std::vector<mfem::FiniteElementCollection*> h1_col, hcurl_col, l2_col;
  std::vector<mfem::FiniteElementSpace*> h1_spaces, hcurl_spaces, l2_spaces;
  std::vector<mfem::Operator*> h1_prol, hcurl_prol, l2_prol;
 
  /// Initialize wih mesh, does not take ownership
  DeRhamFESpaceHierarchy(Mesh2D& mesh_);
  /// Initialize wih mesh, does not take ownership, and refines nref times
  DeRhamFESpaceHierarchy(Mesh2D& mesh_, 
                         unsigned int nref);
  ~DeRhamFESpaceHierarchy();

  /// Add level by \emph{uniform} refinement
  void addLevel();

  /// Prolong from level `lvl` to `lvl + 1` (if lvl < max_lvl)
  void prolong(const mfem::Vector& src0,
               const mfem::Vector& src1,
               const mfem::Vector& src2,
                     mfem::Vector& dest0,
                     mfem::Vector& dest1,
                     mfem::Vector& dest2,
               unsigned int lvl) const;

  /// Restrict from level `lvl` to `lvl - 1` (if lvl > 0)
  void restrict(const mfem::Vector& src0,
                const mfem::Vector& src1,
                const mfem::Vector& src2,
                      mfem::Vector& dest0,
                      mfem::Vector& dest1,
                      mfem::Vector& dest2,
                unsigned int lvl) const;
};

#endif
