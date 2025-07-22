#include "hierarchy.hpp"
#include <mfem/fem/fe_coll.hpp>
#include <mfem/fem/transfer.hpp>
#include <mfem/linalg/handle.hpp>

DeRhamFESpaceHierarchy::DeRhamFESpaceHierarchy(Mesh2D& mesh_) {
  // The mesh
  meshes.push_back(&mesh_);
  // Local, for convenience
  Mesh2D& mesh = *meshes.back();
  
  mfem::FiniteElementCollection* H1_ = 
    new mfem::H1_FECollection(1, 2);
  h1_col.push_back(H1_);
  mfem::FiniteElementSpace* H1 = 
    new mfem::FiniteElementSpace(&mesh, H1_);
  h1_spaces.push_back(H1);

  mfem::FiniteElementCollection* HCurl_ = 
    new mfem::ND_FECollection(1, 2);
  hcurl_col.push_back(HCurl_);
  mfem::FiniteElementSpace* HCurl = 
    new mfem::FiniteElementSpace(&mesh, HCurl_);
  hcurl_spaces.push_back(HCurl);

  mfem::FiniteElementCollection* L2_ = 
    new mfem::L2_FECollection(0, 2);
  l2_col.push_back(L2_);
  mfem::FiniteElementSpace* L2 = 
    new mfem::FiniteElementSpace(&mesh, L2_);
  l2_spaces.push_back(L2);
}

DeRhamFESpaceHierarchy::DeRhamFESpaceHierarchy(Mesh2D& mesh_, 
                                               unsigned int nref):
      DeRhamFESpaceHierarchy(mesh_) {
  for(unsigned int k = 0; k < nref; ++k)
    addLevel();
}

void DeRhamFESpaceHierarchy::addLevel() {
  // Copy mesh
  Mesh2D* new_mesh = new Mesh2D(*meshes.back());
  // Refine
  new_mesh->refine();
  meshes.push_back(new_mesh);
  // Set up new FE spaces
  Mesh2D& mesh = *meshes.back();

  mfem::FiniteElementCollection* H1_ = 
    new mfem::H1_FECollection(1, 2);
  mfem::FiniteElementSpace* H1 = 
    new mfem::FiniteElementSpace(&mesh, H1_);

  mfem::FiniteElementCollection* HCurl_ = 
    new mfem::ND_FECollection(1, 2);
  mfem::FiniteElementSpace* HCurl = 
    new mfem::FiniteElementSpace(&mesh, HCurl_);

  mfem::FiniteElementCollection* L2_ = 
    new mfem::L2_FECollection(0, 2);
  mfem::FiniteElementSpace* L2 = 
    new mfem::FiniteElementSpace(&mesh, L2_);

  // Add prolongation OP
  mfem::OperatorPtr P_H1(mfem::Operator::ANY_TYPE);
  H1->GetTransferOperator(*h1_spaces.back(), P_H1);
  P_H1.SetOperatorOwner(false);
  mfem::Operator* h1_prol_new = P_H1.Ptr();
  // For some reason gives a mfem::PRefinementTransferOperator???
  // We don't want p-refinement!!!
  //mfem::Operator* h1_prol_new = 
    //new mfem::TransferOperator(*h1_spaces.back(), *H1);
  h1_prol.push_back(h1_prol_new);

  mfem::OperatorPtr P_HCurl(mfem::Operator::ANY_TYPE);
  HCurl->GetTransferOperator(*hcurl_spaces.back(), P_HCurl);
  P_HCurl.SetOperatorOwner(false);
  mfem::Operator* hcurl_prol_new = P_HCurl.Ptr();
  //mfem::Operator* hcurl_prol_new = 
  //  new mfem::TransferOperator(*hcurl_spaces.back(), *HCurl);
  hcurl_prol.push_back(hcurl_prol_new);

  mfem::OperatorPtr P_L2(mfem::Operator::ANY_TYPE);
  L2->GetTransferOperator(*l2_spaces.back(), P_L2);
  P_L2.SetOperatorOwner(false);
  mfem::Operator* l2_prol_new = P_L2.Ptr();
  //mfem::Operator* l2_prol_new = 
  //  new mfem::TransferOperator(*l2_spaces.back(), *L2);
  l2_prol.push_back(l2_prol_new);

  // Add the spaces so they can be deleted later
  h1_col.push_back(H1_);
  h1_spaces.push_back(H1);
  hcurl_col.push_back(HCurl_);
  hcurl_spaces.push_back(HCurl);
  l2_col.push_back(L2_);
  l2_spaces.push_back(L2);

  ++max_lvl;
}

void DeRhamFESpaceHierarchy::prolong(const mfem::Vector& src0,
                                     const mfem::Vector& src1,
                                     const mfem::Vector& src2,
                                           mfem::Vector& dest0,
                                           mfem::Vector& dest1,
                                           mfem::Vector& dest2,
                                     unsigned int lvl) const {
  assert(lvl < max_lvl);
  h1_prol[lvl]->Mult(src0, dest0);
  hcurl_prol[lvl]->Mult(src1, dest1);
  l2_prol[lvl]->Mult(src2, dest2);
}

void DeRhamFESpaceHierarchy::restrict(const mfem::Vector& src0,
                                      const mfem::Vector& src1,
                                      const mfem::Vector& src2,
                                            mfem::Vector& dest0,
                                            mfem::Vector& dest1,
                                            mfem::Vector& dest2,
                                      unsigned int lvl) const {
  assert(lvl > 0);
  h1_prol[lvl - 1]->MultTranspose(src0, dest0);
  hcurl_prol[lvl - 1]->MultTranspose(src1, dest1);
  l2_prol[lvl - 1]->MultTranspose(src2, dest2);
}

DeRhamFESpaceHierarchy::~DeRhamFESpaceHierarchy() {
  for(unsigned int k = 0; k < h1_col.size(); ++k) {
    if(k > 0) {
      delete h1_prol[k - 1];
      delete hcurl_prol[k - 1];
      delete l2_prol[k - 1];
      // meshes[0] is given, don't delete
      delete meshes[k];
    }
    delete h1_col[k];
    delete hcurl_col[k];
    delete l2_col[k];
    delete h1_spaces[k];
    delete hcurl_spaces[k];
    delete l2_spaces[k];
  }
}
