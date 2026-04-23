//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#pragma once

#include "MFEMMultiAppTransfer.h"

enum TransferMergePolicy
{
  FIRST_VALID,
  LAST_VALID,
  LOWEST_APP_INDEX,
  SUM,
  AVERAGE,
  MIN,
  MAX
};

/**
 * Shared scaffold for MFEM/libMesh shape-evaluation transfers that will eventually support
 * collecting and merging values from multiple source apps at the same target points.
 *
 * The current implementation intentionally preserves the existing per-source transfer behavior
 * while exposing target-point extraction, source evaluation, and merge hooks in one place.
 */
class MFEMGeneralFieldTransfer : public MFEMMultiAppTransfer
{
public:
  static InputParameters validParams();
  MFEMGeneralFieldTransfer(const InputParameters & params);

protected:
  struct CandidateValue
  {
    mfem::real_t value = std::numeric_limits<mfem::real_t>::infinity();
    unsigned int source_app_id = std::numeric_limits<unsigned int>::max();
    bool valid = false;
    mfem::real_t score = std::numeric_limits<mfem::real_t>::infinity();
  };

  struct MergeState
  {
    mfem::real_t value = std::numeric_limits<mfem::real_t>::infinity();
    mfem::real_t score = std::numeric_limits<mfem::real_t>::infinity();
    unsigned int source_app_id = std::numeric_limits<unsigned int>::max();
    unsigned int count = 0;
    bool initialized = false;
  };

  static MooseEnum getTransferMergePolicyOptions();
  void transferVariable(const unsigned int var_index, bool is_target_local) override;
  void execute() override;
  virtual void extractTargetPoints(const unsigned int var_index,
                                   mfem::Vector & target_pt_coords,
                                  int & n_points) = 0;

  virtual void evaluateActiveSourceAtTargetPoints(const unsigned int var_index,
                                                  const mfem::Vector & target_pt_coords,
                                                  const int & n_points,
                                                  std::vector<CandidateValue> & candidates,
                                                  bool is_target_local) = 0;

  virtual void writeResolvedValues(const unsigned int var_index,
                                   const mfem::Vector & target_pt_coords,
                                   int & n_points,
                                   mfem::Vector & resolved_values) = 0;

  virtual mfem::real_t getInvalidValue() const { return getMFEMOutOfMeshValue(); }

  void initializeMergeStates(const std::size_t n_points, std::vector<MergeState> & states) const;
  void mergeCandidates(const std::vector<CandidateValue> & candidates,
                       std::vector<MergeState> & states) const;
  void finalizeResolvedValues(const std::vector<MergeState> & states,
                              mfem::Vector & resolved_values) const;

  const TransferMergePolicy _merge_policy;
  const bool _error_on_conflicting_overlap;
};

#endif
