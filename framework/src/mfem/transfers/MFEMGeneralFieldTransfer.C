//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "MFEMGeneralFieldTransfer.h"

InputParameters
MFEMGeneralFieldTransfer::validParams()
{
  InputParameters params = MFEMMultiAppTransfer::validParams();
  params.addParam<MooseEnum>(
      "merge_policy",
      getTransferMergePolicyOptions(),
      "How to merge values when multiple source apps contribute to the same target point.");
  params.addParam<bool>("error_on_conflicting_overlap",
                        true,
                        "Whether to error when more than one valid source contributes to a target "
                        "point while using a non-additive merge policy.");
  return params;
}

MFEMGeneralFieldTransfer::MFEMGeneralFieldTransfer(const InputParameters & params)
  : MFEMMultiAppTransfer(params),
    _merge_policy(getParam<MooseEnum>("merge_policy").getEnum<TransferMergePolicy>()),
    _error_on_conflicting_overlap(getParam<bool>("error_on_conflicting_overlap"))
{
}

void
MFEMGeneralFieldTransfer::execute()
{
TIME_SECTION(
    "MFEMMultiAppTransfer::execute", 5, "Perform transfer to and/or from an MFEMProblem.");
// Send from MFEM problem to libMesh problem
for (const auto var_index : make_range(numToVar()))      
  switch (_current_direction)
  {
    case TO_MULTIAPP:
      setActiveFromProblem(getToMultiApp()->problemBase(), 0);
      for (unsigned int i = 0; i < getToMultiApp()->numGlobalApps(); i++)
        if (getToMultiApp()->hasLocalApp(i))
        {
          setActiveToProblem(getToMultiApp()->appProblemBase(i), i);
          transferVariable(var_index, true);
        }
      break;
    case FROM_MULTIAPP:
      setActiveToProblem(getFromMultiApp()->problemBase(), 0);
      for (unsigned int i = 0; i < getFromMultiApp()->numGlobalApps(); i++)
        if (getFromMultiApp()->hasLocalApp(i))
        {
          setActiveFromProblem(getFromMultiApp()->appProblemBase(i), i);
          transferVariable(var_index, true);
        }
      break;
    case BETWEEN_MULTIAPP:
      for (unsigned int to_app_id = 0; to_app_id < getToMultiApp()->numGlobalApps(); to_app_id++)
      {
        const bool is_target_local = getToMultiApp()->hasLocalApp(to_app_id);
        mfem::Vector target_pt_coords;
        int n_points;        
        if (is_target_local)
        {
          setActiveToProblem(getToMultiApp()->appProblemBase(to_app_id), to_app_id);
          extractTargetPoints(var_index, target_pt_coords, n_points);
        }
        std::vector<MergeState> states;
        initializeMergeStates(n_points, states);
        // Evaluate source variable at target points, and update state of best found values across apps
        for (unsigned int from_app_id = 0; from_app_id < getFromMultiApp()->numGlobalApps();
          from_app_id++)
        {
          if (getFromMultiApp()->hasLocalApp(from_app_id))
          {
            setActiveFromProblem(getFromMultiApp()->appProblemBase(from_app_id), from_app_id);          
            std::vector<CandidateValue> candidates;
            evaluateActiveSourceAtTargetPoints(var_index, target_pt_coords, n_points, candidates, is_target_local);
            mergeCandidates(candidates, states);
          }
        }
        if (is_target_local)
        {
          mfem::Vector resolved_values;
          finalizeResolvedValues(states, resolved_values);
          writeResolvedValues(var_index, target_pt_coords, n_points, resolved_values);          
        }
      }
      break;
  }
}
void
MFEMGeneralFieldTransfer::transferVariable(const unsigned int var_index, bool is_target_local)
{
  mfem::Vector target_pt_coords;
  int n_points;
  if (is_target_local)
    extractTargetPoints(var_index, target_pt_coords, n_points);

  std::vector<CandidateValue> candidates;
  evaluateActiveSourceAtTargetPoints(var_index, target_pt_coords, n_points, candidates, is_target_local);

  if (is_target_local)
  {
    std::vector<MergeState> states;
    initializeMergeStates(n_points, states);
    mergeCandidates(candidates, states);

    mfem::Vector resolved_values;
    finalizeResolvedValues(states, resolved_values);
    writeResolvedValues(var_index, target_pt_coords, n_points, resolved_values);
  }
}

MooseEnum
MFEMGeneralFieldTransfer::getTransferMergePolicyOptions()
{
  return MooseEnum("FIRST_VALID LAST_VALID LOWEST_APP_INDEX SUM AVERAGE MIN MAX", "FIRST_VALID");
}

void
MFEMGeneralFieldTransfer::initializeMergeStates(const std::size_t n_points,
                                                std::vector<MergeState> & states) const
{
  states.assign(n_points, MergeState{});
}

void
MFEMGeneralFieldTransfer::mergeCandidates(const std::vector<CandidateValue> & candidates,
                                          std::vector<MergeState> & states) const
{
  mooseAssert(candidates.size() == states.size(), "Mismatched candidate and state counts");

  for (const auto i : index_range(candidates))
  {
    const auto & candidate = candidates[i];
    auto & state = states[i];

    if (!candidate.valid)
      continue;

    if (_error_on_conflicting_overlap &&
        (state.initialized && _merge_policy != TransferMergePolicy::SUM &&
         _merge_policy != TransferMergePolicy::AVERAGE))
      mooseError("Multiple valid sources contributed to the same target point while using merge "
                 "policy '",
                 static_cast<int>(_merge_policy),
                 "'.");

    switch (_merge_policy)
    {
      case TransferMergePolicy::FIRST_VALID:
        if (!state.initialized || state.value == getInvalidValue())
        {
          state.value = candidate.value;
          state.score = candidate.score;
          state.source_app_id = candidate.source_app_id;
          state.count = 1;
          state.initialized = true;
        }
        break;

      case TransferMergePolicy::LAST_VALID:
        state.value = candidate.value;
        state.score = candidate.score;
        state.source_app_id = candidate.source_app_id;
        state.count = 1;
        state.initialized = true;
        break;

      case TransferMergePolicy::LOWEST_APP_INDEX:
        if (!state.initialized || candidate.source_app_id < state.source_app_id)
        {
          state.value = candidate.value;
          state.score = candidate.score;
          state.source_app_id = candidate.source_app_id;
          state.count = 1;
          state.initialized = true;
        }
        break;

      case TransferMergePolicy::SUM:
      case TransferMergePolicy::AVERAGE:
        state.value += candidate.value;
        state.count += 1;
        state.initialized = true;
        break;

      case TransferMergePolicy::MIN:
        if (!state.initialized || candidate.value < state.value)
        {
          state.value = candidate.value;
          state.score = candidate.score;
          state.source_app_id = candidate.source_app_id;
          state.count = 1;
          state.initialized = true;
        }
        break;

      case TransferMergePolicy::MAX:
        if (!state.initialized || candidate.value > state.value)
        {
          state.value = candidate.value;
          state.score = candidate.score;
          state.source_app_id = candidate.source_app_id;
          state.count = 1;
          state.initialized = true;
        }
        break;
    }
  }
}

void
MFEMGeneralFieldTransfer::finalizeResolvedValues(const std::vector<MergeState> & states,
                                                 mfem::Vector & resolved_values) const
{
  resolved_values.SetSize(states.size());
  resolved_values = getInvalidValue();

  for (const auto i : index_range(states))
  {
    const auto & state = states[i];
    if (!state.initialized)
      continue;

    resolved_values(i) = state.value;
    if (_merge_policy == TransferMergePolicy::AVERAGE)
    {
      mooseAssert(state.count, "Average merge should never divide by zero");
      resolved_values /= state.count;
    }
  }
}

#endif
