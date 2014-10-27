
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GLOBAL_SFM_ENGINE_RIG_TIJ_COMPUTATION_H
#define OPENMVG_GLOBAL_SFM_ENGINE_RIG_TIJ_COMPUTATION_H

#include "openMVG/numeric/numeric.h"
#include "openMVG/tracks/tracks.hpp"
#include "software/globalSfM/SfMRigidGlobalEngine.hpp"
#include "software/globalSfM/SfMGlobalEngine_triplet_t_estimator.hpp"

#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"

#undef DYNAMIC
#include "openMVG/bundle_adjustment/problem_data_container.hpp"
#include "software/globalSfM/SfMBundleAdjustmentHelper_tonly.hpp"

#include "openMVG/matching/indexed_sort.hpp"

#include "software/globalSfM/mutexSet.hpp"

namespace openMVG{

// Robust estimation and refinement of a translation and 3D points of an image triplets.
bool estimate_rig_T_triplet(
  const openMVG::tracks::STLMAPTracks & map_tracksCommon,
  const std::map<size_t, std::vector<SIOPointFeature> > & map_feats,
  const std::vector<Mat3> & vec_global_KR_Triplet,
  const Mat3 & K,
  std::vector<Vec3> & vec_tis,
  double & dPrecision, // UpperBound of the precision found by the AContrario estimator
  std::vector<size_t> & vec_inliers,
  const double ThresholdUpperBound, //Threshold used for the trifocal tensor estimation solver used in AContrario Ransac
  const size_t nI,
  const size_t nJ,
  const size_t nK,
  const std::string & sOutDirectory)
{
  bool  bTest=false;

  return bTest;
}

//-- Perform a trifocal estimation of the graph contain in vec_triplets with an
// edge coverage algorithm. It's complexity is sub-linear in term of edges count.
void GlobalRigidReconstructionEngine::computePutativeTranslation_EdgesCoverage(
  const std::map<size_t, Mat3> & map_globalR,
  const std::vector< graphUtils::Triplet > & vec_triplets,
  std::vector<openMVG::relativeInfo > & vec_initialEstimates,
  matching::RigWiseMatches & newpairMatches) const
{
  // The same K matrix is used by all the camera
  const Mat3 _K = Mat3::Identity();

  //-- Prepare global rotations
  std::map<size_t, Mat3> map_global_KR;
  for (std::map<std::size_t, Mat3>::const_iterator iter = map_globalR.begin();
    iter != map_globalR.end(); ++iter)
  {
    map_global_KR[iter->first] = _K * iter->second;
  }

  //-- Prepare tracks count per triplets:
  std::map<size_t, size_t> map_tracksPerTriplets;
#ifdef USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < (int)vec_triplets.size(); ++i)
  {
    const graphUtils::Triplet & triplet = vec_triplets[i];
    const size_t I = triplet.i, J = triplet.j , K = triplet.k;

    RigWiseMatches map_matchesIJK;
    if(_map_Matches_Rig.find(std::make_pair(I,J)) != _map_Matches_Rig.end())
      map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(I,J)));
    else
    if(_map_Matches_Rig.find(std::make_pair(J,I)) != _map_Matches_Rig.end())
      map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(J,I)));

    if(_map_Matches_Rig.find(std::make_pair(I,K)) != _map_Matches_Rig.end())
      map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(I,K)));
    else
    if(_map_Matches_Rig.find(std::make_pair(K,I)) != _map_Matches_Rig.end())
      map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(K,I)));

    if(_map_Matches_Rig.find(std::make_pair(J,K)) != _map_Matches_Rig.end())
      map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(J,K)));
    else
    if(_map_Matches_Rig.find(std::make_pair(K,J)) != _map_Matches_Rig.end())
      map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(K,J)));

    // Compute tracks:
    openMVG::tracks::STLMAPTracks map_tracks;
    TracksBuilder tracksBuilder;
    {
      tracksBuilder.Build(map_matchesIJK);
      tracksBuilder.Filter(3);
      tracksBuilder.ExportToSTL(map_tracks);
    }
#ifdef USE_OPENMP
  #pragma comment omp critical
#endif
    map_tracksPerTriplets[i] = map_tracks.size();
  }
}

} // namespace openMVG

#endif // OPENMVG_GLOBAL_SFM_ENGINE_RIG_TIJ_COMPUTATION_H
