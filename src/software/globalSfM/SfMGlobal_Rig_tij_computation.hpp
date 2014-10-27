
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
  const std::vector<Vec3> & vec_tis,
  std::map<size_t, size_t>  map_RigIdPerImageId,
  const std::vector<openMVG::SfMIO::IntrinsicCameraRigInfo> & vec_intrinsicGroups,
  const double & ThresholdUpperBound)
{
  bool  bTest=false;

  // trifocal tensor evaluation using opengv multi non central relative pose

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

  // compute average focal for trifocal tensor tolerance computation
  double  averageFocal=0.0;

  for(int k(0) ; k < _vec_intrinsicGroups.size(); ++k)
      averageFocal += _vec_intrinsicGroups[k].m_focal ;

  averageFocal /= (double) _vec_intrinsicGroups.size();

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

  typedef std::pair<size_t,size_t> myEdge;

  //-- List all edges
  std::set<myEdge > set_edges;

  for (size_t i = 0; i < vec_triplets.size(); ++i)
  {
    const graphUtils::Triplet & triplet = vec_triplets[i];
    const size_t I = triplet.i, J = triplet.j , K = triplet.k;
    // Add three edges
    set_edges.insert(std::make_pair(std::min(I,J), std::max(I,J)));
    set_edges.insert(std::make_pair(std::min(I,K), std::max(I,K)));
    set_edges.insert(std::make_pair(std::min(J,K), std::max(J,K)));
  }

  // Copy them in vector in order to try to compute them in parallel
  std::vector<myEdge > vec_edges(set_edges.begin(), set_edges.end());

  MutexSet<myEdge> m_mutexSet;

  std::cout << std::endl
    << "Computation of the relative translations over the graph with an edge coverage algorithm" << std::endl;
  #ifdef USE_OPENMP
  #pragma comment omp parallel for schedule(dynamic)
  #endif
  for (int k = 0; k < vec_edges.size(); ++k)
  {
    const myEdge & edge = vec_edges[k];
    //-- If current edge already computed continue
    if (m_mutexSet.isDiscarded(edge) || m_mutexSet.size() == vec_edges.size())
    {
      std::cout << "EDGES WAS PREVIOUSLY COMPUTED" << std::endl;
      continue;
    }

    std::vector<size_t> vec_possibleTriplets;
    // Find the triplet that contain the given edge
    for (size_t i = 0; i < vec_triplets.size(); ++i)
    {
      const graphUtils::Triplet & triplet = vec_triplets[i];
      if (triplet.contain(edge))
      {
        vec_possibleTriplets.push_back(i);
      }
    }

    //-- Sort the triplet according the number of matches they have on their edges
    std::vector<size_t> vec_commonTracksPerTriplets;
    for (size_t i = 0; i < vec_possibleTriplets.size(); ++i)
    {
      vec_commonTracksPerTriplets.push_back(map_tracksPerTriplets[vec_possibleTriplets[i]]);
    }
    //-- If current edge already computed continue
    if (m_mutexSet.isDiscarded(edge))
      continue;

    using namespace indexed_sort;
    std::vector< sort_index_packet_descend < size_t, size_t> > packet_vec(vec_commonTracksPerTriplets.size());
    sort_index_helper(packet_vec, &vec_commonTracksPerTriplets[0]);

    std::vector<size_t> vec_possibleTripletsSorted;
    for (size_t i = 0; i < vec_commonTracksPerTriplets.size(); ++i) {
      vec_possibleTripletsSorted.push_back( vec_possibleTriplets[packet_vec[i].index] );
    }
    vec_possibleTriplets.swap(vec_possibleTripletsSorted);

    // Try to solve the triplets
    // Search the possible triplet:
    for (size_t i = 0; i < vec_possibleTriplets.size(); ++i)
    {
      const graphUtils::Triplet & triplet = vec_triplets[vec_possibleTriplets[i]];
      const size_t I = triplet.i, J = triplet.j , K = triplet.k;
      {
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

        // Select common point:
        STLMAPTracks map_tracksCommon;
        TracksBuilder tracksBuilder;
        {
          tracksBuilder.Build(map_matchesIJK);
          tracksBuilder.Filter(3);
          tracksBuilder.ExportToSTL(map_tracksCommon);
        }

        //--
        // Try to estimate this triplet.
        //--
        // Get rotations:
        std::vector<Mat3> vec_global_KR_Triplet;
        vec_global_KR_Triplet.push_back(map_global_KR[I]);
        vec_global_KR_Triplet.push_back(map_global_KR[J]);
        vec_global_KR_Triplet.push_back(map_global_KR[K]);

        // update precision to have good value for normalized coordinates
        const double ThresholdUpperBound = 0.5 / averageFocal;

        std::vector<Vec3> vec_tis(3);
        std::vector<size_t> vec_inliers;

        if (map_tracksCommon.size() > 50 &&
            estimate_rig_T_triplet(
              map_tracksCommon, _map_feats_normalized, vec_global_KR_Triplet, vec_tis, _map_RigIdPerImageId,
                _vec_intrinsicGroups, ThresholdUpperBound ) ) ;
        {
          std::cout << " triplet estimated."  << endl;
        }
        }
      }
    }
  }
}

} // namespace openMVG

#endif // OPENMVG_GLOBAL_SFM_ENGINE_RIG_TIJ_COMPUTATION_H
