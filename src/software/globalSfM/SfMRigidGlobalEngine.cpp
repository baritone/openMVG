
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/features.hpp"
#include "openMVG/image/image.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching/indexed_sort.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"

#include "software/globalSfM/indexedImageGraph.hpp"
#include "software/globalSfM/indexedImageGraphExport.hpp"
#include "software/globalSfM/SfMRigidGlobalEngine.hpp"
#include "software/globalSfM/SfMGlobal_Rig_tij_computation.hpp"
#include "software/SfM/SfMIOHelper.hpp"
#include "software/SfM/SfMRobust.hpp"
#include "software/SfM/SfMPlyHelper.hpp"

// Rotation averaging
#include "openMVG/multiview/rotation_averaging.hpp"
// Translation averaging
#include "openMVG/linearProgramming/lInfinityCV/global_translations_fromTij.hpp"
#include "openMVG/multiview/translation_averaging_solver.hpp"

// Linear programming solver(s)
#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#ifdef OPENMVG_HAVE_MOSEK
#include "openMVG/linearProgramming/linearProgrammingMOSEK.hpp"
#endif

#include "openMVG/graph/connectedComponent.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"
#include "third_party/stlAddition/stlMap.hpp"
#include "third_party/histogram/histogram.hpp"

#include <opengv/types.hpp>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/NoncentralRelativeAdapter.hpp>
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.hpp>
#include <../test/time_measurement.hpp>

#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

#undef DYNAMIC
#include "openMVG/bundle_adjustment/problem_data_container.hpp"
#include "openMVG/bundle_adjustment/pinhole_ceres_functor.hpp"
#include "openMVG/bundle_adjustment/rig_pinhole_ceres_functor.hpp"
#include "software/globalSfM/SfMBundleAdjustmentHelper_tonly.hpp"

#include "lemon/list_graph.h"
#include <lemon/connectivity.h>

#include "third_party/progress/progress.hpp"
#include "openMVG/system/timer.hpp"

#include <numeric>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <sstream>

#include <stdio.h>
#include <stdlib.h>

using namespace openMVG;
using namespace openMVG::graphUtils;
using namespace svg;
using namespace opengv;

namespace openMVG{

typedef SIOPointFeature FeatureT;
typedef std::vector<FeatureT> featsT;

template<typename T>
void KeepOnlyReferencedElement(
  const std::set<size_t> & Ids,
  T & toFilter)
{
  std::cout << "Must be specialized for your type" << std::endl;
}

// Specialization for Map_RelativeRT
template<>
void KeepOnlyReferencedElement(
  const std::set<size_t> & set_remainingIds,
  GlobalRigidReconstructionEngine::Map_RelativeRT& map_relatives)
{
  GlobalRigidReconstructionEngine::Map_RelativeRT map_relatives_infered;
  for (GlobalRigidReconstructionEngine::Map_RelativeRT::const_iterator
    iter = map_relatives.begin();
    iter != map_relatives.end(); ++iter)
  {
    if (set_remainingIds.find(iter->first.first) != set_remainingIds.end() &&
        set_remainingIds.find(iter->first.second) != set_remainingIds.end())
    {
      map_relatives_infered.insert(*iter);
    }
  }
  map_relatives.swap(map_relatives_infered);
}

// Specialization for RigWiseMatches
template<>
void KeepOnlyReferencedElement(
  const std::set<size_t> & set_remainingIds,
  RigWiseMatches& map_matches)
{
  RigWiseMatches map_matches_Rig_infered;
  for (RigWiseMatches::const_iterator iter = map_matches.begin();
    iter != map_matches.end(); ++iter)
  {
    if (set_remainingIds.find(iter->first.first) != set_remainingIds.end() &&
        set_remainingIds.find(iter->first.second) != set_remainingIds.end())
    {
      map_matches_Rig_infered.insert(*iter);
    }
  }
  map_matches.swap(map_matches_Rig_infered);
}

// Specialization for std::map<size_t,Mat3>
template<>
void KeepOnlyReferencedElement(
  const std::set<size_t> & set_remainingIds,
  std::map<size_t,Mat3>& map_Mat3)
{
  std::map<size_t,Mat3> map_infered;
  for (std::map<size_t,Mat3>::const_iterator iter = map_Mat3.begin();
    iter != map_Mat3.end(); ++iter)
  {
    if (set_remainingIds.find(iter->first) != set_remainingIds.end())
    {
      map_infered.insert(*iter);
    }
  }
  map_Mat3.swap(map_infered);
}

// Specialization for std::vector<openMVG::relativeInfo>
template<>
void KeepOnlyReferencedElement(
  const std::set<size_t> & set_remainingIds,
  std::vector<openMVG::relativeInfo> & map_relativeInfo)
{
  std::vector<openMVG::relativeInfo> map_infered;
  for (std::vector<openMVG::relativeInfo>::const_iterator iter = map_relativeInfo.begin();
    iter != map_relativeInfo.end(); ++iter)
  {
    if (set_remainingIds.find(iter->first.first) != set_remainingIds.end() &&
        set_remainingIds.find(iter->first.second) != set_remainingIds.end())
    {
      map_infered.push_back(*iter);
    }
  }
  map_relativeInfo.swap(map_infered);
}

/// Return imageIds that belongs to the largest bi-edge connected component
template<typename EdgesInterface_T>
std::set<size_t> CleanGraph_Node(
  const EdgesInterface_T & edges,
  const std::vector<std::string> & vec_fileNames,
  const std::string & _sOutDirectory,
  const size_t subCameraNumber)
{
  std::set<size_t> largestBiEdgeCC;

    // Create a graph from pairwise correspondences:
  // - remove not biedge connected component,
  // - keep the largest connected component.

  typedef lemon::ListGraph Graph;
  imageGraph::indexedImageGraph putativeGraph(edges, vec_fileNames, subCameraNumber);

  // Save the graph before cleaning:
  imageGraph::exportToGraphvizData(
    stlplus::create_filespec(_sOutDirectory, "initialGraph"),
    putativeGraph.g);

  // Remove not bi-edge connected edges
  typedef Graph::EdgeMap<bool> EdgeMapAlias;
  EdgeMapAlias cutMap(putativeGraph.g);

  if (lemon::biEdgeConnectedCutEdges(putativeGraph.g, cutMap) > 0)
  {
    // Some edges must be removed because they don't follow the biEdge condition.
    typedef Graph::EdgeIt EdgeIterator;
    EdgeIterator itEdge(putativeGraph.g);
    for (EdgeMapAlias::MapIt it(cutMap); it!=INVALID; ++it, ++itEdge)
    {
      if (*it)
        putativeGraph.g.erase(itEdge); // remove the not bi-edge element
    }
  }

  // Graph is bi-edge connected, but still many connected components can exist
  // Keep only the largest one
  const int connectedComponentCount = lemon::countConnectedComponents(putativeGraph.g);
  std::cout << "\n"
    << "GlobalRigidReconstructionEngine::CleanGraph() :: => connected Component: "
    << connectedComponentCount << std::endl;
  if (connectedComponentCount >= 1)
  {
    // Keep only the largest connected component
    // - list all CC size
    // - if the largest one is meet, keep all the edges that belong to this node

    const std::map<size_t, std::set<lemon::ListGraph::Node> > map_subgraphs = exportGraphToMapSubgraphs(putativeGraph.g);
    size_t count = std::numeric_limits<size_t>::min();
    std::map<size_t, std::set<lemon::ListGraph::Node> >::const_iterator iterLargestCC = map_subgraphs.end();
    for(std::map<size_t, std::set<lemon::ListGraph::Node> >::const_iterator iter = map_subgraphs.begin();
        iter != map_subgraphs.end(); ++iter)
    {
      if (iter->second.size() > count)  {
        count = iter->second.size();
        iterLargestCC = iter;
      }
      std::cout << "Connected component of size: " << iter->second.size() << std::endl;
    }

    //-- Keep only the nodes that are in the largest CC
    for(std::map<size_t, std::set<lemon::ListGraph::Node> >::const_iterator iter = map_subgraphs.begin();
        iter != map_subgraphs.end(); ++iter)
    {
      if (iter == iterLargestCC)
      {
        // list all nodes that belong to the current CC and update the Node Ids list
        const std::set<lemon::ListGraph::Node> & ccSet = iter->second;
        for (std::set<lemon::ListGraph::Node>::const_iterator iter2 = ccSet.begin();
          iter2 != ccSet.end(); ++iter2)
        {
          const size_t Id = (*putativeGraph.map_nodeMapIndex)[*iter2];
          largestBiEdgeCC.insert(Id);

        }
      }
      else
      {
        // remove the edges from the graph
        const std::set<lemon::ListGraph::Node> & ccSet = iter->second;
        for (std::set<lemon::ListGraph::Node>::const_iterator iter2 = ccSet.begin();
          iter2 != ccSet.end(); ++iter2)
        {
          typedef Graph::OutArcIt OutArcIt;
          for (OutArcIt e(putativeGraph.g, *iter2); e!=INVALID; ++e)
          {
            putativeGraph.g.erase(e);
          }
        }
      }
    }
  }

  // Save the graph after cleaning:
  imageGraph::exportToGraphvizData(
    stlplus::create_filespec(_sOutDirectory, "cleanedGraph"),
    putativeGraph.g);

  std::cout << "\n"
    << "Cardinal of nodes: " << lemon::countNodes(putativeGraph.g) << "\n"
    << "Cardinal of edges: " << lemon::countEdges(putativeGraph.g) << std::endl
    << std::endl;

  return largestBiEdgeCC;
}



GlobalRigidReconstructionEngine::GlobalRigidReconstructionEngine(
  const std::string & sImagePath,
  const std::string & sMatchesPath,
  const std::string & sOutDirectory,
  const ERotationAveragingMethod & eRotationAveragingMethod,
  const ETranslationAveragingMethod & eTranslationAveragingMethod,
  bool bHtmlReport)
  : ReconstructionEngine(sImagePath, sMatchesPath, sOutDirectory),
    _eRotationAveragingMethod(eRotationAveragingMethod),
    _eTranslationAveragingMethod(eTranslationAveragingMethod),
    _bRefineIntrinsics(true)
{
  _bHtmlReport = bHtmlReport;
  if (!stlplus::folder_exists(sOutDirectory)) {
    stlplus::folder_create(sOutDirectory);
  }
  if (_bHtmlReport)
  {
    _htmlDocStream = auto_ptr<htmlDocument::htmlDocumentStream>(
      new htmlDocument::htmlDocumentStream("GlobalRigidReconstructionEngine SFM report."));
    _htmlDocStream->pushInfo(
      htmlDocument::htmlMarkup("h1", std::string("Current directory: ") +
      sImagePath));
    _htmlDocStream->pushInfo("<hr>");
  }
}

GlobalRigidReconstructionEngine::~GlobalRigidReconstructionEngine()
{
  ofstream htmlFileStream( string(stlplus::folder_append_separator(_sOutDirectory) +
    "Reconstruction_Report.html").c_str());
  htmlFileStream << _htmlDocStream->getDoc();
}

void GlobalRigidReconstructionEngine::rotationInference(
  Map_RelativeRT & map_relatives)
{
  std::cout
        << "---------------\n"
        << "-- INFERENCE on " << _map_Matches_Rig.size() << " EGs count.\n"
        << "---------------" << std::endl
        << " /!\\  /!\\  /!\\  /!\\  /!\\  /!\\  /!\\ \n"
        << "--- ITERATED BAYESIAN INFERENCE IS NOT RELEASED, SEE C.ZACH CODE FOR MORE INFORMATION" << std::endl
        << " /!\\  /!\\  /!\\  /!\\  /!\\  /!\\  /!\\ \n" << std::endl
        << " A simple inference scheme is used here:" << std::endl
        << "\t only the relative error composition to identity on cycle of length 3 is used." << std::endl;

  //-------------------
  // Triplet inference (test over the composition error)
  //-------------------
  std::vector< graphUtils::Triplet > vec_triplets;
  tripletListing(vec_triplets, _map_ImagesIdPerRigId[0].size());
  //-- Rejection triplet that are 'not' identity rotation (error to identity > 2°)
  tripletRotationRejection(vec_triplets, map_relatives);
}

/// Association of Ids to a contiguous set of Ids
template<typename T>
void reindex(
  const std::vector< std::pair<T,T> > & vec_pairs,
  std::map<T, T> & _reindexForward,
  std::map<T, T> & _reindexBackward)
{
  // get an unique set of Ids
  std::set<size_t> _uniqueId;
  for(typename std::vector< typename std::pair<T,T> >::const_iterator iter = vec_pairs.begin();
        iter != vec_pairs.end(); ++iter)
  {
    _uniqueId.insert(iter->first);
    _uniqueId.insert(iter->second);
  }

  // Build the Forward and Backward mapping
  for(typename std::vector< typename std::pair<T,T> >::const_iterator iter = vec_pairs.begin();
        iter != vec_pairs.end(); ++iter)
  {
    if (_reindexForward.find(iter->first) == _reindexForward.end())
    {
      const size_t dist = std::distance(_uniqueId.begin(), _uniqueId.find(iter->first));
      _reindexForward[iter->first] = dist;
      _reindexBackward[dist] = iter->first;
    }
    if (_reindexForward.find(iter->second) == _reindexForward.end())
    {
      const size_t dist = std::distance(_uniqueId.begin(), _uniqueId.find(iter->second));
      _reindexForward[iter->second] = dist;
      _reindexBackward[dist] = iter->second;
    }
  }
}

bool GlobalRigidReconstructionEngine::computeGlobalRotations(
  ERotationAveragingMethod eRotationAveragingMethod,
  const Map_RelativeRT & map_relatives,
  std::map<size_t, Mat3> & map_globalR) const
{
  // Build relative information for only the largest considered Connected Component
  // - it requires to change the camera indexes, because RotationAveraging is working only with
  //   index ranging in [0 - nbCam]

  std::map<size_t, size_t> _reindexForward, _reindexBackward;

  std::vector< std::pair<size_t,size_t> > vec_pairs;
  for(Map_RelativeRT::const_iterator iter = map_relatives.begin();
        iter != map_relatives.end(); ++iter)
  {
    const openMVG::relativeInfo & rel = *iter;
    vec_pairs.push_back(rel.first);
  }

  reindex(vec_pairs, _reindexForward, _reindexBackward);

  //- A. weight computation
  //- B. solve global rotation computation

  //- A. weight computation: for each pair w = min(1, (#PairMatches/Median(#PairsMatches)))
  std::vector<double> vec_relativeRotWeight;
  {
    vec_relativeRotWeight.reserve(map_relatives.size());
    {
      //-- Compute the median number of matches
      std::vector<double> vec_count;
      for(Map_RelativeRT::const_iterator iter = map_relatives.begin();
        iter != map_relatives.end(); ++iter)
      {
        const openMVG::relativeInfo & rel = *iter;
        // Find the number of support point for this pair
        RigWiseMatches::const_iterator iterMatches = _map_Matches_Rig.find(rel.first);
        if (iterMatches != _map_Matches_Rig.end())
        {
          size_t  supportPointNumber = 0;
          for(PairWiseMatches::const_iterator iterRig = iterMatches->second.begin() ;
              iterRig != iterMatches->second.end(); ++iterRig )
          {
              supportPointNumber += iterRig->second.size();
          }

          vec_count.push_back(supportPointNumber);
        }
      }
      const float thTrustPair = (std::accumulate(vec_count.begin(), vec_count.end(), 0.0f) / vec_count.size()) / 2.0;

      for(Map_RelativeRT::const_iterator iter = map_relatives.begin();
        iter != map_relatives.end(); ++iter)
      {
        const openMVG::relativeInfo & rel = *iter;
        float weight = 1.f; // the relative rotation correspondence point support
        RigWiseMatches::const_iterator iterMatches = _map_Matches_Rig.find(rel.first);
        if (iterMatches != _map_Matches_Rig.end())
        {
          size_t  supportPointNumber = 0;
          for(PairWiseMatches::const_iterator iterRig = iterMatches->second.begin() ;
              iterRig != iterMatches->second.end(); ++iterRig )
          {
              supportPointNumber += iterRig->second.size();
          }

          weight = std::min((float)supportPointNumber/thTrustPair, 1.f);
          vec_relativeRotWeight.push_back(weight);
        }
      }
    }
  }

  using namespace openMVG::rotation_averaging;
  // Setup input data for global rotation computation
  std::vector<RelRotationData> vec_relativeRotEstimate;
  std::vector<double>::const_iterator iterW = vec_relativeRotWeight.begin();
  for(Map_RelativeRT::const_iterator iter = map_relatives.begin();
    iter != map_relatives.end(); ++iter)
  {
    const openMVG::relativeInfo & rel = *iter;
    RigWiseMatches::const_iterator iterMatches = _map_Matches_Rig.find(rel.first);
    if (iterMatches != _map_Matches_Rig.end())
    {
      vec_relativeRotEstimate.push_back(RelRotationData(
        _reindexForward[rel.first.first],
        _reindexForward[rel.first.second],
        rel.second.first, *iterW));
      ++iterW;
    }
  }

  //- B. solve global rotation computation
  bool bSuccess = false;
  std::vector<Mat3> vec_globalR(_reindexForward.size());
  switch(eRotationAveragingMethod)
  {
    case ROTATION_AVERAGING_L2:
    {
      //- Solve the global rotation estimation problem:
      bSuccess = rotation_averaging::l2::L2RotationAveraging(
        _reindexForward.size(),
        vec_relativeRotEstimate,
        vec_globalR);
    }
    break;
    case ROTATION_AVERAGING_L1:
    {
      using namespace openMVG::rotation_averaging::l1;

      //- Solve the global rotation estimation problem:
      const size_t nMainViewID = 0;
      std::vector<bool> vec_inliers;
      bSuccess = rotation_averaging::l1::GlobalRotationsRobust(
        vec_relativeRotEstimate, vec_globalR, nMainViewID, 0.0f, &vec_inliers);

      std::cout << "\ninliers: " << std::endl;
      std::copy(vec_inliers.begin(), vec_inliers.end(), ostream_iterator<bool>(std::cout, " "));
      std::cout << std::endl;
    }
    break;
    default:
    std::cerr << "Unknown rotation averaging method: " << (int) eRotationAveragingMethod << std::endl;
  }

  if (bSuccess)
  {
    //-- Setup the averaged rotation
    for (int i = 0; i < vec_globalR.size(); ++i)  {
      map_globalR[_reindexBackward[i]] = vec_globalR[i];
    }
  }
  else{
    std::cerr << "Global rotation solving failed." << std::endl;
  }

  return bSuccess;
}

bool GlobalRigidReconstructionEngine::Process()
{

//---------------------------------
//-- Global Calibration -----------
//---------------------------------

  //-------------------
  // Load data
  //-------------------

  if(!ReadInputData())  {
    std::cout << "\nError while parsing input data" << std::endl;
    return false;
  }

  //-------------------
  // check input data
  //-------------------

  if(!InputDataIsCorrect())  {
    std::cout << "\nError in your input data" << std::endl;
    return false;
  }

  //--------------------
  // compute liste of pairwise matches per rig
  //--------------------
  ComputeMapMatchesRig();

  //-- Export input graph
  {
    typedef lemon::ListGraph Graph;
    imageGraph::indexedImageGraph putativeGraph(_map_Matches_Rig, _vec_fileNames, _map_ImagesIdPerRigId[0].size());

    // Save the graph before cleaning:
    imageGraph::exportToGraphvizData(
      stlplus::create_filespec(_sOutDirectory, "input_graph"),
      putativeGraph.g);
  }


  openMVG::Timer total_reconstruction_timer;

  //-------------------
  // Only keep the largest biedge connected subgraph
  //-------------------
  {
    const std::set<size_t> set_remainingIds = CleanGraph_Node(_map_Matches_Rig, _vec_fileNames, _sOutDirectory, _map_ImagesIdPerRigId[0].size());
    if(set_remainingIds.empty())
    {
      std::cout << "Invalid input image graph for global SfM" << std::endl;
    return false;
    }
    KeepOnlyReferencedElement(set_remainingIds, _map_Matches_Rig);
  }

  //-------------------
  // Compute relative R|t
  //-------------------

  Map_RelativeRT map_relatives;
  {
    ComputeRelativeRt(map_relatives);
  }

  //-------------------
  // Rotation inference
  //-------------------

  {
    openMVG::Timer timer_Inference;

    rotationInference(map_relatives);

    //-------------------
    // keep the largest biedge connected subgraph
    //-------------------
    const std::set<size_t> set_remainingIds = CleanGraph_Node(_map_Matches_Rig, _vec_fileNames, _sOutDirectory, _map_ImagesIdPerRigId[0].size());
    if(set_remainingIds.empty())
      return false;

    const double time_Inference = timer_Inference.elapsed();

    // Clean Map_RelativeRT and relative matches
    KeepOnlyReferencedElement(set_remainingIds, map_relatives);
    KeepOnlyReferencedElement(set_remainingIds, _map_Matches_Rig);

    //---------------------------------------
    //-- Export geometric filtered matches
    //---------------------------------------
  //  std::ofstream file (string(_sMatchesPath + "/matches.filtered.txt").c_str());
  //  if (file.is_open())
  //    PairedIndMatchToStream(_map_Matches_E, file); // need to modifiy that ?
  //  file.close();

    std::cout << "\n Remaining cameras after inference filter: \n"
      << set_remainingIds.size() << " from a total of " << _vec_fileNames.size() / _vec_intrinsicGroups.size() << std::endl;

    //-- Export statistics about the rotation inference step:
    if (_bHtmlReport)
    {
      using namespace htmlDocument;
      std::ostringstream os;
      os << "Rotation inference.";
      _htmlDocStream->pushInfo("<hr>");
      _htmlDocStream->pushInfo(htmlMarkup("h1",os.str()));

      os.str("");
      os << "-------------------------------" << "<br>"
        << "-- #Camera count: " << set_remainingIds.size() << " remains "
        << "-- from " <<_vec_fileNames.size() << " input images.<br>"
        << "-- timing: " << time_Inference << " second <br>"
        << "-------------------------------" << "<br>";
      _htmlDocStream->pushInfo(os.str());
    }
  }

  //----------------------------
  // Rotation averaging
  //----------------------------

  std::map<std::size_t, Mat3> map_globalR;
  {
    std::set<size_t> set_indeximage;
    for (RigWiseMatches::const_iterator
         iter = _map_Matches_Rig.begin();
         iter != _map_Matches_Rig.end();
         ++iter)
    {
      const size_t I = iter->first.first;
      const size_t J = iter->first.second;
      set_indeximage.insert(I);
      set_indeximage.insert(J);
    }

    std::cout << "\n-------------------------------" << "\n"
      << " Global rotations computation: " << "\n"
      << "   - Ready to compute " << set_indeximage.size() << " global rotations." << "\n"
      << "     from " << map_relatives.size() << " relative rotations\n" << std::endl;

    if (!computeGlobalRotations(
            _eRotationAveragingMethod,
            map_relatives,
            map_globalR))
    {
      std::cerr << "Failed to compute the global rotations." << std::endl;
      return false;
    }
  }

  //-------------------
  // Relative translations estimation (Triplet based translation computation)
  //-------------------
  std::vector<openMVG::relativeInfo > vec_initialRijTijEstimates;
  RigWiseMatches newpairMatches;
  {
    std::cout << "\n-------------------------------" << "\n"
      << " Relative translations computation: " << "\n"
      << "-------------------------------" << std::endl;

    // List putative triplets
    std::vector< graphUtils::Triplet > vec_triplets;
    tripletListing(vec_triplets, _map_ImagesIdPerRigId[0].size());

    // Compute putative translations with an edge coverage algorithm

    openMVG::Timer timerLP_triplet;

    bool  bComputeTrifocal=false;
    if(bComputeTrifocal){
       computePutativeTranslation_EdgesCoverage(map_globalR, vec_triplets, vec_initialRijTijEstimates, newpairMatches);
    }
    else
    {
       vec_initialRijTijEstimates.clear();
       for(Map_RelativeRT::const_iterator  iter = map_relatives.begin();
              iter != map_relatives.end(); ++iter )
       {
           // compute relative rotation of rig pair
           const size_t I = iter->first.first;
           const size_t J = iter->first.second;

           const Mat3  RI = map_globalR.find(I)->second;
           const Mat3  RJ = map_globalR.find(J)->second;

           const Mat3  RIJ = RJ * RI.transpose() ;

           vec_initialRijTijEstimates.push_back( std::make_pair(iter->first , std::make_pair(RIJ, iter->second.second) ) ) ;

       }

    }
    const double timeLP_triplet = timerLP_triplet.elapsed();
    std::cout << "TRIPLET COVERAGE TIMING: " << timeLP_triplet << " seconds" << std::endl;

    //-- Export triplet statistics:
    if (_bHtmlReport)
    {
      using namespace htmlDocument;
      std::ostringstream os;
      os << "Triplet statistics.";
      _htmlDocStream->pushInfo("<hr>");
      _htmlDocStream->pushInfo(htmlMarkup("h1",os.str()));

      os.str("");
      os << "-------------------------------" << "<br>"
        << "-- #Effective translations estimates: " << vec_initialRijTijEstimates.size()/3
        << " from " <<vec_triplets.size() << " triplets.<br>"
        << "-- resulting in " <<vec_initialRijTijEstimates.size() << " translation estimation.<br>"
        << "-- timing to obtain the relative translations: " << timeLP_triplet << " seconds.<br>"
        << "-------------------------------" << "<br>";
      _htmlDocStream->pushInfo(os.str());
    }
  }

  //--Check the relative translation graph:
  //--> Consider only the connected component compound by the translation graph
  //-- Robust translation estimation can perform inference and remove some bad conditioned triplets

  {
    // Build the list of Pairs used by the translations
    std::vector<std::pair<size_t, size_t> > map_pairs_tij;
    for(size_t i = 0; i < vec_initialRijTijEstimates.size(); ++i)
    {
      const openMVG::relativeInfo & rel = vec_initialRijTijEstimates[i];
      map_pairs_tij.push_back(std::make_pair(rel.first.first,rel.first.second));
    }

    const std::set<size_t> set_representedImageIndex = CleanGraph_Node(map_pairs_tij, _vec_fileNames, _sOutDirectory, 1.0);

    std::cout << "\n\n"
      << "We targeting to estimates: " << map_globalR.size()
      << " and we have estimation for: " << set_representedImageIndex.size() << " rigs " << std::endl;

    //-- Clean global rotations that are not in the TRIPLET GRAPH
    KeepOnlyReferencedElement(set_representedImageIndex, map_globalR);
    KeepOnlyReferencedElement(set_representedImageIndex, vec_initialRijTijEstimates);
    KeepOnlyReferencedElement(set_representedImageIndex, newpairMatches);
    // clean _map_matches_E?

    std::cout << "\nRemaining rigs after inference filter: \n"
      << map_globalR.size() << " from a total of " << _vec_fileNames.size() / _vec_intrinsicGroups.size() << std::endl;
  }

  //-------------------
  //-- GLOBAL TRANSLATIONS ESTIMATION from initial triplets t_ij guess
  //-------------------

  {
    const size_t iNRigs = map_globalR.size();
    const size_t iNview = iNRigs * _vec_intrinsicGroups.size() ;

    std::cout << "\n-------------------------------" << "\n"
      << " Global translations computation: " << "\n"
      << "   - Ready to compute " << iNRigs << " global translations." << "\n"
      << "     from " << vec_initialRijTijEstimates.size() << " relative translations\n" << std::endl;

    //-- Update initial estimates in range [0->Nrigs]
    std::map<size_t, size_t> _reindexForward, _reindexBackward;

    std::vector< std::pair<size_t,size_t> > vec_pairs;
        for(size_t i = 0; i < vec_initialRijTijEstimates.size(); ++i)
    {
      const openMVG::relativeInfo & rel = vec_initialRijTijEstimates[i];
      std::pair<size_t,size_t> newPair(rel.first.first, rel.first.second);
      vec_pairs.push_back(newPair);
    }

    reindex( vec_pairs, _reindexForward, _reindexBackward);

    //-- Update initial estimates in range [0->Nrigs]
    for(size_t i = 0; i < vec_initialRijTijEstimates.size(); ++i)
    {
      openMVG::relativeInfo & rel = vec_initialRijTijEstimates[i];
      std::pair<size_t,size_t> newPair(
          _reindexForward[rel.first.first],
          _reindexForward[rel.first.second]);
      rel.first = newPair;
    }


    openMVG::Timer timerLP_translation;

    switch(_eTranslationAveragingMethod)
    {
      case TRANSLATION_AVERAGING_L1:
      {
        double gamma = -1.0;
        std::vector<double> vec_solution;
        {
          vec_solution.resize(iNRigs*3 + vec_initialRijTijEstimates.size() + 1);
          using namespace openMVG::linearProgramming;
          #ifdef OPENMVG_HAVE_MOSEK
            MOSEK_SolveWrapper solverLP(vec_solution.size());
          #else
            OSI_CLP_SolverWrapper solverLP(vec_solution.size());
          #endif

          lInfinityCV::Tifromtij_ConstraintBuilder cstBuilder(vec_initialRijTijEstimates);

          LP_Constraints_Sparse constraint;
          //-- Setup constraint and solver
          cstBuilder.Build(constraint);
          solverLP.setup(constraint);
          //--
          // Solving
          const bool bFeasible = solverLP.solve();
          std::cout << " \n Feasibility " << bFeasible << std::endl;
          //--
          if (bFeasible)  {
            solverLP.getSolution(vec_solution);
            gamma = vec_solution[vec_solution.size()-1];
          }
          else  {
            std::cerr << "Compute global translations: failed" << std::endl;
            return false;
          }
        }

        const double timeLP_translation = timerLP_translation.elapsed();
        //-- Export triplet statistics:
        if (_bHtmlReport)
        {
          using namespace htmlDocument;
          std::ostringstream os;
          os << "Translation fusion statistics.";
          _htmlDocStream->pushInfo("<hr>");
          _htmlDocStream->pushInfo(htmlMarkup("h1",os.str()));

          os.str("");
          os << "-------------------------------" << "<br>"
            << "-- #relative estimates: " << vec_initialRijTijEstimates.size()
            << " converge with gamma: " << gamma << ".<br>"
            << " timing (s): " << timeLP_translation << ".<br>"
            << "-------------------------------" << "<br>";
          _htmlDocStream->pushInfo(os.str());
        }

        std::cout << "Found solution:\n";
        std::copy(vec_solution.begin(), vec_solution.end(), std::ostream_iterator<double>(std::cout, " "));

        std::vector<double> vec_RigTranslation(iNRigs*3,0);
        std::copy(&vec_solution[0], &vec_solution[iNRigs*3], &vec_RigTranslation[0]);

        std::vector<double> vec_rigRelLambdas(&vec_solution[iNRigs*3], &vec_solution[iNRigs*3 + vec_initialRijTijEstimates.size()/3]);
        std::cout << "\nrig position: " << std::endl;
        std::copy(vec_RigTranslation.begin(), vec_RigTranslation.end(), std::ostream_iterator<double>(std::cout, " "));
        std::cout << "\nrig Lambdas: " << std::endl;
        std::copy(vec_rigRelLambdas.begin(), vec_rigRelLambdas.end(), std::ostream_iterator<double>(std::cout, " "));

        // Build a Pinhole camera for each considered rigs
        std::vector<Vec3>  vec_C;
        for (size_t i = 0; i < iNRigs; ++i)
        {
          const size_t rigNum    = _reindexBackward[i];
          const std::vector<size_t> ImageList = _map_ImagesIdPerRigId[rigNum];

          // extract rig pose
          const Mat3 & Ri = map_globalR[rigNum];
          const Vec3 Rigt(vec_RigTranslation[rigNum*3], vec_RigTranslation[rigNum*3+1], vec_RigTranslation[rigNum*3+2]);

          for(size_t j= 0 ; j < ImageList.size(); ++j)
          {
            const size_t I = ImageList[j];
            const size_t subCamId  = _map_IntrinsicIdPerImageId.find(I)->second;

            // compute camera rotation
            const Mat3 & Rcam = _vec_intrinsicGroups[subCamId].m_R;
            const Mat3 & R = Rcam * Ri;

            // compute camera translation
            const Vec3 tCam = -Rcam * _vec_intrinsicGroups[subCamId].m_rigC ;
            const Vec3 t = Rcam * Rigt + tCam;

            const Mat3 & _K = _vec_intrinsicGroups[subCamId].m_K;   // The same K matrix is used by all the camera
            _map_camera[I] = PinholeCamera(_K, R, t);

            //-- Export camera center
            vec_C.push_back(_map_camera[I]._C);
          }

          //update rig map
          _map_rig[rigNum] = std::make_pair(Ri, Rigt);
        }
        plyHelper::exportToPly(vec_C, stlplus::create_filespec(_sOutDirectory, "cameraPath", "ply"));
      }
      break;

      case TRANSLATION_AVERAGING_L2:
      {
          std::cerr << " L2 translation averaging not yet supported " << endl;
          return true;
      }
      break;
    }
  }

  //-------------------
  //-- Initial triangulation of the scene from the computed global motions
  //-------------------
  {
    // keep tracks that are seen in remaining rigs
    RigWiseMatches newpairMatches;
    for(RigWiseMatches::const_iterator iter = _map_Matches_Rig.begin();
          iter != _map_Matches_Rig.end(); ++iter )
    {
      const size_t RigOne = iter->first.first;
      const size_t RigTwo = iter->first.second;

      if( _map_rig.find(RigOne) != _map_rig.end() && _map_rig.find(RigTwo) != _map_rig.end() )
        newpairMatches[iter->first] = iter->second;
    }

    // Compute tracks
    TracksBuilder tracksBuilder;
    {
      tracksBuilder.Build(newpairMatches);
      tracksBuilder.Filter(3);
      tracksBuilder.ExportToSTL(_map_selectedTracks);
    }

    std::vector<double> vec_residuals;

    // Triangulation of all the tracks
    _vec_allScenes.resize(_map_selectedTracks.size());
    {
      std::vector<double> vec_residuals;
      vec_residuals.reserve(_map_selectedTracks.size());

      C_Progress_display my_progress_bar_triangulation( _map_selectedTracks.size(),
      std::cout, "\n\n Initial triangulation:\n");

#ifdef USE_OPENMP
      #pragma omp parallel for schedule(dynamic)
#endif

      for (int idx = 0; idx < _map_selectedTracks.size(); ++idx)
      {
          STLMAPTracks::const_iterator iterTracks = _map_selectedTracks.begin();
          std::advance(iterTracks, idx);

          const submapTrack & subTrack = iterTracks->second;

          // Look to the features required for the triangulation task
          Triangulation trianObj;
          for (submapTrack::const_iterator iterSubTrack = subTrack.begin(); iterSubTrack != subTrack.end(); ++iterSubTrack)
          {
            const size_t imaIndex = iterSubTrack->first;
            const size_t featIndex = iterSubTrack->second;
            const SIOPointFeature & pt = _map_feats[imaIndex][featIndex];

            // Build the P matrix
            trianObj.add(_map_camera[imaIndex]._P, pt.coords().cast<double>());
          }

          // Compute the 3D point and keep point index with negative depth
          const Vec3 Xs = trianObj.compute();
          _vec_allScenes[idx] = Xs;

#ifdef USE_OPENMP
#pragma omp critical
#endif
        //-- Compute residual over all the projections
        {
          for (submapTrack::const_iterator iterSubTrack = subTrack.begin(); iterSubTrack != subTrack.end(); ++iterSubTrack) {
            const size_t imaIndex = iterSubTrack->first;
            const size_t featIndex = iterSubTrack->second;
            const SIOPointFeature & pt = _map_feats[imaIndex][featIndex];
            vec_residuals.push_back(_map_camera[imaIndex].Residual(Xs, pt.coords().cast<double>()));
            // no ordering in vec_residuals since there is parallelism
          }
          ++my_progress_bar_triangulation;
        }
      }

      {
        // Display some statistics of reprojection errors
        std::cout << "\n\nResidual statistics:\n" << std::endl;
        minMaxMeanMedian<double>(vec_residuals.begin(), vec_residuals.end());
        double min, max, mean, median;
        minMaxMeanMedian<double>(vec_residuals.begin(), vec_residuals.end(), min, max, mean, median);

        Histogram<float> histo(0.f, *max_element(vec_residuals.begin(),vec_residuals.end())*1.1f);
        histo.Add(vec_residuals.begin(), vec_residuals.end());
        std::cout << std::endl << "Residual Error pixels: " << std::endl << histo.ToString() << std::endl;

        // Histogram between 0 and 10 pixels
        {
          std::cout << "\n Histogram between 0 and 10 pixels: \n";
          Histogram<float> histo(0.f, 10.f, 20);
          histo.Add(vec_residuals.begin(), vec_residuals.end());
          std::cout << std::endl << "Residual Error pixels: " << std::endl << histo.ToString() << std::endl;
        }

        //-- Export initial triangulation statistics
        if (_bHtmlReport)
        {
          using namespace htmlDocument;
          std::ostringstream os;
          os << "Initial triangulation statistics.";
          _htmlDocStream->pushInfo("<hr>");
          _htmlDocStream->pushInfo(htmlMarkup("h1",os.str()));

          os.str("");
          os << "-------------------------------" << "<br>"
            << "-- #tracks: " << _map_selectedTracks.size() << ".<br>"
            << "-- #observation: " << vec_residuals.size() << ".<br>"
            << "-- residual mean (RMSE): " << std::sqrt(mean) << ".<br>"
            << "-------------------------------" << "<br>";
          _htmlDocStream->pushInfo(os.str());
        }
      }
    }

    plyHelper::exportToPly(_vec_allScenes, stlplus::create_filespec(_sOutDirectory, "raw_pointCloud_LP", "ply"));
  }

  //-------------------
  //-- Bundle Adjustment on translation and structure
  //-------------------

  // Refine only Structure and translations
  bundleAdjustment(_map_rig, _map_camera, _vec_allScenes, _map_selectedTracks, false, true, false);
  plyHelper::exportToPly(_vec_allScenes, stlplus::create_filespec(_sOutDirectory, "raw_pointCloud_BA_T_Xi", "ply"));

  // Refine Structure, rotations and translations
  bundleAdjustment(_map_rig, _map_camera, _vec_allScenes, _map_selectedTracks, true, true, false);
  plyHelper::exportToPly(_vec_allScenes, stlplus::create_filespec(_sOutDirectory, "raw_pointCloud_BA_RT_Xi", "ply"));

  if (_bRefineIntrinsics)
  {
    // Refine Structure, rotations, translations and intrinsics
    bundleAdjustment(_map_rig, _map_camera, _vec_allScenes, _map_selectedTracks, true, true, true);
    plyHelper::exportToPly(_vec_allScenes, stlplus::create_filespec(_sOutDirectory, "raw_pointCloud_BA_KRT_Xi", "ply"));
  }

  //-- Export statistics about the global process
  if (_bHtmlReport)
  {
    using namespace htmlDocument;
    std::ostringstream os;
    os << "Global calibration summary triangulation statistics.";
    _htmlDocStream->pushInfo("<hr>");
    _htmlDocStream->pushInfo(htmlMarkup("h1",os.str()));

    os.str("");
    os << "-------------------------------" << "<br>"
      << "-- Have calibrated: " << _map_camera.size() << " from "
      << _vec_fileNames.size() << " input images.<br>"
      << "-- The scene contains " << _map_selectedTracks.size() << " 3D points.<br>"
      << "-- Total reconstruction time (Inference, global rot, translation fusion, Ba1, Ba2, Ba3): "
      << total_reconstruction_timer.elapsed() << " seconds.<br>"
      << "-------------------------------" << "<br>";
    _htmlDocStream->pushInfo(os.str());
  }

  std::cout << std::endl
    << "-------------------------------" << "\n"
    << "-- Have calibrated: " << _map_camera.size() << " from "
    << _vec_fileNames.size() << " input images.\n"
    << "-- The scene contains " << _map_selectedTracks.size() << " 3D points.\n"
    << "-- Total reconstruction time (Inference, global rot, translation fusion, Ba1, Ba2, Ba3): "
    << total_reconstruction_timer.elapsed() << " seconds.\n"
    << "Relative rotations time was excluded\n"
    << "-------------------------------" << std::endl;

  //-- Export the scene (cameras and structures) to the SfM Data container
  {
    // Cameras
    for (std::map<size_t,PinholeCamera >::const_iterator iter = _map_camera.begin();
      iter != _map_camera.end();  ++iter)
    {
      const PinholeCamera & cam = iter->second;
      _reconstructorData.map_Camera[iter->first] = BrownPinholeCamera(cam._P);
      _reconstructorData.set_imagedId.insert(iter->first);
    }

    // Structure
    size_t i = 0;
    for (std::vector<Vec3>::const_iterator iter = _vec_allScenes.begin();
      iter != _vec_allScenes.end();
      ++iter, ++i)
    {
      const Vec3 & pt3D = *iter;
      _reconstructorData.map_3DPoints[i] = pt3D;
      _reconstructorData.set_trackId.insert(i);
    }
  }
  return true;
}

bool testIntrinsicsEquality(
  SfMIO::IntrinsicCameraRigInfo const &ci1,
  SfMIO::IntrinsicCameraRigInfo const &ci2)
{
  return ci1.m_K == ci2.m_K;
}

bool GlobalRigidReconstructionEngine::ReadInputData()
{
  if (!stlplus::is_folder(_sImagePath) ||
    !stlplus::is_folder(_sMatchesPath) ||
    !stlplus::is_folder(_sOutDirectory))
  {
    std::cerr << std::endl
      << "One of the required directory is not a valid directory" << std::endl;
    return false;
  }

  // a. Read images names
  std::string sListsFile = stlplus::create_filespec(_sMatchesPath,"lists","txt");
  std::string sComputedMatchesFile_E = stlplus::create_filespec(_sMatchesPath,"matches.e","txt");
  if (!stlplus::is_file(sListsFile)||
    !stlplus::is_file(sComputedMatchesFile_E) )
  {
    std::cerr << std::endl
      << "One of the input required file is not a present (lists.txt, matches.e.txt)" << std::endl;
    return false;
  }

  // a. Read images names
  {
    if (!openMVG::SfMIO::loadImageList( _vec_camImageNames,
                                        _vec_intrinsicGroups,
                                        sListsFile) )
    {
      std::cerr << "\nEmpty image list." << std::endl;
      return false;
    }
    else
    {
      // Find to which intrinsic groups each image belong
      for (std::vector<openMVG::SfMIO::CameraRigInfo>::const_iterator iter = _vec_camImageNames.begin();
        iter != _vec_camImageNames.end(); ++iter)
      {
        const openMVG::SfMIO::CameraRigInfo & camInfo = *iter;

        // Find the index of the camera
        const size_t idx = std::distance((std::vector<openMVG::SfMIO::CameraRigInfo>::const_iterator)_vec_camImageNames.begin(), iter);

        // to which intrinsic group each image belongs
        _map_IntrinsicIdPerImageId[idx] = camInfo.m_intrinsicId;

        // to which rigid rig each image belongs
        _map_RigIdPerImageId[idx]       = camInfo.m_rigId;
      }

      for (size_t i = 0; i < _vec_camImageNames.size(); ++i)
      {
        _vec_fileNames.push_back(_vec_camImageNames[i].m_sImageName);
      }
    }
  }

  // b. Read matches (Essential)
  if (!matching::PairedIndMatchImport(sComputedMatchesFile_E, _map_Matches_E)) {
    std::cerr<< "Unable to read the Essential matrix matches" << std::endl;
    return false;
  }

  // Read features:
  for (size_t i = 0; i < _vec_fileNames.size(); ++i)  {
    const size_t camIndex = i;
    if (!loadFeatsFromFile(
      stlplus::create_filespec(_sMatchesPath, stlplus::basename_part(_vec_fileNames[camIndex]), ".feat"),
      _map_feats[camIndex])) {
      std::cerr << "Bad reading of feature files" << std::endl;
      return false;
    }
  }

  // Normalize features:
  for (size_t i = 0; i < _vec_fileNames.size(); ++i)  {
    const size_t camIndex = i;
    const size_t intrinsicId = _map_IntrinsicIdPerImageId[camIndex];

    // load camera matrix
    const Mat3 _K    = _vec_intrinsicGroups[intrinsicId].m_K;
    const Mat3 _Kinv = _K.inverse();

    // array containing normalized feature of image
    std::vector<SIOPointFeature>  normalizedFeatureImageI;
    normalizedFeatureImageI.clear();

    //normalize features
    for( size_t j=0; j < _map_feats[camIndex].size(); ++j)
    {
      const Vec3 x(_map_feats[camIndex][j].x(), _map_feats[camIndex][j].y(), 1.0);
      const Vec3 xBearingVector = _Kinv * x;
      const Vec2 xBearingVectorNormalized = xBearingVector.head(2) / xBearingVector(2);

      const float scale = _map_feats[camIndex][j].scale();
      const float orientation = _map_feats[camIndex][j].orientation();

      const SIOPointFeature  normalized( xBearingVectorNormalized(0), xBearingVectorNormalized(1), scale, orientation);

      normalizedFeatureImageI.push_back(normalized);
    }

    _map_feats_normalized.insert(std::make_pair(camIndex,normalizedFeatureImageI) );

  }

  return true;
}

bool GlobalRigidReconstructionEngine::InputDataIsCorrect()
{
  // check if rotation matrices are rotation matrices
  for(size_t i=0; i < _vec_intrinsicGroups.size() ; ++i )
  {
      Mat3  R = _vec_intrinsicGroups[i].m_R;

      // evaluate difference between R.R^t and I_3
      Mat3  RRt = R.transpose() * R ;
      Mat3  D   = Mat3::Identity() - RRt;

      // R is a rotation matrix if |R| = + 1.0 and R.R^t = I_3
      if( fabs(R.determinant()-1.0) > 1.0e-5 || D.squaredNorm() > 1.0e-8 )
      {
        std::cerr << "Error : Input Rotation Matrix is not a rotation matrix \n";
        return false;
      }
  }

  // generate list of image per rig
  for (std::map<size_t, size_t>::const_iterator iter = _map_RigIdPerImageId.begin();
             iter != _map_RigIdPerImageId.end(); ++ iter)
  {
      const size_t idx = iter->first;
      _map_ImagesIdPerRigId[_map_RigIdPerImageId[idx]].push_back(idx);
  }

  // check that each rig got the same number of subCameras
  size_t  rigSize=0;
  for (std::map<size_t, std::vector<size_t> >::const_iterator iter = _map_ImagesIdPerRigId.begin();
             iter != _map_ImagesIdPerRigId.end(); ++ iter)
  {
      const size_t rigSubCamNumber = iter->second.size();
      if( rigSize == 0)
          rigSize = rigSubCamNumber;

      if( rigSubCamNumber != _vec_intrinsicGroups.size() || rigSubCamNumber != rigSize )
      {
        std::cerr << "Error : The rig does not always have the same number of subcameras or there is more than one rig \n";
        return false;
      }
  }

  return true;
}

void GlobalRigidReconstructionEngine::ComputeMapMatchesRig()
{
  for (PairWiseMatches::const_iterator iter = _map_Matches_E.begin();
    iter != _map_Matches_E.end(); ++iter)
  {
     // extract rig id from image id
     const size_t  I = min(_map_RigIdPerImageId[iter->first.first], _map_RigIdPerImageId[iter->first.second]);
     const size_t  J = max(_map_RigIdPerImageId[iter->first.first], _map_RigIdPerImageId[iter->first.second]);

     const bool  bRigSwap = ( I != _map_RigIdPerImageId[iter->first.first] );

     if( I != J){

       if( !bRigSwap){
            _map_Matches_Rig[make_pair(I,J)][iter->first] = iter->second ;
       }
       else
       {
         // swap I and J matches in order to have matches grouped per rig id
         std::vector <matching::IndMatch>   matches;

         matches.clear();

         for(size_t  k(0); k < _map_Matches_E[iter->first].size(); ++k)
         {
            IndMatch  featurePair;

            featurePair._i = _map_Matches_E[iter->first][k]._j;
            featurePair._j = _map_Matches_E[iter->first][k]._i;

            matches.push_back( featurePair );
         }

         _map_Matches_Rig[make_pair(I,J)][make_pair(iter->first.second, iter->first.first)] = matches;

       }
    }
  }

}

void GlobalRigidReconstructionEngine::ComputeRelativeRt(
  Map_RelativeRT & vec_relatives)
{
  // For each pair, compute the rotation from pairwise point matches:
  // create rig structure using openGV
  translations_t  rigOffsets;
  rotations_t     rigRotations;
  double          averageFocal=0.0;

  for(int k=0; k < _vec_intrinsicGroups.size(); ++k)
  {
      translation_t   t = _vec_intrinsicGroups[k].m_rigC;
      rotation_t      R = _vec_intrinsicGroups[k].m_R.transpose();

      rigOffsets.push_back(t);
      rigRotations.push_back(R);
      averageFocal += _vec_intrinsicGroups[k].m_focal ;
  }

  averageFocal /= (double) _vec_intrinsicGroups.size();

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
  // loop on rigs
  for (int i = 0; i < _map_Matches_Rig.size(); ++i)
  {
    RigWiseMatches::const_iterator iter = _map_Matches_Rig.begin();
    std::advance(iter, i);

    // extract indices of matching rigs
    const size_t R0 = iter->first.first;
    const size_t R1 = iter->first.second;

    // initialize structure used for matching between rigs
    bearingVectors_t bearingVectorsRigOne;
    bearingVectors_t bearingVectorsRigTwo;

    std::vector<int>  camCorrespondencesRigOne;
    std::vector<int>  camCorrespondencesRigTwo;

    // set of subcam that are matched
    std::set<size_t>  setSubCam_rigOne;
    std::set<size_t>  setSubCam_rigTwo;

    // loop on inter-rig correspondences
    for( PairWiseMatches::const_iterator iterMatch =  iter->second.begin() ;
              iterMatch != iter->second.end() ; ++iterMatch ){

      // extract camera id and subchannel number
      const size_t I = iterMatch->first.first;
      const size_t J = iterMatch->first.second;

      const size_t SubI = _map_IntrinsicIdPerImageId[I];
      const size_t SubJ = _map_IntrinsicIdPerImageId[J];

      setSubCam_rigOne.insert(SubI);
      setSubCam_rigTwo.insert(SubJ);

      // extracts features for each pair in order to construct bearing vectors.
      for (size_t l = 0; l < iterMatch->second.size(); ++l)
      {
        bearingVector_t  bearing1;
        bearingVector_t  bearing2;

        // extract normalized keypoints coordinates
        bearing1(0) = _map_feats_normalized[I][iterMatch->second[l]._i].x();
        bearing1(1) = _map_feats_normalized[I][iterMatch->second[l]._i].y();
        bearing1(2) = 1.0;

        bearing2(0) = _map_feats_normalized[J][iterMatch->second[l]._j].x();
        bearing2(1) = _map_feats_normalized[J][iterMatch->second[l]._j].y();
        bearing2(2) = 1.0;

        // normalize bearing vectors
        bearing1 = bearing1 / bearing1.norm();
        bearing2 = bearing2 / bearing2.norm();

        // add bearing vectors to list and update correspondences list
        bearingVectorsRigOne.push_back( bearing1 );
        bearingVectorsRigTwo.push_back( bearing2 );

        camCorrespondencesRigOne.push_back(SubI);
        camCorrespondencesRigTwo.push_back(SubJ);
      }
    }// end loop on inter-rig matches

    //create non-central relative adapter
    relative_pose::NoncentralRelativeAdapter adapter(
          bearingVectorsRigOne,
          bearingVectorsRigTwo,
          camCorrespondencesRigOne,
          camCorrespondencesRigTwo,
          rigOffsets,
          rigRotations);

   //Create a NoncentralRelativePoseSacProblem and Ransac
    sac::Ransac<
        sac_problems::relative_pose::NoncentralRelativePoseSacProblem> ransac;
    boost::shared_ptr<
        sac_problems::relative_pose::NoncentralRelativePoseSacProblem>
        relposeproblem_ptr(
        new sac_problems::relative_pose::NoncentralRelativePoseSacProblem(
        adapter,
        sac_problems::relative_pose::NoncentralRelativePoseSacProblem::SIXPT));
    ransac.sac_model_ = relposeproblem_ptr;
    ransac.threshold_ = 2.0*(1.0 - cos(atan(sqrt(2.0) * 2.5 / averageFocal )));
    ransac.max_iterations_ = 4095;

   //Run the experiment
    struct timeval tic;
    struct timeval toc;
    gettimeofday( &tic, 0 );
    ransac.computeModel();
    gettimeofday( &toc, 0 );
    double ransac_time = TIMETODOUBLE(timeval_minus(toc,tic));

    double inlierProportion = ransac.inliers_.size() / ( (double) camCorrespondencesRigOne.size() );

    if( inlierProportion > 0.65
         &&  setSubCam_rigOne.size() > 0.24 * rigOffsets.size()
         &&  setSubCam_rigTwo.size() > 0.24 * rigOffsets.size() ){

      std::cout << inlierProportion << "  " << camCorrespondencesRigOne.size() << std::endl;

      // retrieve relative rig orientation and translation
      const Mat3  Rrig = ransac.model_coefficients_.block<3,3>(0,0).transpose();
      const Vec3  CRig = ransac.model_coefficients_.col(3);
      const Vec3  tRig = -Rrig * CRig;

      // compute point cloud associated and do BA to refine pose of rigs
      Mat3  K = Mat3::Identity();
      std::vector<Vec3> vec_allScenes;

      // count the number of observations
      size_t nbmeas = 0;

      // compute tracks between rigs
      RigWiseMatches map_matchesR0R1;
      map_matchesR0R1.insert(*_map_Matches_Rig.find(std::make_pair(R0,R1)));

      // Compute tracks:
      openMVG::tracks::STLMAPTracks map_tracks;
      TracksBuilder tracksBuilder;
      {
        tracksBuilder.Build(map_matchesR0R1);
        tracksBuilder.Filter(2);
        tracksBuilder.ExportToSTL(map_tracks);
      }

      // Triangulation of all the tracks
      {
        Map_Camera map_camera;
        std::vector<double> vec_residuals;
        vec_residuals.reserve(map_tracks.size());
        vec_allScenes.resize(map_tracks.size());

         for (int idx = 0; idx < map_tracks.size(); ++idx)
        {
          STLMAPTracks::const_iterator iterTracks = map_tracks.begin();
          std::advance(iterTracks, idx);

          const submapTrack & subTrack = iterTracks->second;

          // Look to the features required for the triangulation task
          Triangulation trianObj;
          for (submapTrack::const_iterator iterSubTrack = subTrack.begin(); iterSubTrack != subTrack.end(); ++iterSubTrack)
          {
            const size_t imaIndex = iterSubTrack->first;
            const size_t featIndex = iterSubTrack->second;
            const SIOPointFeature & pt = _map_feats_normalized[imaIndex][featIndex];

            // extract camera id and subchannel number
            const size_t SubI  = _map_IntrinsicIdPerImageId[imaIndex];
            const size_t rigId = _map_RigIdPerImageId[imaIndex];

            // compute pose of each cameras
            const Mat3 RcamI = _vec_intrinsicGroups[SubI].m_R;
            const Vec3 CI = _vec_intrinsicGroups[SubI].m_rigC;

            Vec3 tI; Mat3 RI;

            if( rigId == R0)
            {
               RI = RcamI;
               tI = -RcamI * CI;
            }
            else
            {
              RI = RcamI * Rrig;
              tI = RI*(-CRig - Rrig.transpose() * CI);
            }

            // build associated camera
            PinholeCamera cam(K, RI, tI);
            map_camera[imaIndex] = cam;

            // Build the P matrix
            trianObj.add(map_camera[imaIndex]._P, pt.coords().cast<double>());
          }

          // Compute the 3D point and keep point index with negative depth
          const Vec3 Xs = trianObj.compute();
          vec_allScenes[idx] = Xs;
        }
      }

      // now do bundle adjustment
      using namespace std;

      const size_t nbRigs = 2;
      const size_t nbCams = _vec_intrinsicGroups.size();
      const size_t nbPoints3D = vec_allScenes.size();

      // Count the number of measurement (sum of the reconstructed track length)
      size_t nbmeasurements = 0;
      for (STLMAPTracks::const_iterator iterTracks = map_tracks.begin();
        iterTracks != map_tracks.end(); ++iterTracks)
      {
        const submapTrack & subTrack = iterTracks->second;
        nbmeasurements += subTrack.size();
      }

      // Setup a BA problem
      using namespace openMVG::bundle_adjustment;
      BA_Problem_data_rigMotionAndIntrinsic<6,6,3> ba_problem; // Will refine [Rotations|Translations] and 3D points

      // Configure the size of the problem
      ba_problem.num_rigs_ = nbRigs;
      ba_problem.num_cameras_ = nbCams;
      ba_problem.num_intrinsics_ = nbCams;
      ba_problem.num_points_ = nbPoints3D;
      ba_problem.num_observations_ = nbmeasurements;

      ba_problem.rig_index_extrinsic.reserve(ba_problem.num_observations_);
      ba_problem.point_index_.reserve(ba_problem.num_observations_);
      ba_problem.camera_index_extrinsic.reserve(ba_problem.num_observations_);
      ba_problem.camera_index_intrinsic.reserve(ba_problem.num_observations_);
      ba_problem.observations_.reserve(2 * ba_problem.num_observations_);

      ba_problem.num_parameters_ =
        6 * ba_problem.num_rigs_         // rigs rotations / translations
        + 6 * ba_problem.num_cameras_    // #[Rotation|translation] = [3x1]|[3x1]
        + 3 * ba_problem.num_intrinsics_ // cameras intrinsics (focal and principal point)
        + 3 * ba_problem.num_points_;    // 3DPoints = [3x1]
      ba_problem.parameters_.reserve(ba_problem.num_parameters_);

      // Fill rigs
      {
        Mat3 R = Mat3::Identity();
        double angleAxis[3];
        ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);

        // translation
        Vec3 t = Vec3::Zero();

        ba_problem.parameters_.push_back(angleAxis[0]);
        ba_problem.parameters_.push_back(angleAxis[1]);
        ba_problem.parameters_.push_back(angleAxis[2]);
        ba_problem.parameters_.push_back(t[0]);
        ba_problem.parameters_.push_back(t[1]);
        ba_problem.parameters_.push_back(t[2]);
      }

      {
        double angleAxis[3];
        ceres::RotationMatrixToAngleAxis((const double*)Rrig.data(), angleAxis);

        // translation
        ba_problem.parameters_.push_back(angleAxis[0]);
        ba_problem.parameters_.push_back(angleAxis[1]);
        ba_problem.parameters_.push_back(angleAxis[2]);
        ba_problem.parameters_.push_back(tRig[0]);
        ba_problem.parameters_.push_back(tRig[1]);
        ba_problem.parameters_.push_back(tRig[2]);
      }

      // Setup rig camera position parameters
      for (size_t iterCam=0; iterCam < ba_problem.num_cameras_ ; ++iterCam )
      {

        const Mat3 R = _vec_intrinsicGroups[iterCam].m_R;
        double angleAxis[3];
        ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
        // translation
        const Vec3 t = -R * _vec_intrinsicGroups[iterCam].m_rigC;
        ba_problem.parameters_.push_back(angleAxis[0]);
        ba_problem.parameters_.push_back(angleAxis[1]);
        ba_problem.parameters_.push_back(angleAxis[2]);
        ba_problem.parameters_.push_back(t[0]);
        ba_problem.parameters_.push_back(t[1]);
        ba_problem.parameters_.push_back(t[2]);
      }

      // Setup rig camera intrinsics parameters
      for (size_t iterCam=0; iterCam < ba_problem.num_cameras_ ; ++iterCam )
      {
        ba_problem.parameters_.push_back( 1.0 );
        ba_problem.parameters_.push_back( 0.0 );
        ba_problem.parameters_.push_back( 0.0 );
      }

      // Fill 3D points
      for (std::vector<Vec3>::const_iterator iterPoint = vec_allScenes.begin();
        iterPoint != vec_allScenes.end();
        ++iterPoint)
      {
        const Vec3 & pt3D = *iterPoint;
        ba_problem.parameters_.push_back(pt3D[0]);
        ba_problem.parameters_.push_back(pt3D[1]);
        ba_problem.parameters_.push_back(pt3D[2]);
      }

      // Fill the measurements
      size_t k = 0;
      for (STLMAPTracks::const_iterator iterTracks = map_tracks.begin();
        iterTracks != map_tracks.end(); ++iterTracks, ++k)
      {
        // Look through the track and add point position
        const tracks::submapTrack & track = iterTracks->second;

        for( tracks::submapTrack::const_iterator iterTrack = track.begin();
          iterTrack != track.end();
          ++iterTrack)
        {
          const size_t imageId = iterTrack->first;
          const size_t featId = iterTrack->second;

          // If imageId reconstructed:
          //  - Add measurements (the feature position)
          //  - Add camidx (map the image number to the camera index)
          //  - Add ptidx (the 3D corresponding point index) (must be increasing)

          //if ( set_camIndex.find(imageId) != set_camIndex.end())
          {
            const std::vector<SIOPointFeature> & vec_feats = _map_feats_normalized[imageId];
            const SIOPointFeature & ptFeat = vec_feats[featId];

            ba_problem.observations_.push_back( ptFeat.x() );
            ba_problem.observations_.push_back( ptFeat.y() );

            ba_problem.point_index_.push_back(k);
            ba_problem.camera_index_extrinsic.push_back(_map_IntrinsicIdPerImageId[imageId]);
            ba_problem.camera_index_intrinsic.push_back(_map_IntrinsicIdPerImageId[imageId]);
            if ( _map_RigIdPerImageId[imageId] == R0 )
               ba_problem.rig_index_extrinsic.push_back(0);
            else
               ba_problem.rig_index_extrinsic.push_back(1);
          }
        }
      }

      // Create residuals for each observation in the bundle adjustment problem. The
      // parameters for cameras and points are added automatically.
      ceres::Problem problem;
      // Set a LossFunction to be less penalized by false measurements
      //  - set it to NULL if you don't want use a lossFunction.
      ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(2.0));
      for (size_t k = 0; k < ba_problem.num_observations(); ++k) {
        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.

        ceres::CostFunction* cost_function =
          new ceres::AutoDiffCostFunction<rig_pinhole_reprojectionError::ErrorFunc_Refine_Rig_Motion_3DPoints, 2, 3, 6, 6, 3>(
            new rig_pinhole_reprojectionError::ErrorFunc_Refine_Rig_Motion_3DPoints(
              &ba_problem.observations()[2 * k]));

        problem.AddResidualBlock(cost_function,
          p_LossFunction,
          ba_problem.mutable_camera_intrinsic_for_observation(k),
          ba_problem.mutable_camera_extrinsic_for_observation(k),
          ba_problem.mutable_rig_extrinsic_for_observation(k),
          ba_problem.mutable_point_for_observation(k));

        // fix intrinsic rig parameters
        problem.SetParameterBlockConstant(
          ba_problem.mutable_camera_extrinsic_for_observation(k) );
        problem.SetParameterBlockConstant(
          ba_problem.mutable_camera_intrinsic_for_observation(k) );
      }

      // fix rig one position
      problem.SetParameterBlockConstant(
        ba_problem.mutable_rig_extrinsic(0) );

      // Configure a BA engine and run it
      //  Make Ceres automatically detect the bundle structure.
      ceres::Solver::Options options;
      // Use a dense back-end since we only consider a two view problem
      options.linear_solver_type = ceres::DENSE_SCHUR;
      options.minimizer_progress_to_stdout = false;
      options.logging_type = ceres::SILENT;

      // Solve BA
      ceres::Solver::Summary summary;
      ceres::Solve(options, &problem, &summary);

      // If no error, get back refined parameters
      Mat3  R = Mat3::Identity();
      Vec3  t = Vec3::Zero();

      if (summary.IsSolutionUsable())
      {
          // Get back 3D points
          size_t k = 0;
          std::vector<Vec3>  finalPoint;

          for (std::vector<Vec3>::iterator iterPoint = vec_allScenes.begin();
                iterPoint != vec_allScenes.end(); ++iterPoint, ++k)
          {
              const double * pt = ba_problem.mutable_points() + k*3;
              Vec3 & pt3D = *iterPoint;
              pt3D = Vec3(pt[0], pt[1], pt[2]);
              finalPoint.push_back(pt3D);
          }

          // retrieve relative translation and rotation of rig
          // Get back rig 1
          Mat3  RotRigOne, RotRigTwo;
          Vec3  tRigOne, tRigTwo;

          {
            const double * cam = ba_problem.mutable_rig_extrinsic() + 0*6;

            // angle axis to rotation matrix
            ceres::AngleAxisToRotationMatrix(cam, RotRigOne.data());

            tRigOne[0] = cam[3]; tRigOne[1] = cam[4]; tRigOne[2] = cam[5];

          }
          // Get back rig 2
          {
            const double * cam = ba_problem.mutable_rig_extrinsic() + 1*6;

            // angle axis to rotation matrix
            ceres::AngleAxisToRotationMatrix(cam, RotRigTwo.data());

            tRigTwo[0] = cam[3]; tRigTwo[1] = cam[4]; tRigTwo[2]=cam[5];
          }

          RelativeCameraMotion(RotRigOne, tRigOne, RotRigTwo, tRigTwo, &R, &t);

      }
      // export rotation for rotation avereging
      vec_relatives[std::make_pair(R0,R1)] = std::make_pair(R,t);
    }
    else
    {
      std::cout << " Pose of rigs " << R0 << " and " << R1 << " is rejected. Matching camera "
      << setSubCam_rigOne.size() << "  " << setSubCam_rigTwo.size()
      << " Inliers proportion  " << inlierProportion << std::endl;
    }
  }
}

void GlobalRigidReconstructionEngine::tripletListing(std::vector< graphUtils::Triplet > & vec_triplets,
   const size_t & numberOfSubcam ) const
{
  vec_triplets.clear();

  imageGraph::indexedImageGraph putativeGraph(_map_Matches_Rig, _vec_fileNames, numberOfSubcam );

  List_Triplets<imageGraph::indexedImageGraph::GraphT>(putativeGraph.g, vec_triplets);

  //Change triplets to ImageIds
  for (size_t i = 0; i < vec_triplets.size(); ++i)
  {
    graphUtils::Triplet & triplet = vec_triplets[i];
    size_t I = triplet.i, J = triplet.j , K = triplet.k;
    I = (*putativeGraph.map_nodeMapIndex)[putativeGraph.g.nodeFromId(I)];
    J = (*putativeGraph.map_nodeMapIndex)[putativeGraph.g.nodeFromId(J)];
    K = (*putativeGraph.map_nodeMapIndex)[putativeGraph.g.nodeFromId(K)];
    size_t triplet_[3] = {I,J,K};
    std::sort(&triplet_[0], &triplet_[3]);
    triplet = graphUtils::Triplet(triplet_[0],triplet_[1],triplet_[2]);
  }
}

void GlobalRigidReconstructionEngine::tripletRotationRejection(
  std::vector< graphUtils::Triplet > & vec_triplets,
  Map_RelativeRT & map_relatives)
{
  Map_RelativeRT map_relatives_validated;

  // DETECTION OF ROTATION OUTLIERS
  std::vector< graphUtils::Triplet > vec_triplets_validated;

  std::vector<float> vec_errToIdentityPerTriplet;
  vec_errToIdentityPerTriplet.reserve(vec_triplets.size());
  // Compute for each length 3 cycles: the composition error
  //  Error to identity rotation.
  for (size_t i = 0; i < vec_triplets.size(); ++i)
  {
    const graphUtils::Triplet & triplet = vec_triplets[i];
    size_t I = triplet.i, J = triplet.j , K = triplet.k;

    //-- Find the three rotations
    const std::pair<size_t,size_t> ij = std::make_pair(I,J);
    const std::pair<size_t,size_t> ji = std::make_pair(J,I);

    Mat3 RIJ;
    if (map_relatives.find(ij) != map_relatives.end())
      RIJ = map_relatives.find(ij)->second.first;
    else
      RIJ = map_relatives.find(ji)->second.first.transpose();

    const std::pair<size_t,size_t> jk = std::make_pair(J,K);
    const std::pair<size_t,size_t> kj = std::make_pair(K,J);

    Mat3 RJK;
    if (map_relatives.find(jk) != map_relatives.end())
      RJK = map_relatives.find(jk)->second.first;
    else
      RJK = map_relatives.find(kj)->second.first.transpose();

    const std::pair<size_t,size_t> ki = std::make_pair(K,I);
    const std::pair<size_t,size_t> ik = std::make_pair(I,K);

    Mat3 RKI;
    if (map_relatives.find(ki) != map_relatives.end())
      RKI = map_relatives.find(ki)->second.first;
    else
      RKI = map_relatives.find(ik)->second.first.transpose();

    Mat3 Rot_To_Identity = RIJ * RJK * RKI; // motion composition
    float angularErrorDegree = static_cast<float>(R2D(getRotationMagnitude(Rot_To_Identity)));
    vec_errToIdentityPerTriplet.push_back(angularErrorDegree);

    if (angularErrorDegree < 2.0f)
    {
      vec_triplets_validated.push_back(triplet);

      if (map_relatives.find(ij) != map_relatives.end())
        map_relatives_validated[ij] = map_relatives.find(ij)->second;
      else
        map_relatives_validated[ji] = map_relatives.find(ji)->second;

      if (map_relatives.find(jk) != map_relatives.end())
        map_relatives_validated[jk] = map_relatives.find(jk)->second;
      else
        map_relatives_validated[kj] = map_relatives.find(kj)->second;

      if (map_relatives.find(ki) != map_relatives.end())
        map_relatives_validated[ki] = map_relatives.find(ki)->second;
      else
        map_relatives_validated[ik] = map_relatives.find(ik)->second;
    }
  }

  map_relatives.swap(map_relatives_validated);

  // Display statistics about rotation triplets error:
  std::cout << "\nStatistics about rotation triplets:" << std::endl;
  minMaxMeanMedian<float>(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end());

  std::sort(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end());

  Histogram<float> histo(0.0f, *max_element(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end()), 180);
  histo.Add(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end());

  svgHisto histosvg;
  histosvg.draw(histo.GetHist(),
                std::make_pair(0.0f, *max_element(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end())),
                stlplus::create_filespec(this->_sOutDirectory, "Triplet_Rotation_Residual_180.svg"),
                600,300);

  histo = Histogram<float>(0.0f, *max_element(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end()), 20);
  histo.Add(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end());

  histosvg.draw(histo.GetHist(),
                std::make_pair(0.0f, *max_element(vec_errToIdentityPerTriplet.begin(), vec_errToIdentityPerTriplet.end())),
                stlplus::create_filespec(this->_sOutDirectory, "Triplet_Rotation_Residual_20.svg"),
                600,300);

  typedef lemon::ListGraph Graph;
  imageGraph::indexedImageGraph putativeGraph(_map_Matches_Rig, _vec_fileNames, _map_ImagesIdPerRigId[0].size() );

  Graph::EdgeMap<bool> edge_filter(putativeGraph.g, false);
  Graph::NodeMap<bool> node_filter(putativeGraph.g, true);

  typedef SubGraph<Graph > subGraphT;
  subGraphT sg(putativeGraph.g, node_filter, edge_filter);

  // Look all edges of the graph and look if exist in one triplet
  for (Graph::EdgeIt iter(putativeGraph.g); iter!=INVALID; ++iter)
  {
    const size_t Idu = (*putativeGraph.map_nodeMapIndex)[sg.u(iter)];
    const size_t Idv = (*putativeGraph.map_nodeMapIndex)[sg.v(iter)];
    //-- Look if the edge Idu,Idv exists in the trifocal tensor list
    for (size_t i = 0; i < vec_triplets_validated.size(); ++i)
    {
      const graphUtils::Triplet & triplet = vec_triplets_validated[i];
      if ( triplet.contain(std::make_pair(Idu, Idv)))
      {
        edge_filter[iter] = true;
        break;
      }
    }
  }

  imageGraph::exportToGraphvizData(
    stlplus::create_filespec(_sOutDirectory, "cleanedGraphTripletRotation"),
    sg);

  {
    std::cout << "\nTriplets filtering based on error on cycles \n";
    std::cout << "Before: " << vec_triplets.size() << " triplets \n"
    << "After: " << vec_triplets_validated.size() << std::endl;
    std::cout << "There is " << lemon::countConnectedComponents (sg)
      << " Connected Component in the filtered graph" << std::endl;
  }

  vec_triplets.clear();
  vec_triplets = vec_triplets_validated;

  size_t removedEdgesCount = 0;

  //-- Remove false edges from the rejected triplets
  {
    for (Graph::EdgeIt iter(putativeGraph.g); iter!=INVALID; ++iter)
    {
      if (!edge_filter[iter])
      {
        removedEdgesCount++;

        size_t Idu = (*putativeGraph.map_nodeMapIndex)[sg.u(iter)];
        size_t Idv = (*putativeGraph.map_nodeMapIndex)[sg.v(iter)];

        //-- Clean relatives matches
        RigWiseMatches::iterator iterF = _map_Matches_Rig.find(std::make_pair(Idu,Idv));
        if (iterF != _map_Matches_Rig.end())
        {
          _map_Matches_Rig.erase(iterF);
        }
        else
        {
          iterF = _map_Matches_Rig.find(std::make_pair(Idv,Idu));
          if (iterF != _map_Matches_Rig.end())
            _map_Matches_Rig.erase(iterF);
        }

        //-- Clean relative motions
        Map_RelativeRT::iterator iterF2 = map_relatives.find(std::make_pair(Idu,Idv));
        if (iterF2 != map_relatives.end())
        {
          map_relatives.erase(iterF2);
        }
        else
        {
          iterF2 = map_relatives.find(std::make_pair(Idv,Idu));
          if (iterF2!= map_relatives.end())
            map_relatives.erase(iterF2);
        }
      }
    }
  }

  std::cout << "\n Relatives edges removed by triplet checking: " << removedEdgesCount << std::endl;
}

void GlobalRigidReconstructionEngine::bundleAdjustment(
    Map_Rig & map_rig,
    Map_Camera & map_camera,
    std::vector<Vec3> & vec_allScenes,
    const STLMAPTracks & map_tracksSelected,
    bool bRefineRotation,
    bool bRefineTranslation,
    bool bRefineIntrinsics)
{
  using namespace std;

  const size_t nbRigs = map_rig.size();
  const size_t nbCams = _vec_intrinsicGroups.size();
  const size_t nbPoints3D = vec_allScenes.size();

  // Count the number of measurement (sum of the reconstructed track length)
  size_t nbmeasurements = 0;
  for (STLMAPTracks::const_iterator iterTracks = map_tracksSelected.begin();
    iterTracks != map_tracksSelected.end(); ++iterTracks)
  {
    const submapTrack & subTrack = iterTracks->second;
    nbmeasurements += subTrack.size();
  }

  // Setup a BA problem
  using namespace openMVG::bundle_adjustment;
  BA_Problem_data_rigMotionAndIntrinsic<6,6,3> ba_problem; // Will refine [Rotations|Translations] and 3D points

  // Configure the size of the problem
  ba_problem.num_rigs_ = nbRigs;
  ba_problem.num_cameras_ = nbCams;
  ba_problem.num_intrinsics_ = nbCams;
  ba_problem.num_points_ = nbPoints3D;
  ba_problem.num_observations_ = nbmeasurements;

  ba_problem.rig_index_extrinsic.reserve(ba_problem.num_observations_);
  ba_problem.point_index_.reserve(ba_problem.num_observations_);
  ba_problem.camera_index_extrinsic.reserve(ba_problem.num_observations_);
  ba_problem.camera_index_intrinsic.reserve(ba_problem.num_observations_);
  ba_problem.observations_.reserve(2 * ba_problem.num_observations_);

  ba_problem.num_parameters_ =
    6 * ba_problem.num_rigs_         // rigs rotations / translations
    + 6 * ba_problem.num_cameras_    // #[Rotation|translation] = [3x1]|[3x1]
    + 3 * ba_problem.num_intrinsics_ // cameras intrinsics (focal and principal point)
    + 3 * ba_problem.num_points_;    // 3DPoints = [3x1]
  ba_problem.parameters_.reserve(ba_problem.num_parameters_);

  // Fill rigs
  std::set<size_t> set_rigIndex;
  std::map<size_t,size_t> map_rigIndexToNumber_extrinsic;
  size_t cpt = 0;
  for (Map_Rig::const_iterator iter = map_rig.begin();
    iter != map_rig.end();  ++iter, ++cpt)
  {
    // in order to map camera index to contiguous number
    set_rigIndex.insert(iter->first);
    map_rigIndexToNumber_extrinsic.insert(std::make_pair(iter->first, cpt));

    const Mat3 R = iter->second.first;
    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // translation
    const Vec3 t = iter->second.second;
    ba_problem.parameters_.push_back(angleAxis[0]);
    ba_problem.parameters_.push_back(angleAxis[1]);
    ba_problem.parameters_.push_back(angleAxis[2]);
    ba_problem.parameters_.push_back(t[0]);
    ba_problem.parameters_.push_back(t[1]);
    ba_problem.parameters_.push_back(t[2]);
  }

  // Setup rig camera position parameters
  for (size_t iter=0; iter < ba_problem.num_cameras_ ; ++iter )
  {
    const Mat3 R = _vec_intrinsicGroups[iter].m_R;
    double angleAxis[3];
    ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
    // translation
    const Vec3 t = -R * _vec_intrinsicGroups[iter].m_rigC;
    ba_problem.parameters_.push_back(angleAxis[0]);
    ba_problem.parameters_.push_back(angleAxis[1]);
    ba_problem.parameters_.push_back(angleAxis[2]);
    ba_problem.parameters_.push_back(t[0]);
    ba_problem.parameters_.push_back(t[1]);
    ba_problem.parameters_.push_back(t[2]);
  }

  // Setup rig camera intrinsics parameters
  for (size_t iter=0; iter < ba_problem.num_cameras_ ; ++iter )
  {
    const Mat3 K = _vec_intrinsicGroups[iter].m_K;
    ba_problem.parameters_.push_back( K(0,0) );
    ba_problem.parameters_.push_back( K(0,2) );
    ba_problem.parameters_.push_back( K(1,2) );
  }

  // Fill 3D points
  for (std::vector<Vec3>::const_iterator iter = vec_allScenes.begin();
    iter != vec_allScenes.end();
    ++iter)
  {
    const Vec3 & pt3D = *iter;
    ba_problem.parameters_.push_back(pt3D[0]);
    ba_problem.parameters_.push_back(pt3D[1]);
    ba_problem.parameters_.push_back(pt3D[2]);
  }

  // Fill the measurements
  size_t k = 0;
  for (STLMAPTracks::const_iterator iterTracks = map_tracksSelected.begin();
    iterTracks != map_tracksSelected.end(); ++iterTracks, ++k)
  {
    // Look through the track and add point position
    const tracks::submapTrack & track = iterTracks->second;

    for( tracks::submapTrack::const_iterator iterTrack = track.begin();
      iterTrack != track.end();
      ++iterTrack)
    {
      const size_t imageId = iterTrack->first;
      const size_t featId = iterTrack->second;
      const size_t rigId = _map_RigIdPerImageId[imageId];

      // If imageId reconstructed:
      //  - Add measurements (the feature position)
      //  - Add camidx (map the image number to the camera index)
      //  - Add ptidx (the 3D corresponding point index) (must be increasing)

      //if ( set_camIndex.find(imageId) != set_camIndex.end())
      {
        const std::vector<SIOPointFeature> & vec_feats = _map_feats[imageId];
        const SIOPointFeature & ptFeat = vec_feats[featId];

        ba_problem.observations_.push_back( ptFeat.x() );
        ba_problem.observations_.push_back( ptFeat.y() );

        ba_problem.point_index_.push_back(k);
        ba_problem.camera_index_extrinsic.push_back(_map_IntrinsicIdPerImageId[imageId]);
        ba_problem.camera_index_intrinsic.push_back(_map_IntrinsicIdPerImageId[imageId]);
        ba_problem.rig_index_extrinsic.push_back(map_rigIndexToNumber_extrinsic[rigId]);

      }
    }
  }

  // Create residuals for each observation in the bundle adjustment problem. The
  // parameters for cameras and points are added automatically.
  ceres::Problem problem;
  // Set a LossFunction to be less penalized by false measurements
  //  - set it to NULL if you don't want use a lossFunction.
  ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(2.0));
  for (size_t k = 0; k < ba_problem.num_observations(); ++k) {
    // Each Residual block takes a point and a camera as input and outputs a 2
    // dimensional residual. Internally, the cost function stores the observed
    // image location and compares the reprojection against the observation.

    ceres::CostFunction* cost_function =
      new ceres::AutoDiffCostFunction<rig_pinhole_reprojectionError::ErrorFunc_Refine_Rig_Motion_3DPoints, 2, 3, 6, 6, 3>(
        new rig_pinhole_reprojectionError::ErrorFunc_Refine_Rig_Motion_3DPoints(
          &ba_problem.observations()[2 * k]));

    problem.AddResidualBlock(cost_function,
      p_LossFunction,
      ba_problem.mutable_camera_intrinsic_for_observation(k),
      ba_problem.mutable_camera_extrinsic_for_observation(k),
      ba_problem.mutable_rig_extrinsic_for_observation(k),
      ba_problem.mutable_point_for_observation(k));
  }

  // Configure a BA engine and run it
  //  Make Ceres automatically detect the bundle structure.
  ceres::Solver::Options options;
  options.preconditioner_type = ceres::JACOBI;
  options.linear_solver_type = ceres::SPARSE_SCHUR;
  if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
    options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
  else
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
      options.sparse_linear_algebra_library_type = ceres::CX_SPARSE;
    else
    {
      // No sparse backend for Ceres.
      // Use dense solving
      options.linear_solver_type = ceres::DENSE_SCHUR;
    }
  options.minimizer_progress_to_stdout = false;
  options.logging_type = ceres::SILENT;
  #ifdef USE_OPENMP
  options.num_threads = omp_get_max_threads();
  options.num_linear_solver_threads = omp_get_max_threads();
  #endif // USE_OPENMP

  // Configure constant parameters (if any)
  {
    std::vector<int> vec_constant_extrinsic; // [R|t]
    if (!bRefineRotation)
    {
      vec_constant_extrinsic.push_back(0);
      vec_constant_extrinsic.push_back(1);
      vec_constant_extrinsic.push_back(2);
    }
    if (!bRefineTranslation)
    {
      vec_constant_extrinsic.push_back(3);
      vec_constant_extrinsic.push_back(4);
      vec_constant_extrinsic.push_back(5);
    }

    for (size_t iExtrinsicId = 0; iExtrinsicId < ba_problem.num_rigs_; ++iExtrinsicId)
    {
      if (!vec_constant_extrinsic.empty())
      {
        ceres::SubsetParameterization *subset_parameterization =
          new ceres::SubsetParameterization(6, vec_constant_extrinsic);
        problem.SetParameterization(ba_problem.mutable_rig_extrinsic(iExtrinsicId),
          subset_parameterization);
      }
    }

    if (!bRefineIntrinsics)
    {
      // No camera intrinsics are being refined,
      // set the whole parameter block as constant for best performance.
      for (size_t iIntrinsicGroupId = 0; iIntrinsicGroupId < ba_problem.num_intrinsics() ; ++iIntrinsicGroupId)
      {
        problem.SetParameterBlockConstant(ba_problem.mutable_cameras_intrinsic(iIntrinsicGroupId));
        problem.SetParameterBlockConstant(ba_problem.mutable_cameras_extrinsic(iIntrinsicGroupId));
      }
    }
  }

  // Solve BA
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  std::cout << summary.FullReport() << std::endl;

  // If no error, get back refined parameters
  if (summary.IsSolutionUsable())
  {
    // Display statistics about the minimization
    std::cout << std::endl
      << "Bundle Adjustment statistics:\n"
      << " Initial RMSE: " << std::sqrt( summary.initial_cost / (ba_problem.num_observations_*2.)) << "\n"
      << " Final RMSE: " << std::sqrt( summary.final_cost / (ba_problem.num_observations_*2.)) << "\n"
      << std::endl;

    // Get back rigs
    size_t i = 0;
    for (Map_Rig::iterator iter = map_rig.begin();
      iter != map_rig.end(); ++iter, ++i)
    {
      // camera motion [R|t]
      const double * cam = ba_problem.mutable_rig_extrinsic(i);
      Mat3 R;
      // angle axis to rotation matrix
      ceres::AngleAxisToRotationMatrix(cam, R.data());
      Vec3 t(cam[3], cam[4], cam[5]);

      std::pair<Mat3, Vec3>  & pRigPose = iter->second ;
      pRigPose = std::make_pair(R, t);
    }

    // Get back 3D points
    i = 0;
    for (std::vector<Vec3>::iterator iter = vec_allScenes.begin();
      iter != vec_allScenes.end(); ++iter, ++i)
    {
      const double * pt = ba_problem.mutable_points() + i*3;
      Vec3 & pt3D = *iter;
      pt3D = Vec3(pt[0], pt[1], pt[2]);
    }

    // Get back camera intrinsics
    i = 0;
    for (i = 0 ; i < _vec_intrinsicGroups.size() ; ++i)
    {
      // camera motion [R|t]
      const double * cam = ba_problem.mutable_cameras_extrinsic(i);
      Mat3 R;
      // angle axis to rotation matrix
      ceres::AngleAxisToRotationMatrix(cam, R.data());
      Vec3 t(cam[3], cam[4], cam[5]);

      // camera intrinsics [f,ppx,ppy]
      double * intrinsics = ba_problem.mutable_cameras_intrinsic(i);
      Mat3 K = Mat3::Identity();
      K << intrinsics[0], 0, intrinsics[1],
           0, intrinsics[0], intrinsics[2],
           0, 0, 1;
      // Update the camera
      _vec_intrinsicGroups[i].m_R    = R;
      _vec_intrinsicGroups[i].m_rigC = -R.transpose() * t ;
      _vec_intrinsicGroups[i].m_K    = K;
    }

    // get back camera positions
    i = 0;
    for (Map_Camera::iterator iter = map_camera.begin();
      iter != map_camera.end(); ++iter, ++i)
    {
      // camera motion [R|t]
      const size_t I = iter->first;
      const size_t intrinsicId = _map_IntrinsicIdPerImageId[I];
      const size_t rigId       = _map_RigIdPerImageId[I];

      // compute camera rotation
      const Mat3 & Ri = map_rig.find(rigId)->second.first ;
      const Mat3 & Rcam = _vec_intrinsicGroups[intrinsicId].m_R;
      const Mat3 & R = Rcam * Ri;

      // compute camera translation
      Vec3 Rigt = map_rig.find(rigId)->second.second ;
      const Vec3 tCam = -Rcam * _vec_intrinsicGroups[intrinsicId].m_rigC ;
      const Vec3 t = Rcam * Rigt + tCam;

      const Mat3 & _K = _vec_intrinsicGroups[intrinsicId].m_K;   // The same K matrix is used by all the camera
      map_camera[I] = PinholeCamera(_K, R, t);

    }

    {
      //-- Export Second bundle adjustment statistics
      if (_bHtmlReport)
      {
        using namespace htmlDocument;
        std::ostringstream os;
        os << "Bundle Adjustment statistics.";
        _htmlDocStream->pushInfo("<hr>");
        _htmlDocStream->pushInfo(htmlMarkup("h1",os.str()));

        os.str("");
        os << "-------------------------------" << "<br>"
          << "-- #observation: " << ba_problem.num_observations_ << ".<br>"
          << "-- residual mean (RMSE): " << std::sqrt(summary.final_cost/ba_problem.num_observations_) << ".<br>"
          << "-- Nb Steps required until convergence: " <<  summary.num_successful_steps + summary.num_unsuccessful_steps << ".<br>"
          << "-------------------------------" << "<br>";
        _htmlDocStream->pushInfo(os.str());
      }

      std::cout << "\n"
        << "-------------------------------" << "\n"
        << "-- #observation: " << ba_problem.num_observations_ << ".\n"
        << "-- residual mean (RMSE): " << std::sqrt(summary.final_cost/ba_problem.num_observations_) << ".\n"
        << "-- Nb Steps required until convergence: " <<  summary.num_successful_steps + summary.num_unsuccessful_steps << ".\n"
        << "-------------------------------" << std::endl;
    }
  }
}

void GlobalRigidReconstructionEngine::ColorizeTracks(
  const STLMAPTracks & map_tracks,
  std::vector<Vec3> & vec_tracksColor) const
{
  // Colorize each track
  //  Start with the most representative image
  //    and iterate to provide a color to each 3D point

  {
    C_Progress_display my_progress_bar(map_tracks.size(),
                                       std::cout,
                                       "\nCompute scene structure color\n");

    vec_tracksColor.resize(map_tracks.size());

    //Build a list of contiguous index for the trackIds
    std::map<size_t, size_t> trackIds_to_contiguousIndexes;
    size_t cpt = 0;
    for (openMVG::tracks::STLMAPTracks::const_iterator it = map_tracks.begin();
      it != map_tracks.end(); ++it, ++cpt)
    {
      trackIds_to_contiguousIndexes[it->first] = cpt;
    }

    // The track list that will be colored (point removed during the process)
    openMVG::tracks::STLMAPTracks mapTrackToColor(map_tracks);
    while( !mapTrackToColor.empty() )
    {
      // Find the most representative image
      //  a. Count the number of visible point for each image
      //  b. Sort to find the most representative image

      std::map<size_t, size_t> map_IndexCardinal; // ImageIndex, Cardinal
      for (openMVG::tracks::STLMAPTracks::const_iterator
        iterT = mapTrackToColor.begin();
        iterT != mapTrackToColor.end();
      ++iterT)
      {
        const size_t trackId = iterT->first;
        const tracks::submapTrack & track = mapTrackToColor[trackId];
        for( tracks::submapTrack::const_iterator iterTrack = track.begin();
          iterTrack != track.end(); ++iterTrack)
        {
          const size_t imageId = iterTrack->first;
          if (map_IndexCardinal.find(imageId) == map_IndexCardinal.end())
            map_IndexCardinal[imageId] = 1;
          else
            ++map_IndexCardinal[imageId];
        }
      }

      // Find the image that is the most represented
      std::vector<size_t> vec_cardinal;
      std::transform(map_IndexCardinal.begin(),
        map_IndexCardinal.end(),
        std::back_inserter(vec_cardinal),
        RetrieveValue());
      using namespace indexed_sort;
      std::vector< sort_index_packet_descend< size_t, size_t> > packet_vec(vec_cardinal.size());
      sort_index_helper(packet_vec, &vec_cardinal[0]);

      //First index is the image with the most of matches
      std::map<size_t, size_t>::const_iterator iterTT = map_IndexCardinal.begin();
      std::advance(iterTT, packet_vec[0].index);
      const size_t indexImage = iterTT->first;
      Image<RGBColor> image;
      ReadImage(
        stlplus::create_filespec(
        _sImagePath,
        stlplus::basename_part(_vec_camImageNames[indexImage].m_sImageName),
        stlplus::extension_part(_vec_camImageNames[indexImage].m_sImageName) ).c_str(), &image);

      // Iterate through the track
      std::set<size_t> set_toRemove;
      for (openMVG::tracks::STLMAPTracks::const_iterator
        iterT = mapTrackToColor.begin();
        iterT != mapTrackToColor.end();
        ++iterT)
      {
        const size_t trackId = iterT->first;
        const tracks::submapTrack & track = mapTrackToColor[trackId];
        tracks::submapTrack::const_iterator iterF = track.find(indexImage);

        if (iterF != track.end())
        {
          // Color the track
          const size_t featId = iterF->second;
          const SIOPointFeature & feat = _map_feats.find(indexImage)->second[featId];
          const RGBColor color = image(feat.y(), feat.x());

          vec_tracksColor[ trackIds_to_contiguousIndexes[trackId] ] = Vec3(color.r(), color.g(), color.b());
          set_toRemove.insert(trackId);
          ++my_progress_bar;
        }
      }
      // Remove colored track
      for (std::set<size_t>::const_iterator iter = set_toRemove.begin();
        iter != set_toRemove.end(); ++iter)
      {
        mapTrackToColor.erase(*iter);
      }
    }
  }
}

} // namespace openMVG
