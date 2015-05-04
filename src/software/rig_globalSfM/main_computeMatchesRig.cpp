
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/image/image.hpp"
#include "openMVG/features/features.hpp"

/// Generic Image Collection image matching
#include "openMVG/matching_image_collection/Matcher_AllInMemory.hpp"
#include "openMVG/matching_image_collection/rig_GeometricFilter.hpp"
#include "openMVG/matching_image_collection/Rig_ACRobust.hpp"
#include "software/SfM/pairwiseAdjacencyDisplay.hpp"
#include "software/SfM/SfMIOHelper.hpp"
#include "openMVG/matching/matcher_brute_force.hpp"
#include "openMVG/matching/matcher_kdtree_flann.hpp"
#include "openMVG/matching/indMatch_utils.hpp"

/// Feature detector and descriptor interface
#include "nonFree/sift/SIFT.hpp"

#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>

using namespace openMVG;
using namespace openMVG::matching;
using namespace openMVG::robust;
using namespace std;

enum ePairMode
{
  PAIR_EXHAUSTIVE = 0,
  PAIR_CONTIGUOUS = 1,
  PAIR_FROM_FILE  = 2
};

// Equality functor to count the number of similar K matrices in the essential matrix case.
bool testIntrinsicsEquality(
  SfMIO::IntrinsicCameraRigInfo const &ci1,
  SfMIO::IntrinsicCameraRigInfo const &ci2)
{
  return ci1.m_K == ci2.m_K && ci1.m_sCameraMaker == ci2.m_sCameraMaker && ci1.m_sCameraModel == ci2.m_sCameraModel;
}

int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sImaDirectory;
  std::string sOutDir = "";
  std::string sGeometricModel = "f";
  float fDistRatio = .6f;
  bool bOctMinus1 = false;
  float dPeakThreshold = 0.04f;
  int iMatchingVideoMode = -1;
  std::string sPredefinedPairList = "";

  cmd.add( make_option('i', sImaDirectory, "imadir") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('r', fDistRatio, "distratio") );
  cmd.add( make_option('s', bOctMinus1, "octminus1") );
  cmd.add( make_option('p', dPeakThreshold, "peakThreshold") );
  cmd.add( make_option('g', sGeometricModel, "geometricModel") );
  cmd.add( make_option('v', iMatchingVideoMode, "videoModeMatching") );
  cmd.add( make_option('l', sPredefinedPairList, "pairList") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--imadir path] \n"
      << "[-o|--outdir path] \n"
      << "\n[Optional]\n"
      << "[-r|--distratio 0.6] \n"
      << "[-s|--octminus1 0 or 1] \n"
      << "[-p|--peakThreshold 0.04 -> 0.01] \n"
      << "[-g|--geometricModel f, e or h] \n"
      << "[-v|--videoModeMatching 2 -> X] \n"
      << "\t sequence matching with an overlap of X images\n"
      << "[-l]--pairList file"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--imadir " << sImaDirectory << std::endl
            << "--outdir " << sOutDir << std::endl
            << "--distratio " << fDistRatio << std::endl
            << "--octminus1 " << bOctMinus1 << std::endl
            << "--peakThreshold " << dPeakThreshold << std::endl
            << "--geometricModel " << sGeometricModel << std::endl
            << "--videoModeMatching " << iMatchingVideoMode << std::endl;

  ePairMode ePairmode = (iMatchingVideoMode == -1 ) ? PAIR_EXHAUSTIVE : PAIR_CONTIGUOUS;

  if (sPredefinedPairList.length()) {
    std::cout << "--pairList " << sPredefinedPairList << std::endl;
    ePairmode = PAIR_FROM_FILE;
    if (iMatchingVideoMode>0) {
      std::cerr << "\nIncompatible options: --videoModeMatching and --pairList" << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  std::string  sGeometricMatchesFilename = "matches.e.txt";

  // -----------------------------
  // a. List images
  // b. Compute features and descriptors
  // c. Compute putative descriptor matches
  // d. Geometric filtering of putative matches
  // e. Export some statistics
  // -----------------------------

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create( sOutDir );

  //---------------------------------------
  // a. List images
  //---------------------------------------
  std::string sListsFile = stlplus::create_filespec(sOutDir, "lists.txt" );
  if (!stlplus::is_file(sListsFile)) {
    std::cerr << std::endl
      << "The input file \""<< sListsFile << "\" is missing" << std::endl;
    return false;
  }

  std::vector<openMVG::SfMIO::CameraRigInfo> vec_camImageName;
  std::vector<openMVG::SfMIO::IntrinsicCameraRigInfo> vec_focalGroup;
  if (!openMVG::SfMIO::loadImageList( vec_camImageName,
                                      vec_focalGroup,
                                      sListsFile) )
  {
    std::cerr << "\nEmpty or invalid image list." << std::endl;
    return false;
  }

  //-- Two alias to ease access to image filenames and image sizes
  std::vector<std::string> vec_fileNames;
  std::vector<std::pair<size_t, size_t> > vec_imagesSize;
  for ( std::vector<openMVG::SfMIO::CameraRigInfo>::const_iterator
    iter_camInfo = vec_camImageName.begin();
    iter_camInfo != vec_camImageName.end();
    iter_camInfo++ )
  {
    vec_imagesSize.push_back( std::make_pair( vec_focalGroup[iter_camInfo->m_intrinsicId].m_w,
                                              vec_focalGroup[iter_camInfo->m_intrinsicId].m_h ) );
    vec_fileNames.push_back( stlplus::create_filespec( sImaDirectory, iter_camInfo->m_sImageName) );
  }

  //---------------------------------------
  // b. Compute features and descriptor
  //    - extract sift features and descriptor
  //    - if keypoints already computed, re-load them
  //    - else save features and descriptors on disk
  //---------------------------------------

  typedef Descriptor<unsigned char, 128> DescriptorT;
  typedef SIOPointFeature FeatureT;
  typedef std::vector<FeatureT> FeatsT;
  typedef vector<DescriptorT > DescsT;
  typedef KeypointSet<FeatsT, DescsT > KeypointSetT;

  {
    Timer timer;
    std::cout << "\n\n - EXTRACT FEATURES - " << std::endl;
    vec_imagesSize.resize(vec_fileNames.size());

    Image<unsigned char> imageGray;

    C_Progress_display my_progress_bar( vec_fileNames.size() );
    for(size_t i=0; i < vec_fileNames.size(); ++i)  {

      std::string sFeat = stlplus::create_filespec(sOutDir,
        stlplus::basename_part(vec_fileNames[i]), "feat");
      std::string sDesc = stlplus::create_filespec(sOutDir,
        stlplus::basename_part(vec_fileNames[i]), "desc");

      //If descriptors or features file are missing, compute them
      if (!stlplus::file_exists(sFeat) || !stlplus::file_exists(sDesc)) {

        if (!ReadImage(vec_fileNames[i].c_str(), &imageGray))
          continue;

        // Compute features and descriptors and export them to files
        KeypointSetT kpSet;
        SIFTDetector(imageGray,
          kpSet.features(), kpSet.descriptors(),
          bOctMinus1, true, dPeakThreshold);
        kpSet.saveToBinFile(sFeat, sDesc);
      }
      ++my_progress_bar;
    }
    std::cout << "Task done in (s): " << timer.elapsed() << std::endl;
  }

  //---------------------------------------
  // c. Compute putative descriptor matches
  //    - L2 descriptor matching
  //    - Keep correspondences only if NearestNeighbor ratio is ok
  //---------------------------------------
  PairWiseMatches map_PutativesMatches;
  // Define the matcher and the used metric (Squared L2)
  // ANN matcher could be defined as follow:
  typedef flann::L2<DescriptorT::bin_type> MetricT;
  typedef ArrayMatcher_Kdtree_Flann<DescriptorT::bin_type, MetricT> MatcherT;
  // Brute force matcher can be defined as following:
  //typedef L2_Vectorized<DescriptorT::bin_type> MetricT;
  //typedef ArrayMatcherBruteForce<DescriptorT::bin_type, MetricT> MatcherT;

  std::cout << std::endl << " - PUTATIVE MATCHES - " << std::endl;
  // If the matches already exists, reload them
  if (stlplus::file_exists(sOutDir + "/matches.putative.txt"))
  {
    PairedIndMatchImport(sOutDir + "/matches.putative.txt", map_PutativesMatches);
    std::cout << "\t PREVIOUS RESULTS LOADED" << std::endl;
  }
  else // Compute the putative matches
  {
    std::cout << "Use: ";
    switch (ePairmode)
    {
      case PAIR_EXHAUSTIVE: std::cout << "exhaustive pairwise matching" << std::endl; break;
      case PAIR_CONTIGUOUS: std::cout << "sequence pairwise matching" << std::endl; break;
      case PAIR_FROM_FILE:  std::cout << "user defined pairwise matching" << std::endl; break;
    }

    Timer timer;
    Matcher_AllInMemory<KeypointSetT, MatcherT> collectionMatcher(fDistRatio);
    if (collectionMatcher.loadData(vec_fileNames, sOutDir))
    {
      // Get pair to match according the matching mode:
      PairsT pairs;
      switch (ePairmode)
      {
        case PAIR_EXHAUSTIVE: pairs = exhaustivePairs(vec_fileNames.size()); break;
        case PAIR_CONTIGUOUS: pairs = contiguousWithOverlap(vec_fileNames.size(), iMatchingVideoMode); break;
        case PAIR_FROM_FILE:
          if(!loadPairs(vec_fileNames.size(), sPredefinedPairList, pairs))
          {
              return EXIT_FAILURE;
          };
          break;
      }
      // Photometric matching of putative pairs
      collectionMatcher.Match(vec_fileNames, pairs, map_PutativesMatches);
      //---------------------------------------
      //-- Export putative matches
      //---------------------------------------
      std::ofstream file (std::string(sOutDir + "/matches.putative.txt").c_str());
      if (file.is_open())
        PairedIndMatchToStream(map_PutativesMatches, file);
      file.close();
    }
    std::cout << "Task done in (s): " << timer.elapsed() << std::endl;
  }
  //-- export putative matches Adjacency matrix
  PairWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
    map_PutativesMatches,
    stlplus::create_filespec(sOutDir, "PutativeAdjacencyMatrix", "svg"));


  //---------------------------------------
  // d. Geometric filtering of putative matches
  //    - AContrario Estimation of the desired geometric model
  //    - Use an upper bound for the a contrario estimated threshold
  //---------------------------------------
  RigWiseMatches map_GeometricMatches;

  ImageCollectionGeometricFilter<FeatureT> collectionGeomFilter;
  if (collectionGeomFilter.loadData(vec_fileNames, sOutDir))
  {
      Timer timer;
      std::cout << std::endl << " - GEOMETRIC FILTERING - " << std::endl;

      // Find to which intrinsic groups each image belong
      std::map < size_t, size_t > map_IntrinsicIdPerImageId;
      std::map < size_t, size_t > map_RigIdPerImageId;
      std::map < size_t, size_t > map_subCamIdPerImageId;

      for (std::vector<openMVG::SfMIO::CameraRigInfo>::const_iterator iter = vec_camImageName.begin();
          iter != vec_camImageName.end(); ++iter)
      {
        const openMVG::SfMIO::CameraRigInfo & camInfo = *iter;

        // Find the index of the camera
        const size_t idx = std::distance((std::vector<openMVG::SfMIO::CameraRigInfo>::const_iterator) vec_camImageName.begin(), iter);

        // to which intrinsic group each image belongs
        map_IntrinsicIdPerImageId[idx] = camInfo.m_intrinsicId;

        // to which rigid rig each image belongs
        map_RigIdPerImageId[idx]       = camInfo.m_rigId;

        // to which subcamera is related image
        map_subCamIdPerImageId[idx]    = camInfo.m_subCameraId;
      }

      // create rig structure using openGV
      translations_t  rigOffsets;
      rotations_t     rigRotations;
      double          averageFocal=1.0e10;

      for(int k=0; k < vec_focalGroup.size(); ++k)
      {
        translation_t   t =  vec_focalGroup[k].m_rigC;
        rotation_t      R =  vec_focalGroup[k].m_R.transpose();

        rigOffsets.push_back(t);
        rigRotations.push_back(R);
        averageFocal = std::min( averageFocal, (double) vec_focalGroup[k].m_focal);
      }

    //  averageFocal /= (double)  vec_focalGroup.size();

      // create the structure with putative match per rigs
      matching::RigWiseMatches  map_Matches_Rig;
      for (PairWiseMatches::const_iterator iter = map_PutativesMatches.begin();
      iter != map_PutativesMatches.end(); ++iter)
      {
        // extract rig id from image id
        const size_t  I = min(map_RigIdPerImageId[iter->first.first], map_RigIdPerImageId[iter->first.second]);
        const size_t  J = max(map_RigIdPerImageId[iter->first.first], map_RigIdPerImageId[iter->first.second]);

        const bool  bRigSwap = ( I != map_RigIdPerImageId[iter->first.first] );

        if( I != J){

          if( !bRigSwap){
            map_Matches_Rig[make_pair(I,J)][iter->first] = iter->second ;
          }
          else
          {
            // swap I and J matches in order to have matches grouped per rig id
            std::vector <matching::IndMatch>   matches;

            matches.clear();

            for(size_t  k(0); k < map_PutativesMatches[iter->first].size(); ++k)
            {
              IndMatch  featurePair;

              featurePair._i = map_PutativesMatches[iter->first][k]._j;
              featurePair._j = map_PutativesMatches[iter->first][k]._i;

              matches.push_back( featurePair );
            }

            map_Matches_Rig[make_pair(I,J)][make_pair(iter->first.second, iter->first.first)] = matches;
          }
        }
      }

      double maxExpectedError = (1.0 - cos(atan(sqrt(2.0) * 4.0 / averageFocal )));

      // Now filter images
      collectionGeomFilter.Filter(
          GeometricFilter_RigEMatrix_AC( maxExpectedError ),
          map_Matches_Rig,
          map_GeometricMatches,
          rigOffsets,
          rigRotations,
          map_IntrinsicIdPerImageId,
          map_subCamIdPerImageId,
          vec_focalGroup
          );

      //-- Perform an additional check to remove pairs with poor overlap
      std::vector<RigWiseMatches::key_type> vec_rigtoRemove;
      for ( size_t i = 0 ; i < map_GeometricMatches.size(); ++i)
      {
          RigWiseMatches::const_iterator iterMap = map_GeometricMatches.begin();
          advance(iterMap, i);

          const PairWiseMatches  map_Matches = iterMap->second;
          std::vector<PairWiseMatches::key_type> vec_toRemove;
          for ( size_t j = 0 ; j < map_Matches.size(); ++j)
          {
            PairWiseMatches::const_iterator iter = map_Matches.begin();
            advance(iter, j);

            const size_t putativePhotometricCount = map_PutativesMatches.find(iter->first)->second.size();
            const size_t putativeGeometricCount = iter->second.size();
            const float ratio = putativeGeometricCount / (float)putativePhotometricCount;
            if (putativeGeometricCount < 50 || ratio < .3f)  {
              // the pair will be removed
              vec_toRemove.push_back(iter->first);
          }

          // -- remove discarded pairs
          for (std::vector<PairWiseMatches::key_type>::const_iterator
            iter =  vec_toRemove.begin(); iter != vec_toRemove.end(); ++iter)
          {
            map_GeometricMatches.at(iterMap->first).erase(*iter);
          }
        }

        // additional check to keep rig matches or not
        if( map_GeometricMatches[iterMap->first].size() < sqrt(rigOffsets.size()) )
          vec_rigtoRemove.push_back(iterMap->first);

      }

      // -- remove discarded pairs
      for (std::vector<RigWiseMatches::key_type>::const_iterator
        iter =  vec_rigtoRemove.begin(); iter != vec_rigtoRemove.end(); ++iter)
      {
        map_GeometricMatches.erase(*iter);
      }

      //---------------------------------------
      //-- Export geometric filtered matches
      //---------------------------------------
      std::ofstream file (string(sOutDir + "/" + sGeometricMatchesFilename).c_str());
      if (file.is_open())
        PairedIndMatchToStream(map_GeometricMatches, file);
      file.close();

      std::cout << "Task done in (s): " << timer.elapsed() << std::endl;

      //-- export Adjacency matrix
      std::cout << "\n Export Adjacency Matrix of the pairwise's geometric matches"
        << std::endl;
      RigWiseMatchingToAdjacencyMatrixSVG(vec_fileNames.size(),
        map_GeometricMatches,
        stlplus::create_filespec(sOutDir, "GeometricAdjacencyMatrix", "svg"));
  }
  return EXIT_SUCCESS;
}
