
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/features/features.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "software/SfM/SfMIOHelper.hpp"

#include <opengv/types.hpp>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/NoncentralRelativeAdapter.hpp>
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.hpp>

using namespace openMVG;
using namespace opengv;

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

#include <vector>
#include <map>

using namespace openMVG::matching;

template <typename FeatureT>
class ImageCollectionGeometricFilter
{
  public:
  ImageCollectionGeometricFilter()
  {
  }

  /// Load all features in memory
  bool loadData(
    const std::vector<std::string> & vec_fileNames, // input filenames
    const std::string & sMatchDir) // where the data are saved
  {
    bool bOk = true;
    for (size_t j = 0; j < vec_fileNames.size(); ++j)  {
      // Load features of Jnth image
      const std::string sFeatJ = stlplus::create_filespec(sMatchDir,
        stlplus::basename_part(vec_fileNames[j]), "feat");
      bOk &= loadFeatsFromFile(sFeatJ, map_Feat[j]);
    }
    return bOk;
  }

  /// Filter all putative correspondences according the templated geometric filter
  template <typename GeometricFilterT>
  void Filter(
    const GeometricFilterT & geometricFilter,  // geometric filter functor
    RigWiseMatches & map_PutativesMatchesPair, // putative correspondences to filter
    RigWiseMatches & map_GeometricMatches,     // filtered putative matches
    translations_t rigOffsets,                 // rig translations
    rotations_t rigRotations,                   // rig rotations
    const std::map < size_t, size_t > map_IntrinsicIdPerImageId,  // map intrinsic image id -> intrinsic id
    const std::map < size_t, size_t > map_subCamIdPerImageId,     // map image id -> subcam Id
    const std::vector<openMVG::SfMIO::IntrinsicCameraRigInfo> vec_focalGroup // set of intrinsic parameters
    ) const
  {
    C_Progress_display my_progress_bar( map_PutativesMatchesPair.size() );

#ifdef USE_OPENMP
  #pragma  omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < (int)map_PutativesMatchesPair.size(); ++i)
    {

      RigWiseMatches::const_iterator iter = map_PutativesMatchesPair.begin();
      advance(iter,i);

      // initialize structure used for matching between rigs
      bearingVectors_t bearingVectorsRigOne;
      bearingVectors_t bearingVectorsRigTwo;

      std::vector<int>  camCorrespondencesRigOne;
      std::vector<int>  camCorrespondencesRigTwo;

      // initialize structure in order to export inliers
      std::vector < std::pair < size_t, std::pair <size_t, size_t > > >  featureAndImages;

      // loop on putative match between the two rigs
      size_t  cpt = 0;
      for( PairWiseMatches::const_iterator iter_pair = iter->second.begin();
            iter_pair != iter->second.end();
            ++iter_pair )
      {
        const size_t iIndex = iter_pair->first.first;
        const size_t jIndex = iter_pair->first.second;
        const std::vector<IndMatch> & vec_PutativeMatches = iter_pair->second;

        // Load features of Inth and Jnth images
        typename std::map<size_t, std::vector<FeatureT> >::const_iterator iterFeatsI = map_Feat.find(iIndex);
        typename std::map<size_t, std::vector<FeatureT> >::const_iterator iterFeatsJ = map_Feat.find(jIndex);
        const std::vector<FeatureT> & kpSetI = iterFeatsI->second;
        const std::vector<FeatureT> & kpSetJ = iterFeatsJ->second;

        //-- extract camera matrix in order to normalize coordinates
        const Mat3  Ki = vec_focalGroup [ map_IntrinsicIdPerImageId.at( iIndex) ].m_K;
        const Mat3  Kj = vec_focalGroup [ map_IntrinsicIdPerImageId.at( jIndex) ].m_K;

        for (size_t i=0; i < vec_PutativeMatches.size(); ++i)
        {
          // intialize bearing vectors
          bearingVector_t  bearingOne;
          bearingVector_t  bearingTwo;

          // compute and store normalized coordinates for image I
          const FeatureT & imaA = kpSetI[vec_PutativeMatches[i]._i];
          bearingOne(0) = imaA.x();
          bearingOne(1) = imaA.y();
          bearingOne(2) = 1.0;

          // compute and store normalized coordinates for image J
          const FeatureT & imaB = kpSetJ[vec_PutativeMatches[i]._j];
          bearingTwo(0) = imaB.x();
          bearingTwo(1) = imaB.y();
          bearingTwo(2) = 1.0;

          // normalize features
          bearingOne = Ki.inverse() * bearingOne;
          bearingTwo = Kj.inverse() * bearingTwo;

          // normalize bearing vectors
          bearingOne = bearingOne / bearingOne.norm();
          bearingTwo = bearingTwo / bearingTwo.norm();

          // add bearing vectors to list and update correspondences list
          bearingVectorsRigOne.push_back( bearingOne  );
          camCorrespondencesRigOne.push_back( map_IntrinsicIdPerImageId.at( iIndex) );

          bearingVectorsRigTwo.push_back( bearingTwo  );
          camCorrespondencesRigTwo.push_back( map_IntrinsicIdPerImageId.at( jIndex) );

          //update datastructure
          featureAndImages.push_back( std::make_pair(i, iter_pair->first) );
          ++cpt;
        }
      }

      //-- Apply the geometric filter
      {
        std::vector<size_t> vec_inliers;
        geometricFilter.Fit(
          bearingVectorsRigOne,
          bearingVectorsRigTwo,
          camCorrespondencesRigOne,
          camCorrespondencesRigTwo,
          rigOffsets,
          rigRotations,
          vec_inliers);

        if(!vec_inliers.empty())
        {
          // export computed matches. Step one, compute PairWiseMatches structure
          PairWiseMatches   rigInliers;

          for (size_t i=0; i < vec_inliers.size(); ++i)
          {
            const std::pair < size_t, size_t >   IJ_pair = featureAndImages[vec_inliers[i]].second;
            const size_t    featIndex = featureAndImages[vec_inliers[i]].first;

            rigInliers[IJ_pair].push_back( iter->second.at(IJ_pair)[featIndex] );
          }
          #ifdef USE_OPENMP
          #pragma omp critical
          #endif
          {
            map_GeometricMatches[iter->first] = rigInliers;
          }
        }
#ifdef USE_OPENMP
#pragma omp critical
#endif
        {
          ++my_progress_bar;
        }
      }
    }
  }

  private:
  // Features per image
  std::map<size_t, std::vector<FeatureT> > map_Feat;
};
