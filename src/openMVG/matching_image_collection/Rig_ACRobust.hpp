
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <opengv/types.hpp>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/NoncentralRelativeAdapter.hpp>
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.hpp>
#include <../test/time_measurement.hpp>

#include "openMVG/multiview/solver_essential_kernel.hpp"
#include "openMVG/multiview/solver_nonCentral_kernel.hpp"
#include "openMVG/multiview/essential.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"
#include <limits>

namespace openMVG {
using namespace openMVG::robust;

//-- A contrario Functor to filter putative corresponding points
//--  thanks estimation of the essential matrix.
struct GeometricFilter_RigEMatrix_AC
{
  GeometricFilter_RigEMatrix_AC(
    const std::map<size_t, Mat3> & K,
    double dPrecision = std::numeric_limits<double>::infinity(),
    size_t iteration = 4096)
    : m_dPrecision(dPrecision), m_stIteration(iteration), m_K(K) {};

  /// Robust fitting of the ESSENTIAL matrix
  void Fit(
    const std::pair<size_t, size_t> pairIndex,
    const Mat & xA,
    const std::pair<size_t, size_t> & imgSizeA,
    const Mat & xB,
    const std::pair<size_t, size_t> & imgSizeB,
    std::vector<size_t> & vec_inliers) const
  {
    vec_inliers.clear();

    std::map<size_t, Mat3>::const_iterator iterK_I = m_K.find(pairIndex.first);
    std::map<size_t, Mat3>::const_iterator iterK_J = m_K.find(pairIndex.second);
    // Check that intrinsic parameters exist for this pair
    if (iterK_I == m_K.end() || iterK_J == m_K.end() )
      return;

    // Define the AContrario adapted Essential matrix solver
    typedef ACKernelAdaptorEssential<
        openMVG::essential::kernel::FivePointKernel,
        openMVG::fundamental::kernel::EpipolarDistanceError,
        UnnormalizerT,
        Mat3>
        KernelType;

    KernelType kernel(xA, imgSizeA.first, imgSizeA.second,
                      xB, imgSizeB.first, imgSizeB.second,
                      iterK_I->second, iterK_J->second);

    // Robustly estimate the Essential matrix with A Contrario ransac
    Mat3 E;
    double upper_bound_precision = m_dPrecision;
    std::pair<double,double> ACRansacOut =
      ACRANSAC(kernel, vec_inliers, m_stIteration, &E, upper_bound_precision);

    if (vec_inliers.size() < KernelType::MINIMUM_SAMPLES *2.5)  {
      vec_inliers.clear();
    }
  }

  double m_dPrecision;  //upper_bound of the precision
  size_t m_stIteration; //maximal number of used iterations
  std::map<size_t, Mat3> m_K; // K intrinsic matrix per image index
};

}; // namespace openMVG
