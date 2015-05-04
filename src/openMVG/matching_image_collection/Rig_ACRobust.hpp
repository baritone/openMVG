
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/multiview/solver_nonCentral_kernel.hpp"
#include "openMVG/multiview/conditioning.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"
#include <limits>

namespace openMVG {
using namespace openMVG::robust;
using namespace opengv;

//-- A contrario Functor to filter putative corresponding points
//--  thanks estimation of the essential matrix.
struct GeometricFilter_RigEMatrix_AC
{
  GeometricFilter_RigEMatrix_AC(
    double dPrecision = std::numeric_limits<double>::infinity(),
    size_t iteration = 4096)
    : m_dPrecision(dPrecision), m_stIteration(iteration) {};

  /// Robust fitting of the rig pose
  void Fit(
    const bearingVectors_t & b1,
    const bearingVectors_t & b2,
    const std::vector<int> & camCorrespondencesRigOne,
    const std::vector<int> & camCorrespondencesRigTwo,
    const translations_t & rigOffsets,
    const rotations_t & rigRotations,
    transformation_t * relativePose,
    std::vector<size_t> & vec_inliers ) const
  {
    vec_inliers.clear();

    // Use the 6 point solver to the pose
    typedef openMVG::noncentral::kernel::GePointSolver SolverType;
    // Define the AContrario adaptor
    typedef ACKernelAdaptorRigPose<  SolverType,
        openMVG::noncentral::kernel::RigAngularError,
        transformation_t>   KernelType;

    KernelType kernel(b1,
                      b2,
                      camCorrespondencesRigOne,
                      camCorrespondencesRigTwo,
                      rigOffsets,
                      rigRotations);

    // Robustly estimation of the Essential matrix and it's precision
    double upper_bound_precision = m_dPrecision;

    std::pair<double,double> acRansacOut = ACRANSAC(kernel, vec_inliers,
        m_stIteration, relativePose, upper_bound_precision, true );

    if (vec_inliers.size() < 2.5 * SolverType::MINIMUM_SAMPLES )  {
        vec_inliers.clear();
    }
  }

  double m_dPrecision;  //upper_bound of the precision
  size_t m_stIteration; //maximal number of used iterations
};

}; // namespace openMVG
