
// Copyright (c) 2010 libmv authors.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.


// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_SOLVER_NONCENTRAL_KERNEL_H_
#define OPENMVG_MULTIVIEW_SOLVER_NONCENTRAL_KERNEL_H_

#include <vector>
#include "openMVG/numeric/numeric.h"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/cameras/PinholeCamera.hpp"
#include <opengv/types.hpp>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/NoncentralRelativeAdapter.hpp>
#include <opengv/triangulation/methods.hpp>
#include <opengv/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.hpp>

namespace openMVG {
namespace noncentral {
namespace kernel {

using namespace std;
using namespace opengv;

/**
 * Six point solver for non central camera system,
 * // [1] "Solutions to minimal generalized relative pose problems".
 * // authors: Stewenius, H., Nister, D., Oskarsson, M., & Astrom K,
 * // Date: 2005:
 * // Conference: Workshop on omnidirectional vision 2005.
 */
struct SixPointSolver {
  enum { MINIMUM_SAMPLES = 9 };
  enum { MAX_MODELS = 1 };
  static void Solve(relative_pose::NoncentralRelativeAdapter & adapter,
                    std::vector<transformation_t> * models,
                    const std::vector<size_t> &indices);
};


// compute reprojection error
struct RigProjError {
  static double Error(size_t sample,
                      const transformation_t & relativePose,
                      relative_pose::NoncentralRelativeAdapter & _adapter)
  {
    // extract pose of cameras
    const Vec3 bearingOne = _adapter.getBearingVector1(sample);
    const Vec3 bearingTwo = _adapter.getBearingVector2(sample);

    const Vec2 x1 = bearingOne.head(2) / bearingOne(2);
    const Vec2 x2 = bearingTwo.head(2) / bearingTwo(2);

    const Mat3 R1 = _adapter.getCamRotation1(sample).transpose();
    const Mat3 R2 = _adapter.getCamRotation2(sample).transpose();

    const Vec3 t1 = - R1 * _adapter.getCamOffset1(sample);
    const Vec3 t2 = - R2 * _adapter.getCamOffset2(sample);

    // retrieve relative pose of rigs
    const translation_t CRig = relativePose.col(3);
    const rotation_t rRig = relativePose.block<3,3>(0,0).transpose();
    const Vec3  tRig = -rRig * CRig;

    // compute relative pose of cameras
    const rotation_t R = R2 * rRig ;
    const translation_t t = R2 * tRig + t2 ;

    // compute 3d point and reprojection error
    const Mat3 K = Mat3::Identity();

    const PinholeCamera cam1(K, R1, t1);
    const PinholeCamera cam2(K, R, t);

    Vec3 X;
    TriangulateDLT(cam1._P, x1, cam2._P, x2, &X);

    return cam1.Residual(X,x1) +  cam2.Residual(X,x2);
  }
};
typedef RigProjError SimpleError;

// compute angular error
struct RigAngularError {
  static double Error(size_t index,
                      const transformation_t & model,
                      relative_pose::NoncentralRelativeAdapter & _adapter)
  {
    // extract rotation and translation from model
    translation_t translation = model.col(3);
    rotation_t rotation = model.block<3,3>(0,0);

    //  initialize variable
    Vec4 p_hom;
    p_hom[3] = 1.0;

    // compute pose
    translation_t cam1Offset = _adapter.getCamOffset1(index);
    rotation_t cam1Rotation = _adapter.getCamRotation1(index);
    translation_t cam2Offset = _adapter.getCamOffset2(index);
    rotation_t cam2Rotation = _adapter.getCamRotation2(index);

    translation_t directTranslation =
        cam1Rotation.transpose() *
        ((translation - cam1Offset) + rotation * cam2Offset);
    rotation_t directRotation =
        cam1Rotation.transpose() * rotation * cam2Rotation;

    _adapter.sett12(directTranslation);
    _adapter.setR12(directRotation);

    transformation_t inverseSolution;
    inverseSolution.block<3,3>(0,0) = directRotation.transpose();
    inverseSolution.col(3) =
        -inverseSolution.block<3,3>(0,0)*directTranslation;

    p_hom.block<3,1>(0,0) =
        opengv::triangulation::triangulate2(_adapter,index);
    bearingVector_t reprojection1 = p_hom.block<3,1>(0,0).normalized();
    bearingVector_t reprojection2 = (inverseSolution * p_hom).normalized();
    bearingVector_t f1 = _adapter.getBearingVector1(index);
    bearingVector_t f2 = _adapter.getBearingVector2(index);

    //bearing-vector based outlier criterium (select threshold accordingly):
    //1-(f1'*f2) = 1-cos(alpha) \in [0:2]
    double reprojError1 = 1.0 - (f1.transpose() * reprojection1);
    double reprojError2 = 1.0 - (f2.transpose() * reprojection2);
    return reprojError1 + reprojError2;
  }
};
typedef RigAngularError SimpleAngularError;

}  // namespace kernel
}  // namespace noncentral
}  // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_SOLVER_NONCENTRAL_KERNEL_H_
