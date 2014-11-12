
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
  enum { MINIMUM_SAMPLES = 6 };
  enum { MAX_MODELS = 1 };
  static void Solve(relative_pose::NoncentralRelativeAdapter & adapter,
                    transformation_t & relativePose,
                    const std::vector<int> &indices);
};


// compute reprojection error
struct RigProjError {
  static double Error(const transformation_t & relativePose,
                      const Vec2 &x1, const Mat3 &R1, const Vec3 &t1,
                      const Vec2 &x2, const Mat3 &R2, const Vec3 &t2)
  {
    // retrieve relative pose of rigs
    const translation_t CRig = relativePose.col(3);
    const rotation_t rRig = relativePose.block<3,3>(0,0);
    const Vec3  tRig = -rRig * CRig;

    // compute relative pose of cameras
    const rotation_t R = R2 * rRig.transpose() * R1.transpose() ;
    const translation_t t = -R * t1 + R2 * tRig + t2;

    // compute 3d point and reprojection error
    const Mat3 K = Mat3::Identity();

    const PinholeCamera cam1(K, Mat3::Identity(), Vec3::Zero() );
    const PinholeCamera cam2(K, R, t);

    Vec3 X;
    TriangulateDLT(cam1._P, x1, cam2._P, x2, &X);

    return cam1.Residual(X,x1) +  cam2.Residual(X,x2);
  }
};
typedef RigProjError SimpleError;

}  // namespace kernel
}  // namespace noncentral
}  // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_SOLVER_NONCENTRAL_KERNEL_H_
