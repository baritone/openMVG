// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_LINFINITY_COMPUTER_VISION_TRANSLATIONANDSTRUCTUREFrom_xi_RI_RIG_H_
#define OPENMVG_LINFINITY_COMPUTER_VISION_TRANSLATIONANDSTRUCTUREFrom_xi_RI_RIG_H_

#include "openMVG/numeric/numeric.h"
#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include <fstream>
#include <utility>
#include <vector>

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

//--
//- Implementation of algorithm from Paper titled :
//- [1] "Multiple-View Geometry under the L_\infty Norm."
//- Author : Fredrik Kahl, Richard Hartley.
//- Date : 9 sept 2008.

//- [2] "Multiple View Geometry and the L_\infty -norm"
//- Author : Fredrik Kahl
//- ICCV 2005.
//--

namespace openMVG   {
namespace lInfinityCV  {

using namespace linearProgramming;

//-- Estimate the translation and the structure
//    from image points coordinates and camera rotations.
//    - Estimation of Ci from Ri and xij
// [1] -> 6.1 Cameras with Known Rotation
//
//    - This implementation Use L1 norm instead of the L2 norm of
//      the paper, it allows to use standard standard LP
//      (simplex) instead of using SOCP (second order cone programming).
//      Implementation by Pierre Moulon
//

/// Encode translation and structure linear program
void EncodeRigTiXi(const Mat & M, //Scene representation
                           const std::vector<Mat3> Ri,
                           const std::vector<Mat3> rigRotation,
                           const std::vector<Vec3> rigOffsets,
                           double sigma, // Start upper bound
                           sRMat & A, Vec & C,
                           std::vector<LP_Constraints::eLP_SIGN> & vec_sign,
                           std::vector<double> & vec_costs,
                           std::vector< std::pair<double,double> > & vec_bounds)
{
  // Build Constraint matrix.
  const size_t Nrig = (size_t) M.row(4).maxCoeff()+1;
  const size_t N3D  = (size_t) M.row(2).maxCoeff()+1;
  const size_t Nobs = M.cols();

  assert(Nrig == Ri.size());

  A.resize(5 * Nobs, 3 * (N3D + Nrig) );

  C.resize(5 * Nobs, 1);
  C.fill(0.0);
  vec_sign.resize(5 * Nobs + 3);

  const size_t transStart  = 0;
  const size_t pointStart  = transStart + 3*Nrig;
# define TVAR(i, el) (0 + 3*(i) + (el))
# define XVAR(j, el) (pointStart + 3*(j) + (el))

  // By default set free variable:
  vec_bounds = std::vector< std::pair<double,double> >(3 * (N3D + Nrig));
  fill( vec_bounds.begin(), vec_bounds.end(), std::make_pair((double)-1e+30, (double)1e+30));
  // Fix the translation ambiguity. (set first cam at (0,0,0))
  vec_bounds[0] = vec_bounds[1] = vec_bounds[2] = std::make_pair(0,0);

  size_t rowPos = 0;
  // Add the cheirality conditions (R_c*R_i*X_j + R_c*T_i + t_c)_3 + Z_ij >= 1
  for (size_t k = 0; k < Nobs; ++k)
  {
    const size_t indexPt3D = M(2,k);
    const size_t indexCam  = M(3,k);
    const size_t indexRig  = M(4,k);

    const Mat3 & R  = Ri[indexRig];
    const Mat3 & Rc = rigRotation[indexCam];
    const Vec3 & tc = -Rc * rigOffsets[indexCam];

    const Mat3 & RcRi = Rc * R;

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = RcRi(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = RcRi(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = RcRi(2,2);
    A.coeffRef(rowPos, TVAR(indexRig, 0)) = Rc(2,0);
    A.coeffRef(rowPos, TVAR(indexRig, 1)) = Rc(2,1);
    A.coeffRef(rowPos, TVAR(indexRig, 2)) = Rc(2,2);
    C(rowPos) = 1.0 - tc(2);
    vec_sign[rowPos] = LP_Constraints::LP_GREATER_OR_EQUAL;
    ++rowPos;

    const Vec2 pt   = M.block<2,1>(0,k);
    const double u = pt(0);
    const double v = pt(1);

    // x-residual =>
    // (R_c*R_i*X_j + R_c*T_i + T_c)_1 / (R_c*R_i*X_j + R_c*T_i + T_c)_3 - u >= -sigma
    // (R_c*R_i*X_j + R_c*T_i + T_c)_1 - u * (R_c*R_i*X_j + R_c*T_i + T_c)_3  + sigma (R_c*R_i*X_j + R_c*T_i + T_c)_3  >= 0.0
    // ((R_c*R_i)_3 * (sigma-u) + (R_c*R_i)_1) * X_j +
    //     + (R_c_3 * (sigma-u) + R_c_1)*t_i + (t_c_1 + t_c_3 * (sigma-u) ) >= 0

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = RcRi(0,0) + (sigma-u) * RcRi(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = RcRi(0,1) + (sigma-u) * RcRi(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = RcRi(0,2) + (sigma-u) * RcRi(2,2);
    A.coeffRef(rowPos, TVAR(indexRig, 0)) = Rc(0,0) + (sigma-u) * Rc(2,0);
    A.coeffRef(rowPos, TVAR(indexRig, 1)) = Rc(0,1) + (sigma-u) * Rc(2,1);
    A.coeffRef(rowPos, TVAR(indexRig, 2)) = Rc(0,2) + (sigma-u) * Rc(2,2);
    C(rowPos) = -tc(0) - (sigma-u) * tc(2);
    vec_sign[rowPos] = LP_Constraints::LP_GREATER_OR_EQUAL;
    ++rowPos;

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = RcRi(0,0) - (sigma+u) * RcRi(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = RcRi(0,1) - (sigma+u) * RcRi(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = RcRi(0,2) - (sigma+u) * RcRi(2,2);
    A.coeffRef(rowPos, TVAR(indexRig, 0)) = Rc(0,0) - (sigma+u) * Rc(2,0);
    A.coeffRef(rowPos, TVAR(indexRig, 1)) = Rc(0,1) - (sigma+u) * Rc(2,1);
    A.coeffRef(rowPos, TVAR(indexRig, 2)) = Rc(0,2) - (sigma+u) * Rc(2,2);
    C(rowPos) = -tc(0) + (sigma+u) * tc(2);
    vec_sign[rowPos] = LP_Constraints::LP_LESS_OR_EQUAL;
    ++rowPos;

    // y-residual =>
    // (R_c*R_i*X_j + R_c*T_i + T_c)_2 / (R_c*R_i*X_j + R_c*T_i + T_c)_3 - v >= -sigma
    // (R_c*R_i*X_j + R_c*T_i + T_c)_2 - v * (R_c*R_i*X_j + R_c*T_i + T_c)_3  + sigma (R_c*R_i*X_j + R_c*T_i + T_c)_3  >= 0.0
    // ((R_c*R_i)_3 * (sigma-v) + (R_c*R_i)_2) * X_j +
    //     + (R_c_3 * (sigma-v) + R_c_2)*t_i + (t_c_2 + t_c_3 * (sigma-u) ) >= 0

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = RcRi(1,0) + (sigma-v) * RcRi(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = RcRi(1,1) + (sigma-v) * RcRi(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = RcRi(1,2) + (sigma-v) * RcRi(2,2);
    A.coeffRef(rowPos, TVAR(indexRig, 0)) = Rc(1,0) + (sigma-v) * Rc(2,0);
    A.coeffRef(rowPos, TVAR(indexRig, 1)) = Rc(1,1) + (sigma-v) * Rc(2,1);
    A.coeffRef(rowPos, TVAR(indexRig, 2)) = Rc(1,2) + (sigma-v) * Rc(2,2);
    C(rowPos) = -tc(1) - (sigma-v) * tc(2);
    vec_sign[rowPos] = LP_Constraints::LP_GREATER_OR_EQUAL;
    ++rowPos;

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = RcRi(1,0) - (sigma+v) * RcRi(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = RcRi(1,1) - (sigma+v) * RcRi(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = RcRi(1,2) - (sigma+v) * RcRi(2,2);
    A.coeffRef(rowPos, TVAR(indexRig, 0)) = Rc(1,0) - (sigma+v) * Rc(2,0);
    A.coeffRef(rowPos, TVAR(indexRig, 1)) = Rc(1,1) - (sigma+v) * Rc(2,1);
    A.coeffRef(rowPos, TVAR(indexRig, 2)) = Rc(1,2) - (sigma+v) * Rc(2,2);
    C(rowPos) = -tc(1) + (sigma+v) * tc(2);
    vec_sign[rowPos] = LP_Constraints::LP_LESS_OR_EQUAL;
    ++rowPos;
  }
# undef TVAR
# undef XVAR
}

/// Kernel that set Linear constraints for the
///   - Translation Registration and Structure Problem.
///  Designed to be used with bisectionLP and LP_Solver interface.
///
/// Implementation of problem of [1] -> 6.1 Cameras with known rotation
///  under a Linear Program form. (With SPARSE constraint matrix).

struct Rig_Translation_Structure_L1_ConstraintBuilder
{
  Rig_Translation_Structure_L1_ConstraintBuilder(
    const std::vector<Mat3> & vec_Ri,
    const Mat & M,
    const std::vector<Mat3> & rigRotation,
    const std::vector<Vec3> & rigOffsets)
  {
    _M = M;
    _vec_Ri = vec_Ri;
    _rigRotation = rigRotation;
    _rigOffsets = rigOffsets;
  }

  /// Setup constraints for the translation and structure problem,
  ///  in the LP_Constraints object.
  bool Build(double gamma, LP_Constraints_Sparse & constraint)
  {
    EncodeRigTiXi(_M, _vec_Ri,
      _rigRotation,
      _rigOffsets,
      gamma,
      constraint._constraintMat,
      constraint._Cst_objective,
      constraint._vec_sign,
      constraint._vec_cost,
      constraint._vec_bounds);

    //-- Setup additional information about the Linear Program constraint
    // We look for nb translations and nb 3D points.
    const size_t N3D  = (size_t) _M.row(2).maxCoeff() + 1;
    const size_t Nrig = (size_t) _M.row(4).maxCoeff() + 1;

    constraint._nbParams = (Nrig + N3D) * 3;

    // to be sure that cost vector is empty
    constraint._vec_cost.clear();

    return true;
  }

  std::vector<Mat3> _vec_Ri;  // Rotation matrix
  Mat _M; // M contains (X,Y,index3dPoint, indexCam)^T
  std::vector<Mat3> _rigRotation; // rotation of rig subcameras
  std::vector<Vec3> _rigOffsets; // optical center of rig subcameras in rig referential frame
};

} // namespace lInfinityCV
} // namespace openMVG

#endif // OPENMVG_LINFINITY_COMPUTER_VISION_TRANSLATIONANDSTRUCTUREFrom_xi_RI_RIG_H_
