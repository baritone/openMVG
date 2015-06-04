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

  A.resize(5 * Nobs, 3 * (N3D + Nrig + rigOffsets.size() ) );

  C.resize(5 * Nobs, 1);
  C.fill(0.0);
  vec_sign.resize(5 * Nobs + 3);

  const size_t transStart  = 0;
  const size_t pointStart  = transStart + 3*Nrig;
  const size_t rigtsStart  = transStart + 3*(Nrig+N3D);

# define TVAR(i, el) (0 + 3*(i) + (el))
# define XVAR(j, el) (pointStart + 3*(j) + (el))
# define RVAR(n, el) (rigtsStart + 3*(n) + (el))

  // By default set free variable:
  vec_bounds = std::vector< std::pair<double,double> >(3 * (N3D + Nrig + rigOffsets.size()));
  fill( vec_bounds.begin(), vec_bounds.end(), std::make_pair((double)-1e+30, (double)1e+30));
  // Fix the translation ambiguity. (set first cam at (0,0,0))
  vec_bounds[0] = vec_bounds[1] = vec_bounds[2] = std::make_pair(0,0);

  // initialize rig translation value and thresholds
  for( size_t l = 0 ; l < rigOffsets.size(); ++l )
  {
      const Vec3 tc = -rigRotation[l] * rigOffsets[l];
      vec_bounds[RVAR(l,0)] = std::make_pair( tc(0), tc(0));
      vec_bounds[RVAR(l,1)] = std::make_pair( tc(1), tc(1));
      vec_bounds[RVAR(l,2)] = std::make_pair( tc(2), tc(2));
  }

  size_t rowPos = 0;
  // Add the cheirality conditions (R_c*R_i*X_j + R_c*T_i + t_c)_3 + Z_ij >= 1
  for (size_t k = 0; k < Nobs; ++k)
  {
    const size_t indexPt3D = M(2,k);
    const size_t indexCam  = M(3,k);
    const size_t indexRig  = M(4,k);

    const Mat3 & R  = Ri[indexRig];
    const Mat3 & Rc = rigRotation[indexCam];

    const Mat3 & RcRi = Rc * R;

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = RcRi(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = RcRi(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = RcRi(2,2);
    A.coeffRef(rowPos, TVAR(indexRig, 0)) = Rc(2,0);
    A.coeffRef(rowPos, TVAR(indexRig, 1)) = Rc(2,1);
    A.coeffRef(rowPos, TVAR(indexRig, 2)) = Rc(2,2);
    A.coeffRef(rowPos, RVAR(indexCam, 2)) = 1.0;
    C(rowPos) = 0.01;
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
    A.coeffRef(rowPos, RVAR(indexCam, 0)) = 1.0;
    A.coeffRef(rowPos, RVAR(indexCam, 2)) = (sigma-u);
    C(rowPos) = 0.0;
    vec_sign[rowPos] = LP_Constraints::LP_GREATER_OR_EQUAL;
    ++rowPos;

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = RcRi(0,0) - (sigma+u) * RcRi(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = RcRi(0,1) - (sigma+u) * RcRi(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = RcRi(0,2) - (sigma+u) * RcRi(2,2);
    A.coeffRef(rowPos, TVAR(indexRig, 0)) = Rc(0,0) - (sigma+u) * Rc(2,0);
    A.coeffRef(rowPos, TVAR(indexRig, 1)) = Rc(0,1) - (sigma+u) * Rc(2,1);
    A.coeffRef(rowPos, TVAR(indexRig, 2)) = Rc(0,2) - (sigma+u) * Rc(2,2);
    A.coeffRef(rowPos, RVAR(indexCam, 0)) = 1.0;
    A.coeffRef(rowPos, RVAR(indexCam, 2)) = -(sigma+u);
    C(rowPos) = 0.0;
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
    A.coeffRef(rowPos, RVAR(indexCam, 1)) = 1.0;
    A.coeffRef(rowPos, RVAR(indexCam, 2)) = (sigma-v);
    C(rowPos) = 0.0;
    vec_sign[rowPos] = LP_Constraints::LP_GREATER_OR_EQUAL;
    ++rowPos;

    A.coeffRef(rowPos, XVAR(indexPt3D, 0)) = RcRi(1,0) - (sigma+v) * RcRi(2,0);
    A.coeffRef(rowPos, XVAR(indexPt3D, 1)) = RcRi(1,1) - (sigma+v) * RcRi(2,1);
    A.coeffRef(rowPos, XVAR(indexPt3D, 2)) = RcRi(1,2) - (sigma+v) * RcRi(2,2);
    A.coeffRef(rowPos, TVAR(indexRig, 0)) = Rc(1,0) - (sigma+v) * Rc(2,0);
    A.coeffRef(rowPos, TVAR(indexRig, 1)) = Rc(1,1) - (sigma+v) * Rc(2,1);
    A.coeffRef(rowPos, TVAR(indexRig, 2)) = Rc(1,2) - (sigma+v) * Rc(2,2);
    A.coeffRef(rowPos, RVAR(indexCam, 1)) = 1.0;
    A.coeffRef(rowPos, RVAR(indexCam, 2)) = -(sigma+v);
    C(rowPos) = 0.0;
    vec_sign[rowPos] = LP_Constraints::LP_LESS_OR_EQUAL;
    ++rowPos;
  }
# undef RVAR
# undef TVAR
# undef XVAR
}

/// Encode translation and structure linear program
void EncodeRigCiXi(const Mat & M, //Scene representation
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

  A.resize(18 * N3D, 3 * (N3D + Nrig) );

  C.resize(18 * N3D, 1);
  C.fill(0.0);
  vec_sign.resize(18 * N3D + 3);

  const size_t transStart  = 0;
  const size_t pointStart  = transStart + 3*Nrig;

# define TVAR(i, el) (0 + 3*(i) + (el))
# define XVAR(j, el) (pointStart + 3*(j) + (el))

  // By default set free variable:
  vec_bounds = std::vector< std::pair<double,double> >(3 * (N3D + Nrig));
  fill( vec_bounds.begin(), vec_bounds.end(), std::make_pair((double)-1e+30, (double)1e+30));
  // Fix the translation ambiguity. (set first cam at (0,0,0))
  vec_bounds[0] = vec_bounds[1] = vec_bounds[2] = std::make_pair(0,0);

  // compute minimal angle between bearing vectors
  double   minAngle = 1.0e10;
  double   maxAngle = 0.0;

  Vec3 b0;
  Vec3 b1;
  Vec3 b2;

  for (size_t k = 0; k < N3D ; ++k)
  {
      // we assume here that each track is of length 3
      // extract bearing vectors
      b0 << M(0,3*k),   M(1,3*k),   1.0;
      b1 << M(0,3*k+1), M(1,3*k+1), 1.0;
      b2 << M(0,3*k+2), M(1,3*k+2), 1.0;

      // extract rotations
      const Mat3  Rc0  = rigRotation[M(3,3*k+0)].transpose();
      const Mat3  Rc1  = rigRotation[M(3,3*k+1)].transpose();
      const Mat3  Rc2  = rigRotation[M(3,3*k+2)].transpose();

      const Mat3  R0  = Ri[M(4,3*k+0)].transpose();
      const Mat3  R1  = Ri[M(4,3*k+1)].transpose();
      const Mat3  R2  = Ri[M(4,3*k+2)].transpose();

      // compute useful quantities
      Vec3  Rb0 = R0 * Rc0 * b0;
      Vec3  Rb1 = R1 * Rc1 * b1;
      Vec3  Rb2 = R2 * Rc2 * b2;

      // normalize bearing vectors
      Rb0 /= Rb0.norm(); Rb1 /= Rb1.norm(); Rb2 /= Rb2.norm();

      // update minAngle
      const double alpha_0 = acos(Rb0.transpose() * Rb1 );
      const double alpha_1 = acos(Rb0.transpose() * Rb2 );
      const double alpha_2 = acos(Rb1.transpose() * Rb2 );

      minAngle = std::min( minAngle, std::min(alpha_0, std::min(alpha_1, alpha_2) ) );
      maxAngle = std::max( maxAngle, std::max(alpha_0, std::max(alpha_1, alpha_2) ) );
  }

  for( size_t l = 0 ; l < N3D; ++l )
  {
        vec_bounds[XVAR(l,0)] = std::make_pair(1.0 / maxAngle, (double)1e+30);
        vec_bounds[XVAR(l,1)] = std::make_pair(1.0 / maxAngle, (double)1e+30);
        vec_bounds[XVAR(l,2)] = std::make_pair(1.0 / maxAngle, (double)1e+30);
  }

  size_t rowPos = 0;
  // Add the cheirality conditions (R_c*R_i*X_j + R_c*T_i + t_c)_3 + Z_ij >= 1
  for (size_t k = 0; k < N3D ; ++k)
  {
      // we assume here that each track is of length 3
      // extract bearing vectors
      b0 << M(0,3*k),   M(1,3*k),   1.0;
      b1 << M(0,3*k+1), M(1,3*k+1), 1.0;
      b2 << M(0,3*k+2), M(1,3*k+2), 1.0;

      // extract rotations
      const Mat3  Rc0  = rigRotation[M(3,3*k+0)].transpose();
      const Mat3  Rc1  = rigRotation[M(3,3*k+1)].transpose();
      const Mat3  Rc2  = rigRotation[M(3,3*k+2)].transpose();

      const Mat3  R0  = Ri[M(4,3*k+0)].transpose();
      const Mat3  R1  = Ri[M(4,3*k+1)].transpose();
      const Mat3  R2  = Ri[M(4,3*k+2)].transpose();

      // compute useful quantities
      const Vec3  Rb0 = R0 * Rc0 * b0;
      const Vec3  Rb1 = R1 * Rc1 * b1;
      const Vec3  Rb2 = R2 * Rc2 * b2;

      Vec3  R_c0 = R0 * rigOffsets[M(3,3*k+0)];
      Vec3  R_c1 = R1 * rigOffsets[M(3,3*k+1)];
      Vec3  R_c2 = R2 * rigOffsets[M(3,3*k+2)];

      // 3D point index
      const size_t  pointIndex = M(2,3*k);

      // encode matrix
      for( int i=0 ; i < 3 ; ++i )
      {
          A.coeffRef(18*k + i, TVAR(0, i)) =  1.0;
          A.coeffRef(18*k + i, TVAR(1, i)) = -1.0;
          A.coeffRef(18*k + i, XVAR(pointIndex, 0)) =  Rb0(i);
          A.coeffRef(18*k + i, XVAR(pointIndex, 1)) = -Rb1(i);
          C(18*k + i) = sigma - R_c0(i) + R_c1(i);
          vec_sign[18*k + i] = LP_Constraints::LP_LESS_OR_EQUAL;

          A.coeffRef(18*k +3+ i, TVAR(0, i)) =  1.0;
          A.coeffRef(18*k +3+ i, TVAR(1, i)) = -1.0;
          A.coeffRef(18*k +3+ i, XVAR(pointIndex, 0)) =  Rb0(i);
          A.coeffRef(18*k +3+ i, XVAR(pointIndex, 1)) = -Rb1(i);
          C(18*k +3+ i) = -sigma - R_c0(i) + R_c1(i);
          vec_sign[18*k +3+ i] = LP_Constraints::LP_GREATER_OR_EQUAL;

          A.coeffRef(18*k +6+ i, TVAR(0, i)) =  1.0;
          A.coeffRef(18*k +6+ i, TVAR(2, i)) = -1.0;
          A.coeffRef(18*k +6+ i, XVAR(pointIndex, 0)) =  Rb0(i);
          A.coeffRef(18*k +6+ i, XVAR(pointIndex, 2)) = -Rb2(i);
          C(18*k +6+ i) = sigma - R_c0(i) + R_c2(i);
          vec_sign[18*k +6+ i] = LP_Constraints::LP_LESS_OR_EQUAL;

          A.coeffRef(18*k +9+ i, TVAR(0, i)) =  1.0;
          A.coeffRef(18*k +9+ i, TVAR(2, i)) = -1.0;
          A.coeffRef(18*k +9+ i, XVAR(pointIndex, 0)) =  Rb0(i);
          A.coeffRef(18*k +9+ i, XVAR(pointIndex, 2)) = -Rb2(i);
          C(18*k +9+ i) = -sigma - R_c0(i) + R_c2(i);
          vec_sign[18*k +9+ i] = LP_Constraints::LP_GREATER_OR_EQUAL;

          A.coeffRef(18*k +12+ i, TVAR(1, i)) =  1.0;
          A.coeffRef(18*k +12+ i, TVAR(2, i)) = -1.0;
          A.coeffRef(18*k +12+ i, XVAR(pointIndex, 1)) =  Rb1(i);
          A.coeffRef(18*k +12+ i, XVAR(pointIndex, 2)) = -Rb2(i);
          C(18*k +12+ i) = sigma - R_c1(i) + R_c2(i);
          vec_sign[18*k +12+ i] = LP_Constraints::LP_LESS_OR_EQUAL;

          A.coeffRef(18*k +15+ i, TVAR(1, i)) =  1.0;
          A.coeffRef(18*k +15+ i, TVAR(2, i)) = -1.0;
          A.coeffRef(18*k +15+ i, XVAR(pointIndex, 1)) =  Rb1(i);
          A.coeffRef(18*k +15+ i, XVAR(pointIndex, 2)) = -Rb2(i);
          C(18*k +15+ i) = -sigma - R_c1(i) + R_c2(i);
          vec_sign[18*k +15+ i] = LP_Constraints::LP_GREATER_OR_EQUAL;
      }
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
    EncodeRigCiXi(_M, _vec_Ri,
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
