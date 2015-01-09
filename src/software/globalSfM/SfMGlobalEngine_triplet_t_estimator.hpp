
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GLOBAL_SFM_ENGINE_TRIPLET_T_ESTIMATOR_H
#define OPENMVG_GLOBAL_SFM_ENGINE_TRIPLET_T_ESTIMATOR_H

#include "openMVG/numeric/numeric.h"

#include "openMVG/multiview/conditioning.hpp"

// Linear programming solver(s)
#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#ifdef OPENMVG_HAVE_MOSEK
#include "openMVG/linearProgramming/linearProgrammingMOSEK.hpp"
#endif

#include "openMVG/linearProgramming/bisectionLP.hpp"
#include "openMVG/linearProgramming/lInfinityCV/tijsAndXis_From_xi_Ri.hpp"

#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"

namespace openMVG {
namespace trifocal {
namespace kernel {

/// A trifocal tensor seen as 3 projective cameras
struct TrifocalTensorModel {
  Mat34 P1, P2, P3;

  static double Error(const TrifocalTensorModel & t, const Vec2 & pt1, const Vec2 & pt2, const Vec2 & pt3)
  {
    // Triangulate and return the reprojection error
    Triangulation triangulationObj;
    triangulationObj.add(t.P1, pt1);
    triangulationObj.add(t.P2, pt2);
    triangulationObj.add(t.P3, pt3);
    const Vec3 X = triangulationObj.compute();

    //- Return max error as a test
    double pt1ReProj = (Project(t.P1, X) - pt1).squaredNorm();
    double pt2ReProj = (Project(t.P2, X) - pt2).squaredNorm();
    double pt3ReProj = (Project(t.P3, X) - pt3).squaredNorm();

    return std::max(pt1ReProj, std::max(pt2ReProj,pt3ReProj));
  }
};

}  // namespace kernel
}  // namespace trifocal
}  // namespace openMVG

namespace openMVG{

using namespace openMVG::trifocal::kernel;

struct tisXisTrifocalSolver {
  enum { MINIMUM_SAMPLES = 4 };
  enum { MAX_MODELS = 1 };
  // Solve the computation of the tensor.
  static void Solve(
    const Mat &pt0, const Mat & pt1, const Mat & pt2,
    const std::vector<Mat3> & vec_KR, const Mat3 & K, std::vector<TrifocalTensorModel> *P,
    const double ThresholdUpperBound)
  {
    //Build the megaMatMatrix
    const int n_obs = pt0.cols();
    Mat4X megaMat(4, n_obs*3);
    {
      size_t cpt = 0;
      for (size_t i = 0; i  < n_obs; ++i) {

        megaMat.col(cpt) << pt0.col(i)(0), pt0.col(i)(1), (double)i, 0.0;
        ++cpt;

        megaMat.col(cpt) << pt1.col(i)(0), pt1.col(i)(1), (double)i, 1.0;
        ++cpt;

        megaMat.col(cpt) << pt2.col(i)(0), pt2.col(i)(1), (double)i, 2.0;
        ++cpt;
        }
    }
    //-- Solve the LInfinity translation and structure from Rotation and points data.
    std::vector<double> vec_solution((3 + MINIMUM_SAMPLES)*3);

    using namespace openMVG::lInfinityCV;

#ifdef OPENMVG_HAVE_MOSEK
    MOSEK_SolveWrapper LPsolver(static_cast<int>(vec_solution.size()));
#else
    OSI_CLP_SolverWrapper LPsolver(static_cast<int>(vec_solution.size()));
#endif

    Translation_Structure_L1_ConstraintBuilder cstBuilder(vec_KR, megaMat);
    double gamma;
    if (BisectionLP<Translation_Structure_L1_ConstraintBuilder, LP_Constraints_Sparse>(
          LPsolver,
          cstBuilder,
          &vec_solution,
          ThresholdUpperBound,//admissibleResidual,
          0.0, 1e-8, 2, &gamma, false))
    {
      std::vector<Vec3> vec_tis(3);
      vec_tis[0] = Vec3(vec_solution[0], vec_solution[1], vec_solution[2]);
      vec_tis[1] = Vec3(vec_solution[3], vec_solution[4], vec_solution[5]);
      vec_tis[2] = Vec3(vec_solution[6], vec_solution[7], vec_solution[8]);

      TrifocalTensorModel PTemp;
      PTemp.P1 = HStack(vec_KR[0], vec_tis[0]);
      PTemp.P2 = HStack(vec_KR[1], vec_tis[1]);
      PTemp.P3 = HStack(vec_KR[2], vec_tis[2]);

      P->push_back(PTemp);
    }
  }

  // Compute the residual of reprojections
  static double Error(const TrifocalTensorModel & Tensor, const Vec2 & pt0, const Vec2 & pt1, const Vec2 & pt2)
  {
    return TrifocalTensorModel::Error(Tensor, pt0, pt1, pt2);
  }
};

template <typename SolverArg,
          typename ErrorArg,
          typename ModelArg>
class TrifocalKernel_ACRansac_N_tisXis
{
public:
  typedef SolverArg Solver;
  typedef ModelArg  Model;


  TrifocalKernel_ACRansac_N_tisXis(const Mat & x1, const Mat & x2, const Mat & x3,
    const std::vector<Mat3> & vec_KRi, const Mat3 & K,
    const double ThresholdUpperBound)
    : x1_(x1), x2_(x2), x3_(x3), vec_KR_(vec_KRi),
      K_(K), ThresholdUpperBound_(ThresholdUpperBound),
      logalpha0_(log10(M_PI)),
      Kinv_(K.inverse())
  {
    // Normalize points by inverse(K)
    ApplyTransformationToPoints(x1_, Kinv_, &x1n_);
    ApplyTransformationToPoints(x2_, Kinv_, &x2n_);
    ApplyTransformationToPoints(x3_, Kinv_, &x3n_);

    vec_KR_[0] = Kinv_ * vec_KR_[0];
    vec_KR_[1] = Kinv_ * vec_KR_[1];
    vec_KR_[2] = Kinv_ * vec_KR_[2];
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<size_t> &samples, std::vector<Model> *models) const {

    // Create a model from the points
    Solver::Solve(
                  ExtractColumns(x1n_, samples),
                  ExtractColumns(x2n_, samples),
                  ExtractColumns(x3n_, samples),
                  vec_KR_, K_, models, ThresholdUpperBound_);
  }

  double Error(size_t sample, const Model &model) const {
    return ErrorArg::Error(model, x1n_.col(sample), x2n_.col(sample), x3n_.col(sample));
  }

  void Errors(const Model &model, std::vector<double> & vec_errors) const {
    for (size_t sample = 0; sample < x1n_.cols(); ++sample)
      vec_errors[sample] = ErrorArg::Error(model, x1n_.col(sample), x2n_.col(sample), x3n_.col(sample));
  }

  size_t NumSamples() const {
    return x1n_.cols();
  }

  void Unnormalize(Model * model) const {
    // Unnormalize model from the computed conditioning.
    model->P1 = K_ * model->P1;
    model->P2 = K_ * model->P2;
    model->P3 = K_ * model->P3;
  }

  double logalpha0() const {return logalpha0_;}

  double multError() const {return 1.0;}

  Mat3 normalizer1() const {return Kinv_;}
  Mat3 normalizer2() const {return Mat3::Identity();}
  double unormalizeError(double val) const { return sqrt(val) / Kinv_(0,0);}

private:
  const Mat & x1_, & x2_, & x3_;
  Mat x1n_, x2n_, x3n_;
  const Mat3 Kinv_, K_;
  const double logalpha0_;
  const double ThresholdUpperBound_;
  std::vector<Mat3> vec_KR_;

};

} // namespace openMVG

namespace openMVG {
namespace trifocal {
namespace kernel {

/// A trifocal tensor seen as 3 projective cameras
struct rigTrifocalTensorModel {
  Mat3  R1, R2, R3;
  Vec3  t1, t2, t3;

  static double Error(const rigTrifocalTensorModel & t, const Vec2 & pt1, const Vec2 & pt2, const Vec2 & pt3,
        const Vec3 & camIndex,  // rig subcamera index for each observation
        const std::vector<Mat3> & rigRotation, // rig subcamera rotation
        const std::vector<Vec3> & rigOffsets ) // rig subcamera translation
  {
    // extract sucamera rotations and translation
    size_t I = (size_t) camIndex[0];
    size_t J = (size_t) camIndex[1];
    size_t K = (size_t) camIndex[2];

    const Mat3 RI = rigRotation[I];  const Vec3 tI = -RI * rigOffsets[I];
    const Mat3 RJ = rigRotation[J];  const Vec3 tJ = -RJ * rigOffsets[J];
    const Mat3 RK = rigRotation[K];  const Vec3 tK = -RK * rigOffsets[K];

    //compute projection matrices
    const Mat34 P1 = HStack(RI * t.R1, RI * t.t1 + tI);
    const Mat34 P2 = HStack(RJ * t.R2, RJ * t.t2 + tJ);
    const Mat34 P3 = HStack(RK * t.R3, RK * t.t3 + tK);

    // Triangulate and return the reprojection error
    Triangulation triangulationObj;
    triangulationObj.add(P1, pt1);
    triangulationObj.add(P2, pt2);
    triangulationObj.add(P3, pt3);
    const Vec3 X = triangulationObj.compute();

    //- Return max error as a test
    double pt1ReProj = (Project(P1, X) - pt1).squaredNorm();
    double pt2ReProj = (Project(P2, X) - pt2).squaredNorm();
    double pt3ReProj = (Project(P3, X) - pt3).squaredNorm();

    return std::max(pt1ReProj, std::max(pt2ReProj,pt3ReProj));
  }
};

}  // namespace kernel
}  // namespace trifocal
}  // namespace openMVG


namespace openMVG{

using namespace openMVG::trifocal::kernel;

struct rigTisXisTrifocalSolver {
  enum { MINIMUM_SAMPLES = 4 };
  enum { MAX_MODELS = 1 };
  // Solve the computation of the tensor.
  static void Solve(
    const Mat &pt0, const Mat & pt1, const Mat & pt2,
    const std::vector<Mat3> & rigRotation, // rotation of subcameras
    const std::vector<Vec3> & rigOffsets, // optical center of rig subcameras in rig referential frame
    const Mat &camIndex, // subcamera index for each rotation
    const std::vector<Mat3> & vec_KR, std::vector<rigTrifocalTensorModel> *P,
    const double ThresholdUpperBound)
  {
    //Build the megaMatMatrix
    const int n_obs = pt0.cols();
    assert(n_obs == camIndex.cols() );

    Mat megaMat(5, n_obs*3);
    {
      size_t cpt = 0;
      for (size_t i = 0; i  < n_obs; ++i) {

        megaMat.col(cpt) << pt0.col(i)(0), pt0.col(i)(1), (double)i, camIndex.col(i)(0), 0.0;
        ++cpt;

        megaMat.col(cpt) << pt1.col(i)(0), pt1.col(i)(1), (double)i, camIndex.col(i)(1), 1.0;
        ++cpt;

        megaMat.col(cpt) << pt2.col(i)(0), pt2.col(i)(1), (double)i, camIndex.col(i)(2), 2.0;
        ++cpt;
        }
    }
    //-- Solve the LInfinity translation and structure from Rotation and points data.
    std::vector<double> vec_solution((3 + MINIMUM_SAMPLES)*3);

    using namespace openMVG::lInfinityCV;

#ifdef OPENMVG_HAVE_MOSEK
    MOSEK_SolveWrapper LPsolver(static_cast<int>(vec_solution.size()));
#else
    OSI_CLP_SolverWrapper LPsolver(static_cast<int>(vec_solution.size()));
#endif

    Rig_Translation_Structure_L1_ConstraintBuilder cstBuilder(vec_KR, megaMat, rigRotation, rigOffsets);
    double gamma;
    if (BisectionLP<Rig_Translation_Structure_L1_ConstraintBuilder, LP_Constraints_Sparse>(
          LPsolver,
          cstBuilder,
          &vec_solution,
          ThresholdUpperBound,//admissibleResidual,
          0.0, 1e-8, 2, &gamma, false))
    {
      std::vector<Vec3> vec_tis(3);
      vec_tis[0] = Vec3(vec_solution[0], vec_solution[1], vec_solution[2]);
      vec_tis[1] = Vec3(vec_solution[3], vec_solution[4], vec_solution[5]);
      vec_tis[2] = Vec3(vec_solution[6], vec_solution[7], vec_solution[8]);

      rigTrifocalTensorModel PTemp;
      PTemp.R1 = vec_KR[0]; PTemp.t1 = vec_tis[0];
      PTemp.R2 = vec_KR[1]; PTemp.t2 = vec_tis[1];
      PTemp.R3 = vec_KR[2]; PTemp.t3 = vec_tis[2];

      P->push_back(PTemp);
    }
  }

  // Compute the residual of reprojections
  static double Error(const rigTrifocalTensorModel & Tensor, const Vec2 & pt0, const Vec2 & pt1, const Vec2 & pt2,
        const Vec3 & camIndex, const std::vector<Mat3> & rigRotation, const std::vector<Vec3> & rigOffsets)
  {
    return rigTrifocalTensorModel::Error(Tensor, pt0, pt1, pt2, camIndex, rigRotation, rigOffsets);
  }
};

template <typename SolverArg,
          typename ErrorArg,
          typename ModelArg>
class rig_TrifocalKernel_ACRansac_N_tisXis
{
public:
  typedef SolverArg Solver;
  typedef ModelArg  Model;


  rig_TrifocalKernel_ACRansac_N_tisXis(const Mat & x1, const Mat & x2, const Mat & x3,
    const std::vector<Mat3> & vec_KRi,
    const std::vector<Mat3> & rigRotation,
    const std::vector<Vec3> & rigOffsets,
    const Mat & camIndex,
    const double ThresholdUpperBound)
    : x1_(x1), x2_(x2), x3_(x3), vec_KR_(vec_KRi),
      vec_rigRotation_(rigRotation),
      vec_rigOffset_(rigOffsets),
      camIndex_(camIndex),
      ThresholdUpperBound_(ThresholdUpperBound),
      logalpha0_(log10(M_PI))
  {
     //initialize normalized coordinates
     // Normalize points by inverse(K)
     const Mat3 Kinv_ = Mat3::Identity();
     ApplyTransformationToPoints(x1_, Kinv_, &x1n_);
     ApplyTransformationToPoints(x2_, Kinv_, &x2n_);
     ApplyTransformationToPoints(x3_, Kinv_, &x3n_);
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<size_t> &samples, std::vector<Model> *models) const {

    // Create a model from the points
    Solver::Solve(
                  ExtractColumns(x1n_, samples),
                  ExtractColumns(x2n_, samples),
                  ExtractColumns(x3n_, samples),
                  vec_rigRotation_, vec_rigOffset_,
                  ExtractColumns(camIndex_, samples),
                  vec_KR_, models, ThresholdUpperBound_);
  }

  void Errors(const Model &model, std::vector<double> & vec_errors) const {
    for (size_t sample = 0; sample < x1n_.cols(); ++sample)
      vec_errors[sample] = ErrorArg::Error(model, x1n_.col(sample), x2n_.col(sample), x3n_.col(sample),
                               camIndex_.col(sample), vec_rigRotation_, vec_rigOffset_);
  }

  double Error(size_t sample, const Model &model) const {
    return ErrorArg::Error(model, x1n_.col(sample), x2n_.col(sample), x3n_.col(sample),
       camIndex_.col(sample), vec_rigRotation_, vec_rigOffset_);
  }

  size_t NumSamples() const {
    return x1n_.cols();
  }

  void Unnormalize(Model * model) const {
    // Unnormalize model from the computed conditioning.
    // Unnormalize model from the computed conditioning.
    const Mat3 K_ = Mat3::Identity();
    model->R1 = K_ * model->R1;
    model->R2 = K_ * model->R2;
    model->R3 = K_ * model->R3;
    model->t1 = K_ * model->t1;
    model->t2 = K_ * model->t2;
    model->t3 = K_ * model->t3;
  }

  Mat3 normalizer1() const {return Mat3::Identity();}
  Mat3 normalizer2() const {return Mat3::Identity();}
  double unormalizeError(double val) const { return sqrt(val);}

  double logalpha0() const {return logalpha0_;}

  double multError() const {return 1.0;}

private:
  const Mat & x1_, & x2_, & x3_;
  Mat x1n_, x2n_, x3n_;
  const Mat & camIndex_;
  const double logalpha0_;
  const double ThresholdUpperBound_;
  std::vector<Mat3> vec_KR_;
  std::vector<Mat3> vec_rigRotation_;
  std::vector<Vec3> vec_rigOffset_;

};

} // namespace openMVG

#endif // OPENMVG_GLOBAL_SFM_ENGINE_TRIPLET_T_ESTIMATOR_H
