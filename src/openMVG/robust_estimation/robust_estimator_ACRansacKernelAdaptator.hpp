
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_ROBUST_ESTIMATOR_ACRANSAC_KERNEL_ADAPTATOR_H_
#define OPENMVG_ROBUST_ESTIMATOR_ACRANSAC_KERNEL_ADAPTATOR_H_

// Here a collection of A contrario Kernel adaptor.
//  - See // [1] "Robust and accurate calibration of camera networks". PhD.
//  - Authors: Pierre MOULON
//
//ACKernelAdaptor
//  i.e. general two view estimation (affine, homography, fundamental)
//ACKernelAdaptorResection
//  For pose/resection estimation
//ACKernelAdaptorResection_K
//  For pose/resection with known camera intrinsic
//ACKernelAdaptorEssential
//  For essential matrix estimation
//ACKernelAdaptor_AngularRadianError
//  i.e. essential matrix estimation between spherical camera
//
// Mainly it add correct data normalization and define the function required
//  by the generic ACRANSAC routine.
//

#include <opengv/types.hpp>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/NoncentralRelativeAdapter.hpp>
#include <opengv/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.hpp>

using namespace opengv;

namespace openMVG {
namespace robust{

/// Two view Kernel adaptator for the A contrario model estimator
/// Handle data normalization and compute the corresponding logalpha 0
///  that depends of the error model (point to line, or point to point)
/// This kernel adaptor is working for affine, homography, fundamental matrix
///  estimation.
template <typename SolverArg,
          typename ErrorArg,
          typename UnnormalizerArg,
          typename ModelArg = Mat3>
class ACKernelAdaptor
{
public:
  typedef SolverArg Solver;
  typedef ModelArg  Model;
  typedef ErrorArg ErrorT;

  ACKernelAdaptor(
    const Mat &x1, int w1, int h1,
    const Mat &x2, int w2, int h2, bool bPointToLine = true)
    : x1_(x1.rows(), x1.cols()), x2_(x2.rows(), x2.cols()),
    N1_(3,3), N2_(3,3), logalpha0_(0.0), bPointToLine_(bPointToLine)
  {
    assert(2 == x1_.rows());
    assert(x1_.rows() == x2_.rows());
    assert(x1_.cols() == x2_.cols());

    NormalizePoints(x1, &x1_, &N1_, w1, h1);
    NormalizePoints(x2, &x2_, &N2_, w2, h2);

    // LogAlpha0 is used to make error data scale invariant
    if(bPointToLine)  {
      // Ratio of containing diagonal image rectangle over image area
      double D = sqrt(w2*(double)w2 + h2*(double)h2); // diameter
      double A = w2*(double)h2; // area
      logalpha0_ = log10(2.0*D/A /N2_(0,0));
    }
    else  {
      // ratio of area : unit circle over image area
      logalpha0_ = log10(M_PI/(w2*(double)h2) /(N2_(0,0)*N2_(0,0)));
    }
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<size_t> &samples, std::vector<Model> *models) const {
    Mat x1 = ExtractColumns(x1_, samples);
    Mat x2 = ExtractColumns(x2_, samples);
    Solver::Solve(x1, x2, models);
  }

  double Error(size_t sample, const Model &model) const {
    return ErrorT::Error(model, x1_.col(sample), x2_.col(sample));
  }

  size_t NumSamples() const {
    return static_cast<size_t>(x1_.cols());
  }

  void Unnormalize(Model * model) const {
    // Unnormalize model from the computed conditioning.
    UnnormalizerArg::Unnormalize(N1_, N2_, &(*model));
  }

  double logalpha0() const {return logalpha0_;}

  double multError() const {return (bPointToLine_)? 0.5 : 1.0;}

  Mat3 normalizer1() const {return N1_;}
  Mat3 normalizer2() const {return N2_;}
  double unormalizeError(double val) const {return sqrt(val) / N2_(0,0);}

private:
  Mat x1_, x2_;       // Normalized input data
  Mat3 N1_, N2_;      // Matrix used to normalize data
  double logalpha0_; // Alpha0 is used to make the error adaptive to the image size
  bool bPointToLine_;// Store if error model is pointToLine or point to point
};

struct UnnormalizerResection {
  static void Unnormalize(const Mat &T, const Mat &U, Mat34 *P){
    *P = T.inverse() * (*P);
  }
};

/// Pose/Resection Kernel adaptator for the A contrario model estimator
template <typename SolverArg,
  typename ErrorArg,
  typename UnnormalizerArg,
  typename ModelArg = Mat34>
class ACKernelAdaptorResection
{
public:
  typedef SolverArg Solver;
  typedef ModelArg  Model;
  typedef ErrorArg ErrorT;

  ACKernelAdaptorResection(const Mat &x2d, int w, int h, const Mat &x3D)
    : x2d_(x2d.rows(), x2d.cols()), x3D_(x3D),
    N1_(3,3), logalpha0_(log10(M_PI))
  {
    assert(2 == x2d_.rows());
    assert(3 == x3D_.rows());
    assert(x2d_.cols() == x3D_.cols());

    NormalizePoints(x2d, &x2d_, &N1_, w, h);
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<size_t> &samples, std::vector<Model> *models) const {
    Mat x1 = ExtractColumns(x2d_, samples);
    Mat x2 = ExtractColumns(x3D_, samples);
    Solver::Solve(x1, x2, models);
  }

  double Error(int sample, const Model &model) const {
    return ErrorT::Error(model, x2d_.col(sample), x3D_.col(sample));
  }

  size_t NumSamples() const { return x2d_.cols(); }

  void Unnormalize(Model * model) const {
    // Unnormalize model from the computed conditioning.
    UnnormalizerArg::Unnormalize(N1_, Mat3::Identity(), &(*model));
  }

  double logalpha0() const {return logalpha0_;}
  double multError() const {return 1.0;} // point to point error
  Mat3 normalizer1() const {return N1_;}
  Mat3 normalizer2() const {return Mat3::Identity();}
  double unormalizeError(double val) const {return sqrt(val) / N1_(0,0);}

private:
  Mat x2d_, x3D_;
  Mat3 N1_;      // Matrix used to normalize data
  double logalpha0_; // Alpha0 is used to make the error adaptive to the image size
};

/// Pose/Resection Kernel adaptator for the A contrario model estimator with
///  known Intrinsic
template <typename SolverArg,
  typename ErrorArg,
  typename UnnormalizerArg,
  typename ModelArg = Mat34>
class ACKernelAdaptorResection_K
{
public:
  typedef SolverArg Solver;
  typedef ModelArg  Model;
  typedef ErrorArg ErrorT;

  ACKernelAdaptorResection_K(const Mat &x2d, const Mat &x3D, const Mat3 & K)
    : x2d_(x2d.rows(), x2d.cols()), x3D_(x3D),
    N1_(K.inverse()),
    logalpha0_(log10(M_PI)), K_(K)
  {
    assert(2 == x2d_.rows());
    assert(3 == x3D_.rows());
    assert(x2d_.cols() == x3D_.cols());

    // Normalize points by inverse(K)
    ApplyTransformationToPoints(x2d, N1_, &x2d_);
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<size_t> &samples, std::vector<Model> *models) const {
    Mat x1 = ExtractColumns(x2d_, samples);
    Mat x2 = ExtractColumns(x3D_, samples);
    Solver::Solve(x1, x2, models);
  }

  double Error(size_t sample, const Model &model) const {
    return ErrorT::Error(model, x2d_.col(sample), x3D_.col(sample));
  }

  size_t NumSamples() const { return x2d_.cols(); }

  void Unnormalize(Model * model) const {
    // Unnormalize model from the computed conditioning.
    (*model) = K_ * (*model);
  }

  double logalpha0() const {return logalpha0_;}
  double multError() const {return 1.0;} // point to point error
  Mat3 normalizer1() const {return N1_;}
  Mat3 normalizer2() const {return Mat3::Identity();}
  double unormalizeError(double val) const {return sqrt(val) / N1_(0,0);}

private:
  Mat x2d_, x3D_;
  Mat3 N1_;      // Matrix used to normalize data
  double logalpha0_; // Alpha0 is used to make the error adaptive to the image size
  Mat3 K_;            // Intrinsic camera parameter
};

/// Essential matrix Kernel adaptator for the A contrario model estimator
template <typename SolverArg,
  typename ErrorArg,
  typename UnnormalizerArg,
  typename ModelArg = Mat3>
class ACKernelAdaptorEssential
{
public:
  typedef SolverArg Solver;
  typedef ModelArg  Model;
  typedef ErrorArg ErrorT;

  ACKernelAdaptorEssential(
    const Mat &x1, int w1, int h1,
    const Mat &x2, int w2, int h2,
    const Mat3 & K1, const Mat3 & K2)
    : x1_(x1), x3D_(x2),
    N1_(Mat3::Identity()), N2_(Mat3::Identity()), logalpha0_(0.0),
    K1_(K1), K2_(K2)
{
    assert(2 == x1_.rows());
    assert(x1_.rows() == x3D_.rows());
    assert(x1_.cols() == x3D_.cols());

    ApplyTransformationToPoints(x1_, K1_.inverse(), &x1k_);
    ApplyTransformationToPoints(x3D_, K2_.inverse(), &x2k_);

    //Point to line probability (line is the epipolar line)
    double D = sqrt(w2*(double)w2 + h2*(double)h2); // diameter
    double A = w2*(double)h2; // area
    logalpha0_ = log10(2.0*D/A * .5);
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<size_t> &samples, std::vector<Model> *models) const {
    Mat x1 = ExtractColumns(x1k_, samples);
    Mat x2 = ExtractColumns(x2k_, samples);
    Solver::Solve(x1, x2, models);
  }

  double Error(size_t sample, const Model &model) const {
    Mat3 F;
    FundamentalFromEssential(model, K1_, K2_, &F);
    return ErrorT::Error(F, this->x1_.col(sample), this->x3D_.col(sample));
  }

  size_t NumSamples() const { return x1_.cols(); }
  void Unnormalize(Model * model) const {}
  double logalpha0() const {return logalpha0_;}
  double multError() const {return 0.5;} // point to line error
  Mat3 normalizer1() const {return N1_;}
  Mat3 normalizer2() const {return N2_;}
  double unormalizeError(double val) const { return val; }

private:
  Mat x1_, x3D_, x1k_, x2k_; // image point and camera plane point.
  Mat3 N1_, N2_;      // Matrix used to normalize data
  double logalpha0_; // Alpha0 is used to make the error adaptive to the image size
  Mat3 K1_, K2_;      // Intrinsic camera parameter
};

/// Essential matrix Kernel adaptator for the A contrario model estimator
template <typename SolverArg,
  typename ErrorArg,
  typename UnnormalizerArg,
  typename ModelArg = Mat34>
class ACKernelAdaptorRigPose
{
public:
  typedef SolverArg Solver;
  typedef ModelArg  Model;
  typedef ErrorArg ErrorT;

  ACKernelAdaptorRigPose(
    const bearingVectors_t & b1,
    const bearingVectors_t & b2,
    const std::vector<int> & scIdOne,
    const std::vector<int> & scIdTwo,
    const translations_t & rigOffsets,
    const rotations_t & rigRotations,
    const size_t & w, const size_t & h )
    : b1_(b1), b2_(b2),
      scIdOne_(scIdOne), scIdTwo_(scIdTwo),
      offSets(rigOffsets), rotations(rigRotations),
      logalpha0_(0.0)
{
    assert(3 == b1_.size());
    assert(b1_.size() == b2_.size());

    //Point to line probability (line is the epipolar line)
    double D = sqrt(w*(double)w + h*(double)h); // diameter
    double A = w*(double)h; // area
    logalpha0_ = log10(2.0*D/A * .5);
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<size_t> &samples,
    relative_pose::NoncentralRelativeAdapter & adapter,
    transformation_t & relativePose)
  const {
    Solver::Solve(adapter, relativePose, samples);
  }

  double Error(size_t sample, const transformation_t & relativePose) const
  {
    const Vec3 bearingOne = b1_[sample];
    const Vec3 bearingTwo = b2_[sample];

    const Vec2 x1 = bearingOne.head(2) / bearingOne(2) ;
    const Vec2 x2 = bearingTwo.head(2) / bearingTwo(2);

    const Mat3 R1 = rotations[scIdOne_[sample]];
    const Mat3 R2 = rotations[scIdTwo_[sample]];

    const Vec3 t1 = - R1 * offSets[scIdOne_[sample]];
    const Vec3 t2 = - R2 * offSets[scIdTwo_[sample]];

    return ErrorT::Error(relativePose, x1, R1, t1, x2, R2, t2);
  }

  size_t NumSamples() const { return b1_.size(); }
  void Unnormalize(Model * model) const {}
  double logalpha0() const {return logalpha0_;}
  double multError() const {return 0.5;} // point to line error
  double unormalizeError(double val) const { return val; }

private:
  bearingVectors_t b1_, b2_; // bearing vectors.
  std::vector<int> scIdOne_, scIdTwo_ ; // cam correspondences
  translations_t offSets ; // camera position in rigid rig frame
  rotations_t  rotations ; // rig subcamera orientation
  double logalpha0_; // Alpha0 is used to make the error adaptive to the image size
};

/// Two view Kernel adaptator for the A contrario model estimator.
/// Specialization to handle radian angular residual error.
template <typename SolverArg,
          typename ErrorArg,
          typename UnnormalizerArg,
          typename ModelArg = Mat3>
class ACKernelAdaptor_AngularRadianError
{
public:
  typedef SolverArg Solver;
  typedef ModelArg  Model;
  typedef ErrorArg ErrorT;

  ACKernelAdaptor_AngularRadianError(
    const Mat & xA,
    const Mat & xB):
    x1_(xA), x2_(xB),
    logalpha0_(log10(1.0/2.0))
  {
    assert(3 == x1_.rows());
    assert(x1_.rows() == x2_.rows());
    assert(x1_.cols() == x2_.cols());
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<size_t> &samples, std::vector<Model> *models) const {
    const Mat x1 = ExtractColumns(x1_, samples);
    const Mat x2 = ExtractColumns(x2_, samples);
    Solver::Solve(x1, x2, models);
  }

  double Error(size_t sample, const Model &model) const {
    return Square(ErrorT::Error(model, x1_.col(sample), x2_.col(sample)));
  }

  size_t NumSamples() const {
    return static_cast<size_t>(x1_.cols());
  }

  void Unnormalize(Model * model) const {
    //-- Do nothing, no normalization in the angular case
  }

  double logalpha0() const {return logalpha0_;}

  double multError() const {return 1./4.;}

  Mat3 normalizer1() const {return Mat3::Identity();}
  Mat3 normalizer2() const {return Mat3::Identity();}
  double unormalizeError(double val) const {return sqrt(val);}

private:
  Mat x1_, x2_;       // Normalized input data
  double logalpha0_; // Alpha0 is used to make the error scale invariant
};

} // namespace robust
} // namespace openMVG
#endif // OPENMVG_ROBUST_ESTIMATOR_ACRANSAC_KERNEL_ADAPTATOR_H_
