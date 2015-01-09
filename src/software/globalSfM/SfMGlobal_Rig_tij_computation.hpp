
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GLOBAL_SFM_ENGINE_RIG_TIJ_COMPUTATION_H
#define OPENMVG_GLOBAL_SFM_ENGINE_RIG_TIJ_COMPUTATION_H

#include "openMVG/numeric/numeric.h"
#include "openMVG/tracks/tracks.hpp"
#include "software/globalSfM/SfMRigidGlobalEngine.hpp"
#include "software/globalSfM/SfMGlobalEngine_triplet_t_estimator.hpp"

#include "openMVG/multiview/essential.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"

#undef DYNAMIC
#include "openMVG/bundle_adjustment/problem_data_container.hpp"
#include "openMVG/bundle_adjustment/rig_pinhole_ceres_functor.hpp"
#include "software/globalSfM/SfMBundleAdjustmentHelper_tonly.hpp"

#include "openMVG/matching/indexed_sort.hpp"

#include "software/globalSfM/mutexSet.hpp"

namespace openMVG{

// Robust estimation and refinement of a translation and 3D points of an image triplets.
bool estimate_T_rig_triplet(
  const openMVG::tracks::STLMAPTracks & map_tracksCommon,
  const std::map<size_t, std::vector<SIOPointFeature> > & map_feats,
  const std::vector<Mat3> & vec_global_KR_Triplet,
  const std::vector<Mat3> & vec_rigRotation,
  const std::vector<Vec3> & vec_rigOffset,
  const std::map<size_t, size_t > & map_intrinsicIdPerImageId,
  const std::vector<Vec3> & subTrackIndex,
  std::vector<Vec3> & vec_tis,
  double & dPrecision, // UpperBound of the precision found by the AContrario estimator
  std::vector<size_t> & vec_inliers,
  const double ThresholdUpperBound, //Threshold used for the trifocal tensor estimation solver used in AContrario Ransac
  const std::string & sOutDirectory,
  const size_t & nI, const size_t & nJ, const size_t & nK)
{
  using namespace linearProgramming;
  using namespace lInfinityCV;

  // Convert data
  Mat x1(2, map_tracksCommon.size());
  Mat x2(2, map_tracksCommon.size());
  Mat x3(2, map_tracksCommon.size());
  Mat camIndex(3, map_tracksCommon.size());

  Mat* xxx[3] = {&x1, &x2, &x3};

  size_t cpt = 0;
  for (STLMAPTracks::const_iterator iterTracks = map_tracksCommon.begin();
    iterTracks != map_tracksCommon.end(); ++iterTracks, ++cpt) {
    const submapTrack & subTrack = iterTracks->second;

    // loop on subtracks
    for (size_t index = 0; index < 3 ; ++index)
    { submapTrack::const_iterator iter = subTrack.begin();
      std::advance(iter, subTrackIndex[cpt][index]);
      const size_t imaIndex = iter->first;
      const size_t featIndex = iter->second;

      // extract features
      const SIOPointFeature & pt = map_feats.find(imaIndex)->second[featIndex];
      xxx[index]->col(cpt)(0) = pt.x();
      xxx[index]->col(cpt)(1)  = pt.y();

      // extract camera indexes
      camIndex.col(cpt)(index) = map_intrinsicIdPerImageId.at(imaIndex);
    }
  }

  using namespace openMVG::trifocal;
  using namespace openMVG::trifocal::kernel;

  typedef rig_TrifocalKernel_ACRansac_N_tisXis<
    rigTisXisTrifocalSolver,
    rigTisXisTrifocalSolver,
    rigTrifocalTensorModel> KernelType;
  KernelType kernel(x1, x2, x3, vec_global_KR_Triplet, vec_rigRotation,
                    vec_rigOffset, camIndex, ThresholdUpperBound);

  const size_t ORSA_ITER = 4096;

  rigTrifocalTensorModel T;
  dPrecision = dPrecision ;//std::numeric_limits<double>::infinity();
  std::pair<double,double> acStat = robust::ACRANSAC(kernel, vec_inliers, ORSA_ITER, &T, dPrecision, false);
  dPrecision = acStat.first;

  //-- Export data in order to have an idea of the precision of the estimates
  vec_tis.resize(3);
  vec_tis[0] = T.t1;
  vec_tis[1] = T.t2;
  vec_tis[2] = T.t3;

  // Fill Xis
  std::vector<double> vec_residuals;
  std::vector<Vec3>   vec_Xis;
  std::vector<size_t> vec_inliers_cleaned;

  for (size_t i = 0; i < vec_inliers.size(); ++i)  {

    // extract subcamera rotations and translation
    size_t I = (size_t) camIndex.col(vec_inliers[i])(0);
    size_t J = (size_t) camIndex.col(vec_inliers[i])(1);
    size_t K = (size_t) camIndex.col(vec_inliers[i])(2);

    const Mat3 RI = vec_rigRotation[I];  const Vec3 tI = -RI * vec_rigOffset[I];
    const Mat3 RJ = vec_rigRotation[J];  const Vec3 tJ = -RJ * vec_rigOffset[J];
    const Mat3 RK = vec_rigRotation[K];  const Vec3 tK = -RK * vec_rigOffset[K];

    //compute projection matrices
    const Mat34 P1 = HStack(RI * T.R1, RI * T.t1 + tI);
    const Mat34 P2 = HStack(RJ * T.R2, RJ * T.t2 + tJ);
    const Mat34 P3 = HStack(RK * T.R3, RK * T.t3 + tK);

    Triangulation triangulation;
    triangulation.add(P1, x1.col(vec_inliers[i]));
    triangulation.add(P2, x2.col(vec_inliers[i]));
    triangulation.add(P3, x3.col(vec_inliers[i]));
    const double residuals = triangulation.error();
    const Vec3 Xs = triangulation.compute();

    // keep only good 3D points for BA
    if (triangulation.minDepth() > 0 && is_finite(Xs[0]) && is_finite(Xs[1])
        || is_finite(Xs[2]) && residuals < ThresholdUpperBound )
    {
        vec_Xis.push_back(Xs);
        vec_residuals.push_back(residuals);
        vec_inliers_cleaned.push_back( vec_inliers[i] );
    }
  }

  vec_inliers.swap(vec_inliers_cleaned);

  double min, max, mean, median;
  minMaxMeanMedian<double>(vec_residuals.begin(), vec_residuals.end(),
    min, max, mean, median);

  bool bTest(vec_inliers.size() > 0.25 * map_tracksCommon.size() );

  if (!bTest)
  {
    std::cout << "Triplet rejected : AC: " << dPrecision
      << " median: " << median
      << " inliers count " << vec_inliers.size()
      << " total putative " << map_tracksCommon.size() << std::endl;
  }

  bool bRefine = true;
  if (bRefine && bTest)
  {
    // BA on tis, Xis
    const size_t nbRigs = 3;
    const size_t nbCams = vec_rigRotation.size();
    const size_t nbPoints3D = vec_Xis.size();

    // Count the number of measurement (sum of the reconstructed track length)
    const size_t nbmeasurements = nbPoints3D * 3;

    // Setup a BA problem
    using namespace openMVG::bundle_adjustment;
    BA_Problem_data_rigMotionAndIntrinsic<6,6,3> ba_problem; // Will refine [Rotations|Translations] and 3D points

    // Configure the size of the problem
    ba_problem.num_rigs_ = nbRigs;
    ba_problem.num_cameras_ = nbCams;
    ba_problem.num_intrinsics_ = nbCams;
    ba_problem.num_points_ = nbPoints3D;
    ba_problem.num_observations_ = nbmeasurements;

    ba_problem.rig_index_extrinsic.reserve(ba_problem.num_observations_);
    ba_problem.point_index_.reserve(ba_problem.num_observations_);
    ba_problem.camera_index_extrinsic.reserve(ba_problem.num_observations_);
    ba_problem.camera_index_intrinsic.reserve(ba_problem.num_observations_);
    ba_problem.observations_.reserve(2 * ba_problem.num_observations_);

    ba_problem.num_parameters_ =
      6 * ba_problem.num_rigs_         // rigs rotations / translations
      + 6 * ba_problem.num_cameras_    // #[Rotation|translation] = [3x1]|[3x1]
      + 3 * ba_problem.num_intrinsics_ // cameras intrinsics (focal and principal point)
      + 3 * ba_problem.num_points_;    // 3DPoints = [3x1]
    ba_problem.parameters_.reserve(ba_problem.num_parameters_);

    // Fill camera
    std::vector<double> vec_Rot(vec_Xis.size()*3, 0.0);
    {
      Mat3 R = vec_global_KR_Triplet[0];
      double angleAxis[3];
      ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
      vec_Rot[0] = angleAxis[0];
      vec_Rot[1] = angleAxis[1];
      vec_Rot[2] = angleAxis[2];

      // translation
      ba_problem.parameters_.push_back(angleAxis[0]);
      ba_problem.parameters_.push_back(angleAxis[1]);
      ba_problem.parameters_.push_back(angleAxis[2]);
      ba_problem.parameters_.push_back(vec_tis[0](0));
      ba_problem.parameters_.push_back(vec_tis[0](1));
      ba_problem.parameters_.push_back(vec_tis[0](2));
    }
    {
      Mat3 R = vec_global_KR_Triplet[1];
      double angleAxis[3];
      ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
      vec_Rot[3] = angleAxis[0];
      vec_Rot[4] = angleAxis[1];
      vec_Rot[5] = angleAxis[2];

      // translation
      ba_problem.parameters_.push_back(angleAxis[0]);
      ba_problem.parameters_.push_back(angleAxis[1]);
      ba_problem.parameters_.push_back(angleAxis[2]);
      ba_problem.parameters_.push_back(vec_tis[1](0));
      ba_problem.parameters_.push_back(vec_tis[1](1));
      ba_problem.parameters_.push_back(vec_tis[1](2));
    }
    {
      Mat3 R = vec_global_KR_Triplet[2];
      double angleAxis[3];
      ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
      vec_Rot[6] = angleAxis[0];
      vec_Rot[7] = angleAxis[1];
      vec_Rot[8] = angleAxis[2];

      // translation
      ba_problem.parameters_.push_back(angleAxis[0]);
      ba_problem.parameters_.push_back(angleAxis[1]);
      ba_problem.parameters_.push_back(angleAxis[2]);
      ba_problem.parameters_.push_back(vec_tis[2](0));
      ba_problem.parameters_.push_back(vec_tis[2](1));
      ba_problem.parameters_.push_back(vec_tis[2](2));
    }

    // Setup rig camera position parameters
    for (size_t iter=0; iter < ba_problem.num_cameras_ ; ++iter )
    {
      const Mat3 R = vec_rigRotation[iter];
      double angleAxis[3];
      ceres::RotationMatrixToAngleAxis((const double*)R.data(), angleAxis);
      // translation
      const Vec3 t = -R * vec_rigOffset[iter];
      ba_problem.parameters_.push_back(angleAxis[0]);
      ba_problem.parameters_.push_back(angleAxis[1]);
      ba_problem.parameters_.push_back(angleAxis[2]);
      ba_problem.parameters_.push_back(t[0]);
      ba_problem.parameters_.push_back(t[1]);
      ba_problem.parameters_.push_back(t[2]);
    }

    // Setup rig camera intrinsics parameters
    for (size_t iterCam=0; iterCam < ba_problem.num_cameras_ ; ++iterCam )
    {
      ba_problem.parameters_.push_back( 1.0 );   // FOCAL LENGTH
      ba_problem.parameters_.push_back( 0.0 );   // PRINCIPAL POINT
      ba_problem.parameters_.push_back( 0.0 );   // PRINCIPAL POINT
    }

    // Fill 3D points
    for (std::vector<Vec3>::const_iterator iter = vec_Xis.begin();
      iter != vec_Xis.end();
      ++iter)
    {
      const Vec3 & pt3D = *iter;
      ba_problem.parameters_.push_back(pt3D[0]);
      ba_problem.parameters_.push_back(pt3D[1]);
      ba_problem.parameters_.push_back(pt3D[2]);
    }

    // Fill the measurements
    for (size_t i = 0; i < vec_inliers.size(); ++i)
    {
      double ppx = 0.0, ppy = 0.0;
      Vec2 ptFeat = x1.col(vec_inliers[i]);
      ba_problem.observations_.push_back( ptFeat.x() - ppx );
      ba_problem.observations_.push_back( ptFeat.y() - ppy );

      ba_problem.point_index_.push_back(i);
      ba_problem.camera_index_extrinsic.push_back(camIndex.col(vec_inliers[i])(0));
      ba_problem.camera_index_intrinsic.push_back(camIndex.col(vec_inliers[i])(0));
      ba_problem.rig_index_extrinsic.push_back(0);

      ptFeat = x2.col(vec_inliers[i]);
      ba_problem.observations_.push_back( ptFeat.x() - ppx );
      ba_problem.observations_.push_back( ptFeat.y() - ppy );

      ba_problem.point_index_.push_back(i);
      ba_problem.camera_index_extrinsic.push_back(camIndex.col(vec_inliers[i])(1));
      ba_problem.camera_index_intrinsic.push_back(camIndex.col(vec_inliers[i])(1));
      ba_problem.rig_index_extrinsic.push_back(1);

      ptFeat = x3.col(vec_inliers[i]);
      ba_problem.observations_.push_back( ptFeat.x() - ppx );
      ba_problem.observations_.push_back( ptFeat.y() - ppy );

      ba_problem.point_index_.push_back(i);
      ba_problem.camera_index_extrinsic.push_back(camIndex.col(vec_inliers[i])(2));
      ba_problem.camera_index_intrinsic.push_back(camIndex.col(vec_inliers[i])(2));
      ba_problem.rig_index_extrinsic.push_back(2);
    }

    // Create residuals for each observation in the bundle adjustment problem. The
    // parameters for cameras and points are added automatically.
    ceres::Problem problem;
    ceres::LossFunction * p_LossFunction = new ceres::HuberLoss(Square(2.0));
    for (size_t i = 0; i < ba_problem.num_observations(); ++i) {
      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.

      ceres::CostFunction* cost_function =
         new ceres::AutoDiffCostFunction<rig_pinhole_reprojectionError::ErrorFunc_Refine_Rig_Motion_3DPoints, 2, 3, 6, 6, 3>(
           new rig_pinhole_reprojectionError::ErrorFunc_Refine_Rig_Motion_3DPoints(
               &ba_problem.observations()[2 * i]));

      problem.AddResidualBlock(cost_function,
        p_LossFunction,
        ba_problem.mutable_camera_intrinsic_for_observation(i),
        ba_problem.mutable_camera_extrinsic_for_observation(i),
        ba_problem.mutable_rig_extrinsic_for_observation(i),
        ba_problem.mutable_point_for_observation(i));

      // fix intrinsic rig parameters
      problem.SetParameterBlockConstant(
          ba_problem.mutable_camera_extrinsic_for_observation(i) );
      problem.SetParameterBlockConstant(
          ba_problem.mutable_camera_intrinsic_for_observation(i) );
    }

    // Configure constant parameters (if any)
    {
      std::vector<int> vec_constant_extrinsic; // [R|t]

      vec_constant_extrinsic.push_back(0);
      vec_constant_extrinsic.push_back(1);
      vec_constant_extrinsic.push_back(2);

      for (size_t iExtrinsicId = 0; iExtrinsicId < ba_problem.num_rigs_; ++iExtrinsicId)
      {
        if (!vec_constant_extrinsic.empty())
        {
          ceres::SubsetParameterization *subset_parameterization =
            new ceres::SubsetParameterization(6, vec_constant_extrinsic);
          problem.SetParameterization(ba_problem.mutable_rig_extrinsic(iExtrinsicId),
            subset_parameterization);
        }
      }
    }

    // fix rig one position
    problem.SetParameterBlockConstant(
      ba_problem.mutable_rig_extrinsic(0) );

    // Configure a BA engine and run it
    //  Make Ceres automatically detect the bundle structure.
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::SPARSE_SCHUR;
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
    options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
    else
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
      options.sparse_linear_algebra_library_type = ceres::CX_SPARSE;
    else
    {
      // No sparse backend for Ceres.
      // Use dense solving
      options.linear_solver_type = ceres::DENSE_SCHUR;
    }

    options.minimizer_progress_to_stdout = false;
    options.logging_type = ceres::SILENT;
#ifdef USE_OPENMP
    options.num_threads = omp_get_max_threads();
    options.num_linear_solver_threads = omp_get_max_threads();
#endif // USE_OPENMP

    // Solve BA
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    // If convergence and no error, get back refined parameters
    if (summary.IsSolutionUsable())
    {
      size_t i = 0;
      // Get back updated cameras
      Vec3 * tt[3] = {&vec_tis[0], &vec_tis[1], &vec_tis[2]};
      for (i=0; i < 3; ++i)
      {
        const double * cam = ba_problem.mutable_rig_extrinsic(i);

        (*tt[i]) = Vec3(cam[3], cam[4], cam[5]);
      }

      // Get back 3D points
      size_t k = 0;
      std::vector<Vec3>  finalPoint;

      for (std::vector<Vec3>::iterator iter = vec_Xis.begin();
        iter != vec_Xis.end(); ++iter, ++k)
      {
        const double * pt = ba_problem.mutable_points() + k*3;
        Vec3 & pt3D = *iter;
        pt3D = Vec3(pt[0], pt[1], pt[2]);
        finalPoint.push_back(pt3D);
      }

      // export point cloud (for debug purpose only)
      // std::ostringstream pairIJK;
      // pairIJK << nI << "_" << nJ << "_" << nK << ".ply";

      // plyHelper::exportToPly(finalPoint, stlplus::create_filespec(sOutDirectory,
      //                 "pointCloud_triplet_t_"+pairIJK.str()) );

    }
  }
  return bTest;
}

//-- Perform a trifocal estimation of the graph contain in vec_triplets with an
// edge coverage algorithm. It's complexity is sub-linear in term of edges count.
void GlobalRigidReconstructionEngine::computePutativeTranslation_EdgesCoverage(
  const std::map<size_t, Mat3> & map_globalR,
  const std::vector< graphUtils::Triplet > & vec_triplets,
  std::vector<openMVG::relativeInfo > & vec_initialEstimates,
  matching::RigWiseMatches & newpairMatches) const
{
  // The same K matrix is used by all the camera
  const Mat3 _K = Mat3::Identity();

  //-- Prepare global rotations
  std::map<size_t, Mat3> map_global_KR;
  for (std::map<std::size_t, Mat3>::const_iterator iter = map_globalR.begin();
    iter != map_globalR.end(); ++iter)
  {
    map_global_KR[iter->first] = _K * iter->second;
  }

  // create rig structure using openGV
  std::vector<Vec3>  rigOffsets;
  std::vector<Mat3>  rigRotations;
  double          averageFocal=0.0;

  for(int k=0; k < _vec_intrinsicGroups.size(); ++k)
  {
      const Vec3 t = _vec_intrinsicGroups[k].m_rigC;
      const Mat3 R = _vec_intrinsicGroups[k].m_R;

      rigOffsets.push_back(t);
      rigRotations.push_back(R);
      averageFocal += _vec_intrinsicGroups[k].m_focal ;
  }

  averageFocal /= (double) _vec_intrinsicGroups.size();

  //-- Prepare tracks count per triplets:
  std::map<size_t, size_t> map_tracksPerTriplets;
#ifdef USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int i = 0; i < (int)vec_triplets.size(); ++i)
  {
    const graphUtils::Triplet & triplet = vec_triplets[i];
    const size_t I = triplet.i, J = triplet.j , K = triplet.k;

    RigWiseMatches map_matchesIJK;
    if(_map_Matches_Rig.find(std::make_pair(I,J)) != _map_Matches_Rig.end())
      map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(I,J)));
    else
    if(_map_Matches_Rig.find(std::make_pair(J,I)) != _map_Matches_Rig.end())
      map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(J,I)));

    if(_map_Matches_Rig.find(std::make_pair(I,K)) != _map_Matches_Rig.end())
      map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(I,K)));
    else
    if(_map_Matches_Rig.find(std::make_pair(K,I)) != _map_Matches_Rig.end())
      map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(K,I)));

    if(_map_Matches_Rig.find(std::make_pair(J,K)) != _map_Matches_Rig.end())
      map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(J,K)));
    else
    if(_map_Matches_Rig.find(std::make_pair(K,J)) != _map_Matches_Rig.end())
      map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(K,J)));

    // Compute tracks:
    openMVG::tracks::STLMAPTracks map_tracks;
    TracksBuilder tracksBuilder;
    {
      tracksBuilder.Build(map_matchesIJK);
      tracksBuilder.Filter(_map_RigIdPerImageId,3);
      tracksBuilder.ExportToSTL(map_tracks);
    }
#ifdef USE_OPENMP
  #pragma omp critical
#endif
   {
    map_tracksPerTriplets[i] = map_tracks.size();
   }
  }

  typedef std::pair<size_t,size_t> myEdge;

  //-- List all edges
  std::set<myEdge > set_edges;

  for (size_t i = 0; i < vec_triplets.size(); ++i)
  {
    const graphUtils::Triplet & triplet = vec_triplets[i];
    const size_t I = triplet.i, J = triplet.j , K = triplet.k;
    // Add three edges
    set_edges.insert(std::make_pair(std::min(I,J), std::max(I,J)));
    set_edges.insert(std::make_pair(std::min(I,K), std::max(I,K)));
    set_edges.insert(std::make_pair(std::min(J,K), std::max(J,K)));
  }

  // Copy them in vector in order to try to compute them in parallel
  std::vector<myEdge > vec_edges(set_edges.begin(), set_edges.end());

  MutexSet<myEdge> m_mutexSet;

  std::cout << std::endl
    << "Computation of the relative translations over the graph with an edge coverage algorithm" << std::endl;
#ifdef USE_OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int k = 0; k < vec_edges.size(); ++k)
  {
    const myEdge & edge = vec_edges[k];
    //-- If current edge already computed continue
    if (m_mutexSet.isDiscarded(edge) || m_mutexSet.size() == vec_edges.size())
    {
      std::cout << "EDGES WAS PREVIOUSLY COMPUTED" << std::endl;
      continue;
    }

    std::vector<size_t> vec_possibleTriplets;
    // Find the triplet that contain the given edge
    for (size_t i = 0; i < vec_triplets.size(); ++i)
    {
      const graphUtils::Triplet & triplet = vec_triplets[i];
      if (triplet.contain(edge))
      {
        vec_possibleTriplets.push_back(i);
      }
    }

    //-- Sort the triplet according the number of matches they have on their edges
    std::vector<size_t> vec_commonTracksPerTriplets;
    for (size_t i = 0; i < vec_possibleTriplets.size(); ++i)
    {
      vec_commonTracksPerTriplets.push_back(map_tracksPerTriplets[vec_possibleTriplets[i]]);
    }
    //-- If current edge already computed continue
    if (m_mutexSet.isDiscarded(edge))
      continue;

    using namespace indexed_sort;
    std::vector< sort_index_packet_descend < size_t, size_t> > packet_vec(vec_commonTracksPerTriplets.size());
    sort_index_helper(packet_vec, &vec_commonTracksPerTriplets[0]);

    std::vector<size_t> vec_possibleTripletsSorted;
    for (size_t i = 0; i < vec_commonTracksPerTriplets.size(); ++i) {
      vec_possibleTripletsSorted.push_back( vec_possibleTriplets[packet_vec[i].index] );
    }
    vec_possibleTriplets.swap(vec_possibleTripletsSorted);

    // Try to solve the triplets
    // Search the possible triplet:
    for (size_t i = 0; i < vec_possibleTriplets.size(); ++i)
    {
      const graphUtils::Triplet & triplet = vec_triplets[vec_possibleTriplets[i]];
      const size_t I = triplet.i, J = triplet.j , K = triplet.k;
      {
        RigWiseMatches map_matchesIJK;
        if(_map_Matches_Rig.find(std::make_pair(I,J)) != _map_Matches_Rig.end())
          map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(I,J)));
        else
        if(_map_Matches_Rig.find(std::make_pair(J,I)) != _map_Matches_Rig.end())
          map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(J,I)));

        if(_map_Matches_Rig.find(std::make_pair(I,K)) != _map_Matches_Rig.end())
          map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(I,K)));
        else
        if(_map_Matches_Rig.find(std::make_pair(K,I)) != _map_Matches_Rig.end())
          map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(K,I)));

        if(_map_Matches_Rig.find(std::make_pair(J,K)) != _map_Matches_Rig.end())
          map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(J,K)));
        else
        if(_map_Matches_Rig.find(std::make_pair(K,J)) != _map_Matches_Rig.end())
          map_matchesIJK.insert(*_map_Matches_Rig.find(std::make_pair(K,J)));

        // Select common point:
        STLMAPTracks map_tracksCommon;
        TracksBuilder tracksBuilder;
        {
          tracksBuilder.Build(map_matchesIJK);
          tracksBuilder.Filter(_map_RigIdPerImageId,3);
          tracksBuilder.ExportToSTL(map_tracksCommon);
        }

        // extract associated subcamera id for each tracks
        std::vector<Vec3> subTrackIndex;

        size_t cpt = 0;
        for (STLMAPTracks::const_iterator iterTracks = map_tracksCommon.begin();
          iterTracks != map_tracksCommon.end(); ++iterTracks, ++cpt) {
          const submapTrack & subTrack = iterTracks->second;
          size_t index = 0;
          size_t subTrackCpt = 0;
          Vec3  rigIndex = -1.0*Vec3::Ones();
          Vec3  subTrackId;
          for (submapTrack::const_iterator iter = subTrack.begin(); iter != subTrack.end(); ++iter, ++subTrackCpt) {
            const size_t imaIndex = iter->first;
            const size_t rigidId = _map_RigIdPerImageId.at(imaIndex);
            if( rigIndex[0] == -1 && index == 0 )
            {
              rigIndex[index]    = rigidId;
              subTrackId[index]  = subTrackCpt;
              ++index;
            }
            if( rigIndex[1] == -1 && rigIndex[0] != rigidId && index == 1 )
            {
              rigIndex[index]    = rigidId;
              subTrackId[index]  = subTrackCpt;
              ++index;
            }
            if( rigIndex[2] == -1 && rigIndex[1] != rigidId && rigIndex[0] != rigidId && index == 2 )
            {
              rigIndex[index]    = rigidId;
              subTrackId[index]  = subTrackCpt;
              ++index;
            }
          }
          subTrackIndex.push_back(subTrackId);
        }

        //--
        // Try to estimate this triplet.
        //--
        // Get rotations:
        std::vector<Mat3> vec_global_KR_Triplet;
        vec_global_KR_Triplet.push_back(map_global_KR.at(I));
        vec_global_KR_Triplet.push_back(map_global_KR.at(J));
        vec_global_KR_Triplet.push_back(map_global_KR.at(K));

        // update precision to have good value for normalized coordinates
        double dPrecision = 4.0 / averageFocal / averageFocal;
        const double ThresholdUpperBound = 1.0 / averageFocal;

        std::vector<Vec3> vec_tis(3);
        std::vector<size_t> vec_inliers;

        if (map_tracksCommon.size() > 50 * rigOffsets.size() &&
            estimate_T_rig_triplet(
                  map_tracksCommon, _map_feats_normalized,  vec_global_KR_Triplet,
                  rigRotations, rigOffsets, _map_IntrinsicIdPerImageId, subTrackIndex,
                  vec_tis, dPrecision, vec_inliers, ThresholdUpperBound, _sOutDirectory, I, J, K) )
        {
          std::cout << dPrecision * averageFocal << "\t" << vec_inliers.size() << std::endl;

          //-- Build the three camera:
          const Mat3 RI = map_globalR.find(I)->second;
          const Mat3 RJ = map_globalR.find(J)->second;
          const Mat3 RK = map_globalR.find(K)->second;
          const Vec3 ti = vec_tis[0];
          const Vec3 tj = vec_tis[1];
          const Vec3 tk = vec_tis[2];

          // Build the 3 relative translations estimations.
          // IJ, JK, IK

          //--- ATOMIC
          #ifdef USE_OPENMP
             #pragma omp critical
          #endif
          {
            Mat3 RijGt;
            Vec3 tij;
            RelativeCameraMotion(RI, ti, RJ, tj, &RijGt, &tij);
            vec_initialEstimates.push_back(
              std::make_pair(std::make_pair(I, J), std::make_pair(RijGt, tij)));

            Mat3 RjkGt;
            Vec3 tjk;
            RelativeCameraMotion(RJ, tj, RK, tk, &RjkGt, &tjk);
            vec_initialEstimates.push_back(
              std::make_pair(std::make_pair(J, K), std::make_pair(RjkGt, tjk)));

            Mat3 RikGt;
            Vec3 tik;
            RelativeCameraMotion(RI, ti, RK, tk, &RikGt, &tik);
            vec_initialEstimates.push_back(
              std::make_pair(std::make_pair(I, K), std::make_pair(RikGt, tik)));

            // Add trifocal inliers as valid 3D points
          }

          //-- Remove the 3 edges validated by the trifocal tensor
          m_mutexSet.discard(std::make_pair(std::min(I,J), std::max(I,J)));
          m_mutexSet.discard(std::make_pair(std::min(I,K), std::max(I,K)));
          m_mutexSet.discard(std::make_pair(std::min(J,K), std::max(J,K)));
          break;
        }
      }
    }
  }
}

} // namespace openMVG

#endif // OPENMVG_GLOBAL_SFM_ENGINE_RIG_TIJ_COMPUTATION_H
