// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/solver_nonCentral_kernel.hpp"
#include "openMVG/numeric/poly.h"
#include <set>

namespace openMVG {
namespace noncentral {
namespace kernel {

using namespace std;
using namespace opengv;

void SixPointSolver::Solve(relative_pose::NoncentralRelativeAdapter & adapter,
                  std::vector<transformation_t> * models,
                  const std::vector<size_t> &indices)
{

  // convert size_t to int for opengv call
  std::vector<int> idx;

  for(size_t i=0; i < indices.size(); ++i)
  {
     idx.push_back( (int) indices[i]);
  }

   // create non central relative sac problem
   sac_problems::relative_pose::NoncentralRelativePoseSacProblem
              problem(adapter,
                      sac_problems::relative_pose::NoncentralRelativePoseSacProblem::SIXPT,
                      false);

   // solve pose problem
   transformation_t relativePose;
   problem.computeModelCoefficients(idx, relativePose);

   models->push_back(relativePose);
}

void GePointSolver::Solve(relative_pose::NoncentralRelativeAdapter & adapter,
std::vector<transformation_t> * models,
const std::vector<size_t> &indices)
{

  // convert size_t to int for opengv call
  std::vector<int> idx;

  for(size_t i=0; i < indices.size(); ++i)
  {
    idx.push_back( (int) indices[i]);
  }

  // create non central relative sac problem
  sac_problems::relative_pose::NoncentralRelativePoseSacProblem
  problem(adapter,
  sac_problems::relative_pose::NoncentralRelativePoseSacProblem::GE,
  false);

  // solve pose problem
  transformation_t relativePose;
  problem.computeModelCoefficients(idx, relativePose);

  models->push_back(relativePose);
}

}  // namespace kernel
}  // namespace noncentral
}  // namespace openMVG
