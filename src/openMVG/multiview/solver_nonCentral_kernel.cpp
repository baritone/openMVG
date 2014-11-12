
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

#include "openMVG/multiview/solver_nonCentral_kernel.hpp"
#include "openMVG/numeric/poly.h"

namespace openMVG {
namespace noncentral {
namespace kernel {

using namespace std;
using namespace opengv;

void SixPointSolver::Solve(adapter_t & adapter,
                  model_t & relativePose,
                  const std::vector<int> &indices)
{
   // create non central relative sac problem
   NoncentralRelativePoseSacProblem
              problem(adapter,
                      sac_problems::relative_pose::NoncentralRelativePoseSacProblem::SIXPT,
                      false);

   // solve pose problem
   problem.computeModelCoefficients(std::vector<int> &indices, relativePose);
}

}  // namespace kernel
}  // namespace noncentral
}  // namespace openMVG
