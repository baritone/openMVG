
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// contributor : Stephane FLOTRON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/features/features.hpp"
#include "openMVG/image/image.hpp"
#include "openMVG/matching/indMatch_utils.hpp"
#include "openMVG/matching/indexed_sort.hpp"
#include "openMVG/multiview/triangulation_nview.hpp"

#include "software/globalSfM/indexedImageGraph.hpp"
#include "software/globalSfM/indexedImageGraphExport.hpp"
#include "software/globalSfM/SfMRigidGlobalEngine.hpp"
#include "software/globalSfM/SfMGlobal_tij_computation.hpp"
#include "software/SfM/SfMIOHelper.hpp"
#include "software/SfM/SfMRobust.hpp"
#include "software/SfM/SfMPlyHelper.hpp"

// Rotation averaging
#include "openMVG/multiview/rotation_averaging.hpp"
// Translation averaging
#include "openMVG/linearProgramming/lInfinityCV/global_translations_fromTriplets.hpp"
#include "openMVG/multiview/translation_averaging_solver.hpp"

// Linear programming solver(s)
#include "openMVG/linearProgramming/linearProgrammingInterface.hpp"
#include "openMVG/linearProgramming/linearProgrammingOSI_X.hpp"
#ifdef OPENMVG_HAVE_MOSEK
#include "openMVG/linearProgramming/linearProgrammingMOSEK.hpp"
#endif

#include "openMVG/graph/connectedComponent.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/vectorGraphics/svgDrawer.hpp"
#include "third_party/stlAddition/stlMap.hpp"
#include "third_party/histogram/histogram.hpp"

#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

#undef DYNAMIC
#include "openMVG/bundle_adjustment/problem_data_container.hpp"
#include "openMVG/bundle_adjustment/pinhole_ceres_functor.hpp"
#include "software/globalSfM/SfMBundleAdjustmentHelper_tonly.hpp"

#include "lemon/list_graph.h"
#include <lemon/connectivity.h>

#include "third_party/progress/progress.hpp"
#include "openMVG/system/timer.hpp"

#include <numeric>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <sstream>

#include <stdio.h>
#include <stdlib.h>

using namespace openMVG;
using namespace openMVG::graphUtils;
using namespace svg;

namespace openMVG{

} // namespace openMVG
