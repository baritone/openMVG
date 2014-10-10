
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// contributor : Stephane FLOTRON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_RIGID_GLOBAL_SFM_ENGINE_H
#define OPENMVG_RIGID_GLOBAL_SFM_ENGINE_H

#include "openMVG/numeric/numeric.h"

#include "openMVG/cameras/PinholeCamera.hpp"
#include "software/SfM/SfMEngine.hpp"
#include "software/SfM/SfMIOHelper.hpp"
#include "software/SfM/SfMReconstructionData.hpp"
#include "software/globalSfM/SfMGlobalEngine.hpp"
class SIOPointFeature;

#include "openMVG/tracks/tracks.hpp"
#include "openMVG/matching/indMatch.hpp"
using namespace openMVG::tracks;

#include "openMVG/graph/triplet_finder.hpp"
#include "openMVG/multiview/translation_averaging_common.hpp"

#include "third_party/htmlDoc/htmlDoc.hpp"

#include <memory>

namespace openMVG{

//------------------
//-- Bibliography --
//------------------
//- [1] "Global Fusion of Relative Motions for Robust, Accurate and Scalable Structure from Motion."
//- Authors: Pierre MOULON, Pascal MONASSE and Renaud MARLET.
//- Date: December 2013.
//- Conference: ICCV.
//- [2] "Disambiguating visual relations using loop constraints."
//- Autors: Christopher ZACH, Manfred KLOPSCHITZ and Marc POLLEFEYS.
//- Date: 2010
//- Conference : CVPR.
//- [3] "Efficient and Robust Large-Scale Rotation Averaging"
//- Authors: Avishek Chatterjee and Venu Madhav Govindu
//- Date: December 2013.
//- Conference: ICCV.

// Implementation of [1]
// Some points differs from the [1] paper to ease open source port:
//-- Relative rotation inference:
//   - only the triplet rejection is performed (in [1] a Bayesian inference on cycle error is performed [2])
//-- Global rotation computation:
//    - a sparse least square,
//    - or, the L1 averaging method of [3].
//-- Linear Programming solver:
//   - in order to have the best performance it is advised to used the MOSEK LP backend.

class GlobalRigidReconstructionEngine : public GlobalReconstructionEngine {


};


} // namespace openMVG

#endif // OPENMVG_RIGID_GLOBAL_SFM_ENGINE_H
