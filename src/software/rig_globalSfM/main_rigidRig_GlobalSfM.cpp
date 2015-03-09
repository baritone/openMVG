
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cstdlib>

#include "software/rig_globalSfM/SfMRigidGlobalEngine.hpp"
using namespace openMVG;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "openMVG/system/timer.hpp"

int main(int argc, char **argv)
{
  using namespace std;
  std::cout << std::endl
    << "-----------------------------------------------------------\n"
    << "Global Structure from Motion:\n"
    << "-----------------------------------------------------------\n"
    << "Open Source implementation of the paper:\n"
    << "\"Global Fusion of Relative Motions for "
    << "Robust, Accurate and Scalable Structure from Motion.\"\n"
    << "Pierre Moulon, Pascal Monasse and Renaud Marlet. "
    << " ICCV 2013." << std::endl
    << "------------------------------------------------------------"
    << std::endl;


  CmdLine cmd;

  std::string sImaDirectory;
  std::string sMatchesDir;
  std::string sOutDir = "";
  bool bColoredPointCloud = false;
  int iRotationAveragingMethod = 2;
  int iTranslationAveragingMethod = 1;
  bool bRefineIntrinsic = false;
  bool bRefineRigStruct = false;

  cmd.add( make_option('i', sImaDirectory, "imadir") );
  cmd.add( make_option('m', sMatchesDir, "matchdir") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('c', bColoredPointCloud, "coloredPointCloud") );
  cmd.add( make_option('r', iRotationAveragingMethod, "rotationAveraging") );
  cmd.add( make_option('t', iTranslationAveragingMethod, "translationAveraging") );
  cmd.add( make_option('f', bRefineIntrinsic, "refineIntrinsic") );
  cmd.add( make_option('e', bRefineRigStruct, "refineRigStructure") );

  try {
    if (argc == 1) throw std::string("Invalid parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--imadir path]\n"
    << "[-m|--matchdir path]\n"
    << "[-o|--outdir path]\n"
    << "[-c|--coloredPointCloud 0(default) or 1]\n"
    << "[-r|--rotationAveraging 2(default L2) or 1 (L1)]\n"
    << "[-t|--translationAveraging 1(default L1) or 2 (L2)]\n"
    << "[-f|--refineRigIntrinsic \n"
    << "\t 0-> keep provided intrinsic,\n"
    << "\t 1-> refine provided intrinsics: (focal, principal point) ] \n"
    << "[-e|--refineRigStructure \n"
    << "\t 0-> keep provided rig rotation and translation,\n"
    << "\t 1-> refine provided rig : (R, t) ] \n"
    << "\n"
    << " ICCV 2013: => -r 2 -t 1"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if (iRotationAveragingMethod < ROTATION_AVERAGING_L1 ||
      iRotationAveragingMethod > ROTATION_AVERAGING_L2 )  {
    std::cerr << "\n Rotation averaging method is invalid" << std::endl;
    return EXIT_FAILURE;
  }

    if (iTranslationAveragingMethod < TRANSLATION_AVERAGING_L1 ||
      iTranslationAveragingMethod > TRANSLATION_AVERAGING_L2 )  {
    std::cerr << "\n Translation averaging method is invalid" << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create(sOutDir);

  //---------------------------------------
  // Incremental reconstruction process
  //---------------------------------------

  openMVG::Timer timer;
  GlobalRigidReconstructionEngine to3DEngine(
    sImaDirectory,
    sMatchesDir,
    sOutDir,
    ERotationAveragingMethod(iRotationAveragingMethod),
    ETranslationAveragingMethod(iTranslationAveragingMethod),
    true);

  to3DEngine.setRefineIntrinsics(bRefineIntrinsic);
  to3DEngine.setRefineRigStruct (bRefineRigStruct);

  if (to3DEngine.Process())
  {
    std::cout << std::endl << " Total Ac-Global-Sfm took (s): " << timer.elapsed() << std::endl;

    //-- Compute color if requested
    const rigReconstructorHelper & reconstructorHelperRef = to3DEngine.refToRigReconstructorHelper();
    std::vector<Vec3> vec_tracksColor;
    if (bColoredPointCloud)
    {
      // Compute the color of each track
      to3DEngine.ColorizeTracks(to3DEngine.getTracks(), vec_tracksColor);
    }

    //-- Export computed data to disk
    reconstructorHelperRef.exportToPly(
      stlplus::create_filespec(sOutDir, "FinalColorized", ".ply"),
      bColoredPointCloud ? &vec_tracksColor : NULL);

    // Export to openMVG format
    std::cout << std::endl << "Export 3D scene to openMVG format" << std::endl
      << " -- Point cloud color: " << (bColoredPointCloud ? "ON" : "OFF") << std::endl;

    if (!reconstructorHelperRef.ExportToOpenMVGFormat(
      stlplus::folder_append_separator(sOutDir) + "SfM_output",
      to3DEngine.getFilenamesVector(),
      sImaDirectory,
      to3DEngine.getImagesSize(),
      to3DEngine.getTracks(),
      bColoredPointCloud ? &vec_tracksColor : NULL,
      false,
      std::string("generated by the Global OpenMVG Calibration Engine")))
    {
      std::cerr << "Error while saving the scene." << std::endl;
    }


    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}
