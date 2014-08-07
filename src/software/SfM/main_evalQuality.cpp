
// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "software/SfMViewer/document.h"
#include "software/SfM/io_readGT.hpp"
#include "software/SfM/tools_precisionEvaluationToGt.hpp"
#include "software/SfM/SfMPlyHelper.hpp"

#include "third_party/htmlDoc/htmlDoc.hpp"

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <cstdlib>

using namespace openMVG;

int main(int argc, char **argv)
{
  using namespace std;

  CmdLine cmd;

  std::string
    sGTDirectory,
    sComputedDirectory,
    sOutDir = "";

  cmd.add( make_option('i', sGTDirectory, "gt") );
  cmd.add( make_option('c', sComputedDirectory, "computed") );
  cmd.add( make_option('o', sOutDir, "outdir") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
      << "[-i|--gt path (where ground truth camera trajectory are saved)] \n"
      << "[-c|--computed path (openMVG SfM_Output directory)] \n"
      << "[-o|--output path (where statistics will be saved)] \n"
      << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create(sOutDir);

  //---------------------------------------
  // Quality evaluation
  //---------------------------------------

  // Load GT:

  std::map< std::string, std::pair<Mat3, Vec3> > map_Rt_gt;
  std::map< size_t, PinholeCamera> map_Cam_gt;
  // READ DATA FROM GT
  {
    std::cout << "\nTry to read data from GT";
    std::vector<std::string> vec_fileNames;
    std::string sGTPath = sGTDirectory;
    readGt( sGTPath, vec_fileNames, map_Rt_gt, map_Cam_gt);
    std::cout << map_Cam_gt.size() << " gt cameras have been found" << std::endl;
  }

  //-- Load the camera of the computed directory
  Document doc;
  if (!doc.load(sComputedDirectory))
  {
    std::cerr << "\nCannot read input calibration." << std::endl;
  }

  // Prepare data for comparison (corresponding camera centers and rotations)
  const std::map<size_t, PinholeCamera> & _map_camera = doc._map_camera;
  std::vector<Vec3> vec_camPosGT, vec_C;
  std::vector<Mat3> vec_camRotGT, vec_camRot;
  for(std::map< size_t, PinholeCamera>::const_iterator iterGT = map_Cam_gt.begin();
    iterGT != map_Cam_gt.end(); ++iterGT)
  {
    size_t index = iterGT->first;
    vec_camPosGT.push_back(iterGT->second._C);
    vec_camRotGT.push_back(iterGT->second._R);

    //-- Computed
    std::map<size_t, PinholeCamera >::const_iterator iterComputed = _map_camera.find(index);
    vec_C.push_back(iterComputed->second._C);
    vec_camRot.push_back(iterComputed->second._R);
  }

  // Visual output of the camera location
  plyHelper::exportToPly(vec_camPosGT, string(stlplus::folder_append_separator(sOutDir) + "camGT.ply").c_str());
  plyHelper::exportToPly(vec_C, string(stlplus::folder_append_separator(sOutDir) + "camComputed.ply").c_str());

  // Evaluation
  htmlDocument::htmlDocumentStream _htmlDocStream("openMVG Quality evaluation.");
  EvaluteToGT(vec_camPosGT, vec_C, vec_camRotGT, vec_camRot, sOutDir, &_htmlDocStream);

  ofstream htmlFileStream( string(stlplus::folder_append_separator(sOutDir) +
    "ExternalCalib_Report.html").c_str());
  htmlFileStream << _htmlDocStream.getDoc();

  return EXIT_SUCCESS;
}

