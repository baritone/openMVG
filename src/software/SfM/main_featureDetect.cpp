
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/image/image.hpp"
#include "openMVG/features/features.hpp"

#include "software/SfM/SfMIOHelper.hpp"

/// Feature detector and descriptor interface
#include "nonFree/sift/SIFT.hpp"

#include "openMVG/system/timer.hpp"

#include "third_party/cmdLine/cmdLine.h"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"
#include "third_party/progress/progress.hpp"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>

using namespace openMVG;
using namespace std;

int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sImaDirectory;
  std::string sOutDir = "";
  bool bOctMinus1 = false;
  float dPeakThreshold = 0.04f;

  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('s', bOctMinus1, "octminus1") );
  cmd.add( make_option('p', dPeakThreshold, "peakThreshold") );

  try {
      if (argc == 1) throw std::string("Invalid command line parameter.");
      cmd.process(argc, argv);
  } catch(const std::string& s) {
      std::cerr << "\nUsage: " << argv[0] << " [OPTIONS...] -o <outdir> <image> ...\n\n"
      << "OPTIONS: \n"
      << "-o|--outdir <path>           Destination directory\n"
      << "-r|--distratio <ratio>       Distance ratio (default=0.6)\n"
      << "-s|--octminus1 [0|1]         When set to 1, SIFT will upscale the image x2\n"
      << "-p|--peakThreshold <thresh>  Peak threshold for SIFT, 0.04 (default) to 0.01\n"
      << std::endl;

      std::cerr << s << std::endl;
      return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
            << argv[0] << std::endl
            << "--outdir " << sOutDir << std::endl
            << "--octminus1 " << bOctMinus1 << std::endl
            << "--peakThreshold " << dPeakThreshold << std::endl;

  if (sOutDir.empty())  {
    std::cerr << "\nIt is an invalid output directory" << std::endl;
    return EXIT_FAILURE;
  }

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create( sOutDir );

  //---------------------------------------
  // Compute features and descriptor
  //    - extract sift features and descriptor
  //    - save features and descriptors on disk
  //---------------------------------------

  typedef Descriptor<unsigned char, 128> DescriptorT;
  typedef SIOPointFeature FeatureT;
  typedef std::vector<FeatureT> FeatsT;
  typedef vector<DescriptorT > DescsT;
  typedef KeypointSet<FeatsT, DescsT > KeypointSetT;

  {
    Timer timer;
    std::cout << "\n\n - EXTRACT FEATURES - " << std::endl;

    Image<unsigned char> imageGray;

    C_Progress_display my_progress_bar( argc );
    for(size_t i=1; i < argc; ++i)  {

      std::string sFeat = stlplus::create_filespec(sOutDir,
        stlplus::basename_part(argv[i]), "feat");
      std::string sDesc = stlplus::create_filespec(sOutDir,
        stlplus::basename_part(argv[i]), "desc");

      //If descriptors or features file are missing, compute them
      if (!stlplus::file_exists(sFeat) || !stlplus::file_exists(sDesc)) {

        if (!ReadImage(argv[i], &imageGray))
          continue;

        // Compute features and descriptors and export them to files
        KeypointSetT kpSet;
        SIFTDetector(imageGray,
          kpSet.features(), kpSet.descriptors(),
          bOctMinus1, true, dPeakThreshold);
        kpSet.saveToBinFile(sFeat, sDesc);
      }
      ++my_progress_bar;
    }
    ++my_progress_bar;
    std::cout << "\nTask done in (s): " << timer.elapsed() << std::endl;
  }
  return EXIT_SUCCESS;
}
