
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_IO_H
#define OPENMVG_SFM_IO_H

#include "openMVG/numeric/numeric.h"
#include "openMVG/split/split.hpp"

#include <fstream>
#include <iterator>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>    // std::lexicographical_compare
#include <cctype>       // std::tolower

namespace openMVG{
namespace SfMIO{

//  basic camera structure
struct CameraInfo
{
  std::string m_sImageName;
  size_t m_intrinsicId;
};

struct IntrinsicCameraInfo
{
  size_t m_w, m_h;
  float m_focal;
  Mat3 m_K;
  double m_k1;
  double m_k2;
  double m_k3;
  bool m_bKnownIntrinsic; // true if 11 or 6, else false
  std::string m_sCameraMaker, m_sCameraModel;

  IntrinsicCameraInfo(): m_w(0), m_h(0), m_K(Mat3::Zero()), m_k1(0), m_k2(0), m_k3(0), m_bKnownIntrinsic(false), m_sCameraModel(""), m_sCameraMaker("")
  {  }

  /// Functor used to tell if two IntrinsicCameraInfo share the same optical properties
  friend bool operator== (IntrinsicCameraInfo const &ci1, IntrinsicCameraInfo const &ci2)
  {
    // Two camera share optical properties if they share the same K matrix (and the same camera name)
    bool bequal = ci1.m_K == ci2.m_K && ci1.m_sCameraMaker == ci2.m_sCameraMaker && ci1.m_sCameraModel == ci2.m_sCameraModel;
    return bequal;
  }
};

// rigid rig camera structure
struct CameraRigInfo
{
  std::string m_sImageName;
  std::string m_sRigName;
  size_t m_intrinsicId;
  size_t m_rigId;
  size_t m_subCameraId;
};

struct IntrinsicCameraRigInfo
{
  size_t m_w, m_h;
  float m_focal;
  Mat3 m_K;
  Mat3 m_R;
  Vec3 m_rigC;
  bool m_bKnownIntrinsic; // true if 11 or 6, else false
  std::string m_sCameraMaker, m_sCameraModel;

  IntrinsicCameraRigInfo(): m_w(0), m_h(0), m_K(Mat3::Zero()), m_bKnownIntrinsic(false), m_sCameraModel(""), m_sCameraMaker("")
      , m_R(Mat3::Zero()), m_rigC(Vec3::Zero())
  {  }

  /// Functor used to tell if two IntrinsicCameraInfo share the same optical properties
  friend bool operator== (IntrinsicCameraRigInfo const &ci1, IntrinsicCameraRigInfo const &ci2)
  {
    // Two camera share optical properties if they share the same K matrix (and the same camera name)
    bool bequal =  ci1.m_K == ci2.m_K && ci1.m_sCameraMaker == ci2.m_sCameraMaker && ci1.m_sCameraModel == ci2.m_sCameraModel
                && ci1.m_R == ci2.m_R && ci1.m_rigC == ci2.m_rigC ;
    return bequal;
  }
};

// Load an image file list
// One basename per line.
// It handle different scenario based on the intrinsic info of the tested image
// - a camera without exif data
// - a camera with exif data found in the database
// - a camera with exif data not found in the database
// - a camera with known intrinsic
static bool loadImageList( std::vector<CameraInfo> & vec_camImageName,
                           std::vector<IntrinsicCameraInfo> & vec_focalGroup,
                           const std::string & sFileName,
                           bool bVerbose = true )
{
  std::ifstream in(sFileName.c_str());
  if(!in.is_open())  {
    std::cerr << std::endl
      << "Impossible to read the specified file." << std::endl;
  }
  std::string sValue;
  std::vector<std::string> vec_str;
  while(getline( in, sValue ) )
  {
    vec_str.clear();
    split( sValue, ";", vec_str );
    if (vec_str.size() == 1)
    {
      std::cerr << "Invalid input file" << std::endl;
      in.close();
      return false;
    }
    std::stringstream oss;
    oss.clear(); oss.str(vec_str[1]);
    size_t width, height;
    oss >> width;
    oss.clear(); oss.str(vec_str[2]);
    oss >> height;

    IntrinsicCameraInfo intrinsicCamInfo;
    intrinsicCamInfo.m_w = width;
    intrinsicCamInfo.m_h = height;

    switch ( vec_str.size() )
    {
      case 3 : // a camera without exif data
      {
         intrinsicCamInfo.m_focal = -1;
         intrinsicCamInfo.m_bKnownIntrinsic = false;
         intrinsicCamInfo.m_sCameraMaker = "";
         intrinsicCamInfo.m_sCameraModel = "";
      }
      break;
      case 5 : // a camera with exif data not found in the database
      {
         intrinsicCamInfo.m_focal = -1;
         intrinsicCamInfo.m_bKnownIntrinsic = false;
         intrinsicCamInfo.m_sCameraMaker = vec_str[3];
         intrinsicCamInfo.m_sCameraModel = vec_str[4];
      }
      break;
      case  6 : // a camera with exif data found in the database
      {
         oss.clear(); oss.str(vec_str[3]);
         float focal;
         oss >> focal;
         intrinsicCamInfo.m_focal = focal;
         intrinsicCamInfo.m_bKnownIntrinsic = true;
         intrinsicCamInfo.m_sCameraMaker = vec_str[4];
         intrinsicCamInfo.m_sCameraModel = vec_str[5];

         Mat3 K;
         K << focal, 0, float(width) / 2.f,
              0, focal, float(height) / 2.f,
              0, 0, 1;
         intrinsicCamInfo.m_K = K;
      }
      break;
      case 12 : // a camera with known intrinsic
      {
        intrinsicCamInfo.m_bKnownIntrinsic = true;
        intrinsicCamInfo.m_sCameraMaker = intrinsicCamInfo.m_sCameraModel = "";

        Mat3 K = Mat3::Identity();

        oss.clear(); oss.str(vec_str[3]);
        oss >> K(0,0);
        oss.clear(); oss.str(vec_str[4]);
        oss >> K(0,1);
        oss.clear(); oss.str(vec_str[5]);
        oss >> K(0,2);
        oss.clear(); oss.str(vec_str[6]);
        oss >> K(1,0);
        oss.clear(); oss.str(vec_str[7]);
        oss >> K(1,1);
        oss.clear(); oss.str(vec_str[8]);
        oss >> K(1,2);
        oss.clear(); oss.str(vec_str[9]);
        oss >> K(2,0);
        oss.clear(); oss.str(vec_str[10]);
        oss >> K(2,1);
        oss.clear(); oss.str(vec_str[11]);
        oss >> K(2,2);

        intrinsicCamInfo.m_K = K;
        intrinsicCamInfo.m_focal = static_cast<float>(K(0,0)); // unknown sensor size;
      }
      break;
      case 26 : // a camera with known intrinsic and rig
      {
        intrinsicCamInfo.m_bKnownIntrinsic = true;
        intrinsicCamInfo.m_sCameraMaker = intrinsicCamInfo.m_sCameraModel = "";

        Mat3 K = Mat3::Identity();

        oss.clear(); oss.str(vec_str[3]);
        oss >> K(0,0);
        oss.clear(); oss.str(vec_str[4]);
        oss >> K(0,1);
        oss.clear(); oss.str(vec_str[5]);
        oss >> K(0,2);
        oss.clear(); oss.str(vec_str[6]);
        oss >> K(1,0);
        oss.clear(); oss.str(vec_str[7]);
        oss >> K(1,1);
        oss.clear(); oss.str(vec_str[8]);
        oss >> K(1,2);
        oss.clear(); oss.str(vec_str[9]);
        oss >> K(2,0);
        oss.clear(); oss.str(vec_str[10]);
        oss >> K(2,1);
        oss.clear(); oss.str(vec_str[11]);
        oss >> K(2,2);

        intrinsicCamInfo.m_K = K;
        intrinsicCamInfo.m_focal = static_cast<float>(K(0,0)); // unknown sensor size;
      }
      break;
      default :
      {
        std::cerr << "Invalid image list line: wrong number of arguments" << std::endl;
        in.close();
        return false;
      }
    }

    std::vector<IntrinsicCameraInfo>::const_iterator iterIntrinsicGroup = find(vec_focalGroup.begin(), vec_focalGroup.end(), intrinsicCamInfo);
    size_t id = -1;
    if ( iterIntrinsicGroup == vec_focalGroup.end())
    {
      vec_focalGroup.push_back(intrinsicCamInfo);
      id = vec_focalGroup.size()-1;
    }
    else
    {
      id = std::distance( std::vector<IntrinsicCameraInfo>::const_iterator(vec_focalGroup.begin()), iterIntrinsicGroup);
    }

    CameraInfo camInfo;
    camInfo.m_sImageName = vec_str[0];
    camInfo.m_intrinsicId = id;
    vec_camImageName.push_back(camInfo);

    vec_str.clear();
  }
  in.close();
  return !(vec_camImageName.empty());
}

//-- Load an image list file but only return camera image names
static bool loadImageList( std::vector<std::string> & vec_camImageName,
                           const std::string & sListFileName,
                           bool bVerbose = true )
{
  vec_camImageName.clear();
  std::vector<openMVG::SfMIO::CameraInfo> vec_camImageIntrinsicInfo;
  std::vector<openMVG::SfMIO::IntrinsicCameraInfo> vec_focalGroup;
  if (loadImageList( vec_camImageIntrinsicInfo,
                      vec_focalGroup,
                      sListFileName,
                      bVerbose) )
  {
    for ( std::vector<openMVG::SfMIO::CameraInfo>::const_iterator
      iter_camInfo = vec_camImageIntrinsicInfo.begin();
      iter_camInfo != vec_camImageIntrinsicInfo.end();
      iter_camInfo++ )
    {
      vec_camImageName.push_back( iter_camInfo->m_sImageName );
    }
  }
  return (!vec_camImageName.empty());
}

// Load an image file list
// One basename per line.
// It handle different scenario based on the intrinsic info of the tested image
// - a camera without exif data
// - a camera with exif data found in the database
// - a camera with exif data not found in the database
// - a camera with known intrinsic
static bool loadImageList(
         std::vector<CameraRigInfo> & vec_camImageName,
         std::vector<IntrinsicCameraRigInfo> & vec_focalGroup,
         const std::string & sFileName,
         bool bVerbose = true )
{
  std::ifstream in(sFileName.c_str());
  if(!in.is_open())  {
    std::cerr << std::endl
      << "Impossible to read the specified file." << std::endl;
  }
  std::string sValue;
  std::vector<std::string> vec_str;
  std::vector<std::string> vec_rigSet;
  while(getline( in, sValue ) )
  {
    vec_str.clear();
    split( sValue, ";", vec_str );
    if (vec_str.size() == 1)
    {
      std::cerr << "Invalid input file" << std::endl;
      in.close();
      return false;
    }
    std::stringstream oss;
    oss.clear(); oss.str(vec_str[1]);
    size_t width, height;
    oss >> width;
    oss.clear(); oss.str(vec_str[2]);
    oss >> height;

    IntrinsicCameraRigInfo intrinsicCamInfo;
    intrinsicCamInfo.m_w = width;
    intrinsicCamInfo.m_h = height;

    switch ( vec_str.size() )
    {
      case 26 : // a camera with known intrinsic and rig structure
      {
        intrinsicCamInfo.m_bKnownIntrinsic = true;
        intrinsicCamInfo.m_sCameraMaker = intrinsicCamInfo.m_sCameraModel = "";

        Mat3 K = Mat3::Identity();

        // get camera matrix
        oss.clear(); oss.str(vec_str[3]);
        oss >> K(0,0);
        oss.clear(); oss.str(vec_str[4]);
        oss >> K(0,1);
        oss.clear(); oss.str(vec_str[5]);
        oss >> K(0,2);
        oss.clear(); oss.str(vec_str[6]);
        oss >> K(1,0);
        oss.clear(); oss.str(vec_str[7]);
        oss >> K(1,1);
        oss.clear(); oss.str(vec_str[8]);
        oss >> K(1,2);
        oss.clear(); oss.str(vec_str[9]);
        oss >> K(2,0);
        oss.clear(); oss.str(vec_str[10]);
        oss >> K(2,1);
        oss.clear(); oss.str(vec_str[11]);
        oss >> K(2,2);

        intrinsicCamInfo.m_K = K;
        intrinsicCamInfo.m_focal = static_cast<float>(K(0,0)); // unknown sensor size;

        // get rotation
        Mat3 R = Mat3::Identity();
        oss.clear(); oss.str(vec_str[14]);
        oss >> R(0,0);
        oss.clear(); oss.str(vec_str[15]);
        oss >> R(0,1);
        oss.clear(); oss.str(vec_str[16]);
        oss >> R(0,2);
        oss.clear(); oss.str(vec_str[17]);
        oss >> R(1,0);
        oss.clear(); oss.str(vec_str[18]);
        oss >> R(1,1);
        oss.clear(); oss.str(vec_str[19]);
        oss >> R(1,2);
        oss.clear(); oss.str(vec_str[20]);
        oss >> R(2,0);
        oss.clear(); oss.str(vec_str[21]);
        oss >> R(2,1);
        oss.clear(); oss.str(vec_str[22]);
        oss >> R(2,2);

        intrinsicCamInfo.m_R = R;

        // get center
        Vec3  C = Vec3::Zero();
        oss.clear(); oss.str(vec_str[23]);
        oss >> C(0);
        oss.clear(); oss.str(vec_str[24]);
        oss >> C(1);
        oss.clear(); oss.str(vec_str[25]);
        oss >> C(2);

        intrinsicCamInfo.m_rigC=C;

      }
      break;
      default :
      {
        std::cerr << "Invalid image list line: wrong number of arguments" << std::endl;
        in.close();
        return false;
      }
    }

    // intrinsic group
    std::vector<IntrinsicCameraRigInfo>::const_iterator iterIntrinsicGroup = find(vec_focalGroup.begin(), vec_focalGroup.end(), intrinsicCamInfo);
    size_t id = -1;
    if ( iterIntrinsicGroup == vec_focalGroup.end())
    {
      vec_focalGroup.push_back(intrinsicCamInfo);
      id = vec_focalGroup.size()-1;
    }
    else
    {
      id = std::distance( std::vector<IntrinsicCameraRigInfo>::const_iterator(vec_focalGroup.begin()), iterIntrinsicGroup);
    }

    // rig Id
    const std::string rigName = vec_str[12];
    std::vector<std::string>::const_iterator iterRigSet = find(vec_rigSet.begin(), vec_rigSet.end(), rigName);
    size_t idr = -1;
    if ( iterRigSet == vec_rigSet.end())
    {
      vec_rigSet.push_back(rigName);
      idr = vec_rigSet.size()-1;
    }
    else
    {
      idr = std::distance( std::vector<std::string>::const_iterator(vec_rigSet.begin()), iterRigSet);
    }


    CameraRigInfo camInfo;
    camInfo.m_sImageName    = vec_str[0];
    camInfo.m_intrinsicId   = id;
    camInfo.m_sRigName      = vec_str[12];
    camInfo.m_rigId         = idr;
    camInfo.m_subCameraId   = atoi(vec_str[13].c_str());
    vec_camImageName.push_back(camInfo);

    vec_str.clear();
  }
  in.close();
  return !(vec_camImageName.empty());
}

} // namespace SfMIO
} // namespace openMVG

#endif // OPENMVG_SFM_IO_H
