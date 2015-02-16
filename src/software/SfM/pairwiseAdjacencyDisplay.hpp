
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_PAIRWISE_ADJACENCY_DISPLAY_H
#define OPENMVG_PAIRWISE_ADJACENCY_DISPLAY_H

#include "third_party/vectorGraphics/svgDrawer.hpp"
using namespace svg;

#include "openMVG/matching/indMatch.hpp"
using namespace openMVG::matching;

namespace openMVG  {

  /// Display pair wises matches as an Adjacency matrix in svg format
  void RigWiseMatchingToAdjacencyMatrixSVG(const size_t NbImages,
  const matching::RigWiseMatches & map_GeometricMatches,
  const std::string & sOutName)
{
  if ( !map_GeometricMatches.empty())
  {
    float scaleFactor = 5.0f;
    svgDrawer svgStream((NbImages+3)*5, (NbImages+3)*5);
    // Parse filtered list
    for ( size_t i = 0 ; i < map_GeometricMatches.size(); ++i)
    {
      RigWiseMatches::const_iterator iterMap = map_GeometricMatches.begin();
      advance(iterMap, i);

      const PairWiseMatches  map_Matches = iterMap->second;

      // loop on intrarig matches
      for ( size_t j = 0 ; j < map_Matches.size(); ++j)
      {
        PairWiseMatches::const_iterator iter = map_Matches.begin();
        advance(iter, j);

        // extract image indexes
        const size_t I = iter->first.first;
        const size_t J = iter->first.second;

        // If the pair have matches display a blue boxes at I,J position.
        if (!iter->second.empty())
        {
          // Display as a tooltip: (IndexI, IndexJ NbMatches)
          std::ostringstream os;
          os << "(" << J << "," << I << " " << iter->second.size() <<")";
          svgStream.drawSquare(J*scaleFactor, I*scaleFactor, scaleFactor/2.0f,
          svgStyle().fill("blue").noStroke());
        } // HINT : THINK ABOUT OPACITY [0.4 -> 1.0] TO EXPRESS MATCH COUNT
      }
    }
    // Display axes with 0 -> NbImages annotation : _|
    ostringstream osNbImages;   osNbImages << NbImages;
    svgStream.drawText((NbImages+1)*scaleFactor, scaleFactor, scaleFactor, "0", "black");
    svgStream.drawText((NbImages+1)*scaleFactor,
    (NbImages)*scaleFactor - scaleFactor, scaleFactor, osNbImages.str(), "black");
    svgStream.drawLine((NbImages+1)*scaleFactor, 2*scaleFactor,
    (NbImages+1)*scaleFactor, (NbImages)*scaleFactor - 2*scaleFactor,
    svgStyle().stroke("black", 1.0));

    svgStream.drawText(scaleFactor, (NbImages+1)*scaleFactor, scaleFactor, "0", "black");
    svgStream.drawText((NbImages)*scaleFactor - scaleFactor,
    (NbImages+1)*scaleFactor, scaleFactor, osNbImages.str(), "black");
    svgStream.drawLine(2*scaleFactor, (NbImages+1)*scaleFactor,
    (NbImages)*scaleFactor - 2*scaleFactor, (NbImages+1)*scaleFactor,
    svgStyle().stroke("black", 1.0));

    ofstream svgFileStream( sOutName.c_str());
    svgFileStream << svgStream.closeSvgFile().str();
  }
}


/// Display pair wises matches as an Adjacency matrix in svg format
void PairWiseMatchingToAdjacencyMatrixSVG(const size_t NbImages,
  const matching::PairWiseMatches & map_Matches,
  const std::string & sOutName)
{
  if ( !map_Matches.empty())
  {
    float scaleFactor = 5.0f;
    svgDrawer svgStream((NbImages+3)*5, (NbImages+3)*5);
    // Go along all possible pair
    for (size_t I = 0; I < NbImages; ++I) {
      for (size_t J = 0; J < NbImages; ++J) {
        // If the pair have matches display a blue boxes at I,J position.
        matching::PairWiseMatches::const_iterator iterSearch =
          map_Matches.find(make_pair(I,J));
        if (iterSearch != map_Matches.end() && !iterSearch->second.empty())
        {
          // Display as a tooltip: (IndexI, IndexJ NbMatches)
          std::ostringstream os;
          os << "(" << J << "," << I << " " << iterSearch->second.size() <<")";
          svgStream.drawSquare(J*scaleFactor, I*scaleFactor, scaleFactor/2.0f,
            svgStyle().fill("blue").noStroke());
        } // HINT : THINK ABOUT OPACITY [0.4 -> 1.0] TO EXPRESS MATCH COUNT
      }
    }
    // Display axes with 0 -> NbImages annotation : _|
    ostringstream osNbImages;   osNbImages << NbImages;
    svgStream.drawText((NbImages+1)*scaleFactor, scaleFactor, scaleFactor, "0", "black");
    svgStream.drawText((NbImages+1)*scaleFactor,
      (NbImages)*scaleFactor - scaleFactor, scaleFactor, osNbImages.str(), "black");
    svgStream.drawLine((NbImages+1)*scaleFactor, 2*scaleFactor,
      (NbImages+1)*scaleFactor, (NbImages)*scaleFactor - 2*scaleFactor,
      svgStyle().stroke("black", 1.0));

    svgStream.drawText(scaleFactor, (NbImages+1)*scaleFactor, scaleFactor, "0", "black");
    svgStream.drawText((NbImages)*scaleFactor - scaleFactor,
      (NbImages+1)*scaleFactor, scaleFactor, osNbImages.str(), "black");
    svgStream.drawLine(2*scaleFactor, (NbImages+1)*scaleFactor,
      (NbImages)*scaleFactor - 2*scaleFactor, (NbImages+1)*scaleFactor,
      svgStyle().stroke("black", 1.0));

    ofstream svgFileStream( sOutName.c_str());
    svgFileStream << svgStream.closeSvgFile().str();
  }
}

} // namespace openMVG

#endif // OPENMVG_PAIRWISE_ADJACENCY_DISPLAY_H
