// Copyright 2017 Google Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


/******************************************************************************
File:        InsetTilespaceIndex.cpp
Description:

Changes:


******************************************************************************/
#include "fusion/autoingest/sysman/InsetTilespaceIndex.h"
#include "common/khException.h"
#include <boost/range/sub_range.hpp>
#include <boost/range/as_literal.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/any_range.hpp>
#include <map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <stdint.h>
#include <assert.h>


QuadtreePath InsetTilespaceIndex::add(const khExtents <double> &extents) {
    int level;
    QuadtreePath quadtreeMbr = getQuadtreeMBR(extents, level, MAX_LEVEL);

    std::vector<khExtents <double> > *mbrExtentsVec;
    std::map<QuadtreePath, std::unordered_set<khExtents <double>>>::iterator quadTreeMbrExtentsVecPair;
    quadTreeMbrExtentsVecPair = _mbrExtentsSetMap.find(quadtreeMbr);

    if (quadTreeMbrExtentsVecPair == _mbrExtentsSetMap.end()) {
        mbrExtentsVec = new std::unordered_set<khExtents <double> >();
        _mbrExtentsSetMap.insert({quadtreeMbr, *mbrExtentsVec});
    } else {
        mbrExtentsVec = &(quadTreeMbrExtentsVecPair->second);
    }
    mbrExtentsVec->insert(extents);
    return quadtreeMbr;
}

QuadtreePath InsetTilespaceIndex::getQuadtreeMBR(const khExtents<double> &extents, int &level, const int max_level) {
    double north = 180;
    double south = -180;
    double west = -180;
    std::string base_path = "";
    char next_qt_node;
    for (level = 0; level < max_level; level += 1) {
        double level_dim_size = pow(2, level);
        double qt_node_size = 180.0 / level_dim_size;
        double north_south_midpoint = (south + north) / 2.0;
        // Get which sub-node the SW vertex is in.
        if (extents.south() <= north_south_midpoint) {
            if (extents.west() <= west + qt_node_size) {
                next_qt_node = '0';
            } else {
                next_qt_node = '1';
            }
        } else {
            if (extents.west() <= west + qt_node_size) {
                next_qt_node = '3';
            } else {
                next_qt_node = '2';
            }
        }
        // Check if NE vertex is in the same sub-node. If
        // not, then break at the node we are at.
        if (extents.north() <= north_south_midpoint) {
            if (extents.east() <= west + qt_node_size) {
                if (next_qt_node != '0') {
                    break;
                }
            } else {
                if (next_qt_node != '1') {
                    break;
                }
                west += qt_node_size;
            }
            north = north_south_midpoint;
        } else {
            if (extents.east() <= west + qt_node_size) {
                if (next_qt_node != '3') {
                    break;
                }
            } else {
                if (next_qt_node != '2') {
                    break;
                }
                west += qt_node_size;
            }
            south = north_south_midpoint;
        }
        // If still contained, descend to the next level of the tree.
        (base_path) += next_qt_node;
    }

    return QuadtreePath(base_path);

}

std::vector <QuadtreePath>
InsetTilespaceIndex::intersectingExtentsQuadtreePaths(QuadtreePath quadtreeMbr, uint32 minLevel, uint32 maxLevel) {
    std::vector <QuadtreePath> mbrHashVec;
    notify(NFY_DEBUG, "intersectingExtentsQuadtreePaths: index has quadtree count of: %lu", _mbrExtentsSetMap.size());
    boost::copy(_mbrExtentsSetMap | boost::adaptors::map_keys,
                std::back_inserter(mbrHashVec));
    notify(NFY_DEBUG, "intersectingExtentsQuadtreePaths: mbrHashVec size is: %lu ", mbrHashVec.size() );                
    std::vector <QuadtreePath> intersectingQuadtrees;

    // TODO - redo this section to use bitwise filtering and partitioning 
    // using the QuadtreePath's internal path_ bits, as this will be most 
    // expeditious.  However, this requires  access to private constructors 
    // and data. BTree lookups in the mbrHashVec could also bring the time 
    // complexity to O(log n)
    notify(NFY_DEBUG, "intersectingExtentsQuadtreePaths: checking levels %d to %d ", minLevel, maxLevel );
    for (uint32 level = minLevel; level <= maxLevel; level++) {
        for (QuadtreePath &otherMbr : mbrHashVec) {
            notify(NFY_DEBUG, "intersectingExtentsQuadtreePaths: Comparing extents from \n\tquery Quadtree:  %s with \n\tOTHER extents Quadtree:  %s",
                quadtreeMbr.AsString().c_str(), otherMbr.AsString().c_str() );
            if (otherMbr.Level() >= minLevel && otherMbr.Level() <= maxLevel) {
                if (QuadtreePath::OverlapsAtLevel(quadtreeMbr, otherMbr, level)) {

                    // TODO - there is probably an issue here with duplicate 
                    // "otherMbr"'s being added to the index - this needs
                    // to be be fixed in order to have the same results as 
                    // the original algo. 
                    intersectingQuadtrees.emplace_back(otherMbr);
                    notify(NFY_DEBUG, "\tOVERLAP at level %d", level );
                    break;
                } else {
                    notify(NFY_DEBUG, "\tno overlap at level %d", level );
                }
            }
            else { 
                notify(NFY_DEBUG, "\tno overlap - OTHER MBR is out of level bounds (%d not in %d->%d)", level, minLevel, maxLevel );
            }
        }
    }
    notify(NFY_DEBUG, "Found %lu intersecting Quadtree indexes with QUERY Quadtree:  %s ",
        intersectingQuadtrees.size(), quadtreeMbr.AsString().c_str() );
    return intersectingQuadtrees;
}


std::vector<khExtents <double> >

InsetTilespaceIndex::intersectingExtents(const QuadtreePath quadtreeMbr, uint32 minLevel, uint32 maxLevel) {
    std::vector <QuadtreePath> intersectingQuadtreeMbrs = intersectingExtentsQuadtreePaths(quadtreeMbr, minLevel,
                                                                                           maxLevel);
    notify(NFY_DEBUG, "intersectingExtents: Converting %lu intersecting Quadtree indexes to extents (MATCHING Quadtrees  %s ",
        intersectingQuadtreeMbrs.size(), quadtreeMbr.AsString().c_str() );
    std::vector<khExtents <double> > intersectingExtentsVec;
    for (auto otherMbr : intersectingQuadtreeMbrs) {
        std::vector<khExtents <double> > extentsVec = _mbrExtentsSetMap[otherMbr];
        notify(NFY_DEBUG, "\tintersectingExtents: Adding %lu intersecting extents from\n\t QUERY Quadtree:  %s OVERLAPs with  \n\tMATCHING extents Quadtree:  %s",
            extentsVec.size(), quadtreeMbr.AsString().c_str(), otherMbr.AsString().c_str() );
        intersectingExtentsVec.insert(intersectingExtentsVec.end(), extentsVec.begin(), extentsVec.end());
    };
    notify(NFY_DEBUG, "intersectingExtents: Found %lu intersecting Extents indexes with QUERY Quadtree:  %s ",
        intersectingExtentsVec.size(), quadtreeMbr.AsString().c_str() );
    
    return intersectingExtentsVec;
}





