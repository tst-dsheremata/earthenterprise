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

// Unit tests for khQtPreindex.

#include "common/khStringUtils.h"
#include "common/quadtreepath.h"
#include "fusion/autoingest/sysman/InsetTilespaceIndex.h"
//#include "fusion/autoingest/sysman/InsetInfo.h"
#include "common/khException.h"
#include "common/khExtents.h"
#include "common/khInsetCoverage.h"
#include "common/khTileAddr.h"
#include "autoingest/plugins/RasterProductAsset.h"
#include "autoingest/plugins/MercatorRasterProductAsset.h"
#include <UnitTest.h>
#include <gtest/gtest.h>
#include <boost/range/sub_range.hpp>
#include <boost/range/as_literal.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/any_range.hpp>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <stdint.h>
#include <assert.h>
#include <khTypes.h>
#include <sstream>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class TestData {
    std::vector <khExtents<double>> extentsVec;

public:
    TestData(std::string fname) {
        char path[354];
        getcwd(path,352) ;
        std::cout << "Current path is " << path << '\n';
        std::cout << "Reading date from file: " << fname;
        std::ifstream input(fname);
        if (input) {
            // First line and all even lines have 4 floats, for decimal degree 
            // representations of the input insets.
            // Every odd line is a QtMBR.  Let's just read them in pairs.
            for (std::string line; std::getline(input, line);)   
            {
                std::istringstream in(line);
                double x1, x2, y1, y2;
                in >> x1 >> x2 >> y1 >> y2;
                notify(NFY_WARN, "Test Data: %f\t%f\t%f\t%f", x1,x2,y1,y2);
                khExtents<double> extents(XYOrder, x1, x2, y1, y2);
                std::getline(input, line);
                // Note -  discarding the MBR for now, we'll recalculate it.
                std::string qtp_txt, qtp_level;
                in >> qtp_txt >> qtp_level;
                // QuadtreePath qtpMbr(qtp_txt);
                extentsVec.push_back(extents);
                notify(NFY_WARN, "Added a new Extents object. extentsVec size is now: %lu", extentsVec.size() );
            }
        }
        else {
            notify(NFY_WARN, "Not a file");
        }
    }

    //std::map<khExtents < double>, QuadtreePath>
    std::vector<khExtents<double>>
    getData() { return extentsVec; }
};

class InsetTilespaceIndexTest : public testing::Test {
protected:
    InsetTilespaceIndex insetTilespaceIndex;
public:

    std::vector <khExtents<double>>
    findInsetsControlAlgo(const khInsetCoverage &queryCoverageArea,
                          const std::vector <khExtents<double>> &testDegExtentsVec) {
        uint beginLevel = queryCoverageArea.beginLevel();
        uint endLevel = queryCoverageArea.endLevel();;
        khExtents<uint32> queryExtents = queryCoverageArea.levelExtents(beginLevel);
        std::vector <uint32> neededIndexes; //This is our return value...
        std::vector <khExtents<double>> matchingExtents;
        std::vector <khExtents<uint32>> tileExtentsVec;
        
        notify(NFY_WARN, "Querying for extents that intersect with query Coverage Area:  \n\t%d, %d, %d, %d ... ",
                queryExtents.beginX(), queryExtents.endX(), queryExtents.beginY(), queryExtents.endY() );

        for (auto degExtents : testDegExtentsVec) {
            khExtents <uint32> tileExtents = DegExtentsToTileExtents(degExtents, 21);
            tileExtentsVec.push_back(tileExtents);
            notify(NFY_WARN, "Converted test extents from decimal degrees \n\t%f, %f, %f, %f ... \nto tilespace: %d, %d, %d, %d \n\t",
                    degExtents.beginX(), degExtents.endX(), degExtents.beginY(), degExtents.endY(),
                    tileExtents.beginX(), tileExtents.endX(), tileExtents.beginY(), tileExtents.endY());
        }
        FindNeededImageryInsets(queryCoverageArea,
                                tileExtentsVec,
                                tileExtentsVec.size(),
                                neededIndexes,
                                beginLevel,
                                endLevel);
        notify(NFY_WARN, "%lu", neededIndexes.size());
        for (uint index : neededIndexes) {
            khExtents<double> ex = testDegExtentsVec[index];
            matchingExtents.push_back(ex);
        };
        return matchingExtents;
    }

    std::vector <khExtents<double>>
    findInsetsExperimentalAlgo(const khInsetCoverage &queryCoverageArea,
                               const std::vector <khExtents<double>> &testDegExtentsVec) {
        std::vector <khExtents<double>> matchingExtentsVec;
        InsetTilespaceIndex qtIndex;

        for (auto extents : testDegExtentsVec) {
            qtIndex.add(extents);
        }
        int level = queryCoverageArea.beginLevel();
        auto queryMBR = qtIndex.getQuadtreeMBR(queryCoverageArea.degreeExtents(), level, queryCoverageArea.endLevel());

        matchingExtentsVec = qtIndex.intersectingExtents(
                queryMBR,
                queryCoverageArea.beginLevel(),
                queryCoverageArea.endLevel());
        notify(NFY_WARN, "Got Here");

        return matchingExtentsVec;
    }


    const bool compareAlgorithmOutputs() {
        TestData dataset("./src/fusion/testdata/TinyTestDataQTPs.txt");
        //std::vector<khExtents<double>> testExtents;
        khExtents<double> queryExtentsArea(XYOrder, 114.032, 114.167, 19.1851, 19.3137);
        khInsetCoverage queryCoverageArea(RasterProductTilespaceFlat, queryExtentsArea, 19, 9, 21);
        std::vector<khExtents<double>> extentsTestDataVec = dataset.getData();
        std::vector<khExtents<double> > requiredExtentsProd = findInsetsControlAlgo(queryCoverageArea, extentsTestDataVec);
        notify(NFY_WARN, "Old Algo Done, %lu insets overlap" , requiredExtentsProd.size());
        std::vector<khExtents<double> > requiredExtentsExp = findInsetsExperimentalAlgo(queryCoverageArea, extentsTestDataVec);
        notify(NFY_WARN, "New Algo Done, %lu insets overlap", requiredExtentsExp.size());
        //std::sort(requiredExtentsProd.begin(), requiredExtentsExp.end()); 
        bool listsMatch = (requiredExtentsProd == requiredExtentsExp);
        return listsMatch;
    }

    InsetTilespaceIndexTest() : insetTilespaceIndex() {};
};

//This test should result in the getExtentMBR method breaking after level 0.

TEST_F(InsetTilespaceIndexTest, BreakBeforeMaxLevelReached
) {
khExtents<double> extent = khExtents<double>(XYOrder, 180, 180, 90, 180);
int level;
insetTilespaceIndex.
getQuadtreeMBR(extent, level, MAX_LEVEL
);
EXPECT_EQ(1, level);
}

//This test should result in the getExtentMBR method looping until the max_level is reached.
TEST_F(InsetTilespaceIndexTest, DecendToMaxLevel
) {
khExtents<double> extent = khExtents<double>(XYOrder, 90, 90, 90, 90);
int level;
insetTilespaceIndex.
getQuadtreeMBR(extent, level, MAX_LEVEL
);
EXPECT_EQ(MAX_LEVEL, level
);
}

TEST_F(InsetTilespaceIndexTest, compareAlgorithmOutputs) {
    EXPECT_EQ(compareAlgorithmOutputs(), true);
}

//This test should result in the quadtree path 202 being returned.
TEST_F(InsetTilespaceIndexTest, ReturnSpecificQTP
) {
khExtents<double> extent = khExtents<double>(XYOrder, 90, 90, 90, 90);
int level;
QuadtreePath qtp = insetTilespaceIndex.getQuadtreeMBR(extent, level, 3);
EXPECT_EQ("202", qtp.

AsString()

);
}

int main(int argc, char *argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

