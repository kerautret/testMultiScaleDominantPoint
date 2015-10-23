#ifndef TESTFUNCTION
#define TESTFUNCTION

#include <iostream>
#include <exception>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/boards/Board2D.h"
#include "myfunctions.h"
#include "testfunctions.h"

using namespace std;
using namespace DGtal;
using namespace Z2i;

#define IS_PLOT true

#define FILENAMESIZE 100

/*** Burred segments decomposition ***/
vector<AlphaThickSegmentComputer2D> blurredSegmentCurveDecomposition(const vector<Point>& aContour, double thickness, const char* filename = NULL, bool verbose = false); /* jump over points */
/*** Burred segments decomposition ***/

/**** Dominant points detections ****/
vector<Point> testDominantPointOnShape(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, bool isSymmetry, bool isClosed, const char* filename = NULL, bool verbose = false);
/**** Dominant points detections ****/

/**** Dominant points selection ****/
vector<Point> testDominantPointSelectionV1(const vector<Point>& DP, const vector<int>& indexDP, int nbDP, const vector<Point>& aContour, bool isClosed, const char* filename = NULL, bool verbose = false); /* by ise*angle */
vector<Point> testDominantPointSelectionV2(const vector<Point>& DP, const vector<int>& indexDP, const vector<Point>& aContour, bool isClosed, const char* filename = NULL, bool verbose = false); /* by max Error crierion FOM1/FOM2/FOM3/FOM_ND */
/**** Dominant points selection ****/

/********** Draw Multi-Thickness Cover **********/
void drawMultiThicknessCover(const vector<Point>& aContour, const vector<double>& thckVect, int nbColor, char* filename);
void drawMultiThicknessCover(const vector<Point>& aContour, const vector<vector<AlphaThickSegmentComputer2D>>& meaningThicknessTangentCover, char* filename);
void drawMultiThicknessCover(const vector<Point>& aContour, const vector<AlphaThickSegmentComputer2D>& meaningThicknessTangentCover, char* filename);
void drawMultiThicknessCover(const vector<Point>& aContour, const vector<vector<AlphaThickSegmentComputer2D>>& meaningThicknessTangentCover, const vector<double>& thckVect, char* filename);
void drawMeaningfulValue(const vector<Point>& aContour, const vector<double> vecMeanVal, const char* filename);
/********** Draw Multi-Thickness Cover **********/

/******** Adaptive tangent cover computation ***********/
vector<AlphaThickSegmentComputer2D> testAdaptiveTangentCover(const vector<Point>& aContour, const vector<double>& vecMT, const char* filename = NULL, bool verbose = false);
/******** Adaptive tangent cover computation ***********/
#endif // TESTFUNCTION

