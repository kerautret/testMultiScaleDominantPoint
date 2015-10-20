#ifndef TESTFUNCTION
#define TESTFUNCTION

#include <iostream>
#include <exception>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/boards/Board2D.h"
#include "myfunctions.h"

using namespace std;
using namespace DGtal;
using namespace Z2i;

#define IS_PLOT true

/*** Burred segments decomposition ***/
vector<AlphaThickSegmentComputer2D> blurredSegmentCurveDecompositionV1(const vector<Point>& aContour, double thickness, const char* filename = NULL, bool verbose = false); /* consider consecutives points */
vector<AlphaThickSegmentComputer2D> blurredSegmentCurveDecompositionV2(const vector<Point>& aContour, double thickness, const char* filename = NULL, bool verbose = false); /* jump over points */
/*** Burred segments decomposition ***/

/*** Curvature estimation ****/
vector<double> curvatureOnCurveV1(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, int type = 0, const char* filename = NULL, bool verbose = false); /* non symmetry */
vector<double> curvatureOnCurveV2(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, int type = 0, const char* filename = NULL, bool verbose = false); /* symmetry */
vector<double> curvatureOnCurveV3(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, int type = 0, const char* filename = NULL, bool verbose = false); /* wrt direct neighbours */
vector<double> curvatureOnCurveCombine(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, int type = 0, const char* filename = NULL, bool verbose = false); /* means of 3 curvatures */
/*** Curvature estimation ****/

/**** Count segments passing by each point of the curves ****/
vector<Point> testCountSegmentsOnShape(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, const char* filename = NULL, bool verbose = false);
/**** Count segments passing by each point of the curves ****/

/**** Dominant points detections ****/
vector<Point> testDominantPointOnShapeV1(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, bool isClosed, const char* filename = NULL, bool verbose = false); /* Phuong's version q = p - 1 */
vector<Point> testDominantPointOnShapeV1(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, const vector<double>& curvature, bool isClosed, const char* filename = NULL, bool verbose = false); /* q = p - 1 + curvature */

vector<Point> testDominantPointOnShapeV2(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, bool isClosed, const char* filename = NULL, bool verbose = false); /* go as far as possible in the monotone sequence */
vector<Point> testDominantPointOnShapeV2(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, const vector<double>& curvature, bool isClosed, const char* filename = NULL, bool verbose = false); /* go as far as possible in the monotone sequence */

vector<Point> testDominantPointOnShapeV3(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, bool isClosed, const char* filename = NULL, bool verbose = false); /* Phuong's version q = p - 1 + with ordering ANGLE */
vector<Point> testDominantPointOnShapeV3(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, const vector<double>& curvature, bool isClosed, const char* filename = NULL, bool verbose = false); /* q = p - 1 + curvature + with ordering ANGLE  */

vector<Point> testDominantPointOnShapeV4(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, bool isClosed, const char* filename = NULL, bool verbose = false); /* Phuong's version q = p - 1 + without ordering ANGLE (only common zone) */
vector<Point> testDominantPointOnShapeV4(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, bool isSymmetry, bool isClosed, const char* filename = NULL, bool verbose = false); /* q = p - 1 + curvature + without ordering ANGLE (only common zone) */
/**** Dominant points detections ****/

/**** Dominant points selection ****/
vector<Point> testDominantPointSelectionV1(const vector<Point>& DP, const vector<int>& indexDP, const vector<double>& curvature, int nbDP, const vector<Point>& aContour, const char* filename = NULL, bool verbose = false); /* by curvature */
vector<Point> testDominantPointSelectionV2(const vector<Point>& DP, const vector<int>& indexDP, int nbDP, const vector<Point>& aContour, bool isClosed, const char* filename = NULL, bool verbose = false); /* by ise */
vector<Point> testDominantPointSelectionV2(const vector<Point>& DP, const vector<int>& indexDP, const vector<double>& curvature, int nbDP, const vector<Point>& aContour, bool isClosed, const char* filename = NULL, bool verbose = false); /* by ise*curvature */
vector<Point> testDominantPointSelectionV3(const vector<Point>& DP, const vector<int>& indexDP, int nbDP, const vector<Point>& aContour, bool isClosed, const char* filename = NULL, bool verbose = false); /* by ise*angle */

vector<Point> testDominantPointSelectionV4(const vector<Point>& DP, const vector<int>& indexDP, const vector<Point>& aContour, bool isClosed, const char* filename = NULL, bool verbose = false); /* by max Error crierion FOM1/FOM2/FOM3/FOM_ND */
/**** Dominant points selection ****/

/**** Tangent space transformation ****/
vector<PointD> tangentspaceV1(const vector<Point>& DP, bool normalized, double* totalLength, const char* filename = NULL, bool verbose = false);
vector<PointD> tangentspaceV2(const vector<Point>& DP, bool normalized, vector<double>& alpha, vector<double>& length, double* totalLength, const char* filename = NULL, bool verbose = false);
vector<PointD> tangentspaceV3(const vector<Point>& DP, bool normalized, vector<PointD>& Ti, vector<double>& alpha, vector<double>& length, double* totalLength, const char* filename = NULL, bool verbose = false);
/**** Tangent space transformation ****/

/****Verification of islolated points ***/
vector<PointD> testIsolatedPointsV1(const vector<PointD>& MP, double alphaMax, bool verbose); /* test isolated point w.r.t. alphaMax */
set<int> testIsolatedPointsV2(const vector<double>& alpha, const vector<double>& length, double alphaMax, double maxLength, bool verbose); /* test isolated point w.r.t. alphaMax + lengthMax */
/****Verification of islolated points ***/

/****Decomposition of Curve into Segments and Arcs ***/
vector<int> testDecompositionSegmentCircleV1(const vector<int>& indexDP, const vector<PointD>& MP, double alphaMax, vector<Point>& segments, vector<Point>& arcs, bool verbose = false); /* Phuong's version with alpha_max */
vector<int> testDecompositionSegmentCircleV1(const vector<int>& indexDP, const vector<PointD>& MP, double alphaMax, double thickness, vector<Point>& segments, vector<Point>& arcs, const char* filename = NULL, bool verbose = false); /* Phuong's version with alpha_max + seg decomposition in TS > 3 */
vector<int> testDecompositionSegmentCircleV2(const vector<Point>& aContour, const vector<int>& indexDP, const vector<PointD>& MP, double alphaMax, double thickness, vector<Point>& segments, vector<Point>& arcs, const char* filename = NULL, bool verbose = false); /* Phuong's version with alpha_max + seg decomposition in TS > 2 + ise_Arc>2*ise_Seg */
vector<int> testDecompositionSegmentCircleV3(const vector<Point>& aContour, const vector<int>& indexDP, const vector<PointD>& MP, double alphaMax, double thickness, vector<Point>& segments, vector<Point>& arcs, int nbPointCir = 3, double iseTol = 2.0, double angleTol = 10.0, const char* filename = NULL, bool verbose = false); /* Phuong's version with alpha_max + seg decomposition in TS >= nbPointCir + ise_Arc>iseTol*ise_Seg */
vector<int> testDecompositionSegmentCircleV3(const vector<Point>& aContour, const vector<int>& indexDP, const vector<PointD>& MP, double alphaMax, double thickness, vector<Point>& segments, vector<Point>& arcs, vector<double>& slope, int nbPointCir = 3, double iseTol = 2.0, double angleTol = 10.0, const char* filename = NULL, bool verbose = false);

void testDecompositionSegmentCircleV4(const vector<Point>& aContour, const vector<int>& indexDP, const vector<PointD>& MP, double thickness, vector<Point>& segments, vector<Point>& arcs, int nbPointCir = 3, double iseTol = 2.0, const char* filename = NULL, bool verbose = false); /* Phuong's version with alpha_max + seg decomposition in TS >= nbPointCir + ise_Arc>iseTol*ise_Seg */

void testDecompositionSegmentCircleV2(const vector<int>& indexDP, const vector<double>& alpha, const vector<double>& length, double alphaMax, double maxLength, vector<Point>& segments, vector<Point>& arcs, bool verbose = false); /* Phuong's version with alpha_max + length_max */
void testDecompositionSegmentCircleV2(const vector<int>& indexDP, const vector<PointD>& MP, const vector<double>& alpha, const vector<double>& length, double alphaMax, double maxLength, double thickness, vector<Point>& segments, vector<Point>& arcs, const char* filename = NULL, bool verbose = false); /* Phuong's version with alpha_max + length_max + seg decomposition in TS */

//new stragegy : arc/circle => monotone step fct
vector<int> testDecompositionSegmentCircleV5(const vector<Point>& aContour, const vector<int>& indexDP, const vector<PointD>& MP, double alphaMax, double thickness, vector<Point>& segments, vector<Point>& arcs, const char* filename = NULL, bool verbose = false);
/****Decomposition of Curve into Segments and Arcs ***/

/**** Draw decomposition from seg of seg and arcs ****/
void drawDecompositionSegmentCircleV1(const vector<Point>& aContour, const vector<Point>& DP, const vector<int> isolatedPoint, const vector<Point>& segments, const vector<Point>& arcs, char* filename); // min ise (= findBestFittingCircle) + min/max radii
void drawDecompositionSegmentCircleV2(const vector<Point>& aContour, const vector<Point>& DP, const vector<int> isolatedPoint, const vector<Point>& segments, const vector<Point>& arcs, const vector<double> radiusTangentSpace, double thickness, double thickness_MP, char* filename); // 2 end points + R in tangent space
void drawDecompositionSegmentCircleV3(const vector<Point>& aContour, const vector<Point>& DP, const vector<int> isolatedPoint, const vector<Point>& segments, const vector<Point>& arcs, double thickness, char* filename); // min ise +/- thickness*sqrt(2)

void drawDecompositionSegmentCircleV4(const vector<Point>& aContour, const vector<Point>& DP, const vector<int> isolatedPoint, const vector<Point>& segments, const vector<Point>& arcs, char* filename); // min/max radius (= findBestChord) + max/min radius

void drawDecompositionSegmentCircleAll(const vector<Point>& aContour, const vector<Point>& DP, const vector<int> isolatedPoint, const vector<Point>& segments, const vector<Point>& arcs, double thickness, double thickness_MP, char* filename); // all = findBestFittingCircle + 2EP and R + sqrt(2)Thickness
/**** Draw decomposition from seg of seg and arcs ****/

/**** Verifying circularity by chord prop ****/
void testVerifyCircularityChrodProp(const vector<Point>& aContour, const vector<Point>& DP, const vector<int>& indexDP, double thickness, vector<Point>& segments, vector<Point>& arcs, bool verbose = false);
/**** Verifying circularity by chord prop ****/

/********** Test circle simple with DP only *************/
int findBestFittingCircle(const vector<Point> aContour, int idStart, int idEnd);
double findMaxRaidusCircle(Point center, double radius, const vector<Point> aContour, int idStart, int idEnd);
double findMinRaidusCircle(Point center, double radius, const vector<Point> aContour, int idStart, int idEnd);

Point findBestChrodCircle(const vector<Point> aContour, int idBegin, int idEnd);
int findBestCircle(const vector<Point> aContour, int idStart, int idEnd);

void testCircleSimpleV1(const vector<Point>& aContour, const vector<Point>& DP, const char* filename = NULL, bool verbose = false);//draw approximated circle
void testCircleSimpleV2(const vector<Point>& aContour, const vector<Point>& DP, const char* filename = NULL, bool verbose = false);//with errors ise/lmax
void testCircleSimpleV3(const vector<Point>& aContour, const vector<Point>& DP, const vector<int>& indexDP, const char* filename = NULL, bool verbose = false);//with approximate radius
void testCircleSimpleV4(const vector<Point>& aContour, const vector<Point>& DP, const vector<int>& indexDP, const char* filename = NULL, bool verbose = false);//with approximate between seg and arc/circle
/********** Test circle simple with DP only *************/
#endif // TESTFUNCTION

