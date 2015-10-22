#include <iostream>
#include "myfunctions.h"
#include "testfunctions.h"

#include "DGtal/geometry/curves/FreemanChain.h"

#define INPUT "DataSet/ContoursTest"
#define OUTPUT "Result/ContoursTest"

//#define INPUT "DataSet/PGMSDP"
//#define OUTPUT "Result/PGMSDP"

using namespace std;

int main()
{
    cout << "Hello World!" << endl;

    /********** read data ***************/
    char fileContour[FILENAMESIZE];
    char file[] = "polygoneBruit2";//flowerNoisePartBis;flowerNoisePartBis;polygoneBruitPart;ContoursTest: polygoneBruit;polygoneBruit2;ellipse; ellipseBruit;flowerSmall2;pentagon50;cercleB2;formeTest6;flower-noise
    sprintf(fileContour,"%s/%s.sdp",INPUT,file);
    vector<Point> aContourTmp = readFile(fileContour,1);
    sprintf(fileContour,"%s/%s.dat",INPUT,file);
    //char file[] = "s01n003";
    if (!std::ifstream(fileContour))
        writeFile(aContourTmp, fileContour, true);
    vector<Point> aContour = readFile(fileContour,true);
    cout<<fileContour<<endl;
    /********** read data ***************/

    /********** convert to freeman code ***************/
    char fileFcContour[FILENAMESIZE];
    sprintf(fileFcContour,"%s/%s.fc",INPUT,file);
    if (!std::ifstream(fileFcContour))
    {
        FreemanChain<Z2i::Integer> fcContour (aContour);
        std::ofstream out(fileFcContour);
        fcContour.selfDisplay(out);
        out<<endl;
        out.close();
    }
    /********** convert to freeman code ***************/

    /********** call multiscale prog ***************/
    char instruction[5*FILENAMESIZE];

    /* meaningful scale *
    char noiseLevelMSFile[FILENAMESIZE];
    sprintf(noiseLevelMSFile,"%s/%sMeanScale.txt",OUTPUT,file);
    //std::system("./../../meaningfulscaleDemo/build/demoIPOL/meaningfulScaleEstim < ../../meaningfulscaleDemo/demoIPOL/Contours/ellipseBruit.fc -printNoiseLevel >  noiseLevel.txt");
    sprintf(instruction,"./../../../ImaGene/ImaGene-forIPOL/build/demoIPOL/meaningfulScaleEstim < %s -printNoiseLevel > %s",fileFcContour,noiseLevelMSFile);
    sprintf(instruction,"./../ImaGene-forIPOL/build/demoIPOL/meaningfulScaleEstim < %s -printNoiseLevel > %s",fileFcContour,noiseLevelMSFile);
    if (!std::ifstream(noiseLevelMSFile))
        std::system(instruction);
    //read file
    vector<double> vecMS = readMeanindfulThicknessFile(noiseLevelMSFile);//*sqrt(2)
    double glNoise = getGlobalNoise(vecMS);
    cout<<"File : "<<file<<" => globalNoise (meaningful thickness) = "<<glNoise<<endl;
    //colorate points by its meaningful thickness
    char fileColorMS[FILENAMESIZE];
    sprintf(fileColorMS,"%s/%s_Color.svg",OUTPUT,file);
    drawMeaningfulValue(aContour,vecMS,fileColorMS);
    * meaningful scale */

    /* meaningful thickness */
    char noiseLevelMTFile[FILENAMESIZE];
    sprintf(noiseLevelMTFile,"%s/%sMeanThickness.txt",OUTPUT,file);
    char smoothContourFile[FILENAMESIZE];
    sprintf(smoothContourFile,"%s/%sSmoothContour.txt",OUTPUT,file);
    ///sprintf(instruction,"./../../../ImaGene/ImaGene-forIPOL/build/tests/TestCompNoiseDetect/displayNoiseBS -srcPolygon %s 0 1 CLOSED -exportNoiseLevel %s -displaySmoothContour %s 0",fileContour,noiseLevelMTFile,smoothContourFile);
    sprintf(instruction,"./../ImaGene-forIPOL/build/tests/TestCompNoiseDetect/displayNoiseBS -srcPolygon %s 0 1 CLOSED -exportNoiseLevel %s -displaySmoothContour %s 0",fileContour,noiseLevelMTFile,smoothContourFile);
    if (!std::ifstream(noiseLevelMTFile))
        std::system(instruction);
    //read file
    vector<double> vecMT = readMeanindfulThicknessFile(noiseLevelMTFile);//*sqrt(2)
    double glNoise = getGlobalNoise(vecMT);
    cout<<"File : "<<file<<" => globalNoise (meaningful thickness) = "<<glNoise<<endl;
    //colorate points by its meaningful thickness
    char fileColorMT[FILENAMESIZE];
    sprintf(fileColorMT,"%s/%s_Color.svg",OUTPUT,file);
    drawMeaningfulValue(aContour,vecMT,fileColorMT);
    /* meaningful thickness */

    if(isRegularNoise(vecMT))
        cout<<"==> Regular noise "<<endl;
    else
        cout<<"==> Not regular noise "<<endl;
    /********** call multiscale prog ***************/

    /******** Adaptive tangent cover ***********/
    vector<AlphaThickSegmentComputer2D> adaptiveTangentCover;
    char fileAdaptMT[FILENAMESIZE];
    sprintf(fileAdaptMT,"%s/%s",OUTPUT,file);
    bool verbose = false;

    adaptiveTangentCover = testAdaptiveTangentCover(aContour,vecMT,fileAdaptMT,verbose);
    /******** Adaptive tangent cover ***********/

    /******** Dominant point detection with the adaptive tangent cover cover ******/
    vector<Point> DP;
    char filenameDP[FILENAMESIZE];
    sprintf(filenameDP,"%s/%s_DP.svg",OUTPUT,file); //DP = Dominant points
    bool isClosed = false;
    bool isSymmetry = false;

    cout<<"======= ATC ========";
    DP = testDominantPointOnShape(adaptiveTangentCover,aContour,isSymmetry,isClosed,filenameDP,verbose);
    cout<<endl<<"===> Num of dominant points is "<<DP.size()<<endl;

    //Find index of DP on the curve
    vector<int> indexDP;
    int indexLastDP=0;
    for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
    {
        indexLastDP = findElement(aContour,*it,indexLastDP);
        indexDP.push_back(indexLastDP);
    }

    if(verbose)
        cout<<"File: "<<fileContour<<endl;
    cout<<"Errors: "; error_All(aContour,DP,indexDP,isClosed);
    /******** Dominant point detection with the adaptive tangent cover cover ******/

    /********** Selection of dominant points ***************/
    char filenameDPnew[FILENAMESIZE];
    sprintf(filenameDPnew,"%s/%s_DPnew.svg",OUTPUT,file);
    vector<Point> newDP = testDominantPointSelectionV2(DP,indexDP,aContour,isClosed,filenameDPnew,verbose); // ISE * ANGLE
    cout<<"===> New num of dominant points is "<<newDP.size()<<endl;

    vector<int> indexNewDP;
    int indexLastNewDP=0;
    for(vector<Point>::const_iterator it = newDP.begin(); it != newDP.end(); it++)
    {
        indexLastNewDP = findElement(aContour,*it,indexLastNewDP);
        indexNewDP.push_back(indexLastNewDP);
    }

    if(verbose)
        cout<<"File: "<<fileContour<<endl;
    cout<<"Errors: "; error_All(aContour,newDP,indexNewDP,isClosed);
    /********** Selection of dominant points ***************/

    cout<<"======= Mid ========";
    double thickness = glNoise;
    vector<AlphaThickSegmentComputer2D> fuzzySegmentSetMid = blurredSegmentCurveDecomposition(aContour,thickness,NULL,false);
    vector<Point> DPM;
    //cout<<"===> Num of seg decomposed is "<<fuzzySegmentSetMid.size()<<endl;
    sprintf(filenameDP,"%s/%s_DPM.svg",OUTPUT,file); //DP = Dominant points
    DPM = testDominantPointOnShape(fuzzySegmentSetMid,aContour,isSymmetry,isClosed,filenameDP,verbose);
    cout<<endl<<"===> Num of dominant points is "<<DPM.size()<<endl;

    //Find index of DP on the curve
    vector<int> indexDPM;
    int indexLastDPM=0;
    for(vector<Point>::const_iterator it = DPM.begin(); it != DPM.end(); it++)
    {
        indexLastDPM = findElement(aContour,*it,indexLastDPM);
        indexDPM.push_back(indexLastDPM);
    }

    if(verbose)
        cout<<"File: "<<fileContour<<endl;
    cout<<"Errors: "; error_All(aContour,DPM,indexDPM,isClosed);
    /******** Dominant point detection with the adaptive tangent cover cover ******/

    /********** Selection of dominant points ***************/
    sprintf(filenameDPnew,"%s/%s_DPMnew_V.svg",OUTPUT,file);
    vector<Point> newDPM = testDominantPointSelectionV2(DPM,indexDPM,aContour,isClosed,filenameDPnew,verbose); // ISE * ANGLE

    cout<<"===> New num of dominant points is "<<newDPM.size()<<endl;
    vector<int> indexNewDPM;
    int indexLastNewDPM=0;
    for(vector<Point>::const_iterator it = newDPM.begin(); it != newDPM.end(); it++)
    {
        indexLastNewDPM = findElement(aContour,*it,indexLastNewDPM);
        indexNewDPM.push_back(indexLastNewDPM);
    }

    if(verbose)
        cout<<"File: "<<fileContour<<endl;
    cout<<"Errors: "; error_All(aContour,newDPM,indexNewDPM,isClosed);
    /********** Selection of dominant points ***************/

    return 0;
}

