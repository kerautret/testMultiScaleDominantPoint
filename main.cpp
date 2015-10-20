#include "testfunctions.h"

#define FILENAMESIZE 100

//#define INPUT "DataSet"
//#define OUTPUT "Result"

//#define INPUT "DataSet/PhuongDataSet"
//#define OUTPUT "Result/PhuongDataSet"

//#define INPUT "DataSet/AvecBruit"
//#define OUTPUT "Result/AvecBruit"

//#define INPUT "DataSet/SansBruit"
//#define OUTPUT "Result/SansBruit"

//#define INPUT "DataSet/PGMSDP"
//#define OUTPUT "Result/PGMSDP"

#define INPUT "DataSet/ContoursTest"
#define OUTPUT "Result/ContoursTest"

//#define INPUT "DataSet/Benmark2D"
//#define OUTPUT "Result/Benmark2D"

//#define INPUT "DataSet/Sketches"
//#define OUTPUT "Result/Sketches"

int main()
{
    /********** read data ***************/
    char fileContour[FILENAMESIZE];

    //char file[] = "courbe";
    //sprintf(fileContour,"%s/%s.txt",INPUT,file);

    //char file[] = "s01n003P3";//s01n003P3;s01n003P3b;courbeTest
    //sprintf(fileContour,"%s/%s.txt",INPUT,file);

    //char file[] = "exemple5";//PhuongDataSet : exemple3;exemple5;exemple7;exemple1-8;symbol004;symbol100;
    //sprintf(fileContour,"%s/%s.sdp",INPUT,file);

    //char file[] = "s01n003";//PGMSDP : s01n003-9;s02n001-11;s03n001-11;s04n001-11;s05n001-11;s06n001-11;s07n001-11;s08n001-11;s09n001-11
    //sprintf(fileContour,"%s/%s.sdp",INPUT,file);

    //char file[] = "flowerNoise";//Benmark2D : s01n003;avion;circle;Cercle20;birdAst;cle;Glas01b;heart;homme;lapin;main;poisson;ray17;courbe.txt
    //sprintf(fileContour,"%s/%s.dat",INPUT,file);//chromosome;leaf;semicircle

    //char file[] = "Part_of_Glas01b";
    //sprintf(fileContour,"%s/%s.txt",INPUT,file);
    //vector<Point> aContour = readFile(fileContour,false);
    //vector<Point> aContour = readFileInverse(fileContour,false);

    char file[] = "formeTest3";//ContoursTest: ellipse; ellipseBruit;flowerSmall2;pentagon50;cercleB2;formeTest6;flower-noise
    sprintf(fileContour,"%s/%s.sdp",INPUT,file);
    vector<Point> aContour = readFile(fileContour,1);

    //char file[] = "s16_01";//Sketches : s01_01; s17_48; s16_01
    //sprintf(fileContour,"%s/%s.txt",INPUT,file);
    //vector<Point> aContour = readFile(fileContour,true);

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
    /* meaningful thickness */
    /********** call multiscale prog ***************/

    double thickness = glNoise;//2.0;
    bool verbose = FALSE;
    bool isClosed = false;
    int version = 2;
    bool isSymmetry = false;
    if(verbose)
    {
        cout<<"File: "<<fileContour<<", Thck: "<<thickness<<", Ver: "<<version<<endl<<endl;
        cout<<"Contour size = "<<aContour.size()<<endl;
    }
    /********** read data ***************/

    if(isSymmetry)
        cout<<"Symmetry angle !"<<endl;
    else
        cout<<"Non symmetry angle !"<<endl;

    char coutFile[FILENAMESIZE];
    //sprintf(coutFile,"%s/%s_Log_T%d_V%d_S%d.txt",OUTPUT,file,int(thickness*10),version,isSymmetry); //BS = blurred segments
    sprintf(coutFile,"%s/%s_Log_T%d_S%d.txt",OUTPUT,file,int(thickness*10),isSymmetry); //BS = blurred segments
    if(verbose)
        freopen(coutFile,"w",stdout);

    cout<<"File: "<<fileContour<<", Thck: "<<thickness<<", Ver: "<<version<<endl<<endl;
    cout<<"Contour size = "<<aContour.size()<<endl;

    /********** calcul the burred-segments corver ***************/
    char filenameBS[FILENAMESIZE];
    sprintf(filenameBS,"%s/%s_BS_T%d.svg",OUTPUT,file,int(thickness*10)); //BS = blurred segments
    //vector<AlphaThickSegmentComputer2D> fuzzySegmentSet = blurredSegmentCurveDecompositionV1(aContour,thickness,filenameBS,verbose);
    vector<AlphaThickSegmentComputer2D> fuzzySegmentSet = blurredSegmentCurveDecompositionV2(aContour,thickness,filenameBS,verbose);
    cout<<"===> Num of seg decomposed is "<<fuzzySegmentSet.size()<<endl<<endl;
    /********** calcul the burred-segments corver ***************/








    /********** tangential cover and domiant point *******/
    char filenameTC[FILENAMESIZE];
    sprintf(filenameTC,"%s/%s_CSP_T%d.svg",OUTPUT,file,int(thickness*10)); //CSP = count segment point
    testCountSegmentsOnShape(fuzzySegmentSet,aContour,filenameTC,verbose);
    /********** tangential cover and domiant point *******/










    /********** calcul the curvature on each point of the contour ***************/
    int type = 1;
    char filenameCur[FILENAMESIZE];
    sprintf(filenameCur,"%s/%s_CP_T%d.svg",OUTPUT,file,int(thickness*10)); //CP = curvature point
    vector<double> curvature;
    if(!isSymmetry)
        curvature = curvatureOnCurveV1(fuzzySegmentSet,aContour,type,NULL,verbose); //non symmetry
    else
        curvature = curvatureOnCurveV2(fuzzySegmentSet,aContour,type,NULL,verbose); //with symmetry
    //vector<double> curvature = curvatureOnCurveCombine(fuzzySegmentSet,aContour,1,NULL,true);
    /********** calcul the curvature on each point of the contour ***************/













    /********** calcul the dominant points ***************/
    ///char filenameDP[FILENAMESIZE];
    ///sprintf(filenameDP,"%s/%s_DP_T%d_V%d_S%d.svg",OUTPUT,file,int(thickness*10),version,isSymmetry); //DP = Dominant points

    //vector<Point> DP = testDominantPointOnShapeV1(fuzzySegmentSet,aContour,filenameDP,verbose);//V1 : mid point
    //vector<Point> DP = testDominantPointOnShapeV1(fuzzySegmentSet,aContour,curvature,filenameDP,verbose);//V2 : curvature
    //vector<Point> DP = testDominantPointOnShapeV2(fuzzySegmentSet,aContour,curvature,filenameDP,verbose);//V3 : new strategy
    //DP.insert(DP.begin(),aContour.front());
    //DP.insert(DP.end(),aContour.back());

    vector<Point> DP;
    char filenameDP[FILENAMESIZE];
    sprintf(filenameDP,"%s/%s_DP_T%d_V%d.svg",OUTPUT,file,int(thickness*10),version); //DP = Dominant points
    writeFile(DP,filenameDP,false);

    /* Version with seq of slope ordering of burred segments*
    if(version==1)
        DP = testDominantPointOnShapeV1(fuzzySegmentSet,aContour,isClosed,filenameDP,verbose);//V1 : mid point
    else if(version==2)
        DP = testDominantPointOnShapeV1(fuzzySegmentSet,aContour,curvature,isClosed,filenameDP,verbose);//V2 : curvature
    * Version with seq of slope ordering */

    /* Version with q starts by last sequence monotone of burred segments *
    if(version==1)
        DP = testDominantPointOnShapeV2(fuzzySegmentSet,aContour,isClosed,filenameDP,verbose);//V1 : mid point
    else if(version==2)
        DP = testDominantPointOnShapeV2(fuzzySegmentSet,aContour,curvature,isClosed,filenameDP,verbose);//V2 : curvature
    * Version with q starts by last sequence monotone of burred segments */

    /* Version with seq of ANGLE ordering of burred segments */
    if(version==1)
        DP = testDominantPointOnShapeV3(fuzzySegmentSet,aContour,isClosed,filenameDP,verbose);//V1 : mid point
    else if(version==2)
        DP = testDominantPointOnShapeV3(fuzzySegmentSet,aContour,curvature,isClosed,filenameDP,verbose);//V2 : curvature
    /* Version with seq of ANGLE ordering */

    /* Version without seq of ANGLE ordering of burred segments *
    if(version==1)
        DP = testDominantPointOnShapeV4(fuzzySegmentSet,aContour,isClosed,filenameDP,verbose);//V1 : mid point
    else if(version==2)
        DP = testDominantPointOnShapeV4(fuzzySegmentSet,aContour,isSymmetry,isClosed,filenameDP,verbose);//V2 : curvature
    * Version without seq of ANGLE ordering */

    cout<<"===> Num of dominant points is "<<DP.size()<<endl<<endl;
    //Find index of DP on the curve
    vector<int> indexDP;
    for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
        indexDP.push_back(findElement(aContour,*it));

    if(verbose)
        cout<<"File: "<<fileContour<<", Thck: "<<thickness<<", Ver: "<<version<<endl;
    ///cout<<"Errors: "; error_All(aContour,DP,indexDP);
    cout<<"Errors: "; error_All(aContour,DP,indexDP,isClosed);
    /********** calcul the dominant points ***************/





    /********** Selection of dominant points ***************/
    int nbDP = DP.size();//int(0.7*DP.size());//DP.size()-4;//int(0.9*DP.size());///DP.size()-9
    char filenameDPnew[FILENAMESIZE];
    sprintf(filenameDPnew,"%s/%s_DPnew_T%d_V%d_N%d.svg",OUTPUT,file,int(thickness*10),version,nbDP);
    //vector<Point> newDP = testDominantPointSelectionV1(DP,indexDP,curvature,nbDP,aContour,filenameDPnew,verbose); // CURVATURE
    //vector<Point> newDP = testDominantPointSelectionV2(DP,indexDP,nbDP,aContour,isClosed,filenameDPnew,verbose); // ISE
    //vector<Point> newDP = testDominantPointSelectionV2(DP,indexDP,curvature,nbDP,aContour,isClosed,filenameDPnew,verbose); // ISE + CURVATURE
    //vector<Point> newDP = testDominantPointSelectionV3(DP,indexDP,nbDP,aContour,isClosed,filenameDPnew,verbose); // ISE * ANGLE
    vector<Point> newDP = testDominantPointSelectionV4(DP,indexDP,aContour,isClosed,filenameDPnew,verbose); // ISE * ANGLE
    //writeFile(newDP,"points.txt",false);

    cout<<"===> New num of dominant points is "<<newDP.size()<<endl;
    vector<int> indexNewDP;
    for(vector<Point>::const_iterator it = newDP.begin(); it != newDP.end(); it++)
        indexNewDP.push_back(findElement(aContour,*it));

    if(verbose)
        cout<<"File: "<<fileContour<<", Thck: "<<thickness<<", Ver: "<<version<<endl;
    ///cout<<"Errors: "; error_All(aContour,newDP,indexNewDP);
    cout<<"Errors: "; error_All(aContour,newDP,indexNewDP,isClosed);
    /********** Selection of dominant points ***************/










    /********* transform DP into the tangent space ********/
    char filenameTSN[FILENAMESIZE];
    sprintf(filenameTSN,"%s/%s_TSN_T%d_V%d.svg",OUTPUT,file,int(thickness*10),version);
    double totalLength=0;
    bool normalized = true;
    //vector<PointD> MP = tangentspaceV1(newDP,normalized,&totalLength,filenameTSN,verbose);
    tangentspaceV1(newDP,normalized,&totalLength,filenameTSN,verbose);

    char filenameTS[FILENAMESIZE];
    sprintf(filenameTS,"%s/%s_TS_T%d_V%d.svg",OUTPUT,file,int(thickness*10),version);
    vector<double> alpha;
    vector<double> length;
    vector<PointD> MP = tangentspaceV2(newDP,false,alpha,length,&totalLength,filenameTS,verbose);
    tangentspaceV2(newDP,false,alpha,length,&totalLength,filenameTS,verbose);
    cout<<"===> Num of middle points in the tangent space is "<<MP.size()<<endl;
    /********* transform DP into the tangent space ********/



    //return 0;





    /********* decomposition segment circle : Phuong version ***************/
    double alpha_max = M_PI/4.0;//(45.0*M_PI)/180.0;
    double thickness_MP = 0.2;//0.5
    vector<Point> DecSegmentsPh;
    vector<Point> DecArcsPh;
    ///testDecompositionSegmentCircleV1(indexDP,MP,alpha_max,DecSegmentsPh,DecArcsPh,false);
    ///cout<<"===> Num of segments (V1) is "<<DecSegmentsPh.size()<<" and of arcs is "<<DecArcsPh.size()<<endl;
    char filenameTSBS[FILENAMESIZE];
    sprintf(filenameTSBS,"%s/%s_TSBS_T%d_V%d.svg",OUTPUT,file,int(thickness*10),version);
    //vector<int> isolatedPoint = testDecompositionSegmentCircleV1(indexNewDP,MP,alpha_max,thickness_MP,DecSegmentsPh,DecArcsPh,filenameTSBS,verbose);
    //vector<int> isolatedPoint = testDecompositionSegmentCircleV2(aContour,indexNewDP,MP,alpha_max,thickness_MP,DecSegmentsPh,DecArcsPh,filenameTSBS,verbose);
    ///vector<int> isolatedPoint = testDecompositionSegmentCircleV3(aContour,indexNewDP,MP,alpha_max,thickness_MP,DecSegmentsPh,DecArcsPh,2,20.0,5.0,filenameTSBS,verbose);
    vector<double> radiusTangentSpace;
    vector<int> isolatedPoint = testDecompositionSegmentCircleV3(aContour,indexNewDP,MP,alpha_max,thickness_MP,DecSegmentsPh,DecArcsPh,radiusTangentSpace,2,3.0,10.0,filenameTSBS,verbose);
    //testDecompositionSegmentCircleV4(aContour,indexNewDP,MP,thickness_MP,DecSegmentsPh,DecArcsPh);
    cout<<"===> Num of segments (V1) is "<<DecSegmentsPh.size()<<" and of arcs is "<<DecArcsPh.size()<<endl;

    char filenameDecSAPh[FILENAMESIZE];
    sprintf(filenameDecSAPh,"%s/%s_DecSAV1_T%d_V%d.svg",OUTPUT,file,int(thickness*10),version);
    drawDecompositionSegmentCircleAll(aContour,newDP,isolatedPoint,DecSegmentsPh,DecArcsPh,thickness,thickness_MP,filenameDecSAPh);
    //drawDecompositionSegmentCircleV4(aContour,newDP,isolatedPoint,DecSegmentsPh,DecArcsPh,filenameDecSAPh);
    /********* decomposition segment circle : Phuong version ***************/












    /********* decomposition segment circle : Phuong version modified ***************
    double alpha_max = M_PI/4.0;
    double thickness_MP = 0.5;
    double shapeDetail = getShapeDetail(aContour); //size of totalLength vs. shape bounding box !
    double length_max = shapeDetail*totalLength/MP.size();//totalLength/10.0;//100000
    vector<Point> DecSegmentsPhM;
    vector<Point> DecArcsPhM;
    ///testDecompositionSegmentCircleV2(indexDP,alpha,length,alpha_max,length_max,DecSegmentsPhM,DecArcsPhM,false);
    ///cout<<"===> Num of segments (V2) is "<<DecSegmentsPhM.size()<<" and of arcs is "<<DecArcsPhM.size()<<endl;
    char filenameTSBSM[FILENAMESIZE];
    sprintf(filenameTSBSM,"%s/%s_TSBSM_T%d_V%d.svg",OUTPUT,file,int(thickness*10),version);
    testDecompositionSegmentCircleV2(indexDP,MP,alpha,length,alpha_max,length_max,thickness_MP,DecSegmentsPhM,DecArcsPhM,filenameTSBSM,true);
    cout<<"===> Num of segments (V2) is "<<DecSegmentsPhM.size()<<" and of arcs is "<<DecArcsPhM.size()<<endl;

    char filenameDecSAPhM[FILENAMESIZE];
    sprintf(filenameDecSAPhM,"%s/%s_DecSAV2_T%d_V%d.svg",OUTPUT,file,int(thickness*10),version);
    Board2D DecSAboardM;
    /// Display contour points
    for (vector<Point>::const_iterator it = aContour.begin();it != aContour.end();it++)
        DecSAboardM << SetMode((*it).className(), "Both") << *it;
    /// Display contour points
    /// Display dominant points
    for (vector<Point>::const_iterator it = DP.begin();it != DP.end();it++)
        DecSAboardM << SetMode("PointVector", "Both")
                   << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,0,142)) )//blue
                   << *it;
    /// Display dominant points
    /// Display decoposition (segments)
    DecSAboardM.setPenColor(Color( 255, 0, 0));
    DecSAboardM.setLineWidth(20.0);
    for (vector<Point>::const_iterator it = DecSegmentsPhM.begin(); it != DecSegmentsPhM.end();it++)
        DecSAboardM.drawLine(aContour.at((*it)[0])[0],aContour.at((*it)[0])[1],aContour.at((*it)[1])[0],aContour.at((*it)[1])[1],1);
    /// Display decoposition (segments)
    /// Display decoposition (cirles)
    DecSAboardM.setPenColor(Color( 0, 255, 0));
    for (vector<Point>::const_iterator it = DecArcsPhM.begin(); it != DecArcsPhM.end();it++)
        DecSAboardM.drawLine(aContour.at((*it)[0])[0],aContour.at((*it)[0])[1],aContour.at((*it)[1])[0],aContour.at((*it)[1])[1],1);
    /// Display decoposition (cirles)
    DecSAboardM.saveSVG(filenameDecSAPhM);
    /********* decomposition segment circle : Phuong version modified ***************/

    return 0;


    /****** simple test of circle : ise/lmax/raidus ********
    char filenameCircle[FILENAMESIZE];
    sprintf(filenameCircle,"%s/%s_Circle_T%d_V%d.svg",OUTPUT,file,int(thickness*10),version);
    testCircleSimpleV1(aContour,DP,filenameCircle,false);
    //testCircleSimpleV2(aContour,DP,filenameCircle,false);
    //testCircleSimpleV3(aContour,DP,indexDP,filenameCircle,false);
    ****** simple test of circle : ise/lmax/raidus ********/

    /* Verifying circularity by chord prop *
    vector<Point> DecSegmentsPhM;
    vector<Point> DecArcsPhM;
    char filenameTSBSM[FILENAMESIZE];
    sprintf(filenameTSBSM,"%s/%s_TSBSM_T%d_V%d.svg",OUTPUT,file,int(thickness*10),version);
    //testVerifyCircularityChrodProp(aContour, DP, indexDP, thickness, DecSegmentsPhM,DecArcsPhM,true);
    testVerifyCircularityChrodProp(aContour, newDP, indexNewDP, thickness, DecSegmentsPhM,DecArcsPhM,true);
    char filenameDecSAPhM[FILENAMESIZE];
    sprintf(filenameDecSAPhM,"%s/%s_DecSAV2_T%d_V%d.svg",OUTPUT,file,int(thickness*10),version);
    Board2D DecSAboardM;
    /// Display contour points
    for (vector<Point>::const_iterator it = aContour.begin();it != aContour.end();it++)
        DecSAboardM << SetMode((*it).className(), "Both") << *it;
    /// Display contour points
    /// Display dominant points
    for (vector<Point>::const_iterator it = DP.begin();it != DP.end();it++)
        DecSAboardM << SetMode("PointVector", "Both")
                   << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,0,142)) )//blue
                   << *it;
    /// Display dominant points
    /// Display decoposition (segments)
    DecSAboardM.setPenColor(Color( 255, 0, 0));
    DecSAboardM.setLineWidth(20.0);
    for (vector<Point>::const_iterator it = DecSegmentsPhM.begin(); it != DecSegmentsPhM.end();it++)
        DecSAboardM.drawLine(aContour.at((*it)[0])[0],aContour.at((*it)[0])[1],aContour.at((*it)[1])[0],aContour.at((*it)[1])[1],1);
    /// Display decoposition (segments)
    /// Display decoposition (cirles)
    DecSAboardM.setPenColor(Color( 0, 255, 0));
    for (vector<Point>::const_iterator it = DecArcsPhM.begin(); it != DecArcsPhM.end();it++)
        DecSAboardM.drawLine(aContour.at((*it)[0])[0],aContour.at((*it)[0])[1],aContour.at((*it)[1])[0],aContour.at((*it)[1])[1],1);
    /// Display decoposition (cirles)
    DecSAboardM.saveSVG(filenameDecSAPhM);
    /* Verifying circularity by chord prop */

    //return 0;


}

