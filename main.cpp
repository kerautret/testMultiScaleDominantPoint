#include <iostream>
#include "myfunctions.h"
#include "testfunctions.h"
#include "DGtal/geometry/curves/FreemanChain.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

namespace po = boost::program_options;
using namespace std;
const std::string version = "1.0.1";



int main(int argc, char *const *argv)
{
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "input contour basename (without sdp extension).")
    ("imageneDir,d", po::value<std::string>(), "specify the imagene dir.")
    ("sourceImageWidth", po::value<unsigned int>(), "source image width")
    ("sourceImageHeight", po::value<unsigned int>(), "source image height")
    ("output,o", po::value<std::string>()->default_value("./"), "output dir (default ./).")
    ("version,v", "display the version num")
    ("eps,e","set output with eps format");

  
  bool parseOK=true;
  po::variables_map vm;
  try{
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }catch(const std::exception& ex){
    trace.info()<< "Error checking program options: "<< ex.what()<< std::endl;
    parseOK=false;
  }
  if(vm.count("version")){
    trace.info() << "version: " << version << std::endl;
    return 0;
  }
  
  po::notify(vm);
  if(vm.count("help")||argc<=1|| !parseOK || !vm.count("input") || !vm.count("imageneDir"))
    {
      trace.info()<< "Apply the adaptative tangential cover and polgonalisation." <<std::endl << "Options: "<<std::endl
		  << general_opt << "\n";
      return 0;
    }
  
  bool eps = vm.count("eps");
  bool displayImageCanvas = vm.count("sourceImageWidth") && vm.count("sourceImageHeight");
  
  unsigned int widthCanvas, heightCanvas = 0;
  if (displayImageCanvas){
    widthCanvas = vm["sourceImageWidth"].as<unsigned int>();
    heightCanvas  = vm["sourceImageHeight"].as<unsigned int>();
  }
  

  /********** read data ***************/
  stringstream fileContour;
  string baseInputName  = vm["input"].as<std::string>(); 
  string outDir  = vm["output"].as<std::string>(); 
  string singleName = baseInputName.substr(baseInputName.find_last_of("/")+1);
  std::string outputExt = baseInputName.substr(baseInputName.find_last_of(".")+1);
  if(outputExt != "sdp"){
    trace.error() << "input file should be sdp file" << std::endl;
    return 1;
  }
  singleName = singleName.substr(0, singleName.find_last_of("."));

  trace.info() << singleName << std::endl;
  string ImaGeneDIR = vm["imageneDir"].as<std::string>();
  fileContour << baseInputName ;
  vector<Point> aContourTmp = readFile(fileContour.str().c_str(),1);
  
    
  fileContour.str("");  
  fileContour << outDir << "/" << singleName << ".dat" ; 
  writeFile(aContourTmp, fileContour.str().c_str(), true);
  vector<Point> aContour = readFile(fileContour.str().c_str(),true);
  cout<<fileContour.str()<<endl;
  /********** read data ***************/



  /********** convert to freeman code ***************/
  stringstream fileContourFC;
  FreemanChain<Z2i::Integer> fcContour (aContour);
  fileContourFC.str("");  
  fileContourFC << outDir << "/" << singleName << ".fc" ; 
  
  std::ofstream out(fileContourFC.str().c_str());
  fcContour.selfDisplay(out);
  out<<endl;
  out.close();
  

    /********** call multiscale prog ***************/
     std::stringstream instruction;
    /* meaningful thickness */
     std::stringstream noiseLevelMTFile;
     noiseLevelMTFile << outDir << "/" << singleName << "MeanThickness.txt";
     std::stringstream  smoothContourFile;
     smoothContourFile << outDir << "/"<< singleName << "SmoothContour.txt";
     instruction << ImaGeneDIR << "/build/tests/TestCompNoiseDetect/displayNoiseBS -srcPolygon " << fileContour.str() 
                 << " 0 1 CLOSED -exportNoiseLevel "<< noiseLevelMTFile.str() 
                 << " -displaySmoothContour " << smoothContourFile.str()  ;
     std::system(instruction.str().c_str());
     //read file
     vector<double> vecMT = readMeanindfulThicknessFile(noiseLevelMTFile.str().c_str());//*sqrt(2)
     double glNoise = getGlobalNoise(vecMT);
     cout<<"File : "<<baseInputName <<" => globalNoise (meaningful thickness) = "<<glNoise<<endl;
    //colorate points by its meaningful thickness
     stringstream fileColorMT; 
     fileColorMT << outDir << "/"<< singleName  << (eps ? "_Color.eps" :"_Color.svg");

     drawMeaningfulValue(aContour,vecMT, fileColorMT.str().c_str(), widthCanvas, heightCanvas);



    // /* meaningful thickness */

    if(isRegularNoise(vecMT))
        cout<<"==> Regular noise "<<endl;
    else
        cout<<"==> Not regular noise "<<endl;
    /********** call multiscale prog ***************/

    /******** Adaptive tangent cover ***********/
    vector<AlphaThickSegmentComputer2D> adaptiveTangentCover;
    stringstream fileAdaptMT;
    fileAdaptMT << outDir<< "/"<< singleName<< "ATC";
    bool verbose = false;

    adaptiveTangentCover = testAdaptiveTangentCover(aContour,vecMT,fileAdaptMT.str().c_str(), (eps? "eps": "svg"), 
                                                    verbose, widthCanvas, heightCanvas);
    /******** Adaptive tangent cover ***********/


    

    /******** Dominant point detection with the adaptive tangent cover cover ******/
    vector<Point> DP;
    stringstream filenameDP;
    filenameDP << outDir << "/" << singleName << (eps? "_DP.eps": "_DP.svg");
    bool isClosed = false;
    bool isSymmetry = false;

     cout<<"======= ATC ========";
     DP = testDominantPointOnShape(adaptiveTangentCover,aContour,
                                   isSymmetry,isClosed,filenameDP.str().c_str(),verbose, widthCanvas,heightCanvas);
     cout<<endl<<"===> Num of dominant points is "<<DP.size()<<endl;

     //Find index of DP on the curve
    vector<int> indexDP;
    int indexLastDP=0;
    for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
    {
        indexLastDP = findElement(aContour,*it,indexLastDP);
        if(indexLastDP!=-1)
            indexDP.push_back(indexLastDP);
    }

    if(verbose)
      cout<<"File: "<<fileContour.str()<<endl;
    cout<<"Errors: "; error_All(aContour,DP,indexDP,isClosed);
    /******** Dominant point detection with the adaptive tangent cover cover ******/


    
    /********** Selection of dominant points ***************/
    stringstream filenameDPnew;
    filenameDPnew << outDir<< "/" << singleName << (eps? "_DPnew.eps":"_DPnew.svg");
    vector<Point> newDP = testDominantPointSelectionV2(DP,indexDP,aContour,
                                                       isClosed,filenameDPnew.str().c_str(),verbose, 
                                                       widthCanvas, heightCanvas); // ISE * ANGLE
    cout<<"===> New num of dominant points is "<<newDP.size()<<endl;

    vector<int> indexNewDP;
    int indexLastNewDP=0;
    for(vector<Point>::const_iterator it = newDP.begin(); it != newDP.end(); it++)
    {
        indexLastNewDP = findElement(aContour,*it,indexLastNewDP);
        if(indexLastNewDP!=-1)
            indexNewDP.push_back(indexLastNewDP);
    }

    if(verbose)
      cout<<"File: "<<fileContour.str()<<endl;
    cout<<"Errors: "; error_All(aContour,newDP,indexNewDP,isClosed);
    /********** Selection of dominant points ***************/


    

    cout<<"======= Mid ========";
    double thickness = glNoise;
    vector<AlphaThickSegmentComputer2D> fuzzySegmentSetMid = blurredSegmentCurveDecomposition(aContour,thickness,NULL,false, widthCanvas, heightCanvas);
    vector<Point> DPM;
    stringstream filenameDPM;
    filenameDPM << outDir << "/" << singleName  << (eps?"_DPM.eps" :"_DPM.svg");
   
    //cout<<"===> Num of seg decomposed is "<<fuzzySegmentSetMid.size()<<endl;
    DPM = testDominantPointOnShape(fuzzySegmentSetMid,aContour,isSymmetry,isClosed,
                                   filenameDPM.str().c_str(),verbose, widthCanvas, heightCanvas);
    cout<<endl<<"===> Num of dominant points is "<<DPM.size()<<endl;

    //Find index of DP on the curve
    vector<int> indexDPM;
    int indexLastDPM=0;
    for(vector<Point>::const_iterator it = DPM.begin(); it != DPM.end(); it++)
    {
        indexLastDPM = findElement(aContour,*it,indexLastDPM);
        //cout<<"indexLastDPM = "<<indexLastDPM<<" Point is "<<(*it)<<endl;
        if(indexLastDPM!=-1)
            indexDPM.push_back(indexLastDPM);
    }

    if(verbose)
      cout<<"File: "<<fileContour.str()<<endl;
    cout<<"Errors: "; error_All(aContour,DPM,indexDPM,isClosed);
    // /******** Dominant point detection with the adaptive tangent cover cover ******/

    /********** Selection of dominant points ***************/
    stringstream filenameDPnewV;
    filenameDPnewV << outDir << "/" << singleName << (eps? "_DPnewV.eps":"_DPnewV.svg");

    vector<Point> newDPM = testDominantPointSelectionV2(DPM,indexDPM,aContour,isClosed,
                                                        filenameDPnewV.str().c_str(),verbose, widthCanvas, heightCanvas); // ISE * ANGLE

    cout<<"===> New num of dominant points is "<<newDPM.size()<<endl;
    vector<int> indexNewDPM;
    int indexLastNewDPM=0;
    for(vector<Point>::const_iterator it = newDPM.begin(); it != newDPM.end(); it++)
    {
        indexLastNewDPM = findElement(aContour,*it,indexLastNewDPM);
        if(indexLastNewDPM!=-1)
            indexNewDPM.push_back(indexLastNewDPM);
    }

    if(verbose)
      cout<<"File: "<<fileContour.str()<<endl;
    cout<<"Errors: "; error_All(aContour,newDPM,indexNewDPM,isClosed);
    /********** Selection of dominant points ***************/

    return 0;
}

