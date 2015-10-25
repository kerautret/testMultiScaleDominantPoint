#include "testfunctions.h"
#include "myfunctions.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"

/*************************************/
/*** Burred segments decomposition ***/
/*************************************/
vector<AlphaThickSegmentComputer2D> blurredSegmentCurveDecomposition(const vector<Point>& aContour, double thickness, const char* filename, bool verbose)
{
    double mu, nu;
    PointD N;
    vector<AlphaThickSegmentComputer2D> fuzzySegmentSet;
    int index=0;
    std::string n (filename);
    std::string outputExt = n.substr(n.find_last_of(".")+1);
    //run over the points on the contours
    for (vector<Point>::const_iterator it = aContour.begin();it != aContour.end();it++)
    {
        //cout<<*it<<endl;
      AlphaThickSegmentComputer2D aSegment(thickness);
        aSegment.init(it);
        /* travel over the contour points and add to the seg */
        //if(aSegment.end() != aContour.end()) aSegment.extendFront();
        while (aSegment.end() != aContour.end() && aSegment.extendFront()){}
        //if(*it == Point(17,12))
        //    cout<<"==> STOP : getStartPoint = "<<findElement(aContour,getStartPoint(aSegment))<<", getEndPoint = "<<findElement(aContour,getEndPoint(aSegment))<<" and EndBack = "<<findElement(aContour,getEndPoint(fuzzySegmentSet.back()))<<endl;
        if(it == aContour.begin())
        {
            if(verbose)
            {
                aSegment.computeParallelStripParams(mu,N,nu);
                cout<<fuzzySegmentSet.size()<<" with slope = "<<getSlope(N[0],N[1])<<" : ("<<N[0]<<" , "<<N[1]<<" , "<<mu<<" , "<<nu<<" ) ==> ";
                cout<<getStartPoint(aSegment)<<" and "<<getEndPoint(aSegment)<<endl;
            }
            fuzzySegmentSet.push_back(aSegment);
        }
        else if(findElement(aContour,getEndPoint(aSegment),index) > findElement(aContour,getEndPoint(fuzzySegmentSet.back()),index))
        {
            if(verbose)
            {
                aSegment.computeParallelStripParams(mu,N,nu);
                cout<<fuzzySegmentSet.size()<<" with slope = "<<getSlope(N[0],N[1])<<" : ("<<N[0]<<" , "<<N[1]<<" , "<<mu<<" , "<<nu<<" ) ==> ";
                cout<<getStartPoint(aSegment)<<" and "<<getEndPoint(aSegment)<<endl;
            }
            fuzzySegmentSet.push_back(aSegment);
        }
        if(getEndPoint(aSegment) == aContour.back())
            break;
        index++;
    }

    /****** write plot file **********/
    if(IS_PLOT)
    {
        /* Display points in common zone */
        vector<int> idBegin;
        vector<int> idEnd;
        int indexBegin=0;
        for (vector<AlphaThickSegmentComputer2D>::const_iterator it = fuzzySegmentSet.begin();it != fuzzySegmentSet.end();it++)
        {
            indexBegin = findElement(aContour,getStartPoint(*it),indexBegin);
            idBegin.push_back(indexBegin);
            idEnd.push_back(findElement(aContour,getEndPoint(*it),indexBegin));
        }
        writeFile(idBegin,"idBeginSeg.txt",false);
        writeFile(idEnd,"idEndSeg.txt",false);
        /* Display points in common zone */
    }
    /****** write plot file **********/

    if(filename != NULL)
    {
        Board2D aBoard;
        //int count=0;
        /* display the boundary */
        aBoard << SetMode("PointVector", "Both");
        for (vector<Point>::const_iterator it = aContour.begin(); it != aContour.end(); it++)
            aBoard << *it;
        /* display the boundary */

        /* Display boundingbox */
        for (vector<AlphaThickSegmentComputer2D>::const_iterator it = fuzzySegmentSet.begin();it != fuzzySegmentSet.end();it++)
            //if(count==0 || count==1 || count==2 || count==3)
        {
            aBoard << SetMode((*it).className(), "BoundingBox") << *it;
        }
        //count++;
        /* Display boundingbox */
        if(outputExt=="svg"){
          aBoard.saveSVG(filename);
        }
        else if (outputExt == "eps"){
           aBoard.saveEPS(filename);
        }
    }
    return fuzzySegmentSet;
}
/*************************************/
/*** Burred segments decomposition ***/
/*************************************/

/************************************/
/**** Dominant points detections ****/
/************************************/
vector<Point> testDominantPointOnShape(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, bool isSymmetry, bool isClosed, const char* filename, bool verbose)
{
  std::string n (filename);
  std::string outputExt = n.substr(n.find_last_of(".")+1);
  int p=1, q=0, Eq=0, Bp=0, m=fuzzySegmentSet.size();
    vector<Point> pile;
    while (p<m && q<m && Bp>=0 && Eq>=0)
    {
        //Bp = findStartElement(aContour,fuzzySegmentSet[p]);
        //Eq = findEndElement(aContour,fuzzySegmentSet[q]);
        Bp = findElement(aContour,getStartPoint(fuzzySegmentSet[p]),Bp);
        Eq = findElement(aContour,getEndPoint(fuzzySegmentSet[q]),Bp);
        //cout<<endl;
        //cout<<"( "<<Bp<<","<<Eq<<" )"<<endl;
        //cout<<endl<<"findStartElement "<<getStartPoint(fuzzySegmentSet[p])<<" findEndElement "<<getEndPoint(fuzzySegmentSet[q])<<endl;
        while(Eq>=Bp && p<m && q<m && Bp>=0 && Eq>=0 && Bp<aContour.size() && Eq<aContour.size())
        {
            //Bp = findStartElement(aContour,fuzzySegmentSet[p]);
            //Eq = findEndElement(aContour,fuzzySegmentSet[q]);
            Bp = findElement(aContour,getStartPoint(fuzzySegmentSet[p]),Bp);
            //Eq = findElement(aContour,getEndPoint(fuzzySegmentSet[q]),Bp);
            //cout<<"findStartElement "<<getStartPoint(fuzzySegmentSet[p])<<" findEndElement "<<getEndPoint(fuzzySegmentSet[q])<<endl;
            //cout<<"( "<<Bp<<","<<Eq<<" )"<<endl;
            if(verbose && Eq>=Bp)
                cout<<"Pair segment : ( "<<q<<","<<p<<" ) and Pair points : ( "<<Bp<<","<<Eq<<" )"<<endl;
            p++;
        }
        if(Eq<Bp)
        {
            p--;//when out of loop => Eq<Bp
            pile.push_back(Point(q,p-1));
            if(verbose)
                cout<<"Pile : ("<<q<<" , "<<p-1<<" ) "<<endl<<endl;
            q=p-1;
        }
        else
        {
            if(p>=m || q>=m)
            {
                pile.push_back(Point(q,p-1));
                if(verbose)
                    cout<<"Pile : ("<<q<<" , "<<p-1<<" ) "<<endl<<endl;
            }
            else
                q = p;
        }
    }
    /* pile construction */

    /* calcul point dominant */
    vector<Point> DP, extremityPoint, commonPoint;
    int indexB = 0, indexE = 0, indexC = 0, indexBi = 0, indexEi = 0;
    for(vector<Point>::const_iterator it = pile.begin(); it != pile.end(); it++)
    {
        q = (*it)[0];
        p = (*it)[1];
        //find point with max (local) curvature
        /*
        int indexB = findStartElement(aContour,fuzzySegmentSet[p]);
        int indexE = findEndElement(aContour,fuzzySegmentSet[q]);
        int indexC = indexB;
        int indexBi = findStartElement(aContour,fuzzySegmentSet[q]);
        int indexEi = findEndElement(aContour,fuzzySegmentSet[p]);
        */
        indexB = findElement(aContour,getStartPoint(fuzzySegmentSet[p]),indexB) ;
        indexE = findElement(aContour,getEndPoint(fuzzySegmentSet[q]),indexB);
        indexC = indexB;
        indexBi = findElement(aContour,getStartPoint(fuzzySegmentSet[q]),indexBi);
        indexEi = findElement(aContour,getEndPoint(fuzzySegmentSet[p]),indexBi);
        if(verbose)
            cout<<"angle at "<<indexC<<" is "<<(acuteAngle(aContour.at(indexBi),aContour.at(indexC),aContour.at(indexEi))*180/M_PI)<<endl;

        for(int i = indexB+1; i<=indexE; i++)
        {
            int minIndex = 0;
            int modIndex = i;
            if(isSymmetry)
            {
                minIndex = fabs(modIndex-indexB) < fabs(modIndex-indexE) ? fabs(modIndex-indexB) : fabs(modIndex-indexE);
                if(fabs(acuteAngle(aContour.at(modIndex),aContour.at(modIndex-minIndex),aContour.at(modIndex+minIndex))) <=
                        fabs(acuteAngle(aContour.at(indexC),aContour.at(modIndex-minIndex),aContour.at(modIndex+minIndex))))
                    indexC = modIndex;
            }
            else
            {
                if(fabs(acuteAngle(aContour.at(indexBi),aContour.at(modIndex),aContour.at(indexEi))) <=
                    fabs(acuteAngle(aContour.at(indexBi),aContour.at(indexC),aContour.at(indexEi))))
                    indexC = modIndex;
            }
            if(verbose)
                cout<<"angle at "<<modIndex<<" is "<<(acuteAngle(aContour.at(indexBi-minIndex),aContour.at(modIndex),aContour.at(indexEi+minIndex))*180/M_PI)<<endl;

        }
        if(verbose)
            cout<<" ==> max local curvature at "<<aContour[indexC]<<endl;

        if(DP.empty() || aContour[indexC] != DP.back())
        {
            DP.push_back(aContour[indexC]);
            extremityPoint.push_back(aContour[indexB]);
            extremityPoint.push_back(aContour[indexE]);
            for(int i=indexB; i<=indexE; i++)
                commonPoint.push_back(aContour[i]);
        }

    }
    DP.insert(DP.begin(),aContour.front());
    if(!isClosed)//donot consider the last elt
        DP.insert(DP.end(),aContour.back());
    /* calcul point dominant */

    /****** write plot file **********/
    if(IS_PLOT)
    {
        /* Display points in common zone */
        vector<int> idCommonPoint;
        int index = 0;
        for(vector<Point>::const_iterator it = commonPoint.begin(); it != commonPoint.end(); it++)
        {
            index = findElement(aContour,*it,index);
            idCommonPoint.push_back(index);
        }
        vector<int> idDominantPoint;
        index = 0;
        for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
        {
            index = findElement(aContour,*it,index);
            idDominantPoint.push_back(index);
        }
        writeFile(idCommonPoint,"idCommonPoint.txt",false);
        writeFile(idDominantPoint,"idDominantPoint.txt",false);
        /* Display points in common zone */
    }
    /****** write plot file **********/

    if(filename != NULL)
    {
        Board2D aBoard;
        /* Display the boundary */
        aBoard << SetMode("PointVector", "Both");
        for (vector<Point>::const_iterator it = aContour.begin(); it != aContour.end(); it++)
            aBoard << *it;
        /* Display the boundary */
        // Display boundingbox
        //for (vector<AlphaThickSegmentComputer2D>::const_iterator it = fuzzySegmentSet.begin();it != fuzzySegmentSet.end();it++)
        //    aBoard << SetMode((*it).className(), "BoundingBox")
        //            << *it;
        // Display boundingbox
        /* Display points in common zone */
        for(vector<Point>::const_iterator it = commonPoint.begin(); it != commonPoint.end(); it++)
            aBoard << SetMode((*it).className(), "Both")
                   << CustomStyle( (*it).className(), new CustomColors( Color(0,255,0), Color(0,92,0)) )
                   <<(*it);
        /* Display points in common zone */
        /* Display dominant point */
        for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
            aBoard << SetMode((*it).className(), "Both")
                   << CustomStyle( (*it).className(), new CustomColors( Color(255,0,0), Color(192,0,0)) )
                   <<(*it);
        /* Display dominant point */
        /* Display the segments by DP */
        aBoard.setPenColor(Color(0, 0, 255));
        aBoard.setLineWidth(10.0);
        //aBoard.drawLine(aContour.front()[0],aContour.front()[1],DP.front()[0],DP.front()[1]);
        for(vector<Point>::const_iterator it = DP.begin(); it+1 != DP.end(); it++)
            aBoard.drawLine((*it)[0],(*it)[1],(*(it+1))[0],(*(it+1))[1]);
        if(isClosed)//NOTE : DP.front() = aContour.front()
            aBoard.drawLine(aContour.front()[0],aContour.front()[1],DP.back()[0],DP.back()[1]);
        /* Display the segments by DP */
        if(outputExt=="svg"){
          aBoard.saveSVG(filename);
        }
        else if (outputExt == "eps"){
           aBoard.saveEPS(filename);
        }
           
    }
    return DP;
}
/************************************/
/**** Dominant points detections ****/
/************************************/


/***********************************/
/**** Dominant points selection ****/
/***********************************/
/* Selection by ise*angle */
vector<Point> testDominantPointSelectionV1(const vector<Point>& DP, const vector<int>& indexDP, int nbDP, const vector<Point>& aContour, bool isClosed, const char* filename, bool verbose)
{
  std::string n (filename);
  std::string outputExt = n.substr(n.find_last_of(".")+1);
    vector<Point> selectedDP;
    for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
        selectedDP.push_back(*it);
    vector<int> indexSelectedDP;
    for(vector<int>::const_iterator it = indexDP.begin(); it != indexDP.end(); it++)
        indexSelectedDP.push_back(*it);
    vector<double> iseDP;
    if(nbDP<DP.size())
    {
        do{
            //cout<<"====> selectedDP.size() = "<<selectedDP.size()<<endl;
            for(int it = 1; it+1 < selectedDP.size(); it++)
            {
                double ise = 0.0, d = 0.0, angle = 0.0;
                int indexStart = indexSelectedDP.at(it-1);
                int indexEnd = indexSelectedDP.at(it+1);
                if(indexStart < 0 || indexEnd < 0)
                    cout<<"Pb in error_ISE index of DP"<<endl;
                for(int index = indexStart+1; index < indexEnd; index++)
                {
                    d = distancePointSegment(aContour.at(index),aContour.at(indexStart), aContour.at(indexEnd));
                    ise += d*d;
                }
                angle = acuteAngle(aContour.at(indexStart),
                                   aContour.at(indexSelectedDP.at(it)),
                                   aContour.at(indexEnd));

                if(verbose)
                    cout<<" DP "<<it<<" : ise = "<<ise<<", and angle = "<<angle<<endl;

                //iseDP.push_back(ise/(indexEnd - indexStart - 1));
                iseDP.push_back(ise/angle);
            }
            iseDP.insert(iseDP.begin(),aContour.size()*aContour.size());//Keep the fist point as DP
            if(!isClosed)//Keep the last point as DP
                iseDP.push_back(aContour.size()*aContour.size());
            else//isClosed => calcule ISE of the last point
            {
                double ise = 0.0, d = 0.0, angle = 0.0;
                int indexStart = indexSelectedDP.at(indexSelectedDP.size()-2);
                int indexEnd = aContour.size();
                //cout<<"it = "<<it<<" : indexStart = "<<indexStart<<", indexEnd = "<<indexEnd<<endl;
                for(int index = indexStart+1; index < indexEnd; index++)
                {
                    d = distancePointSegment(aContour.at(index),aContour.at(indexStart), aContour.front());
                    //cout<<"index "<<index<<" : distance "<<d<<endl;
                    ise += d*d;
                }
                angle = acuteAngle(aContour.at(indexStart),
                                   aContour.at(indexSelectedDP.at(indexSelectedDP.size()-1)),
                                   aContour.front());

                if(verbose)
                    cout<<" DP "<<indexSelectedDP.at(indexSelectedDP.size()-1)<<" : ise = "<<ise<<", and angle = "<<angle<<endl;

                //iseDP.push_back(ise/(indexEnd - indexStart - 1));
                iseDP.push_back(ise/angle);
            }

            vector<int> indexOrderedDP = absSortIndex(iseDP,true);
            if(verbose)
                for(vector<int>::const_iterator it = indexOrderedDP.begin(); it != indexOrderedDP.end(); it++)
                    cout<<"ise of DP "<<*it<< " : "<<iseDP.at(*it)<<endl;

            /* erease the first elt */
            //cout<<"=> min ise at "<<indexOrderedDP.front()<<endl;
            selectedDP.erase(selectedDP.begin()+indexOrderedDP.front());
            indexSelectedDP.erase(indexSelectedDP.begin()+indexOrderedDP.front());
            iseDP.clear();
            /* erease the first elt */
        }while(selectedDP.size()>nbDP);
    }

    if(filename != NULL)
    {
        Board2D aBoard;
        /* Display the boundary */
        aBoard << SetMode("PointVector", "Both");
        for (vector<Point>::const_iterator it = aContour.begin(); it != aContour.end(); it++)
            aBoard << *it;
        /* Display the boundary */
        /* Display old dominant point */
        for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
            aBoard << SetMode((*it).className(), "Both")
                   << CustomStyle( (*it).className(), new CustomColors( Color(255,0,0), Color(192,0,0)) )//red
                   <<(*it);
        /* Display old dominant point */
        /* Display new dominant point */
        for(vector<Point>::const_iterator it = selectedDP.begin(); it != selectedDP.end(); it++)
            aBoard << SetMode((*it).className(), "Both")
                   << CustomStyle( (*it).className(), new CustomColors( Color(0,255,0), Color(0,192,0)) )//green
                   <<(*it);
        /* Display new dominant point */

        /* Display the segments by old DP */
        aBoard.setPenColor(Color(192, 0, 0));
        aBoard.setLineWidth(10.0);
        aBoard.drawLine(aContour.front()[0],aContour.front()[1],DP.front()[0],DP.front()[1]);
        for(vector<Point>::const_iterator it = DP.begin(); it+1 != DP.end(); it++)
            aBoard.drawLine((*it)[0],(*it)[1],(*(it+1))[0],(*(it+1))[1]);
        if(isClosed)//NOTE : DP.front() = aContour.front()
            aBoard.drawLine(aContour.front()[0],aContour.front()[1],DP.back()[0],DP.back()[1]);
        /* Display the segments by old DP */
        /* Display the segments by new DP */
        aBoard.setPenColor(Color(0, 192, 0));
        for(vector<Point>::const_iterator it = selectedDP.begin(); it+1 != selectedDP.end(); it++)
            aBoard.drawLine((*it)[0],(*it)[1],(*(it+1))[0],(*(it+1))[1]);
        if(isClosed)//NOTE : DP.front() = aContour.front()
            aBoard.drawLine(aContour.front()[0],aContour.front()[1],selectedDP.back()[0],selectedDP.back()[1]);
        /* Display the segments by new DP */
        if(outputExt=="svg"){
          aBoard.saveSVG(filename);
        }
        else if (outputExt == "eps"){
           aBoard.saveEPS(filename);
        }

    }
    return selectedDP;
}
/* Selection by ise*angle */

/* Selection by max error criterion */
vector<Point> testDominantPointSelectionV2(const vector<Point>& DP, const vector<int>& indexDP, const vector<Point>& aContour, bool isClosed, const char* filename, bool verbose)
{
  std::string n (filename);
  std::string outputExt = n.substr(n.find_last_of(".")+1);

    vector<Point> selectedDP;
    for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
        selectedDP.push_back(*it);
    vector<int> indexSelectedDP;
    for(vector<int>::const_iterator it = indexDP.begin(); it != indexDP.end(); it++)
        indexSelectedDP.push_back(*it);
    vector<double> iseDP;
    double ise = error_ISE(aContour,selectedDP,indexSelectedDP,isClosed);
    double cr = error_CR(aContour,selectedDP);
    double error_prev, error;
    error_prev = error = (cr*cr)/ise;
    int lastIndex=-1,lastIndexSelectedDP=-1;
    Point lastSelectedDP;

    while(error_prev <= error){
        error_prev = error;
        //cout<<"nd DP = "<<selectedDP.size()<<" errror = "<<error<<endl;
        //cout<<"====> selectedDP.size() = "<<selectedDP.size()<<endl;
        for(int it = 1; it+1 < selectedDP.size(); it++)
        {
            double ise = 0.0, d = 0.0, angle = 0.0;
            int indexStart = indexSelectedDP.at(it-1);
            int indexEnd = indexSelectedDP.at(it+1);
            if(indexStart < 0 || indexEnd < 0)
                cout<<"Pb in error_ISE index of DP"<<endl;
            for(int index = indexStart+1; index < indexEnd; index++)
            {
                d = distancePointSegment(aContour.at(index),aContour.at(indexStart), aContour.at(indexEnd));
                ise += d*d;
            }
            angle = acuteAngle(aContour.at(indexStart),
                               aContour.at(indexSelectedDP.at(it)),
                               aContour.at(indexEnd));

            if(verbose)
                cout<<" DP "<<it<<" : ise = "<<ise<<", and angle = "<<angle<<endl;

            //iseDP.push_back(ise/(indexEnd - indexStart - 1));
            iseDP.push_back(ise/angle);
        }
        iseDP.insert(iseDP.begin(),aContour.size()*aContour.size());//Keep the fist point as DP
        if(!isClosed)//Keep the last point as DP
            iseDP.push_back(aContour.size()*aContour.size());
        else//isClosed => calcule ISE of the last point
        {
            double ise = 0.0, d = 0.0, angle = 0.0;
            int indexStart = indexSelectedDP.at(indexSelectedDP.size()-2);
            int indexEnd = aContour.size();
            //cout<<"it = "<<it<<" : indexStart = "<<indexStart<<", indexEnd = "<<indexEnd<<endl;
            for(int index = indexStart+1; index < indexEnd; index++)
            {
                d = distancePointSegment(aContour.at(index),aContour.at(indexStart), aContour.front());
                //cout<<"index "<<index<<" : distance "<<d<<endl;
                ise += d*d;
            }
            angle = acuteAngle(aContour.at(indexStart),
                               aContour.at(indexSelectedDP.at(indexSelectedDP.size()-1)),
                               aContour.front());

            if(verbose)
                cout<<" DP "<<indexSelectedDP.at(indexSelectedDP.size()-1)<<" : ise = "<<ise<<", and angle = "<<angle<<endl;

            //iseDP.push_back(ise/(indexEnd - indexStart - 1));
            iseDP.push_back(ise/angle);
        }

        vector<int> indexOrderedDP = absSortIndex(iseDP,true);
        if(verbose)
            for(vector<int>::const_iterator it = indexOrderedDP.begin(); it != indexOrderedDP.end(); it++)
                cout<<"ise of DP "<<*it<< " : "<<iseDP.at(*it)<<endl;

        /* erease the first elt */
        lastIndex = indexOrderedDP.front();
        lastSelectedDP = selectedDP.at(lastIndex);
        lastIndexSelectedDP = indexSelectedDP.at(lastIndex);
        //cout<<"=> min ise at "<<indexOrderedDP.front()<<endl;
        selectedDP.erase(selectedDP.begin()+indexOrderedDP.front());
        indexSelectedDP.erase(indexSelectedDP.begin()+indexOrderedDP.front());
        iseDP.clear();
        /* erease the first elt */

        /* update errors */
        ise = error_ISE(aContour,selectedDP,indexSelectedDP,isClosed);
        cr = error_CR(aContour,selectedDP);
        error = (cr*cr)/ise;
        /* update errors */
    }
    /* Push back the last element */
    if(lastIndex !=-1)
    {
        selectedDP.insert(selectedDP.begin()+lastIndex,lastSelectedDP);
        indexSelectedDP.insert(indexSelectedDP.begin()+lastIndex,lastIndexSelectedDP);
    }
    /* Push back the last element */

    if(filename != NULL)
    {
        Board2D aBoard;
        /* Display the boundary */
        aBoard << SetMode("PointVector", "Both");
        for (vector<Point>::const_iterator it = aContour.begin(); it != aContour.end(); it++)
            aBoard << *it;
        /* Display the boundary */
        /* Display old dominant point */
        for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
            aBoard << SetMode((*it).className(), "Both")
                   << CustomStyle( (*it).className(), new CustomColors( Color(255,0,0), Color(192,0,0)) )//red
                   <<(*it);
        /* Display old dominant point */
        /* Display new dominant point */
        for(vector<Point>::const_iterator it = selectedDP.begin(); it != selectedDP.end(); it++)
            aBoard << SetMode((*it).className(), "Both")
                   << CustomStyle( (*it).className(), new CustomColors( Color(0,255,0), Color(0,192,0)) )//green
                   <<(*it);
        /* Display new dominant point */

        /* Display the segments by old DP */
        aBoard.setPenColor(Color(192, 0, 0));
        aBoard.setLineWidth(10.0);
        aBoard.drawLine(aContour.front()[0],aContour.front()[1],DP.front()[0],DP.front()[1]);
        for(vector<Point>::const_iterator it = DP.begin(); it+1 != DP.end(); it++)
            aBoard.drawLine((*it)[0],(*it)[1],(*(it+1))[0],(*(it+1))[1]);
        if(isClosed)//NOTE : DP.front() = aContour.front()
            aBoard.drawLine(aContour.front()[0],aContour.front()[1],DP.back()[0],DP.back()[1]);
        /* Display the segments by old DP */
        /* Display the segments by new DP */
        aBoard.setPenColor(Color(0, 192, 0));
        for(vector<Point>::const_iterator it = selectedDP.begin(); it+1 != selectedDP.end(); it++)
            aBoard.drawLine((*it)[0],(*it)[1],(*(it+1))[0],(*(it+1))[1]);
        if(isClosed)//NOTE : DP.front() = aContour.front()
            aBoard.drawLine(aContour.front()[0],aContour.front()[1],selectedDP.back()[0],selectedDP.back()[1]);
        /* Display the segments by new DP */
        //sprintf(filename,"%s_N%d.svg",filename,selectedDP.size());
        if(outputExt=="svg"){
          aBoard.saveSVG(filename);
        }
        else if (outputExt == "eps"){
           aBoard.saveEPS(filename);
        }
    }
    return selectedDP;
}
/* Selection by max error criterion */
/************************************/
/**** Dominant points selection ****/
/************************************/

/********************************************************/
/**************** Draw Multi-Thickness Cover ************/
/********************************************************/
void drawMultiThicknessCover(const vector<Point>& aContour, const vector<double>& thckVect, int nbColor, char* filename)
{
  std::string n (filename);
  std::string outputExt = n.substr(n.find_last_of(".")+1);

    Board2D aBoard;
    /* display the boundary */
    string className = "PointVector/Both";
    // Creating colormap
    HueShadeColorMap<double> hueMap(0.9, nbColor+1.0);

    //for (vector<Point>::const_iterator it = 0; it != aContour.end(); it++)
    for (int it = 0; it < aContour.size(); it++)
    {
        unsigned int c = (unsigned int)thckVect.at(it);
        aBoard << SetMode("PointVector", "Both")
               << CustomStyle(className, new CustomColors( Color::Gray, hueMap(c) ))
               << aContour.at(it);
    }
    /* display the boundary */
        if(outputExt=="svg"){
          aBoard.saveSVG(filename);
        }
        else if (outputExt == "eps"){
           aBoard.saveEPS(filename);
        }

}

void drawMultiThicknessCover(const vector<Point>& aContour, const vector<vector<AlphaThickSegmentComputer2D> >& meaningThicknessTangentCover, char* filename)
{
  std::string n (filename);
  std::string outputExt = n.substr(n.find_last_of(".")+1);

    Board2D aBoard;
    /* display the boundary */
    aBoard << SetMode("PointVector", "Both");
    for (vector<Point>::const_iterator it = aContour.begin(); it != aContour.end(); it++)
        aBoard << *it;
    /* display the boundary */

    /* Display boundingbox */
    // Draw each segment
    string className = "AlphaThickSegment/BoundingBox";
    unsigned int count = 1;
    // Creating colormap
    HueShadeColorMap<double> hueMap(0.9, meaningThicknessTangentCover.size()+1.0);
    for(vector<vector<AlphaThickSegmentComputer2D> >::const_iterator it = meaningThicknessTangentCover.begin(); it != meaningThicknessTangentCover.end(); it++)
    {
        vector<AlphaThickSegmentComputer2D> fuzzySegmentSet = *it;
        for (vector<AlphaThickSegmentComputer2D>::const_iterator it_bis = fuzzySegmentSet.begin();it_bis != fuzzySegmentSet.end();it_bis++)
        {
            aBoard << SetMode((*it_bis).className(), "BoundingBox")
                   << CustomStyle( className,  new CustomPenColor( hueMap(count) ) )
                   << *it_bis;
        }
        count++;
    }
    /* Display boundingbox */
    if(outputExt=="svg"){
          aBoard.saveSVG(filename);
        }
        else if (outputExt == "eps"){
           aBoard.saveEPS(filename);
        }


}

void drawMultiThicknessCover(const vector<Point>& aContour, const vector<AlphaThickSegmentComputer2D>& meaningThicknessTangentCover, char* filename)
{
  std::string n (filename);
  std::string outputExt = n.substr(n.find_last_of(".")+1);

    Board2D aBoard;
    /* display the boundary */
    aBoard << SetMode("PointVector", "Both");
    for (vector<Point>::const_iterator it = aContour.begin(); it != aContour.end(); it++)
        aBoard << *it;
    /* display the boundary */

    /* Display boundingbox */
    for(vector<AlphaThickSegmentComputer2D>::const_iterator it = meaningThicknessTangentCover.begin(); it != meaningThicknessTangentCover.end(); it++)
    {
        aBoard << SetMode((*it).className(), "BoundingBox") //string className = "AlphaThickSegment/BoundingBox";
               << *it;
    }
    /* Display boundingbox */
        if(outputExt=="svg"){
          aBoard.saveSVG(filename);
        }
        else if (outputExt == "eps"){
           aBoard.saveEPS(filename);
        }

}

void drawMultiThicknessCover(const vector<Point>& aContour, const vector<vector<AlphaThickSegmentComputer2D> >& meaningThicknessTangentCover, const vector<double>& thckVect, char* filename)
{
  std::string n (filename);
  std::string outputExt = n.substr(n.find_last_of(".")+1);

    Board2D aBoard;
    // Creating colormap
    HueShadeColorMap<double> hueMap(0.9, meaningThicknessTangentCover.size()+1.0);
    unsigned int count = 1;
    /* display the boundary */
    string classNamePoint = "PointVector/Both";
    for (int it = 0; it < aContour.size(); it++)
    {
        unsigned int c = (unsigned int)thckVect.at(it);
        aBoard << SetMode("PointVector", "Both")
               << CustomStyle(classNamePoint, new CustomColors( Color::Gray, hueMap(c) ))
               << aContour.at(it);
    }
    /* display the boundary */

    /* Display boundingbox */
    string classNameSeg = "AlphaThickSegment/BoundingBox";
    for(vector<vector<AlphaThickSegmentComputer2D> >::const_iterator it = meaningThicknessTangentCover.begin(); it != meaningThicknessTangentCover.end(); it++)
    {
        vector<AlphaThickSegmentComputer2D> fuzzySegmentSet = *it;
        for (vector<AlphaThickSegmentComputer2D>::const_iterator it_bis = fuzzySegmentSet.begin();it_bis != fuzzySegmentSet.end();it_bis++)
        {
            ///aBoard << SetMode((*it_bis).className(), "BoundingBox") << *it_bis;
            aBoard << SetMode((*it_bis).className(), "BoundingBox")
                   << CustomStyle( classNameSeg,  new CustomPenColor( hueMap( count ) ) )
                   << *it_bis;
        }
        count++;
    }
    /* Display boundingbox */
        if(outputExt=="svg"){
          aBoard.saveSVG(filename);
        }
        else if (outputExt == "eps"){
           aBoard.saveEPS(filename);
        }


}

void drawMeaningfulValue(const vector<Point>& aContour, const vector<double> vecMeanVal, const char* filename)
{
  std::string n (filename);
  std::string outputExt = n.substr(n.find_last_of(".")+1);

    Board2D aBoard;
    aBoard << SetMode("PointVector", "Both");
    int minPos = 0, maxPos = 0;
    for (unsigned int i = 0; i < vecMeanVal.size(); i++)
    {
        if (vecMeanVal.at(i) < vecMeanVal.at(minPos)) // Found a smaller min
            minPos = i;
        if (vecMeanVal.at(i) > vecMeanVal.at(maxPos)) // Found a bigger max
            maxPos = i;
    }
    double min = vecMeanVal.at(minPos);
    double max = vecMeanVal.at(maxPos);
    int a = 100;
    int b = 255;
    //cout<<"min = "<<min<<", max = "<<max<<endl;
    //convert : (min,max) --> (a,b)
    //f(x) = (b-a)(x-min)/(max-min) + a
    double step = (b-a)/(max-min+1);
    vector<int> gradColor;
    for(vector<double>::const_iterator it = vecMeanVal.begin(); it != vecMeanVal.end(); it++)
    {
        int val = round(step*((*it)-min) + a);
        gradColor.push_back(val);
    }

    for (unsigned int i = 0; i < gradColor.size(); i++)
    {
        aBoard << CustomStyle("PointVector", new CustomColors(Color(gradColor.at(i),0,0),Color(gradColor.at(i),0,0)))
               << aContour.at(i);
    }
    /* display the boundary */
        if(outputExt=="svg"){
          aBoard.saveSVG(filename);
        }
        else if (outputExt == "eps"){
           aBoard.saveEPS(filename);
        }


}
/********************************************************/
/**************** Draw Multi-Thickness Cover ************/
/********************************************************/

/****************************************************************/
/*********** Adaptive tangent cover computation *****************/
/****************************************************************/
vector<AlphaThickSegmentComputer2D> testAdaptiveTangentCover(const vector<Point>& aContour, const vector<double>& vecMT, const char* filename, bool verbose)
{
  std::string n (filename);
  std::string outputExt = n.substr(n.find_last_of(".")+1);

    //1. Find vector of thickness element
    vector<double> meaningThicknessElement;
    meaningThicknessElement.push_back(vecMT.front());
    for(vector<double>::const_iterator it = vecMT.begin()+1; it != vecMT.end(); it++)
    {
        double m = (*it);
        if (std::find(meaningThicknessElement.begin(), meaningThicknessElement.end(),m)==meaningThicknessElement.end())
            meaningThicknessElement.push_back(m);
    }
    std::sort(meaningThicknessElement.begin(),meaningThicknessElement.end(),sort_increase);
    for(vector<double>::const_iterator it = meaningThicknessElement.begin(); it != meaningThicknessElement.end(); it++)
        cout<<"meaningThicknessElement : "<<*it<<endl;
    char fileAdaptMT[FILENAMESIZE];
    if(filename != NULL)
        sprintf(fileAdaptMT,"%s_Step1.svg",filename);
    drawMultiThicknessCover(aContour,vecMT,meaningThicknessElement.size(),fileAdaptMT);

    //2. Compute different thickness tangent covers (blurred segments)
    vector<vector<AlphaThickSegmentComputer2D> > meaningThicknessTangentCover(meaningThicknessElement.size());
    int index = 0;
    for(vector<double>::const_iterator it = meaningThicknessElement.begin(); it != meaningThicknessElement.end(); it++)
    {
        double thickness = (*it)*sqrt(2);
        //double thickness = (*it);
        vector<AlphaThickSegmentComputer2D> fuzzySegmentSet = blurredSegmentCurveDecomposition(aContour,thickness,NULL,false);
        cout<<"===> Num of seg decomposed is "<<fuzzySegmentSet.size()<<endl;
        for (vector<AlphaThickSegmentComputer2D>::const_iterator it_bis = fuzzySegmentSet.begin();it_bis != fuzzySegmentSet.end();it_bis++)
            meaningThicknessTangentCover[index].push_back(*it_bis);
        index++;
    }
    if(filename != NULL)
        sprintf(fileAdaptMT,"%s_Step2.svg",filename);
    drawMultiThicknessCover(aContour,meaningThicknessTangentCover,vecMT,fileAdaptMT);//graduate color with vecMT
    ///drawMultiThicknessCover(aContour,meaningThicknessTangentCover,"tmpS2.svg");

    //3. Update thickness of points with tangent covers
    vector<double> vecMTmodified;
    for(vector<double>::const_iterator it = vecMT.begin(); it != vecMT.end(); it++)
        vecMTmodified.push_back(-1);
    int idCover = 0;
    for(vector<vector<AlphaThickSegmentComputer2D> >::const_iterator it = meaningThicknessTangentCover.begin(); it != meaningThicknessTangentCover.end(); it++)
    {
        vector<AlphaThickSegmentComputer2D> fuzzySegmentSet = *it;
        double thickness = meaningThicknessElement.at(idCover);
        int index = 0;
        for (vector<AlphaThickSegmentComputer2D>::const_iterator it_bis = fuzzySegmentSet.begin();it_bis != fuzzySegmentSet.end();it_bis++)
        {
            int idStart = findElement(aContour,getStartPoint(*it_bis),index);
            int idEnd = findElement(aContour,getEndPoint(*it_bis),idStart);
            if(idStart != -1 && idEnd != -1)
            {
                index=idStart;
                double maxThickness = -1;
                for(int i=idStart; i<=idEnd; i++)
                {
                    double thicknessPoint = vecMT.at(i);
                    if(thicknessPoint > maxThickness)
                        maxThickness = thicknessPoint;
                }
                for(int i=idStart; i<=idEnd; i++) // FIXME : idStart + 1 => justify !!!
                {
                    if(vecMTmodified.at(i) < maxThickness && maxThickness==thickness)
                        vecMTmodified.at(i) = maxThickness;
                }
            }
            else
                cout<<"negatif"<<endl;
        }
        idCover++;
    }
    if(filename != NULL)
        sprintf(fileAdaptMT,"%s_Step3.svg",filename);
    drawMultiThicknessCover(aContour,vecMTmodified,meaningThicknessElement.size(),fileAdaptMT);

    //4. Travel over the tangent covers and select the segments w.r.t the associated thickness of points
    vector<vector<AlphaThickSegmentComputer2D> > adaptiveMeaningThicknessTangentCover;
    idCover = 0;
    for(vector<vector<AlphaThickSegmentComputer2D> >::const_iterator it = meaningThicknessTangentCover.begin(); it != meaningThicknessTangentCover.end(); it++)
    {
        adaptiveMeaningThicknessTangentCover.push_back(vector<AlphaThickSegmentComputer2D>());
        vector<AlphaThickSegmentComputer2D> fuzzySegmentSet = *it;
        vector<AlphaThickSegmentComputer2D> adaptiveFuzzySegmentSet;
        int idSeg=0, index=0;
        double thickness = meaningThicknessElement.at(idCover);
        for (vector<AlphaThickSegmentComputer2D>::const_iterator it_bis = fuzzySegmentSet.begin();it_bis != fuzzySegmentSet.end();it_bis++)
        {
            int idStart = findElement(aContour,getStartPoint(*it_bis),index);
            int idEnd = findElement(aContour,getEndPoint(*it_bis),idStart);
            if(idStart != -1 && idEnd != -1)
            {
                index=idStart;
                bool isGoodMT = false, isGoodMTmodif = false;
                for(int i=idStart; i<=idEnd; i++)
                {
                    double thicknessMT = vecMT.at(i); //all elt have same meaningful thickness value (dont contain other meaningful thickness)
                    if(thicknessMT == thickness)
                        isGoodMT = true;
                    double thicknessMTmodif = vecMTmodified.at(i);
                    if(thicknessMTmodif == thickness) //there exist at least one elt in modif having meaningful thickness value
                        isGoodMTmodif = true;
                }
                if(isGoodMTmodif && isGoodMT)
                    adaptiveFuzzySegmentSet.push_back(*it_bis);
            }
            idSeg++;
        }
        for (vector<AlphaThickSegmentComputer2D>::const_iterator it_bis = adaptiveFuzzySegmentSet.begin();it_bis != adaptiveFuzzySegmentSet.end();it_bis++)
            adaptiveMeaningThicknessTangentCover[idCover].push_back(*it_bis);

        idCover++;
    }
    if(filename != NULL)
        sprintf(fileAdaptMT,"%s_Step4.svg",filename);
    drawMultiThicknessCover(aContour,adaptiveMeaningThicknessTangentCover,fileAdaptMT);

    //5. Reorder the multi-thickness tangent cover
    vector<AlphaThickSegmentComputer2D> adaptiveTangentCover;
    int seg=0,nbSeg=0;
    vector<int> idThicknessCover; //stock idSeg of the last seg at idThicknessCover
    for(int it=0; it<meaningThicknessElement.size();it++)
        idThicknessCover.push_back(0);
    for(int it = 0; it < adaptiveMeaningThicknessTangentCover.size(); it++)
        nbSeg += (adaptiveMeaningThicknessTangentCover.at(it)).size();
    cout<<"===> Total num of seg decomposed is "<<nbSeg<<endl;

    while (seg<nbSeg)
    {
        int idMinStart = aContour.size(), idMinEnd = aContour.size();
        int idMin=-1, idSeg=-1;
        //scan adaptiveMeaningThicknessTangentCover
        for(int it = 0; it < adaptiveMeaningThicknessTangentCover.size(); it++)//thickness level = it
        {
            ///vector<Point> extremitySeg;
            //current seg of thickness level idThicknessCover.at(i)
            int idCurrentSeg = idThicknessCover.at(it);
            int index=0;
            if(idCurrentSeg<(adaptiveMeaningThicknessTangentCover.at(it)).size())
            {
                //get idStart and idEnd of seg
                int idCurrentStart = findElement(aContour,getStartPoint((adaptiveMeaningThicknessTangentCover.at(it)).at(idCurrentSeg)),index);
                int idCurrentEnd = findElement(aContour,getEndPoint((adaptiveMeaningThicknessTangentCover.at(it)).at(idCurrentSeg)),idCurrentStart);
                ///extremitySeg.push_back(Point(idStart,idEnd));
                if(idThicknessCover.at(it)<(adaptiveMeaningThicknessTangentCover.at(it)).size())
                {
                    //find min idCurrentStart
                    if(idMinStart==idCurrentStart && idMinEnd<idCurrentEnd)
                    {
                        if(idThicknessCover.at(it)<(adaptiveMeaningThicknessTangentCover.at(it)).size()-1)
                        {
                            idThicknessCover.at(idMin) = idThicknessCover.at(idMin) + 1;
                            seg++;
                        }
                        idSeg = idCurrentSeg;
                        idMin = it;
                        idMinStart = idCurrentStart;
                        idMinEnd = idCurrentEnd;
                        ///extremitySeg.push_back(Point(idCurrentStart,idCurrentEnd));
                    }
                    //if(idMinStart>=idCurrentStart && idMinEnd>=idCurrentEnd)
                    else if(idMinStart>idCurrentStart && idMinEnd>=idCurrentEnd)
                    {
                        idSeg = idCurrentSeg;
                        idMin = it;
                        idMinStart = idCurrentStart;
                        idMinEnd = idCurrentEnd;
                        ///extremitySeg.push_back(Point(idCurrentStart,idCurrentEnd));
                    }
                }
                index=idMinStart;
            }
        }
        adaptiveTangentCover.push_back((adaptiveMeaningThicknessTangentCover.at(idMin)).at(idSeg));
        idThicknessCover.at(idMin) = idThicknessCover.at(idMin) + 1;
        seg++;
    }
    if(filename != NULL)
        sprintf(fileAdaptMT,"%s_Step5.svg",filename);
    drawMultiThicknessCover(aContour,adaptiveTangentCover,fileAdaptMT);

    if(verbose)
    {
        cout<<"==== SEG DECOM ===="<<endl;
        int i = 0,iB=0,iE=0;
        for(vector<AlphaThickSegmentComputer2D>::const_iterator it = adaptiveTangentCover.begin(); it != adaptiveTangentCover.end(); it++)
        {
            iB = findElement(aContour,getStartPoint(*it),iB);
            iE = findElement(aContour,getEndPoint(*it),iB);
            cout<<i<<" : "<<getStartPoint(*it)<<","<<getEndPoint(*it)<<"==> ("<<iB<<") -- "<<"("<<iE<<")"<<endl;
            i++;
        }
        cout<<"==== SEG DECOM of ATC ===="<<endl;
    }

    return adaptiveTangentCover;
}

/****************************************************************/
/************* Adaptive tangent cover computation ***************/
/****************************************************************/
