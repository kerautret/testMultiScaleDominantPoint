#include "testfunctions.h"
#include "myfunctions.h"

/*************************************/
/*** Burred segments decomposition ***/
/*************************************/

/* consider consecutives points */
vector<AlphaThickSegmentComputer2D> blurredSegmentCurveDecompositionV1(const vector<Point>& aContour, double thickness, const char* filename, bool verbose)
{
    int count = 0;
    double mu, nu;
    PointD N;
    vector<AlphaThickSegmentComputer2D> fuzzySegmentSet;
    //run over the points on the contours
    for (vector<Point>::const_iterator it = aContour.begin();it != aContour.end();it++)
    {
        AlphaThickSegmentComputer2D aSegment;
        aSegment.init(it,thickness);
        /* travel over the contour points and add to the seg */
        while (aSegment.end() != aContour.end() && aSegment.extendFront()){}
        if(verbose)
        {
            aSegment.computeParalellStripParams(mu,N,nu);
            cout<<count<<" : ("<<N[0]<<" , "<<N[1]<<" , "<<mu<<" , "<<nu<<" ) ==> ";
            cout<<getStartPoint(aSegment)<<" and "<<getEndPoint(aSegment)<<endl;
        }
        count++;
        fuzzySegmentSet.push_back(aSegment);
        if(aSegment.end() == aContour.end())
            break;
    }

    if(filename != NULL)
    {
        Board2D aBoard;
        /* display the boundary */
        aBoard << SetMode("PointVector", "Both");
        for (vector<Point>::const_iterator it = aContour.begin(); it != aContour.end(); it++)
            aBoard << *it;
        /* display the boundary */

        /* Display boundingbox */
        for (vector<AlphaThickSegmentComputer2D>::const_iterator it = fuzzySegmentSet.begin();it != fuzzySegmentSet.end();it++)
            aBoard << SetMode((*it).className(), "BoundingBox") << *it;
        /* Display boundingbox */
        aBoard.saveSVG(filename);
    }
    return fuzzySegmentSet;
}
/* consider consecutives points */

/* jump over points */
vector<AlphaThickSegmentComputer2D> blurredSegmentCurveDecompositionV2(const vector<Point>& aContour, double thickness, const char* filename, bool verbose)
{
    double mu, nu;
    PointD N;
    vector<AlphaThickSegmentComputer2D> fuzzySegmentSet;
    //run over the points on the contours
    for (vector<Point>::const_iterator it = aContour.begin();it != aContour.end();it++)
    {
        //cout<<*it<<endl;
        AlphaThickSegmentComputer2D aSegment;
        aSegment.init(it,thickness);
        /* travel over the contour points and add to the seg */
        //if(aSegment.end() != aContour.end()) aSegment.extendFront();
        while (aSegment.end() != aContour.end() && aSegment.extendFront()){}
        //if(*it == Point(17,12))
        //    cout<<"==> STOP : getStartPoint = "<<findElement(aContour,getStartPoint(aSegment))<<", getEndPoint = "<<findElement(aContour,getEndPoint(aSegment))<<" and EndBack = "<<findElement(aContour,getEndPoint(fuzzySegmentSet.back()))<<endl;
        if(it == aContour.begin())
        {
            if(verbose)
            {
                aSegment.computeParalellStripParams(mu,N,nu);
                cout<<fuzzySegmentSet.size()<<" with slope = "<<getSlope(N[0],N[1])<<" : ("<<N[0]<<" , "<<N[1]<<" , "<<mu<<" , "<<nu<<" ) ==> ";
                cout<<getStartPoint(aSegment)<<" and "<<getEndPoint(aSegment)<<endl;
            }
            fuzzySegmentSet.push_back(aSegment);
        }
        else if(findElement(aContour,getEndPoint(aSegment)) > findElement(aContour,getEndPoint(fuzzySegmentSet.back())))
        {
            if(verbose)
            {
                aSegment.computeParalellStripParams(mu,N,nu);
                cout<<fuzzySegmentSet.size()<<" with slope = "<<getSlope(N[0],N[1])<<" : ("<<N[0]<<" , "<<N[1]<<" , "<<mu<<" , "<<nu<<" ) ==> ";
                cout<<getStartPoint(aSegment)<<" and "<<getEndPoint(aSegment)<<endl;
            }
            fuzzySegmentSet.push_back(aSegment);
        }
        if(getEndPoint(aSegment) == aContour.back())
            break;
    }

    /****** write plot file **********/
    if(IS_PLOT)
    {
        /* Display points in common zone */
        vector<int> idBegin;
        vector<int> idEnd;
        for (vector<AlphaThickSegmentComputer2D>::const_iterator it = fuzzySegmentSet.begin();it != fuzzySegmentSet.end();it++)
        {
            idBegin.push_back(findElement(aContour,getStartPoint(*it)));
            idEnd.push_back(findElement(aContour,getEndPoint(*it)));
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
        aBoard.saveSVG(filename);
    }
    return fuzzySegmentSet;
}
/* jump over points */

/*************************************/
/*** Burred segments decomposition ***/
/*************************************/


/*********************************/
/****** Curvature estimation *****/
/*********************************/

/* Non symmetry calcul the curvature for every point on the curve */
vector<double> curvatureOnCurveV1(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, int type, const char* filename, bool verbose)
{
    vector<double> curvature;
    int m = fuzzySegmentSet.size(), n = aContour.size(), P_index=0, leftS=0, rightS=0;
    vector<Point> P, P_left, P_right;
    for(vector<Point>::const_iterator it = aContour.begin(); it != aContour.end(); it++)
    {
        //find extremity on the left => min index seg that p belong to
        while(P_index<n && leftS<m && !isBelongTo(P_index,findElement(aContour,getStartPoint(fuzzySegmentSet[leftS])),findElement(aContour,getEndPoint(fuzzySegmentSet[leftS])))) leftS++;
        //find extremity on the right => max index seg that p belong to
        rightS = leftS;
        while(P_index<n && rightS<m && isBelongTo(P_index,findElement(aContour,getStartPoint(fuzzySegmentSet[rightS])),findElement(aContour,getEndPoint(fuzzySegmentSet[rightS])))) rightS++;
        if(rightS > leftS)
        {
            P.push_back(*it);
            Point leftP = getStartPoint(fuzzySegmentSet[leftS]);
            Point rightP = getEndPoint(fuzzySegmentSet[rightS-1]);
            P_left.push_back(leftP);
            P_right.push_back(rightP);
            if(verbose)
                cout<<P_index<<" : "<< *it <<" belong to leftS "<<leftS<<" rightS "<<rightS-1<<endl;
                //cout<<P_index<<" : "<< *it <<" between "<<findElement(aContour,leftP)<<" : "<<leftP<<" and "<<findElement(aContour,rightP)<<" : "<<rightP<<endl;
        }
        P_index++;
    }

    double c;
    for(int it = 0; it< P.size(); it++)
    {
        if(type==1)
            c = getCurvatureAngle(P[it],P_left[it],P_right[it]);
        else if(type==2)
            c = getCurvatureRatio(P[it],P_left[it],P_right[it]);
        else if(type==3)
            c = getCurvatureCosine(P[it],P_left[it],P_right[it]);
        else// if(type==4)
            c = getCurvatureCircle(P[it],P_left[it],P_right[it]);
        curvature.push_back(c);
        if(verbose)
        {
            if(c!=0)
                cout<<"==> curvature at ("<<it<<") "<<P[it]<<" is "<<c<<endl;//<<" and circle radius is "<<fabs(1.0/c)<<endl;
            else
                cout<<"==> curvature at ("<<it<<") "<<P[it]<<" is "<<c<<endl;//<<" and circle radius is INF"<<endl;
        }
    }

    if(filename != NULL)
    {
        Board2D aBoard;
        /* Display the boundary */
        aBoard << SetMode("PointVector", "Both");
        for (vector<Point>::const_iterator it = aContour.begin(); it != aContour.end(); it++)
            aBoard << *it;
        /* Display the boundary */

        /* Display points */
        for(int it = 0; it< P.size(); it++)
        {
            //if(it==4)
            {
                aBoard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,0,142)) )//blue
                       <<P[it];
                aBoard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,192,0)) )//green
                       <<P_left[it]<<P_right[it];
            }
        }
        /* Display points */
        aBoard.saveSVG(filename);
    }

    return curvature;
}
/* Non symmetry calcul the curvature for every point on the curve */

/* Symmetry calcul the curvature for every point on the curve */
vector<double> curvatureOnCurveV2(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, int type, const char* filename, bool verbose)
{
    vector<double> curvature;
    int m = fuzzySegmentSet.size(), n = aContour.size(),P_index=0,leftS=0, rightS=0;
    vector<Point> P, P_left, P_right;
    for(vector<Point>::const_iterator it = aContour.begin(); it != aContour.end(); it++)
    {
        //find extremity on the left => min index seg that p belong to
        while(P_index<n && leftS<m && !isBelongTo(P_index,findElement(aContour,getStartPoint(fuzzySegmentSet[leftS])),findElement(aContour,getEndPoint(fuzzySegmentSet[leftS])))) leftS++;
        //find extremity on the right => max index seg that p belong to
        rightS = leftS;
        while(P_index<n && rightS<m && isBelongTo(P_index,findElement(aContour,getStartPoint(fuzzySegmentSet[rightS])),findElement(aContour,getEndPoint(fuzzySegmentSet[rightS])))) rightS++;
        if(rightS > leftS)
        {
            P.push_back(*it);
            int indexP = findElement(aContour,*it);
            int indexLeft = findElement(aContour,getStartPoint(fuzzySegmentSet[leftS]));
            int indexRight = findElement(aContour,getEndPoint(fuzzySegmentSet[rightS-1]));
            if(fabs(indexP-indexLeft)>fabs(indexRight-indexP))
            {
                P_left.push_back(aContour.at(indexP-fabs(indexRight-indexP)));
                P_right.push_back(aContour.at(indexRight));
                //cout<<" indexP : "<<indexP<<" between "<< indexP-fabs(indexRight-indexP) <<" and "<<indexRight<<endl;
            }
            else
            {
                P_left.push_back(aContour.at(indexLeft));
                P_right.push_back(aContour.at(indexP+fabs(indexP-indexLeft)));
                //cout<<" indexP : "<<indexP<<" between "<< indexLeft <<" and "<<indexP+fabs(indexP-indexLeft)<<endl;
            }
            if(verbose)
                cout<<P_index<<" : "<< *it <<" belong to leftS "<<leftS<<" rightS "<<rightS-1<<endl;
        }
        P_index++;
    }

    double c;
    for(int it = 0; it< P.size(); it++)
    {
        if(type==1)
            c = getCurvatureAngle(P[it],P_left[it],P_right[it]);
        else if(type==2)
            c = getCurvatureRatio(P[it],P_left[it],P_right[it]);
        else if(type==3)
            c = getCurvatureCosine(P[it],P_left[it],P_right[it]);
        else// if(type==4)
            c = getCurvatureCircle(P[it],P_left[it],P_right[it]);
        curvature.push_back(c);
        if(verbose)
        {
            if(c!=0)
                cout<<"==> curvature at ("<<it<<") "<<P[it]<<" is "<<c<<endl;//" and circle radius is "<<fabs(1.0/c)<<endl;
            else
                cout<<"==> curvature at ("<<it<<") "<<P[it]<<" is "<<c<<endl;//" and circle radius is INF"<<endl;
        }
    }

    if(filename != NULL)
    {
        Board2D aBoard;
        /* Display the boundary */
        aBoard << SetMode("PointVector", "Both");
        for (vector<Point>::const_iterator it = aContour.begin(); it != aContour.end(); it++)
            aBoard << *it;
        /* Display the boundary */

        /* Display points */
        for(int it = 0; it< P.size(); it++)
            //if(it==4)
        {
            aBoard << SetMode("PointVector", "Both")
                   << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,0,142)) )//blue
                   <<P[it];
            aBoard << SetMode("PointVector", "Both")
                   << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,192,0)) )//green
                   <<P_left[it]<<P_right[it];
        }
        /* Display points */
        aBoard.saveSVG(filename);
    }

    return curvature;
}
/* Symmetry calcul the curvature for every point on the curve */

/* Direnct neighbours calcul the curvature for every point on the curve */
vector<double> curvatureOnCurveV3(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, int type, const char* filename, bool verbose)
{
    vector<double> curvature;
    vector<Point> P, P_left, P_right;

    P.push_back(aContour.at(0));
    P_left.push_back(aContour.at(0));
    P_right.push_back(aContour.at(1));
    for(vector<Point>::const_iterator it = aContour.begin()+1; it+1 != aContour.end(); it++)
    {
        P.push_back(*it);
        P_left.push_back(*(it-1));
        P_right.push_back(*(it+1));
    }
    P.push_back(aContour.at(aContour.size()-1));
    P_left.push_back(aContour.at(aContour.size()-2));
    P_right.push_back(aContour.at(aContour.size()-1));

    double c;
    for(int it = 0; it< P.size(); it++)
    {
        if(type==1)
            c = getCurvatureAngle(P[it],P_left[it],P_right[it]);
        else if(type==2)
            c = getCurvatureRatio(P[it],P_left[it],P_right[it]);
        else if(type==3)
            c = getCurvatureCosine(P[it],P_left[it],P_right[it]);
        else// if(type==4)
            c = getCurvatureCircle(P[it],P_left[it],P_right[it]);
        curvature.push_back(c);
        if(verbose)
        {
            //cout<<"curvature of point "<<P[it]<<" is "<<P_left[it]<<" and "<<P_right[it]<<endl;
            if(c!=0)
                cout<<"==> curvature at ("<<it<<") "<<P[it]<<" is "<<c<<endl;//<<" and circle radius is "<<fabs(1.0/c)<<endl;
            else
                cout<<"==> curvature at ("<<it<<") "<<P[it]<<" is "<<c<<endl;//" and circle radius is INF"<<endl;
        }
    }

    if(filename != NULL)
    {
        Board2D aBoard;
        /* Display the boundary */
        aBoard << SetMode("PointVector", "Both");
        for (vector<Point>::const_iterator it = aContour.begin(); it != aContour.end(); it++)
            aBoard << *it;
        /* Display the boundary */

        /* Display points */
        for(int it = 0; it< P.size(); it++)
            //if(it==4)
        {
            aBoard << SetMode("PointVector", "Both")
                   << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,0,142)) )//blue
                   <<P[it];
            aBoard << SetMode("PointVector", "Both")
                   << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,192,0)) )//green
                   <<P_left[it]<<P_right[it];
        }
        /* Display points */
        aBoard.saveSVG(filename);
    }

    return curvature;
}
/* Direnct neighbours calcul the curvature for every point on the curve */

/* Means of 3 curvautres calcul the curvature for every point on the curve */
vector<double> curvatureOnCurveCombine(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, int type,const char* filename, bool verbose)
{
    vector<double> curvature1 = curvatureOnCurveV1(fuzzySegmentSet,aContour,type,NULL,false);
    vector<double> curvature2 = curvatureOnCurveV2(fuzzySegmentSet,aContour,type,NULL,false);
    vector<double> curvature3 = curvatureOnCurveV3(fuzzySegmentSet,aContour,type,NULL,false);
    vector<double> curvature;
    for(int it=0; it<aContour.size(); it++ )
        curvature.push_back((curvature1.at(it)+curvature2.at(it)+curvature3.at(it))/3.0);
    return curvature;
}
/* of 3 curvautres calcul the curvature for every point on the curve */

/*********************************/
/****** Curvature estimation *****/
/*********************************/



/**** Count segments passing by each point of the curves ****/
vector<Point> testCountSegmentsOnShape(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, const char* filename, bool verbose)
{    
    vector<int> nbSegPoint(aContour.size(),0);
    vector<int> leftSegPoint(aContour.size(),-1);
    vector<int> rightSegPoint(aContour.size(),-1);
    int idSegCurrent = 0, indexStart=0, indexEnd;
    for(vector<AlphaThickSegmentComputer2D>::const_iterator it = fuzzySegmentSet.begin(); it != fuzzySegmentSet.end(); it++)
    {
        indexStart = findElement(aContour,getStartPoint(*it));
        indexEnd = findElement(aContour,getEndPoint(*it));
        if(indexStart>indexEnd)
            indexEnd += aContour.size();
        //cout<<"idSegCurrent : "<<idSegCurrent<<" indexStart : "<<indexStart<<" indexEnd: "<<indexEnd<<endl;
        for(int i=indexStart; i<=indexEnd; i++)
        {
            //nbSegPoint.at(i%aContour.size()) += 1;
            //if(leftSegPoint.at(i%aContour.size()) == -1) leftSegPoint.at(i%aContour.size()) = idSegCurrent;
            //if(rightSegPoint.at(i%aContour.size())<idSegCurrent) rightSegPoint.at(i%aContour.size()) = idSegCurrent;
            nbSegPoint.at(i) += 1;
            if(leftSegPoint.at(i) == -1) leftSegPoint.at(i) = idSegCurrent;
            if(rightSegPoint.at(i)<idSegCurrent) rightSegPoint.at(i) = idSegCurrent;
        }
        idSegCurrent++;
    }

    vector<Point> maxSegPoint;
    /********** Find points with max seg passing *********
    vector<int> tabZoneB, tabZoneE, tabZoneVal;
    bool isMax = false;
    int idB, idE,val, val_prev, val_succ;
    for(int id = 0; id< nbSegPoint.size()-1;)
    {
        isMax = FALSE;
        idB = idE = id;
        while(isMax == FALSE && id<nbSegPoint.size()-1)
        {
            val = nbSegPoint.at(id);
            val_succ = nbSegPoint.at(id+1);
            if(val < val_succ)
            {
                idB = idE = id+1;
                val_prev = val;
            }
            else if(val > val_succ)
            {
                if(val > val_prev)
                {
                    idE = id;
                    isMax = TRUE;
                    tabZoneB.push_back(idB);
                    tabZoneE.push_back(idE);
                    tabZoneVal.push_back(val);
                }
                val_prev = val;
            }
            else//if(val == val_succ)
                idE = id+1;
            id = id + 1;
        }
    }

    for(int j = 0; j<tabZoneB.size(); j++)
        for(int i = tabZoneB.at(j); i<=tabZoneE.at(j); i++)
            maxSegPoint.push_back(aContour.at(i));
    ********** Find points with max seg passing *********/

    /********** Min angle of points *********
    vector<Point> minAnglePoint;
    for(int j = 0; j<tabZoneB.size(); j++)
    {
        int idMin = tabZoneB.at(j);
        double angle, angleMin = acuteAngle(getStartPoint(fuzzySegmentSet.at(leftSegPoint.at(idMin)%fuzzySegmentSet.size())),
                                     aContour.at(idMin),
                                     getEndPoint(fuzzySegmentSet.at(rightSegPoint.at(idMin)%fuzzySegmentSet.size())));

        for(int i = idMin+1; i<=tabZoneE.at(j); i++)
        {
            angle = acuteAngle(getStartPoint(fuzzySegmentSet.at(leftSegPoint.at(i)%fuzzySegmentSet.size())),
                               aContour.at(i),
                               getEndPoint(fuzzySegmentSet.at(rightSegPoint.at(i)%fuzzySegmentSet.size())));
            if(angle<angleMin)
            {
                idMin = i;
                angleMin = angle;
            }
        }
        minAnglePoint.push_back(aContour.at(idMin));
    }
    ********** Min angle of points *********/

    /******** write plot file ********
    if(IS_PLOT)// && !isClosed
    {
        vector<double> anglePoint;
        vector<double> radiusPoint;
        for(int i=0; i< aContour.size(); i++)
        {
            Point pLeft = getStartPoint(fuzzySegmentSet.at(leftSegPoint.at(i)%fuzzySegmentSet.size()));
            Point pRight = getEndPoint(fuzzySegmentSet.at(rightSegPoint.at(i)%fuzzySegmentSet.size()));
            if(verbose)
                cout<<"Point "<<i<<" has "<<nbSegPoint.at(i)<<
                      " seg passing, with left seg "<<leftSegPoint.at(i)<<
                      " and right seg "<<rightSegPoint.at(i);
            //if(leftSegPoint.at(i%aContour.size()) == rightSegPoint.at(i))
            if(leftSegPoint.at(i) == rightSegPoint.at(i))
            {
                if(verbose)
                    cout<<" => angle = 360"<<endl;
                anglePoint.push_back(2*M_PI);
            }
            else
            {
                //double angle = relativeAngle(pLeft, aContour.at(i%aContour.size()),pRight);
                //double radius = 1.0/getCurvatureCircle(aContour.at(i%aContour.size()),pLeft,pRight);
                double angle = relativeAngle(pLeft, aContour.at(i),pRight);
                double radius = 1.0/getCurvatureCircle(aContour.at(i),pLeft,pRight);

                if(verbose)
                    cout<<" => angle = "<<angle*180/M_PI<<endl;
                anglePoint.push_back(angle);
                radiusPoint.push_back(radius);
            }
        }
        writeFile(nbSegPoint,"nbSeg.txt",false);
        writeFile(anglePoint,"anglePoint.txt",false);
        writeFile(radiusPoint,"radiusPoint.txt",false);
    }
    /******** write plot file ********/

    if(filename != NULL)
    {
        Board2D aBoard;
        int max = *std::max_element(nbSegPoint.begin(),nbSegPoint.end());
        int step_color = 255/(max);
        /* Display the boundary point in function of its nb seg passing */
        aBoard << SetMode("PointVector", "Both");
        int idPoint = 0;
        for (vector<Point>::const_iterator it = aContour.begin(); it != aContour.end(); it++)
        {
            aBoard << SetMode((*it).className(), "Both")
                   << CustomStyle( (*it).className(), new CustomColors( Color(0,255,0), Color(0,255-(nbSegPoint.at(idPoint)*step_color),0)) )
                   <<(*it);
            idPoint++;
        }
        /* Display the boundary point in function of its nb seg passing */
        /* Display points with max nb seg passing *
        aBoard << SetMode("PointVector", "Both");
        for (vector<Point>::const_iterator it = maxSegPoint.begin(); it != maxSegPoint.end(); it++)
        {
            aBoard << SetMode((*it).className(), "Both")
                   << CustomStyle( (*it).className(), new CustomColors( Color(0,0,255), Color(0,0,194)) )
                   <<(*it);
        }
        * Display points with max nb seg passing */
        /* Display points with min angle *
        aBoard << SetMode("PointVector", "Both");
        for (vector<Point>::const_iterator it = minAnglePoint.begin(); it != minAnglePoint.end(); it++)
        {
            aBoard << SetMode((*it).className(), "Both")
                   << CustomStyle( (*it).className(), new CustomColors( Color(255,0,0), Color(194,0,0)) )
                   <<(*it);
        }
        * Display points with min angle */
        aBoard.saveSVG(filename);
    }

    return maxSegPoint;
}
/**** Count segments passing by each point of the curves ****/




/************************************/
/**** Dominant points detections ****/
/************************************/

/* Phuong's version q = p - 1 + mid points */
vector<Point> testDominantPointOnShapeV1(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, bool isClosed, const char* filename, bool verbose)
{    
    /* pile construction */
    int p=1, q=0, m=fuzzySegmentSet.size();
    vector<Point> pile;
    int Eq=0, Bp=0;
    while (p<m && Bp>=0 && Eq>=0)
    {
        Bp = findElement(aContour,getStartPoint(fuzzySegmentSet[p]));
        Eq = findElement(aContour,getEndPoint(fuzzySegmentSet[q]));
        while(Eq>=Bp && p<m && q<m && Bp>=0 && Eq>=0 )
        {
            Bp = findElement(aContour,getStartPoint(fuzzySegmentSet[p]));
            Eq = findElement(aContour,getEndPoint(fuzzySegmentSet[q]));

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
    for(vector<Point>::const_iterator it = pile.begin(); it != pile.end(); it++)
    {
        q = (*it)[0];
        p = (*it)[1];
        if(verbose)
        {
            cout<<endl<<"slope of lines between "<<q<<" and "<<p <<" : "<<endl;
            for(int i = q; i <= p; i++)
            {
                double mu, nu;
                PointD N;
                fuzzySegmentSet[i].computeParalellStripParams(mu,N,nu);
                //normal vector N(u,v), then slope m = -u/v
                cout<<"slope : "<<getSlope(N[0],N[1])<<" and angle is "<<(atan(-N[0]/N[1])*180.0)/M_PI<<endl;
            }
        }

        double mu1, nu1, mu2, nu2;
        PointD N1,N2;
        int i=p, r;
        fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
        fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);
        double slope1 = getSlope(N1[0],N1[1]);
        double slope2 = getSlope(N2[0],N2[1]);
        while(i > q && slope1 == slope2)
        {
            fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
            fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);
            slope1 = getSlope(N1[0],N1[1]);//(N1[1] == 0.0 ? 0.0 : -N1[0]/N1[1]);
            slope2 = getSlope(N2[0],N2[1]);//(N2[1] == 0.0 ? 0.0 : -N2[0]/N2[1]);
            i--;
        }
        if(slope1<slope2)//increasing order
        {
            while(i > q && slope1 <= slope2)
            {
                fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
                fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);
                slope1 = getSlope(N1[0],N1[1]);
                slope2 = getSlope(N2[0],N2[1]);
                i--;
            }
            if(slope1 <= slope2) r = i;
            else r = i+1;
        }
        else //if(slope1>slope2)//decreasing order
        {
            while(i > q && slope1 >= slope2)
            {
                fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
                fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);
                slope1 = getSlope(N1[0],N1[1]);
                slope2 = getSlope(N2[0],N2[1]);
                i--;
            }
            if(slope1 > slope2) r = i;
            else r = i+1;
        }
        if(verbose)
            cout<<"==> monoton seq is ("<<r<<" , "<<p<<")"<<endl;

        //find mid point
        int indexB = findElement(aContour,getStartPoint(fuzzySegmentSet[p]));
        int indexE = findElement(aContour,getEndPoint(fuzzySegmentSet[r]));
        int indexC = indexB + (indexE-indexB)/2;
        if(verbose)
            cout<<" ==> center point is "<<aContour[indexC]<<endl;

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
        for(vector<Point>::const_iterator it = commonPoint.begin(); it != commonPoint.end(); it++)
            idCommonPoint.push_back(findElement(aContour,*it));
        vector<int> idDominantPoint;
        for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
            idDominantPoint.push_back(findElement(aContour,*it));
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
        aBoard.drawLine(aContour.front()[0],aContour.front()[1],DP.front()[0],DP.front()[1]);
        for(vector<Point>::const_iterator it = DP.begin(); it+1 != DP.end(); it++)
            aBoard.drawLine((*it)[0],(*it)[1],(*(it+1))[0],(*(it+1))[1]);
        if(isClosed)//NOTE : DP.front() = aContour.front()
            aBoard.drawLine(aContour.front()[0],aContour.front()[1],DP.back()[0],DP.back()[1]);
        /* Display the segments by DP */
        aBoard.saveSVG(filename);
    }
    return DP;
}
/* Phuong's version q = p - 1 + mid points */

/* Phuong's version q = p - 1 + curvature */
vector<Point> testDominantPointOnShapeV1(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, const vector<double>& curvature, bool isClosed, const char* filename, bool verbose)
{
    /* pile construction */
    int p=1, q=0, m=fuzzySegmentSet.size();
    vector<Point> pile;
    int Eq=0, Bp=0;
    while (p<m && Bp>=0 && Eq>=0)
    {
        Bp = findElement(aContour,getStartPoint(fuzzySegmentSet[p]));
        Eq = findElement(aContour,getEndPoint(fuzzySegmentSet[q]));
        while(Eq>=Bp && p<m && q<m && Bp>=0 && Eq>=0 )
        {
            Bp = findElement(aContour,getStartPoint(fuzzySegmentSet[p]));
            Eq = findElement(aContour,getEndPoint(fuzzySegmentSet[q]));

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
    for(vector<Point>::const_iterator it = pile.begin(); it != pile.end(); it++)
    {
        q = (*it)[0];
        p = (*it)[1];
        if(verbose)
        {
            cout<<endl<<"slope of lines between "<<q<<" and "<<p <<" : "<<endl;
            for(int i = q; i <= p; i++)
            {
                double mu, nu;
                PointD N;
                fuzzySegmentSet[i].computeParalellStripParams(mu,N,nu);
                //normal vector N(u,v), then slope m = -u/v
                cout<<"slope : "<<getSlope(N[0],N[1])<<" and angle is "<<(atan(-N[0]/N[1])*180.0)/M_PI<<endl;
            }
        }

        double mu1, nu1, mu2, nu2;
        PointD N1,N2;
        int i=p, r;
        fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
        fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);
        double slope1 = getSlope(N1[0],N1[1]);
        double slope2 = getSlope(N2[0],N2[1]);
        while(i > q && slope1 == slope2)
        {
            fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
            fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);
            slope1 = getSlope(N1[0],N1[1]);//(N1[1] == 0.0 ? 0.0 : -N1[0]/N1[1]);
            slope2 = getSlope(N2[0],N2[1]);//(N2[1] == 0.0 ? 0.0 : -N2[0]/N2[1]);
            i--;
        }
        if(slope1<slope2)//increasing order
        {
            while(i > q && slope1 <= slope2)
            {
                fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
                fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);
                slope1 = getSlope(N1[0],N1[1]);
                slope2 = getSlope(N2[0],N2[1]);
                i--;
            }
            if(slope1 <= slope2) r = i;
            else r = i+1;
        }
        else //if(slope1>slope2)//decreasing order
        {
            while(i > q && slope1 >= slope2)
            {
                fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
                fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);
                slope1 = getSlope(N1[0],N1[1]);
                slope2 = getSlope(N2[0],N2[1]);
                i--;
            }
            if(slope1 >= slope2) r = i;
            else r = i+1;
        }
        if(verbose)
            cout<<"==> monoton seq is ("<<r<<" , "<<p<<")"<<endl;

        //find point with max (local) curvature
        int indexB = findElement(aContour,getStartPoint(fuzzySegmentSet[p]));
        int indexE = findElement(aContour,getEndPoint(fuzzySegmentSet[r]));
        int indexC = indexB;
        if(verbose)
            cout<<"curvature at "<<indexB<<" is "<<curvature[indexB]<<endl;
        for(int i = indexB+1; i<=indexE && i<curvature.size() ; i++)
        {
            if(fabs(curvature[i]) >= fabs(curvature[indexC]))
                indexC = i;
            if(verbose)
                cout<<"curvature at "<<i<<" is "<<curvature[i]<<endl;
        }
        if(verbose)
            cout<<" ==> max local curvature is "<<curvature[indexC]<<" at "<<aContour[indexC]<<endl;
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
        for(vector<Point>::const_iterator it = commonPoint.begin(); it != commonPoint.end(); it++)
            idCommonPoint.push_back(findElement(aContour,*it));
        vector<int> idDominantPoint;
        for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
            idDominantPoint.push_back(findElement(aContour,*it));
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
        aBoard.drawLine(aContour.front()[0],aContour.front()[1],DP.front()[0],DP.front()[1]);
        for(vector<Point>::const_iterator it = DP.begin(); it+1 != DP.end(); it++)
            aBoard.drawLine((*it)[0],(*it)[1],(*(it+1))[0],(*(it+1))[1]);
        if(isClosed)//NOTE : DP.front() = aContour.front()
                aBoard.drawLine(aContour.front()[0],aContour.front()[1],DP.back()[0],DP.back()[1]);
        /* Display the segments by DP */
        aBoard.saveSVG(filename);
    }
    return DP;
}
/* Phuong's version q = p - 1 + curvature */

/* q = last monotone sequence */
vector<Point> testDominantPointOnShapeV2(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, bool isClosed, const char* filename, bool verbose)
{
    /* pile construction */
    int p=1, q=0, m=fuzzySegmentSet.size();
    vector<Point> pile;
    int Eq=0, Bp=0;
    double mu, nu, mu1, nu1, mu2, nu2,slopeBq,slopeBp,slopeAq,slopeAp;
    PointD N, N1,N2;
    bool doChange = false, isIncreasingB = false, isIncreasingA = false; //flag of increasing/decreasing monton sequence of segment
    while (p<m && Bp>=0 && Eq>=0)
    {
        doChange = false;
        do
        {
            Bp = findElement(aContour,getStartPoint(fuzzySegmentSet[p]));
            Eq = findElement(aContour,getEndPoint(fuzzySegmentSet[q]));
            fuzzySegmentSet[q].computeParalellStripParams(mu1,N1,nu1);
            fuzzySegmentSet[p].computeParalellStripParams(mu2,N2,nu2);
            slopeBq = getSlope(N1[0],N1[1]);
            slopeBp = getSlope(N2[0],N2[1]);
            p++;
        } while (slopeBq==slopeBp && Eq>Bp && p<m);
        p--;
        Bp = findElement(aContour,getStartPoint(fuzzySegmentSet[p]));
        Eq = findElement(aContour,getEndPoint(fuzzySegmentSet[q]));
        fuzzySegmentSet[q].computeParalellStripParams(mu1,N1,nu1);
        fuzzySegmentSet[p].computeParalellStripParams(mu2,N2,nu2);
        slopeBq = getSlope(N1[0],N1[1]);
        slopeBp = getSlope(N2[0],N2[1]);
        if(slopeBq<slopeBp) isIncreasingB = true;
        else isIncreasingB = false;
        if(verbose)
            cout<<"Before, isIncreasing "<<isIncreasingB<<endl;
        while(Eq>=Bp && p<m && q<m && Bp>=0 && Eq>=0 )
        {
            Bp = findElement(aContour,getStartPoint(fuzzySegmentSet[p]));
            Eq = findElement(aContour,getEndPoint(fuzzySegmentSet[q]));
            if(verbose && Eq>=Bp)
                cout<<"Pair segment : ( "<<q<<","<<p<<" ) and Pair points : ( "<<Bp<<","<<Eq<<" )"<<endl;
            if(p>q+1 && Eq>=Bp)
            {
                fuzzySegmentSet[p-1].computeParalellStripParams(mu1,N1,nu1);//q
                fuzzySegmentSet[p].computeParalellStripParams(mu2,N2,nu2);
                slopeAq = getSlope(N1[0],N1[1]);
                slopeAp = getSlope(N2[0],N2[1]);
                if(slopeAq<=slopeAp) isIncreasingA = true;
                else isIncreasingA = false;
                if(verbose)
                    cout<<"After, isIncreasing "<<isIncreasingA<<endl;
                if(isIncreasingA != isIncreasingB)
                {
                    if(!doChange)
                    {
                        if(verbose)
                            cout<<"==> monoton sequence change ! q restart at "<<q+1<<endl;
                        if(p-2 == q)
                            isIncreasingB = isIncreasingA;
                        q = q+1;
                        p = q;
                        doChange = true;
                        if(verbose)
                            cout<<"Before, isIncreasing "<<isIncreasingB<<endl;
                    }
                    else
                        break;
                }
            }
            p++;
        }
        if(Eq<Bp)
        {
            p--;//when out of loop => Eq<Bp
            //Point tmp = Point(q,p-1);
            //if(pile.empty() || tmp != pile.back())
            pile.push_back(Point(q,p-1));
            if(verbose)
                cout<<"Pile : ("<<q<<" , "<<p-1<<" ) "<<endl<<endl;
            q=p-1;
        }
        else//if(Eq>=Bp)
        {
            if(doChange)
            {
                //Point tmp = Point(q,p-1);
                //if(pile.empty() || tmp != pile.back())
                pile.push_back(Point(q,p-1));
                if(verbose)
                    cout<<"Pile : ("<<q<<" , "<<p-1<<" ) "<<endl<<endl;
                q=p-1;
            }
            else
            {
                if(p>=m || q>=m)
                {
                    Point tmp = Point(q,p-1);
                    if(pile.empty() || tmp != pile.back())
                        pile.push_back(Point(q,p-1));
                    if(verbose)
                        cout<<"Pile : ("<<q<<" , "<<p-1<<" ) "<<endl<<endl;
                }
                else
                    q = p;
            }
        }
    }
    /* pile construction */

    /* calcul point dominant */
    vector<Point> DP, extremityPoint, commonPoint;
    for(vector<Point>::const_iterator it = pile.begin(); it != pile.end(); it++)
    {
        q = (*it)[0];
        p = (*it)[1];
        if(verbose)
        {
            cout<<endl<<"slope of lines between "<<q<<" and "<<p <<" : "<<endl;
            for(int i = q; i <= p; i++)
            {
                fuzzySegmentSet[i].computeParalellStripParams(mu,N,nu);
                //normal vector N(u,v), then slope m = -u/v
                cout<<"slope : "<<getSlope(N[0],N[1])<<" and angle is "<<(atan(-N[0]/N[1])*180.0)/M_PI<<endl;
            }
        }
        int i=p, r;
        fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
        fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);
        double slope1 = getSlope(N1[0],N1[1]);
        double slope2 = getSlope(N2[0],N2[1]);
        while(i > q && slope1 == slope2)
        {
            fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
            fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);
            slope1 = getSlope(N1[0],N1[1]);//(N1[1] == 0.0 ? 0.0 : -N1[0]/N1[1]);
            slope2 = getSlope(N2[0],N2[1]);//(N2[1] == 0.0 ? 0.0 : -N2[0]/N2[1]);
            i--;
        }
        if(slope1<slope2)//increasing order
        {
            while(i > q && slope1 <= slope2)
            {
                fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
                fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);
                slope1 = getSlope(N1[0],N1[1]);
                slope2 = getSlope(N2[0],N2[1]);
                i--;
            }
            if(slope1 <= slope2) r = i;
            else r = i+1;
        }
        else //if(slope1>slope2)//decreasing order
        {
            while(i > q && slope1 >= slope2)
            {
                fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
                fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);
                slope1 = getSlope(N1[0],N1[1]);
                slope2 = getSlope(N2[0],N2[1]);
                i--;
            }
            if(slope1 >= slope2) r = i;
            else r = i+1;
        }
        if(verbose)
        {
            cout<<"monoton seq is ("<<r<<" , "<<p<<")"<<endl;
        }
        if(r!=q)
        {
            //cout<<"ERROR at (q = "<<q<<" , p = "<<p<<") "<<endl;
            continue;
        }
        //find point with max (local) curvature
        int indexB = findElement(aContour,getStartPoint(fuzzySegmentSet[p]));
        int indexE = findElement(aContour,getEndPoint(fuzzySegmentSet[r]));
        int indexC = (indexB+indexE)/2;
        if(verbose)
            cout<<" ==> center point is "<<aContour[indexC]<<endl;
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
        for(vector<Point>::const_iterator it = commonPoint.begin(); it != commonPoint.end(); it++)
            idCommonPoint.push_back(findElement(aContour,*it));
        vector<int> idDominantPoint;
        for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
            idDominantPoint.push_back(findElement(aContour,*it));
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
        aBoard.drawLine(aContour.front()[0],aContour.front()[1],DP.front()[0],DP.front()[1]);
        for(vector<Point>::const_iterator it = DP.begin(); it+1 != DP.end(); it++)
            aBoard.drawLine((*it)[0],(*it)[1],(*(it+1))[0],(*(it+1))[1]);
        if(isClosed)//NOTE : DP.front() = aContour.front()
            aBoard.drawLine(aContour.front()[0],aContour.front()[1],DP.back()[0],DP.back()[1]);
        /* Display the segments by DP */
        aBoard.saveSVG(filename);
    }
    return DP;
}

vector<Point> testDominantPointOnShapeV2(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, const vector<double>& curvature, bool isClosed, const char* filename, bool verbose)
{
    /* pile construction */
    int p=1, q=0, m=fuzzySegmentSet.size();
    vector<Point> pile;
    int Eq=0, Bp=0;
    double mu, nu, mu1, nu1, mu2, nu2,slopeBq,slopeBp,slopeAq,slopeAp;
    PointD N, N1,N2;
    bool doChange = false, isIncreasingB = false, isIncreasingA = false; //flag of increasing/decreasing monton sequence of segment
    while (p<m && Bp>=0 && Eq>=0)
    {
        doChange = false;
        do
        {
            Bp = findElement(aContour,getStartPoint(fuzzySegmentSet[p]));
            Eq = findElement(aContour,getEndPoint(fuzzySegmentSet[q]));
            fuzzySegmentSet[q].computeParalellStripParams(mu1,N1,nu1);
            fuzzySegmentSet[p].computeParalellStripParams(mu2,N2,nu2);
            slopeBq = getSlope(N1[0],N1[1]);
            slopeBp = getSlope(N2[0],N2[1]);
            p++;
        } while (slopeBq==slopeBp && Eq>Bp && p<m);
        p--;
        Bp = findElement(aContour,getStartPoint(fuzzySegmentSet[p]));
        Eq = findElement(aContour,getEndPoint(fuzzySegmentSet[q]));
        fuzzySegmentSet[q].computeParalellStripParams(mu1,N1,nu1);
        fuzzySegmentSet[p].computeParalellStripParams(mu2,N2,nu2);
        slopeBq = getSlope(N1[0],N1[1]);
        slopeBp = getSlope(N2[0],N2[1]);
        if(slopeBq<slopeBp) isIncreasingB = true;
        else isIncreasingB = false;
        if(verbose)
            cout<<"Before, isIncreasing "<<isIncreasingB<<endl;
        while(Eq>=Bp && p<m && q<m && Bp>=0 && Eq>=0 )
        {
            Bp = findElement(aContour,getStartPoint(fuzzySegmentSet[p]));
            Eq = findElement(aContour,getEndPoint(fuzzySegmentSet[q]));
            if(verbose && Eq>=Bp)
                cout<<"Pair segment : ( "<<q<<","<<p<<" ) and Pair points : ( "<<Bp<<","<<Eq<<" )"<<endl;
            if(p>q+1 && Eq>=Bp)
            {
                fuzzySegmentSet[p-1].computeParalellStripParams(mu1,N1,nu1);//q
                fuzzySegmentSet[p].computeParalellStripParams(mu2,N2,nu2);
                slopeAq = getSlope(N1[0],N1[1]);
                slopeAp = getSlope(N2[0],N2[1]);
                if(slopeAq<=slopeAp) isIncreasingA = true;
                else isIncreasingA = false;
                if(verbose)
                    cout<<"After, isIncreasing "<<isIncreasingA<<endl;
                if(isIncreasingA != isIncreasingB)
                {
                    if(!doChange)
                    {
                        if(verbose)
                            cout<<"==> monoton sequence change ! q restart at "<<q+1<<endl;
                        if(p-2 == q)
                            isIncreasingB = isIncreasingA;
                        q = q+1;
                        p = q;
                        doChange = true;
                        if(verbose)
                            cout<<"Before, isIncreasing "<<isIncreasingB<<endl;
                    }
                    else
                        break;
                }
            }
            p++;
        }
        if(Eq<Bp)
        {
            p--;//when out of loop => Eq<Bp
            //Point tmp = Point(q,p-1);
            //if(pile.empty() || tmp != pile.back())
            pile.push_back(Point(q,p-1));
            if(verbose)
                cout<<"Pile : ("<<q<<" , "<<p-1<<" ) "<<endl<<endl;
            q=p-1;
        }
        else//if(Eq>=Bp)
        {
            if(doChange)
            {
                //Point tmp = Point(q,p-1);
                //if(pile.empty() || tmp != pile.back())
                pile.push_back(Point(q,p-1));
                if(verbose)
                    cout<<"Pile : ("<<q<<" , "<<p-1<<" ) "<<endl<<endl;
                q=p-1;
            }
            else
            {
                if(p>=m || q>=m)
                {
                    Point tmp = Point(q,p-1);
                    if(pile.empty() || tmp != pile.back())
                        pile.push_back(Point(q,p-1));
                    if(verbose)
                        cout<<"Pile : ("<<q<<" , "<<p-1<<" ) "<<endl<<endl;
                }
                else
                    q = p;
            }
        }
    }
    /* pile construction */

    /* calcul point dominant */
    vector<Point> DP, extremityPoint, commonPoint;
    for(vector<Point>::const_iterator it = pile.begin(); it != pile.end(); it++)
    {
        q = (*it)[0];
        p = (*it)[1];
        if(verbose)
        {
            cout<<endl<<"slope of lines between "<<q<<" and "<<p <<" : "<<endl;
            for(int i = q; i <= p; i++)
            {
                fuzzySegmentSet[i].computeParalellStripParams(mu,N,nu);
                //normal vector N(u,v), then slope m = -u/v
                cout<<"slope : "<<getSlope(N[0],N[1])<<" and angle is "<<(atan(-N[0]/N[1])*180.0)/M_PI<<endl;
            }
        }
        int i=p, r;
        fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
        fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);
        double slope1 = getSlope(N1[0],N1[1]);
        double slope2 = getSlope(N2[0],N2[1]);
        while(i > q && slope1 == slope2)
        {
            fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
            fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);
            slope1 = getSlope(N1[0],N1[1]);//(N1[1] == 0.0 ? 0.0 : -N1[0]/N1[1]);
            slope2 = getSlope(N2[0],N2[1]);//(N2[1] == 0.0 ? 0.0 : -N2[0]/N2[1]);
            i--;
        }
        if(slope1<slope2)//increasing order
        {
            while(i > q && slope1 <= slope2)
            {
                fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
                fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);
                slope1 = getSlope(N1[0],N1[1]);
                slope2 = getSlope(N2[0],N2[1]);
                i--;
            }
            if(slope1 <= slope2) r = i;
            else r = i+1;
        }
        else //if(slope1>slope2)//decreasing order
        {
            while(i > q && slope1 >= slope2)
            {
                fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
                fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);
                slope1 = getSlope(N1[0],N1[1]);
                slope2 = getSlope(N2[0],N2[1]);
                i--;
            }
            if(slope1 >= slope2) r = i;
            else r = i+1;
        }
        if(verbose)
        {
            cout<<"monoton seq is ("<<r<<" , "<<p<<")"<<endl;
        }
        if(r!=q)
        {
            //cout<<"ERROR at (q = "<<q<<" , p = "<<p<<") "<<endl;
            continue;
        }
        //find point with max (local) curvature
        int indexB = findElement(aContour,getStartPoint(fuzzySegmentSet[p]));
        int indexE = findElement(aContour,getEndPoint(fuzzySegmentSet[r]));
        int indexC = indexB;
        if(verbose)
            cout<<"curvature at "<<indexB<<" is "<<curvature[indexB]<<endl;

        for(int i = indexB+1; i<=indexE && i<curvature.size() ; i++)
        {
            if(fabs(curvature[i]) >= fabs(curvature[indexC]))
                indexC = i;
            if(verbose)
                cout<<"curvature at "<<i<<" is "<<curvature[i]<<endl;
        }
        if(verbose)
            cout<<" has max local curvature is "<<curvature[indexC]<<" at "<<aContour[indexC]<<endl;

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
        for(vector<Point>::const_iterator it = commonPoint.begin(); it != commonPoint.end(); it++)
            idCommonPoint.push_back(findElement(aContour,*it));
        vector<int> idDominantPoint;
        for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
            idDominantPoint.push_back(findElement(aContour,*it));
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
        aBoard.drawLine(aContour.front()[0],aContour.front()[1],DP.front()[0],DP.front()[1]);
        for(vector<Point>::const_iterator it = DP.begin(); it+1 != DP.end(); it++)
            aBoard.drawLine((*it)[0],(*it)[1],(*(it+1))[0],(*(it+1))[1]);
        if(isClosed)//NOTE : DP.front() = aContour.front()
            aBoard.drawLine(aContour.front()[0],aContour.front()[1],DP.back()[0],DP.back()[1]);
        /* Display the segments by DP */
        aBoard.saveSVG(filename);
    }
    return DP;
}
/* q = last monotone sequence */

/* q = p-1 with ordering ANGLE of sequence of blurred segments */
vector<Point> testDominantPointOnShapeV3(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, bool isClosed, const char* filename, bool verbose)
{
    int p=1, q=0, Eq=0, Bp=0, m=fuzzySegmentSet.size();
    vector<Point> pile;
    while (p<m && q<m && Bp>=0 && Eq>=0)
    {
        Bp = findStartElement(aContour,fuzzySegmentSet[p]);
        Eq = findEndElement(aContour,fuzzySegmentSet[q]);
        while(Eq>=Bp && p<m && q<m && Bp>=0 && Eq>=0 )
        {
            Bp = findStartElement(aContour,fuzzySegmentSet[p]);
            Eq = findEndElement(aContour,fuzzySegmentSet[q]);
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
    double mu, nu, mu1, nu1, mu2, nu2;
    PointD N, N1, N2;
    Point startS, startS1, startS2, endS, endS1, endS2;
    for(vector<Point>::const_iterator it = pile.begin(); it != pile.end(); it++)
    {
        q = (*it)[0];
        p = (*it)[1];
        if(verbose)
        {
            cout<<endl<<"slope of lines between "<<q<<" and "<<p <<" : "<<endl;
            for(int i = q; i <= p; i++)
            {
                fuzzySegmentSet[i].computeParalellStripParams(mu,N,nu);
                //normal vector N(u,v), then slope m = -u/v
                startS = getStartPoint(fuzzySegmentSet[i]);
                endS = getEndPoint(fuzzySegmentSet[i]);
                N = directionOfVector(N,startS,endS);
                cout<<"slope : "<<getSlope(N[0],N[1])<<" and angle is "<<(atan(-N[0]/N[1])*180.0)/M_PI<<endl;
            }
        }
        int i=p, r;
        fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
        fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);

        startS1 = getStartPoint(fuzzySegmentSet[i]);
        endS1 = getEndPoint(fuzzySegmentSet[i]);
        startS2 = getStartPoint(fuzzySegmentSet[i-1]);
        endS2 = getEndPoint(fuzzySegmentSet[i-1]);
        N1 = directionOfVector(N1,startS1,endS1);
        N2 = directionOfVector(N2,startS2,endS2);

        bool isIncrease = isIncreasingOrder(N1,N2);
        bool isIncrease_prev = isIncrease;

        while(i > q && (isIncrease_prev == isIncrease))
        {
            fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
            fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);

            startS1 = getStartPoint(fuzzySegmentSet[i]);
            endS1 = getEndPoint(fuzzySegmentSet[i]);
            startS2 = getStartPoint(fuzzySegmentSet[i-1]);
            endS2 = getEndPoint(fuzzySegmentSet[i-1]);
            N1 = directionOfVector(N1,startS1,endS1);
            N2 = directionOfVector(N2,startS2,endS2);

            isIncrease_prev = isIncrease;
            isIncrease = isIncreasingOrder(N1,N2);
            i--;
        }
        if(isIncrease_prev == isIncrease) r = i;
        else r = i+1;
        if(verbose)
            cout<<"==> monoton seq is ("<<r<<" , "<<p<<")"<<endl;

        ///r = q;
        //find mid point
        int indexB = findStartElement(aContour,fuzzySegmentSet[p]);
        int indexE = findEndElement(aContour,fuzzySegmentSet[r]);
        int indexC = (indexE+indexB)/2;
        if(verbose)
            cout<<" ==> center point is "<<aContour[indexC]<<endl;

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
        for(vector<Point>::const_iterator it = commonPoint.begin(); it != commonPoint.end(); it++)
            idCommonPoint.push_back(findElement(aContour,*it));
        vector<int> idDominantPoint;
        for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
            idDominantPoint.push_back(findElement(aContour,*it));
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
        aBoard.saveSVG(filename);
    }
    return DP;
}

vector<Point> testDominantPointOnShapeV3(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, const vector<double>& curvature, bool isClosed, const char* filename, bool verbose)
{
    int p=1, q=0, Eq=0, Bp=0, m=fuzzySegmentSet.size();
    vector<Point> pile;
    while (p<m && q<m && Bp>=0 && Eq>=0)
    {
        Bp = findStartElement(aContour,fuzzySegmentSet[p]);
        Eq = findEndElement(aContour,fuzzySegmentSet[q]);
        while(Eq>=Bp && p<m && q<m && Bp>=0 && Eq>=0 )
        {
            Bp = findStartElement(aContour,fuzzySegmentSet[p]);
            Eq = findEndElement(aContour,fuzzySegmentSet[q]);
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
    double mu, nu, mu1, nu1, mu2, nu2;
    PointD N, N1, N2;
    Point startS, startS1, startS2, endS, endS1, endS2;
    for(vector<Point>::const_iterator it = pile.begin(); it != pile.end(); it++)
    {
        q = (*it)[0];
        p = (*it)[1];
        if(verbose)
        {
            cout<<endl<<"slope of lines between "<<q<<" and "<<p <<" : "<<endl;
            for(int i = q; i <= p; i++)
            {
                fuzzySegmentSet[i].computeParalellStripParams(mu,N,nu);
                //normal vector N(u,v), then slope m = -u/v
                startS = getStartPoint(fuzzySegmentSet[i]);
                endS = getEndPoint(fuzzySegmentSet[i]);
                N = directionOfVector(N,startS,endS);
                cout<<"slope : "<<getSlope(N[0],N[1])<<" and angle is "<<(atan(-N[0]/N[1])*180.0)/M_PI<<endl;
            }
        }
        int i=p, r;
        fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
        fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);

        startS1 = getStartPoint(fuzzySegmentSet[i]);
        endS1 = getEndPoint(fuzzySegmentSet[i]);
        startS2 = getStartPoint(fuzzySegmentSet[i-1]);
        endS2 = getEndPoint(fuzzySegmentSet[i-1]);
        N1 = directionOfVector(N1,startS1,endS1);
        N2 = directionOfVector(N2,startS2,endS2);

        bool isIncrease = isIncreasingOrder(N1,N2);
        bool isIncrease_prev = isIncrease;

        while(i > q && (isIncrease_prev == isIncrease))
        {
            fuzzySegmentSet[i].computeParalellStripParams(mu1,N1,nu1);
            fuzzySegmentSet[i-1].computeParalellStripParams(mu2,N2,nu2);

            startS1 = getStartPoint(fuzzySegmentSet[i]);
            endS1 = getEndPoint(fuzzySegmentSet[i]);
            startS2 = getStartPoint(fuzzySegmentSet[i-1]);
            endS2 = getEndPoint(fuzzySegmentSet[i-1]);
            N1 = directionOfVector(N1,startS1,endS1);
            N2 = directionOfVector(N2,startS2,endS2);

            isIncrease_prev = isIncrease;
            isIncrease = isIncreasingOrder(N1,N2);
            i--;
        }
        if(isIncrease_prev == isIncrease) r = i;
        else r = i+1;

        if(verbose)
            cout<<"==> monoton seq is ("<<r<<" , "<<p<<")"<<endl;

        ///r = q;
        //find point with max (local) curvature
        int indexB = findElement(aContour,getStartPoint(fuzzySegmentSet[p]));
        int indexE = findElement(aContour,getEndPoint(fuzzySegmentSet[r]));
        int indexC = indexB;
        if(verbose)
            cout<<"curvature at "<<indexB<<" is "<<curvature[indexB]<<endl;

        for(int i = indexB+1; i<=indexE && i<curvature.size() ; i++)
        {
            if(fabs(curvature[i]) >= fabs(curvature[indexC]))
                indexC = i;
            if(verbose)
                cout<<"curvature at "<<i<<" is "<<curvature[i]<<endl;
        }
        if(verbose)
            cout<<" has max local curvature is "<<curvature[indexC]<<" at "<<aContour[indexC]<<endl;

        if(DP.empty() || (aContour[indexC] != DP.back()))
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
        for(vector<Point>::const_iterator it = commonPoint.begin(); it != commonPoint.end(); it++)
            idCommonPoint.push_back(findElement(aContour,*it));
        vector<int> idDominantPoint;
        for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
            idDominantPoint.push_back(findElement(aContour,*it));
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
        aBoard.saveSVG(filename);
    }
    return DP;
}
/* q = p-1 with ordering ANGLE of sequence of blurred segments */

/* q = p-1 without ordering ANGLE of sequence of blurred segments (only common zone) */
vector<Point> testDominantPointOnShapeV4(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, bool isClosed, const char* filename, bool verbose)
{
    int p=1, q=0, Eq=0, Bp=0, m=fuzzySegmentSet.size();
    vector<Point> pile;
    while (p<m && q<m && Bp>=0 && Eq>=0)
    {
        Bp = findStartElement(aContour,fuzzySegmentSet[p]);
        Eq = findEndElement(aContour,fuzzySegmentSet[q]);
        while(Eq>=Bp && p<m && q<m && Bp>=0 && Eq>=0 )
        {
            Bp = findStartElement(aContour,fuzzySegmentSet[p]);
            Eq = findEndElement(aContour,fuzzySegmentSet[q]);
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
    for(vector<Point>::const_iterator it = pile.begin(); it != pile.end(); it++)
    {
        q = (*it)[0];
        p = (*it)[1];
        //find mid point
        int indexB = findStartElement(aContour,fuzzySegmentSet[p]);
        int indexE = findEndElement(aContour,fuzzySegmentSet[q]);
        int indexC = (indexE+indexB)/2;
        if(verbose)
            cout<<" ==> center point is "<<aContour[indexC]<<endl;

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
        for(vector<Point>::const_iterator it = commonPoint.begin(); it != commonPoint.end(); it++)
            idCommonPoint.push_back(findElement(aContour,*it));
        vector<int> idDominantPoint;
        for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
            idDominantPoint.push_back(findElement(aContour,*it));
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
        aBoard.saveSVG(filename);
    }
    return DP;
}

vector<Point> testDominantPointOnShapeV4(const vector<AlphaThickSegmentComputer2D>& fuzzySegmentSet, const vector<Point>& aContour, bool isSymmetry, bool isClosed, const char* filename, bool verbose)
{
    int p=1, q=0, Eq=0, Bp=0, m=fuzzySegmentSet.size();
    vector<Point> pile;
    while (p<m && q<m && Bp>=0 && Eq>=0)
    {
        Bp = findStartElement(aContour,fuzzySegmentSet[p]);
        Eq = findEndElement(aContour,fuzzySegmentSet[q]);
        while(Eq>=Bp && p<m && q<m && Bp>=0 && Eq>=0 )
        {
            Bp = findStartElement(aContour,fuzzySegmentSet[p]);
            Eq = findEndElement(aContour,fuzzySegmentSet[q]);
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
    for(vector<Point>::const_iterator it = pile.begin(); it != pile.end(); it++)
    {
        q = (*it)[0];
        p = (*it)[1];
        //find point with max (local) curvature
        int indexB = findStartElement(aContour,fuzzySegmentSet[p]);
        int indexE = findEndElement(aContour,fuzzySegmentSet[q]);
        int indexC = indexB;
        int indexBi = findStartElement(aContour,fuzzySegmentSet[q]);
        int indexEi = findEndElement(aContour,fuzzySegmentSet[p]);
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
        for(vector<Point>::const_iterator it = commonPoint.begin(); it != commonPoint.end(); it++)
            idCommonPoint.push_back(findElement(aContour,*it));
        vector<int> idDominantPoint;
        for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
            idDominantPoint.push_back(findElement(aContour,*it));
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
        aBoard.saveSVG(filename);
    }
    return DP;
}
/* q = p-1 without ordering ANGLE of sequence of blurred segments (only common zone) */

/************************************/
/**** Dominant points detections ****/
/************************************/


/***********************************/
/**** Dominant points selection ****/
/***********************************/

/* Selection by curvature value */
vector<Point> testDominantPointSelectionV1(const vector<Point>& DP, const vector<int>& indexDP, const vector<double>& curvature, int nbDP, const vector<Point>& aContour, const char* filename, bool verbose)
{
    vector<double> curvatureDP;
    for(vector<int>::const_iterator it = indexDP.begin(); it != indexDP.end(); it++)
        curvatureDP.push_back(curvature.at(*it));

    vector<int> indexOrderedDP = absSortIndex(curvatureDP,false);
    if(verbose)
        for(vector<int>::const_iterator it = indexOrderedDP.begin(); it != indexOrderedDP.end(); it++)
            cout<<"curvature of DP : "<<curvatureDP.at(*it)<<endl;

    vector<Point> orderedSelectedDP;
    if(nbDP<DP.size())
    {
        vector<Point> selectedDP;
        //for(vector<int>::const_iterator it = indexOrderedDP.begin(); it != indexOrderedDP.end(); it++)
        for(int it=0; it<nbDP; it++)
            selectedDP.push_back(DP.at(indexOrderedDP.at(it)));
        vector<int> indexSelectedDP, indexOrderedSelectedDP;
        for(vector<Point>::const_iterator it = selectedDP.begin(); it != selectedDP.end(); it++)
            indexSelectedDP.push_back(findElement(aContour,*it));
        indexOrderedSelectedDP = sortIndex(indexSelectedDP);
        for(vector<int>::const_iterator it = indexOrderedSelectedDP.begin(); it != indexOrderedSelectedDP.end(); it++)
            orderedSelectedDP.push_back(selectedDP.at(*it));
    }
    else
        for(vector<Point>::const_iterator it = DP.begin(); it != DP.end(); it++)
            orderedSelectedDP.push_back(*it);

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
                   << CustomStyle( (*it).className(), new CustomColors( Color(255,0,0), Color(192,0,0)) )
                   <<(*it);
        /* Display old dominant point */
        /* Display new dominant point */
        for(vector<Point>::const_iterator it = orderedSelectedDP.begin(); it != orderedSelectedDP.end(); it++)
            aBoard << SetMode((*it).className(), "Both")
                   << CustomStyle( (*it).className(), new CustomColors( Color(0,255,0), Color(0,192,0)) )
                   <<(*it);
        /* Display new dominant point */

        /* Display the segments by old DP */
        aBoard.setPenColor(Color(192, 0, 0));
        aBoard.setLineWidth(10.0);
        aBoard.drawLine(aContour.front()[0],aContour.front()[1],DP.front()[0],DP.front()[1]);
        for(vector<Point>::const_iterator it = DP.begin(); it+1 != DP.end(); it++)
            aBoard.drawLine((*it)[0],(*it)[1],(*(it+1))[0],(*(it+1))[1]);
        aBoard.drawLine(aContour.back()[0],aContour.back()[1],DP.back()[0],DP.back()[1]);
        /* Display the segments by old DP */
        /* Display the segments by new DP */
        aBoard.setPenColor(Color(0, 192, 0));
        for(vector<Point>::const_iterator it = orderedSelectedDP.begin(); it+1 != orderedSelectedDP.end(); it++)
            aBoard.drawLine((*it)[0],(*it)[1],(*(it+1))[0],(*(it+1))[1]);
        /* Display the segments by new DP */
        aBoard.saveSVG(filename);
    }
    return orderedSelectedDP;
}
/* Selection by curvature value */

/* Selection by ise value */
vector<Point> testDominantPointSelectionV2(const vector<Point>& DP, const vector<int>& indexDP, int nbDP, const vector<Point>& aContour, bool isClosed, const char* filename, bool verbose)
{
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
                double ise = 0.0, d = 0.0;
                int indexStart = indexSelectedDP.at(it-1);
                int indexEnd = indexSelectedDP.at(it+1);
                if(indexStart < 0 || indexEnd < 0)
                    cout<<"Pb in error_ISE index of DP"<<endl;
                for(int index = indexStart+1; index < indexEnd; index++)
                {
                    d = distancePointSegment(aContour.at(index),aContour.at(indexStart), aContour.at(indexEnd));
                    ise += d*d;
                }
                if(verbose)
                    cout<<" DP "<<it<<" : "<<ise<<endl;
                //iseDP.push_back(ise/(indexEnd - indexStart - 1));
                iseDP.push_back(ise);
            }
            iseDP.insert(iseDP.begin(),aContour.size()*aContour.size());//Keep the fist point as DP
            if(!isClosed)//Keep the last point as DP
                iseDP.push_back(aContour.size()*aContour.size());
            else//isClosed => calcule ISE of the last point
            {
                double ise = 0.0, d = 0.0;
                int indexStart = indexSelectedDP.at(indexSelectedDP.size()-2);
                int indexEnd = aContour.size();
                //cout<<"it = "<<it<<" : indexStart = "<<indexStart<<", indexEnd = "<<indexEnd<<endl;
                for(int index = indexStart+1; index < indexEnd; index++)
                {
                    d = distancePointSegment(aContour.at(index),aContour.at(indexStart), aContour.front());
                    //cout<<"index "<<index<<" : distance "<<d<<endl;
                    ise += d*d;
                }
                if(verbose)
                    cout<<" DP "<<indexSelectedDP.at(indexSelectedDP.size()-1)<<" : "<<ise<<endl;
                //iseDP.push_back(ise/(indexEnd - indexStart - 1));
                iseDP.push_back(ise);
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
                   << CustomStyle( (*it).className(), new CustomColors( Color(255,0,0), Color(192,0,0)) )
                   <<(*it);
        /* Display old dominant point */
        /* Display new dominant point */
        for(vector<Point>::const_iterator it = selectedDP.begin(); it != selectedDP.end(); it++)
            aBoard << SetMode((*it).className(), "Both")
                   << CustomStyle( (*it).className(), new CustomColors( Color(0,255,0), Color(0,192,0)) )
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
        aBoard.saveSVG(filename);
    }
    return selectedDP;
}
/* Selection by ise value */

/* Selection by ise*curvature value */
vector<Point> testDominantPointSelectionV2(const vector<Point>& DP, const vector<int>& indexDP, const vector<double>& curvature, int nbDP, const vector<Point>& aContour, bool isClosed, const char* filename, bool verbose)
{
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
                double ise = 0.0, d = 0.0;
                int indexStart = indexSelectedDP.at(it-1);
                int indexEnd = indexSelectedDP.at(it+1);
                if(indexStart < 0 || indexEnd < 0)
                    cout<<"Pb in error_ISE index of DP"<<endl;
                for(int index = indexStart+1; index < indexEnd; index++)
                {
                    d = distancePointSegment(aContour.at(index),aContour.at(indexStart), aContour.at(indexEnd));
                    ise += d*d;
                }
                if(verbose)
                    cout<<" DP "<<it<<" : "<<ise*curvature.at(indexSelectedDP.at(it))<<endl;
                iseDP.push_back(ise*curvature.at(indexSelectedDP.at(it))*curvature.at(indexSelectedDP.at(it)));
            }
            iseDP.insert(iseDP.begin(),aContour.size()*aContour.size());//Keep the fist point as DP
            if(!isClosed)//Keep the last point as DP
                iseDP.push_back(aContour.size()*aContour.size());
            else//isClosed => calcule ISE of the last point
            {
                double ise = 0.0, d = 0.0;
                int indexStart = indexSelectedDP.at(indexSelectedDP.size()-2);
                int indexEnd = aContour.size();
                //cout<<"it = "<<it<<" : indexStart = "<<indexStart<<", indexEnd = "<<indexEnd<<endl;
                for(int index = indexStart+1; index < indexEnd; index++)
                {
                    d = distancePointSegment(aContour.at(index),aContour.at(indexStart), aContour.front());
                    //cout<<"index "<<index<<" : distance "<<d<<endl;
                    ise += d*d;
                }
                if(verbose)
                    cout<<" DP "<<indexSelectedDP.at(indexSelectedDP.size()-1)<<" : "<<ise*curvature.at(indexSelectedDP.at(indexSelectedDP.size()-1))<<endl;
                //iseDP.push_back(ise/(indexEnd - indexStart - 1));
                iseDP.push_back(ise*curvature.at(indexSelectedDP.at(indexSelectedDP.size()-1)));
            }

            vector<int> indexOrderedDP = absSortIndex(iseDP,true);
            if(verbose)
                for(vector<int>::const_iterator it = indexOrderedDP.begin(); it != indexOrderedDP.end(); it++)
                    cout<<"ise of DP "<<*it<< " : "<<iseDP.at(*it)<<endl;

            /* erease the first elt */
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
                   << CustomStyle( (*it).className(), new CustomColors( Color(255,0,0), Color(192,0,0)) )
                   <<(*it);
        /* Display old dominant point */
        /* Display new dominant point */
        for(vector<Point>::const_iterator it = selectedDP.begin(); it != selectedDP.end(); it++)
            aBoard << SetMode((*it).className(), "Both")
                   << CustomStyle( (*it).className(), new CustomColors( Color(0,255,0), Color(0,192,0)) )
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
        aBoard.saveSVG(filename);
    }
    return selectedDP;
}
/* Selection by ise*curvature value */

/* Selection by ise*angle */
vector<Point> testDominantPointSelectionV3(const vector<Point>& DP, const vector<int>& indexDP, int nbDP, const vector<Point>& aContour, bool isClosed, const char* filename, bool verbose)
{
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
        aBoard.saveSVG(filename);
    }
    return selectedDP;
}
/* Selection by ise*angle */

/* Selection by max error criterion */
vector<Point> testDominantPointSelectionV4(const vector<Point>& DP, const vector<int>& indexDP, const vector<Point>& aContour, bool isClosed, const char* filename, bool verbose)
{
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
        aBoard.saveSVG(filename);
    }
    return selectedDP;
}
/* Selection by max error criterion */

/************************************/
/**** Dominant points selection ****/
/************************************/


/**************************************/
/**** Tangent space transformation ****/
/**************************************/

vector<PointD> tangentspaceV1(const vector<Point>& DP, bool normalized, double* totalLength,const char* filename, bool verbose)
{
    vector<PointD> MidPoints;
    vector<double> alpha;
    vector<double> length;
    double l = 0.0, a = 0.0, totalAngle=0.0;
    *totalLength=0;

    //scan the Dominant Points
    int count = 0;
    for(vector<Point>::const_iterator it = DP.begin(); it+1 != DP.end(); it++)
    {
        //2 consecutive DP = 1 segment => store the length and angle between segments
        l = distancePoints(*it,*(it+1));
        length.push_back(l);
        *totalLength +=l;
        if(it != DP.begin())
        {
            Point v1 = Point((*it)[0] - (*(it-1))[0], (*it)[1] - (*(it-1))[1]);
            Point v2 = Point((*(it+1))[0] - (*it)[0],(*(it+1))[1]-(*it)[1]);
            a = signedAngle(v2,v1);
            alpha.push_back(a);
            totalAngle += a;
        }

        count++;
        if(verbose)
            cout<<count<<" : totalLength = "<<*totalLength<<" and totalAngle = "<<totalAngle<<" (angle = "<<a<<", length = "<<l<<")"<<endl;
    }

    if(length.empty())
        return MidPoints;

    //transform into the tangent space (x,y) = (length,alpha)
    vector<PointD> Ti1,Ti2,Ti;
    PointD p(0.0,0.0);
    if(normalized==true)
    {
        Ti.push_back(p);
        p[0] = length[0]/(*totalLength);
        p[1] = 0.0;
        Ti1.push_back(p);
        Ti.push_back(p);
        if(verbose)
            //cout<<"Ti (1): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<endl;
            cout<<"Tangent space point : "<<Ti.back()<<endl;
        p[0] = Ti1.back()[0] ;
        p[1] = Ti1.back()[1]+alpha[0];
        Ti2.push_back(p);
        Ti.push_back(p);
        if(verbose)
            //cout<<"Ti (2): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
            cout<<"Tangent space point : "<<Ti.back()<<endl;
        for(int it = 1; it < alpha.size(); it++)
        {
            p[0] = Ti2.back()[0] +length[it]/(*totalLength);
            p[1] = Ti2.back()[1];
            Ti1.push_back(p);
            Ti.push_back(p);
            if(verbose)
                //cout<<"Ti (1): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
                cout<<"Tangent space point : "<<Ti.back()<<endl;
            p[0] = Ti1.back()[0];
            p[1] = Ti1.back()[1]+alpha[it];
            Ti2.push_back(p);
            Ti.push_back(p);
            if(verbose)
                //cout<<"Ti (2): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
                cout<<"Tangent space point : "<<Ti.back()<<endl;
        }
        Ti.push_back(PointD(Ti2.back()[0] +length[alpha.size()]/(*totalLength),Ti2.back()[1]));
        if(verbose)
            //cout<<"Ti (3): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
            cout<<"Tangent space point : "<<Ti.back()<<endl;
    }
    else
    {
        Ti.push_back(p);
        p[0] = length[0];
        p[1] = 0.0;
        Ti1.push_back(p);
        Ti.push_back(p);
        if(verbose)
            //cout<<"Ti (1): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<endl;
            cout<<"Tangent space point : "<<Ti.back()<<endl;
        p[0] = Ti1.back()[0] ;
        p[1] = Ti1.back()[1]+alpha[0];
        Ti2.push_back(p);
        Ti.push_back(p);
        if(verbose)
            //cout<<"Ti (2): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
            cout<<"Tangent space point : "<<Ti.back()<<endl;
        for(int it = 1; it < alpha.size(); it++)
        {
            p[0] = Ti2.back()[0] +length[it];
            p[1] = Ti2.back()[1];
            Ti1.push_back(p);
            Ti.push_back(p);
            if(verbose)
                //cout<<"Ti (1): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
                cout<<"Tangent space point : "<<Ti.back()<<endl;
            p[0] = Ti1.back()[0];
            p[1] = Ti1.back()[1]+alpha[it];
            Ti2.push_back(p);
            Ti.push_back(p);
            if(verbose)
                //cout<<"Ti (2): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
                cout<<"Tangent space point : "<<Ti.back()<<endl;
        }
        Ti.push_back(PointD(Ti2.back()[0] +length[alpha.size()],Ti2.back()[1]));
        if(verbose)
            //cout<<"Ti (3): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
            cout<<"Tangent space point : "<<Ti.back()<<endl;
    }

    for(vector<PointD>::const_iterator it = Ti.begin(); it != Ti.end(); it+=2) //it++
    {
        MidPoints.push_back(PointD(((*it)[0]+(*(it+1))[0])/2.0,((*it)[1]+(*(it+1))[1])/2.0));
        if(verbose)
            cout<<"Middle Points : "<<MidPoints.back()<<endl;
    }

    if(filename != NULL)
    {
        Board2D aBoard;
        aBoard.setPenColor(Color( 255, 0, 0 ));
        if(!normalized)
        {
            aBoard.setLineWidth(1.0);
            for(int it = 0; it < Ti.size()-1; it++)
                aBoard.drawLine(Ti[it][0],Ti[it][1],Ti[it+1][0],Ti[it+1][1]);
            for(vector<PointD>::const_iterator it = MidPoints.begin(); it != MidPoints.end(); it++)
                aBoard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,192,0)) )
                       <<*it;
        }
        else
        {
            for(int it = 0; it < Ti.size()-1; it++)
                aBoard.drawLine(Ti[it][0],Ti[it][1],Ti[it+1][0],Ti[it+1][1]);
            aBoard.setPenColor(Color( 0, 255, 0));
            for(vector<PointD>::const_iterator it = MidPoints.begin(); it != MidPoints.end(); it++)
                aBoard.drawDot((*it)[0],(*it)[1],1);
        }
        aBoard.saveSVG(filename);
    }
    return MidPoints;
}

vector<PointD> tangentspaceV2(const vector<Point>& DP, bool normalized, vector<double>& alpha, vector<double>& length, double* totalLength, const char* filename, bool verbose)
{
    vector<PointD> MidPoints;
    double l = 0.0, a = 0.0, totalAngle=0.0;
    *totalLength=0;

    //scan the Dominant Points
    int count = 0;
    alpha.push_back(a);
    for(vector<Point>::const_iterator it = DP.begin(); it+1 != DP.end(); it++)
    {
        //2 consecutive DP = 1 segment => store the length and angle between segments
        l = distancePoints(*it,*(it+1));
        length.push_back(l);
        *totalLength +=l;
        //if(count==42)
        //    cout<<"STOP :"<<*(it-1)<<*(it)<<*(it+1)<<endl;
        if(it != DP.begin())
        {
            Point v1 = Point((*it)[0] - (*(it-1))[0], (*it)[1] - (*(it-1))[1]);
            Point v2 = Point((*(it+1))[0] - (*it)[0],(*(it+1))[1]-(*it)[1]);
            a = signedAngle(v1,v2);
            alpha.push_back(a);
            totalAngle += a;
        }
        count++;
        if(verbose)
            cout<<count<<" : totalLength = "<<*totalLength<<" and totalAngle = "<<totalAngle<<" (angle = "<<a<<", length = "<<l<<")"<<endl;
    }

    if(length.empty())
        return MidPoints;

    //transform into the tangent space (x,y) = (length,alpha)
    vector<PointD> Ti1,Ti2,Ti;
    PointD p(0.0,0.0);
    Ti.push_back(p);
    Ti2.push_back(p);
    if(normalized==true)
    {
        int it;
        for(it = 0; it < alpha.size()-1; it++)
        {
            p[0] = Ti2.back()[0] +length[it]/(*totalLength);
            p[1] = Ti2.back()[1];
            Ti1.push_back(p);
            Ti.push_back(p);
            if(verbose)
                //cout<<"Ti (1): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
                cout<<"Tangent space point : "<<Ti.back()<<endl;
            p[0] = Ti1.back()[0];
            p[1] = Ti1.back()[1]+alpha[it+1];
            Ti2.push_back(p);
            Ti.push_back(p);
            if(verbose)
                //cout<<"Ti (2): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
                cout<<"Tangent space point : "<<Ti.back()<<endl;
        }
        p[0] = Ti2.back()[0] +length[it]/(*totalLength);
        p[1] = Ti2.back()[1];
        Ti.push_back(p);
        if(verbose)
            //cout<<"Ti (1): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
            cout<<"Tangent space point : "<<Ti.back()<<endl;
    }
    else
    {
        int it;
        for(it = 0; it < alpha.size()-1; it++)
        {
            p[0] = Ti2.back()[0]+length[it];
            p[1] = Ti2.back()[1];
            Ti1.push_back(p);
            Ti.push_back(p);
            if(verbose)
                //cout<<"Ti (1): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
                cout<<"Tangent space point : "<<Ti.back()<<endl;
            p[0] = Ti1.back()[0];
            p[1] = Ti1.back()[1]+alpha[it+1];
            Ti2.push_back(p);
            Ti.push_back(p);
            if(verbose)
                //cout<<"Ti (2): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
                cout<<"Tangent space point : "<<Ti.back()<<endl;
        }
        p[0] = Ti2.back()[0]+length[it];
        p[1] = Ti2.back()[1];
        Ti.push_back(p);
        if(verbose)
            //cout<<"Ti (1): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
            cout<<"Tangent space point : "<<Ti.back()<<endl;
    }

    for(vector<PointD>::const_iterator it = Ti.begin(); it != Ti.end(); it+=2) //it++
    {
        MidPoints.push_back(PointD(((*it)[0]+(*(it+1))[0])/2.0,((*it)[1]+(*(it+1))[1])/2.0));
        if(verbose)
            cout<<"Middle Points : "<<MidPoints.back()<<endl;
    }

    if(filename != NULL)
    {
        Board2D aBoard;
        aBoard.setPenColor(Color( 255, 0, 0 ));
        if(!normalized)
        {
            aBoard.setLineWidth(1.0);
            for(int it = 0; it < Ti.size()-1; it++)
                aBoard.drawLine(Ti[it][0],Ti[it][1],Ti[it+1][0],Ti[it+1][1]);
            for(vector<PointD>::const_iterator it = MidPoints.begin(); it != MidPoints.end(); it++)
                aBoard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,192,0)) )
                       <<*it;
        }
        else
        {
            for(int it = 0; it < Ti.size()-1; it++)
                aBoard.drawLine(Ti[it][0],Ti[it][1],Ti[it+1][0],Ti[it+1][1]);
            aBoard.setPenColor(Color( 0, 255, 0));
            for(vector<PointD>::const_iterator it = MidPoints.begin(); it != MidPoints.end(); it++)
                aBoard.drawDot((*it)[0],(*it)[1],1);
        }
        aBoard.saveSVG(filename);
    }

    if(IS_PLOT)
    {
        writeFile(Ti,"Pts_TangentSpace.txt",false);
        writeFile(MidPoints,"MP_TangentSpace.txt",false);
    }
    return MidPoints;
}

vector<PointD> tangentspaceV3(const vector<Point>& DP, bool normalized, vector<PointD>& Ti, vector<double>& alpha, vector<double>& length, double* totalLength, const char* filename, bool verbose)
{
    vector<PointD> MidPoints;
    double l = 0.0, a = 0.0, totalAngle=0.0;
    *totalLength=0;

    //scan the Dominant Points
    int count = 0;
    alpha.push_back(a);
    for(vector<Point>::const_iterator it = DP.begin(); it+1 != DP.end(); it++)
    {
        //2 consecutive DP = 1 segment => store the length and angle between segments
        l = distancePoints(*it,*(it+1));
        length.push_back(l);
        *totalLength +=l;
        if(it != DP.begin())
        {
            Point v1 = Point((*it)[0] - (*(it-1))[0], (*it)[1] - (*(it-1))[1]);
            Point v2 = Point((*(it+1))[0] - (*it)[0],(*(it+1))[1]-(*it)[1]);
            a = signedAngle(v2,v1);
            alpha.push_back(a);
            totalAngle += a;
        }
        count++;
        if(verbose)
            cout<<count<<" : totalLength = "<<*totalLength<<" and totalAngle = "<<totalAngle<<" (angle = "<<a<<", length = "<<l<<")"<<endl;
    }

    if(length.empty())
        return MidPoints;

    //transform into the tangent space (x,y) = (length,alpha)
    vector<PointD> Ti1,Ti2;
    PointD p(0.0,0.0);
    Ti.push_back(p);
    Ti2.push_back(p);
    if(normalized==true)
    {
        int it;
        for(it = 0; it < alpha.size()-1; it++)
        {
            p[0] = Ti2.back()[0] +length[it]/(*totalLength);
            p[1] = Ti2.back()[1];
            Ti1.push_back(p);
            Ti.push_back(p);
            if(verbose)
                //cout<<"Ti (1): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
                cout<<"Tangent space point : "<<Ti.back()<<endl;
            p[0] = Ti1.back()[0];
            p[1] = Ti1.back()[1]+alpha[it+1];
            Ti2.push_back(p);
            Ti.push_back(p);
            if(verbose)
                //cout<<"Ti (2): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
                cout<<"Tangent space point : "<<Ti.back()<<endl;
        }
        p[0] = Ti2.back()[0] +length[it]/(*totalLength);
        p[1] = Ti2.back()[1];
        Ti.push_back(p);
        if(verbose)
            //cout<<"Ti (1): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
            cout<<"Tangent space point : "<<Ti.back()<<endl;
    }
    else
    {
        int it;
        for(it = 0; it < alpha.size()-1; it++)
        {
            p[0] = Ti2.back()[0] +length[it];
            p[1] = Ti2.back()[1];
            Ti1.push_back(p);
            Ti.push_back(p);
            if(verbose)
                //cout<<"Ti (1): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
                cout<<"Tangent space point : "<<Ti.back()<<endl;
            p[0] = Ti1.back()[0];
            p[1] = Ti1.back()[1]+alpha[it+1];
            Ti2.push_back(p);
            Ti.push_back(p);
            if(verbose)
                //cout<<"Ti (2): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
                cout<<"Tangent space point : "<<Ti.back()<<endl;
        }
        p[0] = Ti2.back()[0] +length[it];
        p[1] = Ti2.back()[1];
        Ti.push_back(p);
        if(verbose)
            //cout<<"Ti (1): "<<Ti.back()<<" and Ti1 : "<<Ti1.back()<<" and Ti2 "<<Ti2.back()<<endl;
            cout<<"Tangent space point : "<<Ti.back()<<endl;
    }

    for(vector<PointD>::const_iterator it = Ti.begin(); it != Ti.end(); it+=2) //it++
    {
        MidPoints.push_back(PointD(((*it)[0]+(*(it+1))[0])/2.0,((*it)[1]+(*(it+1))[1])/2.0));
        if(verbose)
            cout<<"Middle Points : "<<MidPoints.back()<<endl;
    }

    if(filename != NULL)
    {
        Board2D aBoard;
        aBoard.setPenColor(Color( 255, 0, 0 ));
        if(!normalized)
        {
            aBoard.setLineWidth(2.0);
            for(int it = 0; it < Ti.size()-1; it++)
                aBoard.drawLine(Ti[it][0],Ti[it][1],Ti[it+1][0],Ti[it+1][1]);
            for(vector<PointD>::const_iterator it = MidPoints.begin(); it != MidPoints.end(); it++)
                aBoard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,192,0)) )
                       <<*it;
        }
        else
        {
            for(int it = 0; it < Ti.size()-1; it++)
                aBoard.drawLine(Ti[it][0],Ti[it][1],Ti[it+1][0],Ti[it+1][1]);
            aBoard.setPenColor(Color( 0, 255, 0));
            for(vector<PointD>::const_iterator it = MidPoints.begin(); it != MidPoints.end(); it++)
                aBoard.drawDot((*it)[0],(*it)[1],1);
        }
        aBoard.saveSVG(filename);
    }
    return MidPoints;
}

/**************************************/
/**** Tangent space transformation ****/
/**************************************/




/****************************************/
/****Verification of islolated points ***/
/****************************************/

/* test isolated point w.r.t. alphaMax */
vector<PointD> testIsolatedPointsV1(const vector<PointD>& MP, double alphaMax, bool verbose)
{
    vector<PointD> isolatedPoints;
    bool res = false;
    vector<PointD>::const_iterator it;
    for (it = MP.begin(); it+1 != MP.end();it++)
    {
        res = isIsololatedPoint(*it, *(it+1),alphaMax);
        if(res==true)
            isolatedPoints.push_back(*it);

        if(verbose)
            cout<<"test isolated points of "<<*(it)<<" is "<<res<<endl;
    }
    return isolatedPoints;
}
/* test isolated point w.r.t. alphaMax */

/* test isolated point w.r.t. alphaMax + lengthMax */
set<int> testIsolatedPointsV2(const vector<double>& alpha, const vector<double>& length, double alphaMax, double maxLength, bool verbose)
{
    set<int> isolatedPoints;
    bool res = false;
    for (int it = 0; it < alpha.size(); it++)
    {
        res = false;
        if((length.at(it) >= maxLength) || (fabs(alpha.at(it)) >= alphaMax))
            res = true;

        if(res==true)
        {
            isolatedPoints.insert(it);
            isolatedPoints.insert(it+1);
        }

        if(verbose)
            cout<<"test isolated of point "<<it<<" is "<<res<<" ( angle = "<<alpha.at(it)<<" , length = "<<length.at(it)<<" )"<<endl;
    }
    if(res==false)
        isolatedPoints.insert(alpha.size());

    return isolatedPoints;
}
/* test isolated point w.r.t. alphaMax + lengthMax */

/****************************************/
/****Verification of islolated points ***/
/****************************************/




/*****************************************************/
/****Decomposition of Curve into Segments and Arcs ***/
/*****************************************************/
/* simple verif of mid points in tangent space => SEG:1,ARC:0,JOINCTION:-1 */
vector<int> testDecompositionSegmentCircleV1(const vector<int>& indexDP, const vector<PointD>& MP, double alphaMax, vector<Point>& segments, vector<Point>& arcs, bool verbose)
{
    /*********** Test of isolated points ***********/
    vector<int> isolatedVector;//SEG:1,ARC:0,JOINCTION:-1
    for (int it = 0; it < MP.size(); it++)
    {
        if(it==0)
        {
            if(fabs(MP.at(it)[1] - MP.at(it+1)[1]) >= alphaMax) //MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG : ";
                segments.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
                isolatedVector.push_back(1);//SEG
            }
            else//==> arc
            {
                if(verbose)
                    cout<<" ARC : ";
                arcs.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
                isolatedVector.push_back(0);//ARC
            }
            if(verbose)
                cout<<"different (V1) of alpha : "<<fabs(MP.at(it)[1] - MP.at(it+1)[1])<<endl;
        }
        else if(it==MP.size()-1)
        {
            if(fabs(MP.at(it)[1] - MP.at(it-1)[1]) >= alphaMax)//MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG : ";
                segments.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
                isolatedVector.push_back(1);//SEG
            }
            else
            {
                if(verbose)
                    cout<<" ARC : ";
                arcs.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
                isolatedVector.push_back(0);//ARC
            }
            if(verbose)
                cout<<"different (V1) of alpha : "<<fabs(MP.at(it)[1] - MP.at(it-1)[1])<<endl;
        }
        else
        {
            if( (fabs(MP.at(it)[1]-MP.at(it-1)[1])>alphaMax) && (fabs(MP.at(it)[1]-MP.at(it+1)[1])>alphaMax) )//MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG : ";
                segments.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
                isolatedVector.push_back(1);//SEG
                cout<<"ISOLATED POINT"<<endl;
            }
            else
            {
                if(verbose)
                    cout<<" ARC : ";
                arcs.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
                if((fabs(MP.at(it)[1]-MP.at(it-1)[1])>alphaMax))// || (fabs(MP.at(it)[1]-MP.at(it+1)[1])>alphaMax)) //MP is an isolated of jonction point
                {
                    isolatedVector.push_back(-1);//JONCTION
                    cout<<"JOINCTION POINT"<<endl;
                }
                else
                    isolatedVector.push_back(0);//ARC
            }
            if(verbose)
                cout<<"different (V1) of alpha : "<<fabs(MP.at(it)[1] - MP.at(it-1)[1])<<" and "<<fabs(MP.at(it)[1] - MP.at(it+1)[1])<<endl;
        }
        ///cout<<"size : "<<aSegment.size()<<" and getNumberSegmentPoints "<<aSegment.getNumberSegmentPoints()<<" it_MP : "<<*(it_MP)<<endl;
    }
    /*********** Test of isolated points ***********/
    return isolatedVector;
}
/* simple verif of mid points in tangent space => SEG:1,ARC:0,JOINCTION:-1 */

/* Decomposition into arc and seg using blurred segments in tangent space */
vector<int> testDecompositionSegmentCircleV1(const vector<int>& indexDP, const vector<PointD>& MP, double alphaMax, double thickness, vector<Point>& segments, vector<Point>& arcs, const char* filename, bool verbose)
{
    /*********** Test of isolated points ***********/
    vector<int> isolatedVector;//SEG:1,ARC:0,JOINCTION:-1
    for (int it = 0; it < MP.size(); it++)
    {
        if(it==0)
        {
            if(fabs(MP.at(it)[1] - MP.at(it+1)[1]) >= alphaMax) //MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG : ";
                isolatedVector.push_back(1);//SEG
            }
            else//==> arc
            {
                if(verbose)
                    cout<<" ARC : ";
                isolatedVector.push_back(0);//ARC
            }
            if(verbose)
                cout<<"different (V1) of alpha : "<<fabs(MP.at(it)[1] - MP.at(it+1)[1])<<endl;
        }
        else if(it==MP.size()-1)
        {
            if(fabs(MP.at(it)[1] - MP.at(it-1)[1]) >= alphaMax)//MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG : ";
                isolatedVector.push_back(1);//SEG
            }
            else
            {
                if(verbose)
                    cout<<" ARC : ";
                isolatedVector.push_back(0);//ARC
            }
            if(verbose)
                cout<<"different (V1) of alpha : "<<fabs(MP.at(it)[1] - MP.at(it-1)[1])<<endl;
        }
        else
        {
            if( (fabs(MP.at(it)[1]-MP.at(it-1)[1])>alphaMax) && (fabs(MP.at(it)[1]-MP.at(it+1)[1])>alphaMax) )//MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG (ISOLATED POINT) ";
                isolatedVector.push_back(1);//SEG
            }
            else
            {
                if(verbose)
                    cout<<" ARC ";
                if((fabs(MP.at(it)[1]-MP.at(it-1)[1])>alphaMax))// || (fabs(MP.at(it)[1]-MP.at(it+1)[1])>alphaMax)) //MP is an isolated of jonction point
                {
                    isolatedVector.push_back(-1);//JONCTION
                    if(verbose)
                        cout<<"(JOINCTION POINT)";
                }
                else
                    isolatedVector.push_back(0);//ARC
            }
            if(verbose)
                cout<<": different (V1) of alpha : "<<fabs(MP.at(it)[1] - MP.at(it-1)[1])<<" and "<<fabs(MP.at(it)[1] - MP.at(it+1)[1])<<endl;
        }
        ///cout<<"size : "<<aSegment.size()<<" and getNumberSegmentPoints "<<aSegment.getNumberSegmentPoints()<<" it_MP : "<<*(it_MP)<<endl;
    }
    /*
    //for(vector<int>::const_iterator it = isolatedVector.begin(); it != isolatedVector.end(); it++)
    for(int it = 0; it < isolatedVector.size(); it++)
        cout<<"Point "<<it<<" is "<<isolatedVector.at(it)<<endl;
    cout<<"-------------------"<<endl;
    */
    /*********** Test of isolated points ***********/

    /*********** Decomposition by blurred segment ****/
    vector<bool> isArc;
    for(int it=0; it<MP.size(); it++)
        isArc.push_back(false);
    vector<AlphaThickSegmentComputer2DD> blurredSegmentTS;
    vector<PointD>::const_iterator it_MP = MP.begin();
    for(int it = 0; it < isolatedVector.size(); it++)
    {
        //cout<<"it = "<<it<<endl;
        int it_start = it;
        if(isolatedVector.at(it)!=1)
        {
            AlphaThickSegmentComputer2DD aSegment;
            aSegment.init(it_MP,thickness);
            while(it<isolatedVector.size() &&  aSegment.end() != MP.end() && aSegment.extendFront()) //isolatedVector.at(it)==0 &&
            {
                //cout<<"it = "<<it<<endl;
                it++;
                it_MP++;
            }
            if(aSegment.getNumberSegmentPoints() > 3)
            {
                //cout<<"ARC : it = "<<it_start<<", isolatedVector = "<<isolatedVector.at(it_start)<<endl;
                blurredSegmentTS.push_back(aSegment);
                for(int i = it_start; i<it; i++)
                {
                    isArc[i]=true;
                    //arcs.push_back(Point(indexDP.at(i),indexDP.at(i+1)));
                }
                arcs.push_back(Point(indexDP.at(it_start),indexDP.at(it)));
                it++;
                it_MP++;
            }
            if(getEndPoint(aSegment) == MP.back())
                break;
            it = it-2;
            it_MP--;
        }
        else
            it_MP++;
    }
    for(int i=0; i<MP.size(); i++)
    {
        if(isArc.at(i)==false)
            segments.push_back(Point(indexDP.at(i),indexDP.at(i+1)));
    }

    if(filename != NULL)
    {
        Board2D TangentSpaceBoard;
        // Draw Middle points in Tangent Space
        for (vector<PointD>::const_iterator it = MP.begin(); it != MP.end();it++)
        {
            TangentSpaceBoard << SetMode((*it).className(), "Both");
            TangentSpaceBoard << *it;
        }

        // Draw Segments in Tangent Space
        int id = 1;
        for (vector<AlphaThickSegmentComputer2DD>::const_iterator it = blurredSegmentTS.begin();it != blurredSegmentTS.end();it++)
        {
            TangentSpaceBoard << SetMode((*it).className(), "BoundingBox");
            TangentSpaceBoard << *it;
            if(verbose)
            {
                double mu, nu;
                PointD N;
                (*it).computeParalellStripParams(mu,N,nu);
                cout<<"Segs "<<id<<" : ("<<N[0]<<" , "<<N[1]<<" , "<<mu<<" , "<<nu<<" ) ==> ";
                cout<<"R = "<<fabs(N[1]/N[0])<<endl;//slope = -N0/N1
                //cout<<getStartPoint(*it)<<" and "<<getEndPoint(*it)<<endl;
                id++;
            }
        }
        TangentSpaceBoard.saveSVG(filename);
    }
    return isolatedVector;
}

vector<int> testDecompositionSegmentCircleV2(const vector<Point>& aContour, const vector<int>& indexDP, const vector<PointD>& MP, double alphaMax, double thickness, vector<Point>& segments, vector<Point>& arcs, const char* filename, bool verbose)
{
    /*********** Test of isolated points ***********/
    vector<int> isolatedVector;//SEG:1,ARC:0,JOINCTION:-1
    for (int it = 0; it < MP.size(); it++)
    {
        if(it==0)
        {
            if(fabs(MP.at(it)[1] - MP.at(it+1)[1]) >= alphaMax) //MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG : ";
                isolatedVector.push_back(1);//SEG
            }
            else//==> arc
            {
                if(verbose)
                    cout<<" ARC : ";
                isolatedVector.push_back(0);//ARC
            }
            if(verbose)
                cout<<"different (V1) of alpha : "<<fabs(MP.at(it)[1] - MP.at(it+1)[1])<<endl;
        }
        else if(it==MP.size()-1)
        {
            if(fabs(MP.at(it)[1] - MP.at(it-1)[1]) >= alphaMax)//MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG : ";
                isolatedVector.push_back(1);//SEG
            }
            else
            {
                if(verbose)
                    cout<<" ARC : ";
                isolatedVector.push_back(0);//ARC
            }
            if(verbose)
                cout<<"different (V1) of alpha : "<<fabs(MP.at(it)[1] - MP.at(it-1)[1])<<endl;
        }
        else
        {
            if( (fabs(MP.at(it)[1]-MP.at(it-1)[1])>alphaMax) && (fabs(MP.at(it)[1]-MP.at(it+1)[1])>alphaMax) )//MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG (ISOLATED POINT) ";
                isolatedVector.push_back(1);//SEG
            }
            else
            {
                if(verbose)
                    cout<<" ARC ";
                if((fabs(MP.at(it)[1]-MP.at(it-1)[1])>alphaMax))// || (fabs(MP.at(it)[1]-MP.at(it+1)[1])>alphaMax)) //MP is an isolated of jonction point
                {
                    isolatedVector.push_back(-1);//JONCTION
                    if(verbose)
                        cout<<"(JOINCTION POINT)";
                }
                else
                    isolatedVector.push_back(0);//ARC
            }
            if(verbose)
                cout<<": different (V1) of alpha : "<<fabs(MP.at(it)[1] - MP.at(it-1)[1])<<" and "<<fabs(MP.at(it)[1] - MP.at(it+1)[1])<<endl;
        }
        ///cout<<"size : "<<aSegment.size()<<" and getNumberSegmentPoints "<<aSegment.getNumberSegmentPoints()<<" it_MP : "<<*(it_MP)<<endl;
    }
    /*********** Test of isolated points ***********/

    /*********** Decomposition by blurred segment ****/
    vector<bool> isArc;
    for(int it=0; it<MP.size(); it++)
        isArc.push_back(false);
    vector<AlphaThickSegmentComputer2DD> blurredSegmentTS;
    vector<PointD>::const_iterator it_MP = MP.begin();
    for(int it = 0; it < isolatedVector.size(); it++)
    {
        //cout<<"it = "<<it<<endl;
        int it_start = it;
        if(isolatedVector.at(it)!=1)
        {
            AlphaThickSegmentComputer2DD aSegment;
            aSegment.init(it_MP,thickness);
            while(it<isolatedVector.size() &&  aSegment.end() != MP.end() && aSegment.extendFront()) //isolatedVector.at(it)==0 &&
            {
                it++;
                it_MP++;
            }
            //if(aSegment.getNumberSegmentPoints() > 2)//at least 4 points on the circle
            if(aSegment.getNumberSegmentPoints() >= 2)//3 points forms a circle
            {
                blurredSegmentTS.push_back(aSegment);
                //verify of ise between arc/circle vs. seg
                int idBegin = indexDP.at(it_start);
                int idEnd = indexDP.at(it);
                int idMid = findBestFittingCircle(aContour,idBegin,idEnd);
                Point center = determineCenter(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd));
                double radius = (determineRadius(center,aContour.at(idBegin)) + determineRadius(center,aContour.at(idMid)) + determineRadius(center,aContour.at(idEnd)))/3.0;
                double ise_Arc = iseContourCircle(aContour,idBegin,idEnd,center,radius);
                double ise_Seg = 0;
                for(int i = it_start; i<it; i++)
                    ise_Seg += iseContourSegment(aContour,indexDP.at(i),indexDP.at(i+1));
                if(verbose)
                    cout<<"ise_Arc = "<<ise_Arc<<" and ise_Seg = "<<ise_Seg<<endl;
                if(ise_Arc>2*ise_Seg)
                {
                    for(int i = it_start; i<it; i++)
                        isArc[i]=false;
                }
                else
                {
                    for(int i = it_start; i<it; i++)
                        isArc[i]=true;
                    arcs.push_back(Point(indexDP.at(it_start),indexDP.at(it)));
                }
                it++;
                it_MP++;
            }
            if(getEndPoint(aSegment) == MP.back())
                break;
            it = it-2;
            it_MP--;
        }
        else
            it_MP++;
    }
    for(int i=0; i<MP.size(); i++)
    {
        if(isArc.at(i)==false)
            segments.push_back(Point(indexDP.at(i),indexDP.at(i+1)));
    }

    if(filename != NULL)
    {
        Board2D TangentSpaceBoard;
        // Draw Middle points in Tangent Space
        for (vector<PointD>::const_iterator it = MP.begin(); it != MP.end();it++)
        {
            TangentSpaceBoard << SetMode((*it).className(), "Both");
            TangentSpaceBoard << *it;
        }

        // Draw Segments in Tangent Space
        int id = 1;
        for (vector<AlphaThickSegmentComputer2DD>::const_iterator it = blurredSegmentTS.begin();it != blurredSegmentTS.end();it++)
        {
            TangentSpaceBoard << SetMode((*it).className(), "BoundingBox");
            TangentSpaceBoard << *it;
            if(verbose)
            {
                double mu, nu;
                PointD N;
                (*it).computeParalellStripParams(mu,N,nu);
                cout<<"Segs "<<id<<" : ("<<N[0]<<" , "<<N[1]<<" , "<<mu<<" , "<<nu<<" ) ==> ";
                cout<<"R = "<<fabs(N[1]/N[0])<<endl;//slope = -N0/N1
                //cout<<getStartPoint(*it)<<" and "<<getEndPoint(*it)<<endl;
                id++;
            }
        }
        TangentSpaceBoard.saveSVG(filename);
    }
    return isolatedVector;
}

vector<int> testDecompositionSegmentCircleV3(const vector<Point>& aContour, const vector<int>& indexDP, const vector<PointD>& MP, double alphaMax, double thickness, vector<Point>& segments, vector<Point>& arcs, int nbPointCir, double iseTol, double angleTol, const char* filename, bool verbose)
{
    /*********** Test of isolated points ***********/
    vector<int> isolatedVector;//SEG:1,ARC:0,JOINCTION:-1
    for (int it = 0; it < MP.size(); it++)
    {
        if(it==0)
        {
            if(fabs(MP.at(it)[1] - MP.at(it+1)[1]) >= alphaMax) //MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG : ";
                isolatedVector.push_back(1);//SEG
            }
            else//==> arc
            {
                if(verbose)
                    cout<<" ARC : ";
                isolatedVector.push_back(0);//ARC
            }
            if(verbose)
                cout<<"different (V1) of alpha : "<<fabs(MP.at(it)[1] - MP.at(it+1)[1])<<endl;
        }
        else if(it==MP.size()-1)
        {
            if(fabs(MP.at(it)[1] - MP.at(it-1)[1]) >= alphaMax)//MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG : ";
                isolatedVector.push_back(1);//SEG
            }
            else
            {
                if(verbose)
                    cout<<" ARC : ";
                isolatedVector.push_back(0);//ARC
            }
            if(verbose)
                cout<<"different (V1) of alpha : "<<fabs(MP.at(it)[1] - MP.at(it-1)[1])<<endl;
        }
        else
        {
            if( (fabs(MP.at(it)[1]-MP.at(it-1)[1])>alphaMax) && (fabs(MP.at(it)[1]-MP.at(it+1)[1])>alphaMax) )//MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG (ISOLATED POINT) ";
                isolatedVector.push_back(1);//SEG
            }
            else
            {
                if(verbose)
                    cout<<" ARC ";
                if((fabs(MP.at(it)[1]-MP.at(it-1)[1])>alphaMax))// || (fabs(MP.at(it)[1]-MP.at(it+1)[1])>alphaMax)) //MP is an isolated of jonction point
                {
                    isolatedVector.push_back(-1);//JONCTION
                    if(verbose)
                        cout<<"(JOINCTION POINT)";
                }
                else
                    isolatedVector.push_back(0);//ARC
            }
            if(verbose)
                cout<<": different (V1) of alpha : "<<fabs(MP.at(it)[1] - MP.at(it-1)[1])<<" and "<<fabs(MP.at(it)[1] - MP.at(it+1)[1])<<endl;
        }
        ///cout<<"size : "<<aSegment.size()<<" and getNumberSegmentPoints "<<aSegment.getNumberSegmentPoints()<<" it_MP : "<<*(it_MP)<<endl;
    }
    /*********** Test of isolated points ***********/

    /*********** Decomposition by blurred segment ****/
    vector<bool> isArc;
    for(int it=0; it<MP.size(); it++)
        isArc.push_back(false);
    vector<AlphaThickSegmentComputer2DD> blurredSegmentTS;
    vector<PointD>::const_iterator it_MP = MP.begin();
    for(int it = 0; it < isolatedVector.size(); it++)
    {
        //cout<<"it = "<<it<<endl;
        int it_start = it;
        if(isolatedVector.at(it)!=1)
        {
            AlphaThickSegmentComputer2DD aSegment;
            aSegment.init(it_MP,thickness);
            while(it<isolatedVector.size() &&  aSegment.end() != MP.end() && aSegment.extendFront()) //isolatedVector.at(it)==0 &&
            {
                it++;
                it_MP++;
            }
            if(aSegment.getNumberSegmentPoints() >= nbPointCir)//at least (nbPointCir+1) points on the circle
            {
                blurredSegmentTS.push_back(aSegment);
                //verify of ise between arc/circle vs. seg
                int idBegin = indexDP.at(it_start);
                int idEnd = indexDP.at(it);
                int idMid = findBestFittingCircle(aContour,idBegin,idEnd);
                double linAngle = relativeAngle(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd))*180/M_PI;
                Point center = determineCenter(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd));
                double radius = (determineRadius(center,aContour.at(idBegin)) + determineRadius(center,aContour.at(idMid)) + determineRadius(center,aContour.at(idEnd)))/3.0;
                double ise_Arc = iseContourCircle(aContour,idBegin,idEnd,center,radius);
                double ise_Seg = 0;
                for(int i = it_start; i<it; i++)
                    ise_Seg += iseContourSegment(aContour,indexDP.at(i),indexDP.at(i+1));
                if(verbose)
                    cout<<"ise_Arc = "<<ise_Arc<<" and ise_Seg = "<<ise_Seg<<" and angle = "<<linAngle<<endl;

                if(ise_Arc>iseTol*ise_Seg || (180-linAngle)<angleTol)
                {
                    //cout<<"==> Ommited, ise_Arc = "<<ise_Arc<<" and ise_Seg = "<<ise_Seg<<" and angle = "<<linAngle<<endl;
                    for(int i = it_start; i<it; i++)
                        isArc[i]=false;
                }
                else
                {
                    for(int i = it_start; i<it; i++)
                        isArc[i]=true;
                    arcs.push_back(Point(indexDP.at(it_start),indexDP.at(it)));
                }
                it++;
                it_MP++;
            }
            if(getEndPoint(aSegment) == MP.back())
                break;
            it = it-2;
            it_MP--;
        }
        else
            it_MP++;
    }
    for(int i=0; i<MP.size(); i++)
    {
        if(isArc.at(i)==false)
            segments.push_back(Point(indexDP.at(i),indexDP.at(i+1)));
    }

    if(filename != NULL)
    {
        Board2D TangentSpaceBoard;
        // Draw Middle points in Tangent Space
        for (vector<PointD>::const_iterator it = MP.begin(); it != MP.end();it++)
        {
            TangentSpaceBoard << SetMode((*it).className(), "Both");
            TangentSpaceBoard << *it;
        }

        // Draw Segments in Tangent Space
        int id = 1;
        for (vector<AlphaThickSegmentComputer2DD>::const_iterator it = blurredSegmentTS.begin();it != blurredSegmentTS.end();it++)
        {
            TangentSpaceBoard << SetMode((*it).className(), "BoundingBox");
            TangentSpaceBoard << *it;
            if(verbose)
            {
                double mu, nu;
                PointD N;
                (*it).computeParalellStripParams(mu,N,nu);
                cout<<"Segs "<<id<<" : ("<<N[0]<<" , "<<N[1]<<" , "<<mu<<" , "<<nu<<" ) ==> ";
                cout<<"R = "<<fabs(N[1]/N[0])<<endl;//slope = -N0/N1
                //cout<<getStartPoint(*it)<<" and "<<getEndPoint(*it)<<endl;
                id++;
            }
        }
        TangentSpaceBoard.saveSVG(filename);
    }
    return isolatedVector;
}

vector<int> testDecompositionSegmentCircleV3(const vector<Point>& aContour, const vector<int>& indexDP, const vector<PointD>& MP, double alphaMax, double thickness, vector<Point>& segments, vector<Point>& arcs, vector<double>& slope, int nbPointCir, double iseTol, double angleTol, const char* filename, bool verbose)
{
    /*********** Test of isolated points ***********/
    vector<int> isolatedVector;//SEG:1,ARC:0,JOINCTION:-1
    for (int it = 0; it < MP.size(); it++)
    {
        if(it==0)
        {
            //if(fabs(MP.at(it)[1] - MP.at(it+1)[1]) >= alphaMax) //MP[0] is an isolated point => segment !
            if(fabs(MP.at(it)[1] - MP.at(it+1)[1]) > alphaMax) //MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG : ";
                isolatedVector.push_back(1);//SEG
            }
            else//==> arc
            {
                if(verbose)
                    cout<<" ARC : ";
                isolatedVector.push_back(0);//ARC
            }
            if(verbose)
                cout<<"different (V1) of alpha : "<<fabs(MP.at(it)[1] - MP.at(it+1)[1])<<endl;
        }
        else if(it==MP.size()-1)
        {
            //if(fabs(MP.at(it)[1] - MP.at(it-1)[1]) >= alphaMax)//MP[0] is an isolated point => segment !
            if(fabs(MP.at(it)[1] - MP.at(it-1)[1]) > alphaMax)//MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG : ";
                isolatedVector.push_back(1);//SEG
            }
            else // ARC or JUNCTION POINT
            {
                if(verbose)
                    cout<<" ARC : ";
                isolatedVector.push_back(0);//ARC
            }
            if(verbose)
                cout<<"different (V1) of alpha : "<<fabs(MP.at(it)[1] - MP.at(it-1)[1])<<endl;
        }
        else
        {
            if( (fabs(MP.at(it)[1]-MP.at(it-1)[1])>alphaMax) && (fabs(MP.at(it)[1]-MP.at(it+1)[1])>alphaMax) )//MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG (ISOLATED POINT) ";
                isolatedVector.push_back(1);//SEG
            }
            else
            {
                if(verbose)
                    cout<<" ARC ";
                if((fabs(MP.at(it)[1]-MP.at(it-1)[1])>alphaMax))// || (fabs(MP.at(it)[1]-MP.at(it+1)[1])>alphaMax)) //MP is an isolated of jonction point
                {
                    isolatedVector.push_back(-1);//JONCTION
                    if(verbose)
                        cout<<"(JOINCTION POINT)";
                }
                else
                    isolatedVector.push_back(0);//ARC
            }
            if(verbose)
                cout<<": different (V1) of alpha : "<<fabs(MP.at(it)[1] - MP.at(it-1)[1])<<" and "<<fabs(MP.at(it)[1] - MP.at(it+1)[1])<<endl;
        }
        ///cout<<"size : "<<aSegment.size()<<" and getNumberSegmentPoints "<<aSegment.getNumberSegmentPoints()<<" it_MP : "<<*(it_MP)<<endl;
    }
    /*********** Test of isolated points ***********/

    /*********** Decomposition by blurred segment ****/
    vector<bool> isArc;
    for(int it=0; it<MP.size(); it++)
        isArc.push_back(false);
    vector<AlphaThickSegmentComputer2DD> blurredSegmentTS;
    vector<PointD>::const_iterator it_MP = MP.begin();
    for(int it = 0; it < isolatedVector.size(); it++)
    {
        //cout<<"it = "<<it<<endl;
        int it_start = it;
        if(isolatedVector.at(it)!=1)
        {
            AlphaThickSegmentComputer2DD aSegment;
            aSegment.init(it_MP,thickness);
            while(it<isolatedVector.size() &&  aSegment.end() != MP.end() && aSegment.extendFront()) //isolatedVector.at(it)==0 &&
            {
                it++;
                it_MP++;
            }
            if(aSegment.getNumberSegmentPoints() >= nbPointCir)//at least (nbPointCir+1) points on the circle
            {
                blurredSegmentTS.push_back(aSegment);
                //verify of ise between arc/circle vs. seg
                int idBegin = indexDP.at(it_start);
                int idEnd = indexDP.at(it);
                int idMid = findBestFittingCircle(aContour,idBegin,idEnd);
                double linAngle = relativeAngle(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd))*180/M_PI;
                Point center = determineCenter(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd));
                double radius = (determineRadius(center,aContour.at(idBegin)) + determineRadius(center,aContour.at(idMid)) + determineRadius(center,aContour.at(idEnd)))/3.0;
                double ise_Arc = iseContourCircle(aContour,idBegin,idEnd,center,radius);
                double ise_Seg = 0;
                for(int i = it_start; i<it; i++)
                    ise_Seg += iseContourSegment(aContour,indexDP.at(i),indexDP.at(i+1));
                if(verbose)
                    cout<<"ise_Arc = "<<ise_Arc<<" and ise_Seg = "<<ise_Seg<<" and angle = "<<linAngle<<endl;

                if(ise_Arc>iseTol*ise_Seg || (180-linAngle)<angleTol)
                {
                    //cout<<"==> Ommited, ise_Arc = "<<ise_Arc<<" and ise_Seg = "<<ise_Seg<<" and angle = "<<linAngle<<endl;
                    for(int i = it_start; i<it; i++)
                        isArc[i]=false;
                }
                else
                {
                    for(int i = it_start; i<it; i++)
                        isArc[i]=true;
                    arcs.push_back(Point(indexDP.at(it_start),indexDP.at(it)));
                    /*
                    double mu, nu;
                    PointD N;
                    aSegment.computeParalellStripParams(mu,N,nu);
                    cout<<"Segs : ("<<N[0]<<" , "<<N[1]<<" , "<<mu<<" , "<<nu<<" ) ==> ";
                    cout<<"R = "<<fabs(N[1]/N[0])<<endl;//slope = -N0/N1
                    slope.push_back(getSlope(N[0],N[1]));
                    */
                    PointD s = getStartPoint(aSegment);
                    PointD e = getEndPoint(aSegment);
                    slope.push_back(fabs((s[0]-e[0])/(s[1]-e[1])));
                }
                it++;
                it_MP++;
            }
            if(getEndPoint(aSegment) == MP.back())
                break;
            it = it-2;
            it_MP--;
        }
        else
            it_MP++;
    }
    for(int i=0; i<MP.size(); i++)
    {
        if(isArc.at(i)==false)
            segments.push_back(Point(indexDP.at(i),indexDP.at(i+1)));
    }

    if(filename != NULL)
    {
        Board2D TangentSpaceBoard;
        // Draw Middle points in Tangent Space
        for (vector<PointD>::const_iterator it = MP.begin(); it != MP.end();it++)
        {
            TangentSpaceBoard << SetMode((*it).className(), "Both");
            TangentSpaceBoard << *it;
        }

        // Draw Segments in Tangent Space
        int id = 1;
        for (vector<AlphaThickSegmentComputer2DD>::const_iterator it = blurredSegmentTS.begin();it != blurredSegmentTS.end();it++)
        {
            TangentSpaceBoard << SetMode((*it).className(), "BoundingBox");
            TangentSpaceBoard << *it;
            if(verbose)
            {
                double mu, nu;
                PointD N;
                (*it).computeParalellStripParams(mu,N,nu);
                cout<<"Segs "<<id<<" : ("<<N[0]<<" , "<<N[1]<<" , "<<mu<<" , "<<nu<<" ) ==> ";
                cout<<"R = "<<fabs(N[1]/N[0])<<endl;//slope = -N0/N1
                //cout<<getStartPoint(*it)<<" and "<<getEndPoint(*it)<<endl;
                id++;
            }
        }
        TangentSpaceBoard.saveSVG(filename);
    }
    return isolatedVector;
}

void testDecompositionSegmentCircleV4(const vector<Point>& aContour, const vector<int>& indexDP, const vector<PointD>& MP, double thickness, vector<Point>& segments, vector<Point>& arcs, int nbPointCir, double iseTol, const char* filename, bool verbose)
{
    /*********** Decomposition by blurred segment ****/
    vector<bool> isArc;
    for(int it=0; it<MP.size(); it++)
        isArc.push_back(false);
    vector<AlphaThickSegmentComputer2DD> blurredSegmentTS;
    vector<PointD>::const_iterator it_MP = MP.begin();
    for(int it = 0; it < MP.size(); it++)
    {
        int it_start = it;
        ///if(isolatedVector.at(it)!=1)
        ///{
            AlphaThickSegmentComputer2DD aSegment;
            aSegment.init(it_MP,thickness);
            while(it<MP.size() &&  aSegment.end() != MP.end() && aSegment.extendFront()) //isolatedVector.at(it)==0 &&
            {
                it++;
                it_MP++;
            }
            if(aSegment.getNumberSegmentPoints() >= nbPointCir)//at least (nbPointCir+1) points on the circle
            {
                blurredSegmentTS.push_back(aSegment);
                //verify of ise between arc/circle vs. seg
                int idBegin = indexDP.at(it_start);
                int idEnd = indexDP.at(it);
                int idMid = findBestFittingCircle(aContour,idBegin,idEnd);
                Point center = determineCenter(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd));
                double radius = (determineRadius(center,aContour.at(idBegin)) + determineRadius(center,aContour.at(idMid)) + determineRadius(center,aContour.at(idEnd)))/3.0;
                double ise_Arc = iseContourCircle(aContour,idBegin,idEnd,center,radius);
                double ise_Seg = 0;
                for(int i = it_start; i<it; i++)
                    ise_Seg += iseContourSegment(aContour,indexDP.at(i),indexDP.at(i+1));
                if(verbose)
                    cout<<"ise_Arc = "<<ise_Arc<<" and ise_Seg = "<<ise_Seg<<endl;
                if(ise_Arc>iseTol*ise_Seg)
                {
                    for(int i = it_start; i<it; i++)
                        isArc[i]=false;
                }
                else
                {
                    for(int i = it_start; i<it; i++)
                        isArc[i]=true;
                    arcs.push_back(Point(indexDP.at(it_start),indexDP.at(it)));
                }
                it++;
                it_MP++;
            }
            if(getEndPoint(aSegment) == MP.back())
                break;
            it = it-2;
            it_MP--;
        ///}
        ///else
        ///    it_MP++;
    }
    for(int i=0; i<MP.size(); i++)
    {
        if(isArc.at(i)==false)
            segments.push_back(Point(indexDP.at(i),indexDP.at(i+1)));
    }

    if(filename != NULL)
    {
        Board2D TangentSpaceBoard;
        // Draw Middle points in Tangent Space
        for (vector<PointD>::const_iterator it = MP.begin(); it != MP.end();it++)
        {
            TangentSpaceBoard << SetMode((*it).className(), "Both");
            TangentSpaceBoard << *it;
        }

        // Draw Segments in Tangent Space
        int id = 1;
        for (vector<AlphaThickSegmentComputer2DD>::const_iterator it = blurredSegmentTS.begin();it != blurredSegmentTS.end();it++)
        {
            TangentSpaceBoard << SetMode((*it).className(), "BoundingBox");
            TangentSpaceBoard << *it;
            if(verbose)
            {
                double mu, nu;
                PointD N;
                (*it).computeParalellStripParams(mu,N,nu);
                cout<<"Segs "<<id<<" : ("<<N[0]<<" , "<<N[1]<<" , "<<mu<<" , "<<nu<<" ) ==> ";
                cout<<"R = "<<fabs(N[1]/N[0])<<endl;//slope = -N0/N1
                //cout<<getStartPoint(*it)<<" and "<<getEndPoint(*it)<<endl;
                id++;
            }
        }
        TangentSpaceBoard.saveSVG(filename);
    }
}
/* Decomposition into arc and seg using blurred segments in tangent space */

void testDecompositionSegmentCircleV2(const vector<int>& indexDP, const vector<double>& alpha, const vector<double>& length, double alphaMax, double maxLength, vector<Point>& segments, vector<Point>& arcs, bool verbose)
{
    int it = 0;
    if((length.at(it) > maxLength) || (fabs(alpha.at(it+1)) >= alphaMax))//MP[0] is an isolated point => segment !
    {
        if(verbose)
            cout<<" SEG : ";
        segments.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
    }
    else//==> arc
    {
        if(verbose)
            cout<<" ARC : ";
        arcs.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
    }
    if(verbose)
        cout<<"different (V2) of angle 1 = "<<fabs(alpha.at(it))<<" and angle 2 = "<< fabs(alpha.at(it+1)) <<" , length = "<<length.at(it)<<endl;

    for (it = 1; it < alpha.size()-1; it++)
    {
        if((length.at(it) >= maxLength) || ((fabs(alpha.at(it)) >= alphaMax) && (fabs(alpha.at(it+1)) >= alphaMax)))//MP[0] is an isolated point => segment !
        {
            if(verbose)
                cout<<" SEG : ";
            segments.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
        }
        else//==> arc
        {
            if(verbose)
                cout<<" ARC : ";
            arcs.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
        }
        if(verbose)
            cout<<"different (V2) of angle 1 = "<<fabs(alpha.at(it))<<" and angle 2 = "<< fabs(alpha.at(it+1)) <<" , length = "<<length.at(it)<<endl;
    }
    if((length.at(it) >= maxLength) || (fabs(alpha.at(it)) >= alphaMax))
    {
        if(verbose)
            cout<<" SEG : ";
        segments.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
    }
    else
    {
        if(verbose)
            cout<<" ARC : ";
        arcs.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
    }
    if(verbose)
        cout<<"different (V2) of angle = "<<fabs(alpha.at(it))<<" , length = "<<length.at(it)<<endl;
}

void testDecompositionSegmentCircleV2(const vector<int>& indexDP, const vector<PointD>& MP, const vector<double>& alpha, const vector<double>& length, double alphaMax, double maxLength, double thickness, vector<Point>& segments, vector<Point>& arcs, const char* filename, bool verbose)
{
    /*********** Test of isolated points ***********/
    vector<bool> isolatedVector;
    for (int it = 0; it < alpha.size(); it++)
    {
        if(it==0)
        {
            if((length.at(it) >= maxLength) || (fabs(alpha.at(it+1)) >= alphaMax))//MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG : ";
                segments.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
                isolatedVector.push_back(true);
            }
            else//==> arc
            {
                if(verbose)
                    cout<<" ARC : ";
                arcs.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
                isolatedVector.push_back(false);
            }
            if(verbose)
                cout<<"different (V2) of angle 1 = "<<fabs(alpha.at(it))<<" and angle 2 = "<< fabs(alpha.at(it+1)) <<" , length = "<<length.at(it)<<endl;
        }
        else if(it==alpha.size()-1)
        {
            if((length.at(it) >= maxLength) || (fabs(alpha.at(it)) >= alphaMax))
            {
                if(verbose)
                    cout<<" SEG : ";
                segments.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
                isolatedVector.push_back(true);
            }
            else
            {
                if(verbose)
                    cout<<" ARC : ";
                arcs.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
                isolatedVector.push_back(false);
            }
            if(verbose)
                cout<<"different (V2) of angle = "<<fabs(alpha.at(it))<<" , length = "<<length.at(it)<<endl;
        }
        else
        {
            if((length.at(it) >= maxLength) || ((fabs(alpha.at(it)) >= alphaMax) && (fabs(alpha.at(it+1)) >= alphaMax)))//MP[0] is an isolated point => segment !
            {
                if(verbose)
                    cout<<" SEG : ";
                segments.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
                isolatedVector.push_back(true);
            }
            else//==> arc
            {
                if(verbose)
                    cout<<" ARC : ";
                arcs.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
                isolatedVector.push_back(false);
            }
            if(verbose)
                cout<<"different (V2) of angle 1 = "<<fabs(alpha.at(it))<<" and angle 2 = "<< fabs(alpha.at(it+1)) <<" , length = "<<length.at(it)<<endl;
        }
    }

    /*
    for(vector<bool>::const_iterator it = isolatedVector.begin(); it!= isolatedVector.end(); it++)
    {
        if(it)
    }
    */
    /*********** Test of isolated points ***********/

    //vector<Point> jPoint;
    vector<AlphaThickSegmentComputer2DD> blurredSegmentTS;
    vector<PointD>::const_iterator it_MP = MP.begin();
    for(vector<bool>::const_iterator it = isolatedVector.begin(); it != isolatedVector.end(); it++)
    {
        if((*it)==true)//a segment
        {
            it_MP++;
        }
        else//an arc
        {
            AlphaThickSegmentComputer2DD aSegment;
            aSegment.init(it_MP,thickness);
            while(it != isolatedVector.end() && (*it)==false &&
                  aSegment.end() != MP.end() && aSegment.extendFront())
            {
                it++;
                it_MP++;
            }
            if(aSegment.getNumberSegmentPoints() > 1) //FIXME !!!
            {
                //if(verbose)
                //    cout<<"====> BP : "<<getStartPoint(aSegment)<<", EP : "<<getEndPoint(aSegment)<<" and MP.back() : "<<MP.back()<<endl;
                blurredSegmentTS.push_back(aSegment);
                if(getEndPoint(aSegment) != MP.back())
                {
                    it_MP--;
                    //jPoint.push_back(Point(indexDP.at(findElement(MP,getEndPoint(aSegment)) );
                }
                else
                    break;
            }
        }

    }

    if(filename != NULL)
    {
        Board2D TangentSpaceBoard;
        // Draw Middle points in Tangent Space
        for (vector<PointD>::const_iterator it = MP.begin(); it != MP.end();it++)
        {
            TangentSpaceBoard << SetMode((*it).className(), "Both");
            TangentSpaceBoard << *it;
        }

        // Draw Segments in Tangent Space
        int id = 1;
        for (vector<AlphaThickSegmentComputer2DD>::const_iterator it = blurredSegmentTS.begin();it != blurredSegmentTS.end();it++)
        {
            TangentSpaceBoard << SetMode((*it).className(), "BoundingBox");
            TangentSpaceBoard << *it;
            if(verbose)
            {
                double mu, nu;
                PointD N;
                (*it).computeParalellStripParams(mu,N,nu);
                cout<<"Segs "<<id<<" : ("<<N[0]<<" , "<<N[1]<<" , "<<mu<<" , "<<nu<<" ) ==> ";
                cout<<"R = "<<fabs(N[1]/N[0])<<endl;//slope = -N0/N1
                //cout<<getStartPoint(*it)<<" and "<<getEndPoint(*it)<<endl;
                id++;
            }
        }
        TangentSpaceBoard.saveSVG(filename);
    }
}



void drawDecompositionSegmentCircleV1(const vector<Point>& aContour, const vector<Point>& DP,const vector<int> isolatedPoint, const vector<Point>& segments, const vector<Point>& arcs, char* filename)
{
    Board2D DecSAboard;
    /// Display contour points
    for (vector<Point>::const_iterator it = aContour.begin();it != aContour.end();it++)
        DecSAboard << SetMode((*it).className(), "Both") << *it;
    /// Display contour points
    /// Display dominant points
    for (vector<Point>::const_iterator it = DP.begin();it != DP.end();it++)
        DecSAboard << SetMode("PointVector", "Both")
                   << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,0,142)) )//blue
                   << *it;
    /// Display dominant points
    /// Display decoposition (segments)
    DecSAboard.setPenColor(Color( 255, 0, 0));
    DecSAboard.setLineWidth(10.0);
    for (vector<Point>::const_iterator it = segments.begin(); it != segments.end();it++)
        DecSAboard.drawLine(aContour.at((*it)[0])[0],aContour.at((*it)[0])[1],aContour.at((*it)[1])[0],aContour.at((*it)[1])[1],1);
    /// Display decoposition (segments)
    /// Display decoposition (cirles)
    for (int it = 0; it < arcs.size();it++)
    {
        int idBegin = arcs.at(it)[0];
        int idEnd = arcs.at(it)[1];
        int idMid = findBestFittingCircle(aContour,idBegin,idEnd);
        ///int idMid = (idBegin + idEnd)/2;
        bool neg = isLeft(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd))<0;
        Point center = determineCenter(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd));
        double radius = (determineRadius(center,aContour.at(idBegin)) + determineRadius(center,aContour.at(idMid)) + determineRadius(center,aContour.at(idEnd)))/3.0;
        double minRadius = findMinRaidusCircle(center,radius,aContour,idBegin,idEnd);
        double maxRadius = findMaxRaidusCircle(center,radius,aContour,idBegin,idEnd);
        cout<<"Radius = "<<radius<<" min radius = "<<minRadius<<" max radius = "<<maxRadius<<" thkness = "<< (maxRadius-minRadius)/2<<endl;
        double startAngle = atan2(aContour.at(idBegin)[1]-center[1], aContour.at(idBegin)[0]-center[0]);
        double endAngle = atan2(aContour.at(idEnd)[1]-center[1], aContour.at(idEnd)[0]-center[0]);
        //DecSAboard.drawArc(center[0],center[1],radius,0,M_PI/4,true);
        DecSAboard.setPenColor(Color( 0, 255, 0));
        DecSAboard.drawArc(center[0],center[1],radius,startAngle,endAngle,neg);
        DecSAboard.drawArc(center[0],center[1],minRadius,startAngle,endAngle,neg);
        DecSAboard.drawArc(center[0],center[1],maxRadius,startAngle,endAngle,neg);
    }
    //FIXME : Post-processing with ise/lmax ...
    /// Display decoposition (cirles)
    /// Display Isolated + Jonction Point
    for(int it = 0; it<isolatedPoint.size(); it++)
    {
        if(isolatedPoint.at(it)==-1)//JOINCTION = ARC
            DecSAboard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(255,0,0), Color(142,0,0)) )//red
                       << DP.at(it);
        if(isolatedPoint.at(it)==1)//ISOLATED = SEG
            DecSAboard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,142,0)) )//green
                       << DP.at(it) << DP.at(it+1);
    }
    /// Display Isolated + Jonction Point
    DecSAboard.saveSVG(filename);
}

void drawDecompositionSegmentCircleV2(const vector<Point>& aContour, const vector<Point>& DP, const vector<int> isolatedPoint, const vector<Point>& segments, const vector<Point>& arcs, const vector<double> radiusTangentSpace, double thickness, double thickness_MP, char* filename)
{
    Board2D DecSAboard;
    /// Display contour points
    for (vector<Point>::const_iterator it = aContour.begin();it != aContour.end();it++)
        DecSAboard << SetMode((*it).className(), "Both") << *it;
    /// Display contour points
    /// Display dominant points
    for (vector<Point>::const_iterator it = DP.begin();it != DP.end();it++)
        DecSAboard << SetMode("PointVector", "Both")
                   << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,0,142)) )//blue
                   << *it;
    /// Display dominant points
    /// Display decoposition (segments)
    DecSAboard.setPenColor(Color( 255, 0, 0));
    DecSAboard.setLineWidth(10.0);
    for (vector<Point>::const_iterator it = segments.begin(); it != segments.end();it++)
        DecSAboard.drawLine(aContour.at((*it)[0])[0],aContour.at((*it)[0])[1],aContour.at((*it)[1])[0],aContour.at((*it)[1])[1],1);
    /// Display decoposition (segments)
    /// Display decoposition (cirles)
    for (int it = 0; it < arcs.size();it++)
    {
        int idBegin = arcs.at(it)[0];
        int idEnd = arcs.at(it)[1];
        int idMid = (idBegin + idEnd)/2;
        bool neg = isLeft(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd))<0;
        double radius = radiusTangentSpace.at(it);
        Point center = determineCenter(aContour.at(idBegin),aContour.at(idEnd),radius,neg);
        double minRadius = findMinRaidusCircle(center,radius,aContour,idBegin,idEnd);
        double maxRadius = findMaxRaidusCircle(center,radius,aContour,idBegin,idEnd);
        cout<<"Radius = "<<radius<<" min radius = "<<minRadius<<" max radius = "<<maxRadius<<" thkness = "<< (maxRadius-minRadius)/2<<endl;
        double startAngle = atan2(aContour.at(idBegin)[1]-center[1], aContour.at(idBegin)[0]-center[0]);
        double endAngle = atan2(aContour.at(idEnd)[1]-center[1], aContour.at(idEnd)[0]-center[0]);
        //DecSAboard.drawArc(center[0],center[1],radius,0,M_PI/4,true);
        DecSAboard.setPenColor(Color( 255, 255, 0));
        DecSAboard.drawArc(center[0],center[1],radius-(radius*thickness_MP+4*thickness)/2.0,startAngle,endAngle,neg);
        DecSAboard.drawArc(center[0],center[1],radius+(radius*thickness_MP+4*thickness)/2.0,startAngle,endAngle,neg);
    }
    //FIXME : Post-processing with ise/lmax ...
    /// Display decoposition (cirles)
    /// Display Isolated + Jonction Point
    for(int it = 0; it<isolatedPoint.size(); it++)
    {
        if(isolatedPoint.at(it)==-1)//JOINCTION = ARC
            DecSAboard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(255,0,0), Color(142,0,0)) )//red
                       << DP.at(it);
        if(isolatedPoint.at(it)==1)//ISOLATED = SEG
            DecSAboard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,142,0)) )//green
                       << DP.at(it) << DP.at(it+1);
    }
    /// Display Isolated + Jonction Point
    DecSAboard.saveSVG(filename);
}

void drawDecompositionSegmentCircleV3(const vector<Point>& aContour, const vector<Point>& DP,const vector<int> isolatedPoint, const vector<Point>& segments, const vector<Point>& arcs, double thickness, char* filename)
{
    Board2D DecSAboard;
    /// Display contour points
    for (vector<Point>::const_iterator it = aContour.begin();it != aContour.end();it++)
        DecSAboard << SetMode((*it).className(), "Both") << *it;
    /// Display contour points
    /// Display dominant points
    for (vector<Point>::const_iterator it = DP.begin();it != DP.end();it++)
        DecSAboard << SetMode("PointVector", "Both")
                   << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,0,142)) )//blue
                   << *it;
    /// Display dominant points
    /// Display decoposition (segments)
    DecSAboard.setPenColor(Color( 255, 0, 0));
    DecSAboard.setLineWidth(10.0);
    for (vector<Point>::const_iterator it = segments.begin(); it != segments.end();it++)
        DecSAboard.drawLine(aContour.at((*it)[0])[0],aContour.at((*it)[0])[1],aContour.at((*it)[1])[0],aContour.at((*it)[1])[1],1);
    /// Display decoposition (segments)
    /// Display decoposition (cirles)
    for (int it = 0; it < arcs.size();it++)
    {
        int idBegin = arcs.at(it)[0];
        int idEnd = arcs.at(it)[1];
        //find best fitting circle
        int idMid = findBestFittingCircle(aContour,idBegin,idEnd);
        bool neg = isLeft(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd))<0;
        //DecSAboard.drawLine(aContour.at(idBegin)[0],aContour.at(idBegin)[1],aContour.at(idEnd)[0],aContour.at(idEnd)[1],1);
        Point center = determineCenter(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd));
        double radius = (determineRadius(center,aContour.at(idBegin)) + determineRadius(center,aContour.at(idMid)) + determineRadius(center,aContour.at(idEnd)))/3.0;
        cout<<"Radius = "<<radius<<" min radius = "<<radius-(sqrt(2)*thickness)<<" max radius = "<<radius+(sqrt(2)*thickness)<<" thkness = "<< sqrt(2)*thickness/2<<endl;
        double startAngle = atan2(aContour.at(idBegin)[1]-center[1], aContour.at(idBegin)[0]-center[0]);
        double endAngle = atan2(aContour.at(idEnd)[1]-center[1], aContour.at(idEnd)[0]-center[0]);
        DecSAboard.setPenColor(Color( 0, 0, 255));
        DecSAboard.drawArc(center[0],center[1],radius-(sqrt(2)*thickness),startAngle,endAngle,neg);
        DecSAboard.drawArc(center[0],center[1],radius+(sqrt(2)*thickness),startAngle,endAngle,neg);
    }
    //FIXME : Post-processing with ise/lmax ...
    /// Display decoposition (cirles)
    /// Display Isolated + Jonction Point
    for(int it = 0; it<isolatedPoint.size(); it++)
    {
        if(isolatedPoint.at(it)==-1)//JOINCTION = ARC
            DecSAboard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(255,0,0), Color(142,0,0)) )//red
                       << DP.at(it);
        if(isolatedPoint.at(it)==1)//ISOLATED = SEG
            DecSAboard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,142,0)) )//green
                       << DP.at(it) << DP.at(it+1);
    }
    /// Display Isolated + Jonction Point
    DecSAboard.saveSVG(filename);
}

void drawDecompositionSegmentCircleV4(const vector<Point>& aContour, const vector<Point>& DP,const vector<int> isolatedPoint, const vector<Point>& segments, const vector<Point>& arcs, char* filename)
{
    Board2D DecSAboard;
    /// Display contour points
    for (vector<Point>::const_iterator it = aContour.begin();it != aContour.end();it++)
        DecSAboard << SetMode((*it).className(), "Both") << *it;
    /// Display contour points
    /// Display dominant points
    for (vector<Point>::const_iterator it = DP.begin();it != DP.end();it++)
        DecSAboard << SetMode("PointVector", "Both")
                   << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,0,142)) )//blue
                   << *it;
    /// Display dominant points
    /// Display decoposition (segments)
    DecSAboard.setPenColor(Color( 255, 0, 0));
    DecSAboard.setLineWidth(10.0);
    for (vector<Point>::const_iterator it = segments.begin(); it != segments.end();it++)
        DecSAboard.drawLine(aContour.at((*it)[0])[0],aContour.at((*it)[0])[1],aContour.at((*it)[1])[0],aContour.at((*it)[1])[1],1);
    /// Display decoposition (segments)
    /// Display decoposition (cirles)
    for (int it = 0; it<arcs.size(); it++)
    {
        int idBegin = arcs.at(it)[0];
        int idEnd = arcs.at(it)[1];
        Point tmp = findBestChrodCircle(aContour, idBegin, idEnd);
        //int newIdBegin = tmp[0];
        //int newIdEnd = tmp[1];
        idBegin = tmp[0];
        idEnd = tmp[1];
        int idMid = findBestFittingCircle(aContour,idBegin,idEnd);
        bool neg = isLeft(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd))<0;
        Point center = determineCenter(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd));
        double radius = (determineRadius(center,aContour.at(idBegin)) + determineRadius(center,aContour.at(idMid)) + determineRadius(center,aContour.at(idEnd)))/3.0;
        double minRadius = findMinRaidusCircle(center,radius,aContour,idBegin,idEnd);
        double maxRadius = findMaxRaidusCircle(center,radius,aContour,idBegin,idEnd);
        cout<<"Radius = "<<radius<<" min radius = "<<minRadius<<" max radius = "<<maxRadius<<" thkness = "<< (maxRadius-minRadius)/2<<endl;
        double startAngle = atan2(aContour.at(idBegin)[1]-center[1], aContour.at(idBegin)[0]-center[0]);
        double endAngle = atan2(aContour.at(idEnd)[1]-center[1], aContour.at(idEnd)[0]-center[0]);
        //DecSAboard.drawArc(center[0],center[1],radius,0,M_PI/4,true);
        DecSAboard.setPenColor(Color( 0, 255, 0));
        DecSAboard.drawArc(center[0],center[1],radius,startAngle,endAngle,neg);
        DecSAboard.drawArc(center[0],center[1],minRadius,startAngle,endAngle,neg);
        DecSAboard.drawArc(center[0],center[1],maxRadius,startAngle,endAngle,neg);
        DecSAboard.setPenColor(Color( 255, 255, 0));
        DecSAboard.drawLine(aContour.at(idBegin)[0],aContour.at(idBegin)[1],aContour.at(idEnd)[0],aContour.at(idEnd)[1],1);
    }
    //FIXME : Post-processing with ise/lmax ...
    /// Display decoposition (cirles)
    /// Display Isolated + Jonction Point
    for(int it = 0; it<isolatedPoint.size(); it++)
    {
        if(isolatedPoint.at(it)==-1)//JOINCTION = ARC
            DecSAboard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(255,0,0), Color(142,0,0)) )//red
                       << DP.at(it);
        if(isolatedPoint.at(it)==1)//ISOLATED = SEG
            DecSAboard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,142,0)) )//green
                       << DP.at(it) << DP.at(it+1);
    }
    /// Display Isolated + Jonction Point
    DecSAboard.saveSVG(filename);
}

void drawDecompositionSegmentCircleAll(const vector<Point>& aContour, const vector<Point>& DP,const vector<int> isolatedPoint, const vector<Point>& segments, const vector<Point>& arcs, double thickness, double thickness_MP, char* filename)
{
    Board2D DecSAboard;
    /// Display contour points
    for (vector<Point>::const_iterator it = aContour.begin();it != aContour.end();it++)
        DecSAboard << SetMode((*it).className(), "Both") << *it;
    /// Display contour points
    /// Display dominant points
    for (vector<Point>::const_iterator it = DP.begin();it != DP.end();it++)
        DecSAboard << SetMode("PointVector", "Both")
                   << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,0,142)) )//blue
                   << *it;
    /// Display dominant points
    /// Display decoposition (segments)
    //DecSAboard.setPenColor(Color( 255, 0, 0));
    DecSAboard.setPenColor(Color( 0, 0, 0));
    DecSAboard.setLineWidth(10.0);
    for (vector<Point>::const_iterator it = segments.begin(); it != segments.end();it++)
        DecSAboard.drawLine(aContour.at((*it)[0])[0],aContour.at((*it)[0])[1],aContour.at((*it)[1])[0],aContour.at((*it)[1])[1],1);
    /// Display decoposition (segments)
    /// Display decoposition (cirles)
    for (int it = 0; it < arcs.size();it++)
    {
        int idBegin = arcs.at(it)[0];
        int idEnd = arcs.at(it)[1];
        //find best fitting circle
        int idMid = findBestFittingCircle(aContour,idBegin,idEnd);
        ///int idMid = (idBegin + idEnd)/2;
        bool neg = isLeft(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd))<0;
        //DecSAboard.drawLine(aContour.at(idBegin)[0],aContour.at(idBegin)[1],aContour.at(idEnd)[0],aContour.at(idEnd)[1],1);
        Point center = determineCenter(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd));
        double radius = (determineRadius(center,aContour.at(idBegin)) + determineRadius(center,aContour.at(idMid)) + determineRadius(center,aContour.at(idEnd)))/3.0;
        ///double radius = radiusTangentSpace.at(it);
        ///Point center = determineCenter(aContour.at(idBegin),aContour.at(idEnd),radius,neg);
        double minRadius = findMinRaidusCircle(center,radius,aContour,idBegin,idEnd);
        double maxRadius = findMaxRaidusCircle(center,radius,aContour,idBegin,idEnd);
        cout<<"Radius = "<<radius<<" min radius = "<<minRadius<<" max radius = "<<maxRadius<<" thkness = "<< (maxRadius-minRadius)/2<<endl;
        double startAngle = atan2(aContour.at(idBegin)[1]-center[1], aContour.at(idBegin)[0]-center[0]);
        double endAngle = atan2(aContour.at(idEnd)[1]-center[1], aContour.at(idEnd)[0]-center[0]);
        //DecSAboard.drawArc(center[0],center[1],radius,0,M_PI/4,true);
        DecSAboard.setPenColor(Color( 0, 255, 0));
        DecSAboard.drawArc(center[0],center[1],radius,startAngle,endAngle,neg);
        DecSAboard.drawArc(center[0],center[1],minRadius,startAngle,endAngle,neg);
        DecSAboard.drawArc(center[0],center[1],maxRadius,startAngle,endAngle,neg);

        DecSAboard.setPenColor(Color( 255, 255, 0));
        DecSAboard.drawArc(center[0],center[1],radius-(radius*thickness_MP+4*thickness)/2.0,startAngle,endAngle,neg);
        DecSAboard.drawArc(center[0],center[1],radius+(radius*thickness_MP+4*thickness)/2.0,startAngle,endAngle,neg);

        DecSAboard.setPenColor(Color( 0, 0, 255));
        DecSAboard.drawArc(center[0],center[1],radius-(sqrt(2)*thickness),startAngle,endAngle,neg);
        DecSAboard.drawArc(center[0],center[1],radius+(sqrt(2)*thickness),startAngle,endAngle,neg);

        //DecSAboard.drawCircle(center[0],center[1],radius);
        //DecSAboard.fillCircle(center[0],center[1],1);//for center
        /*
        DecSAboard.setPenColor(Color( 255, 255, 0));
        DecSAboard.drawLine(aContour.at(idBegin)[0],aContour.at(idBegin)[1],aContour.at(idEnd)[0],aContour.at(idEnd)[1],1);
        */
    }
    //FIXME : Post-processing with ise/lmax ...
    /// Display decoposition (cirles)
    /// Display Isolated + Jonction Point
    for(int it = 0; it<isolatedPoint.size(); it++)
    {
        if(isolatedPoint.at(it)==-1)//JOINCTION = ARC
            DecSAboard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(255,0,0), Color(142,0,0)) )//red
                       << DP.at(it);
        if(isolatedPoint.at(it)==1)//ISOLATED = SEG
            DecSAboard << SetMode("PointVector", "Both")
                       << CustomStyle("PointVector", new CustomColors( Color(0,255,0), Color(0,142,0)) )//green
                       << DP.at(it) << DP.at(it+1);
    }
    /// Display Isolated + Jonction Point
    DecSAboard.saveSVG(filename);

    /*
    if(IS_PLOT)
    {
        // Calcul radius function at points in common zone
        for (int it = 0; it < arcs.size();it++)
        {
            int idBegin = arcs.at(it)[0];
            int idEnd = arcs.at(it)[1];

            int oneThird = (idEnd-idBegin)/3;
            double ise,radius,minIse = -1;
            int idMin = -1;
            Point center;
            for(int idMid = idStart+oneThird ; idMid < (idEnd-oneThird); idMid++)
            {
                center = determineCenter(aContour.at(idStart),aContour.at(idMid),aContour.at(idEnd));
                radius = (determineRadius(center,aContour.at(idStart)) + determineRadius(center,aContour.at(idMid)) + determineRadius(center,aContour.at(idEnd)))/3.0;
                ise = iseContourCircle(aContour,idStart,idEnd,center,radius);
                if((minIse<0) || (ise<minIse))
                {
                    minIse = ise;
                    idMin = idMid;
                }
            }
        }
        for (vector<AlphaThickSegmentComputer2D>::const_iterator it = fuzzySegmentSet.begin();it != fuzzySegmentSet.end();it++)
        {
            idBegin.push_back(findElement(aContour,getStartPoint(*it)));
            idEnd.push_back(findElement(aContour,getEndPoint(*it)));
        }
        writeFile(idBegin,"idBeginSeg.txt",false);
        writeFile(idEnd,"idEndSeg.txt",false);
    }
    */
}


void testVerifyCircularityChrodProp(const vector<Point>& aContour, const vector<Point>& DP, const vector<int>& indexDP, double thickness, vector<Point>& segments, vector<Point>& arcs, bool verbose)
{
    double distanceAB=0, angleLeft, angleRight, angleDiff, angleMax = -1, angleLimite;
    int midPoint,oneThird;
    for(int it = 0; it<DP.size()-2; it++)
    {
        //if(it==3)
        //    cout<<"STOP"<<endl;
        //calcul To_s
        //cout<<"indexDP.at(it) = "<<indexDP.at(it)<<" and indexDP.at(it+1) = "<<indexDP.at(it+1)<<endl;
        distanceAB = distancePoints(DP.at(it),DP.at(it+2));
        angleLimite = 2*atan(2*thickness/distanceAB);
        midPoint = (indexDP.at(it+2) + indexDP.at(it))/2;
        oneThird = (indexDP.at(it+2) - indexDP.at(it))/3;
        for(int idP = indexDP.at(it)+oneThird; idP < indexDP.at(it+2)-oneThird+1; idP++)
        {
            //angleLeft = atan(thickness/distancePoints(DP.at(it),aContour.at(idP)));
            //angleRight = atan(thickness/distancePoints(aContour.at(idP),DP.at(it+1)));
            angleLeft = acuteAngle(DP.at(it),aContour.at(idP),DP.at(it+2));
            angleRight = acuteAngle(DP.at(it),aContour.at(midPoint),DP.at(it+2));
            angleDiff = fabs(angleLeft-angleRight);
            cout<<"angleLeft = "<<angleLeft<<" and angleRight = "<<angleRight<<" and angleDiff = "<<angleDiff<<endl;
            if(angleDiff>angleMax)
                angleMax = angleDiff;
        }
        if(angleMax<=angleLimite)//ARC
        {
            if(verbose)
                cout<<" ARC : ";
            arcs.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
            arcs.push_back(Point(indexDP.at(it+1),indexDP.at(it+2)));
        }
        else//SEG
        {
            if(verbose)
                cout<<" SEG : ";
            segments.push_back(Point(indexDP.at(it),indexDP.at(it+1)));
            segments.push_back(Point(indexDP.at(it+1),indexDP.at(it+2)));
        }
        if(verbose)
            cout<<"angleMax = "<<angleMax<<" angleLimite = "<<angleLimite<<endl;
    }
}

/*****************************************************/
/****Decomposition of Curve into Segments and Arcs ***/
/*****************************************************/

/********************************************************/
/********** Test circle simple with DP only *************/
/********************************************************/
void testCircleSimpleV1(const vector<Point>& aContour, const vector<Point>& DP, const char* filename, bool verbose)//draw approximated circle
{
    //Color white( 255, 255, 255 );
    //Color black( 0, 0, 0 );
    Color red( 255, 0, 0 );
    Color dred( 192, 0, 0 );
    Color green( 0, 255, 0 );
    Color dgreen( 0, 192, 0 );
    //Color blue( 0, 0, 255 );
    //Color dblue( 0, 0, 192 );
    Board2D boardAll;
    boardAll.setPenColor(Color( 255, 0, 0));
    boardAll.setLineWidth(15.0);
    vector<Point> centerVector;
    for (int it = 1; it+1 < DP.size();it++)
    {
        //if(it==42)
        //    cout<<"STOP"<<endl;
        Point center = determineCenter(DP.at(it-1),DP.at(it),DP.at(it+1));
        double radius = (determineRadius(center,DP.at(it-1)) + determineRadius(center,DP.at(it)) + determineRadius(center,DP.at(it+1)))/3.0;
        if(verbose)
        {
            cout<<"center = "<<center<<" and radius = "<<radius<<endl;
            cout<<"length DP segment = "<<distancePoints(DP.at(it-1),DP.at(it)) + distancePoints(DP.at(it),DP.at(it+1))
                <<", length segment = "<<lengthContour(aContour,DP.at(it-1),DP.at(it)) + lengthContour(aContour,DP.at(it),DP.at(it+1))
                <<", and length arc = "<<arcLength(DP.at(it-1),DP.at(it),DP.at(it+1))<<endl;
        }
        //FIXME : real length segment !!!!
        if(filename != NULL)
        {
            //boardAll.drawCircle(center[0],center[1],radius);
            double endAngle = atan2(DP.at(it-1)[1]-center[1], DP.at(it-1)[0]-center[0]);
            double startAngle = atan2(DP.at(it+1)[1]-center[1], DP.at(it+1)[0]-center[0]);
            //boardAll.drawArc(center[0],center[1],radius,0,M_PI/4,true);
            bool neg = isLeft(DP.at(it),DP.at(it-1),DP.at(it+1))<0;
            //if(it==42)
                boardAll.drawArc(center[0],center[1],radius,startAngle,endAngle,neg);
                //boardAll.drawCircle(center[0],center[1],radius);
            //boardAll.fillCircle(center[0],center[1],1);//for center
        }
        centerVector.push_back(center);
    }
    if(filename != NULL)
    {
        for (vector<Point>::const_iterator it = aContour.begin(); it != aContour.end();it++)
            boardAll << CustomStyle( (*it).className(), new MyDrawStyleCustomColor( green, dgreen ) )
                     << *it;
        int count=0;
        for (vector<Point>::const_iterator it = DP.begin(); it != DP.end();it++)
        {
            //if(count>=41 && count<=43)
                boardAll << CustomStyle( (*it).className(), new MyDrawStyleCustomColor( red, dred ) )
                     << *it;
            count++;
        }
        boardAll.saveSVG(filename);
    }
}

void testCircleSimpleV2(const vector<Point>& aContour, const vector<Point>& DP, const char* filename, bool verbose)//with errors ise/lmax
{
    vector<Point> centerVector;
    vector<double> radiusVector;
    Color red( 255, 0, 0 );
    Color dred( 192, 0, 0 );
    Color green( 0, 255, 0 );
    Color dgreen( 0, 192, 0 );
    Board2D boardAll;
    boardAll.setPenColor(Color( 255, 0, 0));
    boardAll.setLineWidth(15.0);

    for (int it = 1; it < DP.size()-1; it++)
    {
        Point center = determineCenter(DP.at(it-1),DP.at(it),DP.at(it+1));
        double radius = (determineRadius(center,DP.at(it-1)) + determineRadius(center,DP.at(it)) + determineRadius(center,DP.at(it+1)))/3.0;
        if(verbose)
        {
            cout<<it-1<<" : center = "<<center<<" and radius = "<<radius<<endl;
            cout<<"length DP segment = "<<distancePoints(DP.at(it-1),DP.at(it)) + distancePoints(DP.at(it),DP.at(it+1))
                <<", length segment = "<< lengthContour(aContour,DP.at(it-1),DP.at(it)) + lengthContour(aContour,DP.at(it),DP.at(it+1))
                <<", and length arc = "<< arcLength(DP.at(it-1),DP.at(it),DP.at(it+1)) <<endl;
        }
        double ise = iseContourCircle(aContour,DP.at(it-1),DP.at(it),center,radius) + iseContourCircle(aContour,DP.at(it),DP.at(it+1),center,radius);
        ///double lmax = lmaxContourCircle(aContour,DP.at(it-1),DP.at(it),center,radius) + lmaxContourCircle(aContour,DP.at(it),DP.at(it+1),center,radius);
        ///if(lmax<4 && radius<50) // && radius>40
        if(ise<100) // && radius>40
        {
            if(filename != NULL)
            {
                //boardAll.drawCircle(center[0],center[1],radius);
                //boardAll.drawArc(center[0],center[1],radius,0,M_PI/4,true);
                bool neg = isLeft(DP.at(it),DP.at(it-1),DP.at(it+1))>0;

                double startAngle = atan2(DP.at(it-1)[1]-center[1], DP.at(it-1)[0]-center[0]);
                double endAngle = atan2(DP.at(it+1)[1]-center[1], DP.at(it+1)[0]-center[0]);
                boardAll.drawArc(center[0],center[1],radius,startAngle,endAngle,neg);
                //boardAll.drawCircle(center[0],center[1],radius);
            }
        }
        //boardAll.fillCircle(center[0],center[1],1);//for center
        centerVector.push_back(center);//of it
        radiusVector.push_back(radius);//of it
    }

    if(filename != NULL)
    {
        for (vector<Point>::const_iterator it = aContour.begin(); it != aContour.end();it++)
            boardAll << CustomStyle( (*it).className(), new MyDrawStyleCustomColor( green, dgreen ) )
                     << *it;
        for (vector<Point>::const_iterator it = DP.begin(); it != DP.end();it++)
            boardAll << CustomStyle( (*it).className(), new MyDrawStyleCustomColor( red, dred ) )
                 << *it;
        boardAll.saveSVG(filename);
    }
}

void testCircleSimpleV3(const vector<Point>& aContour, const vector<Point>& DP, const vector<int>& indexDP, const char* filename, bool verbose)//with approximate radius
{
    Color red( 255, 0, 0 );
    Color dred( 192, 0, 0 );
    Color green( 0, 255, 0 );
    Color dgreen( 0, 192, 0 );
    Board2D boardAll;
    boardAll.setPenColor(Color( 255, 0, 0));
    boardAll.setLineWidth(15.0);
    vector<Point> centerVector;
    vector<double> radiusVector;

    for (int it = 1; it < DP.size()-1; it++)
    {
        Point center = determineCenter(DP.at(it-1),DP.at(it),DP.at(it+1));
        double radius = (determineRadius(center,DP.at(it-1)) + determineRadius(center,DP.at(it)) + determineRadius(center,DP.at(it+1)))/3.0;
        if(verbose)
        {
            cout<<it-1<<" : center = "<<center<<" and radius = "<<radius<<endl;
            cout<<"length DP segment = "<<distancePoints(DP.at(it-1),DP.at(it)) + distancePoints(DP.at(it),DP.at(it+1))
                <<", length segment = "<< lengthContour(aContour,DP.at(it-1),DP.at(it)) + lengthContour(aContour,DP.at(it),DP.at(it+1))
                <<", and length arc = "<< arcLength(DP.at(it-1),DP.at(it),DP.at(it+1)) <<endl;
        }
        centerVector.push_back(center);//of it
        radiusVector.push_back(radius);//of it
    }
    //combine the arcs/circles
    //vector<Point> simpDecSegs;
    vector<Point> simpDecArcs;
    int index = 1, index_prev=1;
    while(index<radiusVector.size())
    {
        index_prev = index;
        while(index<radiusVector.size() && fabs(radiusVector.at(index-1) - radiusVector.at(index))<30)
            index++;
        if(index - index_prev > 3)
        {
            simpDecArcs.push_back(Point(indexDP.at(index_prev-1),indexDP.at(index+1)));
            cout<<"ARC :"<<index_prev-1<<" and "<<index+1<<endl;
        }
        else
        {
            //simpDecSegs.push_back(Point(indexDP.at(index),indexDP.at(index-1)));
            index++;
        }
    }
    if(filename!=NULL)
    {
        for(int i=0; i<simpDecArcs.size(); i++)
        {
            int idBegin = simpDecArcs.at(i)[0];
            int idEnd = simpDecArcs.at(i)[1];
            int idMid = (idBegin+idEnd)/2;
            Point center = determineCenter(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd));
            double radius = (determineRadius(center,aContour.at(idBegin)) + determineRadius(center,aContour.at(idMid)) + determineRadius(center,aContour.at(idEnd)))/3.0;
            double startAngle = atan2(aContour.at(idBegin)[1]-center[1], aContour.at(idBegin)[0]-center[0]);
            double endAngle = atan2(aContour.at(idEnd)[1]-center[1], aContour.at(idEnd)[0]-center[0]);
            bool neg = isLeft(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd))<0;
            boardAll.drawArc(center[0],center[1],radius,startAngle,endAngle,neg);
        }

        for (vector<Point>::const_iterator it = aContour.begin(); it != aContour.end();it++)
            boardAll << CustomStyle( (*it).className(), new MyDrawStyleCustomColor( green, dgreen ) )
                     << *it;
        for (vector<Point>::const_iterator it = DP.begin(); it != DP.end();it++)
            boardAll << CustomStyle( (*it).className(), new MyDrawStyleCustomColor( red, dred ))
                 << *it;
        boardAll.saveSVG(filename);
    }
}

int findBestFittingCircle(const vector<Point> aContour, int idStart, int idEnd)
{
    int oneThird = (idEnd-idStart)/3;
    double ise,radius,minIse = -1;
    int idMin = -1;
    Point center;
    for(int idMid = idStart+oneThird ; idMid < (idEnd-oneThird); idMid++)
    {
        center = determineCenter(aContour.at(idStart),aContour.at(idMid),aContour.at(idEnd));
        radius = (determineRadius(center,aContour.at(idStart)) + determineRadius(center,aContour.at(idMid)) + determineRadius(center,aContour.at(idEnd)))/3.0;
        ise = iseContourCircle(aContour,idStart,idEnd,center,radius);
        if((minIse<0) || (ise<minIse))
        {
            minIse = ise;
            idMin = idMid;
        }
    }
    return idMin;
}

Point findBestChrodCircle(const vector<Point> aContour, int idBegin, int idEnd){
    int oneThird = (idEnd-idBegin)/3;
    bool neg = isLeft(aContour.at(idBegin),aContour.at((idEnd+idBegin)/2),aContour.at(idEnd))<0;
    //verif on the left
    int idLeft = idBegin;
    for(int idMid = idBegin; idMid < idBegin+oneThird; idMid++)
    {
        bool negL = isLeft(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd))<0;
        if(neg != negL)
        {
            idLeft = idMid;
            idBegin = idMid;
        }
    }
    //verif on the right
    int idRight = idEnd;
    for(int idMid = idEnd; idMid > idEnd-oneThird; idMid--)
    {
        bool negR = isLeft(aContour.at(idBegin),aContour.at(idMid),aContour.at(idEnd))<0;
        if(neg != negR)
        {
            idRight = idMid;
            idEnd = idMid;
        }
    }
    return Point(idLeft,idRight);
}

int findBestCircle(const vector<Point> aContour, int idStart, int idEnd)
{
    int oneFith = 1;//(idEnd-idStart)/5;
    double minDistance = distancePoints(aContour.at(idStart),aContour.at(idEnd));
    int idMin = -1;
    ///PointD center = PointD((aContour.at(idStart)[0]+aContour.at(idEnd)[0])/2.0,(aContour.at(idStart)[1]+aContour.at(idEnd)[1])/2.0);
    ///for(int idMid = idStart+oneThird ; idMid < (idEnd-oneThird); idMid++)
    for(int idMid = idStart+oneFith ; idMid < (idEnd-oneFith); idMid++)
    {
        //PointD tmp = PointD(double(aContour.at(idMid)[0]),double(aContour.at(idMid)[1]));
        ///double distance = distancePoints(tmp,center);
        double distance = distancePointSegment(aContour.at(idMid),aContour.at(idStart),aContour.at(idEnd));
        if(distance<minDistance)
        {
            idMin = idMid;
            minDistance = distance;
        }
    }
    return idMin;
}

double findMaxRaidusCircle(Point center, double radius, const vector<Point> aContour, int idStart, int idEnd)
{
    double maxRadius=0,distance=0;
    maxRadius = signedDistancePointCircle(aContour.at(idStart),center,radius);
    for(int idMid = idStart+1; idMid <= idEnd; idMid++)
    {
        distance = signedDistancePointCircle(aContour.at(idMid),center,radius);
        if(distance>0 && distance>maxRadius)
            maxRadius = distance;
    }
    return radius + maxRadius;
}

double findMinRaidusCircle(Point center, double radius, const vector<Point> aContour, int idStart, int idEnd)
{
    double minRadius=0,distance=0;
    minRadius = signedDistancePointCircle(aContour.at(idStart),center,radius);
    for(int idMid = idStart+1; idMid <= idEnd; idMid++)
    {
        distance = signedDistancePointCircle(aContour.at(idMid),center,radius);
        if(distance<0 && distance<minRadius)
            minRadius = distance;
    }
    return fabs(radius + minRadius);
}

/********************************************************/
/********** Test circle simple with DP only *************/
/********************************************************/
