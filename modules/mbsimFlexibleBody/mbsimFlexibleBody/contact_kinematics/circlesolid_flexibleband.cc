/* Copyright (C) 2004-2011 MBSim Development Team
 *
 * This library is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU Lesser General Public 
 * License as published by the Free Software Foundation; either 
 * version 2.1 of the License, or (at your option) any later version. 
 *  
 * This library is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
 * Lesser General Public License for more details. 
 *  
 * You should have received a copy of the GNU Lesser General Public 
 * License along with this library; if not, write to the Free Software 
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA

 * Contact: thorsten.schindler@mytum.de
 *          rzander@users.berlios.de
 
 */

#include<config.h>
#include "mbsimFlexibleBody/contact_kinematics/circlesolid_flexibleband.h"
#include "mbsimFlexibleBody/contours/flexible_band.h"
#include "mbsim/contour.h"
#include "mbsim/contours/circle_solid.h"
#include "mbsim/contours/point.h"
#include "mbsim/functions_contact.h"

using namespace std;
using namespace fmatvec;
using namespace MBSim;

namespace MBSimFlexibleBody { 

  ContactKinematicsCircleSolidFlexibleBand::ContactKinematicsCircleSolidFlexibleBand() : ContactKinematics(),icircle(0),icontour(0),possibleContactsPerNode(1),circle(0),band(0) {}
  ContactKinematicsCircleSolidFlexibleBand::~ContactKinematicsCircleSolidFlexibleBand() {}

  void ContactKinematicsCircleSolidFlexibleBand::assignContours(const vector<Contour*>& contour) {	
    if(dynamic_cast<CircleSolid*>(contour[0])) {
      icircle = 0;
      icontour = 1;
      circle = static_cast<CircleSolid*>(contour[0]);
      band = static_cast<FlexibleBand*>(contour[1]);
    }
    else {
      icircle = 1;
      icontour = 0;
      circle = static_cast<CircleSolid*>(contour[1]);
      band = static_cast<FlexibleBand*>(contour[0]);
    }
    wBand   = band->getWidth();
    hBand   = band->getNormalDistance();
    rCircle = circle->getRadius();

    staticNodes                    =     band->getNodes();
    numberOfPotentialContactPoints = 2 * possibleContactsPerNode*band->getNodes().size()-1;  // dies braeuchte einen eigenen init-Call
    l0     = 1.0   * fabs(band->getAlphaEnd()-band->getAlphaStart())/staticNodes.size(); /* bandwidth of mesh deformer: higher values leads to stronger attraction of last contact points */
    epsTol = 5.e-2 * l0; /* distance, when two contact points should be treated as one */

  }

  //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  void ContactKinematicsCircleSolidFlexibleBand::updateg(std::vector<SingleContact> & contacts) {
    //Vector of lagrange parameters for contact positions
    Vec contactPositions(numberOfPotentialContactPoints,NONINIT);
    int inContact = 0;


    Vec nodes    = staticNodes.copy();
    int nrNodes = nodes.size();

    // always use nodes for contact kinematics
    for(int i=0;i<nrNodes;i++) {
      ContourPointData* cpData = contacts[i].getcpData();
      Vec & g = contacts[i].getg();

      // set lagrange parateter position of current node
      if(cpData[icontour].getLagrangeParameterPosition().size() != 2)
        cpData[icontour].getLagrangeParameterPosition().resize() << VecV(2);
      cpData[icontour].getLagrangeParameterPosition()(0) = nodes(i);
      cpData[icontour].getLagrangeParameterPosition()(1) = 0.0;


      if(cpData[icontour].getLagrangeParameterPosition()(0) < band->getAlphaStart() || cpData[icontour].getLagrangeParameterPosition()(0) > band->getAlphaEnd())
        g(0) = 1.0;
      else {
        band->updateKinematicsForFrame(cpData[icontour],position_cosy);
        Vec Wd = circle->getFrame()->getPosition() - cpData[icontour].getFrameOfReference().getPosition();
        g(0) = nrm2(Wd) - rCircle;

        if(g(0) <= 0.0) {
          contactPositions(inContact) = cpData[icontour].getLagrangeParameterPosition()(0);
          inContact ++;
        }

        {

          Vec Wb = cpData[icontour].getFrameOfReference().getOrientation().col(2);
          cpData[icontour].getLagrangeParameterPosition()(1) = trans(Wb)*Wd; // get contact parameter of second tangential direction

          if(fabs(cpData[icontour].getLagrangeParameterPosition()(1)) > 0.5*wBand)
            g(0) = 1.0;
          else { // calculate the normal distance
            cpData[icontour].getFrameOfReference().getPosition()           += cpData[icontour].getLagrangeParameterPosition()(1)*Wb;
            cpData[icircle] .getFrameOfReference().getOrientation().set(0, -cpData[icontour].getFrameOfReference().getOrientation().col(0));
            cpData[icircle] .getFrameOfReference().getOrientation().set(1, -cpData[icontour].getFrameOfReference().getOrientation().col(1));
            cpData[icircle] .getFrameOfReference().getOrientation().set(2,  cpData[icontour].getFrameOfReference().getOrientation().col(2));
            //cpData[icircle] .getFrameOfReference().getPosition()            =  circle->getFrame()->getPosition() + cpData[icircle].getFrameOfReference().getOrientation().col(0)*rCircle;
            cpData[icircle] .getFrameOfReference().getPosition()            =  circle->getFrame()->getPosition() - rCircle*Wd/nrm2(Wd);

          }
        }
      }
    }

    FuncPairContour1sCircleSolid *func= new FuncPairContour1sCircleSolid(circle,band); // root function for searching contact parameters
    Contact1sSearch search(func);

#if 0
    // deform mesh
    for(int i=0;i<nodes.size();i++) {
      double node_i = nodes(i);
      for(int j=0;j<lastContactPositions.size();j++) {
        nodes(i) += (lastContactPositions(j)-node_i) * exp( - fabs((lastContactPositions(j)-node_i)/ l0));
      }
    }
    if(result.rows()>0) {
      cout << endl <<endl<< "----------------------------------------------\nlast hits: "<< trans(lastContactPositions) << endl;
      Vec tempVec = staticNodes;
      cout << "original Nodes" << endl << trans( tempVec ) << endl;
      cout << "modified Nodes" << endl << trans( nodes   ) << endl;
      cout << trans(result) << endl;
    }
    if ( lastContactPositions.size() ) throw 1;
#endif

    search.setNodes(nodes); // defining search areas for contacts
    Mat result = search.slvAll();
    delete func;

    for(int i=0;i<result.rows();i++) {
      ContourPointData* cpData = contacts[nrNodes + i].getcpData();
      Vec & gi = contacts[i].getg();
      Vec & gNri = contacts[nrNodes + i].getg();

      if(cpData[icontour].getLagrangeParameterPosition().size() != 2)
        cpData[icontour].getLagrangeParameterPosition().resize() << VecV(2);
      cpData[icontour].getLagrangeParameterPosition()(0) = result(i,0);
      cpData[icontour].getLagrangeParameterPosition()(1) = 0.0;

      if(cpData[icontour].getLagrangeParameterPosition()(0) < band->getAlphaStart() || cpData[icontour].getLagrangeParameterPosition()(0) > band->getAlphaEnd())
        gNri(0) = 1.0;
      else {
        bool doubledNode = false;
        for(int j=0;j<staticNodes.size();j++)
          if(fabs(cpData[icontour].getLagrangeParameterPosition()(0)-staticNodes(j))<epsTol)
            doubledNode = true;

        if(doubledNode)
          gNri(0) = 1.0;
        else {
          band->updateKinematicsForFrame(cpData[icontour],position_cosy);
          Vec Wd = circle->getFrame()->getPosition() - cpData[icontour].getFrameOfReference().getPosition();
          gNri(0) = nrm2(Wd) - rCircle;

          if(gi(0) <= 0.0) {
            contactPositions(inContact) = cpData[icontour].getLagrangeParameterPosition()(0);
            inContact ++;
          }

          {
            Vec Wb = cpData[icontour].getFrameOfReference().getOrientation().col(2);
            cpData[icontour].getLagrangeParameterPosition()(1) = trans(Wb)*Wd; // get contact parameter of second tangential direction

            if(fabs(cpData[icontour].getLagrangeParameterPosition()(1)) > 0.5*wBand)
              gi(0) = 1.0;
            else { // calculate the normal distance
              cpData[icontour].getFrameOfReference().getPosition()           += cpData[icontour].getLagrangeParameterPosition()(1)*Wb;
              cpData[icircle] .getFrameOfReference().getOrientation().set(0, -cpData[icontour].getFrameOfReference().getOrientation().col(0));
              cpData[icircle] .getFrameOfReference().getOrientation().set(1, -cpData[icontour].getFrameOfReference().getOrientation().col(1));
              cpData[icircle] .getFrameOfReference().getOrientation().set(2,  cpData[icontour].getFrameOfReference().getOrientation().col(2));
              cpData[icircle] .getFrameOfReference().getPosition()            =  circle->getFrame()->getPosition() + cpData[icircle].getFrameOfReference().getOrientation().col(0)*rCircle;

#if 0
              if(false && (ig[i])(0) < 0.)//1.0e-4)
              {
                cout << endl;
                cout << " g(0)  = " << (ig[i])(0) << "\n------------------------------------------\n" << endl;
                cout << " " << i << " - alpha = " << cpData[icontour].getLagrangeParameterPosition()(0) << endl;
                cout << "     Wb = " << trans(Wb) << "\n Wd = " << trans(Wd) << endl;
                cout << "     rCircle  = " << rCircle << endl;
                cout << " cpData[icontour].r = \n" << trans(cpData[icontour].getFrameOfReference().getPosition()   ) << endl;
                cout << " cpData[icontour].A = \n" <<       cpData[icontour].getFrameOfReference().getOrientation()  << endl;
                cout << " cpData[icircle].r  = \n" << trans(cpData[icircle] .getFrameOfReference().getPosition()   ) << endl;
                cout << " cpData[icircle].A  = \n" <<       cpData[icircle] .getFrameOfReference().getOrientation()  << endl;
                cpData[icontour].getLagrangeParameterPosition()(0) = 0.0;
                band->updateKinematicsForFrame(cpData[icontour],position_cosy);
                cout << " cpData(alpha=0).r = \n" << trans(cpData[icontour].getFrameOfReference().getPosition()   ) << endl;
                cout << " cpData(alpha=0).A = \n" <<       cpData[icontour].getFrameOfReference().getOrientation()  << endl;
                throw 1;
              }
#endif
            }
          }
        }
      }
    }
    lastContactPositions.resize() = contactPositions(Index(0,inContact-1));

    // flush data not being necessary
    for(int i=nrNodes + result.rows();i<numberOfPotentialContactPoints;i++)
    {
      ContourPointData* cpData = contacts[i].getcpData();

      // positions -> dummy and centerOfCircle
      if(cpData[icontour].getLagrangeParameterPosition().size() != 2)
        cpData[icontour].getLagrangeParameterPosition().resize() << VecV(2);
      cpData[icontour].getLagrangeParameterPosition()(0) = 0.0;

      cpData[icontour].getFrameOfReference().setPosition(nodes(nodes.size()-1)*Vec("[1.0;0.0;0.0]"));//Vec(3,INIT,0.0));//
      cpData[icircle ].getFrameOfReference().setPosition(nodes(nodes.size()-1)*Vec("[1.0;0.0;0.0]"));//Vec(3,INIT,0.0));//
      cpData[icontour].getFrameOfReference().getOrientation() = -SqrMat(3,EYE);
      cpData[icircle ].getFrameOfReference().getOrientation() = -SqrMat(3,EYE);
      Vec & gi = contacts[i].getg();
      gi(0) = 1.0;
    }

  }

  void ContactKinematicsCircleSolidFlexibleBand::updateg(std::vector<fmatvec::Vec>::iterator ig, std::vector<ContourPointData*>::iterator icpData) {
    Vec contactPositions(numberOfPotentialContactPoints,NONINIT);
    int inContact = 0;

    Vec nodes    = staticNodes.copy();
    int nrNodes = nodes.size();

    // always use nodes for contact kinematics
    for(int i=0;i<nrNodes;i++) {
      ContourPointData* cpData = icpData[i];

      if(cpData[icontour].getLagrangeParameterPosition().size() != 2)
        cpData[icontour].getLagrangeParameterPosition().resize() = Vec(2,NONINIT);
      cpData[icontour].getLagrangeParameterPosition()(0) = nodes(i);
      cpData[icontour].getLagrangeParameterPosition()(1) = 0.0;

      if(cpData[icontour].getLagrangeParameterPosition()(0) < band->getAlphaStart() || cpData[icontour].getLagrangeParameterPosition()(0) > band->getAlphaEnd())
        (ig[i])(0) = 1.0;
      else {
        band->updateKinematicsForFrame(cpData[icontour],position_cosy);
        Vec Wd = circle->getFrame()->getPosition() - cpData[icontour].getFrameOfReference().getPosition();
        (ig[i])(0) = nrm2(Wd) - rCircle;

        if((ig[i])(0) <= 0.0) {
          contactPositions(inContact) = cpData[icontour].getLagrangeParameterPosition()(0);
          inContact ++;
        }

        {

          Vec Wb = cpData[icontour].getFrameOfReference().getOrientation().col(2);
          cpData[icontour].getLagrangeParameterPosition()(1) = trans(Wb)*Wd; // get contact parameter of second tangential direction

          if(fabs(cpData[icontour].getLagrangeParameterPosition()(1)) > 0.5*wBand)
            (ig[i])(0) = 1.0;
          else { // calculate the normal distance
            cpData[icontour].getFrameOfReference().getPosition()           += cpData[icontour].getLagrangeParameterPosition()(1)*Wb;
            cpData[icircle] .getFrameOfReference().getOrientation().set(0, -cpData[icontour].getFrameOfReference().getOrientation().col(0));
            cpData[icircle] .getFrameOfReference().getOrientation().set(1, -cpData[icontour].getFrameOfReference().getOrientation().col(1));
            cpData[icircle] .getFrameOfReference().getOrientation().set(2,  cpData[icontour].getFrameOfReference().getOrientation().col(2));
            //cpData[icircle] .getFrameOfReference().getPosition()            =  circle->getFrame()->getPosition() + cpData[icircle].getFrameOfReference().getOrientation().col(0)*rCircle;
            cpData[icircle] .getFrameOfReference().getPosition()            =  circle->getFrame()->getPosition() - rCircle*Wd/nrm2(Wd);

          }
        }
      }
    }

    FuncPairContour1sCircleSolid *func= new FuncPairContour1sCircleSolid(circle,band); // root function for searching contact parameters
    Contact1sSearch search(func);

#if 0
    // deform mesh
    for(int i=0;i<nodes.size();i++) {
      double node_i = nodes(i);
      for(int j=0;j<lastContactPositions.size();j++) {
        nodes(i) += (lastContactPositions(j)-node_i) * exp( - fabs((lastContactPositions(j)-node_i)/ l0));
      }
    }
    if(result.rows()>0) {
      cout << endl <<endl<< "----------------------------------------------\nlast hits: "<< trans(lastContactPositions) << endl;
      Vec tempVec = staticNodes;
      cout << "original Nodes" << endl << trans( tempVec ) << endl;
      cout << "modified Nodes" << endl << trans( nodes   ) << endl;
      cout << trans(result) << endl;
    }
    if ( lastContactPositions.size() ) throw 1;
#endif

    search.setNodes(nodes); // defining search areas for contacts
    Mat result = search.slvAll();
    delete func;

    for(int i=0;i<result.rows();i++) {
      ContourPointData* cpData = icpData[nrNodes + i];

      if(cpData[icontour].getLagrangeParameterPosition().size() != 2)
        cpData[icontour].getLagrangeParameterPosition().resize() = Vec(2,NONINIT);
      cpData[icontour].getLagrangeParameterPosition()(0) = result(i,0);
      cpData[icontour].getLagrangeParameterPosition()(1) = 0.0;

      if(cpData[icontour].getLagrangeParameterPosition()(0) < band->getAlphaStart() || cpData[icontour].getLagrangeParameterPosition()(0) > band->getAlphaEnd())
        (ig[nrNodes + i])(0) = 1.0;
      else {
        bool doubledNode = false;
        for(int j=0;j<staticNodes.size();j++)
          if(fabs(cpData[icontour].getLagrangeParameterPosition()(0)-staticNodes(j))<epsTol)
            doubledNode = true;

        if(doubledNode)
          (ig[nrNodes + i])(0) = 1.0;
        else {
          band->updateKinematicsForFrame(cpData[icontour],position_cosy);
          Vec Wd = circle->getFrame()->getPosition() - cpData[icontour].getFrameOfReference().getPosition();
          (ig[nrNodes + i])(0) = nrm2(Wd) - rCircle;

          if((ig[i])(0) <= 0.0) {
            contactPositions(inContact) = cpData[icontour].getLagrangeParameterPosition()(0);
            inContact ++;
          }

          {
            Vec Wb = cpData[icontour].getFrameOfReference().getOrientation().col(2);
            cpData[icontour].getLagrangeParameterPosition()(1) = trans(Wb)*Wd; // get contact parameter of second tangential direction

            if(fabs(cpData[icontour].getLagrangeParameterPosition()(1)) > 0.5*wBand)
              (ig[i])(0) = 1.0;
            else { // calculate the normal distance
              cpData[icontour].getFrameOfReference().getPosition()           += cpData[icontour].getLagrangeParameterPosition()(1)*Wb;
              cpData[icircle] .getFrameOfReference().getOrientation().set(0, -cpData[icontour].getFrameOfReference().getOrientation().col(0));
              cpData[icircle] .getFrameOfReference().getOrientation().set(1, -cpData[icontour].getFrameOfReference().getOrientation().col(1));
              cpData[icircle] .getFrameOfReference().getOrientation().set(2,  cpData[icontour].getFrameOfReference().getOrientation().col(2));
              cpData[icircle] .getFrameOfReference().getPosition()            =  circle->getFrame()->getPosition() + cpData[icircle].getFrameOfReference().getOrientation().col(0)*rCircle;

#if 0
              if(false && (ig[i])(0) < 0.)//1.0e-4)
              {
                cout << endl;
                cout << " g(0)  = " << (ig[i])(0) << "\n------------------------------------------\n" << endl;
                cout << " " << i << " - alpha = " << cpData[icontour].getLagrangeParameterPosition()(0) << endl;
                cout << "     Wb = " << trans(Wb) << "\n Wd = " << trans(Wd) << endl;
                cout << "     rCircle  = " << rCircle << endl;
                cout << " cpData[icontour].r = \n" << trans(cpData[icontour].getFrameOfReference().getPosition()   ) << endl;
                cout << " cpData[icontour].A = \n" <<       cpData[icontour].getFrameOfReference().getOrientation()  << endl;
                cout << " cpData[icircle].r  = \n" << trans(cpData[icircle] .getFrameOfReference().getPosition()   ) << endl;
                cout << " cpData[icircle].A  = \n" <<       cpData[icircle] .getFrameOfReference().getOrientation()  << endl;
                cpData[icontour].getLagrangeParameterPosition()(0) = 0.0;
                band->updateKinematicsForFrame(cpData[icontour],position_cosy);
                cout << " cpData(alpha=0).r = \n" << trans(cpData[icontour].getFrameOfReference().getPosition()   ) << endl;
                cout << " cpData(alpha=0).A = \n" <<       cpData[icontour].getFrameOfReference().getOrientation()  << endl;
                throw 1;
              }
#endif
            }
          }
        }
      }
    }
    lastContactPositions.resize() = contactPositions(Index(0,inContact-1));

    // flush data not being necessary
    for(int i=nrNodes + result.rows();i<numberOfPotentialContactPoints;i++)
    {
      ContourPointData* cpData = icpData[i];

      // positions -> dummy and centerOfCircle
      if(cpData[icontour].getLagrangeParameterPosition().size() != 2)
        cpData[icontour].getLagrangeParameterPosition().resize() = Vec(2,INIT,0.0);
      cpData[icontour].getLagrangeParameterPosition()(0) = 0.0;

      cpData[icontour].getFrameOfReference().setPosition(nodes(nodes.size()-1)*Vec("[1.0;0.0;0.0]"));//Vec(3,INIT,0.0));//
      cpData[icircle ].getFrameOfReference().setPosition(nodes(nodes.size()-1)*Vec("[1.0;0.0;0.0]"));//Vec(3,INIT,0.0));//
      cpData[icontour].getFrameOfReference().getOrientation() = -SqrMat(3,EYE);
      cpData[icircle ].getFrameOfReference().getOrientation() = -SqrMat(3,EYE);
      (ig[i])(0) = 1.0;
    }

  }

}

