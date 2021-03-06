/* Copyright (C) 2004-2009 MBSim Development Team
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
 *
 * Contact: martin.o.foerg@googlemail.com
 *          rzander@users.berlios.de
 */

#include <config.h> 
#include "circlesolid_line.h"
#include "mbsim/contours/line.h"
#include "mbsim/contours/circle_solid.h"


using namespace fmatvec;
using namespace std;

namespace MBSim {

  void ContactKinematicsCircleSolidLine::assignContours(const vector<Contour*> &contour) {
    if(dynamic_cast<CircleSolid*>(contour[0])) {
      icircle = 0; iline = 1;
      circlesolid = static_cast<CircleSolid*>(contour[0]);
      line = static_cast<Line*>(contour[1]);
    } 
    else {
      icircle = 1; iline = 0;
      circlesolid = static_cast<CircleSolid*>(contour[1]);
      line = static_cast<Line*>(contour[0]);
    }
  }

  void ContactKinematicsCircleSolidLine::updateg(Vec &g, ContourPointData *cpData, int index) {

    cpData[iline].getFrameOfReference().setOrientation(line->getFrame()->getOrientation());
    cpData[icircle].getFrameOfReference().getOrientation().set(0, -line->getFrame()->getOrientation().col(0));
    cpData[icircle].getFrameOfReference().getOrientation().set(1, -line->getFrame()->getOrientation().col(1));
    cpData[icircle].getFrameOfReference().getOrientation().set(2, line->getFrame()->getOrientation().col(2));

    Vec3 Wn = cpData[iline].getFrameOfReference().getOrientation().col(0);

    Vec3 Wd = circlesolid->getFrame()->getPosition() - line->getFrame()->getPosition();

    g(0) = Wn.T()*Wd - circlesolid->getRadius();

    cpData[icircle].getFrameOfReference().setPosition(circlesolid->getFrame()->getPosition() - Wn*circlesolid->getRadius());
    cpData[iline].getFrameOfReference().setPosition(cpData[icircle].getFrameOfReference().getPosition() - Wn*g(0));
  }

  void ContactKinematicsCircleSolidLine::updatewb(Vec &wb, const Vec &g, ContourPointData *cpData) {

    Vec3 v2 = cpData[icircle].getFrameOfReference().getOrientation().col(2);
    Vec3 n1 = cpData[iline].getFrameOfReference().getOrientation().col(0);
    // Vec3 n2 = cpData[icircle].getFrameOfReference().getOrientation().col(0);
    Vec3 u1 = cpData[iline].getFrameOfReference().getOrientation().col(1);
    Vec3 u2 = cpData[icircle].getFrameOfReference().getOrientation().col(1);
    Vec3 vC1 = cpData[iline].getFrameOfReference().getVelocity();
    Vec3 vC2 = cpData[icircle].getFrameOfReference().getVelocity();
    Vec3 Om1 = cpData[iline].getFrameOfReference().getAngularVelocity();
    Vec3 Om2 = cpData[icircle].getFrameOfReference().getAngularVelocity();
    double r = circlesolid->getRadius();

    double ad2 = -v2.T()*(Om2-Om1);
    double ad1 = u1.T()*(vC2-vC1) - r*ad2;
    Vec3 s2 = u2*r;

    wb(0) += n1.T()*(-crossProduct(Om1,vC2-vC1) - crossProduct(Om1,u1)*ad1 + crossProduct(Om2,s2)*ad2);
    
    if(wb.size() > 1) 
      wb(1) += u1.T()*(-crossProduct(Om1,vC2-vC1) - crossProduct(Om1,u1)*ad1 + crossProduct(Om2,s2)*ad2);
  }
      
  void ContactKinematicsCircleSolidLine::computeCurvatures(Vec &r, ContourPointData* cpData) {
    r(icircle)=circlesolid->computeCurvature(cpData[icircle]);
    r(iline)=line->computeCurvature(cpData[iline]);
  }

}

