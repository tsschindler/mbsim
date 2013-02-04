/* Copyright (C) 2004-2012 MBSim Development Team
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
 */

#include <config.h>

#include "polynomial_frustum.h"

using namespace std;
using namespace fmatvec;

#ifdef HAVE_OPENMBVCPPINTERFACE
//#include <openmbvcppinterface/frustum.h>
#include <openmbvcppinterface/ivbody.h>
#endif

namespace MBSim {

  PolynomialFrustum::PolynomialFrustum(const std::string & name) :
		        RigidContour(name), parameters(0), height(0.)
#ifdef HAVE_OPENMBVCPPINTERFACE
  , color(LIGHTGRAY), transparency(0.), polynomialPoints(0), circularPoints(25)
#endif
  {

  }

  PolynomialFrustum::~PolynomialFrustum() {
  }

  void PolynomialFrustum::init(InitStage stage) {

    if (stage == MBSim::plot) {
      updatePlotFeatures();

      if (getPlotFeature(plotRecursive) == enabled) {
#ifdef HAVE_OPENMBVCPPINTERFACE
        if(openMBVRigidBody) {
          static_cast<OpenMBV::IvBody*>(openMBVRigidBody)->setIvFileName((this->name + ".iv").c_str());
          static_cast<OpenMBV::IvBody*>(openMBVRigidBody)->setBoundaryEdges(true);
          static_cast<OpenMBV::IvBody*>(openMBVRigidBody)->setInitialTranslation(0.,0.,0.);
          static_cast<OpenMBV::IvBody*>(openMBVRigidBody)->setInitialRotation(0.,0.,0.);

          createInventorFile();
        }
#endif
        RigidContour::init(stage);
      }
    } else
      RigidContour::init(stage);
  }

  void PolynomialFrustum::setPolynomialParameters(
      const vector<double> & parameters_) {
    parameters = parameters_;
  }

  void PolynomialFrustum::setHeight(const double & height_) {
    height = height_;
  }

#ifdef HAVE_OPENMBVCPPINTERFACE
  void PolynomialFrustum::enableOpenMBV(bool enable, int polynomialPoints_, int circularPoints_) {
    if (enable) {
      openMBVRigidBody = new OpenMBV::IvBody;

      if(circularPoints_ <= 0)
        circularPoints = 25;
      else
        circularPoints = circularPoints_;

      if(polynomialPoints_ <= 0)
        polynomialPoints = 4 * parameters.size();
      else
        polynomialPoints = polynomialPoints_;

    } else {
      //TODO: destructor of OpenMBV::RigidBody is private (not public...)
      //      if(openMBVRigidBody)
      //        delete openMBVRigidBody;
      openMBVRigidBody = 0;
    }
  }

  void PolynomialFrustum::setColor(const RGBColor & color_) {
    color = color_;
  }

  void PolynomialFrustum::setTransparency(const double & transparency_) {
    transparency = transparency_;
  }
#endif

  double PolynomialFrustum::getValue(const double & x) {
    double val = 0;
    for (size_t i = 0; i < parameters.size(); i++) {
      val += parameters[i] * pow(x, i);
    }
    return val;
  }

  double PolynomialFrustum::getValueD1(const double & x) {
    double val = 0;
    for (size_t i = 1; i < parameters.size(); i++) {
      val += i * parameters[i] * pow(x, (i - 1));
    }
    return val;
  }

  double PolynomialFrustum::getValueD2(const double & x) {
    double val = 0;
    for (size_t i = 2; i < parameters.size(); i++) {
      val += i * (i - 1) * parameters[i] * pow(x, (i - 2));
    }
    return val;
  }

  double PolynomialFrustum::getXPolyMax() {
    double x = height / 2;
    int iter = 1;
    while (iter < 1000 && fabs(getValueD1(x)) > 1e-6) {
      x = x - getValueD1(x) / getValueD2(x);
      iter++;
    }
    return x;
  }

  double PolynomialFrustum::getRadiSphere() {
    double rad = 0.;
    double fa = getValue(0);
    double fb = getValue(height);
    double fda = getValueD1(0); //f'(a)
    double fdb = getValueD1(height); //f'(b)
    double temp1 = sqrt(pow(fb , 2) + pow(height , 2) / 4);
    double temp2 = sqrt(pow(fa , 2) + pow(height , 2) / 4);
    temp2 = (temp1 > temp2) ? temp1 : temp2;
    if (fda * fdb >= 0) {
      temp1 = 0;
    } else {
      double x = getXPolyMax();
      double fmax = getValue(x);
      temp1 = sqrt(pow(fmax , 2) + pow((x - height / 2) , 2));
    }

    rad = (temp1 > temp2) ? temp1 : temp2;

    return rad;
  }

  void PolynomialFrustum::createInventorFile() {

    //TODO: Use IndexedTriangleSet instead of IndexedFaceSet (should be faster)

    std::ofstream ivFile;

    ivFile.open((this->name + ".iv").c_str());

    /*HEAD*/
    ivFile << "#Inventor V2.1 ascii" << endl << endl;

    ivFile << "Separator" << endl << "{" << endl;

    /*MATERIAL*/
    ivFile << "Material" << endl << "{" << endl;

    ivFile << "diffuseColor " <<  color.red << " " << color.green << " " << color.blue << endl;

    ivFile << "  transparency " << transparency << endl;
    ivFile << "}" << endl;

    /* vertices BEGIN (=Coordinates of the points)*/
    ivFile << "  Coordinate3" << endl << "  {" << endl;
    ivFile << "    point [" << endl;

    for(int i = 0; i < polynomialPoints + 1; i++) {
      double x = (height * i) / polynomialPoints;
      double r = getValue(x);
      for(int j = 0; j < circularPoints; j++) {
        double phi = j * 2. * M_PI / circularPoints;
        double y = r * cos(phi);
        double z = r * sin(phi);
        ivFile << "      " << x << " " << y << " " << z << "," << endl;
      }
    }

    ivFile << "    ]" << endl;

    ivFile << "  }" << endl << endl;

    /*Vertices END*/


    ivFile << "  ShapeHints" << endl << "  {" << endl; // hints BEGIN
    ivFile << "    vertexOrdering COUNTERCLOCKWISE" << endl; //  CLOCKWISE means look inside
    ivFile << "    shapeType UNKNOWN_SHAPE_TYPE" << endl;
    ivFile << "    creaseAngle 0.393" << endl;
    ivFile << "  }" << endl << endl; // hints END

    /*Faces BEGIN (describe the surfaces)*/

    /*Circle Bottom*/
    ivFile << "  IndexedFaceSet" << endl << "  {" << endl;
    ivFile << "    coordIndex [" << endl;

    int firstInd = 0;
    int secondInd = 1;
    int thirdInd = circularPoints - 1;

    for(int i = 0; i < circularPoints - 2; i++) {
      ivFile << firstInd << ", " << secondInd << ", " << thirdInd << "-1" << endl;
      int newThirdInd = secondInd + 1;
      if(thirdInd < secondInd)
        newThirdInd = secondInd - 1;
      firstInd = secondInd;
      secondInd = thirdInd;
      thirdInd = newThirdInd;
    }

    ivFile << "    ]" << endl;
    ivFile << "  }" << endl << endl;

    /*Circle Top*/

    ivFile << "  IndexedFaceSet" << endl << "  {" << endl;
    ivFile << "    coordIndex [" << endl;

    firstInd = circularPoints * polynomialPoints;
    secondInd = circularPoints * polynomialPoints + 1;
    thirdInd = circularPoints * (polynomialPoints +1) - 1;

    for(int i = 0; i < circularPoints - 2; i++) {
      ivFile << firstInd << ", " << secondInd << ", " << thirdInd << "-1" << endl;
      int newThirdInd = secondInd + 1;
      if(thirdInd < secondInd)
        newThirdInd = secondInd - 1;
      firstInd = secondInd;
      secondInd = thirdInd;
      thirdInd = newThirdInd;
    }

    ivFile << "    ]" << endl;
    ivFile << "  }" << endl << endl;

    /*In Between*/
    ivFile << "  IndexedFaceSet" << endl << "  {" << endl;
    ivFile << "    coordIndex [" << endl;

    for(int i = 0; i < circularPoints; i++) {
      int ind = i + 1;

      //ring closure
      if(i >= circularPoints -1) {
        ind = 0;
      }

      //move up current polynomial line
      for(int j = 0; j < polynomialPoints ; j++) {
        ivFile << "      " << i + j * circularPoints << ", " << i + (j+1) * circularPoints << ", " << ind + j * circularPoints << "-1" << endl;
        ivFile << "      " << i + (j+1) * circularPoints << ", " << ind + (j+1) * circularPoints << ", " << ind + j * circularPoints << "-1" << endl;
      }

    }

    ivFile << "    ]" << endl;
    ivFile << "  }" << endl << endl;

    /*Faces END*/

    /*BOTTOM*/
    ivFile << "}" << endl << endl;

    ivFile.close();
  }

}