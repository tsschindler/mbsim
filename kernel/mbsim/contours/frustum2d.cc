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
 */

#include<config.h>
#include "mbsim/contours/frustum2d.h"

#ifdef HAVE_OPENMBVCPPINTERFACE
#include <openmbvcppinterface/frustum.h>
#endif


using namespace std;
using namespace fmatvec;
using namespace boost;

namespace MBSim {

  void Frustum2D::init(InitStage stage) {
    if(stage==plotting) {
      updatePlotFeatures();
  
      if(getPlotFeature(plotRecursive)==enabled) {
  #ifdef HAVE_OPENMBVCPPINTERFACE
        if(getPlotFeature(openMBV)==enabled && openMBVRigidBody) {
          static_pointer_cast<OpenMBV::Frustum>(openMBVRigidBody)->setInitialTranslation(0.,h,0.);
          static_pointer_cast<OpenMBV::Frustum>(openMBVRigidBody)->setInitialRotation(3./2.*M_PI,0,0.);
          static_pointer_cast<OpenMBV::Frustum>(openMBVRigidBody)->setBaseRadius(r(0));
          static_pointer_cast<OpenMBV::Frustum>(openMBVRigidBody)->setTopRadius(r(1));
          static_pointer_cast<OpenMBV::Frustum>(openMBVRigidBody)->setHeight(h);
        }
  #endif
        RigidContour::init(stage);
      }
    }
    else
      RigidContour::init(stage);
  }

}
