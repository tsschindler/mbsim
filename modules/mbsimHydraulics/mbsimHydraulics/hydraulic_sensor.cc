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
 * Contact: markus.ms.schneider@gmail.com
 */

#include <config.h>
#include "mbsimHydraulics/hydraulic_sensor.h"
#include "mbsimHydraulics/hnode.h"
#include "mbsimHydraulics/hline.h"
#include "mbsimHydraulics/environment.h"

using namespace std;
using namespace fmatvec;
using namespace MBSim;
using namespace MBXMLUtils;
using namespace xercesc;

namespace MBSimHydraulics {

  MBSIM_OBJECTFACTORY_REGISTERXMLNAME(FlowSensor, MBSIMHYDRAULICS%"FlowSensor")

  VecV FlowSensor::getSignal() {
    return line->getQIn(); 
  }

  void FlowSensor::initializeUsingXML(DOMElement * element) {
    Sensor::initializeUsingXML(element);
    DOMElement *e;
    e=E(element)->getFirstElementChildNamed(MBSIMHYDRAULICS%"hline");
    lineString=E(e)->getAttribute("ref");
  }

  void FlowSensor::init(InitStage stage) {
    if (stage==resolveXMLPath) {
      if (lineString!="")
        setHLine(getByPath<HLine>(lineString));
      Sensor::init(stage);
    }
    else
      Sensor::init(stage);
  }

  MBSIM_OBJECTFACTORY_REGISTERXMLNAME(PressureSensor, MBSIMHYDRAULICS%"PressureSensor")

  VecV PressureSensor::getSignal() {
    return node->getla(); 
  }

  void PressureSensor::initializeUsingXML(DOMElement * element) {
    Sensor::initializeUsingXML(element);
    DOMElement *e;
    e=E(element)->getFirstElementChildNamed(MBSIMHYDRAULICS%"hnode");
    nodeString=E(e)->getAttribute("ref");
  }

  void PressureSensor::init(InitStage stage) {
    if (stage==resolveXMLPath) {
      if (nodeString!="")
        setHNode(getByPath<HNode>(nodeString));
      Sensor::init(stage);
    }
    else
      Sensor::init(stage);
  }

  MBSIM_OBJECTFACTORY_REGISTERXMLNAME(TemperatureSensor, MBSIMHYDRAULICS%"TemperatureSensor")
  
  void TemperatureSensor::init(InitStage stage) {
    if (stage==preInit) {
      T(0)=HydraulicEnvironment::getInstance()->getTemperature();
      Sensor::init(stage);
    }
    else
      Sensor::init(stage);
  }

  MBSIM_OBJECTFACTORY_REGISTERXMLNAME(KinematicViscositySensor, MBSIMHYDRAULICS%"KinematicViscositySensor")
  
  void KinematicViscositySensor::init(InitStage stage) {
    if (stage==preInit) {
      nu(0)=HydraulicEnvironment::getInstance()->getKinematicViscosity();
      Sensor::init(stage);
    }
    else
      Sensor::init(stage);
  }

}
