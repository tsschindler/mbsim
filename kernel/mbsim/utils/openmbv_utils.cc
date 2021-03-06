/* Copyright (C) 2004-2013 MBSim Development Team
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
 * Contact: martin.o.foerg@gmail.com
 */

#include <config.h>
#include "mbsim/element.h"
#include "mbsim/utils/openmbv_utils.h"
#include <xercesc/dom/DOMProcessingInstruction.hpp>

using namespace std;
using namespace MBXMLUtils;
using namespace fmatvec;
using namespace xercesc;
using namespace boost;

#ifdef HAVE_OPENMBVCPPINTERFACE

namespace MBSim {

  void OpenMBVObject::initializeUsingXML(DOMElement *e) {
    DOMElement *ee;
    ee=E(e)->getFirstElementChildNamed(MBSIM%"diffuseColor");
    if(ee) dc = Element::getVec(ee, 3);
    ee=E(e)->getFirstElementChildNamed(MBSIM%"transparency");
    if(ee) tp = Element::getDouble(ee);

    DOMProcessingInstruction *ID = E(e)->getFirstProcessingInstructionChildNamed("OPENMBV_ID");
    if(ID)
      id = X()%ID->getData();
  }

  void OpenMBVObject::initializeObject(const shared_ptr<OpenMBV::DynamicColoredBody> &object) {
    object->setDiffuseColor(dc(0),dc(1),dc(2));
    object->setTransparency(tp);
    object->setID(id);
  }

  void OpenMBVArrow::initializeUsingXML(DOMElement *e) {
    OpenMBVObject::initializeUsingXML(e);
    DOMElement *ee = E(e)->getFirstElementChildNamed(MBSIM%"scaleLength");
    if(ee) sL = Element::getDouble(ee);
    ee = E(e)->getFirstElementChildNamed(MBSIM%"scaleSize");
    if(ee) sS = Element::getDouble(ee);
    ee = E(e)->getFirstElementChildNamed(MBSIM%"referencePoint");
    if(ee) {
      string rP=string(X()%E(ee)->getFirstTextChild()->getData()).substr(1,string(X()%E(ee)->getFirstTextChild()->getData()).length()-2);
      if(rP=="toPoint")   refPoint=OpenMBV::Arrow::toPoint;
      if(rP=="fromPoint") refPoint=OpenMBV::Arrow::fromPoint;
      if(rP=="midPoint")  refPoint=OpenMBV::Arrow::midPoint;
    }
  }

  shared_ptr<OpenMBV::Arrow> OpenMBVArrow::createOpenMBV(DOMElement *e) {
    shared_ptr<OpenMBV::Arrow> object = OpenMBV::ObjectFactory::create<OpenMBV::Arrow>();
    if(e) initializeUsingXML(e);
    initializeObject(object);
    return object;
  }

  void OpenMBVArrow::initializeObject(const shared_ptr<OpenMBV::Arrow> &object) {
    OpenMBVObject::initializeObject(object);
    object->setDiameter(0.25*sS);
    object->setHeadDiameter(0.5*sS);
    object->setHeadLength(0.75*sS);
    object->setType(type);
    object->setReferencePoint(refPoint);
    object->setScaleLength(sL);
  }

  void OpenMBVFrame::initializeUsingXML(DOMElement *e) {
    OpenMBVObject::initializeUsingXML(e);
    DOMElement *ee = E(e)->getFirstElementChildNamed(MBSIM%"size");
    if(ee) size = Element::getDouble(ee);
    ee = E(e)->getFirstElementChildNamed(MBSIM%"offset");
    if(ee) offset = Element::getDouble(ee);
  }

  shared_ptr<OpenMBV::Frame> OpenMBVFrame::createOpenMBV(DOMElement *e) {
    shared_ptr<OpenMBV::Frame> object = OpenMBV::ObjectFactory::create<OpenMBV::Frame>();
    if(e) initializeUsingXML(e);
    initializeObject(object);
    return object;
  }

  void OpenMBVFrame::initializeObject(const shared_ptr<OpenMBV::Frame> &object) {
    OpenMBVObject::initializeObject(object);
    object->setSize(size);
    object->setOffset(offset);
  }

  void OpenMBVSphere::initializeUsingXML(DOMElement *e) {
    OpenMBVObject::initializeUsingXML(e);
    DOMElement *ee;
    ee = E(e)->getFirstElementChildNamed(MBSIM%xml);
    if(ee) r = Element::getDouble(ee);
  }

  shared_ptr<OpenMBV::Sphere> OpenMBVSphere::createOpenMBV(DOMElement *e) {
    shared_ptr<OpenMBV::Sphere> object = OpenMBV::ObjectFactory::create<OpenMBV::Sphere>();
    if(e) initializeUsingXML(e);
    initializeObject(object);
    return object;
  }

  void OpenMBVSphere::initializeObject(const shared_ptr<OpenMBV::Sphere> &object) {
    OpenMBVObject::initializeObject(object);
    object->setRadius(r);
  }

  void OpenMBVLine::initializeUsingXML(DOMElement *e) {
    OpenMBVObject::initializeUsingXML(e);
    DOMElement *ee;
    ee = E(e)->getFirstElementChildNamed(MBSIM%"length");
    if(ee) l = Element::getDouble(ee);
  }

  shared_ptr<OpenMBV::Cuboid> OpenMBVLine::createOpenMBV(DOMElement *e) {
    shared_ptr<OpenMBV::Cuboid> object = OpenMBV::ObjectFactory::create<OpenMBV::Cuboid>();
    if(e) initializeUsingXML(e);
    initializeObject(object);
    return object;
  }

  void OpenMBVLine::initializeObject(const shared_ptr<OpenMBV::Cuboid> &object) {
    OpenMBVObject::initializeObject(object);
    object->setLength(0,l,0);
  }

  void OpenMBVPlane::initializeUsingXML(DOMElement *e) {
    OpenMBVObject::initializeUsingXML(e);
    DOMElement *ee;
    ee = E(e)->getFirstElementChildNamed(MBSIM%"length");
    if(ee) l = Element::getVec(ee,2);
  }

  shared_ptr<OpenMBV::Cuboid> OpenMBVPlane::createOpenMBV(DOMElement *e) {
    shared_ptr<OpenMBV::Cuboid> object = OpenMBV::ObjectFactory::create<OpenMBV::Cuboid>();
    if(e) initializeUsingXML(e);
    initializeObject(object);
    return object;
  }

  void OpenMBVPlane::initializeObject(const shared_ptr<OpenMBV::Cuboid> &object) {
    OpenMBVObject::initializeObject(object);
    object->setLength(0,l(0),l(1));
  }

  void OpenMBVCuboid::initializeUsingXML(DOMElement *e) {
    OpenMBVObject::initializeUsingXML(e);
    DOMElement *ee;
    ee = E(e)->getFirstElementChildNamed(MBSIM%"length");
    if(ee) l = Element::getVec(ee,3);
  }

  shared_ptr<OpenMBV::Cuboid> OpenMBVCuboid::createOpenMBV(DOMElement *e) {
    shared_ptr<OpenMBV::Cuboid> object = OpenMBV::ObjectFactory::create<OpenMBV::Cuboid>();
    if(e) initializeUsingXML(e);
    initializeObject(object);
    return object;
  }

  void OpenMBVCuboid::initializeObject(const shared_ptr<OpenMBV::Cuboid> &object) {
    OpenMBVObject::initializeObject(object);
    object->setLength(l(0),l(1),l(2));
  }

  void OpenMBVCircle::initializeUsingXML(DOMElement *e) {
    OpenMBVObject::initializeUsingXML(e);
    DOMElement *ee;
    ee = E(e)->getFirstElementChildNamed(MBSIM%"radius");
    if(ee) r = Element::getDouble(ee);
  }

  shared_ptr<OpenMBV::Frustum> OpenMBVCircle::createOpenMBV(DOMElement *e) {
    shared_ptr<OpenMBV::Frustum> object = OpenMBV::ObjectFactory::create<OpenMBV::Frustum>();
    if(e) initializeUsingXML(e);
    initializeObject(object);
    return object;
  }

  void OpenMBVCircle::initializeObject(const shared_ptr<OpenMBV::Frustum> &object) {
    OpenMBVObject::initializeObject(object);
    object->setTopRadius(r);
    object->setBaseRadius(r);
    object->setHeight(0);
  }

  void OpenMBVFrustum::initializeUsingXML(DOMElement *e) {
    OpenMBVObject::initializeUsingXML(e);
    DOMElement *ee;
    ee = E(e)->getFirstElementChildNamed(MBSIM%"topRadius");
    if(ee) t = Element::getDouble(ee);
    ee = E(e)->getFirstElementChildNamed(MBSIM%"baseRadius");
    if(ee) b = Element::getDouble(ee);
    ee = E(e)->getFirstElementChildNamed(MBSIM%"height");
    if(ee) h = Element::getDouble(ee);
  }

  shared_ptr<OpenMBV::Frustum> OpenMBVFrustum::createOpenMBV(DOMElement *e) {
    shared_ptr<OpenMBV::Frustum> object = OpenMBV::ObjectFactory::create<OpenMBV::Frustum>();
    if(e) initializeUsingXML(e);
    initializeObject(object);
    return object;
  }

  void OpenMBVFrustum::initializeObject(const shared_ptr<OpenMBV::Frustum> &object) {
    OpenMBVObject::initializeObject(object);
    object->setTopRadius(t);
    object->setBaseRadius(b);
    object->setHeight(h);
  }

  void OpenMBVExtrusion::initializeUsingXML(DOMElement *e) {
    OpenMBVObject::initializeUsingXML(e);
    DOMElement *ee = E(e)->getFirstElementChildNamed(MBSIM%"height");
    if(ee) h = Element::getDouble(ee);
  }

  shared_ptr<OpenMBV::Extrusion> OpenMBVExtrusion::createOpenMBV(DOMElement *e) {
    shared_ptr<OpenMBV::Extrusion> object = OpenMBV::ObjectFactory::create<OpenMBV::Extrusion>();
    if(e) initializeUsingXML(e);
    initializeObject(object);
    return object;
  }

  void OpenMBVExtrusion::initializeObject(const shared_ptr<OpenMBV::Extrusion> &object) {
    OpenMBVObject::initializeObject(object);
    object->setHeight(h);
  }

  void OpenMBVCoilSpring::initializeUsingXML(DOMElement *e) {
    OpenMBVObject::initializeUsingXML(e);
    DOMElement *ee = E(e)->getFirstElementChildNamed(MBSIM%"numberOfCoils");
    if(ee) n = Element::getInt(ee);
    ee = E(e)->getFirstElementChildNamed(MBSIM%"springRadius");
    if(ee) r = Element::getDouble(ee);
    ee = E(e)->getFirstElementChildNamed(MBSIM%"crossSectionRadius");
    if(ee) cr = Element::getDouble(ee);
    ee = E(e)->getFirstElementChildNamed(MBSIM%"nominalLength");
    if(ee) l = Element::getDouble(ee);
    ee = E(e)->getFirstElementChildNamed(MBSIM%"scaleFactor");
    if(ee) sf = Element::getDouble(ee);
    ee = E(e)->getFirstElementChildNamed(MBSIM%"type");
    if(ee) {
      string typeStr=string(X()%E(ee)->getFirstTextChild()->getData()).substr(1,string(X()%E(ee)->getFirstTextChild()->getData()).length()-2);
      if(typeStr=="tube") type=OpenMBV::CoilSpring::tube;
      if(typeStr=="scaledTube") type=OpenMBV::CoilSpring::scaledTube;
      if(typeStr=="polyline") type=OpenMBV::CoilSpring::polyline;
    }
  }

  shared_ptr<OpenMBV::CoilSpring> OpenMBVCoilSpring::createOpenMBV(DOMElement *e) {
    shared_ptr<OpenMBV::CoilSpring> object = OpenMBV::ObjectFactory::create<OpenMBV::CoilSpring>();
    if(e) initializeUsingXML(e);
    initializeObject(object);
    return object;
  }

  void OpenMBVCoilSpring::initializeObject(const shared_ptr<OpenMBV::CoilSpring> &object) {
    OpenMBVObject::initializeObject(object);
    object->setSpringRadius(r);
    object->setCrossSectionRadius(cr);
    object->setScaleFactor(sf);
    object->setNumberOfCoils(n);
    object->setNominalLength(l);
    object->setType(type);
  }

}

#endif
