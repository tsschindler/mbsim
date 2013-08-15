/*
    MBSimGUI - A fronted for MBSim.
    Copyright (C) 2012 Martin Förg

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _OBSERVER__H_
#define _OBSERVER__H_

#include "element.h"
#include "extended_properties.h"

class Observer : public Element {
  public:
    Observer(const std::string &str, Element *parent);
    ~Observer();
    static Observer* readXMLFile(const std::string &filename, Element *parent);
    virtual int getxSize() {return 0;}
    virtual Element* getByPathSearch(std::string path);
};

class AbsoluteKinematicsObserver : public Observer {
  friend class AbsoluteKinematicsObserverPropertyDialog;
  public:
    AbsoluteKinematicsObserver(const std::string &str, Element *parent);
    ~AbsoluteKinematicsObserver();
    std::string getType() const { return "AbsoluteKinematicsObserver"; }
    virtual void initializeUsingXML(MBXMLUtils::TiXmlElement *element);
    virtual MBXMLUtils::TiXmlElement* writeXMLFile(MBXMLUtils::TiXmlNode *element);
    void initialize();
    ElementPropertyDialog* createPropertyDialog() {return new AbsoluteKinematicsObserverPropertyDialog(this);}
  protected:
    ExtProperty frame, position, velocity, angularVelocity, acceleration, angularAcceleration;
};

#endif