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

#include <config.h>
#include "kinetic_excitation.h"
#include "kinetics_widgets.h"
#include "extended_widgets.h"
#include "ombv_widgets.h"

using namespace std;

KineticExcitation::KineticExcitation(const QString &str, QTreeWidgetItem *parentItem, int ind) : Link(str, parentItem, ind) {

  setText(1,getType());

  properties->addTab("Kinetics");
  //properties->addTab("Constitutive laws");
  properties->addTab("Visualisation");

  forceArrow = new ExtXMLWidget("OpenMBV force arrow",new OMBVArrowWidget("NOTSET"),true);
  ((OMBVArrowWidget*)forceArrow->getWidget())->setID(getID());
  properties->addToTab("Visualisation",forceArrow);

  momentArrow = new ExtXMLWidget("OpenMBV moment arrow",new OMBVArrowWidget("NOTSET"),true);
  ((OMBVArrowWidget*)momentArrow->getWidget())->setID(getID());
  properties->addToTab("Visualisation",momentArrow);

  vector<QWidget*> widget;
  vector<string> name;
  name.push_back("1 frame");
  name.push_back("2 frames");
  widget.push_back(new ConnectFramesWidget(1,this));
  widget.push_back(new ConnectFramesWidget(2,this));

  connections = new ExtXMLWidget("Connections",new XMLWidgetChoiceWidget(name,widget)); 
//  connections = new ExtXMLWidget("Connections",new ConnectFramesWidget(1,this));
  properties->addToTab("Kinetics",connections);

  ForceChoiceWidget *f = new ForceChoiceWidget(MBSIMNS"force", forceArrow);
  force = new ExtXMLWidget("Force",f,true);
  properties->addToTab("Kinetics",force);

  ForceChoiceWidget *m = new ForceChoiceWidget(MBSIMNS"moment", momentArrow);
  moment = new ExtXMLWidget("Moment",m,true);
  properties->addToTab("Kinetics",moment);

  FrameOfReferenceWidget* ref = new FrameOfReferenceWidget(MBSIMNS"frameOfReference",this,0);
  frameOfReference = new ExtXMLWidget("Frame of reference",ref,true);
  properties->addToTab("Kinetics",frameOfReference);

  properties->addStretch();
}

KineticExcitation::~KineticExcitation() {
}

void KineticExcitation::initializeUsingXML(TiXmlElement *element) {
  Link::initializeUsingXML(element);
  frameOfReference->initializeUsingXML(element);
  force->initializeUsingXML(element);
  moment->initializeUsingXML(element);
  connections->initializeUsingXML(element);
}

TiXmlElement* KineticExcitation::writeXMLFile(TiXmlNode *parent) {
  TiXmlElement *ele0 = Link::writeXMLFile(parent);
  frameOfReference->writeXMLFile(ele0);
  force->writeXMLFile(ele0);
  moment->writeXMLFile(ele0);
  connections->writeXMLFile(ele0);
  return ele0;
}
