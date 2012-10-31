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

#include "kinetic_excitation.h"
#include "utils.h"
#include <QtGui/QMenu>
#include "frame.h"

using namespace std;


KineticExcitation::KineticExcitation(const QString &str, QTreeWidgetItem *parentItem, int ind) : Link(str, parentItem, ind) {

  setText(1,getType());

  properties->addTab("Kinetics");
  //properties->addTab("Constitutive laws");

  connections = new ConnectEditor(1, this, properties, Utils::QIconCached("lines.svg"), "Connection", "Kinetics");
  connect(connections->getLineEdit(0),SIGNAL(textChanged(const QString&)),this,SLOT(updateFrameOfReference()));
  force = new ForceLawEditor(properties, Utils::QIconCached("lines.svg"), true);
  moment = new ForceLawEditor(properties, Utils::QIconCached("lines.svg"), false);
  frameOfReference = new FrameOfReferenceEditor(this,properties, Utils::QIconCached("lines.svg"), "Frame of reference", "Kinetics", 0);

  properties->addStretch();
}

KineticExcitation::~KineticExcitation() {
}

void KineticExcitation::updateFrameOfReference() {
  //if(frameOfReference->getFrame()==0)
    frameOfReference->setFrame(connections->getFrame(0));
    frameOfReference->update();
}

void KineticExcitation::initializeUsingXML(TiXmlElement *element) {
  Link::initializeUsingXML(element);
  TiXmlElement *e=element->FirstChildElement(MBSIMNS"frameOfReference");
  if(e)
    saved_frameOfReference=e->Attribute("ref");
  force->initializeUsingXML(element);
  moment->initializeUsingXML(element);
  e=element->FirstChildElement(MBSIMNS"connect");
  saved_ref=e->Attribute("ref");
}

TiXmlElement* KineticExcitation::writeXMLFile(TiXmlNode *parent) {
  TiXmlElement *ele0 = Link::writeXMLFile(parent);
  TiXmlElement *ele1 = new TiXmlElement( MBSIMNS"frameOfReference" );
  if(frameOfReference->getFrame())
    ele1->SetAttribute("ref", frameOfReference->getFrame()->getXMLPath(this,true).toStdString()); // relative path
  ele0->LinkEndChild(ele1);
  force->writeXMLFile(ele0);
  moment->writeXMLFile(ele0);
  ele1 = new TiXmlElement(MBSIMNS"connect");
  if(connections->getFrame(0))
    ele1->SetAttribute("ref", connections->getFrame(0)->getXMLPath(this,true).toStdString()); // relative path
  ele0->LinkEndChild(ele1);
  return ele0;
}

void KineticExcitation::initialize() {
  connections->getLineEdit(0)->blockSignals(true);
  if(saved_ref!="")
    connections->setFrame(0,getByPath<Frame>(saved_ref));
  connections->getLineEdit(0)->blockSignals(false);
  if(saved_frameOfReference!="")
    frameOfReference->setFrame(getByPath<Frame>(saved_frameOfReference));
  else if(connections->getFrame(0))
    frameOfReference->setFrame(connections->getFrame(0));
}