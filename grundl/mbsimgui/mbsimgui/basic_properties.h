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

#ifndef _BASIC_PROPERTIES_H_
#define _BASIC_PROPERTIES_H_

#include <string>
#include "utils.h"
#include "extended_properties.h"

class Element;
class Frame;
class Contour;
class RigidBody;
namespace MBXMLUtils {
  class TiXmlElement;
  class TiXmlNode;
}

class LocalFrameOfReferenceProperty : public Property {
  protected:
    std::string frame;
    Element* element;
    std::string xmlName;
  public:
    LocalFrameOfReferenceProperty(const std::string &frame_="", Element* element_=0, const std::string &xmlName_="") : frame(frame_), element(element_), xmlName(xmlName_) {}
    MBXMLUtils::TiXmlElement* initializeUsingXML(MBXMLUtils::TiXmlElement *element);
    MBXMLUtils::TiXmlElement* writeXMLFile(MBXMLUtils::TiXmlNode *element); 
    void fromWidget(QWidget *widget);
    void toWidget(QWidget *widget);
    void setFrame(const std::string &str) {frame = str;}
    const std::string& getFrame() const {return frame;}
};

class ParentFrameOfReferenceProperty : public Property {
  protected:
    std::string frame;
    Element* element;
    std::string xmlName;
  public:
    ParentFrameOfReferenceProperty(const std::string &frame_="", Element* element_=0, const std::string &xmlName_="") : frame(frame_), element(element_), xmlName(xmlName_) {}
    virtual MBXMLUtils::TiXmlElement* initializeUsingXML(MBXMLUtils::TiXmlElement *element);
    virtual MBXMLUtils::TiXmlElement* writeXMLFile(MBXMLUtils::TiXmlNode *element);
    void fromWidget(QWidget *widget);
    void toWidget(QWidget *widget);
    void setFrame(const std::string &str) {frame = str;}
    const std::string& getFrame() const {return frame;}
};

class FrameOfReferenceProperty : public Property {
  protected:
    std::string frame;
    Element* element;
    std::string xmlName;
  public:
    FrameOfReferenceProperty(const std::string &frame_="", Element* element_=0, const std::string &xmlName_="") : frame(frame_), element(element_), xmlName(xmlName_) {}
    MBXMLUtils::TiXmlElement* initializeUsingXML(MBXMLUtils::TiXmlElement *element);
    MBXMLUtils::TiXmlElement* writeXMLFile(MBXMLUtils::TiXmlNode *element); 
    void fromWidget(QWidget *widget);
    void toWidget(QWidget *widget);
    void setFrame(const std::string &str) {frame = str;}
    const std::string& getFrame() const {return frame;}
};

class ContourOfReferenceProperty : public Property {
  protected:
    std::string contour;
    Element* element;
    std::string xmlName;
  public:
    ContourOfReferenceProperty(const std::string &contour_="", Element* element_=0, const std::string &xmlName_="") : contour(contour_), element(element_), xmlName(xmlName_) {}
    MBXMLUtils::TiXmlElement* initializeUsingXML(MBXMLUtils::TiXmlElement *element);
    MBXMLUtils::TiXmlElement* writeXMLFile(MBXMLUtils::TiXmlNode *element); 
    void fromWidget(QWidget *widget);
    void toWidget(QWidget *widget);
    void setContour(const std::string &str) {contour = str;}
    const std::string& getContour() const {return contour;}
};

class RigidBodyOfReferenceProperty : public Property {
  protected:
    std::string body;
    Element* element;
    std::string xmlName;
  public:
    RigidBodyOfReferenceProperty(const std::string &body_="", Element *element_=0, const std::string &xmlName_="") : body(body_), element(element_), xmlName(xmlName_) {}
    MBXMLUtils::TiXmlElement* initializeUsingXML(MBXMLUtils::TiXmlElement *element);
    MBXMLUtils::TiXmlElement* writeXMLFile(MBXMLUtils::TiXmlNode *element); 
    void fromWidget(QWidget *widget);
    void toWidget(QWidget *widget);
    void setBody(const std::string &str) {body = str;}
    const std::string& getBody() const {return body;}
};

class FileProperty : public Property {

  public:
    FileProperty(const std::string &xmlName_) : xmlName(xmlName_) {}
    virtual MBXMLUtils::TiXmlElement* initializeUsingXML(MBXMLUtils::TiXmlElement *element);
    virtual MBXMLUtils::TiXmlElement* writeXMLFile(MBXMLUtils::TiXmlNode *element);
    void fromWidget(QWidget *widget);
    void toWidget(QWidget *widget);
    const std::string& getFileName() const {return fileName;}
    void setFileName(const std::string &str) {fileName=str;}
    std::string getAbsoluteFilePath() const;
    void setAbsoluteFilePath(const std::string &str);

  protected:
    std::string fileName;
    std::string xmlName;
    std::string absoluteFilePath;
};

class TextProperty : public Property {

  public:
    TextProperty(const std::string &text_, const std::string &xmlName_, int quote_=0) : text(text_), xmlName(xmlName_), quote(quote_) {}
    virtual MBXMLUtils::TiXmlElement* initializeUsingXML(MBXMLUtils::TiXmlElement *element);
    virtual MBXMLUtils::TiXmlElement* writeXMLFile(MBXMLUtils::TiXmlNode *element);
    void fromWidget(QWidget *widget);
    void toWidget(QWidget *widget);
    const std::string& getText() const {return text;}
    void setText(const std::string &text_) {text = text_;}

  protected:
    std::string text;
    std::string xmlName;
    int quote;
};

class DependenciesProperty : public Property {

  public:
    DependenciesProperty(Element* element_, const std::string &xmlName_) : element(element_), xmlName(xmlName_) {}

    void initialize();
    virtual MBXMLUtils::TiXmlElement* initializeUsingXML(MBXMLUtils::TiXmlElement *element);
    virtual MBXMLUtils::TiXmlElement* writeXMLFile(MBXMLUtils::TiXmlNode *element);
    void fromWidget(QWidget *widget);
    void toWidget(QWidget *widget);

  protected:
    Element* element;
    std::string xmlName;
    std::vector<RigidBodyOfReferenceProperty*> refBody;

    void addDependency();
    void updateGeneralizedCoordinatesOfBodies();
};

class ConnectFramesProperty : public Property {

  public:
    ConnectFramesProperty(int n, Element* element);

    void initialize();
    virtual MBXMLUtils::TiXmlElement* initializeUsingXML(MBXMLUtils::TiXmlElement *element);
    virtual MBXMLUtils::TiXmlElement* writeXMLFile(MBXMLUtils::TiXmlNode *element);
    void fromWidget(QWidget *widget);
    void toWidget(QWidget *widget);

  protected:
    std::vector<FrameOfReferenceProperty*> frame;
    Element* element;
};

class ConnectContoursProperty : public Property {

  public:
    ConnectContoursProperty(int n, Element* element);

    void initialize();
    virtual MBXMLUtils::TiXmlElement* initializeUsingXML(MBXMLUtils::TiXmlElement *element);
    virtual MBXMLUtils::TiXmlElement* writeXMLFile(MBXMLUtils::TiXmlNode *element);
    void fromWidget(QWidget *widget);
    void toWidget(QWidget *widget);

  protected:
    std::vector<ContourOfReferenceProperty*> contour;
    Element* element;
};

class SolverTolerancesProperty : public Property {

  public:
    SolverTolerancesProperty();

    virtual MBXMLUtils::TiXmlElement* initializeUsingXML(MBXMLUtils::TiXmlElement *element);
    virtual MBXMLUtils::TiXmlElement* writeXMLFile(MBXMLUtils::TiXmlNode *element);
    void fromWidget(QWidget *widget);
    void toWidget(QWidget *widget);

  protected:
    ExtProperty projection, g, gd, gdd, la, La;
};

class SolverParametersProperty : public Property {

  public:
    SolverParametersProperty();

    virtual MBXMLUtils::TiXmlElement* initializeUsingXML(MBXMLUtils::TiXmlElement *element);
    virtual MBXMLUtils::TiXmlElement* writeXMLFile(MBXMLUtils::TiXmlNode *element);
    void fromWidget(QWidget *widget);
    void toWidget(QWidget *widget);

  protected:
    ExtProperty tolerances;
};

class GearDependencyProperty : public Property {
  public:
    GearDependencyProperty(Element* element);
    void initialize() {refBody.initialize();}
    virtual MBXMLUtils::TiXmlElement* initializeUsingXML(MBXMLUtils::TiXmlElement *element);
    virtual MBXMLUtils::TiXmlElement* writeXMLFile(MBXMLUtils::TiXmlNode *element);
    void fromWidget(QWidget *widget);
    void toWidget(QWidget *widget);
  protected:
    Element* element;
    RigidBodyOfReferenceProperty refBody;
    ExtProperty ratio;
};

class GearDependenciesProperty : public Property {

  public:
    GearDependenciesProperty(Element* element_, const std::string &xmlName_) : element(element_), xmlName(xmlName_) {}

    void initialize();
    virtual MBXMLUtils::TiXmlElement* initializeUsingXML(MBXMLUtils::TiXmlElement *element);
    virtual MBXMLUtils::TiXmlElement* writeXMLFile(MBXMLUtils::TiXmlNode *element);
    void fromWidget(QWidget *widget);
    void toWidget(QWidget *widget);

  protected:
    Element* element;
    std::string xmlName;
    std::vector<GearDependencyProperty*> refBody;

    void addDependency();
    void updateGeneralizedCoordinatesOfBodies();
};

class EmbedProperty : public Property {

  public:
    EmbedProperty(Element *element);
    virtual MBXMLUtils::TiXmlElement* initializeUsingXML(MBXMLUtils::TiXmlElement *element);
    virtual MBXMLUtils::TiXmlElement* writeXMLFile(MBXMLUtils::TiXmlNode *element);
    void fromWidget(QWidget *widget);
    void toWidget(QWidget *widget);
    std::string getFile() const {return static_cast<const FileProperty*>(href.getProperty())->getFileName();}
    bool hasCounter() const {return counterName.isActive();}
    std::string getCounterName() const {return static_cast<const TextProperty*>(counterName.getProperty())->getText();}
    bool hasParameterFile() const {return (parameterList.isActive() && static_cast<const FileProperty*>(parameterList.getProperty())->getFileName()!="");}
    std::string getParameterFile() const {return static_cast<const FileProperty*>(parameterList.getProperty())->getFileName();}

  protected:
    ExtProperty href, count, counterName, parameterList;

};

#endif
