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

#ifndef _INTEGRATOR__H_
#define _INTEGRATOR__H_

#define MBSIMINTNS "{http://mbsim.berlios.de/MBSimIntegrator}"

#include <QtGui/QTreeWidgetItem>
#include "mbsimguitinyxml/tinyxml-src/tinyxml.h"
#include "mbsimguitinyxml/tinyxml-src/tinynamespace.h"
#include "editors.h"
#include "utils.h"
#include <string>
#include <set>

class PropertyDialog;

class Integrator : public QObject, public QTreeWidgetItem {
  Q_OBJECT
  friend class Editor;
  friend class MainWindow;
  protected:
    QString newName(const QString &type);
    bool drawThisPath;
    std::string iconFile;
    bool searchMatched;
    PropertyDialog *properties;
    QMenu *contextMenu;
    DoubleEditor *startTime, *endTime, *plotStepSize;
  public:
    Integrator(const QString &str, QTreeWidgetItem *parentItem, int ind);
    virtual ~Integrator();
    std::string &getIconFile() { return iconFile; }
    virtual QString getInfo();
    virtual void initializeUsingXML(TiXmlElement *element);
    virtual TiXmlElement* writeXMLFile(TiXmlNode *element);
    static Integrator* readXMLFile(const QString &filename, QTreeWidgetItem *parent);
    virtual void writeXMLFile(const QString &name);
    virtual void writeXMLFile() { writeXMLFile(getType()); }
    virtual QString getType() const { return "Integrator"; }
    QMenu* getContextMenu() { return contextMenu; }
    PropertyDialog* getPropertyDialog() { return properties; }
    void setEndTime(double t) {endTime->setValue(t);}
    void setPlotStepSize(double dt) {plotStepSize->setValue(dt);}
    public slots:
      void saveAs();
};

class DOPRI5Integrator : public Integrator {
  public:
    DOPRI5Integrator(const QString &str, QTreeWidgetItem *parentItem, int ind);
    virtual void initializeUsingXML(TiXmlElement *element);
    virtual TiXmlElement* writeXMLFile(TiXmlNode *element);
    virtual QString getType() const { return "DOPRI5Integrator"; }
  protected:
    DoubleEditor *absTol, *relTol, *initialStepSize, *maximalStepSize;
    IntEditor *maxSteps;
};

class LSODEIntegrator : public Integrator {
  public:
    LSODEIntegrator(const QString &str, QTreeWidgetItem *parentItem, int ind);
    virtual void initializeUsingXML(TiXmlElement *element);
    virtual TiXmlElement* writeXMLFile(TiXmlNode *element);
    virtual QString getType() const { return "LSODEIntegrator"; }
  protected:
    DoubleEditor *absTol, *relTol, *initialStepSize, *maximalStepSize, *minimalStepSize;
    IntEditor *maxSteps;
    BoolEditor *stiff;
};

#endif