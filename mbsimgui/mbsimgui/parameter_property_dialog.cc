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
#include "parameter_property_dialog.h"
#include "basic_widgets.h"
#include "variable_widgets.h"
#include "extended_widgets.h"
#include "parameter.h"

using namespace std;

namespace MBSimGUI {

  ParameterPropertyDialog::ParameterPropertyDialog(Parameter *parameter_, QWidget *parent, Qt::WindowFlags f, bool readOnly) : PropertyDialog(parent,f), parameter(parameter_) {
    addTab("General");
    name=new ExtWidget("Name",new TextWidget("",readOnly));
    addToTab("General",name);
  }

  void ParameterPropertyDialog::toWidget(Parameter *parameter) {
    parameter->name.toWidget(name);
  }

  void ParameterPropertyDialog::fromWidget(Parameter *parameter) {
    parameter->name.fromWidget(name);
  }

  StringParameterPropertyDialog::StringParameterPropertyDialog(StringParameter *parameter, QWidget *parent, Qt::WindowFlags f) : ParameterPropertyDialog(parameter,parent,f) {
    value = new ExtWidget("Value",new ChoiceWidget2(new ScalarWidgetFactory("0",vector<QStringList>(2,QStringList()))));
    addToTab("General", value);
  }

  void StringParameterPropertyDialog::toWidget(Parameter *parameter) {
    ParameterPropertyDialog::toWidget(parameter);
    static_cast<StringParameter*>(parameter)->value.toWidget(value);
  }

  void StringParameterPropertyDialog::fromWidget(Parameter *parameter) {
    ParameterPropertyDialog::fromWidget(parameter);
    static_cast<StringParameter*>(parameter)->value.fromWidget(value);
    parameter->setValue(static_cast<PhysicalVariableProperty*>(static_cast<ChoiceProperty2*>(static_cast<StringParameter*>(parameter)->value.getProperty())->getProperty())->getValue());
  }

  ScalarParameterPropertyDialog::ScalarParameterPropertyDialog(ScalarParameter *parameter, QWidget *parent, Qt::WindowFlags f) : ParameterPropertyDialog(parameter,parent,f) {
    value = new ExtWidget("Value",new ChoiceWidget2(new ScalarWidgetFactory("0",vector<QStringList>(2,QStringList()))));
    addToTab("General", value);
  }

  void ScalarParameterPropertyDialog::toWidget(Parameter *parameter) {
    ParameterPropertyDialog::toWidget(parameter);
    static_cast<ScalarParameter*>(parameter)->value.toWidget(value);
  }

  void ScalarParameterPropertyDialog::fromWidget(Parameter *parameter) {
    ParameterPropertyDialog::fromWidget(parameter);
    static_cast<ScalarParameter*>(parameter)->value.fromWidget(value);
    parameter->setValue(static_cast<PhysicalVariableProperty*>(static_cast<ChoiceProperty2*>(static_cast<ScalarParameter*>(parameter)->value.getProperty())->getProperty())->getValue());
  }

  VectorParameterPropertyDialog::VectorParameterPropertyDialog(VectorParameter *parameter, QWidget *parent, Qt::WindowFlags f) : ParameterPropertyDialog(parameter,parent,f) {
    value = new ExtWidget("Value",new ChoiceWidget2(new VecWidgetFactory(3,vector<QStringList>(3,QStringList()))));
    addToTab("General", value);
  }

  void VectorParameterPropertyDialog::toWidget(Parameter *parameter) {
    ParameterPropertyDialog::toWidget(parameter);
    static_cast<VectorParameter*>(parameter)->value.toWidget(value);
  }

  void VectorParameterPropertyDialog::fromWidget(Parameter *parameter) {
    ParameterPropertyDialog::fromWidget(parameter);
    static_cast<VectorParameter*>(parameter)->value.fromWidget(value);
    parameter->setValue(static_cast<PhysicalVariableProperty*>(static_cast<ChoiceProperty2*>(static_cast<VectorParameter*>(parameter)->value.getProperty())->getProperty())->getValue());
  }

  MatrixParameterPropertyDialog::MatrixParameterPropertyDialog(MatrixParameter *parameter,QWidget *parent, Qt::WindowFlags f) : ParameterPropertyDialog(parameter,parent,f) {
    value = new ExtWidget("Value",new ChoiceWidget2(new MatRowsColsVarWidgetFactory(3,3)));
    addToTab("General", value);
  }

  void MatrixParameterPropertyDialog::toWidget(Parameter *parameter) {
    ParameterPropertyDialog::toWidget(parameter);
    static_cast<MatrixParameter*>(parameter)->value.toWidget(value);
  }

  void MatrixParameterPropertyDialog::fromWidget(Parameter *parameter) {
    ParameterPropertyDialog::fromWidget(parameter);
    static_cast<MatrixParameter*>(parameter)->value.fromWidget(value);
    parameter->setValue(static_cast<PhysicalVariableProperty*>(static_cast<ChoiceProperty2*>(static_cast<MatrixParameter*>(parameter)->value.getProperty())->getProperty())->getValue());

  }

  SearchPathParameterPropertyDialog::SearchPathParameterPropertyDialog(SearchPathParameter *parameter,QWidget *parent, Qt::WindowFlags f) : ParameterPropertyDialog(parameter,parent,f,true) {
    value = new ExtWidget("Path",new FileWidget("Path","",3));
    addToTab("General", value);
  }

  void SearchPathParameterPropertyDialog::toWidget(Parameter *parameter) {
    ParameterPropertyDialog::toWidget(parameter);
    static_cast<FileWidget*>(value->getWidget())->setFile(QString::fromStdString(parameter->getValue()));
  }

  void SearchPathParameterPropertyDialog::fromWidget(Parameter *parameter) {
    ParameterPropertyDialog::fromWidget(parameter);
    parameter->setValue(static_cast<FileWidget*>(value->getWidget())->getFile().toStdString());
  }

}
