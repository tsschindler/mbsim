#include "pendulum.h"
#include "body_rigid.h"
#include "objobject.h"

using namespace AMVis;

Pendulum::Pendulum(const string &projectName) : Tree(projectName) {

  double mStab = 0.2;
  double lStab = 0.3;
  double JStab = 1.0/12.0 * mStab * lStab * lStab; 
  double a1 = -0.15*lStab;
  double a2 = 0.15*lStab;

  Vec WrOK(3);
  Vec KrKS(3);
  SymMat Theta(3);

  stab1 = new BodyRigid("Stab1");
  addObject(stab1);
  KrKS(0) = a1;
  SqrMat A(3);
  for(int i=0; i<3; i++)
    A(i,i) = 1;


  stab1->addCoordinateSystem("Ref",-KrKS,A);
  stab1->setRefCoordinateSystem(stab1->getCoordinateSystem("Ref"));

  stab1->setqSize(1);
  stab1->setuSize(1);

  stab1->setMass(mStab);
  Theta(2,2) = JStab;
  stab1->setInertia(Theta);
  stab1->setfAPK(new RotationAxis(Vec("[0;0;1]")));
  ObjObject * obj = new ObjObject(stab1->getFullName(),1,false);
  obj->setObjFilename("objects/pendel1.obj");
  stab1->setAMVisBody(obj);
  obj->setScaleFactor(0.1*0.3);
  obj -> setInitialRotation(0,0,M_PI/2);
  obj->setCalculationOfNormals(3);
  obj->setVertexEPS(1e-5);
  obj-> setNormalEPS(1e-5);
  obj-> setAngleEPS(M_PI*2/9);

  if(1) {
  stab2 = new BodyRigid("Stab2");
  WrOK(0) = lStab/2;
  WrOK(2) = 0.006;
  stab1->addCoordinateSystem("P",WrOK-KrKS,A);
  KrKS(0) = a2;
  stab2->addCoordinateSystem("R",-KrKS,A);
  addObject(stab2);
  stab2->setqSize(1);
  stab2->setuSize(1);
  stab2->setParentCoordinateSystem(stab1->getCoordinateSystem("P"));
  stab2->setRefCoordinateSystem(stab2->getCoordinateSystem("R"));
  stab2->setMass(mStab);
  Theta(2,2) = JStab;
  stab2->setInertia(Theta,stab2->getCoordinateSystem("COG"));
  stab2->setfAPK(new RotationAxis(Vec("[0;0;1]")));
  stab2->setq0(Vec("[-1.6]"));

  obj = new ObjObject(stab2->getFullName(),1,false);
  obj->setObjFilename("objects/pendel2.obj");
  stab2->setAMVisBody(obj);
  obj->setScaleFactor(0.1*0.3);
  obj -> setInitialRotation(0,0,M_PI/2);
  obj->setCalculationOfNormals(3);
  obj->setVertexEPS(1e-5);
  obj-> setNormalEPS(1e-5);
  obj-> setAngleEPS(M_PI*2/9);
  }

}

