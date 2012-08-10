#include "pendulum.h"
#include "mbsim/rigid_body.h"

#ifdef HAVE_OPENMBVCPPINTERFACE
#include <openmbvcppinterface/objbody.h>
#endif

using namespace MBSim;
using namespace fmatvec;
using namespace std;

Pendulum::Pendulum(const string &projectName) : Tree(projectName) {

  double mStab = 0.2;
  double lStab = 0.3;
  double JStab = 1.0/12.0 * mStab * lStab * lStab; 
  double a1 = -0.15*lStab;
  double a2 = 0.15*lStab;

  Vec WrOK(3);
  Vec KrKS(3);
  SymMat Theta(3);

  stab1 = new RigidBody("Stab1");
  Node* node = addObject(0,stab1);
  KrKS(0) = a1;
  SqrMat A(3);
  for(int i=0; i<3; i++)
    A(i,i) = 1;


  stab1->addFrame("Ref",-KrKS,A);
  stab1->setFrameForKinematics(stab1->getFrame("Ref"));
  stab1->setFrameOfReference(getFrame("I"));

  stab1->setqSize(1);
  stab1->setuSize(1);

  stab1->setMass(mStab);
  Theta(2,2) = JStab;
  stab1->setInertiaTensor(Theta);
  stab1->setRotation(new RotationAboutFixedAxis(Vec("[0;0;1]")));
#if HAVE_OPENMBVCPPINTERFACE
  OpenMBV::ObjBody* obj=new OpenMBV::ObjBody;
  obj->setObjFileName("objects/pendel1.obj");
  obj->setScaleFactor(0.1*0.3);
  obj->setInitialRotation(Vec("[0;0;1]")*M_PI/2);
  stab1->setOpenMBVRigidBody(obj);
#endif

  if(1) {
    stab2 = new RigidBody("Stab2");
    WrOK(0) = lStab/2;
    WrOK(2) = 0.006;
    stab1->addFrame("P",WrOK-KrKS,A);
    KrKS(0) = a2;
    stab2->addFrame("R",-KrKS,A);
    addObject(node,stab2);
    stab2->setqSize(1);
    stab2->setuSize(1);
    stab2->setFrameOfReference(stab1->getFrame("P"));
    stab2->setFrameForKinematics(stab2->getFrame("R"));
    stab2->setMass(mStab);
    Theta(2,2) = JStab;
    stab2->setInertiaTensor(Theta,stab2->getFrame("C"));
    stab2->setRotation(new RotationAboutFixedAxis(Vec("[0;0;1]")));
    stab2->setInitialGeneralizedPosition(Vec("[-1.6]"));

#if HAVE_OPENMBVCPPINTERFACE
    OpenMBV::ObjBody* obj=new OpenMBV::ObjBody;
    obj->setObjFileName("objects/pendel2.obj");
    obj->setScaleFactor(0.1*0.3);
    obj->setInitialRotation(Vec("[0;0;1]")*M_PI/2);
    stab2->setOpenMBVRigidBody(obj);
#endif
  }

}
