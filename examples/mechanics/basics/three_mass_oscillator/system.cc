#include "system.h"

#include "mbsim/rigid_body.h"
#include "mbsim/utils/utils.h"
#include "fmatvec/function.h"
#include "mbsim/spring_damper.h"
#include "mbsim/joint.h"
#include "mbsim/constitutive_laws.h"
#include "mbsim/kinetic_excitation.h"
#include "mbsim/functions/kinematic_functions.h"
#include "mbsim/functions/kinetic_functions.h"
#include "mbsim/functions/basic_functions.h"

#ifdef HAVE_OPENMBVCPPINTERFACE
#include "openmbvcppinterface/frustum.h"
#include "openmbvcppinterface/coilspring.h"
#endif

#include <iostream>

using namespace std;
using namespace MBSim;
using namespace fmatvec;

System::System(unsigned int type) : Group("System"+numtostr(int(type))) {

  // input values
  double dDisk=0.05;
  double rhoDisk=800;
  double hDisk=0.01;
  double h0Cylinder=0.1;
  double springStiffness=500;

  // derivated values
  double ADisk=1e-3*round(1e3*M_PI*dDisk*dDisk/4.);
  double mDisk=1e-3*round(1e3*ADisk*hDisk*rhoDisk);
  cout << "ADisk = " << ADisk << endl;
  cout << "mDisk = " << mDisk << endl;

  // single bodies and frames for connection
  addFrame(new FixedRelativeFrame("fK1", 1.*h0Cylinder*Vec("[0;1;0]"), SqrMat(3,EYE)));
  addFrame(new FixedRelativeFrame("fK2", 2.*h0Cylinder*Vec("[0;1;0]"), SqrMat(3,EYE)));
  addFrame(new FixedRelativeFrame("fK3", 3.*h0Cylinder*Vec("[0;1;0]"), SqrMat(3,EYE)));

  RigidBody *k1 = new RigidBody("K1");
  k1->setFrameForKinematics(k1->getFrame("C"));
  k1->setMass(3.*mDisk);
  k1->setInertiaTensor(0.001*SymMat(3,EYE));
  k1->setTranslation(new TranslationAlongAxesXYZ<VecV>);
  k1->addFrame(new FixedRelativeFrame("fK2", 1.*h0Cylinder*Vec("[0;1;0]"), SqrMat(3,EYE)));
  k1->addFrame(new FixedRelativeFrame("fK3", 2.*h0Cylinder*Vec("[0;1;0]"), SqrMat(3,EYE)));
  k1->addFrame(new FixedRelativeFrame("fTop", .5*hDisk*Vec("[0;1;0]"), SqrMat(3,EYE)));
  k1->addFrame(new FixedRelativeFrame("fBottom", -.5*hDisk*Vec("[0;1;0]"), SqrMat(3,EYE)));

  RigidBody *k2 = new RigidBody("K2");
  k2->setFrameForKinematics(k2->getFrame("C"));
  k2->setMass(5.*mDisk);
  k2->setInertiaTensor(0.001*SymMat(3,EYE));
  k2->setTranslation(new TranslationAlongAxesXYZ<VecV>);
  k2->addFrame(new FixedRelativeFrame("fK3", 1.*h0Cylinder*Vec("[0;1;0]"), SqrMat(3,EYE)));
  k2->addFrame(new FixedRelativeFrame("fTop", .5*hDisk*Vec("[0;1;0]"), SqrMat(3,EYE)));
  k2->addFrame(new FixedRelativeFrame("fBottom", -.5*hDisk*Vec("[0;1;0]"), SqrMat(3,EYE)));

  RigidBody *k3 = new RigidBody("K3");
  k3->setFrameForKinematics(k3->getFrame("C"));
  k3->setMass(7.*mDisk);
  k3->setInertiaTensor(0.001*SymMat(3,EYE));
  k3->setTranslation(new TranslationAlongAxesXYZ<VecV>);
  k3->addFrame(new FixedRelativeFrame("fTop", .5*hDisk*Vec("[0;1;0]"), SqrMat(3,EYE)));
  k3->addFrame(new FixedRelativeFrame("fBottom", -.5*hDisk*Vec("[0;1;0]"), SqrMat(3,EYE)));

  // parametrisation of relative kinematics
  if(type==1) { // K2 child of K1, K3 child of K2
    addObject(k1);
    addObject(k2);
    addObject(k3);
    k1->setFrameOfReference(getFrame("fK1"));
    k2->setFrameOfReference(k1->getFrame("fK2"));
    k3->setFrameOfReference(k2->getFrame("fK3"));
  }
  else if(type==2) { // K2 child of K1, K3 child of K1
    addObject(k1);
    addObject(k2);
    addObject(k3);
    k1->setFrameOfReference(getFrame("fK1"));
    k2->setFrameOfReference(k1->getFrame("fK2"));
    k3->setFrameOfReference(k1->getFrame("fK3"));
  }
  else { // absolute parametrisation
    addObject(k1);
    addObject(k2);
    addObject(k3);
    k1->setFrameOfReference(getFrame("fK1"));
    k2->setFrameOfReference(getFrame("fK2"));
    k3->setFrameOfReference(getFrame("fK3"));
  }

  // connection between bodies
  SpringDamper *sp01 = new SpringDamper("S01");
  addLink(sp01);
  sp01->setForceFunction(new LinearSpringDamperForce(springStiffness,0.05*springStiffness,h0Cylinder-.5*hDisk));
  sp01->connect(getFrame("I"), k1->getFrame("fBottom"));

  SpringDamper *sp12 = new SpringDamper("S12");
  addLink(sp12);
  sp12->setForceFunction(new LinearSpringDamperForce(springStiffness,0.05*springStiffness,h0Cylinder-hDisk));
  sp12->connect(k1->getFrame("fTop"), k2->getFrame("fBottom"));

  SpringDamper *sp23 = new SpringDamper("S23");
  addLink(sp23);
  sp23->setForceFunction(new LinearSpringDamperForce(springStiffness,0.05*springStiffness,h0Cylinder-hDisk));
  sp23->connect(k2->getFrame("fTop"), k3->getFrame("fBottom"));

  // joint between the upper two bodies
  Joint *j23;
  j23 = new Joint("xxxJ23");
  j23->setForceDirection(Mat("[0;1;0]"));
  j23->setForceLaw(new BilateralConstraint());
  j23->connect(k2->getFrame("fTop"), k3->getFrame("fBottom"));
  addLink(j23);

  KineticExcitation *a = new KineticExcitation("Anregung");
  addLink(a);
  a->setForceDirection(Vec("[0;1;0]"));
  a->setForceFunction(new ConstantFunction<VecV(double)>(-10));
  a->connect(k3->getFrame("fTop"));

#ifdef HAVE_OPENMBVCPPINTERFACE
  boost::shared_ptr<OpenMBV::Frustum> k1Visu = OpenMBV::ObjectFactory::create<OpenMBV::Frustum>();
  k1Visu->setBaseRadius(dDisk/2.);
  k1Visu->setTopRadius(dDisk/2.);
  k1Visu->setInnerBaseRadius(0);
  k1Visu->setInnerTopRadius(0);
  k1Visu->setHeight(hDisk);
  k1Visu->setInitialRotation(M_PI/2., 0, 0);
  k1Visu->setInitialTranslation(0, -hDisk/2., 0);
  k1Visu->setDiffuseColor(120./360.,1,1);
  k1->setOpenMBVRigidBody(k1Visu);

  boost::shared_ptr<OpenMBV::Frustum> k2Visu = OpenMBV::ObjectFactory::create<OpenMBV::Frustum>();
  k2Visu->setBaseRadius(dDisk/2.);
  k2Visu->setTopRadius(dDisk/2.);
  k2Visu->setInnerBaseRadius(0);
  k2Visu->setInnerTopRadius(0);
  k2Visu->setHeight(hDisk);
  k2Visu->setInitialRotation(M_PI/2., 0, 0);
  k2Visu->setInitialTranslation(0, -hDisk/2., 0);
  k2Visu->setDiffuseColor(120./360.,1,1);
  k2->setOpenMBVRigidBody(k2Visu);

  boost::shared_ptr<OpenMBV::Frustum> k3Visu = OpenMBV::ObjectFactory::create<OpenMBV::Frustum>();
  k3Visu->setBaseRadius(dDisk/2.);
  k3Visu->setTopRadius(dDisk/2.);
  k3Visu->setInnerBaseRadius(0);
  k3Visu->setInnerTopRadius(0);
  k3Visu->setHeight(hDisk);
  k3Visu->setInitialRotation(M_PI/2., 0, 0);
  k3Visu->setInitialTranslation(0, -hDisk/2., 0);
  k3Visu->setDiffuseColor(120./360.,1,1);
  k3->setOpenMBVRigidBody(k3Visu);

  sp01->enableOpenMBVCoilSpring(_springRadius=.75*.5*dDisk,_crossSectionRadius=.1*.25*dDisk,_numberOfCoils=5);

  sp12->enableOpenMBVCoilSpring(_springRadius=.75*.5*dDisk,_crossSectionRadius=.1*.25*dDisk,_numberOfCoils=5);

  sp23->enableOpenMBVCoilSpring(_springRadius=.75*.5*dDisk,_crossSectionRadius=.1*.25*dDisk,_numberOfCoils=5);
#endif

  if(type==1) 
    setPosition(1.5*dDisk*Vec("[1;0;0]")+1.5*dDisk*Vec("[0;0;-1]"));
  else if (type==2)
    setPosition(3.*dDisk*Vec("[1;0;0]")+1.5*dDisk*Vec("[0;0;-1]"));
  else 
    setPosition(0*h0Cylinder*Vec("[1;0;0]")+1.5*dDisk*Vec("[0;0;-1]"));
}

