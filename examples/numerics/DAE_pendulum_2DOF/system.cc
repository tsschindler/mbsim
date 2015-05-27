#include "../DAE_pendulum_2DOF/system.h"

#include "mbsim/rigid_body.h"
#include "mbsim/joint.h"
#include "mbsim/spring_damper.h"
#include "mbsim/environment.h"
#include "mbsim/contours/sphere.h"
#include "mbsim/contact.h"
#include "mbsim/constitutive_laws.h"
//#include "mbsim/utils/symbolic_function.h"

#include "mbsim/functions/basic_functions.h"
#include "mbsim/functions/symbolic_functions.h"
#include "mbsim/functions/kinematic_functions.h"
#include "mbsim/functions/kinetic_functions.h"
#include <casadi/symbolic/sx/sx_tools.hpp>


#ifdef HAVE_OPENMBVCPPINTERFACE
#include "openmbvcppinterface/coilspring.h"
#include <openmbvcppinterface/cuboid.h>
#endif

using namespace MBSim;
using namespace fmatvec;
using namespace std;
using namespace CasADi;


System::System(const string &projectName) :
    DynamicSystemSolver(projectName) {
  // acceleration of gravity
  Vec grav(3);
  grav(1)=-9.81;
  MBSimEnvironment::getInstance()->setAccelerationOfGravity(grav);

  double length_crank = 0.2;
//  double length_spring = 0.1;
//  double radius_mass = length_crank / 10.;
  double mass_crank = 1; //0.038; // m1
//  double mass_mass1 = 2;
  // crank
  RigidBody *crank = new RigidBody("Crank");
  this->addObject(crank);
  // kinematics
  crank->setFrameOfReference(getFrameI());
  Vec3 kinematicsFrameCrankPos;
  kinematicsFrameCrankPos(0) = -length_crank / 2.;
  FixedRelativeFrame * kinematicsFrameCrank = new FixedRelativeFrame("LoadFrame", kinematicsFrameCrankPos);
  FixedRelativeFrame * crankToSpring = new FixedRelativeFrame("crankToSpring", -kinematicsFrameCrankPos);
  crank->addFrame(kinematicsFrameCrank);
  crank->addFrame(crankToSpring);
  //crank->setFrameForKinematics(kinematicsFrameCrank);
  kinematicsFrameCrank->enableOpenMBV(0.5e-1);
  crankToSpring->enableOpenMBV(0.5e-1);
  crank->getFrameC()->enableOpenMBV(0.7e-1);
  
  crank->setMass(mass_crank);
  SymMat inertia_crank(3, INIT, 0.);
  inertia_crank(0, 0) = 1.; // DUMMY
  inertia_crank(1, 1) = 1.; // DUMMY
  inertia_crank(2, 2) = 1; // J1
  crank->setInertiaTensor(inertia_crank);

  crank->setTranslation(new TranslationAlongAxesXY<VecV>);
//  crank->setRotation(new RotationAboutZAxis<VecV>);
//  Vec q0(3);
//  q0(0) = 0.5 * length_crank * cos(M_PI /6);
//  q0(1) = -0.5 * length_crank * sin(M_PI /6);
//  q0(2) = -M_PI /6;
  Vec q0(2);
  q0(0) = 0.5 * length_crank * cos(-M_PI /6);
  q0(1) = 0.5 * length_crank * sin(-M_PI /6);
  crank->setInitialGeneralizedPosition(q0);

  Joint *joint1 = new Joint("Gelenk1");
  addLink(joint1);
  joint1->setForceDirection(Mat("[1, 0; 0, 1; 0, 0]"));
//  joint1->setMomentDirection(Mat("[1, 0; 0, 1; 0, 0]"));
  joint1->connect(getFrame("I"),crank->getFrame("LoadFrame"));
  joint1->setForceLaw(new BilateralConstraint);
//  joint1->setMomentLaw(new BilateralConstraint);

//  SX t("t");
//  SX fexp2 = log(cosh(t));
//  SXFunction foo2(t,fexp2);
//  SymbolicFunction<double(double)> *f2 = new SymbolicFunction<double(double)>(foo2);
//  crank->setRotation(new NestedFunction<RotMat3(double(double))>(new RotationAboutFixedAxis<double>("[0;0;1]"), f2));


  // visualisation
#ifdef HAVE_OPENMBVCPPINTERFACE
  boost::shared_ptr<OpenMBV::Cuboid> openMBVCrank = OpenMBV::ObjectFactory::create<OpenMBV::Cuboid>();
  openMBVCrank->setLength(length_crank, 0.05, 0.05);
  openMBVCrank->setDiffuseColor(0.5, 0, 0);
  openMBVCrank->setTransparency(0.5);
  crank->setOpenMBVRigidBody(openMBVCrank);
#endif


}

