#include "system.h"
#include "mbsim/rigid_body.h"
#include "mbsim/contour.h"
#include "mbsim/constitutive_laws.h"
#include "mbsim/contact.h"
#include "mbsim/load.h"

#ifdef HAVE_OPENMBVCPPINTERFACE
#include <openmbvcppinterface/invisiblebody.h>
#endif

using namespace MBSim;
using namespace fmatvec;
using namespace std;

extern bool rigidContact;

System::System(const string &projectName) : DynamicSystemSolver(projectName) {
  // Gravitation
  Vec grav(3);
  grav(1)=-9.81;
  setAccelerationOfGravity(grav);
  // Parameters
  double l = 0.8;              		
  double h =  0.02;
  double d = 0.1;
  double m = 0.7;
  SymMat Theta(3);
  Theta(1,1) = m*d*d/2.;
  Theta(2,2) = Theta(1,1);
  double alpha = 3.0 * M_PI/180.; 
  double deltax = 0.2;           
  double mu  = 0.3;

  RigidBody* body = new RigidBody("Rod");
  addObject(body);

  body->setFrameOfReference(getFrame("I"));
  body->setFrameForKinematics(body->getFrame("C"));
  body->setMass(m);
  body->setInertiaTensor(Theta);
  body->setTranslation(new LinearTranslation("[1, 0; 0, 1; 0, 0]"));
  body->setRotation(new RotationAboutFixedAxis(Vec("[0;0;1]")));

  // Initial translation and rotation
  Vec q0(3);
  q0(1) = .3;
  body->setq0(q0);
  Vec u0(3);
  u0(2) = -M_PI*10;
  body->setu0(u0);

  CircleSolid *circlecontour=new CircleSolid("Circle",d);
  body->addContour(circlecontour,Vec(3),SqrMat(3,EYE));
#ifdef HAVE_OPENMBVCPPINTERFACE
  circlecontour->enableOpenMBV();
#endif

  double phi = M_PI*0.6; 
  SqrMat A(3);
  A(0,0) = cos(phi);
  A(0,1) = -sin(phi);
  A(1,1) = cos(phi);
  A(1,0) = sin(phi);
  A(2,2) = 1;
  addContour(new Line("Line"),Vec(3),A);

  Contact *rc = new Contact("Contact");
  rc->connect(getContour("Line"),body->getContour("Circle"));
  addLink(rc);
  if(rigidContact) {
    rc->setContactForceLaw(new UnilateralConstraint);
    rc->setContactImpactLaw(new UnilateralNewtonImpact);
    rc->setFrictionForceLaw(new PlanarCoulombFriction(mu));
    rc->setFrictionImpactLaw(new PlanarCoulombImpact(mu));
  } 
  else {
    rc->setContactForceLaw(new LinearRegularizedUnilateralConstraint(1e5,1e4));
    rc->setFrictionForceLaw(new LinearRegularizedPlanarCoulombFriction(mu));
  }

#ifdef HAVE_OPENMBVCPPINTERFACE
  body->getFrame("C")->enableOpenMBV(1.2*d);
  OpenMBV::InvisibleBody* dummy=new OpenMBV::InvisibleBody;
  body->setOpenMBVRigidBody(dummy);
#endif

}

