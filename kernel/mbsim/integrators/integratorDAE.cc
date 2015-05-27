/*
 * IntegratorDAE.cc
 *
 *  Created on: Mar 30, 2015
 *      Author: zwang
 */
#include<config.h>
#include "integratorDAE.h"
#include "mbsim/element.h"

using namespace fmatvec;
using namespace MBSim;
using namespace MBXMLUtils;
using namespace xercesc;

namespace MBSimIntegrator {
  
  IntegratorDAE::IntegratorDAE() :
    aTol(1, INIT, 1e-6), rTol(1, INIT, 1e-4), iTol(0), dt0(0), maxSteps(100000), dtMax(0), DAEIndex(1), FlagErrorTest(1), useExternalJac(0), t(0.), rWork(1, INIT, 0.0), lrWork(0), iWork(1, INIT, 0), liWork(0), rPar(1, INIT, 0.0), iPar(1, INIT, 0), iout(0), idid(0), s0(0.0), time(0.0), FlagPlotIntegrator(true), FlagPlotIntegrationSum(true), FlagCoutInfo(true), Y(1, INIT, 0.0), scaling(1), H(1e-3), rho(1.)  {
    // TODO Auto-generated constructor stub
    
  }
  
  IntegratorDAE::~IntegratorDAE() {
    // TODO Auto-generated destructor stub
  }

  void IntegratorDAE::integrate(DynamicSystemSolver& system_) {
    debugInit();
    preIntegrate(system_);
    subIntegrate(system_, tEnd);
    postIntegrate(system_);
  }

  void IntegratorDAE::initializeUsingXML(DOMElement *element) {
    Integrator::initializeUsingXML(element);
    DOMElement *e;
    e = E(element)->getFirstElementChildNamed(MBSIMINT % "absoluteTolerance");
    if (e)
      setAbsoluteTolerance(Element::getVec(e));
    e = E(element)->getFirstElementChildNamed(MBSIMINT % "absoluteToleranceScalar");
    if (e)
      setAbsoluteTolerance(Element::getDouble(e));
    e = E(element)->getFirstElementChildNamed(MBSIMINT % "relativeTolerance");
    if (e)
      setRelativeTolerance(Element::getVec(e));
    e = E(element)->getFirstElementChildNamed(MBSIMINT % "relativeToleranceScalar");
    if (e)
      setRelativeTolerance(Element::getDouble(e));
    e = E(element)->getFirstElementChildNamed(MBSIMINT % "initialStepSize");
    setInitialStepSize(Element::getDouble(e));
    e = E(element)->getFirstElementChildNamed(MBSIMINT % "maximalStepSize");
    setMaximalStepSize(Element::getDouble(e));
    e = E(element)->getFirstElementChildNamed(MBSIMINT % "maximalNumberOfSteps");
    if (e)
      setMaxStepNumber(Element::getInt(e));
  }
} /* namespace MBSimIntegrator */
