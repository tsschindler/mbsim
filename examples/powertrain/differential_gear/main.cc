#include "system.h"
#include <mbsim/integrators/integrators.h>

using namespace std;
using namespace MBSim;
using namespace MBSimIntegrator;

int main (int argc, char* argv[])
{
  // Einzelne Bausteine des MKS erschaffen
  DynamicSystemSolver *sys = new System("TS");

  // Bausteine zum Gesamtsystem zusammenfuegen (zu einem DGL-System) 
  sys->initialize();

  // LSODEIntegrator integrator;
  // RKSuiteIntegrator integrator;
  // RADAU5Integrator integrator;
  //TimeSteppingIntegrator integrator;
  //integrator.setdt(1e-4);
  //
  DOPRI5Integrator integrator;

  integrator.setEndTime(10.0);
  integrator.setPlotStepSize(1e-2);

  // Vec z(4);
  // z(0) = 0.13;
  // z(1) = -1.13;
  // z(2) = -1.1;
  // z(3) = 3.1;

  // cout << sys->zdot(z,0)<< endl;
  // cout << sys->getM()<<endl;
  // cout << sys->geth()<<endl;
  // throw 5;

  integrator.integrate(*sys);
  cout << "finished"<<endl;
  delete sys;

  return 0;

}

