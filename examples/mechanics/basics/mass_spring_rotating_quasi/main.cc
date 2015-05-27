#include <mbsim/integrators/integrators.h>
#include "system.h"

using namespace std;
using namespace MBSim;
using namespace MBSimIntegrator;

int main(int argc, char* argv[]) {

  string inputFileName("config.dat");
  ifstream configFile(inputFileName.c_str());
  if (!configFile.is_open()) {
    cout << "Can not open file " << "config.dat" << endl;
    throw 1;
  }

  string IntegratorType;
  getline(configFile, IntegratorType);
  configFile.close();

  if (IntegratorType == "TS") {
    // build single modules
    System *sys = new System("TS");

    // add modules to overall dynamical system
    sys->initialize();
    TimeSteppingIntegrator integrator;
    cout << "TimeSteppingIntegrator. \n";
    integrator.setStepSize(1e-3);
    integrator.setEndTime(4.0);
    integrator.setPlotStepSize(1e-3);
    integrator.integrate(*sys);
    cout << "finished" << endl;
    delete sys;
  }
  else if (IntegratorType == "QS") {
    // build single modules
    System *sys = new System("QS");
    // add modules to overall dynamical system
    sys->initialize();
    QuasiStaticIntegrator integrator;
    cout << "QuasiStaticIntegrator. \n";
    integrator.setStepSize(1e-3);
    integrator.setEndTime(4.0);
    integrator.setPlotStepSize(1e-3);
    integrator.integrate(*sys);
    integrator.sethTolerance(1e-6);
    integrator.setgTolerance(1e-6);
    cout << "finished" << endl;
    delete sys;
  }
  else {
    cout << "Please assign the integrator type in the config.dat file." << endl;
    throw 1;
  }

  return 0;
}

