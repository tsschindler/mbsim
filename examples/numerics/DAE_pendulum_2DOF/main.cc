#include <mbsim/integrators/integrators.h>

#include "../DAE_pendulum_2DOF/system.h"

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
  double tEnd = 4;
  string IntegratorType;
  getline(configFile, IntegratorType);
  configFile.close();
  // if any argument is given, then use time stepping integrator, otherwise daskr integrator is used.
  if (IntegratorType == "TS") {
    // build single modules
    System *sys = new System("TS");

    // add modules to overall dynamical system
    sys->initialize();
    TimeSteppingIntegrator integrator;
    cout << "TimeSteppingIntegrator. \n";
    integrator.setStepSize(1e-3);
    integrator.setEndTime(tEnd);
    integrator.setPlotStepSize(1e-3);
    integrator.integrate(*sys);
    cout << "finished" << endl;
    delete sys;
  }
  else if (IntegratorType == "DASKR") {
    // build single modules
    System *sys = new System("DASKR");
    // add modules to overall dynamical system
    sys->initialize();
    DASKRIntegrator integrator;
    cout << "DASKR Integrator. \n";
    integrator.setEndTime(tEnd);
    integrator.setPlotStepSize(1e-3);
    integrator.integrate(*sys);
    cout << "finished" << endl;
    delete sys;
  }
  else if (IntegratorType == "RADAU5") {
    // build single modules
    System *sys = new System("RADAU5");
    // add modules to overall dynamical system
    sys->initialize();
    RADAU5Integrator integrator;
    cout << "RADAU5 Integrator. \n";
    integrator.setEndTime(tEnd);
    integrator.setPlotStepSize(1e-3);
    integrator.integrate(*sys);
    cout << "finished" << endl;
    delete sys;
  }
  else if (IntegratorType == "PHEM56") {
    // build single modules
    System *sys = new System("PHEM56");
    // add modules to overall dynamical system
    sys->initialize();
    PHEM56Integrator integrator;
    cout << "PHEM56 Integrator. \n";
    integrator.setEndTime(tEnd);
    integrator.setPlotStepSize(1e-3);
    integrator.integrate(*sys);
    cout << "finished" << endl;
    delete sys;
  }
  else {
    cout << "Please assign the integrator type in the config.dat file." << endl;
    throw 1;
  }

  return 0;
}

