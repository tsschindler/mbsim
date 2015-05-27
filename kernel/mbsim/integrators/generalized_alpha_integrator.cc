/* Copyright (C) 2004-2014 MBSim Development Team
 *
 * This library is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU Lesser General Public 
 * License as published by the Free Software Foundation; either 
 * version 2.1 of the License, or (at your option) any later version. 
 *  
 * This library is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
 * Lesser General Public License for more details. 
 *  
 * You should have received a copy of the GNU Lesser General Public 
 * License along with this library; if not, write to the Free Software 
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 *
 * Contact: martin.o.foerg@googlemail.com
 */

#include<config.h>
#include<mbsim/dynamic_system_solver.h>
#include "generalized_alpha_integrator.h"
#include <mbsim/utils/nonlinear_algebra.h>
#include <mbsim/numerics/nonlinear_algebra/multi_dimensional_newton_method.h>

#include <time.h>
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>

#ifndef NO_ISO_14882
using namespace std;
#endif

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;

using namespace fmatvec;
using namespace MBSim;
using namespace MBXMLUtils;
using namespace xercesc;

namespace MBSimIntegrator {

  MBSIM_OBJECTFACTORY_REGISTERXMLNAME(GeneralizedAlphaIntegrator, MBSIMINT % "GeneralizedAlphaIntegrator")

  GeneralizedAlphaIntegrator::GeneralizedAlphaIntegrator() :
      dt(1e-3), t(0.), tPlot(0.), gTol(1e-6), hTol(1e-4), iter(0), step(0), integrationSteps(0), maxIter(0), sumIter(0), updateJacobianEvery(1), s0(0.), time(0.), stepPlot(0), rho_inf(0.8) {
  }

  void GeneralizedAlphaIntegrator::integratorSetups(DynamicSystemSolver& system_) {
    zSize = system_.getzSize();  // zSize = qSize + uSize[0] + xSize;
    qSize = system_.getqSize();
    uSize = system_.getuSize();
    laBiSize = system_.getlaSize(); // todo: check whether need to different laBiSize and laUniSize
    YSize = zSize + laBiSize;
    //    if (DAEIndex==21) YSize += laBiSize;
    t = tStart;

    Y.resize(YSize, INIT, 0.0);
    ud.resize(uSize, INIT, 0.0);
    a.resize(uSize, INIT, 0.0);

    Index Iq(0, qSize - 1);
    Index Iu(qSize, qSize + uSize - 1);
    Index Ix(qSize + uSize, zSize - 1);
    Index Iz(0, zSize - 1);
    q >> Y(Iq);
    u >> Y(Iu);
    x >> Y(Ix);
    z >> Y(Iz);

    if (z0.size())
      z = z0;
    else
      system_.initz(z);

    //todo: initialize ud, a ?? M^(-1)(h + W * la);

    // aTol and rTol can be chosen both scalars or else both arrays of length NEQ.
    assert(aTol.size() == rTol.size());

    // alpha parameters
    alpha_m = (2 * rho_inf - 1) / (rho_inf + 1);
    alpha_f = rho_inf / (rho_inf + 1);
    gamma = 0.5 - alpha_m + alpha_f;
    beta = 0.25 * (1 - alpha_m + alpha_f) * (1 - alpha_m + alpha_f);
    beta_stroke = (1 - alpha_m) / (dt * dt * beta * (1 - alpha_f));
    gamma_stroke = gamma / (dt * beta);

  }
  void GeneralizedAlphaIntegrator::preIntegrate(DynamicSystemSolver& system_) {
    // initialisation
    system = &system_;
    DAEIndex = 3; // default index
    assert(dtPlot >= dt);

    integratorSetups(system_);
    if (FlagCoutInfo) {
      cout << "DAE index " << DAEIndex;
      if (useExternalJac)
        cout << " external analytical Jacobian." << endl;
      else
        cout << " internally numerical computed Jacobian." << endl;
    }

    cout.setf(ios::scientific, ios::floatfield);

    if (FlagPlotIntegrator) {
      integPlot.open((name + ".plt").c_str());
      integPlot << "#1 t [s]:" << endl;
      integPlot << "#2 dt [s]:" << endl;
      integPlot << "#3 calculation time [s]:" << endl;
    }

    stepPlot = (int) (dtPlot / dt + 0.5);
    if (fabs(stepPlot * dt - dtPlot) > dt * dt) {
      cout << "WARNING: Due to the plot-Step settings it is not possible to plot exactly at the correct times." << endl;
    }

    s0 = clock();
  }

  void GeneralizedAlphaIntegrator::subIntegrate(DynamicSystemSolver& system, double tEnd) {
    /* FIND EQUILIBRIUM*/
    alphaFun fun_hg(&system, this, gamma_stroke, beta_stroke);

    qla.resize(qSize + laBiSize);
    la.resize(laBiSize);
    Index qInd = Index(0, qSize - 1);
    Index laInd = Index(qSize, qla.size() - 1);

    qlaStart.resize(qla.size());

    /* use MultiDimNewtonMethod*/
//    MultiDimNewtonMethod newton(&fun_hg);
//    newton.setLinearAlgebra(1);
//    newton.setTolerance(gTol);

    /* use  MultiDimensionalNewtonMethod */
    NewtonJacobianFunction * jac = new NumericalNewtonJacobianFunction();
    MultiDimensionalNewtonMethod newton;
    map<Index, double> tolerances;
    tolerances.insert(pair<Index, double>(qInd, hTol));
    tolerances.insert(pair<Index, double>(laInd, gTol));
    LocalResidualCriteriaFunction cfunc(tolerances);
    GlobalResidualCriteriaFunction cfuncGlob(gTol);
    StandardDampingFunction dfunc(30, 0.2);
    newton.setFunction(&fun_hg);
    newton.setJacobianFunction(jac);
    newton.setCriteriaFunction(&cfunc);
    newton.setDampingFunction(&dfunc);
    newton.setJacobianUpdateFreq(updateJacobianEvery);

    static_cast<DynamicSystem&>(system).plot(t, dt);

    while (t < tEnd) { // time loop
      integrationSteps++;

      fun_hg.setT(t);

      AlphaStep(newton);

      system.setLa(la); // only for plotting?

      // todo: check whether the plot make the openMBV unstable.
      if ((step * stepPlot - integrationSteps) < 0) {
        /* WRITE OUTPUT */
        step++;
        static_cast<DynamicSystem&>(system).plot(t, dt);
        double s1 = clock();
        time += (s1 - s0) / CLOCKS_PER_SEC;
        s0 = s1;
        integPlot << t << " " << dt << " " << iter << " " << time << " " << system.getlaSize() << endl;
        if (output)
          cout << "   t = " << t << ",\tdt = " << dt << ",\titer = " << setw(5) << setiosflags(ios::left) << iter << "\r" << flush;
        tPlot += dtPlot;
      }

      /* UPDATE SYSTEM FOR NEXT STEP*/
      t += dt;   // step 0: update time, go into new time step.
//      x += system.deltax(z, t, dt);  // todo: framework is not ready for x.
    }
  }

  void GeneralizedAlphaIntegrator::AlphaStep(MultiDimensionalNewtonMethod& newton) {
    // preparing for the starting values: q, u, ud, la
    q = q + dt * u + dt * dt * (0.5 - beta) * a;
    u = u + dt * (1 - gamma) * a;
    la.init(0.0);
    a = 1 / (1 - alpha_m) * (alpha_f * ud - alpha_m * a);

    q = q + dt * dt * beta * a;
    u = u + dt * gamma * a;
    ud.init(0.0);

    // copy q, la to qlaStart
    qlaStart.set(Index(0, qSize - 1), q);
    qlaStart.set(Index(qSize, qlaStart.size() - 1), la);

    qla = newton.solve(qlaStart);  // u, ud are updated inside function operators.
    if (newton.getInfo() != 0)
      throw MBSimError("ERROR (GeneralizedAlphaIntegrator::subIntegrate): No convergence of Newton method for the new time step");
    iter = newton.getNumberOfIterations();

    for (int i = 0; i < q.size(); i++)
      q(i) = qla(i);
    for (int i = 0; i < la.size(); i++)
      la(i) = qla(q.size() + i);

    a = a + (1 - alpha_f) / (1 - alpha_m) * ud;

  }

  void GeneralizedAlphaIntegrator::postIntegrate(DynamicSystemSolver& system) {
    integPlot.close();

    typedef tee_device<ostream, ofstream> TeeDevice;
    typedef stream<TeeDevice> TeeStream;
    ofstream integSum((name + ".sum").c_str());
    TeeDevice ts_tee(cout, integSum);
    TeeStream ts_split(ts_tee);

    ts_split << endl << endl << "******************************" << endl;
    ts_split << "INTEGRATION SUMMARY: " << endl;
    ts_split << "End time [s]: " << tEnd << endl;
    ts_split << "Integration time [s]: " << time << endl;
    ts_split << "Integration steps: " << integrationSteps << endl;
    ts_split << "Maximum number of iterations: " << maxIter << endl;
    ts_split << "Average number of iterations: " << double(sumIter) / integrationSteps << endl;
    ts_split << "******************************" << endl;
    ts_split.flush();
    ts_split.close();

    cout.unsetf(ios::scientific);
    cout << endl;
  }

  fmatvec::Vec alphaFun::operator()(const fmatvec::Vec& qla) {
    int qSize = sys->getqSize();
    int laSize = sys->getlaSize();
    int qlaSize = qSize + laSize;

    Index Iq(0, qSize - 1);
    Index Ila(qSize, qlaSize - 1);

    Vec delta_q = qla(Iq) - sys->getq();

    sys->setu(sys->getu() + gamma_stroke * delta_q);  // u is needed for updateh();
    inte->setUd(inte->getUd() + beta_stroke * delta_q);  // ud is directly used in hg function

    // set the q of system to be the input q
    sys->setq(qla(0, qSize - 1));
    sys->setLa(qla(qSize, qlaSize - 1));

    sys->updateStateDependentVariables(t);
    sys->updateg(t); // for joint, gap distance need also be updated.
//    sys->checkActive(1);   // todo: flag = 1, gap distance level
    sys->updategd(t);  //todo:  gd is needed for updating h,
//    sys->updateT(t);
    sys->updateJacobians(t); // needed for calculating W.
    sys->updateh(t);
    sys->updateM(t);
    if ((qlaSize - qSize) > 0)
      sys->updateW(t);
//    sys->updateV(t);  //not needed, because updateV is only needed for contact calculation
//    sys->updateG(t);  //  todo: not needed, because la is not solved by contact iteration

    // get the new h vector
    Vec hg;
    hg.resize(qla.size());

    // need new M, ud, h, W, g, input la (not from the system, but from the input value),
    hg(0, qSize - 1) = sys->getM() * inte->getUd() - sys->geth() - sys->getW() * qla(qSize, qlaSize - 1);
    hg(qSize, qlaSize - 1) = sys->getg().copy();

//    cout << "t = "  << t << "\n";
//    cout << "sys.geth() = "  <<  sys->geth().T() << "\n";
//    cout << "la= " << sys->getla().T() << endl;
//    cout << "W= " << sys->getW().T() << endl;
//    cout << "w * la = "  << Vec(sys->getW() * sys->getla()).T()  << "\n \n";
//    cout << "hg = "  << hg.T()  << "\n \n";

    return hg;
  }

}

