/*
 * IntegratorDAE.h
 *
 *  Created on: Mar 30, 2015
 *      Author: zwang
 */

#ifndef KERNEL_MBSIM_INTEGRATORS_INTEGRATORDAE_H_
#define KERNEL_MBSIM_INTEGRATORS_INTEGRATORDAE_H_

#include <fmatvec/fmatvec.h>
#include "integrator.h"

using namespace fmatvec;

namespace MBSimIntegrator {
  
  class IntegratorDAE : public Integrator {
    public:
      IntegratorDAE();
      virtual ~IntegratorDAE();
      void setAbsoluteTolerance(const fmatvec::Vec &aTol_) {aTol.resize() = aTol_;}
      void setAbsoluteTolerance(double aTol_) {aTol.resize() = Vec(1, INIT, aTol_);}
      void setRelativeTolerance(const fmatvec::Vec &rTol_) {rTol.resize() = rTol_;}
      void setRelativeTolerance(double rTol_) {rTol.resize() = Vec(1, INIT, rTol_);}
      void setInitialStepSize(double dt0_) {dt0 = dt0_;}
      void setMaximalStepSize(double dtMax_) {dtMax = dtMax_;}
      void setMaxStepNumber(int maxSteps_) {maxSteps = maxSteps_;}
      void useExternalJacobian(bool flag = true) {useExternalJac = flag;}

      virtual void initializeUsingXML(xercesc::DOMElement *element);

      virtual void integratorSetups(MBSim::DynamicSystemSolver& system) = 0;
      virtual void integrate(MBSim::DynamicSystemSolver& system);

      double getH() const {
        return H;
      }
      
      void setH(double h) {
        H = h;
      }
      
      double getRho() const {
        return rho;
      }
      
      void setRho(double rho) {
        this->rho = rho;
      }
      
      bool isScaling() const {
        return scaling;
      }
      
      void setScaling(bool scaling) {
        this->scaling = scaling;
      }
      
      int getDAEIndex() const {
        return DAEIndex;
      }
      
      void setDAEIndex(int daeIndex) {
        DAEIndex = daeIndex;
      }

    protected:
      /**
       * \brief Absolute Tolerance
       */
      fmatvec::Vec aTol;

      /**
       * \brief Relative Tolerance
       */
      fmatvec::Vec rTol;

      /**
       * \brief Flag of scalar or vector error setting
       */
      int iTol;

      /**
       * \brief Step size for the first step
       */
      double dt0;

      /**
       * \brief Maximum number of steps
       */
      int maxSteps;

      /**
       * \brief Maximal step size
       */
      double dtMax;

      /**
       * \brief Index of DAE integrator
       */
      int DAEIndex;

      /**
       * \brief include or exclude algebraic variables or scale variables with stepsize
       */
      int FlagErrorTest;

      /**
       * \brief Use external Jacobian (true) or numerical calculation by integrator (false=default)
       */
      bool useExternalJac;

      /**
       * \brief Time (current value of the independent variable t)
       */
      double t;

      /**
       * \brief Plot time point
       */
       double tPlot;

      /**
       * \brief Real work array
       */
      Vec rWork;

      /**
       * \brief Length of the real work array
       */
      int lrWork;

      /**
       * \brief Integer work array
       */
      Vector<Var, int> iWork;

      /**
       * \brief Length of the integer work array
       */
      int liWork;

      /**
       * \brief Real parameter array, which can be used for communication  between the calling program and subroutines
       */
      Vec rPar;

      /**
       * \brief Integer parameter array, which can be used for communication  between the calling program and subroutines
       */
      Vector<Var, int> iPar;

      /**
       * Flag of outputMode of external integrator during the integration
       */
      int iout;

      /**
       * \brief Error code of the external integrator
       */
      int idid;

      /**
       * \brief Computational time for integration
       */
      double s0, time;

      /**
       * \brief Filestream for integrator info at each step
       */
      std::ofstream integPlot;

      /**
       * \brief Flag: filestream for integrator info at each step
       */
      bool FlagPlotIntegrator;

      /**
       * \brief Flag: write integration summary into a file
       */
      bool FlagPlotIntegrationSum;

      /**
       * \brief Flag: write integration info (Index, status of calc. init. values, integ. summary) to cout
       */
      bool FlagCoutInfo;

      /**
       * \brief System state vectors
       */
      Vec z, Y;

      /**
       * \brief System state vector sizes
       */
      int qSize, uSize, xSize, laBiSize, zSize, YSize;

      /**
       * \brief Flag of whether scaling the equation;
       */
      bool scaling;

      /**
       * \brief Initial time step size when scaling is not applied;
       */
      double H;

      /**
       * \brief rho * s = penalty factor;
       */
      double rho;

  };

} /* namespace MBSimIntegrator */

#endif /* KERNEL_MBSIM_INTEGRATORS_INTEGRATORDAE_H_ */
