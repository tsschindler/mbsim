/* Copyright (C) 2007 Zhan Wang
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
 * Contact:
 *   wangzhanrock@gmail.com
 *
 */

#ifndef _DASKR_INTEGRATOR_H_
#define _DASKR_INTEGRATOR_H_

using namespace fmatvec;

#include "integratorDAE.h"
#include <time.h>

namespace MBSimIntegrator {

  /** \brief DAE-Integrator DASKR (index 1) with root function
   *   based on original the DAE-Integrator DASKR (up to index 2) implemented by rhuber.
   *  \author Zhan Wang
   */
  class DASKRIntegrator : public IntegratorDAE {

      friend class DAETSIntegrator;

    private:
      static void residuum(double* t, double* Y_, double* Ydot_, double* CJ, double* res_, int* ires, double* rpar, int* ipar);
      static void jac(double* t, double* Y_, double* Ydot_, double* PD, double* CJ, double* rpar, int* ipar);
      static void rt(int* neq, double* t, double* Y_, double* Ydot_, int* nrt, double *rval, double* rpar, int* ipar);
      /**
       * \brief Flag: use Krylov Method to solve the linear system inside DASKR (false=default)
       */
      bool KrylovMethod;

      /**
       * \brief Derivatives of state variables at time t
       */
      Vec Ydot;

      /**
       * \brief Maximum method order (one to five)
       */
      int maxOrder;

      /**
       * \brief Integer array to communicate with integrator
       */
      Vector<Var, int> info;

      /**
       * \brief Integer array for output to indicate where one or more roots were found
       */
      Vector<Var, int> jroot;

      /**
       * \brief Number of root functions
       */
      int nrt;

      /**
       * \brief Integration statistics variables
       */
      int IntSteps, RefusedSteps, FuncEvals, JacEval, RootEvals;

      /**
       * \brief Work space required by JAC and/or PSOL, and space for data to
       *        be communicated from JAC to PSOL is made available in the form
       *        of arrays WP and IWP, which are parts of the RWORK and IWORK
       *        arrays, respectively.  The lengths of these real and integer
       *        work spaces WP and IWP must be supplied in LENWP and LENIWP,
       *        respectively, as follows..
       *        IWORK(27) = LENWP = length of real work space WP
       *        IWORK(28) = LENIWP = length of integer work space IWP
       *        LENWP and LENIWP should be provided by the user when Krylov
       *        method is used.
       */
      int LENWP, LENIWP;

    public:

      DASKRIntegrator();
      ~DASKRIntegrator() {
      }

      void setLENWP(int LENWP_) {
        LENWP = LENWP_;
      }
      void setLENIWP(int LENIWP_) {
        LENIWP = LENIWP_;
      }
      void setMaxOrder(int order_);

//      void setGGLStabilisation() {DAEIndex=21;};
      void integrate(MBSim::DynamicSystemSolver& system);
      int computeInitialConditions(MBSim::DynamicSystemSolver& system, bool FlagPlot, bool FlagCoutInfo_);
      void integratorSetups(MBSim::DynamicSystemSolver& system);
      void setFlagErrorTest(int Flag);

      /** subroutines used by integrate */
      void preIntegrate(MBSim::DynamicSystemSolver& system);
      void subIntegrate(MBSim::DynamicSystemSolver& system, double tEnd);
      void postIntegrate(MBSim::DynamicSystemSolver& system);


  };

}

#endif
