/* Copyright (C) 2004-2010 MBSim Development Team
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

#ifndef _GENERALIZED_ALPHA_INTEGRATOR_H_
#define _GENERALIZED_ALPHA_INTEGRATOR_H_

#include "integratorDAE.h"
#include "mbsim/functions/function.h"
#include<mbsim/dynamic_system_solver.h>
#include<mbsim/integrators/integrators.h>
#include <mbsim/numerics/nonlinear_algebra/multi_dimensional_newton_method.h>

namespace MBSimIntegrator {

  /** 
   * brief quasi static time-stepping integrator
   * \author Zhan Wang
   * \date 2015-05-05 initial commit (Zhan Wang)
   */
  class GeneralizedAlphaIntegrator : public IntegratorDAE {
    public:
      /**
       * \brief constructor
       */
      GeneralizedAlphaIntegrator();

      /**
       * \brief destructor
       */
      virtual ~GeneralizedAlphaIntegrator() {
      }

      void integratorSetups(MBSim::DynamicSystemSolver& system_);
      void preIntegrate(MBSim::DynamicSystemSolver& system);
      void subIntegrate(MBSim::DynamicSystemSolver& system, double tEnd);
      void postIntegrate(MBSim::DynamicSystemSolver& system);
      void AlphaStep(MBSim::MultiDimensionalNewtonMethod& newton);


      /* GETTER / SETTER */
      void setStepSize(double dt_) {
        dt = dt_;
      }

      void setgTolerance(double tolerance_) {
        gTol = tolerance_;
      }

      void sethTolerance(double tolerance_) {
        hTol = tolerance_;
      }

      void setupdateJacobianEvery(int value) {
        updateJacobianEvery = value;
      }
      
      const fmatvec::Vec& getU() const {
        return u;
      }
      
      void setU(const fmatvec::Vec& u) {
        this->u = u;
      }
      
      const fmatvec::Vec& getUd() const {
        return ud;
      }
      
      void setUd(const fmatvec::Vec& ud) {
        this->ud = ud;
      }

      double getRhoInf() const {
        return rho_inf;
      }

      void setRhoInf(double rhoInf) {
        rho_inf = rhoInf;
      }

      /***************************************************/

    private:
      /**
       * \brief step size
       */
      double dt;

      /**
       * \brief time and plot time
       */
      double t, tPlot;

      /*!
       * \brief tolerance for the newton iteration for distances
       */
      double gTol;

      /*!
       * \brief tolerance for newton iteration for forces
       */
      double hTol;

      /**
       * \brief iteration counter for constraints, plots, integration, maximum constraints, cummulation constraint
       */
      int iter, step, integrationSteps, maxIter, sumIter;

      /*!
       * \brief value of how often the Jacobian should be updated every step
       */
      int updateJacobianEvery;

      /**
       * \brief computing time counter
       */
      double s0, time;

      /**
       * \brief plot step difference
       */
      int stepPlot;

      /**
       * \brief state, position, velocity, order coordinate of dynamical system
       */
      fmatvec::Vec z, q, u, x, la;

      /**
       * \brief acceleration vector
       */
      fmatvec::Vec ud;

      /**
       * \brief solution vector and starting vector of the Newton method
       */

      VecV qla, qlaStart;

      /**
       * \brief artificial acceleration
       */
      fmatvec::Vec a;

      /**
       * \brief file stream for integration information
       */
      std::ofstream integPlot;

      /**
       * \brief alpla parameters
       */
      double alpha_m, alpha_f;

      /**
       * \brief
       */
      double gamma, beta;

      /**
       * \brief
       */
      double gamma_stroke, beta_stroke;

      /**
       * \brief user-specified value of the spectral radius in the high-frequency limit
       */
      double rho_inf;
  };

  /*!
    * \brief calculate h vector according the new q and system boundary conditions
    * \author Zhan Wang
    * \date 2015-05-05 initial commit (Zhan Wang)
    */
   class alphaFun : public MBSim::Function<fmatvec::Vec(fmatvec::Vec)> {
     public:
       /*!
        * \brief constructor
        */
       alphaFun(MBSim::DynamicSystemSolver* sys_, MBSimIntegrator::GeneralizedAlphaIntegrator* inte_, double gamma_stroke_, double beta_stroke_) :
           sys(sys_), inte(inte_),t(0), gamma_stroke(gamma_stroke_), beta_stroke(beta_stroke_) {
       }

       /*!
        * \brief destructor
        */
       virtual ~alphaFun() {
       }
       ;

       /* INHERITED INTERFACE */
       fmatvec::Vec operator()(const fmatvec::Vec& qla);

       void setT(double t) {
         this->t = t;
       }

       /*******************************************************/

     private:

       MBSim::DynamicSystemSolver* sys;
       MBSimIntegrator::GeneralizedAlphaIntegrator* inte;
       double t;
       double gamma_stroke, beta_stroke;

   };
}

#endif /* _QUASI_STATIC_INTEGRATOR_H_ */

