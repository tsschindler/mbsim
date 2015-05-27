/* Copyright (C) 2004-2014 MBSim Development Team
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

#ifndef _DYNAMIC_SYSTEM_H_
#define _DYNAMIC_SYSTEM_H_

#include "mbsim/element.h"
#include "mbsim/mbsim_event.h"

namespace H5 {
  class Group;
}

namespace MBSim {
  class Frame;
  class FixedRelativeFrame;
  class Contour;
  class Object;
  class Link;
  class ModellingInterface;
  class Contact;
  class InverseKineticsJoint;
  class Observer;

  /**
   * \brief dynamic system as topmost hierarchical level
   * \author Martin Foerg
   * \date 2009-03-26 some comments (Thorsten Schindler)
   * \date 2009-04-06 ExtraDynamicInterface included (Thorsten Schindler)
   * \date 2009-06-14 OpenMP (Thorsten Schindler)
   * \date 2009-07-08 relative dynamic system location (Thorsten Schindler)
   * \date 2009-07-16 splitted link / object right hand side (Thorsten Schindler)
   * \date 2009-07-27 implicit integration improvement (Thorsten Schindler)
   * \date 2009-07-28 splitted interfaces (Thorsten Schindler)
   * \date 2009-12-14 revised inverse kinetics (Martin Foerg)
   * \date 2010-07-06 modifications for timestepper ssc, e.g LinkStatus and buildListOfSetValuedLinks (Robert Huber)
   * \date 2012-05-08 OpenMP completely removed; will be inserted again soon (Jan Clauberg)
   * \date 2012-05-08 modifications for AutoTimeSteppingSSCIntegrator (Jan Clauberg)
   * \date 2013-06-23 removed extra dynamics (Martin Foerg)
   * \date 2014-09-16 contact forces are calculated on acceleration level (Thorsten Schindler)
   */
  class DynamicSystem : public Element {
    public:
      /** 
       * \brief constructor
       */
      DynamicSystem(const std::string &name);

      /**
       * \brief destructor
       */
      virtual ~DynamicSystem();

      /* INTERFACE FOR DERIVED CLASSES */
      virtual void updateT(double t); 
      virtual void updateh(double t, int i=0); 
      virtual void updateStateDependentVariables(double t) = 0;
      virtual void updateStateDerivativeDependentVariables(double t);
      virtual void updateM(double t, int i=0);
      virtual void updateJacobians(double t, int j=0) = 0;
      virtual void updatedq(double t, double dt); 
      virtual void updateud(double t, int i=0) { THROW_MBSIMERROR("(DynamicSystem::updateud): Not implemented!"); }
      virtual void updatezd(double t) = 0;
      virtual void updatedu(double t, double dt) = 0;
      virtual void updateqd(double t);
      virtual void sethSize(int hSize_, int i=0);
      virtual int gethSize(int i=0) const { return hSize[i]; }
      virtual int getqSize() const { return qSize; }
      virtual int getuSize(int i=0) const { return uSize[i]; }
      virtual void calcqSize();
      virtual void calcuSize(int j=0);
      //virtual int getqInd(DynamicSystem* sys);
      virtual int getuInd(int i=0) { return uInd[i]; }
      //virtual int getuInd(DynamicSystem* sys, int i=0);
      //virtual void setqInd(int qInd_) { qInd = qInd_; }
      //virtual void setuInd(int uInd_, int i=0) { uInd[i] = uInd_; }
      virtual void setqInd(int qInd_);
      virtual void setuInd(int uInd_, int i=0);
      virtual void sethInd(int hInd_, int i=0);
      virtual void setxInd(int xInd_);
      //virtual int gethInd(DynamicSystem* sys, int i=0); 
      virtual const fmatvec::Vec& getq() const { return q; };
      virtual fmatvec::Vec& getq() { return q; };
      virtual const fmatvec::Vec& getu() const { return u; };
      virtual fmatvec::Vec& getu() { return u; };
      virtual H5::GroupBase *getPlotGroup() { return plotGroup; }
      virtual PlotFeatureStatus getPlotFeature(PlotFeature fp) { return Element::getPlotFeature(fp); };
      virtual PlotFeatureStatus getPlotFeatureForChildren(PlotFeature fp) { return Element::getPlotFeatureForChildren(fp); };
#ifdef HAVE_OPENMBVCPPINTERFACE
      virtual boost::shared_ptr<OpenMBV::Group> getOpenMBVGrp();
#endif

      virtual void updatewb(double t, int j=0); 
      virtual void updateW(double t, int j=0); 
      virtual void updateV(double t, int j=0); 
      virtual void updateg(double t);
      virtual void updategd(double t);
      virtual void updateStopVector(double t); 
      virtual void updateLinkStatus(double t);
      virtual void updateLinkStatusReg(double t);

      virtual void updategInverseKinetics(double t); 
      virtual void updategdInverseKinetics(double t);
      virtual void updateWInverseKinetics(double t, int j=0); 
      virtual void updatehInverseKinetics(double t, int j=0); 
      virtual void updateJacobiansInverseKinetics(double t, int j=0); 
      virtual void updatebInverseKinetics(double t); 

      virtual void updatedx(double t, double dt); 
      virtual void updatexd(double t);
      virtual void calcxSize();
      const fmatvec::Vec& getx() const { return x; };
      fmatvec::Vec& getx() { return x; };
      int getxSize() const { return xSize; }
      void updatexRef(const fmatvec::Vec &ref);
      void updatexdRef(const fmatvec::Vec &ref);
      virtual void init(InitStage stage);
      virtual void initz();
      virtual void writez(H5::GroupBase *group);
      virtual void readz0(H5::GroupBase *parent);
      /*****************************************************/

      /* INHERITED INTERFACE OF ELEMENT */
      /** DEPRECATED */
      virtual std::string getType() const { return "DynamicSystem"; }
      virtual void setDynamicSystemSolver(DynamicSystemSolver* sys);
      virtual void plot(double t, double dt);
      virtual void plotAtSpecialEvent(double t, double dt=1.);
      virtual void closePlot();
      /*****************************************************/

      /* INTERFACE FOR DERIVED CLASSES */
      /**
       * \brief compute Cholesky decomposition of mass matrix TODO necessary?
       */
      virtual void facLLM(int i=0) = 0;

      /**
       * \brief solve contact equations with single step fixed point scheme
       * \return iterations of solver
       */
      virtual int solveConstraintsFixpointSingle();

      /**
       * \brief solve impact equations with single step fixed point scheme on velocity level 
       * \param time step-size
       * \return iterations of solver
       */
      virtual int solveImpactsFixpointSingle(double dt);

      /**
       * \brief solve contact equations with Gauss-Seidel scheme
       * \return iterations of solver
       */
      virtual int solveConstraintsGaussSeidel();

      /**
       * \brief solve impact equations with Gauss-Seidel scheme on velocity level 
       * \param time step-size
       * \return iterations of solver
       */
      virtual int solveImpactsGaussSeidel(double dt);

      /**
       * \brief solve contact equations with Newton scheme
       * \return iterations of solver
       */
      virtual int solveConstraintsRootFinding();

      /**
       * \brief solve impact equations with Newton scheme on velocity level 
       * \param time step-size
       * \return iterations of solver
       */
      virtual int solveImpactsRootFinding(double dt);

      /**
       * \brief compute JACOBIAN of contact equations
       */
      virtual int jacobianConstraints();

      /**
       * \brief compute JACOBIAN of contact equations on velocity level
       */
      virtual int jacobianImpacts();

      /**
       * \brief validate force laws concerning given tolerances
       */
      virtual void checkConstraintsForTermination();

      /**
       * \brief validate force laws concerning given tolerances on velocity level
       */
      virtual void checkImpactsForTermination(double dt);

      /**
       * \brief update relaxation factors for contact equations
       */
      virtual void updaterFactors();

      /**
       * \param name of the frame
       * \param check for existence of frame
       * \return frame
       */
      virtual Frame* getFrame(const std::string &name, bool check=true) const;

      /**
       * \param name of the contour
       * \param check for existence of contour
       * \return contour
       */
      virtual Contour* getContour(const std::string &name, bool check=true) const;
      /*****************************************************/

      /* GETTER / SETTER */
      void setPosition(const fmatvec::Vec3& PrPF_) { PrPF = PrPF_; }
      void setOrientation(const fmatvec::SqrMat3& APF_) { APF = APF_; }
      void setFrameOfReference(Frame *frame) { R = frame; };
      const fmatvec::Vec3& getPosition() const { return PrPF; }
      const fmatvec::SqrMat3& getOrientation() const { return APF; }
      const Frame* getFrameOfReference() const { return R; };

      const fmatvec::Vec& getqd() const {return qd;}
      const fmatvec::Vec& getud(int i = 0) const {return ud[i];}
      const fmatvec::Vec& getxd() const { return xd; };
      fmatvec::Vec& getxd() { return xd; };
      const fmatvec::Vec& getx0() const { return x0; };
      fmatvec::Vec& getx0() { return x0; };

      const fmatvec::Mat& getT() const { return T; };
      const fmatvec::SymMat& getM(int i=0) const { return M[i]; };
      const fmatvec::SymMat& getLLM(int i=0) const { return LLM[i]; };
      fmatvec::SymMat& getLLM(int i=0) { return LLM[i]; };
      const fmatvec::Vec& geth(int j=0) const { return h[j]; };
      fmatvec::Vec& geth(int j=0) { return h[j]; };
      const fmatvec::Vec& getf() const { return f; };
      fmatvec::Vec& getf() { return f; };

      const fmatvec::Mat& getW(int i=0) const { return W[i]; }
      fmatvec::Mat& getW(int i=0) { return W[i]; }
      const fmatvec::Vec& getwb() const { return wb; }
      fmatvec::Vec& getwb() { return wb; }
      const fmatvec::Mat& getV(int i=0) const { return V[i]; }
      fmatvec::Mat& getV(int i=0) { return V[i]; }

      const fmatvec::Vec& getla() const { return la; }
      fmatvec::Vec& getla() { return la; }
      const fmatvec::Vec& getg() const { return g; }
      fmatvec::Vec& getg() { return g; }
      const fmatvec::Vec& getgd() const { return gd; }
      fmatvec::Vec& getgd() { return gd; }
      const fmatvec::Vec& getgdd() const { return gdd; }
      fmatvec::Vec& getgdd() { return gdd; }
      const fmatvec::Vec& getrFactor() const { return rFactor; }
      fmatvec::Vec& getrFactor() { return rFactor; }
      fmatvec::Vec& getsv() { return sv; }
      const fmatvec::Vec& getsv() const { return sv; }
      fmatvec::VecInt& getjsv() { return jsv; }
      const fmatvec::VecInt& getjsv() const { return jsv; }
      fmatvec::VecInt& getLinkStatus() { return LinkStatus; }
      fmatvec::VecInt& getLinkStatusReg() { return LinkStatusReg; }
      const fmatvec::VecInt& getLinkStatus() const { return LinkStatus; }
      const fmatvec::VecInt& getLinkStatusReg() const { return LinkStatusReg; }
      const fmatvec::Vec& getres() const { return res; }
      fmatvec::Vec& getres() { return res; }
      const fmatvec::Vec& getcorr() const { return corr; };
      fmatvec::Vec& getcorr() { return corr; };

      void setx(const fmatvec::Vec& x_) { x = x_; }
      void setx0(const fmatvec::Vec &x0_) { x0 = x0_; }
      void setx0(double x0_) { x0 = fmatvec::Vec(1,fmatvec::INIT,x0_); }

      int getxInd() { return xInd; }
      int getlaInd() const { return laInd; } 

      int gethInd(int i=0) { return hInd[i]; }
      void setlaInd(int ind) { laInd = ind; }
      void setgInd(int ind) { gInd = ind; }
      void setgdInd(int ind) { gdInd = ind; }
      void setrFactorInd(int ind) { rFactorInd = ind; }
      virtual void setsvInd(int svInd_);
      void setLinkStatusInd(int LinkStatusInd_) {LinkStatusInd = LinkStatusInd_;};      
      void setLinkStatusRegInd(int LinkStatusRegInd_) {LinkStatusRegInd = LinkStatusRegInd_;};

      int getzSize() const { return qSize + uSize[0] + xSize; }

      void setqSize(int qSize_) { qSize = qSize_; }
      void setuSize(int uSize_, int i=0) { uSize[i] = uSize_; }
      void setxSize(int xSize_) { xSize = xSize_; }

      int getlaSize() const { return laSize; } 
      int getgSize() const { return gSize; } 
      int getgdSize() const { return gdSize; } 
      int getrFactorSize() const { return rFactorSize; } 
      int getsvSize() const { return svSize; }
      int getLinkStatusSize() const { return LinkStatusSize; }
      int getLinkStatusRegSize() const { return LinkStatusRegSize; }
      /*****************************************************/

      const std::vector<Object*>& getObjects() const { return object; }
      const std::vector<Link*>& getLinks() const { return link; }
      const std::vector<DynamicSystem*>& getDynamicSystems() const { return dynamicsystem; }
      const std::vector<Frame*>& getFrames() const { return frame; }
      const std::vector<Contour*>& getContours() const { return contour; }

      /**
       * \brief references to positions of dynamic system parent
       * \param vector to be referenced
       */
      void updateqRef(const fmatvec::Vec &ref); 

      /**
       * \brief references to differentiated positions of dynamic system parent
       * \param vector to be referenced
       */
      void updateqdRef(const fmatvec::Vec &ref);

      /**
       * \brief references to velocities of dynamic system parent
       * \param vector to be referenced
       */
      void updateuRef(const fmatvec::Vec &ref);

      /**
       * \brief references to velocities of dynamic system parent
       * \param vector to be referenced
       */
      void updateuallRef(const fmatvec::Vec &ref);

      /**
       * \brief references to differentiated velocities of dynamic system parent
       * \param vector to be referenced
       */
      void updateudRef(const fmatvec::Vec &ref, int i=0);

      /**
       * \brief references to velocities of dynamic system parent
       * \param vector to be referenced
       */
      void updateudallRef(const fmatvec::Vec &ref, int i=0);

      /**
       * \brief references to smooth right hand side of dynamic system parent
       * \param complete vector to be referenced
       * \param vector concerning objects to be referenced
       * \param vector concerning links to be referenced
       * \param index of normal usage and inverse kinetics
       */
      void updatehRef(const fmatvec::Vec &hRef, int i=0);

      /**
       * \brief references to nonsmooth right hand side of dynamic system parent
       * \param vector to be referenced
       */
      void updaterRef(const fmatvec::Vec &ref, int j=0);

      /**
       * \brief references to linear transformation matrix between differentiated positions and velocities of dynamic system parent
       * \param matrix to be referenced
       */
      void updateTRef(const fmatvec::Mat &ref);

      /**
       * \brief references to mass matrix of dynamic system parent
       * \param matrix to be referenced
       * \param index of normal usage and inverse kinetics
       */
      void updateMRef(const fmatvec::SymMat &ref, int i=0);

      /**
       * \brief references to Cholesky decomposition of dynamic system parent
       * \param matrix to be referenced
       * \param index of normal usage and inverse kinetics
       */
      void updateLLMRef(const fmatvec::SymMat &ref, int i=0);

      /**
       * \brief references to relative distances of dynamic system parent
       * \param vector to be referenced
       */
      void updategRef(const fmatvec::Vec &ref);

      /**
       * \brief references to relative velocities of dynamic system parent
       * \param vector to be referenced
       */
      void updategdRef(const fmatvec::Vec &ref);

      /**
       * \brief references to contact force parameters of dynamic system parent
       * \param vector to be referenced
       */
      void updatelaRef(const fmatvec::Vec &ref);

      void updatelaInverseKineticsRef(const fmatvec::Vec &ref);
      void updatebInverseKineticsRef(const fmatvec::Mat &ref);

      /**
       * \brief references to TODO of dynamic system parent
       * \param vector to be referenced
       */      
      void updatewbRef(const fmatvec::Vec &ref);

      /**
       * \brief references to contact force direction matrix of dynamic system parent
       * \param matrix to be referenced
       * \param index of normal usage and inverse kinetics
       */
      void updateWRef(const fmatvec::Mat &ref, int i=0);

      /**
       * \brief references to contact force direction matrix of dynamic system parent
       * \param matrix to be referenced
       * \param index of normal usage and inverse kinetics
       */
      void updateWInverseKineticsRef(const fmatvec::Mat &ref, int i=0);

      /**
       * \brief references to condensed contact force direction matrix of dynamic system parent
       * \param matrix to be referenced
       * \param index of normal usage and inverse kinetics
       */
      void updateVRef(const fmatvec::Mat &ref, int i=0);

      /**
       * \brief references to stopvector (rootfunction for event driven integrator) of dynamic system parent
       * \param vector to be referenced
       */
      void updatesvRef(const fmatvec::Vec& ref);

      /**
       * \brief references to boolean evaluation of stopvector concerning roots of dynamic system parent
       * \param vector to be referenced
       */
      void updatejsvRef(const fmatvec::VecInt &ref);

      /**
       * \brief references to status vector of set valued links with piecewise link equations (which piece is valid)
       * \param vector to be referenced 
       */
      void updateLinkStatusRef(const fmatvec::VecInt &LinkStatusParent);

      /**
       * \brief references to status vector of single valued links
       * \param vector to be referenced 
       */
      void updateLinkStatusRegRef(const fmatvec::VecInt &LinkStatusRegParent);      
            
      /**
       * \brief references to residuum of contact equations of dynamic system parent
       * \param vector to be referenced
       */
      void updateresRef(const fmatvec::Vec &ref);

      /**
       * \brief references to relaxation factors for contact equations of dynamic system parent
       * \param vector to be referenced
       */
      void updaterFactorRef(const fmatvec::Vec &ref);

      void clearElementLists();

      /**
       * \brief build flat list of dynamic systems
       * \param list of dynamic systems
       */
      void buildListOfDynamicSystems(std::vector<DynamicSystem*> &sys);

      /**
       * \brief build flat list of objects
       * \param list of objects
       */
      void buildListOfObjects(std::vector<Object*> &obj);

      /**
       * \brief build flat list of links
       * \param list of links
       */
      void buildListOfLinks(std::vector<Link*> &lnk);

      /**
       * \brief build flat list of all setvalued links
       * \param list of links
       */
      void buildListOfSetValuedLinks(std::vector<Link*> &lnk);

      /**
       * \brief build flat list of frames
       * \param list of frames
       */
      void buildListOfFrames(std::vector<Frame*> &frm);

      /**
       * \brief build flat list of contours
       * \param list of contours
       */
      void buildListOfContours(std::vector<Contour*> &cnt);

      /**
       * \brief build flat list of models
       * \param list of models
       */
      void buildListOfModels(std::vector<ModellingInterface*> &model);

      /**
       * \brief build flat list of inverse kinetics links
       * \param list of inverse kinetics links
       */
      void buildListOfInverseKineticsLinks(std::vector<Link*> &lnk);

      /**
       * \brief build flat list of observers
       * \param list of observers
       */
      void buildListOfObservers(std::vector<Observer*> &obsrv);

      /**
       * \brief analyse constraints of dynamic systems for usage in inverse kinetics
       */
      void setUpInverseKinetics();

      /**
       * \brief distribute links to set- and single valued container
       */
      void setUpLinks();

      /**
       * \return flag, if vector of active relative distances has changed (set-valued laws)
       */
      bool gActiveChanged();
      
      /**
       * \return flag, if vector of active relative distances has changed (single-valued laws)
       */
      bool gActiveChangedReg();
      
      /**
       * \return flag, if an impact occured in the system
       */
      bool detectImpact();

      /**
       * \brief calculates size of stop vector
       */
      void calcsvSize();

      /**
       * \brief calculates size of contact force parameters
       */
      void calclaSize(int j);

      /**
       * \brief calculates size of set-valued link status vector
       */
      void calcLinkStatusSize();

      /**
       * \brief calculates size of single-valued link status vector
       */
      void calcLinkStatusRegSize();

      /**
       * \brief calculates size of contact force parameters
       */
      void calclaInverseKineticsSize();

      /**
       * \brief calculates size of contact force parameters
       */
      void calcbInverseKineticsSize();

      /**
       * \brief calculates size of relative distances
       */
      void calcgSize(int j);

      /**
       * \brief calculates size of gap velocities
       * \param flag to decide which contacts are included in the calculation
       * 0 = all contacts
       * 1 = closed contacts
       * 2 = contacts which stay closed
       * 3 = sticking contacts
       *
       * see SingleContact for the implementation
       */
      void calcgdSize(int j);

      /**
       * \brief calculates size of relaxation factors for contact equations
       */
      void calcrFactorSize(int j);

      /** 
       * \brief rearrange vector of active setvalued links
       */
      void setUpActiveLinks();

      /**
       * \brief check if set-valued contacts are active and set corresponding attributes
       * \param flag to decide which criteria are used to define 'activity'
       * 1 = position level activity (gActive = normal gap not larger than tolerance)
       * 2 = velocity level activity (if gActive, gdActive[0] = normal gap velocity not larger than tolerance ; if gdActive[0], gdActive[1] = tangential gap velocity not larger than tolerance, i.e., sticking)
       * 3 = velocity level activity (if gActive, gdActive[0] = new normal gap velocity AFTER impact not larger than tolerance ; if gdActive[0], gdActive[1] = new tangential gap velocity AFTER impact not larger than tolerance, i.e., sticking)
       * 4 = acceleration level activity (if gActive and gdActive[0], gddActive[0] = normal gap acceleration not larger than tolerance ; if gddActive[0] and gdActive[1], gddActive[1] = tangential gap acceleration not larger than tolerance, i.e., stay sticking)
       * 5 = activity clean-up: if there is no activity on acceleration or velocity level, also more basic levels are set to non-active
       * 6 = closing sets all activities
       * 7 = slip-stick transition sets tangential velocity and acceleration activity
       * 8 = opening and stick-slip transition set corresponding acceleration activity to non-active
       *
       * see SingleContact for the implementation
       */
      void checkActive(int i);

      /**
       * \brief check if single-valued contacts are active
       */
      void checkActiveReg(int i);

      /**
       * \param tolerance for relative velocity
       */
      virtual void setgTol(double tol);

      /**
       * \param tolerance for relative velocity
       */
      virtual void setgdTol(double tol);

      /**
       * \param tolerance for relative acceleration
       */
      virtual void setgddTol(double tol);

      /**
       * \param tolerance for contact force
       */
      virtual void setlaTol(double tol);

      /**
       * \param tolerance for impact
       */
      virtual void setLaTol(double tol);

      /**
       * \param maximum relaxation factor for contact equations
       */
      void setrMax(double rMax);

      void addFrame(FixedRelativeFrame *frame);

      /**
       * \param frame
       * \return index of frame TODO renaming
       */
      int frameIndex(const Frame *frame_) const;

      /**
       * \param dynamic system to add
       */
      void addGroup(DynamicSystem *dynamicsystem);

      /**
       * \param name of the dynamic system
       * \param check for existence of dynamic system
       * \return dynamic system
       */
      DynamicSystem* getGroup(const std::string &name,bool check=true) const;

      /**
       * \param object to add
       */
      void addObject(Object *object);

      /**
       * \param name of the object
       * \param check for existence of object
       * \return object
       */
      Object* getObject(const std::string &name,bool check=true) const;

      /**
       * \param link to add
       */
      void addLink(Link *link);

      /**
       * \param add link for inverse kinetics
       */
      void addInverseKineticsLink(Link *link);

      Observer* getObserver(const std::string &name,bool check=true) const;
      void addObserver(Observer *element);

      /**
       * \param name of the link
       * \param check for existence of link
       * \return link
       */
      Link* getLink(const std::string &name,bool check=true) const;

      /**
       * \param modell to add
       */
      void addModel(ModellingInterface *modell);

      /**
       * \param name of the model
       * \param check for existence of model
       * \return modelling interface
       */
      ModellingInterface* getModel(const std::string &name, bool check=true) const;

      /** Return frame "I" */
      FixedRelativeFrame *getFrameI() { return I; }

      virtual Element *getChildByContainerAndName(const std::string &container, const std::string &name) const;

      virtual void updatecorr(int j);
      void updatecorrRef(const fmatvec::Vec &ref);
      void calccorrSize(int j);

      void checkRoot();

    protected:
      /**
       * \brief parent frame
       */
      Frame *R;

      /**
       * \brief relative translation with respect to parent frame
       */
      fmatvec::Vec3 PrPF;

      /**
       * \brief relative rotation with respect to parent frame
       */
      fmatvec::SqrMat3 APF;

      /** 
       * \brief container for possible ingredients
       */
      std::vector<Object*> object;
      std::vector<Link*> link;
      std::vector<Link*> linkSingleValued;
      std::vector<Link*> linkSetValued;
      std::vector<Link*> linkSetValuedActive;
      std::vector<ModellingInterface*> model;
      std::vector<DynamicSystem*> dynamicsystem;
      std::vector<Link*> inverseKineticsLink;
      std::vector<Observer*> observer;
      std::vector< std::vector<Link*> > linkOrdered;

      /** 
       * \brief linear relation matrix of position and velocity parameters
       */
      fmatvec::Mat T;

      /**
       * \brief mass matrix
       */
      fmatvec::SymMat M[2];

      /** 
       * \brief Cholesky decomposition of mass matrix
       */
      fmatvec::SymMat LLM[2];

      /**
       * \brief positions, differentiated positions, initial positions
       */
      fmatvec::Vec q, qd, q0;

      /**
       * \brief velocities, differentiated velocities, initial velocities
       */
      fmatvec::Vec u, ud[2], u0;

      /**
       * \brief order one parameters, differentiated order one parameters, initial order one parameters
       */
      fmatvec::Vec x, xd, x0;

      /**
       * \brief smooth, smooth with respect to objects, smooth with respect to links, nonsmooth and order one right hand side
       */
      fmatvec::Vec h[2], r[2], f;

      /**
       * \brief 
       */
      fmatvec::Mat W[2], V[2];

      /**
       * \brief contact force parameters
       */
      fmatvec::Vec la;

      /** 
       * \brief relative distances and velocities
       */
      fmatvec::Vec g, gd, gdd;

      /**
       * \brief TODO
       */
      fmatvec::Vec wb;

      /**
       * \brief residuum of nonlinear contact equations for Newton scheme
       */
      fmatvec::Vec res;

      /**
       * \brief rfactors for relaxation nonlinear contact equations
       */
      fmatvec::Vec rFactor;

      /**
       * \brief stop vector (root functions for event driven integration
       */
      fmatvec::Vec sv;

      /**
       * \brief boolean evaluation of stop vector concerning roots
       */
      fmatvec::VecInt jsv;

      /**
       * \brief status of set-valued links 
       */
      fmatvec::VecInt LinkStatus;

      /**
       * \brief status of single-valued links 
       */
      fmatvec::VecInt LinkStatusReg;

      /** 
       * \brief size and local start index of positions relative to parent
       */
      int qSize, qInd;

      /** 
       * \brief size and local start index of velocities relative to parent
       */
      int uSize[2], uInd[2];

      /** 
       * \brief size and local start index of order one parameters relative to parent
       */
      int xSize, xInd;

      /** 
       * \brief size and local start index of order smooth right hand side relative to parent
       */
      int hSize[2], hInd[2];

      /** 
       * \brief size and local start index of relative distances relative to parent
       */
      int gSize, gInd;

      /** 
       * \brief size and local start index of relative velocities relative to parent
       */
      int gdSize, gdInd;

      /** 
       * \brief size and local start index of contact force parameters relative to parent
       */
      int laSize, laInd;

      /** 
       * \brief size and local start index of rfactors relative to parent
       */
      int rFactorSize, rFactorInd;

      /** 
       * \brief size and local start index of stop vector relative to parent
       */
      int svSize, svInd;

      /**
       * \brief size and local start index of set-valued link status vector relative to parent
       */
      int LinkStatusSize, LinkStatusInd;
      
      /**
       * \brief size and local start index of single-valued link status vector relative to parent
       */
      int LinkStatusRegSize, LinkStatusRegInd;
      
      /**
       * \brief vector of frames and contours
       */
      std::vector<Frame*> frame;
      std::vector<Contour*> contour;

#ifdef HAVE_OPENMBVCPPINTERFACE
      boost::shared_ptr<OpenMBV::Group> openMBVGrp;
#endif
      boost::shared_ptr<H5::File> hdf5File;

      /**
       * \param contour to add
       */
      void addContour(Contour* contour);

      /** A pointer to frame "I" */
      FixedRelativeFrame *I;

      /** 
       * \brief size of contact force parameters of special links relative to parent
       */
      int laInverseKineticsSize, bInverseKineticsSize;

      fmatvec::Mat WInverseKinetics[2], bInverseKinetics;
      fmatvec::Vec laInverseKinetics;

      int corrSize, corrInd;
      fmatvec::Vec corr;

      std::string saved_frameOfReference;
  };
}

#endif

