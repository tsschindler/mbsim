/* Copyright (C) 2004-2012 MBSim Development Team
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

#ifndef POLYNOMIAL_FRUSTUM_H_
#define POLYNOMIAL_FRUSTUM_H_
#include <iostream>
#include <mbsim/contour.h>
#include <mbsim/contours/rectangle.h>
#include <mbsim/utils/colors.h>
#include <fmatvec/fmatvec.h>
#include <mbsim/functions/function.h>
#ifdef HAVE_OPENMBVCPPINTERFACE
#include <mbsim/utils/boost_parameters.h>
#endif

namespace MBSim {

  /*!
   * \brief Frustum contour with a polynomial radius over height course
   * \author Kilian Grundl, Tingting Sun
   * \date  09.10.2012
   */

  class PolynomialFrustum : public MBSim::RigidContour {
    public:
      /*!
       * \brief Constructor
       * \param name   Name of the contour
       * \param param_
       */
      PolynomialFrustum(const std::string & name, const fmatvec::Vec & param_);

      /*!
       * \brief Destructor
       */
      virtual ~PolynomialFrustum();

      /* INHERITED INTERFACE OF ELEMENT */
      std::string getType() const {
        return "PolynomialFrustum";
      }
      virtual void init(InitStage stage);
      /***************************************************/

      /* INHERITED INTERFACE OF CONTOUR*/
      /*!
       * \brief returns the two lagrange parameters
       *    index 0: x,   height coordinate of frustum
       *    index 1: phi, angle of point
       */
      virtual fmatvec::Vec2 computeLagrangeParameter(const fmatvec::Vec3 &WrPoint);
      /*********************************/

      /*!
       * \brief set height of frustum
       * \param height height of the frustum
       */
      void setHeight(const double & height_);

      /*!
       * \brief return height of frustum
       */
      double getHeight();

      /*!
       * \brief return height of frustum
       */
      double getHeight() const;

#ifdef HAVE_OPENMBVCPPINTERFACE
      /*!
       * \brief enable visualisation
       * \param enable           enable or disable openmbv
       * \param polynomialPoints how fine should the grid be in polynomial direction
       * \param circulaPoints    how fine should the grid be in circular direction
       */
      BOOST_PARAMETER_MEMBER_FUNCTION( (void), enableOpenMBV, tag, (optional (polynomialPoints,(int),0)(circularPoints,(int),25)(diffuseColor,(const fmatvec::Vec3&),"[-1;1;1]")(transparency,(double),0))) { enableOpenMBV_(diffuseColor,transparency,polynomialPoints,circularPoints); }

      /*!
       * \brief set color of body
       */
      void setColor(const RGBColor & color);

      /*!
       * \brief set transparency of body
       */
      void setTransparency(const double & transparency);
#endif

      /*!
       * \brief get value at position of 0 derivative
       */
      double getValue(const double & x);

      /*!
       * \brief get value at position of 0 derivative
       */
      double getValue(const double & x) const;

      /*!
       * \brief get value at position of first derivative
       */
      double getValueD1(const double & x);

      /*!
       * \brief get value at position of first derivative
       */
      double getValueD1(const double & x) const;

      /*!
       * \brief get value at position of second derivative
       */
      double getValueD2(const double & x);

      /*!
       * \brief get value at position of second derivative
       */
      double getValueD2(const double & x) const;

      /*!
       * \brief get inner peak value of the polynomial
       */
      double getXPolyMax();

      /*!
       * \brief take the largest one among the distance from ((a+b)/2,0) to (a,f(a)),(b,f(b)) and (x,f(x)) where f(x) achieves maximum at x within[a,b]
       * \return the radius of the circumsphere
       */
      double getEnclosingSphereRadius();

      /*!
       * \brief return the center of the enclosing sphere
       */
      fmatvec::Vec3 getEnclosingSphereCenter();

      /*!
       * \brief get coefficient vector of the polynomial
       */
      const fmatvec::Vec & getPolynomialParameters();

      /*!
       * \brief returns the point in local coordinates of the frustum at the position (x, phi)
       */
      fmatvec::Vec3 computePoint(const double & x, const double & phi);

      /*!
       * \brief returns the normal in local coordinates of the frustum at the position (x, phi)
       */
      fmatvec::Vec3 computeNormal(const double & x, const double & phi);

      /*!
       * \brief returns the tangent in radial direction in local coordinates of the frustum at the position (x, phi)
       */
      fmatvec::Vec3 computeTangentRadial(const double & x, const double & phi);

      /*!
       * \brief returns the tangent in azimuthal direction in local coordinates of the frustum at the point x, phi
       */
      fmatvec::Vec3 computeTangentAzimuthal(const double & x, const double & phi);

      /*!
       * \brief in 2D plane, given a point outside a polynomial curve, search for the closest point on the curve to the point
       * \param x_0: starting point of the polynomial domain
       * \param x_end: end point of the polynomial domain
       * \param P: the given point
       * \return: a dimension 3 vector, vec[0] storing the dis, vec[1] and vec[2] for the 2D position of the closest point on the curve
       */
      fmatvec::Vec3 CP_toP_onPolycurve2D(double x_0, double x_end, fmatvec::Vec2 P);

    protected:
      /*!
       * \computes a sphere that totally encloses the frustum
       */
      void computeEnclosingSphere();

      /*!
       * \brief vector holding our parameters of the polynom describing the frustum (radius over height)
       * parameters=[a0,a1,a2,...an]
       */
      fmatvec::Vec parameters;

      /*!
       * \brief height of the frustum
       */
      double height;

      /*!
       * \brief radius of the enclosing sphere
       */
      double sphereRadius;

#ifdef HAVE_OPENMBVCPPINTERFACE
      /*!
       * \brief color values for the iv-body
       */
      RGBColor color;

      /*!
       * \brief transparency value for the body
       */
      double transparency;

      /*!
       * \brief grid points in polynomial direction
       */
      int polynomialPoints;

      /*!
       * \brief grid points in polynomial direction
       */
      int circularPoints;

      /*!
       * \brief create inventor file for visualisation
       */
      void createInventorFile();

      void enableOpenMBV_(const fmatvec::Vec3 &dc, double tp, int polynomialPoints, int circularPoints);
#endif
  };

  /*!
   * \brief this class denotes polynomial equation like this:a0+a1*x+a2*x^2+...+an*x^n=rhs
   * \param para: coefficient vector of the left side
   */

  class ContactPolyfun : public Function<double(double)> {
    public:
      /*!
       * \brief Constructor
       */
      ContactPolyfun(const double & rhs, const PolynomialFrustum * frustum);
      virtual ~ContactPolyfun() {
      }

      virtual double operator()(const double &x);
      void initializeUsingXML() {
      }

    protected:
      double rhs;
      const PolynomialFrustum * frustum;

  };

}

#endif /* POLYNOMIAL_FRUSTUM_H_ */
