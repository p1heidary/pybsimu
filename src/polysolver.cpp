/*! \file polysolver.cpp
 *  \brief Polynomial solver
 */

/* Copyright (c) 2005-2009,2012 Taneli Kalvas. All rights reserved.
 *
 * You can redistribute this software and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * 
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this library (file "COPYING" included in the package);
 * if not, write to the Free Software Foundation, Inc., 51 Franklin
 * Street, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Technology Transfer
 * Department at TTD@lbl.gov. Other questions, comments and bug
 * reports should be sent directly to the author via email at
 * taneli.kalvas@jyu.fi.
 * 
 * NOTICE. This software was developed under partial funding from the
 * U.S.  Department of Energy.  As such, the U.S. Government has been
 * granted for itself and others acting on its behalf a paid-up,
 * nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and
 * display publicly.  Beginning five (5) years after the date
 * permission to assert copyright is obtained from the U.S. Department
 * of Energy, and subject to any subsequent five (5) year renewals,
 * the U.S. Government is granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in
 * the Software to reproduce, prepare derivative works, distribute
 * copies to the public, perform publicly and display publicly, and to
 * permit others to do so.
 */

/*! \file polysolver.cpp
 *  \brief Polynomial solvers
 * 
 *  Copyright (C) 2003 CERN, K. S. K\"{o}lbig and T. Kalvas
 *
 *  Quadratic and quartic solvers converted to C and implemented into
 *  the GSL-extras library by Andrew W. Steiner and Andy Buckley
 *
 *  Cubic solver taken from gsl-1.12 by T. Kalvas
 * 
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or (at
 *  your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


#include <iostream>
#include <cmath>
#include "polysolver.hpp"


//#define DEBUG_POLYSOLVER 1


#define SWAPD(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)


inline double MAXD( double a, double b )
{
    return( a > b ? a : b );
}


uint32_t solve_quadratic( double a, double b, double c, 
			  double *x0, double *x1 )
{
#ifdef DEBUG_POLYSOLVER
    std::cout << "solve_quadratic( " 
	      << a << ", "
	      << b << ", "
	      << c << " )\n";
#endif

    // Handle linear case
    if( a == 0.0 ) {
	if( b == 0.0 ) {
#ifdef DEBUG_POLYSOLVER
	    std::cout << "linear with no roots\n";
#endif
	    return( 0 );
        } else {
	    *x0 = -c / b;
#ifdef DEBUG_POLYSOLVER
	    std::cout << "linear with one root at x = " << *x0 << "\n";
#endif
	    return( 1 );
        }
    }
    
    double disc = b*b - 4.0*a*c;
    if( disc > 0 ) {
	if( b == 0 ) {
	    double r = fabs( 0.5 * sqrt (disc) / a );
	    *x0 = -r;
	    *x1 =  r;
        } else {
	    double sgnb = (b > 0 ? 1 : -1);
	    double temp = -0.5 * (b + sgnb * sqrt (disc));
	    double r1 = temp / a ;
	    double r2 = c / temp ;
	    
	    if( r1 < r2 ) {
		*x0 = r1 ;
		*x1 = r2 ;
            } else {
		*x0 = r2 ;
		*x1 = r1 ;
            }
        }
	return( 2 );
    } else if( disc == 0 ) {
	*x0 = -0.5 * b / a ;
	*x1 = -0.5 * b / a ;
	return( 2 );
    } else {
	return( 0 );
    }
}


/* Finds the real roots of a x^3 + b x^2 + c x + d = 0
 */
uint32_t solve_cubic( double a, double b, double c, double d, 
		      double *x0, double *x1, double *x2 )
{
#ifdef DEBUG_POLYSOLVER
    std::cout << "solve_cubic( " 
	      << a << ", "
	      << b << ", "
	      << c << ", "
	      << d << " )\n";
#endif

    // Handle quadratic/linear case
    if( a == 0.0 )
	return( solve_quadratic( b, c, d, x0, x1 ) );

    double inva = 1.0/a;
    a = inva*b;
    b = inva*c;
    c = inva*d;

    double q = (a * a - 3 * b);
    double r = (2 * a * a * a - 9 * a * b + 27 * c);

    double Q = q / 9;
    double R = r / 54;

    double Q3 = Q * Q * Q;
    double R2 = R * R;

    double CR2 = 729 * r * r;
    double CQ3 = 2916 * q * q * q;

    if( R == 0 && Q == 0 ) {

	*x0 = - a / 3;
	*x1 = - a / 3;
	*x2 = - a / 3;
	return( 3 );

    } else if( CR2 == CQ3 ) {

	double sqrtQ = sqrt( Q );
	if( R > 0 ) {
	    *x0 = -2 * sqrtQ  - a / 3;
	    *x1 = sqrtQ - a / 3;
	    *x2 = sqrtQ - a / 3;
        } else {
	    *x0 = - sqrtQ  - a / 3;
	    *x1 = - sqrtQ - a / 3;
	    *x2 = 2 * sqrtQ - a / 3;
        }
	return( 3 );

    } else if( CR2 < CQ3 ) {

	double sqrtQ = sqrt( Q );
	double sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
	double theta = acos( R / sqrtQ3 );
	double norm = -2 * sqrtQ;
	*x0 = norm * cos( theta / 3 ) - a / 3;
	*x1 = norm * cos( (theta + 2.0 * M_PI) / 3 ) - a / 3;
	*x2 = norm * cos( (theta - 2.0 * M_PI) / 3 ) - a / 3;
      
	// Sort *x0, *x1, *x2 into increasing order
	if( *x0 > *x1 )
	    SWAPD( *x0, *x1 );
	if( *x1 > *x2 ) {
	    SWAPD( *x1, *x2 );          
	    if( *x0 > *x1 )
		SWAPD( *x0, *x1 );
        }
	return( 3 );

    } else {

	double sgnR = (R >= 0 ? 1 : -1);
	double A = -sgnR * pow( fabs( R ) + sqrt( R2 - Q3 ), 1.0/3.0 );
	double B = Q / A ;
	*x0 = A + B - a / 3;
	return( 1 );
    }
}


/* Finds the real roots of x^4 + a x^3 + b x^2 + c x + d = 0
 */
uint32_t solve_quartic( double a, double b, double c, double d,
		      double *x0, double *x1, double *x2, double *x3 )
{
    double u[3];
    double aa, pp, qq, rr, rc, sc, tc, mt;
    double w1r, w1i, w2r, w2i, w3r;
    double v[3], v1, v2, arg, theta;
    double disc, h;
    int k1 = 0, k2 = 0;
    double zarr[4];

    /* Deal easily with the cases where the quartic is degenerate. The
     * ordering of solutions is done explicitly. */
    if (0 == b && 0 == c)
    {
	if (0 == d)
        {
	    if (a > 0)
            {
		*x0 = -a;
		*x1 = 0.0;
		*x2 = 0.0;
		*x3 = 0.0;
            }
	    else
            {
		*x0 = 0.0;
		*x1 = 0.0;
		*x2 = 0.0;
		*x3 = -a;
            }
	    return 4;
        }
	else if (0 == a)
        {
	    if (d > 0)
            {
		return 0;
            }
	    else
            {
		*x1 = sqrt (sqrt (-d));
		*x0 = -(*x1);
		return 2;
            }
        }
    }

    if (0.0 == c && 0.0 == d)
    {
	*x0=0.0;
	*x1=0.0;
	if( solve_quadratic( 1.0, a, b, x2, x3 ) == 0 ) {
	    mt=3;
	} else {
	    mt=1;
	}
    }
    else 
    {
	/* For non-degenerate solutions, proceed by constructing and
	 * solving the resolvent cubic */
	aa = a * a;
	pp = b - (3.0/8.0) * aa;
	qq = c - (1.0/2.0) * a * (b - (1.0/4.0) * aa);
	rr = d - (1.0/4.0) * (a * c - (1.0/4.0) * aa * (b - (3.0/16.0) * aa));
	rc = (1.0/2.0) * pp;
	sc = (1.0/4.0) * ((1.0/4.0) * pp * pp - rr);
	tc = -((1.0/8.0) * qq * (1.0/8.0) * qq);

	/* This code solves the resolvent cubic in a convenient fashion
	 * for this implementation of the quartic. If there are three real
	 * roots, then they are placed directly into u[].  If two are
	 * complex, then the real root is put into u[0] and the real
	 * and imaginary part of the complex roots are placed into
	 * u[1] and u[2], respectively. Additionally, this
	 * calculates the discriminant of the cubic and puts it into the
	 * variable disc. */
	{
	    double qcub = (rc * rc - 3 * sc);
	    double rcub = (2 * rc * rc * rc - 9 * rc * sc + 27 * tc);

	    double Q = qcub / 9;
	    double R = rcub / 54;

	    double Q3 = Q * Q * Q;
	    double R2 = R * R;

	    double CR2 = 729 * rcub * rcub;
	    double CQ3 = 2916 * qcub * qcub * qcub;

	    disc = (CR2 - CQ3) / 2125764.0;

	    if (0 == R && 0 == Q)
	    {
		u[0] = -rc / 3;
		u[1] = -rc / 3;
		u[2] = -rc / 3;
	    }
	    else if (CR2 == CQ3)
	    {
		double sqrtQ = sqrt (Q);
		if (R > 0)
		{
		    u[0] = -2 * sqrtQ - rc / 3;
		    u[1] = sqrtQ - rc / 3;
		    u[2] = sqrtQ - rc / 3;
		}
		else
		{
		    u[0] = -sqrtQ - rc / 3;
		    u[1] = -sqrtQ - rc / 3;
		    u[2] = 2 * sqrtQ - rc / 3;
		}
	    }
	    else if (CR2 < CQ3)
	    {
		double sqrtQ = sqrt (Q);
		double sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
		double theta = acos (R / sqrtQ3);
		if (R / sqrtQ3 >= 1.0) theta = 0.0;
		{
		    double norm = -2 * sqrtQ;
	  
		    u[0] = norm * cos (theta / 3) - rc / 3;
		    u[1] = norm * cos ((theta + 2.0 * M_PI) / 3) - rc / 3;
		    u[2] = norm * cos ((theta - 2.0 * M_PI) / 3) - rc / 3;
		}
	    }
	    else
	    {
		double sgnR = (R >= 0 ? 1 : -1);
		double modR = fabs (R);
		double sqrt_disc = sqrt (R2 - Q3);
		double A = -sgnR * pow (modR + sqrt_disc, 1.0 / 3.0);
		double B = Q / A;
		double mod_diffAB = fabs (A - B);

		u[0] = A + B - rc / 3;
		u[1] = -0.5 * (A + B) - rc / 3;
		u[2] = -(sqrt (3.0) / 2.0) * mod_diffAB;
	    }
	}
	/* End of solution to resolvent cubic */

	/* Combine the square roots of the roots of the cubic 
	 * resolvent appropriately. Also, calculate 'mt' which 
	 * designates the nature of the roots:
	 * mt=1 : 4 real roots (disc == 0)
	 * mt=2 : 0 real roots (disc < 0)
	 * mt=3 : 2 real roots (disc > 0)
	 */

	if (0.0 == disc) 
	    u[2] = u[1];

	if (0 >= disc)
	{
	    mt = 2; 

	    /* One would think that we could return 0 here and exit,
	     * since mt=2. However, this assignment is temporary and
	     * changes to mt=1 under certain conditions below.  
	     */
	  
	    v[0] = fabs (u[0]);
	    v[1] = fabs (u[1]);
	    v[2] = fabs (u[2]);

	    v1 = MAXD( MAXD( v[0], v[1] ), v[2] );
	    /* Work out which two roots have the largest moduli */
	    k1 = 0, k2 = 0;
	    if (v1 == v[0])
	    {
		k1 = 0;
		v2 = MAXD( v[1], v[2] );
	    }
	    else if (v1 == v[1])
	    {
		k1 = 1;
		v2 = MAXD( v[0], v[2] );
	    }
	    else
	    {
		k1 = 2;
		v2 = MAXD( v[0], v[1] );
	    }

	    if (v2 == v[0])
	    {
		k2 = 0;
	    }
	    else if (v2 == v[1])
	    {
		k2 = 1;
	    }
	    else
	    {
		k2 = 2;
	    }
	  
	    if (0.0 <= u[k1]) 
	    {
		w1r=sqrt(u[k1]);
		w1i=0.0;
	    } 
	    else 
	    {
		w1r=0.0;
		w1i=sqrt(-u[k1]);
	    }
	    if (0.0 <= u[k2]) 
	    {
		w2r=sqrt(u[k2]);
		w2i=0.0;
	    } 
	    else 
	    {
		w2r=0.0;
		w2i=sqrt(-u[k2]);
	    }
	}
	else
	{
	    mt = 3;

	    if (0.0 == u[1] && 0.0 == u[2]) 
	    {
		arg = 0.0;
	    } 
	    else 
	    {
		arg = sqrt(sqrt(u[1] * u[1] + u[2] * u[2]));
	    }
	    theta = atan2(u[2], u[1]);
	  
	    w1r = arg * cos(theta / 2.0);
	    w1i = arg * sin(theta / 2.0);
	    w2r = w1r;
	    w2i = -w1i;
	}
  
	/* Solve the quadratic to obtain the roots to the quartic */
	w3r = qq / 8.0 * (w1i * w2i - w1r * w2r) / 
	    (w1i * w1i + w1r * w1r) / (w2i * w2i + w2r * w2r);
	h = a / 4.0;

	zarr[0] = w1r + w2r + w3r - h;
	zarr[1] = -w1r - w2r + w3r - h;
	zarr[2] = -w1r + w2r - w3r - h;
	zarr[3] = w1r - w2r - w3r - h;
      
	/* Arrange the roots into the variables z0, z1, z2, z3 */
	if (2 == mt)
        {
	    if (u[k1] >= 0 && u[k2] >= 0)
            {
		mt = 1;
		*x0 = zarr[0];
		*x1 = zarr[1];
		*x2 = zarr[2];
		*x3 = zarr[3];
            }
	    else
	    {
		return 0;
	    }
	}
	else 
        {
	    *x0 = zarr[0];
	    *x1 = zarr[1];
        }
    }
  
    /* Sort the roots as usual */
    if (1 == mt)
    {
	/* Roots are all real, sort them by the real part */
	if (*x0 > *x1)
	    SWAPD (*x0, *x1);
	if (*x0 > *x2)
	    SWAPD (*x0, *x2);
	if (*x0 > *x3)
	    SWAPD (*x0, *x3);

	if (*x1 > *x2)
	    SWAPD (*x1, *x2);
	if (*x2 > *x3)
        {
	    SWAPD (*x2, *x3);
	    if (*x1 > *x2)
		SWAPD (*x1, *x2);
        }
	return 4;
    }
    else
    {
	/* 2 real roots */
	if (*x0 > *x1)
	    SWAPD (*x0, *x1);
    }

    return 2;
}
