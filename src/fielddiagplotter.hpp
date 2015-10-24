/*! \file fielddiagplotter.hpp
 *  \brief Non-interactive field diagnostic plotter
 */

/* Copyright (c) 2005-2011 Taneli Kalvas. All rights reserved.
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

#ifndef FIELDDIAGPLOTTER_HPP
#define FIELDDIAGPLOTTER_HPP 1


#include "geometry.hpp"

#include "fielddiagplot.hpp"
#include "plotter.hpp"



/*! \brief Non-interactive fielddiag plotter.
 *
 *  This class ties together Plotter, which provides basic graphics
 *  features including frame control, and FieldDiagPlot, which makes the
 *  fielddiag plot.
 */
class FieldDiagPlotter : public Plotter, public FieldDiagPlot {

    virtual void build_plot( void );

public:

    /*! \brief Constructor for field diagnostic plotter.
     */
    FieldDiagPlotter( const Geometry &geom );

    /*! \brief Destructor for field diagnostic plotter.
     */
    virtual ~FieldDiagPlotter();


};


#endif
