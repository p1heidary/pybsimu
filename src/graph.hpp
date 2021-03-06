/*! \file graph.hpp
 *  \brief Base for plottable graphs
 */

/* Copyright (c) 2005-2009,2011 Taneli Kalvas. All rights reserved.
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

#ifndef GRAPH_HPP
#define GRAPH_HPP 1


#include <cairo.h>
#include "coordmapper.hpp"


/*! \brief Abstract base class for drawable plots.
 *
 *  Graph type objects are used in Frame type plots to provide the
 *  graphical view of the data content and legend sample.
 */
class Graph {

public:

    /*! \brief Virtual destructor.
     */
    virtual ~Graph() {}

    /*! \brief Plot graph with cairo.
     *
     *  Plot the graph using \a cairo and coordinate mapper \a cm. The
     *  visible range of plot is given in array \a range in order \a
     *  xmin, \a ymin, \a xmax, \a ymax. The graph should be able to
     *  handle any range values. Also \a min > \a max.
     *
     *  Called by Frame during drawing.
     */
    virtual void plot( cairo_t *cairo, const Coordmapper *cm, const double range[4] ) = 0;

    /*! \brief Plot sample for legend.
     *
     *  Plot graph sample for legend at cairo coordinates \a (x,y).
     */
    virtual void plot_sample( cairo_t *cairo, double x, double y, double width, double height ) = 0;

    /*! \brief Get bounding box of graph.
     *
     *  Returns the bounding box of the graph in array \a bbox in
     *  order \a xmin, \a ymin, \a xmax, \a ymax.
     */
    virtual void get_bbox( double bbox[4] ) = 0;
};


#endif


