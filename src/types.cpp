/*! \file types.cpp
 *  \brief Base types
 */

/* Copyright (c) 2005-2012 Taneli Kalvas. All rights reserved.
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

#include "types.hpp"


const char *coordinate_axis_string[] = {
    "x", 
    "y", 
    "r", 
    "z" 
};


const char *coordinate_axis_string_with_unit[] = {
    "x (m)", 
    "y (m)", 
    "r (m)", 
    "z (m)" 
};


const char *trajectory_diagnostic_string[] = {
    "none",
    "t",
    "x",
    "v_x",
    "y",
    "r",
    "v_y",
    "v_r",
    "\\omega",
    "v_\\theta",
    "z",
    "v_z",
    "o",
    "v_o",
    "p",
    "v_p",
    "q",
    "v_q",
    "x\'",
    "y\'",
    "r\'",
    "\\alpha\'",
    "z\'",
    "o\'",
    "p\'",
    "I",
    "E_k",
    "q/m",
    "charge",
    "mass",
    "particle no"
};


const char *trajectory_diagnostic_string_with_unit[] = {
    "none ()",
    "t (s)",
    "x (m)",
    "v_x (m/s)",
    "y (m)",
    "r (m)",
    "v_y (m/s)",
    "v_r (m/s)",
    "\\omega (rad/s)",
    "v_\\theta (m/s)",
    "z (m)",
    "v_z (m/s)",
    "o (m)",
    "v_o (m/s)",
    "p (m)",
    "v_p (m/s)",
    "q (m)",
    "v_q (m/s)",
    "x\' (rad)",
    "y\' (rad)",
    "r\' (rad)",
    "\\alpha\' (rad)",
    "z\' (rad)",
    "o\' (rad)",
    "p\' (rad)",
    "I (A)",
    "E_k (eV)",
    "q/m (e/u)",
    "charge (e)",
    "mass (u)",
    "particle no"
};


const char *trajectory_diagnostic_string_unit[] = {
    "",
    "s",
    "m",
    "m/s",
    "m",
    "m",
    "m/s",
    "m/s",
    "rad/s",
    "m/s",
    "m",
    "m/s",
    "m",
    "m/s",
    "m",
    "m/s",
    "m",
    "m/s",
    "rad",
    "rad",
    "rad",
    "rad",
    "rad",
    "rad",
    "rad",
    "A",
    "eV",
    "e/u",
    "e",
    "u",
    ""
};


