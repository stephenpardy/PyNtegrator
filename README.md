MCOrbits
========

Orbital integration and likelihood testing for MCMC parameter searches of tidal fields.

Python code acts as a wrapper for integration in C and for the emcee package which searches parameter space.

Naive setup and running:

python setup.py build_ext --inplace

For a single temperature:
python mcorbit_serial.py some_name_for_output some_IC_file

For multiple temperatures:
python mcorbit_parallel.py some_name_for_output some_IC_file

See ICS_Example for example of initial conditions file structure

/***************************************************************************
 *   Copyright (C) 2013, 2014 by Stephen Pardy, Ana Bonaca                 * 
 *                               &  Andreas Kuepper                        *
 *   contact: spardy@astro.wisc.edu                                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
