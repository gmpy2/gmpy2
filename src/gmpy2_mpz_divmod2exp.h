/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpz_divmod2exp.h                                                  *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP, MPFR, and MPC multiple precision           *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000 - 2009 Alex Martelli                                     *
 *                                                                         *
 * Copyright 2008 - 2025 Case Van Horsen                                   *
 *                                                                         *
 * This file is part of GMPY2.                                             *
 *                                                                         *
 * GMPY2 is free software: you can redistribute it and/or modify it under  *
 * the terms of the GNU Lesser General Public License as published by the  *
 * Free Software Foundation, either version 3 of the License, or (at your  *
 * option) any later version.                                              *
 *                                                                         *
 * GMPY2 is distributed in the hope that it will be useful, but WITHOUT    *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or   *
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public    *
 * License for more details.                                               *
 *                                                                         *
 * You should have received a copy of the GNU Lesser General Public        *
 * License along with GMPY2; if not, see <http://www.gnu.org/licenses/>    *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef GMPY_MPZ_DIVMOD2EXP_H
#define GMPY_MPZ_DIVMOD2EXP_H

#ifdef __cplusplus
extern "C" {
#endif

static PyObject * GMPy_MPZ_c_divmod_2exp(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_c_div_2exp(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_c_mod_2exp(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_f_divmod_2exp(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_f_div_2exp(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_f_mod_2exp(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_t_divmod_2exp(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_t_div_2exp(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_t_mod_2exp(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif
