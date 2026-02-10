/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_pow.h                                                             *
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

#ifndef GMPY2_POW_H
#define GMPY2_POW_H

#ifdef __cplusplus
extern "C" {
#endif

/* Private API */

static PyObject * GMPy_Number_Pow_Slot(PyObject *base, PyObject *exp, PyObject *mod);

static PyObject * GMPy_Integer_PowWithType(PyObject *base, int btype, PyObject *exp, int etype, PyObject *mod, CTXT_Object *context);
static PyObject * GMPy_Rational_PowWithType(PyObject *base, int btype, PyObject *exp, int etype, PyObject *mod, CTXT_Object *context);
static PyObject * GMPy_Real_PowWithType(PyObject *base, int btype, PyObject *exp, int etype, PyObject *mod, CTXT_Object *context);
static PyObject * GMPy_Complex_PowWithType(PyObject *base, int btype, PyObject *exp, int etype, PyObject *mod, CTXT_Object *context);
static PyObject * GMPy_Integer_PowMod(PyObject *self, PyObject *args);
static PyObject * GMPy_Integer_PowMod_Sec(PyObject *self, PyObject *args);

static PyObject * GMPy_Context_Pow(PyObject *self, PyObject *args);
static PyObject * GMPy_Number_Pow(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context);

#ifdef __cplusplus
}
#endif
#endif
