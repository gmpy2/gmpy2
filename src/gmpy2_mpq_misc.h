/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpq_misc.h                                                        *
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

#ifndef GMPY_MPQ_MISC_H
#define GMPY_MPQ_MISC_H

#ifdef __cplusplus
extern "C" {
#endif

static PyObject * GMPy_MPQ_Attrib_GetNumer(MPQ_Object *self, void *closure);
static PyObject * GMPy_MPQ_Attrib_GetDenom(MPQ_Object *self, void *closure);
static PyObject * GMPy_MPQ_Attrib_GetReal(MPQ_Object *self, void *closure);
static PyObject * GMPy_MPQ_Attrib_GetImag(MPQ_Object *self, void *closure);
static PyObject * GMPy_MPQ_Method_As_Integer_Ratio(PyObject *self, PyObject *Py_UNUSED(ignored));
static PyObject * GMPy_MPQ_Method_From_As_Integer_Ratio(PyTypeObject *type, PyObject *const *args, Py_ssize_t nargs);

static PyObject * GMPy_MPQ_Function_Numer(PyObject *self, PyObject *other);
static PyObject * GMPy_MPQ_Function_Denom(PyObject *self, PyObject *other);
static PyObject * GMPy_MPQ_Function_Qdiv(PyObject *self, PyObject *args);
static PyObject * GMPy_MPQ_Method_Ceil(PyObject *self, PyObject *other);
static PyObject * GMPy_MPQ_Method_Floor(PyObject *self, PyObject *other);
static PyObject * GMPy_MPQ_Method_Trunc(PyObject *self, PyObject *other);
static PyObject * GMPy_MPQ_Method_Round(PyObject *self, PyObject *other);
static int        GMPy_MPQ_NonZero_Slot(MPQ_Object *x);

#ifdef __cplusplus
}
#endif
#endif
