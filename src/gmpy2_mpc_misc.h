/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpc_misc.h                                                        *
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

#ifndef GMPY_MPC_MISC_H
#define GMPY_MPC_MISC_H

#ifdef __cplusplus
extern "C" {
#endif

static PyObject * GMPy_Complex_Phase(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Phase(PyObject *self, PyObject *other);

#ifdef MPC_110
static PyObject * GMPy_Complex_Root_Of_Unity(PyObject *n, PyObject *k, CTXT_Object *context);
static PyObject * GMPy_Context_Root_Of_Unity(PyObject *self, PyObject *args);
#endif

static PyObject * GMPy_Complex_Norm(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Norm(PyObject *self, PyObject *other);

static PyObject * GMPy_Complex_Polar(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Polar(PyObject *self, PyObject *other);

static PyObject * GMPy_Complex_Rect(PyObject *r, PyObject *phi, CTXT_Object *context);
static PyObject * GMPy_Context_Rect(PyObject *self, PyObject *args);

static PyObject * GMPy_Complex_Proj(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Proj(PyObject *self, PyObject *other);

static PyObject * GMPy_MPC_Conjugate_Method(PyObject *self, PyObject *args);
static PyObject * GMPy_MPC_GetPrec_Attrib(MPC_Object *self, void *closure);
static PyObject * GMPy_MPC_GetRc_Attrib(MPC_Object *self, void *closure);
static PyObject * GMPy_MPC_GetImag_Attrib(MPC_Object *self, void *closure);
static PyObject * GMPy_MPC_GetReal_Attrib(MPC_Object *self, void *closure);
static int        GMPy_MPC_NonZero_Slot(MPC_Object *self);
static PyObject * GMPy_MPC_SizeOf_Method(PyObject *self, PyObject *other);

static PyObject * GMPy_MPC_Conjugate_Method(PyObject *self, PyObject *args);
static PyObject * GMPy_MPC_GetPrec_Attrib(MPC_Object *self, void *closure);
static PyObject * GMPy_MPC_GetRc_Attrib(MPC_Object *self, void *closure);
static PyObject * GMPy_MPC_GetImag_Attrib(MPC_Object *self, void *closure);
static PyObject * GMPy_MPC_GetReal_Attrib(MPC_Object *self, void *closure);
static int        GMPy_MPC_NonZero_Slot(MPC_Object *self);

#ifdef __cplusplus
}
#endif
#endif
