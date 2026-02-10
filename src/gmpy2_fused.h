/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_fused.h                                                           *
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

#ifndef GMPY_FUSED_H
#define GMPY_FUSED_H

#ifdef __cplusplus
extern "C" {
#endif

static PyObject * GMPy_IntegerWithType_FMA(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, CTXT_Object *context);
static PyObject * GMPy_RationalWithType_FMA(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, CTXT_Object *context);
static PyObject * GMPy_RealWithType_FMA(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, CTXT_Object *context);
static PyObject * GMPy_ComplexWithType_FMA(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, CTXT_Object *context);
static PyObject * GMPy_Number_FMA(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context);
static PyObject * GMPy_Context_FMA(PyObject *self, PyObject *args);

static PyObject * GMPy_IntegerWithType_FMS(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, CTXT_Object *context);
static PyObject * GMPy_RationalWithType_FMS(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, CTXT_Object *context);
static PyObject * GMPy_RealWithType_FMS(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, CTXT_Object *context);
static PyObject * GMPy_ComplexWithType_FMS(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, CTXT_Object *context);
static PyObject * GMPy_Number_FMS(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context);
static PyObject * GMPy_Context_FMS(PyObject *self, PyObject *args);

#if MPFR_VERSION_MAJOR > 3
static PyObject * GMPy_IntegerWithType_FMMA(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, PyObject *t, int ttype, CTXT_Object *context);
static PyObject * GMPy_RationalWithType_FMMA(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, PyObject *t, int ttype,CTXT_Object *context);
static PyObject * GMPy_RealWithType_FMMA(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, PyObject *t, int ttype,CTXT_Object *context);
static PyObject * GMPy_Number_FMMA(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context);
static PyObject * GMPy_Context_FMMA(PyObject *self, PyObject *args);

static PyObject * GMPy_IntegerWithType_FMMS(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, PyObject *t, int ttype, CTXT_Object *context);
static PyObject * GMPy_RationalWithType_FMMS(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, PyObject *t, int ttype,CTXT_Object *context);
static PyObject * GMPy_RealWithType_FMMS(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, PyObject *t, int ttype,CTXT_Object *context);
static PyObject * GMPy_Number_FMMS(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context);
static PyObject * GMPy_Context_FMMS(PyObject *self, PyObject *args);
#endif

#ifdef __cplusplus
}
#endif
#endif
