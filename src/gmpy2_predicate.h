/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_predicate.h                                                       *
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

#ifndef GMPY_PREDICATE_H
#define GMPY_PREDICATE_H

#ifdef __cplusplus
extern "C" {
#endif
static PyObject * GMPy_RealWithType_Is_NAN(PyObject *x, int xtype, CTXT_Object *context);
static PyObject * GMPy_ComplexWithType_Is_NAN(PyObject *x, int xtype, CTXT_Object *context);
static PyObject * GMPy_Number_Method_Is_NAN(PyObject *self, PyObject *args);
static PyObject * GMPy_Number_Is_NAN(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Is_NAN(PyObject *self, PyObject *other);

static PyObject * GMPy_RealWithType_Is_Infinite(PyObject *x, int xtype, CTXT_Object *context);
static PyObject * GMPy_ComplexWithType_Is_Infinite(PyObject *x, int xtype, CTXT_Object *context);
static PyObject * GMPy_Number_Method_Is_Infinite(PyObject *self, PyObject *args);
static PyObject * GMPy_Number_Is_Infinite(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Is_Infinite(PyObject *self, PyObject *other);

static PyObject * GMPy_RealWithType_Is_Finite(PyObject *x, int xtype, CTXT_Object *context);
static PyObject * GMPy_ComplexWithType_Is_Finite(PyObject *x, int xtype, CTXT_Object *context);
static PyObject * GMPy_Number_Method_Is_Finite(PyObject *self, PyObject *args);
static PyObject * GMPy_Number_Is_Finite(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Is_Finite(PyObject *self, PyObject *other);

static PyObject * GMPy_RealWithType_Is_Zero(PyObject *x, int xtype, CTXT_Object *context);
static PyObject * GMPy_ComplexWithType_Is_Zero(PyObject *x, int xtype, CTXT_Object *context);
static PyObject * GMPy_Number_Method_Is_Zero(PyObject *self, PyObject *args);
static PyObject * GMPy_Number_Is_Zero(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Is_Zero(PyObject *self, PyObject *other);

static PyObject * GMPy_RealWithType_Is_Signed(PyObject *x, int xtype, CTXT_Object *context);
static PyObject * GMPy_MPFR_Is_Signed_Method(PyObject *self, PyObject *args);
static PyObject * GMPy_Number_Is_Signed(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Is_Signed(PyObject *self, PyObject *other);

static PyObject * GMPy_RealWithType_Is_Regular(PyObject *x, int xtype, CTXT_Object *context);
static PyObject * GMPy_MPFR_Is_Regular_Method(PyObject *self, PyObject *args);
static PyObject * GMPy_Number_Is_Regular(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Is_Regular(PyObject *self, PyObject *other);

static PyObject * GMPy_RealWithType_Is_Integer(PyObject *x, int xtype, CTXT_Object *context);
static PyObject * GMPy_MPFR_Is_Integer_Method(PyObject *self, PyObject *args);
static PyObject * GMPy_Number_Is_Integer(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Is_Integer(PyObject *self, PyObject *other);

static PyObject * GMPy_RealWithType_Is_LessGreater(PyObject *x, int xtype, PyObject *y, int ytype, CTXT_Object *context);
static PyObject * GMPy_Number_Is_LessGreater(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_Is_LessGreater(PyObject *self, PyObject *args);

static PyObject * GMPy_RealWithType_Is_Unordered(PyObject *x, int xtype, PyObject *y, int ytype, CTXT_Object *context);
static PyObject * GMPy_Number_Is_Unordered(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_Is_Unordered(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif
