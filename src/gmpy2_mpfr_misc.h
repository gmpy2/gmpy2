/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpfr_misc.h                                                       *
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

#ifndef GMPY_MPFR_MISC_H
#define GMPY_MPFR_MISC_H

#ifdef __cplusplus
extern "C" {
#endif

static PyObject * GMPy_Real_F2Q(PyObject *x, PyObject *y, CTXT_Object *context);

static PyObject * GMPy_MPFR_Free_Cache(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_Can_Round(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_get_emax_max(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_get_max_precision(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_get_exp(PyObject *self, PyObject *other);
static PyObject * GMPy_MPFR_set_exp(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_set_sign(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_copy_sign(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_Integer_Ratio_Method(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_Mantissa_Exp_Method(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_Simple_Fraction_Method(PyObject *self, PyObject *args, PyObject *keywds);
static PyObject * GMPy_MPFR_set_nan(PyObject *self, PyObject *other);
static PyObject * GMPy_MPFR_set_inf(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_set_zero(PyObject *self, PyObject *args);

static PyObject * GMPy_MPFR_GetPrec_Attrib(MPFR_Object *self, void *closure);
static PyObject * GMPy_MPFR_GetRc_Attrib(MPFR_Object *self, void *closure);
static PyObject * GMPy_MPFR_GetImag_Attrib(MPFR_Object *self, void *closure);
static PyObject * GMPy_MPFR_GetReal_Attrib(MPFR_Object *self, void *closure);
static PyObject * GMPy_MPFR_Get_Mpmath_MPF_Tuple(MPFR_Object *self, void *closure);
static int        GMPy_MPFR_NonZero_Slot(MPFR_Object *self);
static PyObject * GMPy_MPFR_SizeOf_Method(PyObject *self, PyObject *other);

static PyObject * GMPy_MPFR_CheckRange(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_CheckRange(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_CheckRange(PyObject *self, PyObject *other);

static PyObject * GMPy_MPFR_Method_Round10(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif
