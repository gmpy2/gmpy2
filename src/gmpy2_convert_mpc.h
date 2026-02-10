/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_convert_mpc.h                                                     *
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

#ifndef GMPY2_CONVERT_MPC_H
#define GMPY2_CONVERT_MPC_H

#ifdef __cplusplus
extern "C" {
#endif

/* Conversions with Pympc */

static MPC_Object *   GMPy_MPC_From_MPC(MPC_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_PyComplex(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_MPFR(MPFR_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_PyFloat(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_MPZ(MPZ_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_MPQ(MPQ_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_Fraction(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_PyLong(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_PyStr(PyObject *s, int base, mpfr_prec_t rbits, mpfr_prec_t ibits, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_Complex(PyObject* obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
#if 0
static MPC_Object *   GMPy_MPC_From_ComplexAndCopy(PyObject* obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
#endif

static PyObject *     GMPy_MPC_Float_Slot(PyObject *self);
static PyObject *     GMPy_MPC_Int_Slot(PyObject *self);
static PyObject *     GMPy_PyStr_From_MPC(MPC_Object *self, int base, int digits, CTXT_Object *context);
static PyObject *     GMPy_PyComplex_From_MPC(PyObject *self, PyObject *other);
/* support str() and repr() */
static PyObject *     GMPy_MPC_Str_Slot(MPC_Object *self);
static PyObject *     GMPy_MPC_Repr_Slot(MPC_Object *self);

/* Miscellaneous */
#ifdef SHARED
/* static int            GMPy_MPC_convert_arg(PyObject *arg, PyObject **ptr); */
static GMPy_MPC_ConvertArg_RETURN GMPy_MPC_ConvertArg GMPy_MPC_ConvertArg_PROTO;

#endif

#ifdef __cplusplus
}
#endif
#endif
