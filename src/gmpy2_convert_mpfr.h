/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_convert_mpfr.h                                                    *
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

#ifndef GMPY2_CONVERT_MPFR_H
#define GMPY2_CONVERT_MPFR_H

#ifdef __cplusplus
extern "C" {
#endif

/********* Pympfr Conversions *********/

/* Conversions with Pympfr */

static MPFR_Object *    GMPy_MPFR_From_MPFR(MPFR_Object *obj, mpfr_prec_t prec, CTXT_Object *context);
static MPFR_Object *    GMPy_MPFR_From_PyLong(PyObject *obj, mpfr_prec_t prec, CTXT_Object *context);
static MPFR_Object *    GMPy_MPFR_From_PyFloat(PyObject *obj, mpfr_prec_t prec, CTXT_Object *context);
static MPFR_Object *    GMPy_MPFR_From_MPZ(MPZ_Object *obj, mpfr_prec_t prec, CTXT_Object *context);
static MPFR_Object *    GMPy_MPFR_From_MPQ(MPQ_Object *obj, mpfr_prec_t prec, CTXT_Object *context);
static MPFR_Object *    GMPy_MPFR_From_Fraction(PyObject *obj, mpfr_prec_t prec, CTXT_Object *context);
static MPFR_Object *    GMPy_MPFR_From_PyStr(PyObject *s, int base, mpfr_prec_t prec, CTXT_Object *context);
static MPFR_Object *    GMPy_MPFR_From_Real(PyObject* obj, mpfr_prec_t prec, CTXT_Object *context);
static MPFR_Object *    GMPy_MPFR_From_RealWithTypeAndCopy(PyObject* obj, int xtype, mpfr_prec_t prec, CTXT_Object *context);

static PyObject *       GMPy_PyLong_From_MPFR(MPFR_Object *obj, CTXT_Object *context);
static MPZ_Object *     GMPy_MPZ_From_MPFR(MPFR_Object *obj, CTXT_Object *context);
static XMPZ_Object *    GMPy_XMPZ_From_MPFR(MPFR_Object *self, CTXT_Object *context);
static MPQ_Object  *    GMPy_MPQ_From_MPFR(MPFR_Object *self, CTXT_Object *context);
static PyObject *       GMPy_PyFloat_From_MPFR(MPFR_Object *self, CTXT_Object *context);
static PyObject *       GMPy_PyStr_From_MPFR(MPFR_Object *self, int base, int digits, CTXT_Object *context);

static PyObject *       GMPy_MPFR_Str_Slot(MPFR_Object *self);
static PyObject *       GMPy_MPFR_Repr_Slot(MPFR_Object *self);
static PyObject *       GMPy_MPFR_Int_Slot(MPFR_Object *self);
static PyObject *       GMPy_MPFR_Float_Slot(MPFR_Object *self);

/* Miscellaneous */
#ifdef SHARED
/* static int           GMPy_MPFR_ConvertArg(PyObject *arg, PyObject **ptr); */
static GMPy_MPFR_ConvertArg_RETURN GMPy_MPFR_ConvertArg GMPy_MPFR_ConvertArg_PROTO;
#endif

static PyObject *       stern_brocot(MPFR_Object* self, MPFR_Object *err, mpfr_prec_t prec, int mayz, CTXT_Object *context);
static PyObject *       mpfr_ascii(mpfr_t self, int base, int digits, int round);

#ifdef __cplusplus
}
#endif
#endif
