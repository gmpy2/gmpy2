/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_convert_gmp.h                                                      *
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

#ifndef GMPY2_CONVERT_GMP_H
#define GMPY2_CONVERT_GMP_H

#ifdef __cplusplus
extern "C" {
#endif

/* ======================================================================== *
 * Conversion between native Python objects and MPZ.                        *
 * ======================================================================== */

static MPZ_Object *    GMPy_MPZ_From_PyLong(PyObject *obj, CTXT_Object *context);
static MPZ_Object *    GMPy_MPZ_From_PyStr(PyObject *s, int base, CTXT_Object *context);
static MPZ_Object *    GMPy_MPZ_From_PyFloat(PyObject *obj, CTXT_Object *context);

static MPZ_Object *    GMPy_MPZ_From_Integer(PyObject *obj, CTXT_Object *context);
static MPZ_Object *    GMPy_MPZ_From_IntegerAndCopy(PyObject *obj, CTXT_Object *context);
static MPZ_Object *    GMPy_MPZ_From_IntegerWithType(PyObject *obj, int xtype, CTXT_Object *context);
static MPZ_Object *    GMPy_MPZ_From_IntegerWithTypeAndCopy(PyObject *obj, int xtype, CTXT_Object *context);

static PyObject *      GMPy_MPZ_Str_Slot(MPZ_Object *self);
static PyObject *      GMPy_MPZ_Repr_Slot(MPZ_Object *self);

static PyObject *      GMPy_PyLong_From_MPZ(MPZ_Object *obj, CTXT_Object *context);
static PyObject *      GMPy_PyFloat_From_MPZ(MPZ_Object *obj, CTXT_Object *context);
static PyObject *      GMPy_PyStr_From_MPZ(MPZ_Object *obj, int base, int option, CTXT_Object *context);

#ifdef SHARED
/* static int             GMPy_MPZ_ConvertArg(PyObject *arg, PyObject **ptr); */
static GMPy_MPZ_ConvertArg_RETURN GMPy_MPZ_ConvertArg GMPy_MPZ_ConvertArg_PROTO;
#endif


/* ======================================================================== *
 * Conversion between native Python objects/MPZ and XMPZ.                   *
 * ======================================================================== */

static XMPZ_Object *   GMPy_XMPZ_From_PyLong(PyObject *obj, CTXT_Object *context);
static XMPZ_Object *   GMPy_XMPZ_From_PyStr(PyObject *s, int base, CTXT_Object *context);
static XMPZ_Object *   GMPy_XMPZ_From_PyFloat(PyObject *self, CTXT_Object *context);
static XMPZ_Object *   GMPy_XMPZ_From_MPZ(MPZ_Object *obj, CTXT_Object *context);
static XMPZ_Object *   GMPy_XMPZ_From_XMPZ(XMPZ_Object *obj, CTXT_Object *context);

static PyObject *      GMPy_XMPZ_Str_Slot(XMPZ_Object *self);
static PyObject *      GMPy_XMPZ_Repr_Slot(XMPZ_Object *self);

static PyObject *      GMPy_PyStr_From_XMPZ(XMPZ_Object *self, int base, int option, CTXT_Object *context);
static MPZ_Object *    GMPy_MPZ_From_XMPZ(XMPZ_Object *obj, CTXT_Object *context);

/* ======================================================================== *
 * Conversion between native Python objects/MPZ/XMPZ and MPQ.               *
 * ======================================================================== */

static MPQ_Object *    GMPy_MPQ_From_PyLong(PyObject *obj, CTXT_Object *context);
static MPQ_Object *    GMPy_MPQ_From_PyStr(PyObject *s, int base, CTXT_Object *context);
static MPQ_Object *    GMPy_MPQ_From_PyFloat(PyObject *obj, CTXT_Object *context);
static MPQ_Object *    GMPy_MPQ_From_Fraction(PyObject *obj, CTXT_Object *context);
static MPQ_Object *    GMPy_MPQ_From_MPZ(MPZ_Object *obj, CTXT_Object *context);
static MPQ_Object *    GMPy_MPQ_From_XMPZ(XMPZ_Object *obj, CTXT_Object *context);

static MPQ_Object *    GMPy_MPQ_From_Rational(PyObject* obj, CTXT_Object *context);
static MPQ_Object *    GMPy_MPQ_From_RationalAndCopy(PyObject* obj, CTXT_Object *context);
static MPQ_Object *    GMPy_MPQ_From_NumberWithType(PyObject* obj, int xtype, CTXT_Object *context);
static MPQ_Object *    GMPy_MPQ_From_Number(PyObject* obj, CTXT_Object *context);
static MPQ_Object *    GMPy_MPQ_From_RationalWithTypeAndCopy(PyObject* obj, int xtype, CTXT_Object *context);
static MPQ_Object *    GMPy_MPQ_From_RationalWithType(PyObject *obj, int xtype, CTXT_Object *context);
static PyObject *      GMPy_PyLong_From_MPQ(MPQ_Object *obj, CTXT_Object *context);
static PyObject *      GMPy_PyStr_From_MPQ(MPQ_Object *obj, int base, int option, CTXT_Object *context);
static PyObject *      GMPy_PyFloat_From_MPQ(MPQ_Object *obj, CTXT_Object *context);
static MPZ_Object *    GMPy_MPZ_From_MPQ(MPQ_Object *obj, CTXT_Object *context);
static XMPZ_Object *   GMPy_XMPZ_From_MPQ(MPQ_Object *obj, CTXT_Object *context);

/* support str() and repr() */
static PyObject *      GMPy_MPQ_Str_Slot(MPQ_Object *obj);
static PyObject *      GMPy_MPQ_Repr_Slot(MPQ_Object *obj);
static PyObject *      GMPy_MPQ_Float_Slot(MPQ_Object *obj);
static PyObject *      GMPy_MPQ_Int_Slot(MPQ_Object *obj);

#ifdef SHARED
/* int GMPy_MPQ_convert_arg(PyObject *arg, PyObject **ptr); */
static GMPy_MPQ_ConvertArg_RETURN GMPy_MPQ_ConvertArg GMPy_MPQ_ConvertArg_PROTO;
#endif

#ifdef __cplusplus
}
#endif
#endif
