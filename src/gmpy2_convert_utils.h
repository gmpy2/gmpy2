/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_convert_utils.h                                                    *
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

#ifndef GMPY2_CONVERT_UTILS_H
#define GMPY2_CONVERT_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

/* ======================================================================== *
 * Conversion between Integer objects and C types.                          *
 * ======================================================================== */

static long            GMPy_Integer_AsLongWithType(PyObject *x, int xtype);
static long            GMPy_Integer_AsLong(PyObject *x);
static unsigned long   GMPy_Integer_AsUnsignedLongWithType(PyObject *x, int xtype);
static unsigned long   GMPy_Integer_AsUnsignedLong(PyObject *x);
static long            GMPy_Integer_AsUnsignedLongOrLong(PyObject *x, int *is_signed);

static PY_LONG_LONG    GMPy_Integer_AsLongLongWithType(PyObject *x, int xtype);
static PY_LONG_LONG    GMPy_Integer_AsLongLong(PyObject *x);

static mp_bitcnt_t     GMPy_Integer_AsMpBitCnt(PyObject *x);
static mp_bitcnt_t     GMPy_PyLong_AsMpBitCnt(PyObject *x);
static PyObject *      GMPy_PyLong_FromMpBitCnt(mp_bitcnt_t n);

static unsigned PY_LONG_LONG GMPy_Integer_AsUnsignedLongLongWithType(PyObject *x, int xtype);
static unsigned PY_LONG_LONG GMPy_Integer_AsUnsignedLongLong(PyObject *x);

#define GMPy_Integer_AsSsize_t (Py_ssize_t)GMPy_Integer_AsLongLong

static mp_bitcnt_t GMPy_Integer_AsMpBitCnt(PyObject *x);

#ifdef __cplusplus
}
#endif
#endif
