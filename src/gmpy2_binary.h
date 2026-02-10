/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_binary.c                                                          *
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

#ifndef GMPY_BINARY_H
#define GMPY_BINARY_H

#ifdef __cplusplus
extern "C" {
#endif

/* Conversion routines between GMPY2 objects and a compact, portable
 * binary representation. The binary format of GMPY2 is not compatible
 * with GMPY 1.x. Methods to read the old format are provided.
 */

static PyObject * GMPy_MPZ_From_Old_Binary(PyObject *self, PyObject *other);
static PyObject * GMPy_MPQ_From_Old_Binary(PyObject *self, PyObject *other);
static PyObject * GMPy_MPFR_From_Old_Binary(PyObject *self, PyObject *other);

static PyObject * GMPy_MPANY_From_Binary(PyObject *self, PyObject *other);
static PyObject * GMPy_MPANY_To_Binary(PyObject *self, PyObject *other);

static PyObject * GMPy_MPZ_To_Binary(MPZ_Object *self);
static PyObject * GMPy_XMPZ_To_Binary(XMPZ_Object *self);
static PyObject * GMPy_MPQ_To_Binary(MPQ_Object *self);
static PyObject * GMPy_MPFR_To_Binary(MPFR_Object *self);
static PyObject * GMPy_MPC_To_Binary(MPC_Object *self);

#ifdef __cplusplus
}
#endif
#endif
