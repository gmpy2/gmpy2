/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_xmpz_inplace.c                                                    *
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

#ifndef GMPY_XMPZ_INPLACE_H
#define GMPY_XMPZ_INPLACE_H

#ifdef __cplusplus
extern "C" {
#endif

static PyObject * GMPy_XMPZ_IAdd_Slot(PyObject *self, PyObject *other);
static PyObject * GMPy_XMPZ_ISub_Slot(PyObject *self, PyObject *other);
static PyObject * GMPy_XMPZ_IMul_Slot(PyObject *self, PyObject *other);
static PyObject * GMPy_XMPZ_IFloorDiv_Slot(PyObject *self, PyObject *other);
static PyObject * GMPy_XMPZ_IRem_Slot(PyObject *self, PyObject *other);
static PyObject * GMPy_XMPZ_IRshift_Slot(PyObject *self, PyObject *other);
static PyObject * GMPy_XMPZ_ILshift_Slot(PyObject *self, PyObject *other);
static PyObject * GMPy_XMPZ_IPow_Slot(PyObject *self, PyObject *other, PyObject *mod);
static PyObject * GMPy_XMPZ_IAnd_Slot(PyObject *self, PyObject *other);
static PyObject * GMPy_XMPZ_IXor_Slot(PyObject *self, PyObject *other);
static PyObject * GMPy_XMPZ_IIor_Slot(PyObject *self, PyObject *other);

#ifdef __cplusplus
}
#endif
#endif
