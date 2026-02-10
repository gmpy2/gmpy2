/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpz_bitops.h                                                      *
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

#ifndef GMPY_MPZ_BITOPS_H
#define GMPY_MPZ_BITOPS_H

#ifdef __cplusplus
extern "C" {
#endif

static PyObject * GMPy_MPZ_bit_length_function(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_bit_length_method(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_bit_mask(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_bit_scan0_function(PyObject *self, PyObject *const *args, Py_ssize_t nargs);
static PyObject * GMPy_MPZ_bit_scan0_method(PyObject *self, PyObject *const *args, Py_ssize_t nargs);
static PyObject * GMPy_MPZ_bit_scan1_function(PyObject *self, PyObject *const *args, Py_ssize_t nargs);
static PyObject * GMPy_MPZ_bit_scan1_method(PyObject *self, PyObject *const *args, Py_ssize_t nargs);

static PyObject * GMPy_MPZ_bit_test_function(PyObject *self, PyObject *const *args, Py_ssize_t nargs);
static PyObject * GMPy_MPZ_bit_test_method(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_bit_clear_function(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_bit_clear_method(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_bit_set_function(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_bit_set_method(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_bit_flip_function(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_bit_flip_method(PyObject *self, PyObject *other);

static PyObject * GMPy_MPZ_popcount(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_hamdist(PyObject *self, PyObject *args);

static PyObject * GMPy_MPZ_Invert_Slot(MPZ_Object *self);
static PyObject * GMPy_MPZ_And_Slot(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Ior_Slot(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Xor_Slot(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Rshift_Slot(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Lshift_Slot(PyObject *self, PyObject *other);

#ifdef __cplusplus
}
#endif
#endif
