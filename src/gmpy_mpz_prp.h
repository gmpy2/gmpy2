/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz_prp.h                                                          *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP, MPFR, and MPC multiple precision           *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2012 - 2025 Case Van Horsen                                   *
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

#ifndef GMPY_PRP_H
#define GMPY_PRP_H

#ifdef __cplusplus
extern "C" {
#endif

static PyObject * GMPY_mpz_is_fermat_prp(PyObject *self, PyObject *args);
static PyObject * GMPY_mpz_is_euler_prp(PyObject *self, PyObject *args);
static PyObject * GMPY_mpz_is_strong_prp(PyObject *self, PyObject *args);
static PyObject * GMPY_mpz_is_fibonacci_prp(PyObject *self, PyObject *args);
static PyObject * GMPY_mpz_is_lucas_prp(PyObject *self, PyObject *args);
static PyObject * GMPY_mpz_is_stronglucas_prp(PyObject *self, PyObject *args);
static PyObject * GMPY_mpz_is_extrastronglucas_prp(PyObject *self, PyObject *args);
static PyObject * GMPY_mpz_is_selfridge_prp(PyObject *self, PyObject *args);
static PyObject * GMPY_mpz_is_strongselfridge_prp(PyObject *self, PyObject *args);
static PyObject * GMPY_mpz_is_bpsw_prp(PyObject *self, PyObject *args);
static PyObject * GMPY_mpz_is_strongbpsw_prp(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif
