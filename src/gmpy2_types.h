/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_types.h                                                           *
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

#ifndef GMPY_TYPES_H
#define GMPY_TYPES_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    PyObject_HEAD
    mpz_t z;
    Py_hash_t hash_cache;
} MPZ_Object;

typedef struct {
    PyObject_HEAD
    mpz_t z;
} XMPZ_Object;

typedef struct {
    PyObject_HEAD
    mpq_t q;
    Py_hash_t  hash_cache;
} MPQ_Object;

typedef struct {
    PyObject_HEAD
    mpfr_t f;
    Py_hash_t hash_cache;
    int rc;
} MPFR_Object;

typedef struct {
    PyObject_HEAD
    mpc_t c;
    Py_hash_t hash_cache;
    int rc;
} MPC_Object;

typedef struct {
    PyObject_HEAD
    gmp_randstate_t state;
} RandomState_Object;

typedef struct {
    mpfr_prec_t mpfr_prec;   /* current precision in bits, for MPFR */
    mpfr_rnd_t mpfr_round;   /* current rounding mode for float (MPFR) */
    mpfr_exp_t emax;         /* maximum exponent */
    mpfr_exp_t emin;         /* minimum exponent */
    int subnormalize;        /* if 1, subnormalization is performed */
    int underflow;           /* did an underflow occur? */
    int overflow;            /* did an overflow occur? */
    int inexact;             /* was the result inexact? */
    int invalid;             /* invalid operation (i.e. NaN)? */
    int erange;              /* did a range error occur? */
    int divzero;             /* divided by zero? */
    int traps;               /* if 0, do not trap any exceptions */
                             /* if not 0, then raise traps per bits above  */
    mpfr_prec_t real_prec;   /* current precision in bits, for Re(MPC) */
    mpfr_prec_t imag_prec;   /* current precision in bits, for Im(MPC) */
    mpfr_rnd_t real_round;   /* current rounding mode for Re(MPC) */
    mpfr_rnd_t imag_round;   /* current rounding mode for Im(MPC) */
    int allow_complex;       /* if 1, allow mpfr functions to return an mpc */
    int rational_division;   /* if 1, mpz/mpz returns an mpq result */
    int allow_release_gil;   /* if 1, release GIL for mpz operations */
} gmpy_context;

typedef struct {
    PyObject_HEAD
    gmpy_context ctx;
} CTXT_Object;


#ifdef __cplusplus
}
#endif

#endif
