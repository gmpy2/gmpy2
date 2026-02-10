/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_convert.h                                                          *
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

#ifndef GMPY2_CONVERT_H
#define GMPY2_CONVERT_H

#ifdef __cplusplus
extern "C" {
#endif

/* The following macros classify the numeric types that are supported by
 * gmpy2.
 */
#define HAS_MPZ_CONVERSION(x) (PyObject_HasAttrString(x, "__mpz__"))
#define HAS_MPQ_CONVERSION(x) (PyObject_HasAttrString(x, "__mpq__"))
#define HAS_MPFR_CONVERSION(x) (PyObject_HasAttrString(x, "__mpfr__"))
#define HAS_MPC_CONVERSION(x) (PyObject_HasAttrString(x, "__mpc__"))

#define HAS_STRICT_MPZ_CONVERSION(x) (HAS_MPZ_CONVERSION(x) && \
                                     !HAS_MPQ_CONVERSION(x))
#define HAS_STRICT_MPFR_CONVERSION(x) (HAS_MPFR_CONVERSION(x) && \
                                      !HAS_MPC_CONVERSION(x))

#define IS_FRACTION(x) (!strcmp(Py_TYPE(x)->tp_name, "Fraction"))

#define IS_RATIONAL_ONLY(x) (MPQ_Check(x) || IS_FRACTION(x) || \
                             HAS_MPQ_CONVERSION(x))

#define IS_INTEGER(x) (MPZ_Check(x) || PyLong_Check(x) || \
                       XMPZ_Check(x) || HAS_STRICT_MPZ_CONVERSION(x))
#define IS_RATIONAL(x) (MPQ_Check(x) || IS_FRACTION(x) || \
                        MPZ_Check(x) || PyLong_Check(x) || \
                        XMPZ_Check(x) || HAS_MPQ_CONVERSION(x) || \
                        HAS_MPZ_CONVERSION(x))
#define IS_DECIMAL(x) (!strcmp(Py_TYPE(x)->tp_name, "decimal.Decimal") || \
                       !strcmp(Py_TYPE(x)->tp_name, "Decimal"))
#define IS_REAL_ONLY(x) (MPFR_Check(x) || PyFloat_Check(x) || \
                         HAS_STRICT_MPFR_CONVERSION(x) || IS_DECIMAL(x))
#define IS_REAL(x) (IS_RATIONAL(x) || IS_REAL_ONLY(x))

#define IS_COMPLEX_ONLY(x) (MPC_Check(x) || PyComplex_Check(x) || \
                            HAS_MPC_CONVERSION(x))
#define IS_COMPLEX(x) (IS_REAL(x) || IS_COMPLEX_ONLY(x))

/* Define constants used in gmpy2_convert.c->GMPy_ObjectType. */

#define OBJ_TYPE_UNKNOWN        0
#define OBJ_TYPE_MPZ            1
#define OBJ_TYPE_XMPZ           2
#define OBJ_TYPE_PyInteger      3
#define OBJ_TYPE_HAS_MPZ        4
/* 5 TO 14 reserved for additional integer types. */
#define OBJ_TYPE_INTEGER        15

#define OBJ_TYPE_MPQ            16
#define OBJ_TYPE_PyFraction     17
#define OBJ_TYPE_HAS_MPQ        18
/* 19 to 30 reserved for additional rational types. */
#define OBJ_TYPE_RATIONAL       31

#define OBJ_TYPE_MPFR           32
#define OBJ_TYPE_PyFloat        33
#define OBJ_TYPE_HAS_MPFR       34
/* 35 to 46 reserved for additional real types. */
#define OBJ_TYPE_REAL           47

#define OBJ_TYPE_MPC            48
#define OBJ_TYPE_PyComplex      49
#define OBJ_TYPE_HAS_MPC        50
/* 50 to 62 reserved for additional complex types. */
#define OBJ_TYPE_COMPLEX        63

#define OBJ_TYPE_MAX            64

/* The following macros are the recommended method to check the result of the
 * object type check.
 */

#define IS_TYPE_UNKNOWN(x)          (!OBJ_TYPE_UNKNOWN)
#define IS_TYPE_MPZ(x)              (x == OBJ_TYPE_MPZ)
#define IS_TYPE_XMPZ(x)             (x == OBJ_TYPE_XMPZ)
#define IS_TYPE_MPZANY(x)           ((x == OBJ_TYPE_MPZ) || \
                                     (x == OBJ_TYPE_XMPZ))
#define IS_TYPE_PyInteger(x)        (x == OBJ_TYPE_PyInteger)
#define IS_TYPE_HAS_MPZ(x)          (x == OBJ_TYPE_HAS_MPZ)
#define IS_TYPE_INTEGER(x)          ((x > OBJ_TYPE_UNKNOWN) &&  \
                                     (x < OBJ_TYPE_INTEGER))

#define IS_TYPE_MPQ(x)              (x == OBJ_TYPE_MPQ)
#define IS_TYPE_PyFraction(x)       (x == OBJ_TYPE_PyFraction)
#define IS_TYPE_HAS_MPQ(x)          (x == OBJ_TYPE_HAS_MPQ)
#define IS_TYPE_RATIONAL(x)         ((x > OBJ_TYPE_UNKNOWN) && \
                                     (x < OBJ_TYPE_RATIONAL))
#define IS_TYPE_RATIONAL_ONLY(x)    ((x > OBJ_TYPE_INTEGER) && \
                                     (x < OBJ_TYPE_RATIONAL))

#define IS_TYPE_MPFR(x)             (x == OBJ_TYPE_MPFR)
#define IS_TYPE_PyFloat(x)          (x == OBJ_TYPE_PyFloat)
#define IS_TYPE_HAS_MPFR(x)         (x == OBJ_TYPE_HAS_MPFR)
#define IS_TYPE_REAL(x)             ((x > OBJ_TYPE_UNKNOWN) && \
                                     (x < OBJ_TYPE_REAL))
#define IS_TYPE_REAL_ONLY(x)        ((x > OBJ_TYPE_RATIONAL) && \
                                     (x < OBJ_TYPE_REAL))

#define IS_TYPE_MPC(x)              (x == OBJ_TYPE_MPC)
#define IS_TYPE_PyComplex(x)        (x == OBJ_TYPE_PyComplex)
#define IS_TYPE_HAS_MPC(x)          (x == OBJ_TYPE_HAS_MPC)
#define IS_TYPE_COMPLEX(x)          ((x > OBJ_TYPE_UNKNOWN) && \
                                     (x < OBJ_TYPE_COMPLEX))
#define IS_TYPE_COMPLEX_ONLY(x)     ((x > OBJ_TYPE_REAL) && \
                                     (x < OBJ_TYPE_COMPLEX))

/* Since the macros are used in gmpy2's codebase, these functions are skipped
 * until they are needed for the C API in the future.
 */

#if 0
/* Checks for mpz, xmpz, and the integer types included with Python. */
static int GMPy_isInteger(PyObject *obj);

/* Checks for the Fraction type included with Python. */
static int GMPy_isFraction(PyObject *obj);

/* Combined mpq, isInteger() and isFraction() check. */
static int GMPy_isRational(PyObject *obj);

/* Combined mpfr, PyFloat, and isRational() check. */
static int GMPy_isReal(PyObject *obj);

/* Combined mpc, PyComplex, and isReal() check. */
static int GMPy_isComplex(PyObject *obj);
#endif

/* ======== C helper routines ======== */
static int             mpz_set_PyStr(mpz_t z, PyObject *s, int base);
static PyObject *      mpz_ascii(mpz_t z, int base, int option, int which);

#ifdef __cplusplus
}
#endif
#endif
