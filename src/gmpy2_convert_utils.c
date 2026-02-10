/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_convert_utils.c                                                   *
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

/* ======================================================================== *
 * Conversion between Integer objects and C types.                          *
 * ======================================================================== *
 *
 * Optimized routines for converting an Integer object (Python's integer
 * type(s), mpz, plus types defining __mpz__) to various C types.
 */


static long
GMPy_Integer_AsLongWithType(PyObject *x, int xtype)
{
    if IS_TYPE_PyInteger(xtype) {
        return PyLong_AsLong(x);
    }

    if (IS_TYPE_MPZANY(xtype)) {
        if (mpz_fits_slong_p(MPZ(x))) {
            return (long) mpz_get_si(MPZ(x));
        }
        else {
            OVERFLOW_ERROR("value could not be converted to C long");
            return -1;
        }
    }

    if (IS_TYPE_HAS_MPZ(xtype)) {
        long result = 0;
        PyObject *temp_mpz = PyObject_CallMethod(x, "__mpz__", NULL);

        if ((temp_mpz != NULL) && MPZ_Check(temp_mpz)) {
            if (mpz_fits_slong_p(MPZ(temp_mpz))) {
                result = (long) mpz_get_si(MPZ(temp_mpz));
            }
            else {
                OVERFLOW_ERROR("value could not be converted to C long");
                result = -1;
            }
        }
        Py_XDECREF(temp_mpz);
        return result;
    }
    TYPE_ERROR("could not convert object to integer");
    return (long)-1;
}

static long
GMPy_Integer_AsLong(PyObject *x)
{
    return GMPy_Integer_AsLongWithType(x, GMPy_ObjectType(x));
}

static unsigned long
GMPy_Integer_AsUnsignedLongWithType(PyObject *x, int xtype)
{
    if IS_TYPE_PyInteger(xtype) {
        return PyLong_AsUnsignedLong(x);
    }

    if (IS_TYPE_MPZANY(xtype)) {
        if (mpz_fits_ulong_p(MPZ(x))) {
            return (unsigned long) mpz_get_ui(MPZ(x));
        }
        else {
            OVERFLOW_ERROR("value could not be converted to C long");
            return (unsigned long)-1;
        }
    }

    if (IS_TYPE_HAS_MPZ(xtype)) {
        unsigned long result = 0;
        PyObject *temp_mpz = PyObject_CallMethod(x, "__mpz__", NULL);

        if ((temp_mpz != NULL) && MPZ_Check(temp_mpz)) {
            if (mpz_fits_ulong_p(MPZ(temp_mpz))) {
                result = (unsigned long) mpz_get_ui(MPZ(temp_mpz));
            }
            else {
                OVERFLOW_ERROR("value could not be converted to C long");
                result = (unsigned long)-1;
            }
        }
        Py_XDECREF(temp_mpz);
        return result;
    }
    TYPE_ERROR("could not convert object to integer");
    return (unsigned long)-1;
}

static unsigned long
GMPy_Integer_AsUnsignedLong(PyObject *x)
{
    return GMPy_Integer_AsUnsignedLongWithType(x, GMPy_ObjectType(x));
}

static long
GMPy_Integer_AsUnsignedLongOrLong(PyObject *x, int *is_signed)
{
    long result = (long)GMPy_Integer_AsUnsignedLong(x);
    if (result == -1 && PyErr_Occurred()) {
        *is_signed = 1;
        PyErr_Clear();
        result = GMPy_Integer_AsLong(x);
        if (result == -1 && PyErr_Occurred()) {
            return -1;
        }
    }
    return result;
}

#define PY_ABS_LLONG_MIN (0-(unsigned PY_LONG_LONG)PY_LLONG_MIN)

static PY_LONG_LONG
GMPy_Integer_AsLongLongWithType(PyObject *x, int xtype)
{
    int sign = 0;

    if IS_TYPE_PyInteger(xtype) {
        return PyLong_AsLongLong(x);
    }

    if (IS_TYPE_MPZANY(xtype)) {
        unsigned PY_LONG_LONG xtemp = 0;
        PY_LONG_LONG result = 0;
        sign = mpz_sgn(MPZ(x));
        if (sign) {
            if (mpz_sizeinbase(MPZ(x), 256) <= sizeof(xtemp)) {
                mpz_export(&xtemp, NULL, 1, sizeof(xtemp), 0, 0, MPZ(x));
            }
            if (xtemp <= (unsigned PY_LONG_LONG)PY_LLONG_MAX) {
                result = (PY_LONG_LONG)xtemp * sign;
            }
            else if (sign < 0 && xtemp == PY_ABS_LLONG_MIN) {
                result = PY_LLONG_MIN;
            }
            else {
                OVERFLOW_ERROR("value could not be converted to C long long");
                result = -1;
            }
        }
        return result;
    }

    if (IS_TYPE_HAS_MPZ(xtype)) {
        unsigned PY_LONG_LONG xtemp = 0;
        PY_LONG_LONG result = 0;
        PyObject *temp_mpz = PyObject_CallMethod(x, "__mpz__", NULL);

        if (temp_mpz != NULL && MPZ_Check(temp_mpz)) {
            sign = mpz_sgn(MPZ(temp_mpz));
            if (sign) {
                if (mpz_sizeinbase(MPZ(x), 256) <= sizeof(xtemp)) {
                    mpz_export(&xtemp, NULL, 1, sizeof(xtemp), 0, 0, MPZ(x));
                }
                if (xtemp <= (unsigned PY_LONG_LONG)PY_LLONG_MAX) {
                    result = (PY_LONG_LONG)xtemp * sign;
                }
                else if (sign < 0 && xtemp == PY_ABS_LLONG_MIN) {
                    result = PY_LLONG_MIN;
                }
                else {
                    OVERFLOW_ERROR("value could not be converted to C long long");
                    result = -1;
                }
            }
        }
        Py_XDECREF(temp_mpz);
        return result;
    }
    TYPE_ERROR("could not convert object to integer");
    return -1;
}
static PY_LONG_LONG
GMPy_Integer_AsLongLong(PyObject *x)
{
    return GMPy_Integer_AsLongLongWithType(x, GMPy_ObjectType(x));
}
static mp_bitcnt_t
GMPy_Integer_AsMpBitCnt(PyObject *x)
{
    if(sizeof(mp_bitcnt_t) == sizeof(unsigned long long))
        return GMPy_Integer_AsUnsignedLongLong(x);
    if(sizeof(mp_bitcnt_t) == sizeof(unsigned long))
        return GMPy_Integer_AsUnsignedLong(x);
    SYSTEM_ERROR("Unexpected sizes in GMPy_Integer_AsMpBitCnt.");
}

static mp_bitcnt_t
GMPy_PyLong_AsMpBitCnt(PyObject *x)
{
    if(sizeof(mp_bitcnt_t) == sizeof(unsigned long long))
        return PyLong_AsUnsignedLongLong(x);
    if(sizeof(mp_bitcnt_t) == sizeof(unsigned long))
        return PyLong_AsUnsignedLong(x);
    SYSTEM_ERROR("Unexpected sizes in GMPy_PyLong_AsMpBitCnt.");
}

static PyObject *
GMPy_PyLong_FromMpBitCnt(mp_bitcnt_t n)
{
    if(sizeof(mp_bitcnt_t) == sizeof(unsigned long long))
        return PyLong_FromUnsignedLongLong((unsigned long long)n);
    if(sizeof(mp_bitcnt_t) == sizeof(unsigned long))
        return PyLong_FromUnsignedLong((unsigned long)n);
    SYSTEM_ERROR("Unexpected sizes in GMPy_PyLong_FromMpBitCnt.");
}
static unsigned PY_LONG_LONG
GMPy_Integer_AsUnsignedLongLongWithType(PyObject *x, int xtype)
{
    if IS_TYPE_PyInteger(xtype) {
        return PyLong_AsUnsignedLongLong(x);
    }

    if (IS_TYPE_MPZANY(xtype)) {
        unsigned PY_LONG_LONG xtemp = 0;
        int sign = mpz_sgn(MPZ(x));

        if (sign >= 0) {
            if (mpz_sizeinbase(MPZ(x), 256) <= sizeof(xtemp)) {
                mpz_export(&xtemp, NULL, 1, sizeof(xtemp), 0, 0, MPZ(x));
            }
            else {
                OVERFLOW_ERROR("value could not be converted to C long long");
                xtemp = (unsigned PY_LONG_LONG)-1;
            }
        }
        return xtemp;
    }

    if (IS_TYPE_HAS_MPZ(xtype)) {
        unsigned PY_LONG_LONG xtemp = 0;
        PyObject *temp_mpz = PyObject_CallMethod(x, "__mpz__", NULL);

        if (temp_mpz != NULL && MPZ_Check(temp_mpz)) {
            int sign = mpz_sgn(MPZ(temp_mpz));

            if (sign >= 0) {
                if (mpz_sizeinbase(MPZ(x), 256) <= sizeof(xtemp)) {
                    mpz_export(&xtemp, NULL, 1, sizeof(xtemp), 0, 0, MPZ(x));
                }
                else {
                    OVERFLOW_ERROR("value could not be converted to C long long");
                    xtemp = (unsigned PY_LONG_LONG)-1;
                }
            }
        }
        Py_XDECREF((PyObject*)temp_mpz);
        return xtemp;
    }
    TYPE_ERROR("could not convert object to integer");
    return -1;
}

static unsigned PY_LONG_LONG
GMPy_Integer_AsUnsignedLongLong(PyObject *x)
{
    return GMPy_Integer_AsUnsignedLongLongWithType(x, GMPy_ObjectType(x));
}
