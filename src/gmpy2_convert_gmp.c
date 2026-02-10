/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_convert_gmp.c                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP, MPFR, and MPC multiple precision           *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000 - 2009 Alex Martelli                                     *
 *                                                                         *
 * Copyright 2008 - 2025 Case Van Horsen                                   *
 *                                                                         *
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

/* This file contains all the conversion functions for GMP data types.
 *
 * Overview
 * --------
 * gmpy2 tries to optimize the performance and accuracy of conversions from
 * other numeric types. gmpy2 uses a LBYL (Look Before You Leap) approach and
 * identifies the numeric type before conversion before conversion to a gmpy2
 * type. The basic operations (+, -, *, /) are optimized to directly work with
 * some basic types such as C longs or doubles.
 */

#include "pythoncapi_compat.h"

/* ======================================================================== *
 * Conversion between native Python objects and MPZ.                        *
 * ======================================================================== */

/* To support creation of temporary mpz objects. */
static int
mpz_set_PyLong(mpz_t z, PyObject *obj)
{
#ifndef PYPY_VERSION
    const PyLongLayout *layout = PyLong_GetNativeLayout();
    PyLongExport long_export = {0, 0, 0, 0, 0};

    if (PyLong_Export(obj, &long_export) < 0) {
        /* LCOV_EXCL_START */
        return -1;
        /* LCOV_EXCL_STOP */
    }
    if (long_export.digits) {
        mpz_import(z, long_export.ndigits, layout->digits_order,
                   layout->digit_size, layout->digit_endianness,
                   layout->digit_size*8 - layout->bits_per_digit,
                   long_export.digits);
        if (long_export.negative) {
            mpz_neg(z, z);
        }
        PyLong_FreeExport(&long_export);
    }
    else {
        const int64_t value = long_export.value;

        if (LONG_MIN <= value && value <= LONG_MAX) {
            mpz_set_si(z, value);
        }
        else {
            mpz_import(z, 1, -1, sizeof(int64_t), 0, 0, &value);
            if (value < 0) {
                mpz_t tmp;
                mpz_init(tmp);
                mpz_ui_pow_ui(tmp, 2, 64);
                mpz_sub(z, z, tmp);
                mpz_clear(tmp);
            }
        }
    }
    return 0;
#else
    int overflow;
    long value = PyLong_AsLongAndOverflow(obj, &overflow);
    if (!overflow) {
        mpz_set_si(z, value);
        return 0;
    }

    PyObject *s = PyNumber_ToBase(obj, 16);

    if (!s) {
        /* LCOV_EXCL_START */
        return -1;
        /* LCOV_EXCL_STOP */
    }

    const char *str = PyUnicode_AsUTF8(s), *p = str;

    if (!str) {
        /* LCOV_EXCL_START */
        Py_DECREF(s);
        return -1;
        /* LCOV_EXCL_STOP */
    }

    int negative = (str[0] == '-');

    p += 2;
    if (negative) {
        p++;
    }
    mpz_init_set_str(z, p, 16);
    Py_DECREF(s);
    if (negative) {
        mpz_neg(z, z);
    }
    return 0;
#endif
}

static MPZ_Object *
GMPy_MPZ_From_PyLong(PyObject *obj, CTXT_Object *context)
{
    MPZ_Object *result;

    if(!(result = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (mpz_set_PyLong(MPZ(result), obj)) {
        /* LCOV_EXCL_START */
        Py_DECREF((PyObject*)result);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    return result;
}

static MPZ_Object *
GMPy_MPZ_From_PyStr(PyObject *s, int base, CTXT_Object *context)
{
    MPZ_Object *result;

    if (!(result = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (mpz_set_PyStr(result->z, s, base) == -1) {
        Py_DECREF((PyObject*)result);
        return NULL;
    }
    return result;
}

static MPZ_Object *
GMPy_MPZ_From_PyFloat(PyObject *obj, CTXT_Object *context)
{
    MPZ_Object *result;

    if ((result = GMPy_MPZ_New(context))) {
        double d = PyFloat_AsDouble(obj);

        if (isnan(d)) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'mpz' does not support NaN");
            return NULL;
        }
        if (isinf(d)) {
            Py_DECREF((PyObject*)result);
            OVERFLOW_ERROR("'mpz' does not support Infinity");
            return NULL;
        }
        mpz_set_d(result->z, d);
    }
    return result;
}

static PyObject *
GMPy_PyLong_From_MPZ(MPZ_Object *obj, CTXT_Object *context)
{
    if (mpz_fits_slong_p(obj->z)) {
        return PyLong_FromLong(mpz_get_si(obj->z));
    }

#ifndef PYPY_VERSION
    const PyLongLayout *layout = PyLong_GetNativeLayout();
    size_t size = (mpz_sizeinbase(obj->z, 2) +
                   layout->bits_per_digit - 1)/layout->bits_per_digit;
    void *digits;
    PyLongWriter *writer = PyLongWriter_Create(mpz_sgn(obj->z) < 0, size,
                                               &digits);
    if (writer == NULL) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_export(digits, NULL, layout->digits_order, layout->digit_size,
               layout->digit_endianness,
               layout->digit_size*8 - layout->bits_per_digit, obj->z);
    return PyLongWriter_Finish(writer);
#else
    PyObject *str = GMPy_PyStr_From_MPZ(obj, 16, 0, NULL);

    if (!str) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    PyObject *res = PyLong_FromUnicodeObject(str, 16);

    Py_DECREF(str);
    return res;
#endif
}

static PyObject *
GMPy_MPZ_Int_Slot(MPZ_Object *self)
{
    return GMPy_PyLong_From_MPZ(self, NULL);
}

static CTXT_Object *
_get_ieee_context(long bitwidth)
{
    assert(bitwidth != 16 || bitwidth%32 == 0);

    CTXT_Object *result = (CTXT_Object*)GMPy_CTXT_New();

    if (!result) {
        return NULL;
    }
    if (bitwidth == 16) {
        result->ctx.mpfr_prec = 11;
        result->ctx.emax = 16;
    }
    else if (bitwidth == 32) {
        result->ctx.mpfr_prec = 24;
        result->ctx.emax = 128;
    }
    else if (bitwidth == 64) {
        result->ctx.mpfr_prec = 53;
        result->ctx.emax = 1024;
    }
    else if (bitwidth == 128) {
        result->ctx.mpfr_prec = 113;
        result->ctx.emax = 16384;
    }
    else {
        double bitlog2 = floor((4 * log(bitwidth) / log(2.0)) + 0.5);
        result->ctx.mpfr_prec = bitwidth - (long)(bitlog2) + 13;
        result->ctx.emax = 1 << (bitwidth - result->ctx.mpfr_prec - 1);
    }
    result->ctx.subnormalize = 1;
    result->ctx.emin = 4 - result->ctx.emax - result->ctx.mpfr_prec;
    return result;
}

static PyObject *
GMPy_PyFloat_From_MPZ(MPZ_Object *obj, CTXT_Object *unused)
{
    MPFR_Object *tmp;
    CTXT_Object *context = _get_ieee_context(64);

    CHECK_CONTEXT(context);

    if (!(tmp = GMPy_MPFR_New(53, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_clear_flags();
    tmp->rc = mpfr_set_z(tmp->f, obj->z, MPFR_RNDN);
    GMPY_MPFR_CHECK_RANGE(tmp, context);
    GMPY_MPFR_SUBNORMALIZE(tmp, context);
    GMPY_MPFR_EXCEPTIONS(tmp, context);

    PyObject *res = GMPy_PyFloat_From_MPFR(tmp, context);
    Py_DECREF(tmp);
    Py_DECREF(context);
    if (isinf(PyFloat_AS_DOUBLE(res))) {
        OVERFLOW_ERROR("'mpz' too large to convert to float");
        Py_CLEAR(res);
    }
    return res;
}

static PyObject *
GMPy_MPZ_Float_Slot(MPZ_Object *self)
{
    return GMPy_PyFloat_From_MPZ(self, NULL);
}

static PyObject *
GMPy_PyStr_From_MPZ(MPZ_Object *obj, int base, int option, CTXT_Object *context)
{
    return mpz_ascii(obj->z, base, option, 0);
}

static MPZ_Object *
GMPy_MPZ_From_Integer(PyObject *obj, CTXT_Object *context)
{
    MPZ_Object *result = NULL;

    if (MPZ_Check(obj)) {
        Py_INCREF(obj);
        return (MPZ_Object*)obj;
    }

    if (PyLong_Check(obj))
        return GMPy_MPZ_From_PyLong(obj, context);

    if (XMPZ_Check(obj))
        return GMPy_MPZ_From_XMPZ((XMPZ_Object*)obj, context);

    if (HAS_STRICT_MPZ_CONVERSION(obj)) {
        result = (MPZ_Object *) PyObject_CallMethod(obj, "__mpz__", NULL);

        if (result != NULL && MPZ_Check(result)) {
            return result;
        }
        else {
            Py_XDECREF((PyObject*)result);
            goto error;
        }
    }

  error:
    TYPE_ERROR("cannot convert object to mpz");
    return NULL;
}

static MPZ_Object *
GMPy_MPZ_From_IntegerAndCopy(PyObject *obj, CTXT_Object *context)
{
    MPZ_Object *result = NULL;

    if (MPZ_Check(obj)) {
        if (!(result = GMPy_MPZ_New(context))) {
            /* LCOV_EXCL_START */
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        mpz_set(result->z, MPZ(obj));
        return result;
    }

    if (PyLong_Check(obj))
        return GMPy_MPZ_From_PyLong(obj, context);

    if (XMPZ_Check(obj))
        return GMPy_MPZ_From_XMPZ((XMPZ_Object*)obj, context);

    if (HAS_STRICT_MPZ_CONVERSION(obj)) {
        result = (MPZ_Object *) PyObject_CallMethod(obj, "__mpz__", NULL);

        if (result != NULL && MPZ_Check(result)) {
            return result;
        }
        else {
            Py_XDECREF((PyObject*)result);
            goto error;
        }
    }

  error:
    TYPE_ERROR("cannot convert object to mpz");
    return NULL;
}

static MPZ_Object *
GMPy_MPZ_From_IntegerWithType(PyObject *obj, int xtype, CTXT_Object *context)
{
    MPZ_Object *result = NULL;

    if (IS_TYPE_MPZ(xtype)) {
        Py_INCREF(obj);
        return (MPZ_Object*)obj;
    }

    if (IS_TYPE_PyInteger(xtype))
        return GMPy_MPZ_From_PyLong(obj, context);

    if (IS_TYPE_XMPZ(xtype))
        return GMPy_MPZ_From_XMPZ((XMPZ_Object*)obj, context);

    if (IS_TYPE_HAS_MPZ(xtype)) {
        result = (MPZ_Object *) PyObject_CallMethod(obj, "__mpz__", NULL);

        if (result != NULL && MPZ_Check(result)) {
            return result;
        }
        else {
            Py_XDECREF((PyObject*)result);
            goto error;
        }
    }

  error:
    TYPE_ERROR("cannot convert object to mpz");
    return NULL;
}

static MPZ_Object *
GMPy_MPZ_From_IntegerWithTypeAndCopy(PyObject *obj, int xtype, CTXT_Object *context)
{
    MPZ_Object *result = NULL, *temp = NULL;

    result = GMPy_MPZ_From_IntegerWithType(obj, xtype, context);

    if (result == NULL)
        return result;

    if (Py_REFCNT(result) == 1)
        return result;

    if (!(temp = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_set(temp->z, result->z);
    Py_DECREF((PyObject*)result);
    return temp;
}

/* str and repr implementations for mpz */
static PyObject *
GMPy_MPZ_Str_Slot(MPZ_Object *self)
{
    /* base-10, no tag */
    return GMPy_PyStr_From_MPZ(self, 10, 0, NULL);
}

static PyObject *
GMPy_MPZ_Repr_Slot(MPZ_Object *self)
{
    /* base-10, with tag */
    return GMPy_PyStr_From_MPZ(self, 10, 1, NULL);
}

#ifdef SHARED
/* Helper function for argument parsing. Not used in static build. */

static int
GMPy_MPZ_ConvertArg(PyObject *arg, PyObject **ptr)
{
    MPZ_Object *result = GMPy_MPZ_From_IntegerWithType(arg, GMPy_ObjectType(arg), NULL);

    if (result) {
        *ptr = (PyObject*)result;
        return 1;
    }
    else {
        TYPE_ERROR("argument can not be converted to 'mpz'");
        return 0;
    }
}
#endif

/* ======================================================================== *
 * Conversion between native Python objects/MPZ and XMPZ.                   *
 * ======================================================================== */

static XMPZ_Object *
GMPy_XMPZ_From_PyLong(PyObject *obj, CTXT_Object *context)
{
    XMPZ_Object *result;

    if(!(result = GMPy_XMPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (mpz_set_PyLong(result->z, obj)) {
        /* LCOV_EXCL_START */
        Py_DECREF((PyObject*)result);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    return result;
}

static XMPZ_Object *
GMPy_XMPZ_From_PyStr(PyObject *s, int base, CTXT_Object *context)
{
    XMPZ_Object *result;

    if (!(result = GMPy_XMPZ_New(context)))
        return NULL;

    if (mpz_set_PyStr(result->z, s, base) == -1) {
        Py_DECREF((PyObject*)result);
        return NULL;
    }
    return result;
}

static XMPZ_Object *
GMPy_XMPZ_From_PyFloat(PyObject *obj, CTXT_Object *context)
{
    XMPZ_Object *result;

    if ((result = GMPy_XMPZ_New(context))) {
        double d = PyFloat_AsDouble(obj);

        if (isnan(d)) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'xmpz' does not support NaN");
            return NULL;
        }
        if (isinf(d)) {
            Py_DECREF((PyObject*)result);
            OVERFLOW_ERROR("'xmpz' does not support Infinity");
            return NULL;
        }
        mpz_set_d(result->z, d);
    }
    return result;
}

static XMPZ_Object *
GMPy_XMPZ_From_MPZ(MPZ_Object *obj, CTXT_Object *context)
{
    XMPZ_Object *result;

    if ((result = GMPy_XMPZ_New(context)))
        mpz_set(result->z, obj->z);

    return result;
}

static XMPZ_Object *
GMPy_XMPZ_From_XMPZ(XMPZ_Object *obj, CTXT_Object *context)
{
    XMPZ_Object *result;

    if ((result = GMPy_XMPZ_New(context)))
        mpz_set(result->z, obj->z);

    return result;
}

static PyObject *
GMPy_PyStr_From_XMPZ(XMPZ_Object *obj, int base, int option, CTXT_Object *context)
{
    return mpz_ascii(obj->z, base, option, 1);
}

static MPZ_Object *
GMPy_MPZ_From_XMPZ(XMPZ_Object *obj, CTXT_Object *context)
{
    MPZ_Object *result;

    if ((result = GMPy_MPZ_New(context)))
        mpz_set(result->z, obj->z);

    return result;
}

/* str and repr implementations for xmpz */
static PyObject *
GMPy_XMPZ_Str_Slot(XMPZ_Object *self)
{
    /* base-10, no tag */
    return GMPy_PyStr_From_XMPZ(self, 10, 0, NULL);
}

static PyObject *
GMPy_XMPZ_Repr_Slot(XMPZ_Object *self)
{
    /* base-10, with tag */
    return GMPy_PyStr_From_XMPZ(self, 10, 1, NULL);
}

/* ======================================================================== *
 * Conversion between native Python objects/MPZ/XMPZ and MPQ.               *
 * ======================================================================== */

static MPQ_Object *
GMPy_MPQ_From_PyLong(PyObject *obj, CTXT_Object *context)
{
    MPQ_Object *result;
    MPZ_Object *temp;

    temp = GMPy_MPZ_From_PyLong(obj, context);

    if (!temp)
        return NULL;

    if ((result = GMPy_MPQ_New(context))) {
        mpq_set_z(result->q, temp->z);
        Py_DECREF((PyObject*)temp);
    }
    return result;
}

static MPQ_Object *
GMPy_MPQ_From_PyStr(PyObject *s, int base, CTXT_Object *context)
{
    MPQ_Object *result;
    char *cp;
    char exp_char = 'E';
    long expt = 0;
    PyObject *ascii_str = ascii_str = GMPy_RemoveIgnoredASCII(s);

    if (!ascii_str) return NULL;

    if (!(result = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    cp = PyBytes_AsString(ascii_str);

    {
        char *whereslash = strchr((char*)cp, '/');
        char *wheredot = strchr((char*)cp, '.');
        char *whereexp = strchr((char*)cp, 'E');

        if (!whereexp) {
            whereexp = strchr((char*)cp, 'e');
            exp_char = 'e';
        }

        if (whereslash && wheredot) {
            VALUE_ERROR("illegal string: both . and / found");
            goto error;
        }

        if (wheredot && (base != 10)) {
            VALUE_ERROR("illegal string: embedded . requires base=10");
            goto error;
        }

        /* If base=10, no slash is found, and an exponent symbol is found, then
         * assume we have decimal number in scientific format.
         */
        if (whereexp && !whereslash && (base == 10)) {
            /* Temporarily shorten the string and continue processing as
             * normal. We'll deal with the exponent later.
             */
            *whereexp = '\0';
            expt = atol(whereexp+1);
        }

        /* Handle the case where an embedded decimal point exists. An optional
         * exponent is also allowed. We take advantage of the fact that
         * mpz_set_str() ignores embedded space characters. We count the
         * number of digits after the decimal point and then replace the '.'
         * with a ' '. We've already inserted a NUL byte to terminate the
         * string in the case when an exponent was entered. The value for
         * the exponent has already been read.
         */

        if (wheredot) {
            unsigned char *counter;
            long digits = 0;

            counter = (unsigned char*)wheredot;
            digits = 0;
            *wheredot = ' ';
            while (*++counter != '\0') {
                if (isdigit(*counter))
                    digits++;
            }
            if (-1 == mpz_set_str(mpq_numref(result->q), (char*)cp, base)) {
                if (wheredot)
                    *wheredot = '.';
                /* Restore the exponent! */
                if (whereexp && !whereslash && (base == 10))
                    *whereexp = exp_char;
                VALUE_ERROR("invalid digits");
                goto error;
            }
            /* Process the exponent. */
            digits = expt - digits;
            if (digits < 0) {
                mpz_ui_pow_ui(mpq_denref(result->q), 10, (unsigned long)(-digits));
            }
            else {
                /* Use the denominator instead of a temporary variable. */
                mpz_ui_pow_ui(mpq_denref(result->q), 10, (unsigned long)(digits));
                mpz_mul(mpq_numref(result->q), mpq_numref(result->q), mpq_denref(result->q));
                mpz_set_ui(mpq_denref(result->q), 1);
            }
            mpq_canonicalize(result->q);

            /* Restore the decimal point. */
            *wheredot = '.';

            /* Restore the exponent! */
            if (whereexp && (base == 10))
                *whereexp = exp_char;

            goto finish;
        }

        if (whereslash)
            *whereslash = '\0';

        /* Read the numerator. */
        if (-1 == mpz_set_str(mpq_numref(result->q), (char*)cp, base)) {
            if (whereslash)
                *whereslash = '/';
            VALUE_ERROR("invalid digits");
            goto error;
        }

        /* If a slash was present, read the denominator. */
        if (whereslash) {
            *whereslash = '/';
            if (-1 == mpz_set_str(mpq_denref(result->q), whereslash+1, base)) {
                VALUE_ERROR("invalid digits");
                goto error;
            }
            if (0 == mpz_sgn(mpq_denref(result->q))) {
                ZERO_ERROR("zero denominator in mpq()");
                goto error;
            }
            mpq_canonicalize(result->q);
        }
        else {
            /* Since no slash was present, either the denominator is 1 or a
             * power of 10.
             */

            if (expt <= 0) {
                mpz_ui_pow_ui(mpq_denref(result->q), 10, (unsigned long)(-expt));
            }
            else {
                mpz_ui_pow_ui(mpq_denref(result->q), 10, (unsigned long)(expt));
                mpz_mul(mpq_numref(result->q), mpq_numref(result->q), mpq_denref(result->q));
                mpz_set_ui(mpq_denref(result->q), 1);
            }
            mpq_canonicalize(result->q);
            if (whereexp && (base == 10))
                *whereexp = exp_char;
        }
    }

  finish:
    Py_XDECREF(ascii_str);
    return result;

  error:
    Py_DECREF((PyObject*)result);
    Py_XDECREF(ascii_str);
    return NULL;
}

static MPQ_Object *
GMPy_MPQ_From_PyFloat(PyObject *obj, CTXT_Object *context)
{
    MPQ_Object *result;

    if ((result = GMPy_MPQ_New(context))) {
        double d = PyFloat_AsDouble(obj);

        if (isnan(d)) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'mpq' does not support NaN");
            return NULL;
        }
        if (isinf(d)) {
            Py_DECREF((PyObject*)result);
            OVERFLOW_ERROR("'mpq' does not support Infinity");
            return NULL;
        }
        mpq_set_d(result->q, d);
    }

    return result;
}

static MPQ_Object *
GMPy_MPQ_From_MPZ(MPZ_Object *obj, CTXT_Object *context)
{
    MPQ_Object *result;

    if ((result = GMPy_MPQ_New(context)))
        mpq_set_z(result->q, obj->z);

    return result;
}

static MPQ_Object *
GMPy_MPQ_From_XMPZ(XMPZ_Object *obj, CTXT_Object *context)
{
    MPQ_Object *result;

    if ((result = GMPy_MPQ_New(context)))
        mpq_set_z(result->q, obj->z);

    return result;
}

static MPZ_Object *
GMPy_MPZ_From_MPQ(MPQ_Object *obj, CTXT_Object *context)
{
    MPZ_Object *result;

    if ((result = GMPy_MPZ_New(context)))
        mpz_set_q(result->z, obj->q);

    return result;
}

static XMPZ_Object *
GMPy_XMPZ_From_MPQ(MPQ_Object *obj, CTXT_Object *context)
{
    XMPZ_Object *result;

    if ((result = GMPy_XMPZ_New(context)))
        mpz_set_q(result->z, obj->q);

    return result;
}

static PyObject *
GMPy_PyLong_From_MPQ(MPQ_Object *obj, CTXT_Object *context)
{
    PyObject *result;
    MPZ_Object *temp;

    temp = GMPy_MPZ_From_MPQ(obj, context);

    if (!temp) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = GMPy_PyLong_From_MPZ(temp, context);
    Py_DECREF((PyObject*)temp);

    return result;
}

static PyObject *
GMPy_MPQ_Int_Slot(MPQ_Object *self)
{
    return GMPy_PyLong_From_MPQ(self, NULL);
}

static char* _qtag = "mpq(";

static PyObject *
GMPy_PyStr_From_MPQ(MPQ_Object *obj, int base, int option, CTXT_Object *context)
{
    PyObject *result = NULL, *numstr = NULL, *denstr = NULL;
    char buffer[50], *p;

    numstr = mpz_ascii(mpq_numref(obj->q), base, 0, 0);
    if (!numstr) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    /* Check if denominator is 1 and no tag is requested. If so, just
     * return the numerator.
     */
    if (!(option & 1) && (0 == mpz_cmp_ui(mpq_denref(obj->q),1)))
        return numstr;

    denstr = mpz_ascii(mpq_denref(obj->q), base, 0, 0);
    if (!denstr) {
        /* LCOV_EXCL_START */
        Py_DECREF(numstr);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    /* Build the format string. */
    p = buffer;
    if (option & 1) {
        strcpy(p, _qtag);
        p += strlen(p);
    }

    *(p++) = '%';
    *(p++) = 'U';
    if (option & 1)
        *(p++) = ',';
    else
        *(p++) = '/';
    *(p++) = '%';
    *(p++) = 'U';
    if (option & 1)
        *(p++) = ')';
    *(p++) = '\00';
    result = PyUnicode_FromFormat(buffer, numstr, denstr);
    Py_DECREF(numstr);
    Py_DECREF(denstr);
    return result;
}

static PyObject *
GMPy_PyFloat_From_MPQ(MPQ_Object *obj, CTXT_Object *unused)
{
    MPFR_Object *tmp;
    CTXT_Object *context = _get_ieee_context(64);

    CHECK_CONTEXT(context);

    if (!(tmp = GMPy_MPFR_New(53, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_clear_flags();
    tmp->rc = mpfr_set_q(tmp->f, obj->q, MPFR_RNDN);
    GMPY_MPFR_CHECK_RANGE(tmp, context);
    GMPY_MPFR_SUBNORMALIZE(tmp, context);
    GMPY_MPFR_EXCEPTIONS(tmp, context);

    PyObject *res = GMPy_PyFloat_From_MPFR(tmp, context);
    Py_DECREF(tmp);
    Py_DECREF(context);
    if (isinf(PyFloat_AS_DOUBLE(res))) {
        OVERFLOW_ERROR("'mpq' too large to convert to float");
        Py_CLEAR(res);
    }
    return res;
}

static PyObject *
GMPy_MPQ_Float_Slot(MPQ_Object *self)
{
    return GMPy_PyFloat_From_MPQ(self, NULL);
}

static MPQ_Object*
GMPy_MPQ_From_Fraction(PyObject* obj, CTXT_Object *context)
{
    MPQ_Object *result;
    PyObject *num, *den;

    if (!(result = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }
    mpq_set_si(result->q, 0, 1);

    num = PyObject_GetAttrString(obj, "numerator");
    den = PyObject_GetAttrString(obj, "denominator");
    if (!num || !PyLong_Check(num) || !den || !PyLong_Check(den)) {
        SYSTEM_ERROR("Object does not appear to be Fraction");
        goto error;
    }
    if (mpz_set_PyLong(mpq_numref(result->q), num)) {
        /* LCOV_EXCL_START */
        goto error;
        /* LCOV_EXCL_STOP */
    }
    if (mpz_set_PyLong(mpq_denref(result->q), den)) {
        /* LCOV_EXCL_START */
        goto error;
        /* LCOV_EXCL_STOP */
    }
    Py_DECREF(num);
    Py_DECREF(den);
    return result;
error:
    Py_XDECREF(num);
    Py_XDECREF(den);
    Py_DECREF((PyObject*)result);
    return NULL;
}

static MPQ_Object*
GMPy_MPQ_From_Number(PyObject *obj, CTXT_Object *context)
{
    if (MPQ_Check(obj)) {
        Py_INCREF(obj);
        return (MPQ_Object*)obj;
    }

    if (MPZ_Check(obj))
        return GMPy_MPQ_From_MPZ((MPZ_Object*)obj, context);

    if (MPFR_Check(obj))
        return GMPy_MPQ_From_MPFR((MPFR_Object*)obj, context);

    if (PyFloat_Check(obj))
        return GMPy_MPQ_From_PyFloat(obj, context);

    if (PyLong_Check(obj))
        return GMPy_MPQ_From_PyLong(obj, context);

    if (XMPZ_Check(obj))
        return GMPy_MPQ_From_XMPZ((XMPZ_Object*)obj, context);

    if (IS_FRACTION(obj))
        return GMPy_MPQ_From_Fraction(obj, context);

    PyObject *pair = PyObject_CallMethod(obj, "as_integer_ratio", NULL);
    if (pair != NULL) {
         MPQ_Object *res = (MPQ_Object*)GMPy_MPQ_NewInit(&MPQ_Type, pair, NULL);
         Py_DECREF(pair);
         return res;
    }
    PyErr_Clear();

    if (HAS_MPQ_CONVERSION(obj)) {
        MPQ_Object * res = (MPQ_Object *) PyObject_CallMethod(obj, "__mpq__", NULL);

        if (res != NULL && MPQ_Check(res)) {
            return res;
        }
        else {
            Py_XDECREF((PyObject*)res);
            goto error;
        }
    }

    if (HAS_MPZ_CONVERSION(obj)) {
        MPZ_Object * res = (MPZ_Object *) PyObject_CallMethod(obj, "__mpz__", NULL);

        if (res != NULL && MPZ_Check(res)) {
            MPQ_Object * temp = GMPy_MPQ_From_MPZ(res, context);
            Py_DECREF(res);
            return temp;
        }
        else {
            Py_XDECREF((PyObject*)res);
            goto error;
        }
    }

  error:
    TYPE_ERROR("cannot convert object to mpq");
    return NULL;
}

static MPQ_Object*
GMPy_MPQ_From_NumberWithType(PyObject *obj, int xtype, CTXT_Object *context)
{
    if (IS_TYPE_MPQ(xtype)) {
        Py_INCREF(obj);
        return (MPQ_Object*)obj;
    }

    if (IS_TYPE_MPZ(xtype))
        return GMPy_MPQ_From_MPZ((MPZ_Object*)obj, context);

    if (IS_TYPE_MPFR(xtype))
        return GMPy_MPQ_From_MPFR((MPFR_Object*)obj, context);

    if (IS_TYPE_PyFloat(xtype))
        return GMPy_MPQ_From_PyFloat(obj, context);

    if (IS_TYPE_PyInteger(xtype))
        return GMPy_MPQ_From_PyLong(obj, context);

    if (IS_TYPE_XMPZ(xtype))
        return GMPy_MPQ_From_XMPZ((XMPZ_Object*)obj, context);

    if (IS_TYPE_PyFraction(xtype))
        return GMPy_MPQ_From_Fraction(obj, context);

    if (IS_TYPE_HAS_MPQ(xtype)) {
        MPQ_Object * res = (MPQ_Object *) PyObject_CallMethod(obj, "__mpq__", NULL);

        if (res != NULL && MPQ_Check(res)) {
            return res;
        }
        else {
            Py_XDECREF((PyObject*)res);
            goto error;
        }
    }

    if (IS_TYPE_HAS_MPZ(xtype)) {
        MPZ_Object * res = (MPZ_Object *) PyObject_CallMethod(obj, "__mpz__", NULL);

        if (res != NULL && MPZ_Check(res)) {
            MPQ_Object * temp = GMPy_MPQ_From_MPZ(res, context);
            Py_DECREF(res);
            return temp;
        }
        else {
            Py_XDECREF((PyObject*)res);
            goto error;
        }
    }

  error:
    TYPE_ERROR("cannot convert object to mpq");
    return NULL;
}

static MPQ_Object*
GMPy_MPQ_From_RationalWithTypeAndCopy(PyObject *obj, int xtype, CTXT_Object *context)
{
    MPQ_Object *result = NULL, *temp = NULL;

    result = GMPy_MPQ_From_RationalWithType(obj, xtype, context);

    if (result == NULL)
        return result;

    if (Py_REFCNT(result) == 1)
        return result;

    if (!(temp = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpq_set(temp->q, result->q);
    Py_DECREF((PyObject*)result);
    return temp;
}

static MPQ_Object*
GMPy_MPQ_From_Rational(PyObject *obj, CTXT_Object *context)
{
    if (MPQ_Check(obj)) {
        Py_INCREF(obj);
        return (MPQ_Object*)obj;
    }

    if (MPZ_Check(obj))
        return GMPy_MPQ_From_MPZ((MPZ_Object*)obj, context);

    if (PyLong_Check(obj))
        return GMPy_MPQ_From_PyLong(obj, context);

    if (XMPZ_Check(obj))
        return GMPy_MPQ_From_XMPZ((XMPZ_Object*)obj, context);

    if (IS_FRACTION(obj))
        return GMPy_MPQ_From_Fraction(obj, context);

    if (HAS_MPQ_CONVERSION(obj)) {
        MPQ_Object * res = (MPQ_Object *) PyObject_CallMethod(obj, "__mpq__", NULL);

        if (res != NULL && MPQ_Check(res)) {
            return res;
        }
        else {
            Py_XDECREF((PyObject*)res);
            goto error;
        }
    }

    if (HAS_MPZ_CONVERSION(obj)) {
        MPZ_Object * res = (MPZ_Object *) PyObject_CallMethod(obj, "__mpz__", NULL);

        if (res != NULL && MPZ_Check(res)) {
            MPQ_Object * temp = GMPy_MPQ_From_MPZ(res, context);
            Py_DECREF(res);
            return temp;
        }
        else {
            Py_XDECREF((PyObject*)res);
            goto error;
        }
    }

  error:
    TYPE_ERROR("cannot convert object to mpq");
    return NULL;
}

static MPQ_Object*
GMPy_MPQ_From_RationalAndCopy(PyObject *obj, CTXT_Object *context)
{
    MPQ_Object *result = NULL, *temp = NULL;

    result = GMPy_MPQ_From_Rational(obj, context);

    if (result == NULL)
        return result;

    if (Py_REFCNT(result) == 1)
        return result;

    if (!(temp = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpq_set(temp->q, result->q);
    Py_DECREF((PyObject*)result);
    return temp;
}


static MPQ_Object*
GMPy_MPQ_From_RationalWithType(PyObject *obj, int xtype, CTXT_Object *context)
{
    if (IS_TYPE_MPQ(xtype)) {
        Py_INCREF(obj);
        return (MPQ_Object*)obj;
    }

    if (IS_TYPE_MPZ(xtype))
        return GMPy_MPQ_From_MPZ((MPZ_Object*)obj, context);

    if (IS_TYPE_PyInteger(xtype))
        return GMPy_MPQ_From_PyLong(obj, context);

    if (IS_TYPE_XMPZ(xtype))
        return GMPy_MPQ_From_XMPZ((XMPZ_Object*)obj, context);

    if (IS_TYPE_PyFraction(xtype))
        return GMPy_MPQ_From_Fraction(obj, context);

    if (IS_TYPE_HAS_MPQ(xtype)) {
        MPQ_Object * res = (MPQ_Object *) PyObject_CallMethod(obj, "__mpq__", NULL);

        if (res != NULL && MPQ_Check(res)) {
            return res;
        }
        else {
            Py_XDECREF((PyObject*)res);
            goto error;
        }
    }

    if (IS_TYPE_HAS_MPZ(xtype)) {
        MPZ_Object * res = (MPZ_Object *) PyObject_CallMethod(obj, "__mpz__", NULL);

        if (res != NULL && MPZ_Check(res)) {
            MPQ_Object * temp = GMPy_MPQ_From_MPZ(res, context);
            Py_DECREF(res);
            return temp;
        }
        else {
            Py_XDECREF((PyObject*)res);
            goto error;
        }
    }

  error:
    TYPE_ERROR("cannot convert object to mpq");
    return NULL;
}

/*
 * coerce any number to a mpq
 */

#ifdef SHARED
/* Helper function for argument parsing. Not used in static build. */

int
GMPy_MPQ_ConvertArg(PyObject *arg, PyObject **ptr)
{
    MPQ_Object* result = GMPy_MPQ_From_NumberWithType(arg, GMPy_ObjectType(arg), NULL);

    if (result) {
        *ptr = (PyObject*)result;
        return 1;
    }
    else {
        if (!PyErr_Occurred()) {
            TYPE_ERROR("argument can not be converted to 'mpq'");
        }
        return 0;
    }
}

#endif

/* str and repr implementations for mpq */
static PyObject *
GMPy_MPQ_Str_Slot(MPQ_Object *self)
{
    /* base-10, no tag */
    return GMPy_PyStr_From_MPQ(self, 10, 0, NULL);
}

static PyObject *
GMPy_MPQ_Repr_Slot(MPQ_Object *self)
{
    /* base-10, with tag */
    return GMPy_PyStr_From_MPQ(self, 10, 1, NULL);
}

