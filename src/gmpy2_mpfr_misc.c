/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpfr_misc.c                                                       *
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

PyDoc_STRVAR(GMPy_doc_function_f2q,
"f2q(x, err=0, /) -> mpz | mpq\n\n"
"Return the 'best' `mpq` approximating x to within relative error err.\n"
"Default is the precision of x. Uses Stern-Brocot tree to find the\n"
"'best' approximation. An `mpz` object is returned if the denominator\n"
"is 1. If err<0, relative error is 2.0 ** err.");

static PyObject *
GMPy_Real_F2Q(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *tempx, *tempy = NULL;
    PyObject *result;

    CHECK_CONTEXT(context);

    if (y) {
        if (!(tempy = GMPy_MPFR_From_Real(y, 1, context))) {
            /* LCOV_EXCL_START */
            return NULL;
            /* LCOV_EXCL_STOP */
        }
    }

    if (!(tempx = GMPy_MPFR_From_Real(x, 1, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)tempy);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    /* See gmpy2_convert_mpfr for stern_brocot(). */

    result = stern_brocot(tempx, tempy, 0, 1, context);
    Py_DECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    return result;
}

static PyObject *
GMPy_Number_F2Q(PyObject *x, PyObject *y, CTXT_Object *context)
{
    if (IS_REAL(x) && (!y || IS_REAL(y)))
        return GMPy_Real_F2Q(x, y, context);

    TYPE_ERROR("f2q() argument types not supported");
    return NULL;
}

static PyObject *
GMPy_Context_F2Q(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) < 1 || PyTuple_GET_SIZE(args) > 2) {
        TYPE_ERROR("f2q() requires 1 or 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        /* LCOV_EXCL_START */
        context = (CTXT_Object*)self;
        /* LCOV_EXCL_STOP */
    }
    else {
        CHECK_CONTEXT(context);
    }

    if (PyTuple_GET_SIZE(args) == 1) {
        return GMPy_Number_F2Q(PyTuple_GET_ITEM(args, 0), NULL, context);
    }
    else {
        return GMPy_Number_F2Q(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context);
    }
}

PyDoc_STRVAR(GMPy_doc_mpfr_free_cache,
"free_cache() -> None\n\n"
"Free the internal cache of constants maintained by MPFR.");

static PyObject *
GMPy_MPFR_Free_Cache(PyObject *self, PyObject *args)
{
    mpfr_free_cache();
    Py_RETURN_NONE;
}

PyDoc_STRVAR(GMPy_doc_mpfr_can_round,
"can_round(b, err, rnd1, rnd2, prec, /) -> bool\n\n"
"Let b be an approximation to an unknown number x that is rounded\n"
"according to rnd1. Assume the b has an error at most two to the power\n"
"of E(b)-err where E(b) is the exponent of b. Then return `True` if x\n"
"can be rounded correctly to prec bits with rounding mode rnd2.");

static PyObject *
GMPy_MPFR_Can_Round(PyObject *self, PyObject *args)
{
    int rnd1, rnd2;
    long err, prec;
    PyObject *b;

    if (!PyArg_ParseTuple(args, "O!liil", &MPFR_Type, &b, &err, &rnd1, &rnd2, &prec)) {
        return NULL;
    }

    if (!(rnd1 == MPFR_RNDN || rnd1 == MPFR_RNDZ ||
          rnd1 == MPFR_RNDU || rnd1 == MPFR_RNDD ||
          rnd1 == MPFR_RNDA)) {
        VALUE_ERROR("invalid value for rounding mode");
        return NULL;
    }

    if (!(rnd2 == MPFR_RNDN || rnd2 == MPFR_RNDZ ||
          rnd2 == MPFR_RNDU || rnd2 == MPFR_RNDD ||
          rnd2 == MPFR_RNDA)) {
        VALUE_ERROR("invalid value for rounding mode");
        return NULL;
    }

    if (prec < MPFR_PREC_MIN || prec > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }

    if (mpfr_can_round(MPFR(b), err, rnd1, rnd2, prec))
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpfr_get_emin_min,
"get_emin_min() -> int\n\n"
"Return the minimum possible exponent that can be set for `mpfr`.");

static PyObject *
GMPy_MPFR_get_emin_min(PyObject *self, PyObject *args)
{
    return PyLong_FromSsize_t((Py_ssize_t)mpfr_get_emin_min());
}

PyDoc_STRVAR(GMPy_doc_mpfr_get_emax_max,
"get_emax_max() -> int\n\n"
"Return the maximum possible exponent that can be set for `mpfr`.");

static PyObject *
GMPy_MPFR_get_emax_max(PyObject *self, PyObject *args)
{
    return PyLong_FromSsize_t((Py_ssize_t)mpfr_get_emax_max());
}

PyDoc_STRVAR(GMPy_doc_mpfr_get_max_precision,
"get_max_precision() -> int\n\n"
"Return the maximum bits of precision that can be used for calculations.\n"
"Note: to allow extra precision for intermediate calculations, avoid\n"
"setting precision close the maximum precision.");

static PyObject *
GMPy_MPFR_get_max_precision(PyObject *self, PyObject *args)
{
    return PyLong_FromSsize_t((Py_ssize_t)MPFR_PREC_MAX);
}

PyDoc_STRVAR(GMPy_doc_mpfr_get_exp,
"get_exp(x, /) -> int\n\n"
"Return the exponent of x. Returns 0 for NaN or Infinity and\n"
"sets the `context.erange` flag of the current context and will\n"
"raise an exception if `context.trap_erange` is set.");

static PyObject *
GMPy_MPFR_get_exp(PyObject *self, PyObject *other)
{
    PyObject *result = NULL;
    Py_ssize_t exp;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (!(MPFR_Check(other))) {
        TYPE_ERROR("get_exp() requires 'mpfr' argument");
        return NULL;
    }

    if (mpfr_regular_p(MPFR(other))) {
        exp = (Py_ssize_t)mpfr_get_exp(MPFR(other));
        result = PyLong_FromSsize_t(exp);
    }
    else if (mpfr_zero_p(MPFR(other))) {
        result = PyLong_FromSsize_t(0);
    }
    else {
        context->ctx.erange = 1;
        if (context->ctx.traps & TRAP_ERANGE) {
            GMPY_ERANGE("Can not get exponent from NaN or Infinity.");
        }
        else {
            result = PyLong_FromSsize_t(0);
        }
    }
    return result;
}

PyDoc_STRVAR(GMPy_doc_mpfr_set_exp,
"set_exp(x, n, /) -> mpfr\n\n"
"Set the exponent of x to n. If n is outside the range of\n"
"valid exponents, `set_exp()` will set the `context.erange` flag\n"
"of the current context and either return the original value or raise\n"
"an exception if `context.trap_erange` is set.");

static PyObject *
GMPy_MPFR_set_exp(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    PyObject *temp;
    mpfr_exp_t _oldemin, _oldemax, exp;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_GET_SIZE(args) != 2 ||
        !MPFR_Check(PyTuple_GET_ITEM(args, 0)) ||
        !PyLong_Check(PyTuple_GET_ITEM(args, 1))) {
        TYPE_ERROR("set_exp() requires 'mpfr', 'integer' arguments");
        return NULL;
    }

    temp = PyTuple_GET_ITEM(args, 0);
    exp = (mpfr_exp_t)PyLong_AsLong(PyTuple_GET_ITEM(args, 1));
    if (exp == -1 && PyErr_Occurred()) {
        VALUE_ERROR("exponent too large");
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(mpfr_get_prec(MPFR(temp)), context))) {
        return NULL;
    }

    _oldemin = mpfr_get_emin();
    _oldemax = mpfr_get_emax();
    mpfr_set_emin(context->ctx.emin);
    mpfr_set_emax(context->ctx.emax);

    mpfr_set(MPFR(result), MPFR(temp), GET_MPFR_ROUND(context));
    result->rc = mpfr_set_exp(MPFR(result), exp);

    mpfr_set_emin(_oldemin);
    mpfr_set_emax(_oldemax);

    if (result->rc) {
        context->ctx.erange = 1;
        if (context->ctx.traps & TRAP_ERANGE) {
            GMPY_ERANGE("new exponent is out-of-bounds");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpfr_set_sign,
"set_sign(x, s, /) -> mpfr\n\n"
"If s is `True`, then return x with the sign bit set.");

static PyObject *
GMPy_MPFR_set_sign(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_GET_SIZE(args) != 2 ||
        !MPFR_Check(PyTuple_GET_ITEM(args, 0)) ||
        !PyLong_Check(PyTuple_GET_ITEM(args, 1))) {
        TYPE_ERROR("set_sign() requires 'mpfr', 'boolean' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    result->rc = mpfr_setsign(MPFR(result), MPFR(PyTuple_GET_ITEM(args, 0)),
                              PyObject_IsTrue(PyTuple_GET_ITEM(args, 1)),
                              GET_MPFR_ROUND(context));

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpfr_copy_sign,
"copy_sign(x, y, /) -> mpfr\n\n"
"Return an `mpfr` composed of x with the sign of y.");

static PyObject *
GMPy_MPFR_copy_sign(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_GET_SIZE(args) != 2 ||
        !MPFR_Check(PyTuple_GET_ITEM(args, 0)) ||
        !MPFR_Check(PyTuple_GET_ITEM(args, 1))) {
        TYPE_ERROR("copy_sign() requires 'mpfr', 'mpfr' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    result->rc = mpfr_copysign(MPFR(result), MPFR(PyTuple_GET_ITEM(args, 0)),
                               MPFR(PyTuple_GET_ITEM(args, 1)),
                               GET_MPFR_ROUND(context));

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpfr_set_nan,
"nan() -> mpfr\n\n"
"Return an `mpfr` initialized to NaN (Not-A-Number).");

static PyObject *
GMPy_MPFR_set_nan(PyObject *self, PyObject *other)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if ((result = GMPy_MPFR_New(0, context))) {
        mpfr_set_nan(result->f);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpfr_set_inf,
"inf(n, /) -> mpfr\n\n"
"Return an `mpfr` initialized to Infinity with the same sign as n.\n"
"If n is not given, +Infinity is returned.");

static PyObject *
GMPy_MPFR_set_inf(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    long s = 1;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_Size(args) == 1) {
        s = GMPy_Integer_AsLongWithType(PyTuple_GET_ITEM(args, 0), 
                                        GMPy_ObjectType(PyTuple_GET_ITEM(args, 0)));
        if (s == -1 && PyErr_Occurred()) {
            return NULL;
        }
    }

    if ((result = GMPy_MPFR_New(0, context))) {
        mpfr_set_inf(result->f, s < 0 ? -1 : 1);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpfr_set_zero,
"zero(n, /) -> mpfr\n\n"
"Return an `mpfr` initialized to 0.0 with the same sign as n.\n"
"If n is not given, +0.0 is returned.");

static PyObject *
GMPy_MPFR_set_zero(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    long s = 1;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_Size(args) == 1) {
        s = GMPy_Integer_AsLongWithType(PyTuple_GET_ITEM(args, 0), 
                                        GMPy_ObjectType(PyTuple_GET_ITEM(args, 0)));
        if (s == -1 && PyErr_Occurred()) {
            return NULL;
        }
    }

    if ((result = GMPy_MPFR_New(0, context))) {
        mpfr_set_zero(result->f, s < 0 ? -1 : 1);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_method_integer_ratio,
"x.as_integer_ratio() -> tuple[mpz, mpz]\n\n"
"Return the exact rational equivalent of an `mpfr`. Value is a `tuple`\n"
"for compatibility with Python's `float.as_integer_ratio`.");

/* Note: almost identical code exists in gmpy2_convert_mpfr.c as the
 * function GMPy_MPQ_From_MPFR. They should be refactored.
 */

static PyObject *
GMPy_MPFR_Integer_Ratio_Method(PyObject *self, PyObject *args)
{
    MPZ_Object *num, *den;
    mpfr_exp_t temp, twocount;
    PyObject *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (mpfr_nan_p(MPFR(self))) {
        VALUE_ERROR("Cannot pass NaN to mpfr.as_integer_ratio.");
        return NULL;
    }

    if (mpfr_inf_p(MPFR(self))) {
        OVERFLOW_ERROR("Cannot pass Infinity to mpfr.as_integer_ratio.");
        return NULL;
    }

    num = GMPy_MPZ_New(context);
    den = GMPy_MPZ_New(context);
    if (!num || !den) {
        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)num);
        Py_XDECREF((PyObject*)den);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (mpfr_zero_p(MPFR(self))) {
        mpz_set_ui(num->z, 0);
        mpz_set_ui(den->z, 1);
    }
    else {
        temp = mpfr_get_z_2exp(num->z, MPFR(self));
        twocount = (mpfr_exp_t)mpz_scan1(num->z, 0);
        if (twocount) {
            temp += twocount;
            mpz_div_2exp(num->z, num->z, twocount);
        }
        mpz_set_ui(den->z, 1);
        if (temp > 0)
            mpz_mul_2exp(num->z, num->z, temp);
        else if (temp < 0)
            mpz_mul_2exp(den->z, den->z, -temp);
    }
    result = Py_BuildValue("(NN)", (PyObject*)num, (PyObject*)den);
    if (!result) {
        /* LCOV_EXCL_START */
        Py_DECREF((PyObject*)num);
        Py_DECREF((PyObject*)den);
        /* LCOV_EXCL_STOP */
    }
    return result;
}

PyDoc_STRVAR(GMPy_doc_method_mantissa_exp,
"x.as_mantissa_exp() -> tuple[mpz, mpz]\n\n"
"Return the mantissa and exponent of an `mpfr`.");

static PyObject *
GMPy_MPFR_Mantissa_Exp_Method(PyObject *self, PyObject *args)
{
    MPZ_Object *mantissa , *exponent;
    mpfr_exp_t temp;
    PyObject *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (mpfr_nan_p(MPFR(self))) {
        VALUE_ERROR("Cannot pass NaN to mpfr.as_mantissa_exp.");
        return NULL;
    }

    if (mpfr_inf_p(MPFR(self))) {
        OVERFLOW_ERROR("Cannot pass Infinity to mpfr.as_mantissa_exp.");
        return NULL;
    }

    mantissa = GMPy_MPZ_New(context);
    exponent = GMPy_MPZ_New(context);
    if (!mantissa || !exponent) {
        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)mantissa);
        Py_XDECREF((PyObject*)exponent);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (mpfr_zero_p(MPFR(self))) {
        mpz_set_ui(mantissa->z, 0);
        mpz_set_ui(exponent->z, 1);
    }
    else {
        temp = mpfr_get_z_2exp(mantissa->z, MPFR(self));
        mpz_set_si(exponent->z, temp);
    }
    result = Py_BuildValue("(NN)", (PyObject*)mantissa, (PyObject*)exponent);
    if (!result) {
        /* LCOV_EXCL_START */
        Py_DECREF((PyObject*)mantissa);
        Py_DECREF((PyObject*)exponent);
        /* LCOV_EXCL_STOP */
    }
    return result;
}

PyDoc_STRVAR(GMPy_doc_method_simple_fraction,
"x.as_simple_fraction(precision=0) -> mpq\n\n"
"Return a simple rational approximation to x. The result will be\n"
"accurate to 'precision' bits. If 'precision' is 0, the precision\n"
"of 'x' will be used.");

static PyObject *
GMPy_MPFR_Simple_Fraction_Method(PyObject *self, PyObject *args, PyObject *keywds)
{
    mpfr_prec_t prec = 0;
    static char *kwlist[] = {"precision", NULL};
    CTXT_Object *context = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "|l", kwlist, &prec)) {
        return NULL;
    }

    return (PyObject*)stern_brocot((MPFR_Object*)self, 0, prec, 0, context);
}

/* Implement the .precision attribute of an mpfr. */

static PyObject *
GMPy_MPFR_GetPrec_Attrib(MPFR_Object *self, void *closure)
{
    return PyLong_FromSsize_t((Py_ssize_t)mpfr_get_prec(self->f));
}

/* Implement the .rc attribute of an mpfr. */

static PyObject *
GMPy_MPFR_GetRc_Attrib(MPFR_Object *self, void *closure)
{
    return PyLong_FromLong((long)self->rc);
}

/* Implement the .imag attribute of an mpfr. */

static PyObject *
GMPy_MPFR_GetImag_Attrib(MPFR_Object *self, void *closure)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if ((result = GMPy_MPFR_New(0, context)))
        mpfr_set_zero(result->f, 1);
    return (PyObject*)result;
}

/* Implement the .real attribute of an mpfr. */

static PyObject *
GMPy_MPFR_GetReal_Attrib(MPFR_Object *self, void *closure)
{
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
}

static PyObject *
GMPy_MPFR_Get_Mpmath_MPF_Tuple(MPFR_Object *self, void *closure)
{
    long sign;
    mp_bitcnt_t bc;
    MPZ_Object *mantissa;
    MPZ_Object *exponent;
    mpfr_exp_t temp;
    PyObject *tuple;

    if (!(tuple = PyTuple_New(4))) {
        return NULL;
    }

    mantissa = GMPy_MPZ_New(NULL);
    exponent = GMPy_MPZ_New(NULL);
    if (!mantissa || !exponent) {
        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)mantissa);
        Py_XDECREF((PyObject*)exponent);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (mpfr_zero_p(MPFR(self))) {
        mpz_set_ui(MPZ(mantissa), 0);
        mpz_set_ui(MPZ(exponent), 1);
    }
    else {
        temp = mpfr_get_z_2exp(MPZ(mantissa), MPFR(self));
        mpz_set_si(MPZ(exponent), temp);
    }

    sign = mpz_sgn(MPZ(mantissa)) < 0;
    mpz_abs(MPZ(mantissa), MPZ(mantissa));
    bc = mpz_sizeinbase(MPZ(mantissa), 2);

    PyTuple_SET_ITEM(tuple, 0, PyLong_FromLong(sign));
    PyTuple_SET_ITEM(tuple, 1, (PyObject*)mantissa);
    PyTuple_SET_ITEM(tuple, 2, GMPy_PyLong_From_MPZ(exponent, NULL));
    PyTuple_SET_ITEM(tuple, 3, GMPy_PyLong_FromMpBitCnt(bc));

    return tuple;
}

/* Implement the nb_bool slot. */

static int
GMPy_MPFR_NonZero_Slot(MPFR_Object *self)
{
    return !mpfr_zero_p(self->f);
}

PyDoc_STRVAR(GMPy_doc_function_check_range,
"check_range(x, /) -> mpfr\n\n"
"Return a new `mpfr` with exponent that lies within the current range\n"
"of emin and emax.");

PyDoc_STRVAR(GMPy_doc_context_check_range,
"context.check_range(x, /) -> mpfr\n\n"
"Return a new `mpfr` with exponent that lies within the range of emin\n"
"and emax specified by context.");

static PyObject *
GMPy_MPFR_CheckRange(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result;

    CHECK_CONTEXT(context);

    if ((result = GMPy_MPFR_New(mpfr_get_prec(MPFR(x)), context))) {
        mpfr_set(result->f, MPFR(x), GET_MPFR_ROUND(context));
        mpfr_clear_flags();
        _GMPy_MPFR_Cleanup(&result, context);
    }
    return (PyObject*)result;
}

static PyObject *
GMPy_Number_CheckRange(PyObject *x, CTXT_Object *context)
{
    if (MPFR_Check(x))
        return GMPy_MPFR_CheckRange(x, context);

    TYPE_ERROR("check_range() argument types not supported");
    return NULL;
}

static PyObject *
GMPy_Context_CheckRange(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_CheckRange(other, context);
}

PyDoc_STRVAR(GMPy_doc_mpfr_sizeof_method,
"x.__sizeof__()\n\n"
"Returns the amount of memory consumed by x.");

static PyObject *
GMPy_MPFR_SizeOf_Method(PyObject *self, PyObject *other)
{
    return PyLong_FromSize_t(sizeof(MPFR_Object) + \
        (((MPFR(self))->_mpfr_prec + mp_bits_per_limb - 1) / \
        mp_bits_per_limb) * sizeof(mp_limb_t));
}

PyDoc_STRVAR(GMPy_doc_method_round10,
"x.__round__(n = 0, /) -> mpfr\n\n"
"Return x rounded to n decimal digits before (n < 0) or after (n > 0)\n"
"the decimal point. Rounds to an integer if n is not specified.");

static PyObject *
GMPy_MPFR_Method_Round10(PyObject *self, PyObject *args)
{
    long digits = 0;
    mpz_t temp;
    MPFR_Object *resultf = 0;
    MPZ_Object *resultz;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    /* If the size of args is 0, we just return an mpz. */

    if (PyTuple_GET_SIZE(args) == 0) {
        if ((resultz = GMPy_MPZ_New(context))) {
            if (mpfr_nan_p(MPFR(self))) {
                Py_DECREF((PyObject*)resultz);
                VALUE_ERROR("'mpz' does not support NaN");
                return NULL;
            }
            if (mpfr_inf_p(MPFR(self))) {
                Py_DECREF((PyObject*)resultz);
                OVERFLOW_ERROR("'mpz' does not support Infinity");
                return NULL;
            }
            /* return code is ignored */
            mpfr_get_z(resultz->z, MPFR(self), MPFR_RNDN);
        }
        return (PyObject*)resultz;
    }

    /* Now we need to return an mpfr, so handle the simple cases. */

    if (!mpfr_regular_p(MPFR(self))) {
        Py_INCREF(self);
        return self;
    }

    if (PyTuple_GET_SIZE(args) >= 1) {
        digits = PyLong_AsLong(PyTuple_GET_ITEM(args, 0));
        if (digits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("__round__() requires 'int' argument");
            return NULL;
        }
    }

    if (!(resultf = GMPy_MPFR_New(mpfr_get_prec(MPFR(self))+100, context))) {
        return NULL;
    }

    mpz_init(temp);
    mpz_ui_pow_ui(temp, 10, digits > 0 ? digits : -digits);
    if (digits >= 0) {
        mpfr_mul_z(resultf->f, MPFR(self), temp, MPFR_RNDN);
    }
    else {
        mpfr_div_z(resultf->f, MPFR(self), temp, MPFR_RNDN);
    }

    mpfr_rint(resultf->f, resultf->f, MPFR_RNDN);

    if (digits >= 0) {
        mpfr_div_z(resultf->f, resultf->f, temp, MPFR_RNDN);
    }
    else {
        mpfr_mul_z(resultf->f, resultf->f, temp, MPFR_RNDN);
    }
    mpfr_prec_round(resultf->f, mpfr_get_prec(MPFR(self)), MPFR_RNDN);

    mpz_clear(temp);
    return((PyObject*)resultf);
}
