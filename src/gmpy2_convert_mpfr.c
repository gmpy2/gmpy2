/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_convert_mpfr.c                                                     *
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

/* This file contains all the conversion functions for MPFR data types.
 *
 * Overview
 * --------
 * gmpy2 tries to optimize the performance and accuracy of conversions from
 * other numeric types. The basic operations (+, -, *, /) are optimized to
 * directly work with the basic types such as C longs or doubles.
 *
 * gmpy2 supports two different strategies for creating new references to
 * an mpfr instance. If bits (or prec) is set to 0, the precison of the
 * result exactly matches the precision of the context. The exponent range
 * is also limited to the exponents defined in the contest. This is the
 * default behavior of the mpfr() function.
 *
 * If bits (or prec) is set to 1, the precision of the result depends on the
 * type of the source.
 *
 *   If the source number is already a radix-2 floating-point number,
 *   the precision is not changed. In practical terms, this only applies
 *   to sources operands that are either an mpfr or Python double.
 *
 * In addition, the exponent range is taken from the global emin/emax values
 * set in the MPFR library.
 *
 */

/* If prec == 0:
 *   Return an reference to a new mpfr instance. The context argument
 *   specifies the precision, exponent range, and whether or not subnormals
 *   are allowed.
 *
 *   This is the default behavior of the mpfr() constructor.
 *
 * If prec == 1:
 *   Return an additional reference to an existing mpfr. The precision is not
 *   changed. The exponent is not checked and may be outside the bounds
 *   specified in the context argument.
 *
 *   This is used internally for parsing arguments to functions.
 *
 * If prec >= 2:
 *   Return a reference to a new mpfr instance. The precision is taken from
 *   the argument. n instance with the specified precision and the
 *   exponent is valid in the current context.
 *
 *   This is used by mpfr() with an optional precision keyword.
 *
 * Since references to existing objects may be returned, the result should not
 * modified in-place.
 */

static MPFR_Object *
GMPy_MPFR_From_MPFR(MPFR_Object *obj, mpfr_prec_t prec, CTXT_Object *context)
{
    MPFR_Object *result = NULL;

    /* Optimize the critical case when prec==1 or obj is NaN or Inf. */

    if (prec == 1 || !mpfr_number_p(obj->f)) {
        Py_INCREF((PyObject*)obj);
        return obj;
    }

    CHECK_CONTEXT(context);

    if (prec == 0)
        prec = GET_MPFR_PREC(context);

    /* Try to identify when an additional reference to existing instance can
     * be returned. It is possible when (1) the precision matches, (2) the
     * exponent is valid and not in the range that might require subnormal-
     * ization, and (3) subnormalize is not enabled.
     */

    if (prec == mpfr_get_prec(obj->f) &&
        !context->ctx.subnormalize &&
        obj->f->_mpfr_exp >= (context->ctx.emin + mpfr_get_prec(obj->f) - 1) &&
        obj->f->_mpfr_exp <= context->ctx.emax
        ) {

        Py_INCREF((PyObject*)obj);
        return obj;
    }

    if ((result = GMPy_MPFR_New(prec, context))) {
        mpfr_clear_flags();
        result->rc = mpfr_set(result->f, obj->f, GET_MPFR_ROUND(context));
        _GMPy_MPFR_Cleanup(&result, context);
    }
    return result;
}

static MPFR_Object *
GMPy_MPFR_From_PyLong(PyObject *obj, mpfr_prec_t prec, CTXT_Object *context)
{
    MPFR_Object *result = NULL;
    MPZ_Object *tempx = NULL;
    int was_one = 0;
    long temp;

    CHECK_CONTEXT(context);

    if (prec == 0)
        prec = GET_MPFR_PREC(context);

    if (prec == 1) {
        /* This should be made accurate, but is it worth the overhead? */
        /* It is only used if the value fits in a C long. */
        prec = 64;
        was_one = 1;
    }

    temp = GMPy_Integer_AsLongWithType(obj, GMPy_ObjectType(obj));

    if ((temp == -1) && PyErr_Occurred()) {
        PyErr_Clear();
        if (!(tempx = GMPy_MPZ_From_PyLong(obj, context))) {
            return NULL;
        }
        if (was_one)
            prec = 1;
        result = GMPy_MPFR_From_MPZ(tempx, prec, context);
        Py_DECREF((PyObject*)tempx);
        return result;
    }
    else {
        if (!(result = GMPy_MPFR_New(prec, context)))
            return NULL;
        mpfr_clear_flags();
        result->rc = mpfr_set_si(result->f, temp, GET_MPFR_ROUND(context));
        if (!was_one) {
            GMPY_MPFR_CHECK_RANGE(result, context);
        }
        GMPY_MPFR_EXCEPTIONS(result, context);
    }

    return result;
}

/* If prec==0, then the precision of the current context is used.
 *
 * If prec==1, then the precision of Python float type is used (typically 53).
 *
 * If prec>=2, then the specified precision is used.
 */

static MPFR_Object *
GMPy_MPFR_From_PyFloat(PyObject *obj, mpfr_prec_t prec, CTXT_Object *context)
{
    MPFR_Object *result;

    CHECK_CONTEXT(context);

    if (prec == 0)
        prec = GET_MPFR_PREC(context);
    else if (prec == 1)
        prec = DBL_MANT_DIG;

    if ((result = GMPy_MPFR_New(prec, context))) {
        mpfr_clear_flags();
        result->rc = mpfr_set_d(result->f, PyFloat_AS_DOUBLE(obj), GET_MPFR_ROUND(context));
        if (prec != 1) {
            GMPY_MPFR_CHECK_RANGE(result, context);
        }
        GMPY_MPFR_SUBNORMALIZE(result, context);
        GMPY_MPFR_EXCEPTIONS(result, context);
    }
    return result;
}

/* If prec==0, then the precision of the current context is used.
 *
 * If prec==1, then an exact conversion is done by using the bit length of the
 * argument as the precision.
 *
 * If prec>=2, then the specified precision is used.
 */

static MPFR_Object *
GMPy_MPFR_From_MPZ(MPZ_Object *obj, mpfr_prec_t prec, CTXT_Object *context)
{
    MPFR_Object *result;
    int was_one = 0;
    size_t bitlen;

    CHECK_CONTEXT(context);

    if (prec == 0)
        prec = GET_MPFR_PREC(context);

    if (prec == 1) {
        bitlen = mpz_sizeinbase(obj->z, 2);
        if (bitlen < MPFR_PREC_MIN) {
            bitlen = MPFR_PREC_MIN;
	}
	if (bitlen > MPFR_PREC_MAX) {
	    OVERFLOW_ERROR("'mpz' to large to convert to 'mpfr'\n");
	    return NULL;
	}
	prec = (mpfr_prec_t)bitlen;
        was_one = 1;
    }

    if ((result = GMPy_MPFR_New(prec, context))) {
        mpfr_clear_flags();
        result->rc = mpfr_set_z(result->f, obj->z, GET_MPFR_ROUND(context));
        if (!was_one) {
            GMPY_MPFR_CHECK_RANGE(result, context);
        }
        GMPY_MPFR_EXCEPTIONS(result, context);
    }

    return result;
}

/* If prec<2, then the precision of the current context is used.
 *
 * If prec>=2, then the specified precision is used.
 */

static MPFR_Object *
GMPy_MPFR_From_MPQ(MPQ_Object *obj, mpfr_prec_t prec, CTXT_Object *context)
{
    MPFR_Object *result;

    CHECK_CONTEXT(context);

    if (prec < 2)
        prec = GET_MPFR_PREC(context);

    if ((result = GMPy_MPFR_New(prec, context))) {
        mpfr_clear_flags();
        result->rc = mpfr_set_q(result->f, obj->q, GET_MPFR_ROUND(context));
        GMPY_MPFR_CHECK_RANGE(result, context);
        GMPY_MPFR_SUBNORMALIZE(result, context);
        GMPY_MPFR_EXCEPTIONS(result, context);
    }

    return result;
}

/* If prec<2, then the precision of the current context is used.
 *
 * If prec>=2, then the specified precision is used.
 */

static MPFR_Object *
GMPy_MPFR_From_Fraction(PyObject *obj, mpfr_prec_t prec, CTXT_Object *context)
{
    MPFR_Object *result = NULL;
    MPQ_Object *tempq;

    CHECK_CONTEXT(context);

    if ((tempq = GMPy_MPQ_From_Fraction(obj, context))) {
        result = GMPy_MPFR_From_MPQ(tempq, prec, context);
        Py_DECREF((PyObject*)tempq);
    }

    return result;
}

static MPFR_Object *
GMPy_MPFR_From_PyStr(PyObject *s, int base, mpfr_prec_t prec, CTXT_Object *context)
{
    MPFR_Object *result;
    MPQ_Object *tempq;
    char *cp, *endptr;
    Py_ssize_t len;
    PyObject *ascii_str = ascii_str = GMPy_RemoveIgnoredASCII(s);

    if (!ascii_str) return NULL;

    CHECK_CONTEXT(context);

    if (prec < 2)
        prec = GET_MPFR_PREC(context);

    len = PyBytes_Size(ascii_str);
    cp = PyBytes_AsString(ascii_str);

    if (!(result = GMPy_MPFR_New(prec, context))) {
        Py_XDECREF(ascii_str);
        return NULL;
    }

    /* delegate the rest to MPFR */
    mpfr_clear_flags();
    result->rc = mpfr_strtofr(result->f, cp, &endptr, base, GET_MPFR_ROUND(context));
    Py_XDECREF(ascii_str);

    if (len != (Py_ssize_t)(endptr - cp)) {
        VALUE_ERROR("invalid digits");
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    /* If the context requests subnormals and the result is in the range for subnormals,
     * we use exact conversion via intermediate conversion to an mpq.
     *
     * The sticky bit returned by MPFR's string conversion appears to only reflect the
     * portion of the string needed to compute the correctly rounded result. It does not
     * accurately reflect whether or not the result is larger or smaller than the entire
     * input string. A correct sticky bit is needed by mfpr_subnormalize. Converting the
     * string to an mpq and then converting the mpq to an mpfr does properly set the
     * sticky bit.
     */

    if (mpfr_regular_p(result->f) &&
        base == 10 &&
        context->ctx.subnormalize &&
        result->f->_mpfr_exp <= context->ctx.emin + mpfr_get_prec(result->f) - 1)
        {

        if (!(tempq = GMPy_MPQ_From_PyStr(s, base, context))) {
            Py_DECREF((PyObject*)result);
            return NULL;
        }

        mpfr_clear_flags();
        result->rc = mpfr_set_q(result->f, tempq->q, GET_MPFR_ROUND(context));
        Py_DECREF((PyObject*)tempq);
    }

    GMPY_MPFR_CHECK_RANGE(result, context);
    GMPY_MPFR_SUBNORMALIZE(result, context);
    GMPY_MPFR_EXCEPTIONS(result, context);

    return result;
}

/* GMPy_MPFR_From_Real() converts a real number (see IS_REAL()) to an mpfr.
 *
 * If prec==0, then the result has the precision of the current context.
 *
 * If prec==1 and the value can be converted exactly (i.e. the input value is
 * a floating-point number using radix-2 representation or an integer), then
 * the conversion is done with the maximum possible precision. If the input
 * value can't be converted exactly, then the context precision is used.
 *
 * If prec >=2, then the specified precision is used.
 *
 * The return value is guaranteed to have a valid exponent.
 */

static MPFR_Object *
GMPy_MPFR_From_RealWithType(PyObject *obj, int xtype, mp_prec_t prec, CTXT_Object *context)
{
    CHECK_CONTEXT(context);

    if (IS_TYPE_MPFR(xtype))
        return GMPy_MPFR_From_MPFR((MPFR_Object*)obj, prec, context);

    if (IS_TYPE_PyFloat(xtype))
        return GMPy_MPFR_From_PyFloat(obj, prec, context);

    if (IS_TYPE_MPQ(xtype))
        return GMPy_MPFR_From_MPQ((MPQ_Object*)obj, prec, context);

    if (IS_TYPE_MPZANY(xtype))
        return GMPy_MPFR_From_MPZ((MPZ_Object*)obj, prec, context);

    if (IS_TYPE_PyInteger(xtype))
        return GMPy_MPFR_From_PyLong(obj, prec, context);

    if (IS_TYPE_PyFraction(xtype))
        return GMPy_MPFR_From_Fraction(obj, prec, context);

    if (IS_TYPE_HAS_MPFR(xtype)) {
        MPFR_Object *res = (MPFR_Object *) PyObject_CallMethod(obj, "__mpfr__", NULL);

        if (res != NULL && MPFR_Check(res)) {
            return res;
        }
        else {
            Py_XDECREF((PyObject*)res);
            goto error;
        }
    }

    if (IS_TYPE_HAS_MPQ(xtype)) {
        MPQ_Object *res = (MPQ_Object *) PyObject_CallMethod(obj, "__mpq__", NULL);

        if (res != NULL && MPQ_Check(res)) {
            MPFR_Object * temp =  GMPy_MPFR_From_MPQ(res, prec, context);
            Py_DECREF(res);
            return temp;
        }
        else {
            Py_XDECREF((PyObject*)res);
            goto error;
        }
    }

    if (IS_TYPE_HAS_MPZ(xtype)) {
        MPZ_Object *res = (MPZ_Object *) PyObject_CallMethod(obj, "__mpz__", NULL);

        if (res != NULL && MPZ_Check(res)) {
            MPFR_Object * temp =  GMPy_MPFR_From_MPZ(res, prec, context);
            Py_DECREF(res);
            return temp;
        }
        else {
            Py_XDECREF((PyObject*)res);
            goto error;
        }
    }

  error:
    TYPE_ERROR("object could not be converted to 'mpfr'");
    return NULL;
}

static MPFR_Object *
GMPy_MPFR_From_Real(PyObject *obj, mp_prec_t prec, CTXT_Object *context)
{
    return GMPy_MPFR_From_RealWithType(obj, GMPy_ObjectType(obj),
                                      prec, context);
}

static MPFR_Object *
GMPy_MPFR_From_RealWithTypeAndCopy(PyObject *obj, int xtype, mp_prec_t prec, CTXT_Object *context)
{
    MPFR_Object *result = NULL, *temp = NULL;

    result = GMPy_MPFR_From_RealWithType(obj, xtype, prec, context);

    if (result == NULL)
        return result;

    if (Py_REFCNT(result) == 1)
        return result;

    if (!(temp = GMPy_MPFR_New(mpfr_get_prec(result->f), context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    /* Since the precision of temp is the same as the precision of result,
     * there shouldn't be any rounding.
     */

    mpfr_set(temp->f, result->f, MPFR_RNDN);
    Py_DECREF((PyObject*)result);
    return temp;
}

static MPZ_Object *
GMPy_MPZ_From_MPFR(MPFR_Object *obj, CTXT_Object *context)
{
    MPZ_Object *result;

    CHECK_CONTEXT(context);

    if ((result = GMPy_MPZ_New(context))) {
        if (mpfr_nan_p(MPFR(obj))) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'mpz' does not support NaN");
            return NULL;
        }
        if (mpfr_inf_p(MPFR(obj))) {
            Py_DECREF((PyObject*)result);
            OVERFLOW_ERROR("'mpz' does not support Infinity");
            return NULL;
        }
        /* return code is ignored */
        mpfr_get_z(result->z, MPFR(obj), GET_MPFR_ROUND(context));
    }

    return result;
}

static XMPZ_Object *
GMPy_XMPZ_From_MPFR(MPFR_Object *self, CTXT_Object *context)
{
    XMPZ_Object *result;

    CHECK_CONTEXT(context);

    if ((result = GMPy_XMPZ_New(context))) {
        if (mpfr_nan_p(MPFR(self))) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'xmpz' does not support NaN");
            return NULL;
        }
        if (mpfr_inf_p(MPFR(self))) {
            Py_DECREF((PyObject*)result);
            OVERFLOW_ERROR("'xmpz' does not support Infinity");
            return NULL;
        }
        /* return code is ignored */
        mpfr_get_z(result->z, MPFR(self), GET_MPFR_ROUND(context));
    }

    return result;
}

/* Return the simpliest rational number that approximates 'self' to the
 * requested error bound.
 *
 * 'bits' is used as the working precision for the calculations.
 *   - If (bits == 0) then the precision of 'self' is used.
 *   - If (bits > precision of 'self') then use precision of 'self'.
 *
 * 'err' is the maximum error in the rational approximation.
 *   - If (err > 0) then prec
 *   - If (err == NULL), then the requested precision is 1/(2**prec).
 *       This should return the smallest fraction that returns the
 *       same interval that includes 'self'.
 *  'err'. If 'err' is negative, then the requested
 * precision is -2**abs(int(err)). If 'err' is NULL, then the requested
 * precision is -2**prec. If 'prec' is 0, then the requested precision is
 * the precision of 'self'.
 */


static PyObject *
stern_brocot(MPFR_Object* self, MPFR_Object *err, mpfr_prec_t bits, int mayz, CTXT_Object *context)
{
    PyObject* result = NULL;
    MPQ_Object *mpqres = NULL;
    MPZ_Object *mpzres = NULL;

    int i, negative, errsign;
    mpfr_t f, al, a, r1[3], r2[3], temperr, minerr, curerr, newerr, temp;

    if (mpfr_nan_p(self->f)) {
        VALUE_ERROR("Cannot convert NaN to a number.");
        return NULL;
    }

    if (mpfr_inf_p(self->f)) {
        OVERFLOW_ERROR("Cannot convert Infinity to a number.");
        return NULL;
    }

    if (!bits) {
        bits = mpfr_get_prec(self->f);
    }

    if (err) {
        /* Make a copy. */
        mpfr_init2(temperr, mpfr_get_prec(err->f));
        mpfr_set(temperr, err->f, MPFR_RNDN);
    }
    else {
        mpfr_init2(temperr, mpfr_get_prec(self->f));
        mpfr_set_ui(temperr, 0, MPFR_RNDN);
    }

    errsign = mpfr_sgn(temperr);
    if (errsign <= 0 && (bits < 2 || bits > mpfr_get_prec(self->f))) {
        VALUE_ERROR("Requested precision out-of-bounds.");
        mpfr_clear(temperr);
        return NULL;
    }

    if (errsign == 0) {
        mpfr_set_si(temperr, 1, MPFR_RNDN);
        mpfr_div_2exp(temperr, temperr, bits, MPFR_RNDN);
    }
    else if (errsign < 0) {
        long ubits;
        mpfr_abs(temperr, temperr, MPFR_RNDN);
        mpfr_floor(temperr, temperr);
        ubits = mpfr_get_si(temperr, MPFR_RNDN);
        if (ubits < 2 || ubits > mpfr_get_prec(self->f)) {
            VALUE_ERROR("Requested precision out-of-bounds.");
            mpfr_clear(temperr);
            return NULL;
        }
        mpfr_set_si(temperr, 1, MPFR_RNDN);
        mpfr_div_2exp(temperr, temperr, ubits, MPFR_RNDN);
    }

    if (!(mpqres = GMPy_MPQ_New(context)) ||
        !(mpzres = GMPy_MPZ_New(context))) {
        Py_XDECREF((PyObject*)mpqres);
        Py_XDECREF((PyObject*)mpzres);
        mpfr_clear(temperr);
        return NULL;
    }

    mpfr_init2(minerr, bits);
    mpfr_set(minerr, temperr, MPFR_RNDN);
    mpfr_clear(temperr);

    mpfr_init2(f, bits);
    if (mpfr_sgn(self->f) < 0) {
        negative = 1;
        mpfr_abs(f, self->f, MPFR_RNDN);
    }
    else {
        negative = 0;
        mpfr_set(f, self->f, MPFR_RNDN);
    }

    mpfr_init2(al, bits);
    mpfr_set(al, f, MPFR_RNDN);
    mpfr_init2(a, bits);
    mpfr_floor(a, al);
    mpfr_init2(temp, bits);
    for (i=0; i<3; ++i) {
        mpfr_init2(r1[i], bits);
        mpfr_init2(r2[i], bits);
    }

    mpfr_set_si(r1[0], 0, MPFR_RNDN);
    mpfr_set_si(r1[1], 0, MPFR_RNDN);
    mpfr_set_si(r1[2], 1, MPFR_RNDN);
    mpfr_set_si(r2[0], 0, MPFR_RNDN);
    mpfr_set_si(r2[1], 1, MPFR_RNDN);
    mpfr_set(r2[2], a, MPFR_RNDN);
    mpfr_init2(curerr, bits);
    mpfr_init2(newerr, bits);
    mpfr_reldiff(curerr, f, a, MPFR_RNDN);

    while (mpfr_cmp(curerr, minerr) > 0) {
        mpfr_sub(temp, al, a, MPFR_RNDN);
        mpfr_ui_div(al, 1, temp, MPFR_RNDN);
        mpfr_floor(a, al);
        mpfr_swap(r1[0], r1[1]);
        mpfr_swap(r1[1], r1[2]);
        //mpfr_fma(r1[2], r1[1], a, r1[0], MPFR_RNDN);
        mpfr_mul(r1[2], r1[1], a, MPFR_RNDN);
        mpfr_add(r1[2], r1[2], r1[0], MPFR_RNDN);
        mpfr_swap(r2[0], r2[1]);
        mpfr_swap(r2[1], r2[2]);
        //mpfr_fma(r2[2], r2[1], a, r2[0], MPFR_RNDN);
        mpfr_mul(r2[2], r2[1], a, MPFR_RNDN);
        mpfr_add(r2[2], r2[2], r2[0], MPFR_RNDN);
        mpfr_div(temp, r2[2], r1[2], MPFR_RNDN);
        mpfr_reldiff(newerr, f, temp, MPFR_RNDN);
        if(mpfr_cmp(curerr, newerr) <= 0) {
            mpfr_swap(r1[1],r1[2]);
            mpfr_swap(r2[1],r2[2]);
            break;
        }
        mpfr_swap(curerr, newerr);
    }

    /* Note: both mpqres and mpzrec have been created. Remember to delete the
     * one you don't need.
     */

    if (mayz && (mpfr_cmp_ui(r1[2],1) == 0)) {
        mpfr_get_z(mpzres->z, r2[2], MPFR_RNDN);
        if (negative) {
            mpz_neg(mpzres->z, mpzres->z);
        }
        result = (PyObject*)mpzres;
        Py_DECREF(mpqres);
    }
    else {
        mpfr_get_z(mpq_numref(mpqres->q), r2[2], MPFR_RNDN);
        mpfr_get_z(mpq_denref(mpqres->q), r1[2], MPFR_RNDN);
        if (negative) {
            mpz_neg(mpq_numref(mpqres->q), mpq_numref(mpqres->q));
        }
        result = (PyObject*)mpqres;
        Py_DECREF(mpzres);
    }

    mpfr_clear(minerr);
    mpfr_clear(al);
    mpfr_clear(a);
    mpfr_clear(f);
    for (i=0; i<3; ++i) {
        mpfr_clear(r1[i]);
        mpfr_clear(r2[i]);
    }
    mpfr_clear(curerr);
    mpfr_clear(newerr);
    mpfr_clear(temp);
    return result;
}

static MPQ_Object *
GMPy_MPQ_From_MPFR(MPFR_Object *self, CTXT_Object *context)
{
    mpfr_exp_t temp, twocount;
    MPQ_Object *result;

    CHECK_CONTEXT(context);

    if (mpfr_nan_p(self->f)) {
        VALUE_ERROR("can not convert NaN to MPQ");
        return NULL;
    }

    if (mpfr_inf_p(self->f)) {
        OVERFLOW_ERROR("can not convert Infinity to MPQ");
        return NULL;
    }

    if (!(result = GMPy_MPQ_New(context))) {
        return NULL;
    }

    if (mpfr_zero_p(self->f)) {
        mpz_set_ui(mpq_numref(result->q), 0);
        mpz_set_ui(mpq_denref(result->q), 1);
    }
    else {
        temp = mpfr_get_z_2exp(mpq_numref(result->q), self->f);
        twocount = (mpfr_exp_t)mpz_scan1(mpq_numref(result->q), 0);
        if (twocount) {
            temp += twocount;
            mpz_div_2exp(mpq_numref(result->q), mpq_numref(result->q), twocount);
        }
        mpz_set_ui(mpq_denref(result->q), 1);
        if (temp > 0)
            mpz_mul_2exp(mpq_numref(result->q), mpq_numref(result->q), temp);
        else if (temp < 0)
            mpz_mul_2exp(mpq_denref(result->q), mpq_denref(result->q), -temp);
    }
    return result;
}

static PyObject *
GMPy_PyLong_From_MPFR(MPFR_Object *obj, CTXT_Object *context)
{
    PyObject *result;
    MPZ_Object *tempz;

    CHECK_CONTEXT(context);

    if (!(tempz = GMPy_MPZ_From_MPFR(obj, context)))
        return NULL;

    result = GMPy_PyLong_From_MPZ(tempz, context);
    Py_DECREF((PyObject*)tempz);

    return result;
}

static PyObject *
GMPy_MPFR_Int_Slot(MPFR_Object *self) {
    return GMPy_PyLong_From_MPFR(self, NULL);
}

static PyObject *
GMPy_PyFloat_From_MPFR(MPFR_Object *self, CTXT_Object *context)
{
    CHECK_CONTEXT(context);

    double res = mpfr_get_d(self->f, GET_MPFR_ROUND(context));

    return PyFloat_FromDouble(res);
}

static PyObject *
GMPy_MPFR_Float_Slot(MPFR_Object *self) {
    return GMPy_PyFloat_From_MPFR(self, NULL);
}

static PyObject*
GMPy_PyStr_From_MPFR(MPFR_Object *self, int base, int digits, CTXT_Object *context)
{
    PyObject *result;
    char *buffer;
    mpfr_exp_t the_exp;

    CHECK_CONTEXT(context);

    /* check arguments are valid */
    if (!((base >= 2) && (base <= 62))) {
        VALUE_ERROR("base must be in the interval [2,62]");
        return NULL;
    }
    if ((digits < 0) || (digits == 1)) {
        VALUE_ERROR("digits must be 0 or >= 2");
        return NULL;
    }

    /* Process special cases first */
    if (!(mpfr_regular_p(self->f))) {
        if (mpfr_nan_p(self->f)) {
            return Py_BuildValue("(sii)", "nan", 0, 0);
        }
        else if (mpfr_inf_p(self->f) && !mpfr_signbit(self->f)) {
            return Py_BuildValue("(sii)", "inf", 0, 0);
        }
        else if (mpfr_inf_p(self->f) && mpfr_signbit(self->f)) {
            return Py_BuildValue("(sii)", "-inf", 0, 0);
        }
        /* 0 is not considered a 'regular" number */
        else if (mpfr_signbit(self->f)) {
            return Py_BuildValue("(sii)", "-0", 0, mpfr_get_prec(self->f));
        }
        else {
            return Py_BuildValue("(sii)", "0", 0, mpfr_get_prec(self->f));
        }
    }

    /* obtain digits-string and exponent */
    buffer = mpfr_get_str(0, &the_exp, base, digits, self->f, GET_MPFR_ROUND(context));
    if (!*buffer) {
        SYSTEM_ERROR("Internal error in Pympfr_To_PyStr");
        return NULL;
    }

    result = Py_BuildValue("(sii)", buffer, the_exp, mpfr_get_prec(self->f));
    mpfr_free_str(buffer);
    return result;
}

#ifdef SHARED
/*
 * coerce any number to a mpf
 */

int
GMPy_MPFR_ConvertArg(PyObject *arg, PyObject **ptr)
{
    MPFR_Object* newob = GMPy_MPFR_From_RealWithType(arg, GMPy_ObjectType(arg),
                                                     1, NULL);

    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        TYPE_ERROR("argument can not be converted to 'mpfr'");
        return 0;
    }
}
#endif

/* str and repr implementations for mpfr */
static PyObject *
GMPy_MPFR_Str_Slot(MPFR_Object *self)
{
    PyObject *result, *temp;
    long precision;
    char fmtstr[60];

    precision = (long)(log10(2) * (double)mpfr_get_prec(MPFR(self))) + 2;

    sprintf(fmtstr, "{0:.%ldg}", precision);

    temp = PyUnicode_FromString(fmtstr);
    if (!temp)
        return NULL;
    result = PyObject_CallMethod(temp, "format", "O", self);
    Py_DECREF(temp);
    return result;
}

static PyObject *
GMPy_MPFR_Repr_Slot(MPFR_Object *self)
{
    PyObject *result, *temp;
    long precision, bits;
    char fmtstr[60];

    bits = mpfr_get_prec(MPFR(self));
    precision = (long)(log10(2) * (double)bits) + 2;

    if (mpfr_number_p(MPFR(self)) && bits != DBL_MANT_DIG)
        sprintf(fmtstr, "mpfr('{0:.%ldg}',%ld)", precision, bits);
    else
        sprintf(fmtstr, "mpfr('{0:.%ldg}')", precision);

    temp = PyUnicode_FromString(fmtstr);
    if (!temp)
        return NULL;
    result = PyObject_CallMethod(temp, "format", "O", self);
    Py_DECREF(temp);
    return result;
}

static PyObject *
mpfr_ascii(mpfr_t self, int base, int digits, int round)
{
    PyObject *result;
    char *buffer;
    mpfr_exp_t the_exp;

    /* Process special cases first */
    if (!(mpfr_regular_p(self))) {
        if (mpfr_nan_p(self)) {
            return Py_BuildValue("(sii)", "nan", 0, 0);
        }
        else if (mpfr_inf_p(self) && !mpfr_signbit(self)) {
            return Py_BuildValue("(sii)", "inf", 0, 0);
        }
        else if (mpfr_inf_p(self) && mpfr_signbit(self)) {
            return Py_BuildValue("(sii)", "-inf", 0, 0);
        }
        /* 0 is not considered a 'regular" number */
        else if (mpfr_signbit(self)) {
            return Py_BuildValue("(sii)", "-0", 0, mpfr_get_prec(self));
        }
        else {
            return Py_BuildValue("(sii)", "0", 0, mpfr_get_prec(self));
        }
    }

    /* obtain digits-string and exponent */
    buffer = mpfr_get_str(0, &the_exp, base, digits, self, round);
    if (!*buffer) {
        SYSTEM_ERROR("Internal error in mpfr_ascii");
        return NULL;
    }

    result = Py_BuildValue("(sii)", buffer, the_exp, mpfr_get_prec(self));
    mpfr_free_str(buffer);
    return result;
}
