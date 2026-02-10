/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpq_misc.c                                                        *
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

static PyObject *
GMPy_MPQ_Attrib_GetNumer(MPQ_Object *self, void *closure)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    if ((result = GMPy_MPZ_New(context)))
        mpz_set(result->z, mpq_numref(self->q));
    return (PyObject*)result;
}

static PyObject *
GMPy_MPQ_Attrib_GetReal(MPQ_Object *self, void *closure)
{
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
}

static PyObject *
GMPy_MPQ_Attrib_GetDenom(MPQ_Object *self, void *closure)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    if ((result = GMPy_MPZ_New(context)))
        mpz_set(result->z, mpq_denref(self->q));
    return (PyObject*)result;
}

static PyObject *
GMPy_MPQ_Attrib_GetImag(MPQ_Object *self, void *closure)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    if ((result = GMPy_MPZ_New(context)))
        mpz_set_ui(result->z, 0);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpq_function_numer,
"numer($module, x, /)\n--\n\n"
"Return the numerator of x.");

static PyObject *
GMPy_MPQ_Function_Numer(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    MPQ_Object *tempq;
    CTXT_Object *context = NULL;

    if (!(result = (MPZ_Object*)GMPy_MPZ_New(context)))
        return NULL;

    if (!(tempq = GMPy_MPQ_From_Rational(other, context))) {
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    mpz_set(result->z, mpq_numref(tempq->q));
    Py_DECREF((PyObject*)tempq);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpq_function_denom,
"denom($module, x, /)\n--\n\n"
"Return the denominator of x.");

static PyObject *
GMPy_MPQ_Function_Denom(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    MPQ_Object *tempq;
    CTXT_Object *context = NULL;

    if (!(result = (MPZ_Object*)GMPy_MPZ_New(context)))
        return NULL;

    if (!(tempq = GMPy_MPQ_From_Rational(other, context))) {
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    mpz_set(result->z, mpq_denref(tempq->q));
    Py_DECREF((PyObject*)tempq);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpq_method_as_integer_ratio,
"as_integer_ratio($self)\n--\n\n\
Return a pair of integers, whose ratio is exactly equal to the\n\
original number.  The ratio is in lowest terms and has a\n\
positive denominator.");
static PyObject *
GMPy_MPQ_Method_As_Integer_Ratio(PyObject *self, PyObject *Py_UNUSED(ignored))
{
    return PyTuple_Pack(2, GMPy_MPQ_Attrib_GetNumer((MPQ_Object*)self, NULL),
                        GMPy_MPQ_Attrib_GetDenom((MPQ_Object*)self, NULL));
}

PyDoc_STRVAR(GMPy_doc_mpq_method_from_float,
"from_float($self, f, /)\n--\n\n\
Converts a finite float to a rational number, exactly.");

PyDoc_STRVAR(GMPy_doc_mpq_method_from_decimal,
"from_decimal($self, dec, /)\n--\n\n\
Converts a finite `decimal.Decimal` instance to a rational number, exactly.");

static PyObject *
GMPy_MPQ_Method_From_As_Integer_Ratio(PyTypeObject *type, PyObject *const *args, Py_ssize_t nargs)
{
    PyObject *pair, *result;

    if (nargs != 1) {
        TYPE_ERROR("missing 1 required positional argument");
        return NULL;
    }

    pair = PyObject_CallMethod(args[0], "as_integer_ratio", NULL);
    if (pair == NULL) {
        return NULL;
    }

    result = GMPy_MPQ_NewInit(type, pair, NULL);
    Py_DECREF(pair);

    return result;
}

PyDoc_STRVAR(GMPy_doc_function_qdiv,
"qdiv($module, x, y=1, /)\n--\n\n"
"Return x/y as `mpz` if possible, or as `mpq` if x is not exactly\n"
"divisible by y.");

static PyObject *
GMPy_MPQ_Function_Qdiv(PyObject *self, PyObject *args)
{
    Py_ssize_t argc;
    PyObject *result = NULL, *x, *y;
    MPQ_Object *tempx = NULL, *tempy = NULL, *tempr = NULL;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    /* Validate the argument(s). */

    /* If there is only one argument, it should be either an integer or
     * rational type. If it is an integer, then immediately return an mpz().
     * If it is a rational type, convert it to an mpq and check the denominator.
     */

    argc = PyTuple_GET_SIZE(args);
    if (argc == 1) {
        x = PyTuple_GET_ITEM(args, 0);

        if (!IS_RATIONAL(x)) {
            goto arg_error;
        }

        if (IS_INTEGER(x)) {
            return (PyObject*)GMPy_MPZ_From_Integer(x, context);
        }

        if (!(tempx = GMPy_MPQ_From_Rational(x, context))) {
            return NULL;
        }

        if (mpz_cmp_ui(mpq_denref(tempx->q), 1) == 0) {
            if ((result = (PyObject*)GMPy_MPZ_New(context))) {
                mpz_set(MPZ(result), mpq_numref(tempx->q));
            }
            Py_DECREF((PyObject*)tempx);
            return result;
        }
        else {
            return (PyObject*)tempx;
        }
    }

    /* If there are two rational arguments, just convert them both to mpq,
     * divide, and then check the denominator.
     */

    if (argc == 2) {
        x = PyTuple_GET_ITEM(args, 0);
        y = PyTuple_GET_ITEM(args, 1);

        if (!IS_RATIONAL(x) || !IS_RATIONAL(y)) {
            goto arg_error;
        }

        if (!(tempx = GMPy_MPQ_From_Rational(x, context)) ||
            !(tempy = GMPy_MPQ_From_Rational(y, context))) {
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
        }

        if (mpq_sgn(tempy->q) == 0) {
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            ZERO_ERROR("qdiv() division by zero");
            return NULL;
        }

        /* tempr contains the result of the division. */

        if (!(tempr = GMPy_MPQ_New(context))) {
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            return NULL;
        }

        mpq_div(tempr->q, tempx->q, tempy->q);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);

        if (mpz_cmp_ui(mpq_denref(tempr->q), 1) == 0) {
            if ((result = (PyObject*)GMPy_MPZ_New(context))) {
                mpz_set(MPZ(result), mpq_numref(tempr->q));
            }
            Py_DECREF((PyObject*)tempr);
            return result;
        }
        else {
            return (PyObject*)tempr;
        };
    }

  arg_error:
    TYPE_ERROR("qdiv() requires 1 or 2 integer or rational arguments");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_mpq_method_floor,
"Return greatest integer less than or equal to an mpq.");

static PyObject *
GMPy_MPQ_Method_Floor(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if ((result = GMPy_MPZ_New(context))) {
        mpz_fdiv_q(result->z, mpq_numref(MPQ(self)), mpq_denref(MPQ(self)));
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpq_method_ceil,
"Return least integer greater than or equal to an mpq.");

static PyObject *
GMPy_MPQ_Method_Ceil(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if ((result = GMPy_MPZ_New(context))) {
        mpz_cdiv_q(result->z, mpq_numref(MPQ(self)), mpq_denref(MPQ(self)));
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpq_method_trunc,
"Return integer portion of an mpq.");

static PyObject *
GMPy_MPQ_Method_Trunc(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if ((result = GMPy_MPZ_New(context))) {
        mpz_tdiv_q(result->z, mpq_numref(MPQ(self)), mpq_denref(MPQ(self)));
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpq_method_round,
"__round__($self, ndigits=0, /)\n--\n\nRound an mpq to power of 10.");

static PyObject *
GMPy_MPQ_Method_Round(PyObject *self, PyObject *args)
{
    Py_ssize_t round_digits = 0;
    MPQ_Object *resultq;
    MPZ_Object *resultz;
    mpz_t temp;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    /* If args is NULL or the size of args is 0, we just return an mpz. */

    if (!args || PyTuple_GET_SIZE(args) == 0) {
        if (!(resultz = GMPy_MPZ_New(context))) {
            return NULL;
        }

        CTXT_Object *context = NULL;
        mpfr_rnd_t ctx_round;

        CHECK_CONTEXT(context);
        ctx_round = GET_MPFR_ROUND(context);
        if (ctx_round == MPFR_RNDA) {
            if (mpq_sgn(MPQ(self)) < 0) {
                ctx_round = MPFR_RNDD;
            }
            else {
                ctx_round = MPFR_RNDU;
            }
        }
        switch (ctx_round) {
            case MPFR_RNDZ:
                mpz_tdiv_q(resultz->z, mpq_numref(MPQ(self)),
                           mpq_denref(MPQ(self)));
                break;
            case MPFR_RNDD:
                mpz_fdiv_q(resultz->z, mpq_numref(MPQ(self)),
                           mpq_denref(MPQ(self)));
                break;
            case MPFR_RNDU:
                mpz_cdiv_q(resultz->z, mpq_numref(MPQ(self)),
                           mpq_denref(MPQ(self)));
                break;
            case MPFR_RNDN:
            default:
                {
                    mpz_t rem;

                    mpz_init(rem);
                    mpz_fdiv_qr(resultz->z, rem, mpq_numref(MPQ(self)),
                                mpq_denref(MPQ(self)));
                    mpz_mul_2exp(rem, rem, 1);
                    if (mpz_cmp(rem, mpq_denref(MPQ(self))) > 0) {
                        mpz_add_ui(resultz->z, resultz->z, 1);
                    }
                    else if (mpz_cmp(rem, mpq_denref(MPQ(self))) == 0) {
                        if (mpz_odd_p(resultz->z)) {
                            mpz_add_ui(resultz->z, resultz->z, 1);
                        }
                    }
                    mpz_clear(rem);
                }
            }
        return (PyObject*)resultz;
    }

    if (PyTuple_GET_SIZE(args) > 1) {
        TYPE_ERROR("Too many arguments for __round__()");
        return NULL;
    }

    if (PyTuple_GET_SIZE(args) == 1) {
        round_digits = PyLong_AsSsize_t(PyTuple_GET_ITEM(args, 0));
        if (round_digits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("__round__() requires 'int' argument");
            return NULL;
        }
    }

    if (!(resultq = GMPy_MPQ_New(context))) {
        return NULL;
    }

    mpz_init(temp);
    mpz_ui_pow_ui(temp, 10, round_digits > 0 ? (unsigned long)round_digits : (unsigned long)-round_digits);

    mpq_set(resultq->q, MPQ(self));
    if (round_digits > 0) {
        mpz_mul(mpq_numref(resultq->q), mpq_numref(resultq->q), temp);
        mpq_canonicalize(resultq->q);
        if (!(resultz = (MPZ_Object*)GMPy_MPQ_Method_Round((PyObject*)resultq, NULL))) {
            mpz_clear(temp);
            return NULL;
        }
        mpz_set(mpq_numref(resultq->q), resultz->z);
        Py_DECREF((PyObject*)resultz);
        mpz_set(mpq_denref(resultq->q), temp);
        mpz_clear(temp);
        mpq_canonicalize(resultq->q);
    }
    else {
        mpz_mul(mpq_denref(resultq->q), mpq_denref(resultq->q), temp);
        mpq_canonicalize(resultq->q);
        if (!(resultz = (MPZ_Object*)GMPy_MPQ_Method_Round((PyObject*)resultq, NULL))) {
            mpz_clear(temp);
            return NULL;
        }
        mpq_set_ui(resultq->q, 0, 1);
        mpz_mul(mpq_numref(resultq->q), resultz->z, temp);
        Py_DECREF((PyObject*)resultz);
        mpz_clear(temp);
        mpq_canonicalize(resultq->q);
    }
    return (PyObject*)resultq;
}

static int
GMPy_MPQ_NonZero_Slot(MPQ_Object *self)
{
    return mpq_sgn(self->q) != 0;
}

PyDoc_STRVAR(GMPy_doc_mpq_method_limit_denominator,
"limit_denominator($self, max_denominator=1000000)\n--\n\n"
"Closest fraction to self with denominator at most max_denominator.");

static PyObject *
GMPy_MPQ_Method_Limit_Denominator(PyObject *self, PyObject *const *args,
                         Py_ssize_t nargs, PyObject *kwnames)
{
    Py_ssize_t i, nkws = 0;
    PyObject *arg;
    int argidx[2] = {-1, -1};
    const char *kwname;

    if (nargs > 1) {
        TYPE_ERROR("limit_denominator() takes at most 1 positional arguments");
        return NULL;
    }
    if (nargs >= 1) {
        argidx[0] = 0;
    }

    if (kwnames) {
        nkws = PyTuple_GET_SIZE(kwnames);
    }
    if (nkws > 1) {
        TYPE_ERROR("limit_denominator() takes at most 1 keyword arguments");
        return NULL;
    }
    for (i = 0; i < nkws; i++) {
        kwname = PyUnicode_AsUTF8(PyTuple_GET_ITEM(kwnames, i));
        if (strcmp(kwname, "max_denominator") == 0) {
            if (nargs == 0) {
                argidx[0] = (int)(nargs + i);
            }
            else {
                TYPE_ERROR("argument for limit_denominator() given by name ('max_denominator') and position (1)");
                return NULL;
            }
        }
        else {
            TYPE_ERROR("got an invalid keyword argument for limit_denominator()");
            return NULL;
        }
    }

    mpz_t max_denominator;

    if (argidx[0] >= 0) {
        arg = args[argidx[0]];
        if (PyLong_Check(arg)) {
            arg = (PyObject *)GMPy_MPZ_From_PyLong(arg, NULL);
            if (!arg) {
                return NULL;
            }
            mpz_init_set(max_denominator, MPZ(arg));
            Py_DECREF(arg);
        }
        else {
            TYPE_ERROR("limit_denominator() takes an integer argument 'max_denominator'");
            return NULL;
        }
    }
    else {
        mpz_init_set_ui(max_denominator, 1000000);
    }
    if (mpz_cmp_ui(max_denominator, 1) <= 0) {
        VALUE_ERROR("max_denominator should be at least 1");
        return NULL;
    }

    /* Literal translation of the CPython's
     * Fraction.limit_denominator() follows.  See
     * https://github.com/python/cpython/issues/95723 for
     * correctness proof.
     */

    mpz_t n, d;

    mpz_init_set(n, mpq_numref(MPQ(self)));
    mpz_init_set(d, mpq_denref(MPQ(self)));

    if (mpz_cmp(d, max_denominator) <= 0) {
        mpz_clears(max_denominator, n, d, NULL);
        return self;
    }

    mpz_t p0, q0, p1, q1, a, q2, t;

    mpz_init_set_ui(p0, 0);
    mpz_init_set_ui(q0, 1);
    mpz_init_set_ui(p1, 1);
    mpz_init_set_ui(q1, 0);
    mpz_inits(a, q2, t, NULL);
    while (1) {
        mpz_fdiv_q(a, n, d);
        mpz_mul(t, a, q1);
        mpz_add(q2, q0, t);
        if (mpz_cmp(q2, max_denominator) > 0) {
            break;
        }
        mpz_set(q0, p0);
        mpz_set(p0, p1);
        mpz_mul(t, a, p1);
        mpz_add(p1, q0, t);
        mpz_set(q0, q1);
        mpz_set(q1, q2);

        mpz_mul(a, a, d);
        mpz_set(t, n);
        mpz_set(n, d);
        mpz_sub(d, t, a);
    }
    mpz_sub(t, max_denominator, q0);
    mpz_fdiv_q(t, t, q1);
    mpz_mul(a, t, q1);
    mpz_add(q0, q0, a);
    mpz_mul(a, q0, d);
    mpz_mul_ui(a, a, 2);

    MPQ_Object *result = GMPy_MPQ_New(NULL);

    if (!result) {
        /* LCOV_EXCL_START */
        mpz_clears(n, d, p0, p1, q0, q1, a, q2, t, NULL);
        return NULL;
        /* LCOV_EXCL_STOP */
    }
    if (mpz_cmp(a, mpq_denref(MPQ(self))) <= 0) {
        mpq_set_num(MPQ(result), p1);
        mpq_set_den(MPQ(result), q1);
    }
    else {
        mpz_mul(t, t, p1);
        mpz_add(p0, p0, t);
        mpq_set_num(MPQ(result), p0);
        mpq_set_den(MPQ(result), q0);
    }
    mpz_clears(n, d, p0, p1, q0, q1, a, q2, t, NULL);
    return (PyObject *)result;
}

PyDoc_STRVAR(GMPy_doc_mpq_method_sizeof,
"__sizeof__($self)\n--\n\n"
"Returns the amount of memory consumed by self. Note: deleted mpq objects\n"
"are reused and may or may not be resized when a new value is assigned.");

static PyObject *
GMPy_MPQ_Method_Sizeof(PyObject *self, PyObject *other)
{
    return PyLong_FromSize_t(sizeof(MPQ_Object) + \
        (mpq_numref(MPQ(self))->_mp_alloc * sizeof(mp_limb_t)) + \
        (mpq_denref(MPQ(self))->_mp_alloc * sizeof(mp_limb_t)));
}

PyDoc_STRVAR(GMPy_doc_mpq_method_is_integer,
"is_integer($self)\n--\n\n"
"Return `True` if self is an integer.");

static PyObject *
GMPy_MPQ_Method_IsInteger(PyObject *self, PyObject *other)
{
    return PyBool_FromLong(!mpz_cmp_ui(mpq_denref(MPQ(self)), 1));
}
