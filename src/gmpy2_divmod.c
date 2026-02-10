/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_divmod.c                                                          *
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

/* This file implements __divmod__ and context.divmod().
 */

static PyObject *
GMPy_Integer_DivModWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                            CTXT_Object *context)
{
    PyObject *result = NULL;
    MPZ_Object *tempx = NULL, *tempy = NULL, *rem = NULL, *quo = NULL;

    CHECK_CONTEXT(context);

    if (!(result = PyTuple_New(2)) ||
        !(rem = GMPy_MPZ_New(context)) ||
        !(quo = GMPy_MPZ_New(context))) {

        /* LCOV_EXCL_START */
        goto error;
        /* LCOV_EXCL_STOP */
    }

    if (IS_TYPE_MPZANY(xtype)) {
        if (IS_TYPE_MPZANY(ytype)) {
            if (mpz_sgn(MPZ(y)) == 0) {
                ZERO_ERROR("division or modulo by zero");
                goto error;
            }
            GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
            mpz_fdiv_qr(quo->z, rem->z, MPZ(x), MPZ(y));
            GMPY_MAYBE_END_ALLOW_THREADS(context);
            PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
            PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
            return result;
        }

        if (IS_TYPE_PyInteger(ytype)) {
            int error;
            long temp = PyLong_AsLongAndOverflow(y, &error);

            if (error) {
                /* Use quo->z as a temporary variable. */
                if (mpz_set_PyLong(quo->z, y)) {
                    /* LCOV_EXCL_START */
                    Py_DECREF((PyObject*)rem);
                    Py_DECREF((PyObject*)quo);
                    Py_DECREF((PyObject*)result);
                    return NULL;
                    /* LCOV_EXCL_STOP */
                }
                GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
                mpz_fdiv_qr(quo->z, rem->z, MPZ(x), quo->z);
                GMPY_MAYBE_END_ALLOW_THREADS(context);
            }
            else if (temp > 0) {
                mpz_fdiv_qr_ui(quo->z, rem->z, MPZ(x), temp);
            }
            else if (temp == 0) {
                ZERO_ERROR("division or modulo by zero");
                goto error;
            }
            else {
                mpz_cdiv_qr_ui(quo->z, rem->z, MPZ(x), -temp);
                mpz_neg(quo->z, quo->z);
            }
            PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
            PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
            return result;
        }
    }

    if (IS_TYPE_MPZANY(ytype)) {
        if (IS_TYPE_PyInteger(xtype)) {
            if (mpz_sgn(MPZ(y)) == 0) {
                ZERO_ERROR("division or modulo by zero");
                goto error;
            }
            else {
                if (mpz_set_PyLong(quo->z, x)) {
                    /* LCOV_EXCL_START */
                    Py_DECREF((PyObject*)rem);
                    Py_DECREF((PyObject*)quo);
                    Py_DECREF((PyObject*)result);
                    return NULL;
                    /* LCOV_EXCL_STOP */
                }
                GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
                mpz_fdiv_qr(quo->z, rem->z, quo->z, MPZ(y));
                GMPY_MAYBE_END_ALLOW_THREADS(context);
                PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
                PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
                return (PyObject*)result;
            }
        }
    }

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype)) {

        if (!(tempx = GMPy_MPZ_From_IntegerWithType(x, xtype, context)) ||
            !(tempy = GMPy_MPZ_From_IntegerWithType(y, ytype, context))) {

            /* LCOV_EXCL_START */
            goto error;
            /* LCOV_EXCL_STOP */
        }
        if (mpz_sgn(tempy->z) == 0) {
            ZERO_ERROR("division or modulo by zero");
            goto error;
        }
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_fdiv_qr(quo->z, rem->z, tempx->z, tempy->z);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
        PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
        return result;
    }

    /* LCOV_EXCL_START */
    TYPE_ERROR("divmod() arguments not supported");
    /* LCOV_EXCL_STOP */
  error:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_XDECREF((PyObject*)rem);
    Py_XDECREF((PyObject*)quo);
    Py_XDECREF(result);
    return NULL;
}

static PyObject *
GMPy_Rational_DivModWithType(PyObject *x, int xtype, PyObject *y, int ytype, 
                             CTXT_Object *context)
{
    MPQ_Object *tempx = NULL, *tempy = NULL, *rem = NULL;
    MPZ_Object *quo = NULL;
    PyObject *result = NULL;

    CHECK_CONTEXT(context);

    if (!(result = PyTuple_New(2)) ||
        !(rem = GMPy_MPQ_New(context)) ||
        !(quo = GMPy_MPZ_New(context))) {

        /* LCOV_EXCL_START */
        goto error;
        /* LCOV_EXCL_STOP */
    }

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype)) {

        if (!(tempx = GMPy_MPQ_From_RationalWithType(x, xtype, context)) ||
            !(tempy = GMPy_MPQ_From_RationalWithType(y, ytype, context))) {

            /* LCOV_EXCL_START */
            goto error;
            /* LCOV_EXCL_STOP */
        }
        if (mpq_sgn(tempy->q) == 0) {
            ZERO_ERROR("division or modulo by zero");
            goto error;
        }

        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpq_div(rem->q, tempx->q, tempy->q);
        mpz_fdiv_q(quo->z, mpq_numref(rem->q), mpq_denref(rem->q));
        /* Need to calculate x - quo * y. */
        mpq_set_z(rem->q, quo->z);
        mpq_mul(rem->q, rem->q, tempy->q);
        mpq_sub(rem->q, tempx->q, rem->q);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
        PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
        return result;
    }

    /* LCOV_EXCL_START */
    TYPE_ERROR("divmod() arguments not supported");
  error:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_XDECREF((PyObject*)rem);
    Py_XDECREF((PyObject*)quo);
    Py_XDECREF(result);
    return NULL;
    /* LCOV_EXCL_STOP */
}

static PyObject *
GMPy_Real_DivModWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                         CTXT_Object *context)
{
    MPFR_Object *tempx = NULL, *tempy = NULL, *quo = NULL, *rem = NULL, *temp;
    PyObject *result = NULL;

    CHECK_CONTEXT(context);

    if (!(result = PyTuple_New(2)) ||
        !(rem = GMPy_MPFR_New(0, context)) ||
        !(quo = GMPy_MPFR_New(0, context))) {

        /* LCOV_EXCL_START */
        goto error;
        /* LCOV_EXCL_STOP */
    }

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype)) {

        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context)) ||
            !(tempy = GMPy_MPFR_From_RealWithType(y, ytype, 1, context))) {

            /* LCOV_EXCL_START */
            goto error;
            /* LCOV_EXCL_STOP */
        }
        if (mpfr_zero_p(tempy->f)) {
            context->ctx.divzero = 1;
            if (context->ctx.traps & TRAP_DIVZERO) {
                GMPY_DIVZERO("divmod() division by zero");
                goto error;
            }
            mpfr_set_nan(quo->f);
            mpfr_set_nan(rem->f);
            goto okay;
        }

        if (mpfr_nan_p(tempx->f) || mpfr_nan_p(tempy->f) || mpfr_inf_p(tempx->f)) {
            context->ctx.invalid = 1;
            if (context->ctx.traps & TRAP_INVALID) {
                GMPY_INVALID("divmod() invalid operation");
                goto error;
            }
            mpfr_set_nan(quo->f);
            mpfr_set_nan(rem->f);
            goto okay;
        }

        if (mpfr_inf_p(tempy->f)) {
            context->ctx.invalid = 1;
            if (context->ctx.traps & TRAP_INVALID) {
                GMPY_INVALID("divmod() invalid operation");
                goto error;
            }
            if (mpfr_zero_p(tempx->f)) {
                mpfr_set_zero(quo->f, mpfr_sgn(tempy->f));
                mpfr_set_zero(rem->f, mpfr_sgn(tempy->f));
            }
            else if ((mpfr_signbit(tempx->f)) != (mpfr_signbit(tempy->f))) {
                mpfr_set_si(quo->f, -1, MPFR_RNDN);
                mpfr_set_inf(rem->f, mpfr_sgn(tempy->f));
            }
            else {
                mpfr_set_si(quo->f, 0, MPFR_RNDN);
                rem->rc = mpfr_set(rem->f, tempx->f, MPFR_RNDN);
            }
            goto okay;
        }

        if (!(temp = GMPy_MPFR_New(0, context))) {
            /* LCOV_EXCL_START */
            goto error;
            /* LCOV_EXCL_STOP */
        }
        mpfr_fmod(rem->f, tempx->f, tempy->f, MPFR_RNDN);
        mpfr_sub(temp->f, tempx->f, rem->f, MPFR_RNDN);
        mpfr_div(quo->f, temp->f, tempy->f, MPFR_RNDN);
        Py_DECREF((PyObject*)temp);

        if (!mpfr_zero_p(rem->f)) {
            if ((mpfr_sgn(tempy->f) < 0) != (mpfr_sgn(rem->f) < 0)) {
                mpfr_add(rem->f, rem->f, tempy->f, MPFR_RNDN);
                mpfr_sub_ui(quo->f, quo->f, 1, MPFR_RNDN);
            }
        }
        else {
            mpfr_copysign(rem->f, rem->f, tempy->f, MPFR_RNDN);
        }

        if (!mpfr_zero_p(quo->f)) {
            mpfr_round(quo->f, quo->f);
        }
        else {
            mpfr_setsign(quo->f, quo->f, mpfr_sgn(tempx->f) * mpfr_sgn(tempy->f) - 1, MPFR_RNDN);
        }

        GMPY_MPFR_CHECK_RANGE(quo, context);
        GMPY_MPFR_CHECK_RANGE(rem, context);
        GMPY_MPFR_SUBNORMALIZE(quo, context);
        GMPY_MPFR_SUBNORMALIZE(rem, context);
      okay:
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
        PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    TYPE_ERROR("divmod() arguments not supported");
  error:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_XDECREF((PyObject*)rem);
    Py_XDECREF((PyObject*)quo);
    Py_XDECREF(result);
    return NULL;
    /* LCOV_EXCL_STOP */
}

static PyObject *
GMPy_Complex_DivModWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                            CTXT_Object *context)
{
    TYPE_ERROR("can't take floor or mod of complex number.");
    return NULL;
}

static PyObject *
GMPy_Number_DivMod_Slot(PyObject *x, PyObject *y)
{
    int xtype = GMPy_ObjectType(x);
    int ytype = GMPy_ObjectType(y);

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype))
        return GMPy_Integer_DivModWithType(x, xtype, y, ytype, NULL);

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype))
        return GMPy_Rational_DivModWithType(x, xtype, y, ytype, NULL);

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype))
        return GMPy_Real_DivModWithType(x, xtype, y, ytype, NULL);

    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype))
        return GMPy_Complex_DivModWithType(x, xtype, y, ytype, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_Number_DivMod(PyObject *x, PyObject *y, CTXT_Object *context)
{
    int xtype = GMPy_ObjectType(x);
    int ytype = GMPy_ObjectType(y);

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype))
        return GMPy_Integer_DivModWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype))
        return GMPy_Rational_DivModWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype))
        return GMPy_Real_DivModWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype))
        return GMPy_Complex_DivModWithType(x, xtype, y, ytype, context);

    TYPE_ERROR("divmod() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_context_divmod,
"divmod($self, x, y, /)\n--\n\n"
"Return divmod(x, y).\n\n"
"Note: overflow, underflow, and inexact exceptions are not supported for\n"
"`mpfr` arguments.");

static PyObject *
GMPy_Context_DivMod(PyObject *self, PyObject *args)
{
    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("div_mod() requires 2 arguments");
        return NULL;
    }

    return GMPy_Number_DivMod(PyTuple_GET_ITEM(args, 0),
                              PyTuple_GET_ITEM(args, 1),
                              (CTXT_Object*)self);
}
