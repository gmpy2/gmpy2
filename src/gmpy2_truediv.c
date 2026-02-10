/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_truediv.c                                                         *
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

#include "fr_div.c"

/* This file implements the / operator, gmpy2.div() and context.div().
 *
 * Public API
 * ==========
 * The following function is available as part of GMPY2's C API. A NULL value
 * for context implies the function should use the currently active context.
 */

/* Divide two Integer objects (see convert.c/isInteger) using true division.
 * If an error occurs, NULL is returned and an exception is set. If either x
 * or y can't be converted into an mpz, Py_NotImplemented is returned. */

static PyObject *
GMPy_Integer_TrueDivWithType(PyObject *x, int xtype, PyObject *y, int ytype, 
                             CTXT_Object *context)
{
    MPZ_Object *tempx = NULL, *tempy = NULL;
    mpq_t tempq;
    MPFR_Object *result = NULL;

    CHECK_CONTEXT(context);

    if (GET_DIV_MODE(context))
        return GMPy_Rational_TrueDivWithType(x, xtype, y, ytype, context);

    if (!(result = GMPy_MPFR_New(0, context)))
        return NULL;

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype)) {
        if (!(tempx = GMPy_MPZ_From_IntegerWithType(x, xtype, context)) ||
            !(tempy = GMPy_MPZ_From_IntegerWithType(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        if (mpz_sgn(tempy->z) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }

        mpq_init(tempq);
        mpq_set_num(tempq, tempx->z);
        mpq_set_den(tempq, tempy->z);
        mpq_canonicalize(tempq);

        mpfr_clear_flags();

        result->rc = mpfr_set_q(result->f, tempq, GET_MPFR_ROUND(context));

        mpq_clear(tempq);

        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        _GMPy_MPFR_Cleanup(&result, context);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    TYPE_ERROR("div() argument type not supported");
    return NULL;
    /* LCOV_EXCL_STOP */
}
static PyObject *
GMPy_Rational_TrueDivWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                              CTXT_Object *context)
{
    MPQ_Object *result = NULL, *tempx = NULL, *tempy = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPQ_New(context)))
        return NULL;

    if (IS_TYPE_MPQ(xtype) && IS_TYPE_MPQ(ytype)) {
        if (mpq_sgn(MPQ(y)) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpq_div(result->q, MPQ(x), MPQ(y));
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        return (PyObject*)result;
    }

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype)) {
        if (!(tempx = GMPy_MPQ_From_RationalWithType(x, xtype, context)) ||
            !(tempy = GMPy_MPQ_From_RationalWithType(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        if (mpq_sgn(tempy->q) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }

        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpq_div(result->q, tempx->q, tempy->q);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    TYPE_ERROR("div() argument type not supported");
    return NULL;
    /* LCOV_EXCL_STOP */
}

/* Attempt true division of two numbers and return an mpfr. The code path is
 * optimized by checking for mpfr objects first. Returns Py_NotImplemented if
 * both objects are not valid reals.  */

static PyObject *
GMPy_Real_TrueDivWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                          CTXT_Object *context)
{
    MPFR_Object *result = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (IS_TYPE_MPFR(xtype) && IS_TYPE_MPFR(ytype)) {
        mpfr_clear_flags();

        result->rc = mpfr_div(result->f, MPFR(x), MPFR(y), GET_MPFR_ROUND(context));
        _GMPy_MPFR_Cleanup(&result, context);
        return (PyObject*)result;
    }

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype)) {
        MPFR_Object *tempx = NULL, *tempy = NULL;

        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context)) ||
            !(tempy = GMPy_MPFR_From_RealWithType(y, ytype, 1, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpfr_clear_flags();
        result->rc = mpfr_div(result->f, tempx->f, tempy->f, GET_MPFR_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        _GMPy_MPFR_Cleanup(&result, context);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    TYPE_ERROR("div() argument type not supported");
    return NULL;
    /* LCOV_EXCL_STOP */
}

static PyObject *
GMPy_Complex_TrueDivWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                             CTXT_Object *context)
{
    MPC_Object *result = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (IS_TYPE_MPC(xtype) && IS_TYPE_MPC(ytype)) {

        if (MPC_IS_ZERO_P(y)) {
            context->ctx.divzero = 1;
            if (context->ctx.traps & TRAP_DIVZERO) {
                GMPY_DIVZERO("'mpc' division by zero");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
        }
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        result->rc = mpc_div(result->c, MPC(x), MPC(y), GET_MPC_ROUND(context));
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        _GMPy_MPC_Cleanup(&result, context);
        return (PyObject*)result;
    }

    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_REAL(ytype)) {
        MPC_Object *tempx = NULL;
        MPFR_Object *tempy = NULL;

        if (mpfr_zero_p(MPFR(y))) {
            context->ctx.divzero = 1;
            if (context->ctx.traps & TRAP_DIVZERO) {
                GMPY_DIVZERO("'mpc' division by zero");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
        }

        if (!(tempx = GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context)) ||
            !(tempy = GMPy_MPFR_From_RealWithType(y, ytype, 1, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        result->rc = mpc_div_fr(result->c, tempx->c, MPFR(tempy), GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        _GMPy_MPC_Cleanup(&result, context);
        return (PyObject*)result;
    }

    if (IS_TYPE_COMPLEX(ytype) && IS_TYPE_REAL(xtype)) {
        MPC_Object *tempy = NULL;
        MPFR_Object *tempx = NULL;

        if (MPC_IS_ZERO_P(y)) {
            context->ctx.divzero = 1;
            if (context->ctx.traps & TRAP_DIVZERO) {
                GMPY_DIVZERO("'mpc' division by zero");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
        }

        if (!(tempy = GMPy_MPC_From_ComplexWithType(y, ytype, 1, 1, context)) ||
            !(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        result->rc = mpc_fr_div(result->c, MPFR(tempx), tempy->c, GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        _GMPy_MPC_Cleanup(&result, context);
        return (PyObject*)result;
    }

    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype)) {
        MPC_Object *tempx = NULL, *tempy = NULL;

        if (MPC_IS_ZERO_P(y)) {
            context->ctx.divzero = 1;
            if (context->ctx.traps & TRAP_DIVZERO) {
                GMPY_DIVZERO("'mpc' division by zero");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
        }

        if (!(tempx = GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context)) ||
            !(tempy = GMPy_MPC_From_ComplexWithType(y, ytype, 1, 1, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        result->rc = mpc_div(result->c, tempx->c, tempy->c, GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        _GMPy_MPC_Cleanup(&result, context);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    TYPE_ERROR("div() argument type not supported");
    return NULL;
    /* LCOV_EXCL_STOP */
}

/* Implement all the slot methods here. */

static PyObject *
GMPy_Number_TrueDiv_Slot(PyObject *x, PyObject *y)
{
    int xtype = GMPy_ObjectType(x);
    int ytype = GMPy_ObjectType(y);

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype))
        return GMPy_Integer_TrueDivWithType(x, xtype, y, ytype, NULL);

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype))
        return GMPy_Rational_TrueDivWithType(x, xtype, y, ytype, NULL);

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype))
        return GMPy_Real_TrueDivWithType(x, xtype, y, ytype, NULL);

    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype))
        return GMPy_Complex_TrueDivWithType(x, xtype, y, ytype, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

PyDoc_STRVAR(GMPy_doc_truediv,
"div($module, x, y, /)\n--\n\n"
"Return x / y; uses true division.");

static PyObject *
GMPy_Number_TrueDiv(PyObject *x, PyObject *y, CTXT_Object *context)
{
    int xtype = GMPy_ObjectType(x);
    int ytype = GMPy_ObjectType(y);

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype))
        return GMPy_Integer_TrueDivWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype))
        return GMPy_Rational_TrueDivWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype))
        return GMPy_Real_TrueDivWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype))
        return GMPy_Complex_TrueDivWithType(x, xtype, y, ytype, context);

    TYPE_ERROR("div() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_context_truediv,
"div($self, x, y, /)\n--\n\n"
"Return x / y; uses true division.");

static PyObject *
GMPy_Context_TrueDiv(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("div() requires 2 arguments.");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_TrueDiv(PyTuple_GET_ITEM(args, 0),
                               PyTuple_GET_ITEM(args, 1),
                               context);
}
