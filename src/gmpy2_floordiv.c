/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_floordiv.c                                                        *
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

/* This file implements the // operator, gmpy2.floor_div() and
 * context.floor_div().
 */

static PyObject *
GMPy_Integer_FloorDivWithType(PyObject *x, int xtype, PyObject *y, int ytype, 
                              CTXT_Object *context)
{
    MPZ_Object *result = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (IS_TYPE_MPZANY(xtype)) {
        if (IS_TYPE_MPZANY(ytype)) {
            if (mpz_sgn(MPZ(y)) == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
            mpz_fdiv_q(result->z, MPZ(x), MPZ(y));
            GMPY_MAYBE_END_ALLOW_THREADS(context);
            return (PyObject*)result;
        }

        if (IS_TYPE_PyInteger(ytype)) {
            int error;
            long temp = PyLong_AsLongAndOverflow(y, &error);

            if (!error) {
                if (temp > 0) {
                    mpz_fdiv_q_ui(result->z, MPZ(x), temp);
                }
                else if (temp == 0) {
                    ZERO_ERROR("division or modulo by zero");
                    Py_DECREF((PyObject*)result);
                    return NULL;
                }
                else {
                    mpz_cdiv_q_ui(result->z, MPZ(x), -temp);
                    mpz_neg(result->z, result->z);
                }
            }
            else {
                if (mpz_set_PyLong(result->z, y)) {
                    /* LCOV_EXCL_START */
                    Py_DECREF((PyObject*)result);
                    return NULL;
                    /* LCOV_EXCL_STOP */
                }
                GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
                mpz_fdiv_q(result->z, MPZ(x), result->z);
                GMPY_MAYBE_END_ALLOW_THREADS(context);
            }
            return (PyObject*)result;
        }
    }

    if (IS_TYPE_MPZANY(ytype)) {
        if (mpz_sgn(MPZ(y)) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_DECREF((PyObject*)result);
            return NULL;
        }

        if (IS_TYPE_PyInteger(xtype)) {
            if (mpz_set_PyLong(result->z, x)) {
                /* LCOV_EXCL_START */
                Py_DECREF((PyObject*)result);
                return NULL;
                /* LCOV_EXCL_STOP */
            }
            GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
            mpz_fdiv_q(result->z, result->z, MPZ(y));
            GMPY_MAYBE_END_ALLOW_THREADS(context);
            return (PyObject*)result;
        }
    }

    if (IS_TYPE_INTEGER(xtype) && (IS_TYPE_INTEGER(ytype))) {
        MPZ_Object *tempx = NULL, *tempy = NULL;

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

        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_fdiv_q(result->z, tempx->z, tempy->z);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    TYPE_ERROR("floor_div() argument type not supported");
    return NULL;
    /* LCOV_EXCL_STOP */
}

static PyObject *
GMPy_Rational_FloorDivWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                                CTXT_Object *context)
{
    MPZ_Object *result = NULL;
    MPQ_Object *tempq = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPZ_New(context)) ||
        !(tempq = GMPy_MPQ_New(context))) {
         /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempq);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (IS_TYPE_MPQ(xtype) && IS_TYPE_MPQ(ytype)) {
        if (mpq_sgn(MPQ(y)) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_DECREF((PyObject*)result);
            Py_DECREF((PyObject*)tempq);
            return NULL;
        }
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpq_div(tempq->q, MPQ(x), MPQ(y));
        mpz_fdiv_q(result->z, mpq_numref(tempq->q), mpq_denref(tempq->q));
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_DECREF((PyObject*)tempq);
        return (PyObject*)result;
    }

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype)) {
        MPQ_Object *tempx = NULL, *tempy = NULL;

        if (!(tempx = GMPy_MPQ_From_RationalWithType(x, xtype, context)) ||
            !(tempy = GMPy_MPQ_From_RationalWithType(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            Py_DECREF((PyObject*)tempq);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        if (mpq_sgn(tempy->q) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            Py_DECREF((PyObject*)tempq);
            return NULL;
        }

        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpq_div(tempq->q, tempx->q, tempy->q);
        mpz_fdiv_q(result->z, mpq_numref(tempq->q), mpq_denref(tempq->q));
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        Py_DECREF((PyObject*)tempq);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)tempq);
    Py_DECREF((PyObject*)result);
    TYPE_ERROR("floor_div() argument type not supported");
    return NULL;
    /* LCOV_EXCL_STOP */
}

/* Attempt floor division of two numbers and return an mpfr. The code path is
 * optimized by checking for mpfr objects first. Returns Py_NotImplemented if
 * both objects are not valid reals.  */

static PyObject *
GMPy_Real_FloorDivWithType(PyObject *x, int xtype, PyObject *y, int ytype,
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
            result->rc = mpfr_floor(result->f, result->f);
            _GMPy_MPFR_Cleanup(&result, context);
            return (PyObject*)result;
        }

    /* Handle the remaining cases.
     * Note: verify that MPZ if converted at full precision! */

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
        result->rc = mpfr_div(result->f, MPFR(tempx), MPFR(tempy), GET_MPFR_ROUND(context));
        result->rc = mpfr_floor(result->f, result->f);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        _GMPy_MPFR_Cleanup(&result, context);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    TYPE_ERROR("floor_div() argument type not supported");
    return NULL;
    /* LCOV_EXCL_STOP */
}

static PyObject *
GMPy_Complex_FloorDivWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                         CTXT_Object *context)
{
    TYPE_ERROR("can't take floor of complex number");
    return NULL;
}

/* Implement all the slot methods here. */

static PyObject *
GMPy_Number_FloorDiv_Slot(PyObject *x, PyObject *y)
{
    int xtype = GMPy_ObjectType(x);
    int ytype = GMPy_ObjectType(y);

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype))
        return GMPy_Integer_FloorDivWithType(x, xtype, y, ytype, NULL);

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype))
        return GMPy_Rational_FloorDivWithType(x, xtype, y, ytype, NULL);

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype))
        return GMPy_Real_FloorDivWithType(x, xtype, y, ytype, NULL);

    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype))
        return GMPy_Complex_FloorDivWithType(x, xtype, y, ytype, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

PyDoc_STRVAR(GMPy_doc_floordiv,
"floor_div($module, x, y, /)\n--\n\n"
"Return x // y; uses floor division.");

static PyObject *
GMPy_Number_FloorDiv(PyObject *x, PyObject *y, CTXT_Object *context)
{
    int xtype = GMPy_ObjectType(x);
    int ytype = GMPy_ObjectType(y);

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype))
        return GMPy_Integer_FloorDivWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype))
        return GMPy_Rational_FloorDivWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype))
        return GMPy_Real_FloorDivWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype))
        return GMPy_Complex_FloorDivWithType(x, xtype, y, ytype, context);

    TYPE_ERROR("floor_div() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_context_floordiv,
"floor_div($self, x, y, /)\n--\n\n"
"Return x // y; uses floor division.");

static PyObject *
GMPy_Context_FloorDiv(PyObject *self, PyObject *args)
{
    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("floor_div() requires 2 arguments");
        return NULL;
    }

    return GMPy_Number_FloorDiv(PyTuple_GET_ITEM(args, 0),
                                PyTuple_GET_ITEM(args, 1),
                                (CTXT_Object*)self);
}
