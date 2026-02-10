/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mod.c                                                             *
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

/* This file implements __mod__, gmpy2.mod(), and context.mod().
 */

static PyObject *
GMPy_Integer_ModWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                         CTXT_Object *context)
{
    MPZ_Object *result;

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
            mpz_fdiv_r(result->z, MPZ(x), MPZ(y));
            GMPY_MAYBE_END_ALLOW_THREADS(context);
            return (PyObject*)result;
        }

        if (IS_TYPE_PyInteger(ytype)) {
            int error;
            long temp = PyLong_AsLongAndOverflow(y, &error);

            if (!error) {
                if (temp > 0) {
                    mpz_fdiv_r_ui(result->z, MPZ(x), temp);
                }
                else if (temp == 0) {
                    ZERO_ERROR("division or modulo by zero");
                    Py_DECREF((PyObject*)result);
                    return NULL;
                }
                else {
                    mpz_cdiv_r_ui(result->z, MPZ(x), -temp);
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
                mpz_fdiv_r(result->z, MPZ(x), result->z);
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

        if (PyLong_Check(x)) {
            if (mpz_set_PyLong(result->z, x)) {
                /* LCOV_EXCL_START */
                Py_DECREF((PyObject*)result);
                return NULL;
                /* LCOV_EXCL_STOP */
            }
            mpz_fdiv_r(result->z, result->z, MPZ(y));
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
        mpz_fdiv_r(result->z, tempx->z, tempy->z);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    TYPE_ERROR("mod() argument type not supported");
    return NULL;
    /* LCOV_EXCL_STOP */
}

static PyObject *
GMPy_Rational_ModWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                          CTXT_Object *context)
{
    MPQ_Object *tempx = NULL, *tempy = NULL, *result = NULL;
    MPZ_Object *tempz = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPQ_New(context)) ||
        !(tempz = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)tempz);
        Py_XDECREF((PyObject*)result);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype)) {
        if (!(tempx = GMPy_MPQ_From_RationalWithType(x, xtype, context)) ||
            !(tempy = GMPy_MPQ_From_RationalWithType(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)tempz);
            Py_DECREF((PyObject*)result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        if (mpq_sgn(tempy->q) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)tempz);
            Py_DECREF((PyObject*)result);
            return NULL;
        }

        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpq_div(result->q, tempx->q, tempy->q);
        mpz_fdiv_q(tempz->z, mpq_numref(result->q), mpq_denref(result->q));
        /* Need to calculate x - tempz * y. */
        mpq_set_z(result->q, tempz->z);
        mpq_mul(result->q, result->q, tempy->q);
        mpq_sub(result->q, tempx->q, result->q);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        Py_DECREF((PyObject*)tempz);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    TYPE_ERROR("mod() argument type not supported");
    return NULL;
    /* LCOV_EXCL_STOP */
}

static PyObject *
GMPy_Real_ModWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                      CTXT_Object *context)
{
    MPFR_Object *tempx = NULL, *tempy = NULL, *result = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
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
                GMPY_DIVZERO("mod() modulo by zero");
                goto error;
            }
        }

        mpfr_clear_flags();

        if (mpfr_nan_p(tempx->f) || mpfr_nan_p(tempy->f) || mpfr_inf_p(tempx->f)) {
            context->ctx.invalid = 1;
            if (context->ctx.traps & TRAP_INVALID) {
                GMPY_INVALID("mod() invalid operation");
                goto error;
            }
            mpfr_set_nan(result->f);
        }
        else if (mpfr_inf_p(tempy->f)) {
            context->ctx.invalid = 1;
            if (context->ctx.traps & TRAP_INVALID) {
                GMPY_INVALID("mod() invalid operation");
                goto error;
            }
            if (mpfr_signbit(tempy->f)) {
                mpfr_set_inf(result->f, -1);
            }
            else {
                result->rc = mpfr_set(result->f, tempx->f,
                                      GET_MPFR_ROUND(context));
            }
        }
        else {
            mpfr_fmod(result->f, tempx->f, tempy->f, GET_MPFR_ROUND(context));

            if (!mpfr_zero_p(result->f)) {
                if ((mpfr_sgn(tempy->f) < 0) != (mpfr_sgn(result->f) < 0)) {
                    mpfr_add(result->f, result->f, tempy->f, GET_MPFR_ROUND(context));
                }
            }
            else {
                mpfr_copysign(result->f, result->f, tempy->f, GET_MPFR_ROUND(context));
            }
        }
        _GMPy_MPFR_Cleanup(&result, context);

        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    TYPE_ERROR("mod() argument type not supported");
    return NULL;
    /* LCOV_EXCL_STOP */

  error:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)result);
    return NULL;
}

static PyObject *
GMPy_Complex_ModWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                         CTXT_Object *context)
{
    TYPE_ERROR("can't take mod of complex number");
    return NULL;
}

static PyObject *
GMPy_Number_Mod_Slot(PyObject *x, PyObject *y)
{
    int xtype = GMPy_ObjectType(x);
    int ytype = GMPy_ObjectType(y);

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype))
        return GMPy_Integer_ModWithType(x, xtype, y, ytype, NULL);

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype))
        return GMPy_Rational_ModWithType(x, xtype, y, ytype, NULL);

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype))
        return GMPy_Real_ModWithType(x, xtype, y, ytype, NULL);

    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype))
        return GMPy_Complex_ModWithType(x, xtype, y, ytype, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

PyDoc_STRVAR(GMPy_doc_mod,
"mod($module, x, y, /)\n--\n\n"
"Return mod(x, y).\n"
"Note: overflow, underflow, and inexact exceptions are not supported for\n"
"mpfr arguments to mod().");

static PyObject *
GMPy_Number_Mod(PyObject *x, PyObject *y, CTXT_Object *context)
{
    int xtype = GMPy_ObjectType(x);
    int ytype = GMPy_ObjectType(y);

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype))
        return GMPy_Integer_ModWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype))
        return GMPy_Rational_ModWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype))
        return GMPy_Real_ModWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype))
        return GMPy_Complex_ModWithType(x, xtype, y, ytype, context);

    TYPE_ERROR("mod() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_context_mod,
"mod($self, x, y, /)\n--\n\n"
"Return mod(x, y).\n"
"Note: overflow, underflow, and inexact exceptions are not supported for\n"
"`mpfr` arguments.");

static PyObject *
GMPy_Context_Mod(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("mod() requires 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Mod(PyTuple_GET_ITEM(args, 0),
                           PyTuple_GET_ITEM(args, 1),
                           context);
}
