/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_sub.c                                                             *
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

/* This file implements the - operator, gmpy2.sub(), and context.sub().
 */

/* Subtract two Integer objects (see gmpy2_convert.c). If an error occurs,
 * NULL is returned and an exception is set. If either x or y can't be
 * converted into an mpz, Py_NotImplemented is returned. */

static PyObject *
GMPy_Integer_SubWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                         CTXT_Object *context)
{
    MPZ_Object *result = NULL;

    if (!(result = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (IS_TYPE_MPZANY(xtype)) {
        if (IS_TYPE_MPZANY(ytype)) {
            GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
            mpz_sub(result->z, MPZ(x), MPZ(y));
            GMPY_MAYBE_END_ALLOW_THREADS(context);
            return (PyObject*)result;
        }
        if (IS_TYPE_PyInteger(ytype)) {
            int error;
            long temp = PyLong_AsLongAndOverflow(y, &error);

            if (!error) {
                if (temp >= 0) {
                    mpz_sub_ui(result->z, MPZ(x), temp);
                }
                else {
                    mpz_add_ui(result->z, MPZ(x), -temp);
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
                mpz_sub(result->z, MPZ(x), result->z);
                GMPY_MAYBE_END_ALLOW_THREADS(context);
            }
            return (PyObject*)result;
        }
    }

    if (IS_TYPE_MPZANY(ytype)) {
        if (IS_TYPE_PyInteger(xtype)) {
            int error;
            long temp = PyLong_AsLongAndOverflow(x, &error);

            if (!error) {
                if (temp >= 0) {
                    mpz_ui_sub(result->z, temp, MPZ(y));
                }
                else {
                    mpz_add_ui(result->z, MPZ(y), -temp);
                    mpz_neg(result->z, result->z);
                }
            }
            else {
                if (mpz_set_PyLong(result->z, x)) {
                    /* LCOV_EXCL_START */
                    Py_DECREF((PyObject*)result);
                    return NULL;
                    /* LCOV_EXCL_STOP */
                }
                GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
                mpz_sub(result->z, result->z, MPZ(y));
                GMPY_MAYBE_END_ALLOW_THREADS(context);
            }
            return (PyObject*)result;
        }
    }

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype)) {
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

        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_sub(result->z, tempx->z, tempy->z);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    TYPE_ERROR("sub() argument type not supported");
    return NULL;
    /* LCOV_EXCL_STOP */
}

/* Subtract two Rational objects (see gmpy2_convert.h). Returns None and
 * raises TypeError if both objects are not valid rationals. Pympq_Sub_Rational
 * is intended to be called from GMPy_Number_Sub(). */

static PyObject *
GMPy_Rational_SubWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                          CTXT_Object *context)
{
    MPQ_Object *result = NULL;

    if (!(result = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (IS_TYPE_MPQ(xtype) && IS_TYPE_MPQ(ytype)) {
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpq_sub(result->q, MPQ(x), MPQ(y));
        GMPY_MAYBE_END_ALLOW_THREADS(context);
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
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpq_sub(result->q, tempx->q, tempy->q);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    TYPE_ERROR("sub() argument type not supported");
    return NULL;
    /* LCOV_EXCL_STOP */
}

/*   GMPy_Real_Sub(x, y, context) returns x-y using the provided context. If
 *   provided context is NULL, then the current context is used. If an error
 *   occurs, NULL is returned and an exception is set. If either x or y can't
 *   be converted to an mpfr, then Py_NotImplemented is returned.
 *   GMPy_Real_Sub() will not try to promote the result to a different type
 *   (i.e. mpc).
*/

static PyObject *
GMPy_Real_SubWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                      CTXT_Object *context)
{
    MPFR_Object *result = NULL;

    if (!(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (IS_TYPE_MPFR(xtype) && IS_TYPE_MPFR(ytype)) {
        mpfr_clear_flags();
        result->rc = mpfr_sub(result->f, MPFR(x), MPFR(y), GET_MPFR_ROUND(context));
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
        result->rc = mpfr_sub(result->f, MPFR(tempx), MPFR(tempy), GET_MPFR_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        goto done;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    TYPE_ERROR("sub() argument type not supported");
    return NULL;
    /* LCOV_EXCL_STOP */

  done:
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

/* GMPy_Complex_Sub(x, y, context) returns x-y using the provided context. If
 * context is NULL, then the current context is used. If an error occurs, NULL
 * is returned and an exception is set. If either x or y can't be converted to
 * an mpc, then Py_NotImplemented is returned. */

static PyObject *
GMPy_Complex_SubWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                         CTXT_Object *context)
{
    MPC_Object *result = NULL;

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (IS_TYPE_MPC(xtype) && IS_TYPE_MPC(ytype)) {
        result->rc = mpc_sub(result->c, MPC(x), MPC(y), GET_MPC_ROUND(context));
        _GMPy_MPC_Cleanup(&result, context);
        return (PyObject*)result;
    }

    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_REAL(ytype)) {
        MPC_Object *tempx = NULL;
        MPFR_Object *tempy = NULL;

        if (!(tempx = GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context)) ||
            !(tempy = GMPy_MPFR_From_RealWithType(y, ytype, 1, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        result->rc = mpc_sub_fr(result->c, tempx->c, MPFR(tempy), GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        _GMPy_MPC_Cleanup(&result, context);
        return (PyObject*)result;
    }

    if (IS_TYPE_COMPLEX(ytype) && IS_TYPE_REAL(xtype)) {
        MPC_Object *tempy = NULL;
        MPFR_Object *tempx = NULL;

        if (!(tempy = GMPy_MPC_From_ComplexWithType(y, ytype, 1, 1, context)) ||
            !(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        result->rc = mpc_fr_sub(result->c, MPFR(tempx), tempy->c, GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        _GMPy_MPC_Cleanup(&result, context);
        return (PyObject*)result;
    }

    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype)) {
        MPC_Object *tempx = NULL, *tempy = NULL;

        if (!(tempx = GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context)) ||
            !(tempy = GMPy_MPC_From_ComplexWithType(y, ytype, 1, 1, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        result->rc = mpc_sub(result->c, tempx->c, tempy->c, GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        _GMPy_MPC_Cleanup(&result, context);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    TYPE_ERROR("sub() argument type not supported");
    return NULL;
    /* LCOV_EXCL_STOP */
}

static PyObject *
GMPy_Number_Sub(PyObject *x, PyObject *y, CTXT_Object *context)
{
    CHECK_CONTEXT(context);

    int xtype = GMPy_ObjectType(x);
    int ytype = GMPy_ObjectType(y);

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype))
        return GMPy_Integer_SubWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype))
        return GMPy_Rational_SubWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype))
        return GMPy_Real_SubWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype))
        return GMPy_Complex_SubWithType(x, xtype, y, ytype, context);

    TYPE_ERROR("sub() argument type not supported");
    return NULL;
}

/* Implement all the slot methods here. */

static PyObject *
GMPy_Number_Sub_Slot(PyObject *x, PyObject *y)
{
    CTXT_Object *context = NULL;
    CHECK_CONTEXT(context);

    int xtype = GMPy_ObjectType(x);
    int ytype = GMPy_ObjectType(y);

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype))
        return GMPy_Integer_SubWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype))
        return GMPy_Rational_SubWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype))
        return GMPy_Real_SubWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype))
        return GMPy_Complex_SubWithType(x, xtype, y, ytype, context);

    Py_RETURN_NOTIMPLEMENTED;
}

/* Implement context.sub() and gmpy2.sub(). */

PyDoc_STRVAR(GMPy_doc_sub,
"sub($module, x, y, /)\n--\n\n"
"Return x - y.");

PyDoc_STRVAR(GMPy_doc_context_sub,
"sub($self, x, y, /)\n--\n\n"
"Return x - y.");

static PyObject *
GMPy_Context_Sub(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("sub() requires 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Sub(PyTuple_GET_ITEM(args, 0),
                           PyTuple_GET_ITEM(args, 1),
                           context);
}
