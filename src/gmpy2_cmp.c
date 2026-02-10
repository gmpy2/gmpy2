/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_cmp.c                                                             *
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

PyDoc_STRVAR(GMPy_doc_mpany_cmp,
"cmp($module, x, y, /)\n--\n\n"
"Return -1 if x < y; 0 if x = y; or 1 if x > y. Both x and y must be\n"
"integer, rational or real. Note: 0 is returned (and exception flag set)\n"
"if either argument is NaN.");

static PyObject * _return_cmp(int c)
{
    if (c < 0) return PyLong_FromLong(-1);
    if (c > 0) return PyLong_FromLong(1);
    return PyLong_FromLong(0);
}

static PyObject * _return_negated_cmp(int c)
{
    if (c < 0) return PyLong_FromLong(1);
    if (c > 0) return PyLong_FromLong(-1);
    return PyLong_FromLong(0);
}

static PyObject *
GMPy_MPANY_cmp(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *result = NULL;
    int xtype, ytype;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("cmp() requires 2 arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);

    xtype = GMPy_ObjectType(x);
    ytype = GMPy_ObjectType(y);

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype)) {
        MPZ_Object *tempx = NULL, *tempy = NULL;

        if (!(tempx = GMPy_MPZ_From_IntegerWithType(x, xtype, context)) ||
            !(tempy = GMPy_MPZ_From_IntegerWithType(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        result = _return_cmp(mpz_cmp(tempx->z, tempy->z));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return result;
    }

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_INTEGER(ytype)) {
        MPQ_Object *tempx = NULL;
        MPZ_Object *tempy = NULL;

        if (!(tempx = GMPy_MPQ_From_RationalWithType(x, xtype, context)) ||
            !(tempy = GMPy_MPZ_From_IntegerWithType(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        result = _return_cmp(mpq_cmp_z(tempx->q, tempy->z));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return result;
    }

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_RATIONAL(ytype)) {
        MPZ_Object *tempx = NULL;
        MPQ_Object *tempy = NULL;

        if (!(tempx = GMPy_MPZ_From_IntegerWithType(x, xtype, context)) ||
            !(tempy = GMPy_MPQ_From_RationalWithType(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        result = _return_negated_cmp(mpq_cmp_z(tempy->q, tempx->z));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return result;
    }

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype)) {
        MPQ_Object *tempx = NULL, *tempy = NULL;

        if (!(tempx = GMPy_MPQ_From_RationalWithType(x, xtype, context)) ||
            !(tempy = GMPy_MPQ_From_RationalWithType(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        result = _return_cmp(mpq_cmp(tempx->q, tempy->q));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return result;
    }

    /* We perform exact comparisons between the mpz, mpq, and mpfr types.
     */

    if (IS_TYPE_REAL(xtype) && IS_TYPE_INTEGER(ytype)) {
            MPFR_Object *tempx = NULL;
            MPZ_Object *tempy = NULL;

        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context)) ||
            !(tempy = GMPy_MPZ_From_IntegerWithType(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpfr_clear_flags();
        result =_return_cmp(mpfr_cmp_z(tempx->f, tempy->z));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        GMPY_CHECK_ERANGE(result, context, "invalid comparison with NaN");
        return result;
    }

    if (IS_TYPE_REAL(xtype) && IS_TYPE_RATIONAL(ytype)) {
            MPFR_Object *tempx = NULL;
            MPQ_Object *tempy = NULL;

        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context)) ||
            !(tempy = GMPy_MPQ_From_RationalWithType(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpfr_clear_flags();
        result =_return_cmp(mpfr_cmp_q(tempx->f, tempy->q));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        GMPY_CHECK_ERANGE(result, context, "invalid comparison with NaN");
        return result;
    }

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype)) {
            MPFR_Object *tempx = NULL;
            MPFR_Object *tempy = NULL;

        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context)) ||
            !(tempy = GMPy_MPFR_From_RealWithType(y, ytype, 1, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpfr_clear_flags();
        result =_return_cmp(mpfr_cmp(tempx->f, tempy->f));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        GMPY_CHECK_ERANGE(result, context, "invalid comparison with NaN");
        return result;
    }

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_REAL(ytype)) {
            MPZ_Object *tempx = NULL;
            MPFR_Object *tempy = NULL;

        if (!(tempx = GMPy_MPZ_From_IntegerWithType(x, xtype, context)) ||
            !(tempy = GMPy_MPFR_From_RealWithType(y, ytype, 1, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpfr_clear_flags();
        result =_return_negated_cmp(mpfr_cmp_z(tempy->f, tempx->z));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        GMPY_CHECK_ERANGE(result, context, "invalid comparison with NaN");
        return result;
    }

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_REAL(ytype)) {
            MPQ_Object *tempx = NULL;
            MPFR_Object *tempy = NULL;

        if (!(tempx = GMPy_MPQ_From_RationalWithType(x, xtype, context)) ||
            !(tempy = GMPy_MPFR_From_RealWithType(y, ytype, 1, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpfr_clear_flags();
        result =_return_negated_cmp(mpfr_cmp_q(tempy->f, tempx->q));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        GMPY_CHECK_ERANGE(result, context, "invalid comparison with NaN");
        return result;
    }

    TYPE_ERROR("cmp() requires integer, rational, or real arguments");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_mpany_cmp_abs,
"cmp_abs($module, x, y, /)\n--\n\n"
"Return -1 if abs(x) < abs(y); 0 if abs(x) = abs(y); or 1 else.");

static PyObject *
GMPy_MPANY_cmp_abs(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *result = NULL;
    int xtype, ytype;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("cmp() requires 2 arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);

    xtype = GMPy_ObjectType(x);
    ytype = GMPy_ObjectType(y);

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype)) {
        MPZ_Object *tempx = NULL, *tempy = NULL;

        if (!(tempx = GMPy_MPZ_From_IntegerWithType(x, xtype, context)) ||
            !(tempy = GMPy_MPZ_From_IntegerWithType(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        result = _return_cmp(mpz_cmpabs(tempx->z, tempy->z));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return result;
    }

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_INTEGER(ytype)) {
        MPQ_Object *tempx = NULL;
        MPZ_Object *tempy = NULL;

        if (!(tempx = GMPy_MPQ_From_RationalWithTypeAndCopy(x, xtype, context)) ||
            !(tempy = GMPy_MPZ_From_IntegerWithTypeAndCopy(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpq_abs(tempx->q, tempx->q);
        mpz_abs(tempy->z, tempy->z);

        result = _return_cmp(mpq_cmp_z(tempx->q, tempy->z));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return result;
    }

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_RATIONAL(ytype)) {
        MPZ_Object *tempx = NULL;
        MPQ_Object *tempy = NULL;

        if (!(tempx = GMPy_MPZ_From_IntegerWithTypeAndCopy(x,xtype, context)) ||
            !(tempy = GMPy_MPQ_From_RationalWithTypeAndCopy(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpz_abs(tempx->z, tempx->z);
        mpq_abs(tempy->q, tempy->q);

        result = _return_negated_cmp(mpq_cmp_z(tempy->q, tempx->z));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return result;
    }

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype)) {
        MPQ_Object *tempx = NULL, *tempy = NULL;

        if (!(tempx = GMPy_MPQ_From_RationalWithTypeAndCopy(x, xtype, context)) ||
            !(tempy = GMPy_MPQ_From_RationalWithTypeAndCopy(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpq_abs(tempx->q, tempx->q);
        mpq_abs(tempy->q, tempy->q);

        result = _return_cmp(mpq_cmp(tempx->q, tempy->q));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return result;
    }

    /* We perform exact comparisons between the mpz, mpq, and mpfr types.
     */

    if (IS_TYPE_REAL(xtype) && IS_TYPE_INTEGER(ytype)) {
            MPFR_Object *tempx = NULL;
            MPZ_Object *tempy = NULL;

        if (!(tempx = GMPy_MPFR_From_RealWithTypeAndCopy(x, xtype, 1, context)) ||
            !(tempy = GMPy_MPZ_From_IntegerWithTypeAndCopy(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpfr_clear_flags();
        mpfr_abs(tempx->f, tempx->f, MPFR_RNDN);
        mpz_abs(tempy->z, tempy->z);

        result =_return_cmp(mpfr_cmp_z(tempx->f, tempy->z));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        GMPY_CHECK_ERANGE(result, context, "invalid comparison with NaN");
        return result;
    }

    if (IS_TYPE_REAL(xtype) && IS_TYPE_RATIONAL(ytype)) {
            MPFR_Object *tempx = NULL;
            MPQ_Object *tempy = NULL;

        if (!(tempx = GMPy_MPFR_From_RealWithTypeAndCopy(x, xtype, 1, context)) ||
            !(tempy = GMPy_MPQ_From_RationalWithTypeAndCopy(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpfr_clear_flags();
        mpfr_abs(tempx->f, tempx->f, MPFR_RNDN);
        mpq_abs(tempy->q, tempy->q);

        result =_return_cmp(mpfr_cmp_q(tempx->f, tempy->q));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        GMPY_CHECK_ERANGE(result, context, "invalid comparison with NaN");
        return result;
    }

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype)) {
            MPFR_Object *tempx = NULL;
            MPFR_Object *tempy = NULL;

        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context)) ||
            !(tempy = GMPy_MPFR_From_RealWithType(y, ytype, 1, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpfr_clear_flags();
        result =_return_cmp(mpfr_cmpabs(tempx->f, tempy->f));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        GMPY_CHECK_ERANGE(result, context, "invalid comparison with NaN");
        return result;
    }

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_REAL(ytype)) {
            MPZ_Object *tempx = NULL;
            MPFR_Object *tempy = NULL;

        if (!(tempx = GMPy_MPZ_From_IntegerWithTypeAndCopy(x, xtype, context)) ||
            !(tempy = GMPy_MPFR_From_RealWithTypeAndCopy(y, ytype, 1, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpfr_clear_flags();
        mpz_abs(tempx->z, tempx->z);
        mpfr_abs(tempy->f, tempy->f, MPFR_RNDN);

        result =_return_negated_cmp(mpfr_cmp_z(tempy->f, tempx->z));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        GMPY_CHECK_ERANGE(result, context, "invalid comparison with NaN");
        return result;
    }

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_REAL(ytype)) {
            MPQ_Object *tempx = NULL;
            MPFR_Object *tempy = NULL;

        if (!(tempx = GMPy_MPQ_From_RationalWithTypeAndCopy(x, xtype, context)) ||
            !(tempy = GMPy_MPFR_From_RealWithTypeAndCopy(y, xtype, 1, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpfr_clear_flags();
        mpq_abs(tempx->q, tempx->q);
        mpfr_abs(tempy->f, tempy->f, MPFR_RNDN);

        result =_return_negated_cmp(mpfr_cmp_q(tempy->f, tempx->q));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        GMPY_CHECK_ERANGE(result, context, "invalid comparison with NaN");
        return result;
    }


    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype)) {
            MPC_Object *tempx = NULL;
            MPC_Object *tempy = NULL;

        if (!(tempx = GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context)) ||
            !(tempy = GMPy_MPC_From_ComplexWithType(y, ytype, 1, 1, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpfr_clear_flags();
        result =_return_cmp(mpc_cmp_abs(tempx->c, tempy->c));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        GMPY_CHECK_ERANGE(result, context, "invalid comparison with NaN");
        return result;
    }

    TYPE_ERROR("cmp_abs() requires integer, rational, real, or complex arguments");
    return NULL;
}

