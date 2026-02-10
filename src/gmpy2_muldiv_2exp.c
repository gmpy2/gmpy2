/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_muldiv_2exp.c                                                     *
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
GMPy_Real_Mul_2exp(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *result, *tempx;
    long exp;
    int is_signed = 0;

    CHECK_CONTEXT(context);

    exp = GMPy_Integer_AsUnsignedLongOrLong(y, &is_signed);
    if (exp == -1 && PyErr_Occurred()) {
        return NULL;
    }

    result = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPFR_From_Real(x, 1, context);
    if (!result || !tempx) {
        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_clear_flags();

    if (is_signed) {
        result->rc = mpfr_mul_2si(result->f, tempx->f, exp,
                                  GET_MPFR_ROUND(context));
    }
    else {
        result->rc = mpfr_mul_2ui(result->f, tempx->f, (unsigned long)exp,
                                  GET_MPFR_ROUND(context));
    }
    Py_DECREF((PyObject*)tempx);
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Complex_Mul_2exp(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPC_Object *result, *tempx;
    long exp;
    int is_signed = 0;

    CHECK_CONTEXT(context);

    exp = GMPy_Integer_AsUnsignedLongOrLong(y, &is_signed);
    if (exp == -1 && PyErr_Occurred()) {
        return NULL;
    }

    result = GMPy_MPC_New(0, 0, context);
    tempx = GMPy_MPC_From_Complex(x, 1, 1, context);
    if (!result || !tempx) {
        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (is_signed) {
        result->rc = mpc_mul_2si(result->c, tempx->c, exp,
                                 GET_MPC_ROUND(context));
    }
    else {
        result->rc = mpc_mul_2ui(result->c, tempx->c, exp,
                                 GET_MPC_ROUND(context));
    }
    Py_DECREF((PyObject*)tempx);
    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_context_mul_2exp,
"mul_2exp($self, x, n, /)\n--\n\n"
"Return `mpfr` or `mpc` multiplied by 2**n.");

PyDoc_STRVAR(GMPy_doc_function_mul_2exp,
"mul_2exp($module, x, n, /)\n--\n\n"
"Return x multiplied by 2**n.");

static PyObject *
GMPy_Number_Mul_2exp(PyObject *x, PyObject *y, CTXT_Object *context)
{
    if (IS_REAL(x))
        return GMPy_Real_Mul_2exp(x, y, context);

    if (IS_COMPLEX(x))
        return GMPy_Complex_Mul_2exp(x, y, context);

    TYPE_ERROR("mul_2exp() argument type not supported");
    return NULL;
}

static PyObject *
GMPy_Context_Mul_2exp(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("mul_2exp() requires 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Mul_2exp(PyTuple_GET_ITEM(args, 0),
                                PyTuple_GET_ITEM(args, 1),
                                context);
}

/* ======================================================================= */

static PyObject *
GMPy_Real_Div_2exp(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *result, *tempx;
    long exp;
    int is_signed = 0;

    CHECK_CONTEXT(context);

    exp = GMPy_Integer_AsUnsignedLongOrLong(y, &is_signed);
    if (exp == -1 && PyErr_Occurred()) {
        return NULL;
    }

    result = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPFR_From_Real(x, 1, context);
    if (!result || !tempx) {
        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_clear_flags();

    if (is_signed) {
        result->rc = mpfr_div_2si(result->f, tempx->f, exp,
                                  GET_MPFR_ROUND(context));
    }
    else  {
        result->rc = mpfr_div_2ui(result->f, tempx->f, exp,
                                  GET_MPFR_ROUND(context));
    }
    Py_DECREF((PyObject*)tempx);
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Complex_Div_2exp(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPC_Object *result, *tempx;
    long exp;
    int is_signed = 0;

    CHECK_CONTEXT(context);

    exp = GMPy_Integer_AsUnsignedLongOrLong(y, &is_signed);
    if (exp == -1 && PyErr_Occurred()) {
        return NULL;
    }

    result = GMPy_MPC_New(0, 0, context);
    tempx = GMPy_MPC_From_Complex(x, 1, 1, context);
    if (!result || !tempx) {
        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (is_signed) {
        result->rc = mpc_div_2si(result->c, tempx->c, exp,
                                 GET_MPC_ROUND(context));
    }
    else {
        result->rc = mpc_div_2ui(result->c, tempx->c, exp,
                                 GET_MPC_ROUND(context));
    }
    Py_DECREF((PyObject*)tempx);
    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_context_div_2exp,
"div_2exp($self, x, n, /)\n--\n\n"
"Return `mpfr` or `mpc` divided by 2**n.");

PyDoc_STRVAR(GMPy_doc_function_div_2exp,
"div_2exp($module, x, n, /)\n--\n\n"
"Return x divided by 2**n.");

static PyObject *
GMPy_Number_Div_2exp(PyObject *x, PyObject *y, CTXT_Object *context)
{
    if (IS_REAL(x))
        return GMPy_Real_Div_2exp(x, y, context);

    if (IS_COMPLEX(x))
        return GMPy_Complex_Div_2exp(x, y, context);

    TYPE_ERROR("div_2exp() argument type not supported");
    return NULL;
}

static PyObject *
GMPy_Context_Div_2exp(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("div_2exp() requires 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Div_2exp(PyTuple_GET_ITEM(args, 0),
                                PyTuple_GET_ITEM(args, 1),
                                context);
}
