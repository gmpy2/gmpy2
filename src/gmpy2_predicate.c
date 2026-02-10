/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_predicate.c                                                       *
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


/* To work with the MPC_IS_ macros, NAN, INF, FINITE, and ZERO are
 * all upper-case.
 */

PyDoc_STRVAR(GMPy_doc_function_is_nan,
"is_nan($module, x, /)\n--\n\n"
"Return `True` if x is NaN (Not-A-Number) else `False`.");

PyDoc_STRVAR(GMPy_doc_context_is_nan,
"is_nan($self, x, /)\n--\n\n"
"Return `True` if x is NaN (Not-A-Number) else `False`.");

PyDoc_STRVAR(GMPy_doc_method_is_nan,
"is_nan($self)\n--\n\n"
"Return `True` if self is NaN (Not-A-Number) else `False`.");

static PyObject *
GMPy_RealWithType_Is_NAN(PyObject *x, int xtype, CTXT_Object *context)
{
    MPFR_Object *tempx = NULL;
    int res = 0;

    if (IS_TYPE_MPFR(xtype)) {
        res = mpfr_nan_p(MPFR(x));
    }
    else {
        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context)))
            return NULL;
        res = mpfr_nan_p(tempx->f);
        Py_DECREF((PyObject*)tempx);
    }

    if (res)
        Py_RETURN_TRUE;
    else 
        Py_RETURN_FALSE;
}

static PyObject *
GMPy_ComplexWithType_Is_NAN(PyObject *x, int xtype, CTXT_Object *context)
{
    MPC_Object *tempx = NULL;
    int res = 0;

    if (IS_TYPE_MPC(xtype)) {
        res = MPC_IS_NAN_P(x);
    }
    else {
        if (!(tempx = GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context)))
            return NULL;
        res = MPC_IS_NAN_P(tempx);
        Py_DECREF(tempx);
    }

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

/* This creates the following functions:
 *   GMPy_Number_Is_NAN(x, context)
 *   GMPy_Context_Is_NAN(self, other)
 *   GMPy_Number_Method_Is_NAN(self, args)
 *
 * They assume the following functions have been defined above:
 *   GMPy_RealWithType_Is_NAN(x, xtype, context)
 *   GMPy_ComplexWithType_Is_NAN(x, xtype, context)
 */

GMPY_MPFR_MPC_UNIOP_TEMPLATEWT(Is_NAN, is_nan)

PyDoc_STRVAR(GMPy_doc_function_is_infinite,
"is_infinite($module, x, /)\n--\n\n"
"Return `True` if x is +Infinity or -Infinity. If x is an `mpc`, return `True`\n"
"if either x.real or x.imag is infinite. Otherwise return `False`.");

PyDoc_STRVAR(GMPy_doc_context_is_infinite,
"is_infinite($self, x, /)\n--\n\n"
"Return `True` if x is +Infinity or -Infinity. If x is an `mpc`, return `True`\n"
"if either x.real or x.imag is infinite. Otherwise return `False`.");

PyDoc_STRVAR(GMPy_doc_method_is_infinite,
"is_infinite($self)\n--\n\n"
"Return `True` if self is +Infinity or -Infinity. If self is an `mpc`, return `True`\n"
"if either self.real or self.imag is infinite. Otherwise return `False`.");

static PyObject *
GMPy_RealWithType_Is_Infinite(PyObject *x, int xtype, CTXT_Object *context)
{
    MPFR_Object *tempx = NULL;
    int res = 0;

    if (IS_TYPE_MPFR(xtype)) {
        res = mpfr_inf_p(MPFR(x));
    }
    else {
        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context)))
            return NULL;
        res = mpfr_inf_p(tempx->f);
        Py_DECREF(tempx);
    }

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

static PyObject *
GMPy_ComplexWithType_Is_Infinite(PyObject *x, int xtype, CTXT_Object *context)
{
    MPC_Object *tempx = NULL;
    int res = 0;

    if (IS_TYPE_MPC(xtype)) {
        res = MPC_IS_INF_P(x);
    }
    else {
        if (!(tempx = GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context)))
            return NULL;
        res = MPC_IS_INF_P(tempx);
        Py_DECREF(tempx);
    }

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

GMPY_MPFR_MPC_UNIOP_TEMPLATEWT(Is_Infinite, is_infinite)

PyDoc_STRVAR(GMPy_doc_function_is_finite,
"is_finite($module, x, /)\n--\n\n"
"Return `True` if x is an actual number (i.e. non NaN or Infinity). If x is\n"
"an `mpc`, return `True` if both x.real and x.imag are finite.");

PyDoc_STRVAR(GMPy_doc_context_is_finite,
"is_finite($self, x, /)\n--\n\n"
"Return `True` if x is an actual number (i.e. non NaN or Infinity). If x is\n"
"an `mpc`, return `True` if both x.real and x.imag are finite.");

PyDoc_STRVAR(GMPy_doc_method_is_finite,
"is_finite($self)\n--\n\n"
"Return `True` if self is an actual number (i.e. non NaN or Infinity). If self is\n"
"an `mpc`, return `True` if both self.real and self.imag are finite.");

static PyObject *
GMPy_RealWithType_Is_Finite(PyObject *x, int xtype, CTXT_Object *context)
{
    MPFR_Object *tempx = NULL;
    int res = 0;

    if (IS_TYPE_MPFR(xtype)) {
        res = mpfr_number_p(MPFR(x));
    }
    else {
        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context)))
            return NULL;
        res = mpfr_number_p(tempx->f);
        Py_DECREF(tempx);
    }

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

static PyObject *
GMPy_ComplexWithType_Is_Finite(PyObject *x, int xtype, CTXT_Object *context)
{
    MPC_Object *tempx = NULL;
    int res = 0;

    if (IS_TYPE_MPC(xtype)) {
        res = MPC_IS_FINITE_P(x);
    }
    else {
        if (!(tempx = GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context)))
            return NULL;
        res = MPC_IS_FINITE_P(tempx);
        Py_DECREF(tempx);
    }

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

GMPY_MPFR_MPC_UNIOP_TEMPLATEWT(Is_Finite, is_finite)

PyDoc_STRVAR(GMPy_doc_function_is_zero,
"is_zero($module, x, /)\n--\n\n"
"Return `True` if x is equal to 0. If x is an `mpc`, return `True` if both x.real\n"
"and x.imag are equal to 0.");

PyDoc_STRVAR(GMPy_doc_context_is_zero,
"is_zero($self, x, /)\n--\n\n"
"Return `True` if x is equal to 0. If x is an `mpc`, return `True` if both x.real\n"
"and x.imag are equal to 0.");

PyDoc_STRVAR(GMPy_doc_method_is_zero,
"is_zero($self)\n--\n\n"
"Return `True` if self is equal to 0. If self is an `mpc`, return `True` if both self.real\n"
"and self.imag are equal to 0.");

static PyObject *
GMPy_RealWithType_Is_Zero(PyObject *x, int xtype, CTXT_Object *context)
{
    MPFR_Object *tempx = NULL;
    int res = 0;

    if (IS_TYPE_MPFR(xtype)) {
        res = mpfr_zero_p(MPFR(x));
    }
    else {
        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context)))
            return NULL;
        res = mpfr_zero_p(tempx->f);
        Py_DECREF(tempx);
    }

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

static PyObject *
GMPy_ComplexWithType_Is_Zero(PyObject *x, int xtype, CTXT_Object *context)
{
    MPC_Object *tempx = NULL;
    int res = 0;

    if (IS_TYPE_MPC(xtype)) {
        res = MPC_IS_ZERO_P(x);
    }
    else {
        if (!(tempx = GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context)))
            return NULL;
        res = MPC_IS_ZERO_P(tempx);
        Py_DECREF(tempx);
    }

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

GMPY_MPFR_MPC_UNIOP_TEMPLATEWT(Is_Zero, is_zero)

PyDoc_STRVAR(GMPy_doc_function_is_signed,
"is_signed($module, x, /)\n--\n\n"
"Return `True` if the sign bit of x is set.");

PyDoc_STRVAR(GMPy_doc_context_is_signed,
"is_signed($self, x, /)\n--\n\n"
"Return `True` if the sign bit of x is set.");

PyDoc_STRVAR(GMPy_doc_method_is_signed,
"is_signed($self)\n--\n\n"
"Return `True` if the sign bit of self is set.");

static PyObject *
GMPy_RealWithType_Is_Signed(PyObject *x, int xtype, CTXT_Object *context)
{
    MPFR_Object *tempx = NULL;
    int res;

    if (IS_TYPE_MPFR(xtype)) {
        res = mpfr_signbit(MPFR(x));
    }
    else {
        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context))) {
            return NULL;
        }
        res = mpfr_signbit(tempx->f);
        Py_DECREF((PyObject*)tempx);
    }

    if (res) {
        Py_RETURN_TRUE;
    }
    else {
        Py_RETURN_FALSE;
    }
}

static PyObject *
GMPy_MPFR_Is_Signed_Method(PyObject *self, PyObject *args)
{
    return GMPy_Number_Is_Signed(self, NULL);
}

GMPY_MPFR_UNIOP_TEMPLATEWT(Is_Signed, is_signed)

PyDoc_STRVAR(GMPy_doc_function_is_regular,
"is_regular($module, x, /)\n--\n\n"
"Return `True` if x is not zero, NaN, or Infinity; `False` otherwise.");

PyDoc_STRVAR(GMPy_doc_context_is_regular,
"is_regular($self, x, /)\n--\n\n"
"Return `True` if x is not zero, NaN, or Infinity; `False` otherwise.");

PyDoc_STRVAR(GMPy_doc_method_is_regular,
"is_regular($self)\n--\n\n"
"Return `True` if self is not zero, NaN, or Infinity; `False` otherwise.");

static PyObject *
GMPy_RealWithType_Is_Regular(PyObject *x, int xtype, CTXT_Object *context)
{
    MPFR_Object *tempx = NULL;
    int res;

    if (IS_TYPE_MPFR(xtype)) {
        res = mpfr_regular_p(MPFR(x));
    }
    else {
        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context))) {
            return NULL;
        }
        res = mpfr_regular_p(tempx->f);
        Py_DECREF((PyObject*)tempx);
    }

    if (res) {
        Py_RETURN_TRUE;
    }
    else {
        Py_RETURN_FALSE;
    }
}

static PyObject *
GMPy_MPFR_Is_Regular_Method(PyObject *self, PyObject *args)
{
    return GMPy_Number_Is_Regular(self, NULL);
}

GMPY_MPFR_UNIOP_TEMPLATEWT(Is_Regular, is_regular)

PyDoc_STRVAR(GMPy_doc_function_is_integer,
"is_integer($module, x, /)\n--\n\n"
"Return `True` if x is an integer; `False` otherwise.");

PyDoc_STRVAR(GMPy_doc_context_is_integer,
"is_integer($self, x, /)\n--\n\n"
"Return `True` if x is an integer; `False` otherwise.");

PyDoc_STRVAR(GMPy_doc_method_is_integer,
"is_integer($self)\n--\n\n"
"Return `True` if self is an integer; `False` otherwise.");

static PyObject *
GMPy_RealWithType_Is_Integer(PyObject *x, int xtype, CTXT_Object *context)
{
    MPFR_Object *tempx = NULL;
    int res;

    if (IS_TYPE_MPFR(xtype)) {
        res = mpfr_integer_p(MPFR(x));
    }
    else {
        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype,  1, context))) {
            return NULL;
        }
        res = mpfr_integer_p(tempx->f);
        Py_DECREF((PyObject*)tempx);
    }

    if (res) {
        Py_RETURN_TRUE;
    }
    else {
        Py_RETURN_FALSE;
    }
}

static PyObject *
GMPy_MPFR_Is_Integer_Method(PyObject *self, PyObject *args)
{
    return GMPy_Number_Is_Integer(self, NULL);
}

GMPY_MPFR_UNIOP_TEMPLATEWT(Is_Integer, is_integer)

PyDoc_STRVAR(GMPy_doc_function_is_lessgreater,
"is_lessgreater($module, x, y, /)\n--\n\n"
"Return `True` if x > y or x < y. Return `False` if x == y or either x\n"
"and/or y is NaN.");

static PyObject *
GMPy_RealWithType_Is_LessGreater(PyObject *x, int xtype, PyObject *y, int ytype, CTXT_Object *context)
{
    MPFR_Object *tempx = NULL, *tempy = NULL;
    int res;

    tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context);
    tempy = GMPy_MPFR_From_RealWithType(y, ytype, 1, context);
    if (!tempx || !tempy) {
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        return NULL;
    }

    res = mpfr_lessgreater_p(tempx->f, tempy->f);
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

GMPY_MPFR_BINOP_TEMPLATEWT(Is_LessGreater, is_lessgreater)

PyDoc_STRVAR(GMPy_doc_function_is_unordered,
"is_unordered($module, x, y, /)\n--\n\n"
"Return `True` if either x and/or y is NaN.");

static PyObject *
GMPy_RealWithType_Is_Unordered(PyObject *x, int xtype, PyObject *y, int ytype, CTXT_Object *context)
{
    MPFR_Object *tempx = NULL, *tempy = NULL;
    int res;

    tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context);
    tempy = GMPy_MPFR_From_RealWithType(y, ytype, 1, context);
    if (!tempx || !tempy) {
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        return NULL;
    }

    res = mpfr_unordered_p(tempx->f, tempy->f);
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

GMPY_MPFR_BINOP_TEMPLATEWT(Is_Unordered, is_unordered)
