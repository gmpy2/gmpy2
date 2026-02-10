/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_math.c                                                            *
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

PyDoc_STRVAR(GMPy_doc_context_sin,
"sin($self, x, /)\n--\n\n"
"Return sine of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_sin,
"sin($module, x, /)\n--\n\n"
"Return sine of x; x in radians.");

GMPY_MPFR_MPC_UNIOP_EXWT(Sin, sin)

PyDoc_STRVAR(GMPy_doc_context_cos,
"cos($self, x, /)\n--\n\n"
"Return cosine of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_cos,
"cos($module, x, /)\n--\n\n"
"Return cosine of x; x in radians.");

GMPY_MPFR_MPC_UNIOP_EXWT(Cos, cos)

PyDoc_STRVAR(GMPy_doc_context_tan,
"tan($self, x, /)\n--\n\n"
"Return tangent of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_tan,
"tan($module, x, /)\n--\n\n"
"Return tangent of x; x in radians.");

GMPY_MPFR_MPC_UNIOP_EXWT(Tan, tan)

PyDoc_STRVAR(GMPy_doc_context_atan,
"atan($self, x, /)\n--\n\n"
"Return inverse tangent of x; result in radians.");

PyDoc_STRVAR(GMPy_doc_function_atan,
"atan($module, x, /)\n--\n\n"
"Return inverse tangent of x; result in radians.");

GMPY_MPFR_MPC_UNIOP_EXWT(Atan, atan)

PyDoc_STRVAR(GMPy_doc_context_sinh,
"sinh($self, x, /)\n--\n\n"
"Return hyperbolic sine of x.");

PyDoc_STRVAR(GMPy_doc_function_sinh,
"sinh($module, x, /)\n--\n\n"
"Return hyperbolic sine of x.");

GMPY_MPFR_MPC_UNIOP_EXWT(Sinh, sinh)

PyDoc_STRVAR(GMPy_doc_context_cosh,
"cosh($self, x, /)\n--\n\n"
"Return hyperbolic cosine of x.");

PyDoc_STRVAR(GMPy_doc_function_cosh,
"cosh($module, x, /)\n--\n\n"
"Return hyperbolic cosine of x.");

GMPY_MPFR_MPC_UNIOP_EXWT(Cosh, cosh)

PyDoc_STRVAR(GMPy_doc_context_tanh,
"tanh($self, x, /)\n--\n\n"
"Return hyperbolic tangent of x.");

PyDoc_STRVAR(GMPy_doc_function_tanh,
"tanh($module, x, /)\n--\n\n"
"Return hyperbolic tangent of x.");

GMPY_MPFR_MPC_UNIOP_EXWT(Tanh, tanh)

PyDoc_STRVAR(GMPy_doc_context_asinh,
"asinh($self, x, /)\n--\n\n"
"Return inverse hyperbolic sine of x.");

PyDoc_STRVAR(GMPy_doc_function_asinh,
"asinh($module, x, /)\n--\n\n"
"Return inverse hyperbolic sine of x.");

GMPY_MPFR_MPC_UNIOP_EXWT(Asinh, asinh)

PyDoc_STRVAR(GMPy_doc_context_acosh,
"acosh($self, x, /)\n--\n\n"
"Return inverse hyperbolic cosine of x.");

PyDoc_STRVAR(GMPy_doc_function_acosh,
"acosh($module, x, /)\n--\n\n"
"Return inverse hyperbolic cosine of x.");

GMPY_MPFR_MPC_UNIOP_EXWT(Acosh, acosh)

/* Section 2:
 * These functions accept a single argument and return an mpfr result.
 *
 * GMPY_MPFR_UNIOP(NAME, FUNC) creates the following functions:
 *     GMPy_Real_NAME(x, context)
 *     GMPy_Number_NAME(x, context)
 *     GMPy_Context_NAME(self, other)
 *     - called with METH_O
 */

PyDoc_STRVAR(GMPy_doc_context_sec,
"sec($self, x, /)\n--\n\n"
"Return secant of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_sec,
"sec($module, x, /)\n--\n\n"
"Return secant of x; x in radians.");

GMPY_MPFR_UNIOP_EXWT(Sec, sec)

PyDoc_STRVAR(GMPy_doc_context_csc,
"csc($self, x, /)\n--\n\n"
"Return cosecant of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_csc,
"csc($module, x, /)\n--\n\n"
"Return cosecant of x; x in radians.");

GMPY_MPFR_UNIOP_EXWT(Csc, csc)

PyDoc_STRVAR(GMPy_doc_context_cot,
"cot($self, x, /)\n--\n\n"
"Return cotangent of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_cot,
"cot($module, x, /)\n--\n\n"
"Return cotangent of x; x in radians.");

GMPY_MPFR_UNIOP_EXWT(Cot, cot)

PyDoc_STRVAR(GMPy_doc_context_sech,
"sech($self, x, /)\n--\n\n"
"Return hyperbolic secant of x.");

PyDoc_STRVAR(GMPy_doc_function_sech,
"sech($module, x, /)\n--\n\n"
"Return hyperbolic secant of x.");

GMPY_MPFR_UNIOP_EXWT(Sech, sech)

PyDoc_STRVAR(GMPy_doc_context_csch,
"csch($self, x, /)\n--\n\n"
"Return hyperbolic cosecant of x.");

PyDoc_STRVAR(GMPy_doc_function_csch,
"csch($module, x, /)\n--\n\n"
"Return hyperbolic cosecant of x.");

GMPY_MPFR_UNIOP_EXWT(Csch, csch)

PyDoc_STRVAR(GMPy_doc_context_coth,
"coth($self, x, /)\n--\n\n"
"Return hyperbolic cotangent of x.");

PyDoc_STRVAR(GMPy_doc_function_coth,
"coth($module, x, /)\n--\n\n"
"Return hyperbolic cotangent of x.");

GMPY_MPFR_UNIOP_EXWT(Coth, coth)

PyDoc_STRVAR(GMPy_doc_context_rec_sqrt,
"rec_sqrt($self, x, /)\n--\n\n"
"Return the reciprocal of the square root of x.");

PyDoc_STRVAR(GMPy_doc_function_rec_sqrt,
"rec_sqrt($module, x, /)\n--\n\n"
"Return the reciprocal of the square root of x.");

GMPY_MPFR_UNIOP_EXWT(RecSqrt, rec_sqrt)

PyDoc_STRVAR(GMPy_doc_context_rint,
"rint($self, x, /)\n--\n\n"
"Return x rounded to the nearest integer using the context rounding\n"
"mode.");

PyDoc_STRVAR(GMPy_doc_function_rint,
"rint($module, x, /)\n--\n\n"
"Return x rounded to the nearest integer using the current rounding\n"
"mode.");

GMPY_MPFR_UNIOP_EXWT(Rint, rint)

PyDoc_STRVAR(GMPy_doc_context_rint_ceil,
"rint_ceil($self, x, /)\n--\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"next higher or equal integer and then, if needed, using the context\n"
"rounding mode.");

PyDoc_STRVAR(GMPy_doc_function_rint_ceil,
"rint_ceil($module, x, /)\n--\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"next higher or equal integer and then, if needed, using the current\n"
"rounding mode.");

GMPY_MPFR_UNIOP_EXWT(RintCeil, rint_ceil)

PyDoc_STRVAR(GMPy_doc_context_rint_floor,
"rint_floor($self, x, /)\n--\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"next lower or equal integer and then, if needed, using the context\n"
"rounding mode.");

PyDoc_STRVAR(GMPy_doc_function_rint_floor,
"rint_floor($module, x, /)\n--\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"next lower or equal integer and then, if needed, using the current\n"
"rounding mode.");

GMPY_MPFR_UNIOP_EXWT(RintFloor, rint_floor)

PyDoc_STRVAR(GMPy_doc_context_rint_round,
"rint_round($self, x, /)\n--\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"nearest integer (ties away from 0) and then, if needed, using\n"
"the context rounding mode.");

PyDoc_STRVAR(GMPy_doc_function_rint_round,
"rint_round($module, x, /)\n--\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"nearest integer (ties away from 0) and then, if needed, using\n"
"the current rounding mode.");

GMPY_MPFR_UNIOP_EXWT(RintRound, rint_round)

PyDoc_STRVAR(GMPy_doc_context_rint_trunc,
"rint_trunc($self, x, /)\n--\n\n"
"Return x rounded to the nearest integer by first rounding towards\n"
"zero and then, if needed, using the context rounding mode.");

PyDoc_STRVAR(GMPy_doc_function_rint_trunc,
"rint_trunc($module, x, /)\n--\n\n"
"Return x rounded to the nearest integer by first rounding towards\n"
"zero and then, if needed, using the current rounding mode.");

GMPY_MPFR_UNIOP_EXWT(RintTrunc, rint_trunc)

PyDoc_STRVAR(GMPy_doc_context_frac,
"frac($self, x, /)\n--\n\n"
"Return fractional part of x.");

PyDoc_STRVAR(GMPy_doc_function_frac,
"frac($module, x, /)\n--\n\n"
"Return fractional part of x.");

GMPY_MPFR_UNIOP_EXWT(Frac, frac)

PyDoc_STRVAR(GMPy_doc_context_cbrt,
"cbrt($self, x, /)\n--\n\n"
"Return the cube root of x.");

PyDoc_STRVAR(GMPy_doc_function_cbrt,
"cbrt($module, x, /)\n--\n\n"
"Return the cube root of x.");

GMPY_MPFR_UNIOP_EXWT(Cbrt, cbrt)

PyDoc_STRVAR(GMPy_doc_context_log2,
"log2($self, x, /)\n--\n\n"
"Return base-2 logarithm of x.");

PyDoc_STRVAR(GMPy_doc_function_log2,
"log2($module, x, /)\n--\n\n"
"Return base-2 logarithm of x.");

GMPY_MPFR_UNIOP_EXWT(Log2, log2)

PyDoc_STRVAR(GMPy_doc_context_exp2,
"exp2($self, x, /)\n--\n\n"
"Return 2**x.");

PyDoc_STRVAR(GMPy_doc_function_exp2,
"exp2($module, x, /)\n--\n\n"
"Return 2**x.");

GMPY_MPFR_UNIOP_EXWT(Exp2, exp2)

PyDoc_STRVAR(GMPy_doc_context_exp10,
"exp10($self, x, /)\n--\n\n"
"Return 10**x.");

PyDoc_STRVAR(GMPy_doc_function_exp10,
"exp10($module, x, /)\n--\n\n"
"Return 10**x.");

GMPY_MPFR_UNIOP_EXWT(Exp10, exp10)

PyDoc_STRVAR(GMPy_doc_context_log1p,
"log1p($self, x, /)\n--\n\n"
"Return natural logarithm of (1+x).");

PyDoc_STRVAR(GMPy_doc_function_log1p,
"log1p($module, x, /)\n--\n\n"
"Return natural logarithm of (1+x).");

GMPY_MPFR_UNIOP_EXWT(Log1p, log1p)

PyDoc_STRVAR(GMPy_doc_context_expm1,
"expm1($self, x, /)\n--\n\n"
"Return exp(x) - 1.");

PyDoc_STRVAR(GMPy_doc_function_expm1,
"expm1($module, x, /)\n--\n\n"
"Return exp(x) - 1.");

GMPY_MPFR_UNIOP_EXWT(Expm1, expm1)

PyDoc_STRVAR(GMPy_doc_context_eint,
"eint($self, x, /)\n--\n\n"
"Return exponential integral of x.");

PyDoc_STRVAR(GMPy_doc_function_eint,
"eint($module, x, /)\n--\n\n"
"Return exponential integral of x.");

GMPY_MPFR_UNIOP_EXWT(Eint, eint)

PyDoc_STRVAR(GMPy_doc_context_li2,
"li2($self, x, /)\n--\n\n"
"Return real part of dilogarithm of x.");

PyDoc_STRVAR(GMPy_doc_function_li2,
"li2($module, x, /)\n--\n\n"
"Return real part of dilogarithm of x.");

GMPY_MPFR_UNIOP_EXWT(Li2, li2)

PyDoc_STRVAR(GMPy_doc_context_gamma,
"gamma($self, x, /)\n--\n\n"
"Return gamma of x.");

PyDoc_STRVAR(GMPy_doc_function_gamma,
"gamma($module, x, /)\n--\n\n"
"Return gamma of x.");

GMPY_MPFR_UNIOP_EXWT(Gamma, gamma)

PyDoc_STRVAR(GMPy_doc_function_gamma_inc,
"gamma_inc($module, a, x, /)\n--\n\n\
Return (upper) incomplete gamma of a and x.");

PyDoc_STRVAR(GMPy_doc_context_gamma_inc,
"gamma_inc($self, a, x, /)\n--\n\n\
Return (upper) incomplete gamma of a and x.");

GMPY_MPFR_BINOP_EXWT(Gamma_Inc, gamma_inc)

PyDoc_STRVAR(GMPy_doc_context_lngamma,
"lngamma($self, x, /)\n--\n\n"
"Return natural logarithm of gamma(x).");

PyDoc_STRVAR(GMPy_doc_function_lngamma,
"lngamma($module, x, /)\n--\n\n"
"Return natural logarithm of gamma(x).");

GMPY_MPFR_UNIOP_EXWT(Lngamma, lngamma)

PyDoc_STRVAR(GMPy_doc_context_digamma,
"digamma($self, x, /)\n--\n\n"
"Return digamma of x.");

PyDoc_STRVAR(GMPy_doc_function_digamma,
"digamma($module, x, /)\n--\n\n"
"Return digamma of x.");

GMPY_MPFR_UNIOP_EXWT(Digamma, digamma)

PyDoc_STRVAR(GMPy_doc_context_zeta,
"zeta($self, x, /)\n--\n\n"
"Return Riemann zeta of x.");

PyDoc_STRVAR(GMPy_doc_function_zeta,
"zeta($module, x, /)\n--\n\n"
"Return Riemann zeta of x.");

GMPY_MPFR_UNIOP_EXWT(Zeta, zeta)

PyDoc_STRVAR(GMPy_doc_context_erf,
"erf($self, x, /)\n--\n\n"
"Return error function of x.");

PyDoc_STRVAR(GMPy_doc_function_erf,
"erf($module, x, /)\n--\n\n"
"Return error function of x.");

GMPY_MPFR_UNIOP_EXWT(Erf, erf)

PyDoc_STRVAR(GMPy_doc_context_erfc,
"erfc($self, x, /)\n--\n\n"
"Return complementary error function of x.");

PyDoc_STRVAR(GMPy_doc_function_erfc,
"erfc($module, x, /)\n--\n\n"
"Return complementary error function of x.");

GMPY_MPFR_UNIOP_EXWT(Erfc, erfc)

PyDoc_STRVAR(GMPy_doc_context_j0,
"j0($self, x, /)\n--\n\n"
"Return first kind Bessel function of order 0 of x.");

PyDoc_STRVAR(GMPy_doc_function_j0,
"j0($module, x, /)\n--\n\n"
"Return first kind Bessel function of order 0 of x.");

GMPY_MPFR_UNIOP_EXWT(J0, j0)

PyDoc_STRVAR(GMPy_doc_context_j1,
"j1($self, x, /)\n--\n\n"
"Return first kind Bessel function of order 1 of x.");

PyDoc_STRVAR(GMPy_doc_function_j1,
"j1($module, x, /)\n--\n\n"
"Return first kind Bessel function of order 1 of x.");

GMPY_MPFR_UNIOP_EXWT(J1, j1)

PyDoc_STRVAR(GMPy_doc_context_y0,
"y0($self, x, /)\n--\n\n"
"Return second kind Bessel function of order 0 of x.");

PyDoc_STRVAR(GMPy_doc_function_y0,
"y0($module, x, /)\n--\n\n"
"Return second kind Bessel function of order 0 of x.");

GMPY_MPFR_UNIOP_EXWT(Y0, y0)

PyDoc_STRVAR(GMPy_doc_context_y1,
"y1($self, x, /)\n--\n\n"
"Return second kind Bessel function of order 1 of x.");

PyDoc_STRVAR(GMPy_doc_function_y1,
"y1($module, x, /)\n--\n\n"
"Return second kind Bessel function of order 1 of x.");

GMPY_MPFR_UNIOP_EXWT(Y1, y1)

PyDoc_STRVAR(GMPy_doc_context_ai,
"ai($self, x, /)\n--\n\n"
"Return Airy function of x.");

PyDoc_STRVAR(GMPy_doc_function_ai,
"ai($module, x, /)\n--\n\n"
"Return Airy function of x.");

GMPY_MPFR_UNIOP_EXWT(Ai, ai)

/* Section 3:
 * The following functions may return an mpc result for certain mpfr arguments.
 * Since the expectional values vary between functions, the 'Real' and 'Complex'
 * functions do not use macros. However, they do use a macro to create the
 * higher-level functions.
 *
 * GMPY_MPFR_MPC_UNIOP_TEMPLATE_EXWT(NAME, FUNC) creates the following functions:
 *     GMPy_Number_NAME(x, context)
 *     - assumes GMPy_RealWithType_NAME & GMPy_ComplexWithType_NAME exist
 *     GMPy_Context_NAME(self, other)
 *     - called with METH_O
 */

PyDoc_STRVAR(GMPy_doc_context_acos,
"acos($self, x, /)\n--\n\n"
"Return inverse cosine of x; result in radians.");

PyDoc_STRVAR(GMPy_doc_function_acos,
"acos($module, x, /)\n--\n\n"
"Return inverse cosine of x; result in radians.");

/* Helper function assumes x is of type mpfr. */
static PyObject *
_GMPy_MPFR_Acos(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result = NULL;

    if (!mpfr_nan_p(MPFR(x)) &&
            (mpfr_cmp_si(MPFR(x), 1) > 0 || mpfr_cmp_si(MPFR(x), -1) < 0) &&
            context->ctx.allow_complex
       ) {
        return GMPy_ComplexWithType_Acos(x, OBJ_TYPE_MPFR, context);
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    mpfr_clear_flags();

    result->rc = mpfr_acos(result->f, MPFR(x), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

/* Helper function assumes x is of type mpc. */
static PyObject *
_GMPy_MPC_Acos(PyObject *x, CTXT_Object *context)
{
    MPC_Object *result = NULL;

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        return NULL;
    }

    result->rc = mpc_acos(result->c, MPC(x), GET_MPC_ROUND(context));
    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_RealWithType_Acos(PyObject *x, int xtype, CTXT_Object *context)
{
    MPFR_Object *tempx = NULL;
    PyObject *result = NULL;

    if (IS_TYPE_MPFR(xtype)) {
        return _GMPy_MPFR_Acos(x, context);
    }
    else {
        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context))) {
            return NULL;
        }
        result = _GMPy_MPFR_Acos((PyObject*)tempx, context);
        Py_DECREF(tempx);
        return result;
    }
}

static PyObject *
GMPy_ComplexWithType_Acos(PyObject *x, int xtype, CTXT_Object *context)
{
    MPC_Object *tempx = NULL;
    PyObject *result = NULL;

    if (IS_TYPE_MPC(xtype)) {
        return _GMPy_MPC_Acos(x, context);
    }
    else {
        if (!(tempx = GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context))) {
            return NULL;
        }
        result = _GMPy_MPC_Acos((PyObject*)tempx, context);
        Py_DECREF(tempx);
        return result;
    }
}

GMPY_MPFR_MPC_UNIOP_TEMPLATE_EXWT(Acos, acos)

PyDoc_STRVAR(GMPy_doc_context_asin,
"asin($self, x, /)\n--\n\n"
"Return inverse sine of x; result in radians.");

PyDoc_STRVAR(GMPy_doc_function_asin,
"asin($module, x, /)\n--\n\n"
"Return inverse sine of x; result in radians.");

static PyObject *
_GMPy_MPFR_Asin(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result = NULL;

    if (!mpfr_nan_p(MPFR(x)) &&
            (mpfr_cmp_si(MPFR(x), 1) > 0 || mpfr_cmp_si(MPFR(x), -1) < 0) &&
            context->ctx.allow_complex
       ) {
        return GMPy_ComplexWithType_Asin(x, OBJ_TYPE_MPFR, context);
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    mpfr_clear_flags();

    result->rc = mpfr_asin(result->f, MPFR(x), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
_GMPy_MPC_Asin(PyObject *x, CTXT_Object *context)
{
    MPC_Object *result = NULL;

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        return NULL;
    }

    result->rc = mpc_asin(result->c, MPC(x), GET_MPC_ROUND(context));
    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_RealWithType_Asin(PyObject *x, int xtype, CTXT_Object *context)
{
    MPFR_Object *tempx = NULL;
    PyObject *result = NULL;

    if (IS_TYPE_MPFR(xtype)) {
        return _GMPy_MPFR_Asin(x, context);
    }
    else {
        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context))) {
            return NULL;
        }
        result = _GMPy_MPFR_Asin((PyObject*)tempx, context);
        Py_DECREF(tempx);
        return result;
    }
}

static PyObject *
GMPy_ComplexWithType_Asin(PyObject *x, int xtype, CTXT_Object *context)
{
    MPC_Object *tempx = NULL;
    PyObject *result = NULL;

    if (IS_TYPE_MPC(xtype)) {
        return _GMPy_MPC_Asin(x, context);
    }
    else {
        if (!(tempx = GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context))) {
            return NULL;
        }
        result = _GMPy_MPC_Asin((PyObject*)tempx, context);
        Py_DECREF(tempx);
        return result;
    }
}

GMPY_MPFR_MPC_UNIOP_TEMPLATE_EXWT(Asin, asin)

PyDoc_STRVAR(GMPy_doc_context_atanh,
"atanh($self, x, /)\n--\n\n"
"Return inverse hyperbolic tanget of x.");

PyDoc_STRVAR(GMPy_doc_function_atanh,
"atanh($module, x, /)\n--\n\n"
"Return inverse hyperbolic tangent of x.");

static PyObject *
_GMPy_MPFR_Atanh(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result = NULL;

    if (!mpfr_nan_p(MPFR(x)) &&
            (mpfr_cmp_si(MPFR(x), 1) > 0 || mpfr_cmp_si(MPFR(x), -1) < 0) &&
            context->ctx.allow_complex
       ) {
        return GMPy_ComplexWithType_Atanh(x, OBJ_TYPE_MPFR, context);
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    mpfr_clear_flags();

    result->rc = mpfr_atanh(result->f, MPFR(x), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
_GMPy_MPC_Atanh(PyObject *x, CTXT_Object *context)
{
    MPC_Object *result = NULL;

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        return NULL;
    }

    result->rc = mpc_atanh(result->c, MPC(x), GET_MPC_ROUND(context));
    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_RealWithType_Atanh(PyObject *x, int xtype, CTXT_Object *context)
{
    MPFR_Object *tempx = NULL;
    PyObject *result = NULL;

    if (IS_TYPE_MPFR(xtype)) {
        return _GMPy_MPFR_Atanh(x, context);
    }
    else {
        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context))) {
            return NULL;
        }
        result = _GMPy_MPFR_Atanh((PyObject*)tempx, context);
        Py_DECREF(tempx);
        return result;
    }
}

static PyObject *
GMPy_ComplexWithType_Atanh(PyObject *x, int xtype, CTXT_Object *context)
{
    MPC_Object *tempx = NULL;
    PyObject *result = NULL;

    if (IS_TYPE_MPC(xtype)) {
        return _GMPy_MPC_Atanh(x, context);
    }
    else {
        if (!(tempx = GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context))) {
            return NULL;
        }
        result = _GMPy_MPC_Atanh((PyObject*)tempx, context);
        Py_DECREF(tempx);
        return result;
    }
}

GMPY_MPFR_MPC_UNIOP_TEMPLATE_EXWT(Atanh, atanh)

PyDoc_STRVAR(GMPy_doc_function_atan2,
"atan2($module, y, x, /)\n--\n\n"
"Return arc-tangent of (y/x); result in radians.");

PyDoc_STRVAR(GMPy_doc_context_atan2,
"atan2($self, y, x, /)\n--\n\n"
"Return arc-tangent of (y/x); result in radians.");

GMPY_MPFR_BINOP_EXWT(Atan2, atan2)

PyDoc_STRVAR(GMPy_doc_function_hypot,
"hypot($module, x, y, /)\n--\n\n"
"Return square root of (x**2 + y**2).");

PyDoc_STRVAR(GMPy_doc_context_hypot,
"hypot($self, x, y, /)\n--\n\n"
"Return square root of (x**2 + y**2).");

GMPY_MPFR_BINOP_EXWT(Hypot, hypot)

static PyObject *
GMPy_RealWithType_Sin_Cos(PyObject *x, int xtype, CTXT_Object *context)
{
    MPFR_Object *s = NULL, *c = NULL, *tempx = NULL;
    PyObject *result = NULL;
    int code;

    s = GMPy_MPFR_New(0, context);
    c = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context);
    result = PyTuple_New(2);
    if (!s || !c || !tempx || !result) {
        Py_XDECREF(s);
        Py_XDECREF(c);
        Py_XDECREF(tempx);
        Py_XDECREF(result);
        return NULL;
    }

    mpfr_clear_flags();

    code = mpfr_sin_cos(s->f, c->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF(tempx);

    s->rc = code & 0x03;
    c->rc = code >> 2;
    if (s->rc == 2) s->rc = -1;
    if (c->rc == 2) c->rc = -1;

    _GMPy_MPFR_Cleanup(&s, context);
    _GMPy_MPFR_Cleanup(&c, context);

    if (!s || !c) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_DECREF(result);
        return NULL;
    }

    PyTuple_SET_ITEM(result, 0, (PyObject*)s);
    PyTuple_SET_ITEM(result, 1, (PyObject*)c);
    return result;
}

static PyObject *
_GMPy_MPC_Sin_Cos(PyObject *x, CTXT_Object *context)
{
    MPC_Object *s = NULL, *c = NULL;
    PyObject *result = NULL;
    int code;

    s = GMPy_MPC_New(0, 0, context);
    c = GMPy_MPC_New(0, 0, context);
    result = PyTuple_New(2);
    if (!s || !c || !result) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        return NULL;
    }

    code = mpc_sin_cos(s->c, c->c, MPC(x), GET_MPC_ROUND(context), GET_MPC_ROUND(context));

    s->rc = MPC_INEX1(code);
    c->rc = MPC_INEX2(code);

    _GMPy_MPC_Cleanup(&s, context);
    _GMPy_MPC_Cleanup(&c, context);

    if (!s || !c) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        return NULL;
    }

    PyTuple_SET_ITEM(result, 0, (PyObject*)s);
    PyTuple_SET_ITEM(result, 1, (PyObject*)c);
    return result;
}

static PyObject *
GMPy_ComplexWithType_Sin_Cos(PyObject *x, int xtype, CTXT_Object *context)
{
    PyObject *result, *tempx;

    if (!(tempx = (PyObject*)GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context))) {
        return NULL;
    }

    result = _GMPy_MPC_Sin_Cos(tempx, context);
    Py_DECREF(tempx);
    return result;
}

PyDoc_STRVAR(GMPy_doc_context_sin_cos,
"sin_cos($self, x, /)\n--\n\n"
"Return a tuple containing the sine and cosine of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_sin_cos,
"sin_cos($module, x, /)\n--\n\n"
"Return a tuple containing the sine and cosine of x; x in radians.");

GMPY_MPFR_MPC_UNIOP_TEMPLATE_EXWT(Sin_Cos, sin_cos)

static PyObject *
GMPy_RealWithType_Sinh_Cosh(PyObject *x, int xtype, CTXT_Object *context)
{
    MPFR_Object *s = NULL, *c = NULL, *tempx = NULL;
    PyObject *result = NULL;
    int code;

    s = GMPy_MPFR_New(0, context);
    c = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context);
    result = PyTuple_New(2);
    if (!s || !c || !tempx || !result) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        return NULL;
    }

    mpfr_clear_flags();
    code = mpfr_sinh_cosh(s->f, c->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF(tempx);

    s->rc = code & 0x03;
    c->rc = code >> 2;
    if (s->rc == 2) s->rc = -1;
    if (c->rc == 2) c->rc = -1;

    _GMPy_MPFR_Cleanup(&s, context);
    _GMPy_MPFR_Cleanup(&c, context);

    if (!s || !c) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_DECREF(result);
        return NULL;
    }

    PyTuple_SET_ITEM(result, 0, (PyObject*)s);
    PyTuple_SET_ITEM(result, 1, (PyObject*)c);
    return result;
}

PyDoc_STRVAR(GMPy_doc_context_sinh_cosh,
"sinh_cosh($self, x, /)\n--\n\n"
"Return a `tuple` containing the hyperbolic sine and cosine of x.");

PyDoc_STRVAR(GMPy_doc_function_sinh_cosh,
"sinh_cosh($module, x, /)\n--\n\n"
"Return a `tuple` containing the hyperbolic sine and cosine of x.");

GMPY_MPFR_UNIOP_TEMPLATEWT(Sinh_Cosh, sinh_cosh)

PyDoc_STRVAR(GMPy_doc_function_degrees,
"degrees($module, x, /)\n--\n\n"
"Convert angle x from radians to degrees.\n"
"Note: In rare cases the result may not be correctly rounded.");

PyDoc_STRVAR(GMPy_doc_context_degrees,
"degrees($self, x, /)\n--\n\n"
"Convert angle x from radians to degrees.\n"
"Note: In rare cases the result may not be correctly rounded.");

static PyObject *
GMPy_Context_Degrees(PyObject *self, PyObject *other)
{
    MPFR_Object *result, *tempx, *temp;
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    result = GMPy_MPFR_New(0, context);
    temp = GMPy_MPFR_New(context->ctx.mpfr_prec + 100, context);
    tempx = GMPy_MPFR_From_Real(other, 1, context);
    if (!result || !temp || !tempx) {
        Py_XDECREF((PyObject*)temp);
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    mpfr_const_pi(temp->f, MPFR_RNDN);
    mpfr_ui_div(temp->f, 180, temp->f, MPFR_RNDN);

    mpfr_clear_flags();

    mpfr_mul(result->f, temp->f, tempx->f, MPFR_RNDN);

    Py_DECREF((PyObject*)temp);
    Py_DECREF((PyObject*)tempx);
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_function_radians,
"radians($module, x, /)\n--\n\n"
"Convert angle x from degrees to radians.\n"
"Note: In rare cases the result may not be correctly rounded.");

PyDoc_STRVAR(GMPy_doc_context_radians,
"radians($self, x, /)\n--\n\n"
"Convert angle x from degrees to radians.\n"
"Note: In rare cases the result may not be correctly rounded.");

static PyObject *
GMPy_Context_Radians(PyObject *self, PyObject *other)
{
    MPFR_Object *result, *tempx, *temp;
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    result = GMPy_MPFR_New(0, context);
    temp = GMPy_MPFR_New(context->ctx.mpfr_prec + 100, context);
    tempx = GMPy_MPFR_From_Real(other, 1, context);
    if (!result || !temp || !tempx) {
        Py_XDECREF((PyObject*)temp);
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    mpfr_const_pi(temp->f, MPFR_RNDN);
    mpfr_div_ui(temp->f, temp->f, 180, MPFR_RNDN);

    mpfr_clear_flags();

    mpfr_mul(result->f, tempx->f, temp->f, MPFR_RNDN);

    Py_DECREF((PyObject*)temp);
    Py_DECREF((PyObject*)tempx);
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_context_log10,
"log10($self, x, /)\n--\n\n"
"Return the base-10 logarithm of x.");

PyDoc_STRVAR(GMPy_doc_function_log10,
"log10($module, x, /)\n--\n\n"
"Return the base-10 logarithm of x.");

GMPY_MPFR_MPC_UNIOP_EXWT(Log10, log10)

PyDoc_STRVAR(GMPy_doc_context_log,
"log($self, x, /)\n--\n\n"
"Return the natural logarithm of x.");

PyDoc_STRVAR(GMPy_doc_function_log,
"log($module, x, /)\n--\n\n"
"Return the natural logarithm of x.");

GMPY_MPFR_MPC_UNIOP_EXWT(Log, log)

PyDoc_STRVAR(GMPy_doc_context_exp,
"exp($self, x, /)\n--\n\n"
"Return the exponential of x.");

PyDoc_STRVAR(GMPy_doc_function_exp,
"exp($module, x, /)\n--\n\n"
"Return the exponential of x.");

GMPY_MPFR_MPC_UNIOP_EXWT(Exp, exp)

PyDoc_STRVAR(GMPy_doc_context_sqrt,
"sqrt($self, x, /)\n--\n\n"
"Return the square root of x.");

PyDoc_STRVAR(GMPy_doc_function_sqrt,
"sqrt($module, x, /)\n--\n\n"
"Return the square root of x.");

static PyObject *
GMPy_RealWithType_Sqrt(PyObject *x, int xtype, CTXT_Object *context)
{
    MPFR_Object *result = NULL;

    CHECK_CONTEXT(context);

    if (IS_TYPE_MPFR(xtype)) {
        if (mpfr_sgn(MPFR(x)) < 0 && context->ctx.allow_complex) {
            return GMPy_ComplexWithType_Sqrt(x, xtype, context);
        }

        if (!(result = GMPy_MPFR_New(0, context))) {
            return NULL;
        }

        mpfr_clear_flags();
        result->rc = mpfr_sqrt(result->f, MPFR(x), GET_MPFR_ROUND(context));
        _GMPy_MPFR_Cleanup(&result, context);
        return (PyObject*)result;
    }

    if (IS_TYPE_REAL(xtype)) {
        MPFR_Object *tempx = NULL;

        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context))) {
            return NULL;
        }

        if (mpfr_sgn(MPFR(tempx)) < 0 && context->ctx.allow_complex) {
            PyObject *res = NULL;

            res = GMPy_ComplexWithType_Sqrt((PyObject*)tempx, OBJ_TYPE_MPFR, context);
            Py_DECREF(tempx);
            return res;
        }
        if (!(result = GMPy_MPFR_New(0, context))) {
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }

        mpfr_clear_flags();
        result->rc = mpfr_sqrt(result->f, MPFR(tempx), GET_MPFR_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        _GMPy_MPFR_Cleanup(&result, context);
        return (PyObject*)result;
    }

    TYPE_ERROR("sqrt() argument type not supported");
    return NULL;
}

static PyObject *
GMPy_ComplexWithType_Sqrt(PyObject *x, int xtype, CTXT_Object *context)
{
    MPC_Object *result = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        return NULL;
    }

    if (IS_TYPE_MPC(xtype)) {
        result->rc = mpc_sqrt(result->c, MPC(x), GET_MPFR_ROUND(context));
        _GMPy_MPC_Cleanup(&result, context);
        return (PyObject*)result;
    }

    if (IS_TYPE_COMPLEX(xtype)) {
        MPC_Object *tempx = NULL;

        if (!(tempx = GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context))) {
            Py_DECREF(result);
            return NULL;
        }

        result->rc = mpc_sqrt(result->c, MPC(tempx), GET_MPFR_ROUND(context));
        Py_DECREF(tempx);
        _GMPy_MPC_Cleanup(&result, context);
        return (PyObject*)result;
    }

    TYPE_ERROR("sqrt() argument type not supported");
    return NULL;
}

GMPY_MPFR_MPC_UNIOP_TEMPLATE_EXWT(Sqrt, sqrt)

PyDoc_STRVAR(GMPy_doc_function_root,
"root($module, x, n, /)\n--\n\n"
"Return n-th root of x. The result always an `mpfr`.\n"
"Note: not IEEE 754-2008 compliant; result differs when\n"
"x = -0 and n is even. See `rootn()`.");

PyDoc_STRVAR(GMPy_doc_context_root,
"root($self, x, n, /)\n--\n\n"
"Return n-th root of x. The result always an `mpfr`.\n"
"Note: not IEEE 754-2008 compliant; result differs when\n"
"x = -0 and n is even. See `context.rootn()`.");

PyDoc_STRVAR(GMPy_doc_function_rootn,
"rootn($module, x, n, /)\n--\n\n"
"Return n-th root of x. The result always an `mpfr`.\n"
"Note: this is IEEE 754-2008 compliant version of `root()`.");

PyDoc_STRVAR(GMPy_doc_context_rootn,
"rootn($self, x, n, /)\n--\n\n"
"Return n-th root of x. The result always an `mpfr`.\n"
"Note: this is IEEE 754-2008 compliant version of `context.root()`.");

#if MPFR_VERSION_MAJOR > 3

/* Since mpfr_root is deprecated in MPFR 4, we use mpfr_rootn_ui to
 * mimic the behavior of mpfr_root. And the converse will be true when
 * using MPFR 3.
 */

static PyObject *
GMPy_Real_Rootn(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *result = NULL, *tempx = NULL;
    unsigned long n;

    CHECK_CONTEXT(context);

    result = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPFR_From_Real(x, 1, context);
    n = GMPy_Integer_AsUnsignedLong(y);

    if (!result || !tempx || (n == (unsigned long)(-1) && PyErr_Occurred())) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    mpfr_clear_flags();
    result->rc = mpfr_rootn_ui(result->f, tempx->f, n, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Real_Root(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *result = NULL, *tempx = NULL;
    unsigned long n;

    CHECK_CONTEXT(context);

    result = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPFR_From_Real(x, 1, context);
    n = GMPy_Integer_AsUnsignedLong(y);

    if (!result || !tempx || (n == (unsigned long)(-1) && PyErr_Occurred())) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    mpfr_clear_flags();

    /* Mimic the non-compliant IEEE 752-2008 behavior. */

    if (mpfr_zero_p(tempx->f)) {
        mpfr_set(result->f, tempx->f, GET_MPFR_ROUND(context));
    }
    else {
        result->rc = mpfr_rootn_ui(result->f, tempx->f, n, GET_MPFR_ROUND(context));
    }

    Py_DECREF((PyObject*)tempx);
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}
#else
static PyObject *
GMPy_Real_Rootn(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *result = NULL, *tempx = NULL;
    unsigned long n;

    CHECK_CONTEXT(context);

    result = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPFR_From_Real(x, 1, context);
    n = GMPy_Integer_AsUnsignedLong(y);

    if (!result || !tempx || (n == (unsigned long)(-1) && PyErr_Occurred())) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    mpfr_clear_flags();

    /* Mimic the compliant IEEE 752-2008 behavior. */

    if (mpfr_zero_p(tempx->f) && mpfr_signbit(tempx->f)) {
        if ((n & 1)) {
            /* Odd, so result is -0. */
            mpfr_set_zero(result->f, -1);
        }
        else {
            /* Even, so result is 0. */
            mpfr_set_zero(result->f, 1);
        }
    }
    else {
        result->rc = mpfr_root(result->f, tempx->f, n, GET_MPFR_ROUND(context));
    }
    Py_DECREF((PyObject*)tempx);
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Real_Root(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *result = NULL, *tempx = NULL;
    unsigned long n;

    CHECK_CONTEXT(context);

    result = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPFR_From_Real(x, 1, context);
    n = GMPy_Integer_AsUnsignedLong(y);

    if (!result || !tempx || (n == (unsigned long)(-1) && PyErr_Occurred())) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    mpfr_clear_flags();
    result->rc = mpfr_root(result->f, tempx->f, n, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}
#endif

static PyObject *
GMPy_Number_Rootn(PyObject *x, PyObject *y, CTXT_Object *context)
{
    if (IS_REAL(x) && PyLong_Check(y))
        return GMPy_Real_Rootn(x, y, context);
    TYPE_ERROR("rootn() argument type not supported");
    return NULL;
}

static PyObject *
GMPy_Context_Rootn(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;
    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("rootn() requires 2 arguments");
        return NULL;
    }
    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }
    return GMPy_Number_Rootn(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context);
}

static PyObject *
GMPy_Number_Root(PyObject *x, PyObject *y, CTXT_Object *context)
{
    if (IS_REAL(x) && PyLong_Check(y))
        return GMPy_Real_Root(x, y, context);
    TYPE_ERROR("root() argument type not supported");
    return NULL;
}

static PyObject *
GMPy_Context_Root(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;
    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("root() requires 2 arguments");
        return NULL;
    }
    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }
    return GMPy_Number_Root(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context);
}

PyDoc_STRVAR(GMPy_doc_function_jn,
"jn($module, n, x, /)\n--\n\n"
"Return the first kind Bessel function of order n of x.\n"
"Note: the order of the arguments changed in gmpy2 2.2.0a2");

PyDoc_STRVAR(GMPy_doc_context_jn,
"jn($self, n, x, /)\n--\n\n"
"Return the first kind Bessel function of order n of x.\n"
"Note: the order of the arguments changed in gmpy2 2.2.0a2");

GMPY_MPFR_BINOP_LONG_REALWT(Jn, jn)

PyDoc_STRVAR(GMPy_doc_function_yn,
"yn($module, n, x, /)\n--\n\n"
"Return the second kind Bessel function of order n of x.\n"
"Note: the order of the arguments changed in gmpy2 2.2.0a2");

PyDoc_STRVAR(GMPy_doc_context_yn,
"yn($self, n, x, /)\n--\n\n"
"Return the second kind Bessel function of order n of x.\n"
"Note: the order of the arguments changed in gmpy2 2.2.0a2");

GMPY_MPFR_BINOP_LONG_REALWT(Yn, yn)

PyDoc_STRVAR(GMPy_doc_function_agm,
"agm($module, x, y, /)\n--\n\n"
"Return arithmetic-geometric mean of x and y.");

PyDoc_STRVAR(GMPy_doc_context_agm,
"agm($self, x, y, /)\n--\n\n"
"Return arithmetic-geometric mean of x and y.");

GMPY_MPFR_BINOPWT(AGM, agm)

PyDoc_STRVAR(GMPy_doc_function_maxnum,
"maxnum($module, x, y, /)\n--\n\n"
"Return the maximum number of x and y. If x and y are not `mpfr`, they are\n"
"converted to `mpfr`. The result is rounded to match the current context.\n"
"If only one of x or y is a number, then that number is returned.");

PyDoc_STRVAR(GMPy_doc_context_maxnum,
"maxnum($self, x, y, /)\n--\n\n"
"Return the maximum number of x and y. If x and y are not `mpfr`, they are\n"
"converted to `mpfr`. The result is rounded to match the specified context.\n"
"If only one of x or y is a number, then that number is returned.");

GMPY_MPFR_BINOPWT(Maxnum, max)

PyDoc_STRVAR(GMPy_doc_function_minnum,
"minnum($module, x, y, /)\n--\n\n"
"Return the minimum number of x and y. If x and y are not `mpfr`, they are\n"
"converted to `mpfr`. The result is rounded to match the current context.\n"
"If only one of x or y is a number, then that number is returned.");

PyDoc_STRVAR(GMPy_doc_context_minnum,
"minnum($self, x, y, /)\n--\n\n"
"Return the minimum number of x and y. If x and y are not `mpfr`, they are\n"
"converted to `mpfr`. The result is rounded to match the specified context.\n"
"If only one of x or y is a number, then that number is returned.");

GMPY_MPFR_BINOPWT(Minnum, min)

PyDoc_STRVAR(GMPy_doc_function_remainder,
"remainder($module, x, y, /)\n--\n\n"
"Return x - n*y where n is the integer quotient of x/y, rounded to\n"
"the nearest integer and ties rounded to even.");

PyDoc_STRVAR(GMPy_doc_context_remainder,
"remainder($self, x, y, /)\n--\n\n"
"Return x - n*y where n is the integer quotient of x/y, rounded to\n"
"the nearest integer and ties rounded to even.");

GMPY_MPFR_BINOPWT(Remainder, remainder)

PyDoc_STRVAR(GMPy_doc_function_fmod,
"fmod($module, x, y, /)\n--\n\n"
"Return x - n*y where n is the integer quotient of x/y, rounded to 0.");

PyDoc_STRVAR(GMPy_doc_context_fmod,
"fmod($self, x, y, /)\n--\n\n"
"Return x - n*y where n is the integer quotient of x/y, rounded to 0.");

GMPY_MPFR_BINOPWT(Fmod, fmod)

PyDoc_STRVAR(GMPy_doc_function_round2,
"round2($module, x, n=0, /)\n--\n\n"
"Return x rounded to n bits.  Uses default precision if n is not specified.\n"
"See `round_away()` to access the mpfr_round() function of the MPFR.");

PyDoc_STRVAR(GMPy_doc_context_round2,
"round2($self, x, n=0, /)\n--\n\n"
"Return x rounded to n bits.  Uses default precision if n is not specified.\n"
"See `context.round_away()` to access the mpfr_round() function of the MPFR.");

static PyObject *
GMPy_Real_Round2(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *result, *tempx;
    long n = 0;

    CHECK_CONTEXT(context);

    if (y) {
        n = PyLong_AsLong(y);
        if (n && ((n == -1 && PyErr_Occurred()) || n < MPFR_PREC_MIN || n > MPFR_PREC_MAX)) {
            VALUE_ERROR("invalid precision");
            return NULL;
        }
    }

    if (!n) {
        n = GET_MPFR_PREC(context);
    }

    if (!(tempx = GMPy_MPFR_From_Real(x, 1, context))) {
        return NULL;
    }
    if (!(result = GMPy_MPFR_New(mpfr_get_prec(tempx->f), context))) {
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    mpfr_set(result->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    mpfr_clear_flags();

    result->rc = mpfr_prec_round(result->f, n, GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Number_Round2(PyObject *x, PyObject *y, CTXT_Object *context)
{
    if (IS_REAL(x) && (!y || PyLong_Check(y)))
        return GMPy_Real_Round2(x, y, context);

    TYPE_ERROR("round2() argument type not supported");
    return NULL;
}

static PyObject *
GMPy_Context_Round2(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) < 1 || PyTuple_GET_SIZE(args) > 2) {
        TYPE_ERROR("round2() requires 1 or 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    if (PyTuple_GET_SIZE(args) == 1) {
        return GMPy_Number_Round2(PyTuple_GET_ITEM(args, 0), NULL, context);
    }
    else {
        return GMPy_Number_Round2(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context);
    }
}

PyDoc_STRVAR(GMPy_doc_function_reldiff,
"reldiff($module, x, y, /)\n--\n\n"
"Return the relative difference between x and y. Result is equal to\n"
"abs(x-y)/x.");

PyDoc_STRVAR(GMPy_doc_context_reldiff,
"reldiff($self, x, y, /)\n--\n\n"
"Return the relative difference between x and y. Result is equal to\n"
"abs(x-y)/x.");

static PyObject *
GMPy_RealWithType_RelDiff(PyObject *x, int xtype, PyObject *y, int ytype, CTXT_Object *context)
{
    MPFR_Object *tempx = NULL, *tempy = NULL, *result = NULL;

    CHECK_CONTEXT(context);

    result = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context);
    tempy = GMPy_MPFR_From_RealWithType(y, ytype, 1, context);
    if (!result || !tempx || !tempy) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)tempy);
        return NULL;
    }

    mpfr_clear_flags();

    mpfr_reldiff(result->f, tempx->f, tempy->f, GET_MPFR_ROUND(context));
    result->rc = 0;
    _GMPy_MPFR_Cleanup(&result, context);
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    return (PyObject*)result;
}

GMPY_MPFR_BINOP_TEMPLATEWT(RelDiff, reldiff)

PyDoc_STRVAR(GMPy_doc_mpfr_ceil_method,
"__ceil__($self)\n--\n\n"
"Return an 'mpfr' that is the smallest integer >= x.");

PyDoc_STRVAR(GMPy_doc_function_ceil,
"ceil($module, x, /)\n--\n\n"
"Return an 'mpfr' that is the smallest integer >= x.");

PyDoc_STRVAR(GMPy_doc_context_ceil,
"ceil($self, x, /)\n--\n\n"
"Return an `mpfr` that is the smallest integer >= x.");

GMPY_MPFR_UNIOP_NOROUNDWT(Ceil, ceil)

PyDoc_STRVAR(GMPy_doc_mpfr_floor_method,
"__floor__($self)\n--\n\n"
"Return an 'mpfr' that is the largest integer <= x.");

PyDoc_STRVAR(GMPy_doc_function_floor,
"floor($module, x, /)\n--\n\n"
"Return an `mpfr` that is the largest integer <= x.");

PyDoc_STRVAR(GMPy_doc_context_floor,
"floor($self, x, /)\n--\n\n"
"Return an `mpfr` that is the largest integer <= x.");

GMPY_MPFR_UNIOP_NOROUNDWT(Floor, floor)

PyDoc_STRVAR(GMPy_doc_mpfr_trunc_method,
"__trunc__($self)\n--\n\n"
"Return an 'mpfr' that is truncated towards 0. Same as\n"
"x.floor() if x>=0 or x.ceil() if x<0.");

PyDoc_STRVAR(GMPy_doc_function_trunc,
"trunc($module, x, /)\n\n"
"Return an `mpfr` that is x truncated towards 0. Same as\n"
"x.floor() if x>=0 or x.ceil() if x<0.");

PyDoc_STRVAR(GMPy_doc_context_trunc,
"trunc($self, x, /)\n--\n\n"
"Return an `mpfr` that is x truncated towards 0. Same as\n"
"x.floor() if x>=0 or x.ceil() if x<0.");

GMPY_MPFR_UNIOP_NOROUNDWT(Trunc, trunc)

PyDoc_STRVAR(GMPy_doc_function_round_away,
"round_away($module, x, /)\n--\n\n"
"Return an `mpfr` that is x rounded to the nearest integer,\n"
"with ties rounded away from 0.");

PyDoc_STRVAR(GMPy_doc_context_round_away,
"round_away($self, x, /)\n--\n\n"
"Return an `mpfr` that is x rounded to the nearest integer,\n"
"with ties rounded away from 0.");

GMPY_MPFR_UNIOP_NOROUND_NOMETHODWT(RoundAway, round)

PyDoc_STRVAR(GMPy_doc_function_modf,
"modf($module, x, /)\n--\n\n"
"Return a `tuple` containing the integer and fractional portions\n"
"of x.");

PyDoc_STRVAR(GMPy_doc_context_modf,
"modf($self, x, /)\n--\n\n"
"Return a `tuple` containing the integer and fractional portions\n"
"of x.");

static PyObject *
GMPy_RealWithType_Modf(PyObject *x, int xtype, CTXT_Object *context)
{
    MPFR_Object *s = NULL, *c = NULL, *tempx = NULL;
    PyObject *result = NULL;
    int code;

    tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context);
    s = GMPy_MPFR_New(0, context);
    c = GMPy_MPFR_New(0, context);
    result = PyTuple_New(2);
    if (! tempx || !s || !c || !result) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        return NULL;
    }

    mpfr_clear_flags();

    code = mpfr_modf(s->f, c->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    s->rc = code & 0x03;
    c->rc = code >> 2;
    if (s->rc == 2) s->rc = -1;
    if (c->rc == 2) c->rc = -1;

    _GMPy_MPFR_Cleanup(&s, context);
    _GMPy_MPFR_Cleanup(&c, context);

    if (!s || !c) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_DECREF(result);
        return NULL;
    }

    PyTuple_SET_ITEM(result, 0, (PyObject*)s);
    PyTuple_SET_ITEM(result, 1, (PyObject*)c);
    return result;
}

GMPY_MPFR_UNIOP_TEMPLATEWT(Modf, modf)

PyDoc_STRVAR(GMPy_doc_function_lgamma,
"lgamma($module, x, /)\n--\n\n"
"Return a `tuple` containing the logarithm of the absolute value of\n"
"gamma(x) and the sign of gamma(x)");

PyDoc_STRVAR(GMPy_doc_context_lgamma,
"lgamma($self, x, /)\n--\n\n"
"Return a `tuple` containing the logarithm of the absolute value of\n"
"gamma(x) and the sign of gamma(x)");

static PyObject *
GMPy_RealWithType_Lgamma(PyObject *x, int xtype, CTXT_Object *context)
{
    PyObject *result = NULL;
    MPFR_Object *value = NULL, *tempx = NULL;
    int signp = 0;

    tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context);
    value = GMPy_MPFR_New(0, context);
    result = PyTuple_New(2);
    if (!tempx || !value || !result) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)value);
        Py_XDECREF(result);
        return NULL;
    }

    mpfr_clear_flags();

    value->rc = mpfr_lgamma(value->f, &signp, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    _GMPy_MPFR_Cleanup(&value, context);

    if (!value) {
        Py_DECREF(result);
        return NULL;
    }

    PyTuple_SET_ITEM(result, 0, (PyObject*)value);
    PyTuple_SET_ITEM(result, 1, PyLong_FromLong((long)signp));
    return result;
}

GMPY_MPFR_UNIOP_TEMPLATEWT(Lgamma, lgamma)

PyDoc_STRVAR(GMPy_doc_function_remquo,
"remquo($module, x, y, /)\n--\n\n"
"Return a `tuple` containing the remainder(x,y) and the low bits of the\n"
"quotient.");

PyDoc_STRVAR(GMPy_doc_context_remquo,
"remquo($self, x, y, /)\n--\n\n"
"Return a `tuple` containing the remainder(x,y) and the low bits of the\n"
"quotient.");

static PyObject *
GMPy_RealWithType_RemQuo(PyObject *x, int xtype, PyObject *y, int ytype, CTXT_Object *context)
{
    PyObject *result;
    MPFR_Object *value, *tempx, *tempy;
    long quobits = 0;

    CHECK_CONTEXT(context);

    value = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context);
    tempy = GMPy_MPFR_From_RealWithType(y, ytype, 1, context);
    result = PyTuple_New(2);
    if (!value || !tempx || !tempx || !result) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)tempy);
        Py_XDECREF((PyObject*)value);
        Py_XDECREF(result);
        return NULL;
    }

    mpfr_clear_flags();

    value->rc = mpfr_remquo(value->f, &quobits, tempx->f, tempy->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    _GMPy_MPFR_Cleanup(&value, context);

    PyTuple_SET_ITEM(result, 0, (PyObject*)value);
    PyTuple_SET_ITEM(result, 1, PyLong_FromLong(quobits));
    return result;
}

GMPY_MPFR_BINOP_TEMPLATEWT(RemQuo, remquo)

PyDoc_STRVAR(GMPy_doc_function_frexp,
"frexp($module, x, /)\n--\n\n"
"Return a `tuple` containing the exponent and mantissa of x.");

PyDoc_STRVAR(GMPy_doc_context_frexp,
"frexp($self, x, /)\n--\n\n"
"Return a `tuple` containing the exponent and mantissa of x.");

static PyObject *
GMPy_RealWithType_Frexp(PyObject *x, int xtype, CTXT_Object *context)
{
    PyObject *result = NULL;
    MPFR_Object *value = NULL, *tempx = NULL;
    mpfr_exp_t exp = 0;

    value = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context);
    result = PyTuple_New(2);
    if (!value || !result || !tempx) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)value);
        Py_XDECREF(result);
        return NULL;
    }

    mpfr_clear_flags();
    value->rc = mpfr_frexp(&exp, value->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    _GMPy_MPFR_Cleanup(&value, context);

    PyTuple_SET_ITEM(result, 0, PyLong_FromSsize_t((Py_ssize_t)exp));
    PyTuple_SET_ITEM(result, 1, (PyObject*)value);
    return result;
}

GMPY_MPFR_UNIOP_TEMPLATEWT(Frexp, frexp)

PyDoc_STRVAR(GMPy_doc_function_next_toward,
"next_toward($module, x, y, /)\n--\n\n"
"Return the next `mpfr` from x in the direction of y. The result has\n"
"the same precision as x.");

PyDoc_STRVAR(GMPy_doc_context_next_toward,
"next_toward($self, x, y, /)\n--\n\n"
"Return the next `mpfr` from x in the direction of y. The result has\n"
"the same precision as x.");

static PyObject *
GMPy_Context_NextToward(PyObject *self, PyObject *args)
{
    MPFR_Object *result, *tempx, *tempy;
    CTXT_Object *context = NULL;
    int direction;
    mpfr_rnd_t temp_round;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("next_toward() requires 2 arguments");
        return NULL;
    }

    tempx = GMPy_MPFR_From_Real(PyTuple_GET_ITEM(args, 0), 1, context);
    tempy = GMPy_MPFR_From_Real(PyTuple_GET_ITEM(args, 1), 1, context);
    if (!tempx || !tempy) {
        TYPE_ERROR("next_toward() argument type not supported");
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)tempy);
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(mpfr_get_prec(tempx->f), context))) {
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return NULL;
    }

    mpfr_clear_flags();

    mpfr_set(result->f, tempx->f, GET_MPFR_ROUND(context));
    mpfr_nexttoward(result->f, tempy->f);
    result->rc = 0;
    direction = mpfr_signbit(tempy->f);
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    temp_round = GET_MPFR_ROUND(context);
    if (direction)
        context->ctx.mpfr_round = MPFR_RNDD;
    else
         context->ctx.mpfr_round = MPFR_RNDU;
    _GMPy_MPFR_Cleanup(&result, context);
    context->ctx.mpfr_round = temp_round;
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_function_next_above,
"next_above($module, x, /)\n--\n\n"
"Return the next `mpfr` from x toward +Infinity.");

PyDoc_STRVAR(GMPy_doc_context_next_above,
"next_above($self, x, /)\n--\n\n"
"Return the next `mpfr` from x toward +Infinity.");

static PyObject *
GMPy_Context_NextAbove(PyObject *self, PyObject *other)
{
    MPFR_Object *result, *tempx;
    CTXT_Object *context = NULL;
    mpfr_rnd_t temp_round;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    if (!(tempx = GMPy_MPFR_From_Real(other, 1, context))) {
        TYPE_ERROR("next_above() argument type not supported");
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(mpfr_get_prec(tempx->f), context))) {
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    mpfr_clear_flags();

    mpfr_set(result->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    mpfr_nextabove(result->f);
    result->rc = 0;
    temp_round = GET_MPFR_ROUND(context);
    context->ctx.mpfr_round = MPFR_RNDU;
    _GMPy_MPFR_Cleanup(&result, context);
    context->ctx.mpfr_round = temp_round;
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_function_next_below,
"next_below($module, x, /)\n--\n\n"
"Return the next `mpfr` from x toward -Infinity.");

PyDoc_STRVAR(GMPy_doc_context_next_below,
"next_below($self, x, /)\n--\n\n"
"Return the next `mpfr` from x toward -Infinity.");

static PyObject *
GMPy_Context_NextBelow(PyObject *self, PyObject *other)
{
    MPFR_Object *result, *tempx;
    CTXT_Object *context = NULL;
    mpfr_rnd_t temp_round;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    if (!(tempx = GMPy_MPFR_From_Real(other, 1, context))) {
        TYPE_ERROR("next_below() argument type not supported");
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(mpfr_get_prec(tempx->f), context))) {
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    mpfr_clear_flags();

    mpfr_set(result->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    mpfr_nextbelow(result->f);
    result->rc = 0;
    temp_round = GET_MPFR_ROUND(context);
    context->ctx.mpfr_round = MPFR_RNDD;
    _GMPy_MPFR_Cleanup(&result, context);
    context->ctx.mpfr_round = temp_round;
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_function_factorial,
"factorial($module, n, /)\n--\n\n"
"Return the floating-point approximation to the factorial of n.\n\n"
"See `fac()` to get the exact integer result.");

PyDoc_STRVAR(GMPy_doc_context_factorial,
"factorial($self, n, /)\n--\n\n"
"Return the floating-point approximation to the factorial of n.\n\n"
"See `fac()` to get the exact integer result.");

static PyObject *
GMPy_Context_Factorial(PyObject *self, PyObject *other)
{
    MPFR_Object *result;
    unsigned long n;
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    n = GMPy_Integer_AsUnsignedLong(other);
    if ((n == (unsigned long)-1) && PyErr_Occurred()) {
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    /* Force result to be 'inf' if n >= 44787928. The matches the behavior of
     * MPFR compiled on 64-bit Linux. MPFR compiled on 32-bit Linux and both 32
     * and 64-bit versions of Windows will enter an infinite loop. The constant
     * 44787929 occurs in the MPFR file gamma.c.
     */

    mpfr_clear_flags();

    if (n >= 44787928) {
        mpfr_set_inf(result->f, 1);
        mpfr_set_overflow();
    }
    else {
        mpfr_fac_ui(result->f, n, GET_MPFR_ROUND(context));
    }

    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_function_fsum,
"fsum($module, iterable, /)\n--\n\n"
"Return an accurate sum of the values in the iterable.");

PyDoc_STRVAR(GMPy_doc_context_fsum,
"fsum($self, iterable, /)\n--\n\n"
"Return an accurate sum of the values in the iterable.");

static PyObject *
GMPy_Context_Fsum(PyObject *self, PyObject *other)
{
    MPFR_Object *temp, *result;
    mpfr_ptr *tab;
    int errcode;
    Py_ssize_t i, seq_length = 0;
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    if (!(other = PySequence_List(other))) {
        Py_DECREF((PyObject*)result);
        TYPE_ERROR("argument must be an iterable");
        return NULL;
    }

    /* other contains a new list containing all the values from the
     * iterable. Now make sure each item in the list is an mpfr.
     */

    seq_length = PyList_GET_SIZE(other);
    if (seq_length > LONG_MAX) {
        OVERFLOW_ERROR("temporary array is too large");
	Py_DECREF(other);
	Py_DECREF((PyObject*)result);
	return NULL;
    }
    for (i=0; i < seq_length; i++) {
        if (!(temp = GMPy_MPFR_From_Real(PyList_GET_ITEM(other, i), 1, context))) {
            Py_DECREF(other);
            Py_DECREF((PyObject*)result);
            TYPE_ERROR("all items in iterable must be real numbers");
            return NULL;
        }

        errcode = PyList_SetItem(other, i,(PyObject*)temp);
        if (errcode < 0) {
            Py_DECREF(other);
            Py_DECREF((PyObject*)result);
            TYPE_ERROR("all items in iterable must be real numbers");
            return NULL;
        }
    }

    /* create an array of pointers to the mpfr_t field of a Pympfr object */

    if (!(tab = (mpfr_ptr *)malloc((sizeof(mpfr_srcptr) * seq_length)))) {
        Py_DECREF(other);
        Py_DECREF((PyObject*)result);
        return PyErr_NoMemory();
    }
    for (i=0; i < seq_length; i++) {
        temp = (MPFR_Object*)PyList_GET_ITEM(other, i);
        tab[i] = temp->f;
    }

    mpfr_clear_flags();

    /* The cast is safe since we have compared seq_length to LONG_MAX. */
    result->rc = mpfr_sum(result->f, tab, (unsigned long)seq_length, GET_MPFR_ROUND(context));
    Py_DECREF(other);
    free(tab);

    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}
