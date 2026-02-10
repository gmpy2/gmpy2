/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_macros.h                                                          *
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

/* This file contains a collection of macros that can be used to reduce
 * repetitive code. As new macros are written to support of refactoring,
 * they should be placed here.
 */

/* NAME is used as part of the GMPy function name. It usually uses an upper-
 *      case first character.
 * FUNC is the component of the actual name used by MPFR and MPC.
 *
 * Note: the following macro can release the GIL.
 * GMPY_MPFR_MPC_UNIOP_EXWT(NAME, FUNC) creates the following functions:
 *     GMPy_RealWithType_NAME(x, xtype, context)
 *     GMPy_ComplexWithType_NAME(x, xtype, context)
 *     GMPy_Number_NAME(x, context)
 *     GMPy_Context_NAME(self, other)
 *     - called with METH_O
 *
 * Note: the following macro is only used for is_xxx tests so it does 
 *       not release the GIL.
 * GMPY_MPFR_MPC_UNIOP_TEMPLATEWT(NAME, FUNC) creates the following functions:
 *     GMPy_Number_NAME(x, context)
 *     - assumes GMPy_RealWithType_NAME & GMPy_ComplexWithType_NAME exist
 *     GMPy_Context_NAME(self, other)
 *     - called with METH_O
 *
 * GMPY_MPFR_MPC_TRIOP_EX(NAME, FUNC) creates the following functions:
 *     GMPy_Real_NAME(x, y, Z, context)
 *     GMPy_Complex_NAME(x, y, Z, context)
 *     GMPy_Number_NAME(x, y, Z, context)
 *     - assumes GMPy_Integer_NAME & GMPy_Rational_NAME also exist
 *     GMPy_Context_NAME(self, args)
 *     - called with METH_VARARGS
 *
 * GMPY_MPFR_MPC_UNIOP_TEMPLATE(NAME, FUNC) creates the following functions:
 *     GMPy_Number_NAME(x, context)
 *     - assumes GMPy_Real_NAME & GMPy_Complex_NAME exist
 *     GMPy_Context_NAME(self, other)
 *     - called with METH_O
 *
 * GMPY_MPFR_UNIOP(NAME, FUNC) creates the following functions:
 *     GMPy_Real_NAME(x, context)
 *     GMPy_Number_NAME(x, context)
 *     GMPy_Context_NAME(self, other)
 *     - called with METH_O
 *
 * GMPY_MPFR_UNIOP_TEMPLATE(NAME, FUNC) creates the following functions:
 *     GMPy_Number_NAME(x, context)
 *     - assumes GMPy_Real_NAME exists
 *     GMPy_Context_NAME(self, other)
 *     - called with METH_O
 *
 * GMPY_MPFR_BINOP(NAME, FUNC) creates the following functions:
 *     GMPy_Real_NAME(x, y, context)
 *     GMPy_Number_NAME(x, y, context)
 *     GMPy_Context_NAME(self, args)
 *     - called with METH_VARARGS
 *
 */

#define GMPY_MPFR_MPC_UNIOP_EXWT(NAME, FUNC) \
static PyObject * \
GMPy_RealWithType_##NAME(PyObject *x, int xtype, CTXT_Object *context) \
{ \
    MPFR_Object *result = NULL, *tempx = NULL; \
    if (IS_TYPE_MPFR(xtype)) { \
        if (!(result = GMPy_MPFR_New(0, context))) return NULL; \
        mpfr_clear_flags(); \
        result->rc = mpfr_##FUNC(result->f, MPFR(x), GET_MPFR_ROUND(context)); \
        _GMPy_MPFR_Cleanup(&result, context); \
        return (PyObject*)result; \
    } \
    if (IS_TYPE_REAL(xtype)) { \
        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context))) return NULL; \
        if (!(result = GMPy_MPFR_New(0, context))) { \
            Py_DECREF(tempx); \
            return NULL; \
        } \
        mpfr_clear_flags(); \
        result->rc = mpfr_##FUNC(result->f, MPFR(tempx), GET_MPFR_ROUND(context)); \
        _GMPy_MPFR_Cleanup(&result, context); \
        Py_DECREF(tempx); \
        return (PyObject*)result; \
    } \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_ComplexWithType_##NAME(PyObject *x, int xtype, CTXT_Object *context) \
{ \
    MPC_Object *result = NULL, *tempx = NULL; \
    if (IS_TYPE_MPC(xtype)) { \
        if (!(result = GMPy_MPC_New(0, 0, context))) return NULL; \
        result->rc = mpc_##FUNC(result->c, MPC(x), GET_MPC_ROUND(context)); \
        _GMPy_MPC_Cleanup(&result, context); \
        return (PyObject*)result; \
    } \
    if (IS_TYPE_COMPLEX(xtype)) { \
        if (!(tempx = GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context))) return NULL; \
        if (!(result = GMPy_MPC_New(0, 0, context))) { \
            Py_DECREF(tempx); \
            return NULL; \
        } \
        result->rc = mpc_##FUNC(result->c, MPC(tempx), GET_MPC_ROUND(context)); \
        _GMPy_MPC_Cleanup(&result, context); \
        Py_DECREF(tempx); \
        return (PyObject*)result; \
    } \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    CHECK_CONTEXT(context); \
    int xtype = GMPy_ObjectType(x);\
    if (IS_TYPE_REAL(xtype)) \
        return GMPy_RealWithType_##NAME(x, xtype, context); \
    if (IS_TYPE_COMPLEX(xtype)) \
        return GMPy_ComplexWithType_##NAME(x, xtype, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(other, context); \
}

/*********************************************************************/

#define GMPY_MPFR_MPC_UNIOP_TEMPLATEWT(NAME, FUNC) \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    CHECK_CONTEXT(context); \
    int xtype = GMPy_ObjectType(x); \
    if (IS_TYPE_REAL(xtype)) \
        return GMPy_RealWithType_##NAME(x, xtype, context); \
    if (IS_TYPE_COMPLEX(xtype)) \
        return GMPy_ComplexWithType_##NAME(x, xtype, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(other, context); \
} \
static PyObject * \
GMPy_Number_Method_##NAME(PyObject *self, PyObject *args) \
{ \
    return GMPy_Number_##NAME(self, NULL); \
} \

/*********************************************************************/

#define GMPY_MPFR_MPC_UNIOP_TEMPLATE_EXWT(NAME, FUNC) \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    CHECK_CONTEXT(context); \
    int xtype = GMPy_ObjectType(x); \
    if (IS_TYPE_REAL(xtype)) \
        return GMPy_RealWithType_##NAME(x, xtype, context); \
    if (IS_TYPE_COMPLEX(xtype)) \
        return GMPy_ComplexWithType_##NAME(x, xtype, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(other, context); \
} \

/*********************************************************************/

#define GMPY_MPFR_MPC_TRIOP_TEMPLATEWT(NAME, FUNC) \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context) \
{ \
    CHECK_CONTEXT(context); \
    int xtype = GMPy_ObjectType(x); \
    int ytype = GMPy_ObjectType(y); \
    int ztype = GMPy_ObjectType(z); \
    if (IS_TYPE_MPZ(xtype) && IS_TYPE_MPZ(ytype) && IS_TYPE_MPZ(ztype)) \
        return _GMPy_MPZ_##NAME(x, y, z, context); \
    if (IS_TYPE_MPQ(xtype) && IS_TYPE_MPQ(ytype) && IS_TYPE_MPQ(ztype)) \
        return _GMPy_MPQ_##NAME(x, y,  z, context); \
    if (IS_TYPE_MPFR(xtype) && IS_TYPE_MPFR(ytype) && IS_TYPE_MPFR(ztype)) \
        return _GMPy_MPFR_##NAME(x, y, z, context); \
    if (IS_TYPE_MPC(xtype) && IS_TYPE_MPC(ytype) && IS_TYPE_MPC(ztype)) \
        return _GMPy_MPC_##NAME(x, y, z, context); \
    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype) && IS_TYPE_INTEGER(ztype)) \
        return GMPy_IntegerWithType_##NAME(x, xtype, y, ytype, z, ztype, context); \
    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype) && IS_TYPE_RATIONAL(ztype)) \
        return GMPy_RationalWithType_##NAME(x, xtype, y, ytype, z, ztype, context); \
    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype) && IS_TYPE_REAL(ztype)) \
        return GMPy_RealWithType_##NAME(x, xtype, y, ytype, z, ztype, context); \
    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype) && IS_TYPE_COMPLEX(ztype)) \
        return GMPy_ComplexWithType_##NAME(x, xtype, y, ytype, z, ztype, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *args) \
{ \
    CTXT_Object *context = NULL; \
    if (PyTuple_GET_SIZE(args) != 3) { \
        TYPE_ERROR(#FUNC"() requires 3 arguments"); \
        return NULL; \
    } \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(PyTuple_GET_ITEM(args, 0), \
                              PyTuple_GET_ITEM(args, 1), \
                              PyTuple_GET_ITEM(args, 2), context); \
}

/*********************************************************************/

#define GMPY_MPFR_QUADOP_TEMPLATEWT(NAME, FUNC) \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context) \
{ \
    CHECK_CONTEXT(context); \
    int xtype = GMPy_ObjectType(x); \
    int ytype = GMPy_ObjectType(y); \
    int ztype = GMPy_ObjectType(z); \
    int ttype = GMPy_ObjectType(t); \
    if (IS_TYPE_MPZ(xtype) && IS_TYPE_MPZ(ytype) && IS_TYPE_MPZ(ztype) && IS_TYPE_MPZ(ttype)) \
        return _GMPy_MPZ_##NAME(x, y, z, t, context); \
    if (IS_TYPE_MPQ(xtype) && IS_TYPE_MPQ(ytype) && IS_TYPE_MPQ(ztype) && IS_TYPE_MPQ(ttype)) \
        return _GMPy_MPQ_##NAME(x, y,  z, t, context); \
    if (IS_TYPE_MPFR(xtype) && IS_TYPE_MPFR(ytype) && IS_TYPE_MPFR(ztype) && IS_TYPE_MPFR(ttype)) \
        return _GMPy_MPFR_##NAME(x, y, z, t, context); \
    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype) && IS_TYPE_INTEGER(ztype) && IS_TYPE_INTEGER(ttype)) \
        return GMPy_IntegerWithType_##NAME(x, xtype, y, ytype, z, ztype, t, ttype, context); \
    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype) && IS_TYPE_RATIONAL(ztype) && IS_TYPE_RATIONAL(ttype)) \
        return GMPy_RationalWithType_##NAME(x, xtype, y, ytype, z, ztype, t, ttype, context); \
    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype) && IS_TYPE_REAL(ztype) && IS_TYPE_REAL(ttype)) \
        return GMPy_RealWithType_##NAME(x, xtype, y, ytype, z, ztype, t, ttype, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *args) \
{ \
    CTXT_Object *context = NULL; \
    if (PyTuple_GET_SIZE(args) != 4) { \
        TYPE_ERROR(#FUNC"() requires 4 arguments"); \
        return NULL; \
    } \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(PyTuple_GET_ITEM(args, 0), \
                              PyTuple_GET_ITEM(args, 1), \
                              PyTuple_GET_ITEM(args, 2), \
                              PyTuple_GET_ITEM(args, 3), context); \
}

/*********************************************************************/

#define GMPY_MPFR_UNIOP_NOROUNDWT(NAME, FUNC) \
static PyObject * \
GMPy_RealWithType_##NAME(PyObject *x, int xtype, CTXT_Object *context) \
{ \
    MPFR_Object *result = NULL, *tempx = NULL; \
    result = GMPy_MPFR_New(0, context); \
    tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context); \
    if (!result || !tempx) { \
        Py_XDECREF((PyObject*)result); \
        Py_XDECREF((PyObject*)tempx); \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, tempx->f); \
    Py_DECREF((PyObject*)tempx); \
    _GMPy_MPFR_Cleanup(&result, context); \
    MPZ_Object *mpz_result = GMPy_MPZ_From_MPFR(result, context); \
    Py_DECREF((PyObject*)result); \
    return (PyObject*)mpz_result; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    int xtype = GMPy_ObjectType(x); \
    if (IS_TYPE_REAL(xtype)) \
        return GMPy_RealWithType_##NAME(x, xtype, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_MPFR_Method_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    CHECK_CONTEXT(context); \
    return GMPy_Number_##NAME(self, context); \
}\
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(other, context); \
}

/*********************************************************************/

#define GMPY_MPFR_UNIOP_NOROUND_NOMETHODWT(NAME, FUNC) \
static PyObject * \
GMPy_RealWithType_##NAME(PyObject *x, int xtype, CTXT_Object *context) \
{ \
    MPFR_Object *result = NULL, *tempx = NULL; \
    result = GMPy_MPFR_New(0, context); \
    tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context); \
    if (!result || !tempx) { \
        Py_XDECREF((PyObject*)result); \
        Py_XDECREF((PyObject*)tempx); \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, tempx->f); \
    Py_DECREF((PyObject*)tempx); \
    _GMPy_MPFR_Cleanup(&result, context); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    int xtype = GMPy_ObjectType(x); \
    if (IS_TYPE_REAL(xtype)) \
        return GMPy_RealWithType_##NAME(x, xtype, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(other, context); \
}

/*********************************************************************/

#define GMPY_MPFR_UNIOP_EXWT(NAME, FUNC) \
static PyObject * \
GMPy_RealWithType_##NAME(PyObject *x, int xtype, CTXT_Object *context) \
{ \
    MPFR_Object *result = NULL, *tempx = NULL; \
    if (IS_TYPE_MPFR(xtype)) { \
        if (!(result = GMPy_MPFR_New(0, context))) return NULL; \
        mpfr_clear_flags(); \
        result->rc = mpfr_##FUNC(result->f, MPFR(x), GET_MPFR_ROUND(context)); \
        _GMPy_MPFR_Cleanup(&result, context); \
        return (PyObject*)result; \
    } \
    if (IS_TYPE_REAL(xtype)) { \
        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context))) return NULL; \
        if (!(result = GMPy_MPFR_New(0, context))) { \
            Py_DECREF(tempx); \
            return NULL; \
        } \
        mpfr_clear_flags(); \
        result->rc = mpfr_##FUNC(result->f, MPFR(tempx), GET_MPFR_ROUND(context)); \
        Py_DECREF(tempx); \
        _GMPy_MPFR_Cleanup(&result, context); \
        return (PyObject*)result; \
    } \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    CHECK_CONTEXT(context); \
    int xtype = GMPy_ObjectType(x);\
    if (IS_TYPE_REAL(xtype)) \
        return GMPy_RealWithType_##NAME(x, xtype, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(other, context); \
}

/*********************************************************************/

#define GMPY_MPFR_UNIOP_TEMPLATEWT(NAME, FUNC) \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    CHECK_CONTEXT(context); \
    int xtype = GMPy_ObjectType(x);\
    if (IS_TYPE_REAL(xtype)) \
        return GMPy_RealWithType_##NAME(x, xtype, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(other, context); \
}

/*********************************************************************/

#define GMPY_MPFR_BINOPWT(NAME, FUNC) \
static PyObject * \
GMPy_RealWithType_##NAME(PyObject *x, int xtype, PyObject *y, int ytype, CTXT_Object *context) \
{ \
    MPFR_Object *result = NULL, *tempx = NULL, *tempy = NULL; \
    tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context); \
    tempy = GMPy_MPFR_From_RealWithType(y, ytype, 1, context); \
    result = GMPy_MPFR_New(0, context); \
    if (!result || !tempx || !tempy) { \
        Py_XDECREF((PyObject*)tempx); \
        Py_XDECREF((PyObject*)tempy); \
        Py_XDECREF((PyObject*)result); \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, tempx->f, tempy->f, GET_MPFR_ROUND(context)); \
    Py_DECREF((PyObject*)tempx); \
    Py_DECREF((PyObject*)tempy); \
    _GMPy_MPFR_Cleanup(&result, context); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, PyObject *y, CTXT_Object *context) \
{ \
    int xtype = GMPy_ObjectType(x); \
    int ytype = GMPy_ObjectType(y); \
    CHECK_CONTEXT(context); \
    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype)) \
        return GMPy_RealWithType_##NAME(x, xtype, y, ytype, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *args) \
{ \
    CTXT_Object *context = NULL; \
    if (PyTuple_GET_SIZE(args) != 2) { \
        TYPE_ERROR(#FUNC"() requires 2 arguments"); \
        return NULL; \
    } \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context); \
} \

/* Macro to support functions that require ('mpfr', 'int').
 * More precisely, the first argument must pass IS_REAL() and the second
 * argument must pass IS_INTEGER(). The calling sequence passes n first
 * to the MPFR library.*/

/*********************************************************************/

#define GMPY_MPFR_BINOP_REAL_LONGWT(NAME, FUNC) \
static PyObject * \
GMPy_RealWithType_##NAME(PyObject *x, int xtype, PyObject *y, int ytype, CTXT_Object *context) \
{ \
    MPFR_Object *result = NULL, *tempx = NULL; \
    long n; \
    result = GMPy_MPFR_New(0, context); \
    tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context); \
    n = GMPy_Integer_AsLongWithType(y, ytype); \
    if (!result || !tempx || (n == -1 && PyErr_Occurred())) { \
        Py_XDECREF((PyObject*)tempx); \
        Py_XDECREF((PyObject*)result); \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, n, tempx->f, GET_MPFR_ROUND(context)); \
    Py_DECREF((PyObject*)tempx); \
    _GMPy_MPFR_Cleanup(&result, context); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, PyObject *y, CTXT_Object *context) \
{ \
    int xtype = GMPy_ObjectType(x); \
    int ytype = GMPy_ObjectType(y); \
    CHECK_CONTEXT(context); \
    if (IS_TYPE_REAL(xtype) && IS_TYPE_INTEGER(ytype)) \
        return GMPy_RealWithType_##NAME(x, xtype, y, ytype, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *args) \
{ \
    CTXT_Object *context = NULL; \
    if (PyTuple_GET_SIZE(args) != 2) { \
        TYPE_ERROR(#FUNC"() requires 2 arguments"); \
        return NULL; \
    } \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context); \
} \

/*********************************************************************/

/* Support standard calling convention for jn() and yn().*/

#define GMPY_MPFR_BINOP_LONG_REALWT(NAME, FUNC) \
static PyObject * \
GMPy_RealWithType_##NAME(PyObject *x, int xtype, PyObject *y, int ytype, CTXT_Object *context) \
{ \
    MPFR_Object *result = NULL, *temp = NULL; \
    long n; \
    result = GMPy_MPFR_New(0, context); \
    temp = GMPy_MPFR_From_RealWithType(y, ytype, 1, context); \
    n = GMPy_Integer_AsLongWithType(x, xtype); \
    if (!result || !temp || (n == -1 && PyErr_Occurred())) { \
        Py_XDECREF((PyObject*)temp); \
        Py_XDECREF((PyObject*)result); \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, n, temp->f, GET_MPFR_ROUND(context)); \
    Py_DECREF((PyObject*)temp); \
    _GMPy_MPFR_Cleanup(&result, context); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, PyObject *y, CTXT_Object *context) \
{ \
    int xtype = GMPy_ObjectType(x); \
    int ytype = GMPy_ObjectType(y); \
    CHECK_CONTEXT(context); \
    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_REAL(ytype)) \
        return GMPy_RealWithType_##NAME(x, xtype, y, ytype, context); \
    TYPE_ERROR(#FUNC"() argument types not supported. Note that the argument order " \
                    "for jn() and yn() has changed to (int, mpfr) instead of (mpfr,int)."); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *args) \
{ \
    CTXT_Object *context = NULL; \
    if (PyTuple_GET_SIZE(args) != 2) { \
        TYPE_ERROR(#FUNC"() requires 2 arguments"); \
        return NULL; \
    } \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context); \
} \

/*********************************************************************/

#define GMPY_MPFR_BINOP_TEMPLATEWT(NAME, FUNC) \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, PyObject *y, CTXT_Object *context) \
{ \
    int xtype = GMPy_ObjectType(x); \
    int ytype = GMPy_ObjectType(y); \
    CHECK_CONTEXT(context); \
    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype)) \
        return GMPy_RealWithType_##NAME(x, xtype, y, ytype, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *args) \
{ \
    CTXT_Object *context = NULL; \
    if (PyTuple_GET_SIZE(args) != 2) { \
        TYPE_ERROR(#FUNC"() requires 2 arguments"); \
        return NULL; \
    } \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context); \
} \

/*********************************************************************/

#define GMPY_MPFR_BINOP_EXWT(NAME, FUNC) \
static PyObject * \
GMPy_RealWithType_##NAME(PyObject *x, int xtype, PyObject *y, int ytype, CTXT_Object *context) \
{ \
    MPFR_Object *result = NULL, *tempx = NULL, *tempy = NULL; \
    CHECK_CONTEXT(context); \
    tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context); \
    tempy = GMPy_MPFR_From_RealWithType(y, ytype, 1, context); \
    result = GMPy_MPFR_New(0, context); \
    if (!tempx || !tempy || !result) { \
        Py_XDECREF(tempx); \
        Py_XDECREF(tempy); \
        Py_XDECREF(result); \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, tempx->f, tempy->f, GET_MPFR_ROUND(context)); \
    Py_DECREF(tempx); \
    Py_DECREF(tempy); \
    _GMPy_MPFR_Cleanup(&result, context); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, PyObject *y, CTXT_Object *context) \
{ \
    int xtype = GMPy_ObjectType(x); \
    int ytype = GMPy_ObjectType(y); \
    CHECK_CONTEXT(context); \
    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype)) \
        return GMPy_RealWithType_##NAME(x, xtype, y, ytype, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *args) \
{ \
    CTXT_Object *context = NULL; \
    if (PyTuple_GET_SIZE(args) != 2) { \
        TYPE_ERROR(#FUNC"() requires 2 arguments"); \
        return NULL; \
    } \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context); \
}

#define SWAP(T, a, b)  \
    do {               \
        T tmp = a;     \
        a = b;         \
        b = tmp;       \
    } while (0);

static inline void
revstr(char *s, size_t l, size_t r)
{
    while (l < r) {
        SWAP(char, s[l], s[r]);
        l++;
        r--;
    }
}
