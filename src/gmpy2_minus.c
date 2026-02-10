/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_minus.c                                                           *
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

/* This file implements __neg__ and context.minus().
 *
 * Public API
 * ==========
 * The following function is available as part of GMPY2's C API. If the value
 * of context is NULL, then the function should use the currently active
 * context.
 *
 *   GMPy_Number_Minus(Number, context)
 *
 * Private API
 * ===========
 *   GMPy_MPZ_Minus_Slot
 *   GMPy_MPQ_Minus_Slot
 *   GMPy_MPFR_Minus_Slot
 *   GMPy_MPC_Minus_Slot
 *
 *   GMPy_Integer_Minus(Integer, context|NULL)
 *   GMPy_Rational_Minus(Rational, context|NULL)
 *   GMPy_Real_Minus(Real, context|NULL)
 *   GMPy_Complex_Minus(Complex, context|NULL)
 *
 *   GMPy_Context_Minus(context, args)
 */

static PyObject *
_GMPy_MPZ_Minus(PyObject *x, CTXT_Object *context)
{
    MPZ_Object *result;

    if (!(result = GMPy_MPZ_New(context))) {
        return NULL;
    }

    mpz_neg(result->z, MPZ(x));
    return (PyObject*)result;
}
static PyObject *
GMPy_Integer_MinusWithType(PyObject *x, int xtype, CTXT_Object *context)
{
    PyObject *result, *tempx;

    if (!(tempx = (PyObject*)GMPy_MPZ_From_IntegerWithType(x, xtype, context))) {
        return NULL;
    }

    result = _GMPy_MPZ_Minus(tempx, context);
    Py_DECREF(tempx);
    return result;
}

static PyObject *
GMPy_MPZ_Minus_Slot(MPZ_Object *x)
{
    return _GMPy_MPZ_Minus((PyObject*)x, NULL);
}

static PyObject *
_GMPy_MPQ_Minus(PyObject *x, CTXT_Object *context)
{
    MPQ_Object *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPQ_New(context))) {
        return NULL;
    }

    mpq_neg(result->q, MPQ(x));
    return (PyObject*)result;
}

static PyObject *
GMPy_Rational_MinusWithType(PyObject *x, int xtype, CTXT_Object *context)
{
    PyObject *result, *tempx;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPQ_From_RationalWithType(x, xtype, context))) {
        return NULL;
    }

    result = _GMPy_MPQ_Minus(tempx, context);
    Py_DECREF(tempx);
    return result;
}

static PyObject *
GMPy_MPQ_Minus_Slot(MPQ_Object *x)
{
    return _GMPy_MPQ_Minus((PyObject*)x, NULL);
}

static PyObject *
_GMPy_MPFR_Minus(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    mpfr_clear_flags();

    result->rc = mpfr_neg(result->f, MPFR(x), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Real_MinusWithType(PyObject *x, int xtype, CTXT_Object *context)
{
    PyObject *result, *tempx;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPFR_From_RealWithType(x, xtype, 1, context))) {
        return NULL;
    }

    result = _GMPy_MPFR_Minus(tempx, context);
    Py_DECREF(tempx);
    return result;
}

static PyObject *
GMPy_MPFR_Minus_Slot(MPFR_Object *x)
{
    return _GMPy_MPFR_Minus((PyObject*)x, NULL);
}

static PyObject *
_GMPy_MPC_Minus(PyObject *x, CTXT_Object *context)
{
    MPC_Object *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        return NULL;
    }

    result->rc = mpc_neg(result->c, MPC(x), GET_MPC_ROUND(context));
    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Complex_MinusWithType(PyObject *x, int xtype, CTXT_Object *context)
{
    PyObject *result, *tempx;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context))) {
        return NULL;
    }

    result = _GMPy_MPC_Minus(tempx, context);
    Py_DECREF(tempx);
    return result;
}

static PyObject *
GMPy_MPC_Minus_Slot(MPC_Object *x)
{
    return _GMPy_MPC_Minus((PyObject*)x, NULL);
}

static PyObject *
GMPy_Number_Minus(PyObject *x, CTXT_Object *context)
{
    int xtype = GMPy_ObjectType(x);

    if (IS_TYPE_MPZ(xtype))
        return _GMPy_MPZ_Minus(x, context);

    if (IS_TYPE_MPQ(xtype))
        return _GMPy_MPQ_Minus(x, context);

    if (IS_TYPE_MPFR(xtype))
        return _GMPy_MPFR_Minus(x, context);

    if (IS_TYPE_MPC(xtype))
        return _GMPy_MPC_Minus(x, context);

    if (IS_TYPE_INTEGER(xtype))
        return GMPy_Integer_MinusWithType(x, xtype, context);

    if (IS_TYPE_RATIONAL(xtype))
        return GMPy_Rational_MinusWithType(x, xtype, context);

    if (IS_TYPE_REAL(xtype))
        return GMPy_Real_MinusWithType(x, xtype, context);

    if (IS_TYPE_COMPLEX(xtype))
        return GMPy_Complex_MinusWithType(x, xtype, context);

    TYPE_ERROR("minus() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_context_minus,
"minus($self, x, /)\n--\n\n"
"Return -x. The context is applied to the result.");

static PyObject *
GMPy_Context_Minus(PyObject *self, PyObject *args)
{
    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("minus() requires 1 argument.");
        return NULL;
    }

    return GMPy_Number_Minus(PyTuple_GET_ITEM(args, 0), (CTXT_Object*)self);
}
