/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_square.c                                                             *
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

/* Public API
 * ==========
 * The following function is available as part of GMPY2's C API. A NULL value
 * for context implies the function should use the currently active context.
 *
 *   GMPy_Number_Square(Number, Number, context|NULL)
 *
 * Private API
 * ===========
 *   GMPy_Integer_Square(Integer, Integer, context|NULL)
 *   GMPy_Rational_Square(Rational, Rational, context|NULL)
 *   GMPy_Real_Square(Real, Real, context|NULL)
 *   GMPy_Complex_Square(Complex, Complex, context|NULL)
 *
 *   GMPy_Context_Square(context, args)
 *
 */

static PyObject *
_GMPy_MPZ_Square(PyObject *x, CTXT_Object *context)
{
    MPZ_Object *result = NULL;

    if (!(result = GMPy_MPZ_New(context))) {
        return NULL;
    }

    mpz_mul(result->z, MPZ(x), MPZ(x));
    return (PyObject*)result;
}

static PyObject *
GMPy_Integer_Square(PyObject *x, CTXT_Object *context)
{
    PyObject *result, *tempx;

    if (!(tempx = (PyObject*)GMPy_MPZ_From_Integer(x, context))) {
        return NULL;
    }

    result = _GMPy_MPZ_Square(tempx, context);
    Py_DECREF(tempx);
    return result;
}

static PyObject *
_GMPy_MPQ_Square(PyObject *x, CTXT_Object *context)
{
    MPQ_Object *result;

    if (!(result = GMPy_MPQ_New(context))) {
        return NULL;
    }

    mpq_mul(result->q, MPQ(x), MPQ(x));
    return (PyObject*)result;
}

static PyObject *
GMPy_Rational_Square(PyObject *x, CTXT_Object *context)
{
    PyObject *result, *tempx;

    if (!(tempx = (PyObject*)GMPy_MPQ_From_Rational(x, context))) {
        return NULL;
    }

    result = _GMPy_MPQ_Square(tempx, context);
    Py_DECREF(tempx);
    return result;
}

static PyObject *
_GMPy_MPFR_Square(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    mpfr_clear_flags();

    mpfr_sqr(result->f, MPFR(x), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Real_Square(PyObject *x, CTXT_Object *context)
{
    PyObject *result, *tempx;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPFR_From_Real(x, 1, context))) {
        return NULL;
    }

    result = _GMPy_MPFR_Square(tempx, context);
    Py_DECREF(tempx);
    return result;
}

static PyObject *
_GMPy_MPC_Square(PyObject *x, CTXT_Object *context)
{
    MPC_Object *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        return NULL;
    }

    mpc_sqr(result->c, MPC(x), GET_MPC_ROUND(context));
    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}
static PyObject *
GMPy_Complex_Square(PyObject *x, CTXT_Object *context)
{
    PyObject *result, *tempx;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPC_From_Complex(x, 1, 1, context))) {
        return NULL;
    }

    result = _GMPy_MPC_Square(tempx, context);
    Py_DECREF(tempx);
    return result;
}

PyDoc_STRVAR(GMPy_doc_function_square,
"square($module, x, /)\n--\n\n"
"Return x * x.");

PyDoc_STRVAR(GMPy_doc_context_square,
"square($self, x, /)\n--\n\n"
"Return x * x.");

static PyObject *
GMPy_Number_Square(PyObject *x, CTXT_Object *context)
{
    if (MPZ_Check(x))
        return _GMPy_MPZ_Square(x, context);

    if (MPQ_Check(x))
        return _GMPy_MPQ_Square(x, context);

    if (MPFR_Check(x))
        return _GMPy_MPFR_Square(x, context);

    if (MPC_Check(x))
        return _GMPy_MPC_Square(x, context);

    if (IS_INTEGER(x))
        return GMPy_Integer_Square(x, context);

    if (IS_RATIONAL(x))
        return GMPy_Rational_Square(x, context);

    if (IS_REAL(x))
        return GMPy_Real_Square(x, context);

    if (IS_COMPLEX(x))
        return GMPy_Complex_Square(x, context);

    TYPE_ERROR("square() argument type not supported");
    return NULL;
}

static PyObject *
GMPy_Context_Square(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Square(other, context);
}
