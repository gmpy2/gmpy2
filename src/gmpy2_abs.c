/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_abs.c                                                             *
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

/* This file implements __abs__, gmpy2.abs(), and context.abs().
 *
 * Public API
 * ==========
 * The following function is available as part of GMPY2's C API. If the value
 * of context is NULL, then the function should use the currently active
 * context.
 *
 *   GMPy_Number_Abs(Number, context)
 *
 * Private API
 * ===========
 *   GMPy_MPZ_Abs_Slot
 *   GMPy_MPQ_Abs_Slot
 *   GMPy_MPFR_Abs_Slot
 *   GMPy_MPC_Abs_Slot
 *
 *   GMPy_Integer_AbsWithType(Integer, xtype, context|NULL)
 *   GMPy_Rational_AbsWithType(Rational, xtype, context|NULL)
 *   GMPy_Real_AbsWithType(Real, xtype, context|NULL)
 *   GMPy_Complex_AbsWithType(Complex, xtype, context|NULL)
 *
 *   GMPy_Context_Abs(context, obj)
 */

static PyObject *
GMPy_Integer_AbsWithType(PyObject *x, int xtype, CTXT_Object *context)
{
    MPZ_Object *result = NULL;

    if (IS_TYPE_MPZ(xtype)) {
        if (mpz_sgn(MPZ(x)) >= 0) {
            Py_INCREF(x);
            return x;
        }
        else {
            if ((result = GMPy_MPZ_New(context)))
                mpz_abs(result->z, MPZ(x));
            return (PyObject*)result;
        }
    }

    /* This is safe because result is not an incremented reference to an
     * existing value. Why?
     *   1) No values are interned like Python's integers.
     *   2) MPZ is already handled so GMPy_MPZ_From_IntegerWithType() can't
     *      return an incremented reference to an existing value (which it
     *      would do if passed an MPZ).
     */

    if ((result = GMPy_MPZ_From_IntegerWithType(x, xtype, context))) {
        mpz_abs(result->z, result->z);
    }

    return (PyObject*)result;
}

static PyObject *
GMPy_MPZ_Abs_Slot(MPZ_Object *x)
{
    return GMPy_Integer_AbsWithType((PyObject*)x, OBJ_TYPE_MPZ, NULL);
}

static PyObject *
GMPy_Rational_AbsWithType(PyObject *x, int xtype, CTXT_Object *context)
{
    MPQ_Object *result = NULL;

    if (IS_TYPE_MPQ(xtype)) {
        if (mpz_sgn(mpq_numref(MPQ(x))) >= 0) {
            Py_INCREF(x);
            return x;
        }
        else {
            if ((result = GMPy_MPQ_New(context))) {
                mpq_set(result->q, MPQ(x));
                mpz_abs(mpq_numref(result->q), mpq_numref(result->q));
            }
            return (PyObject*)result;
        }
    }

    /* This is safe because result is not an incremented reference to an
     * existing value. MPQ is already handled so GMPy_MPQ_From_Rational()
     * can't return an incremented reference to an existing value (which it
     * would do if passed an MPQ).
     */

    if ((result = GMPy_MPQ_From_RationalWithType(x, xtype, context))) {
        mpz_abs(mpq_numref(result->q), mpq_numref(result->q));
    }

    return (PyObject*)result;
}

static PyObject *
GMPy_MPQ_Abs_Slot(MPQ_Object *x)
{
    return GMPy_Rational_AbsWithType((PyObject*)x, OBJ_TYPE_MPQ, NULL);
}

static PyObject *
GMPy_Real_AbsWithType(PyObject *x, int xtype, CTXT_Object *context)
{
    MPFR_Object *result = NULL, *tempx = NULL;
    mpfr_exp_t _oldemin, _oldemax;

    CHECK_CONTEXT(context);
    
    if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context)) ||
        !(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)result);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_clear_flags();
    
    _oldemin = mpfr_get_emin();
    _oldemax = mpfr_get_emax();
	mpfr_set_emin(GET_EMIN(context));
	mpfr_set_emax(GET_EMAX(context));

    result->rc = mpfr_abs(result->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);

	mpfr_set_emin(_oldemin);
	mpfr_set_emax(_oldemax);
	
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_MPFR_Abs_Slot(MPFR_Object *x)
{
    return GMPy_Real_AbsWithType((PyObject*)x, OBJ_TYPE_MPFR, NULL);
}

static PyObject *
GMPy_Complex_AbsWithType(PyObject *x, int xtype, CTXT_Object *context)
{
    MPFR_Object *result = NULL;
    MPC_Object *tempx = NULL;

    CHECK_CONTEXT(context);

    if (!(tempx = GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context)) ||
        !(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)result);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_clear_flags();

    result->rc = mpc_abs(result->f, tempx->c, GET_MPC_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_MPC_Abs_Slot(MPC_Object *x)
{
    return GMPy_Complex_AbsWithType((PyObject*)x, OBJ_TYPE_MPC, NULL);
}

static PyObject *
GMPy_Number_Abs(PyObject *x, CTXT_Object *context)
{
    int xtype = GMPy_ObjectType(x);
    
    if (IS_TYPE_INTEGER(xtype))
        return GMPy_Integer_AbsWithType(x, xtype, context);

    if (IS_TYPE_RATIONAL(xtype))
        return GMPy_Rational_AbsWithType(x, xtype, context);

    if (IS_TYPE_REAL(xtype))
        return GMPy_Real_AbsWithType(x, xtype, context);

    if (IS_TYPE_COMPLEX(xtype))
        return GMPy_Complex_AbsWithType(x, xtype, context);

    TYPE_ERROR("abs() argument type not supported");
    return NULL;
}

/* Implement context.abs(). The following code assumes it used a as method of
 * a context. */

PyDoc_STRVAR(GMPy_doc_context_abs,
"abs($self, x, /)\n--\n\n"
"Return abs(x), the context is applied to the result.");

static PyObject *
GMPy_Context_Abs(PyObject *self, PyObject *other)
{
    return GMPy_Number_Abs(other, (CTXT_Object*)self);
}

