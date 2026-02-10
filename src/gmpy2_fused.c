/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_fused.c                                                           *
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
_GMPy_MPZ_FMA(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    MPZ_Object *result = NULL;

    if (!(result = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
    mpz_mul(result->z, MPZ(x), MPZ(y));
    mpz_add(result->z, result->z, MPZ(z));
    GMPY_MAYBE_END_ALLOW_THREADS(context);
    return (PyObject*)result;
}

static PyObject *
GMPy_IntegerWithType_FMA(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL;

    if (!(tempx = (PyObject*)GMPy_MPZ_From_IntegerWithType(x, xtype, context)) ||
        !(tempy = (PyObject*)GMPy_MPZ_From_IntegerWithType(y, ytype, context)) ||
        !(tempz = (PyObject*)GMPy_MPZ_From_IntegerWithType(z, ztype, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPZ_FMA(tempx, tempy, tempz, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    return result;
}

static PyObject *
_GMPy_MPQ_FMA(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    MPQ_Object *result;

    if (!(result = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
    mpq_mul(result->q, MPQ(x), MPQ(y));
    mpq_add(result->q, result->q, MPQ(z));
    GMPY_MAYBE_END_ALLOW_THREADS(context);
    return (PyObject*)result;
}

static PyObject *
GMPy_RationalWithType_FMA(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL;

    if (!(tempx = (PyObject*)GMPy_MPQ_From_RationalWithType(x, xtype, context)) ||
        !(tempy = (PyObject*)GMPy_MPQ_From_RationalWithType(y, ytype, context)) ||
        !(tempz = (PyObject*)GMPy_MPQ_From_RationalWithType(z, ztype, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPQ_FMA(tempx, tempy, tempz, context);
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)tempz);
    return (result);
}

static PyObject *
_GMPy_MPFR_FMA(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    MPFR_Object *result;

    if (!(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_clear_flags();

    result->rc = mpfr_fma(result->f, MPFR(x), MPFR(y), MPFR(z), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_RealWithType_FMA(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL;

    if (!(tempx = (PyObject*)GMPy_MPFR_From_RealWithType(x, xtype, 1, context)) ||
        !(tempy = (PyObject*)GMPy_MPFR_From_RealWithType(y, ytype, 1, context)) ||
        !(tempz = (PyObject*)GMPy_MPFR_From_RealWithType(z, ztype, 1, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPFR_FMA(tempx, tempy, tempz, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    return result;
}

static PyObject *
_GMPy_MPC_FMA(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    MPC_Object *result;

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result->rc = mpc_fma(result->c, MPC(x), MPC(y), MPC(z), GET_MPC_ROUND(context));
    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_ComplexWithType_FMA(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL;

    if (!(tempx = (PyObject*)GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context)) ||
        !(tempy = (PyObject*)GMPy_MPC_From_ComplexWithType(y, ytype, 1, 1, context)) ||
        !(tempz = (PyObject*)GMPy_MPC_From_ComplexWithType(z, ztype, 1, 1, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPC_FMA(tempx, tempy, tempz, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    return result;
}

PyDoc_STRVAR(GMPy_doc_context_fma,
"fma($self, x, y, z, /)\n--\n\n"
"Return correctly rounded result of (x * y) + z.");

PyDoc_STRVAR(GMPy_doc_function_fma,
"fma($module, x, y, z, /)\n--\n\n"
"Return correctly rounded result of (x * y) + z.");

GMPY_MPFR_MPC_TRIOP_TEMPLATEWT(FMA, fma)

static PyObject *
_GMPy_MPZ_FMS(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    MPZ_Object *result;

    if (!(result = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
    mpz_mul(result->z, MPZ(x), MPZ(y));
    mpz_sub(result->z, result->z, MPZ(z));
    GMPY_MAYBE_END_ALLOW_THREADS(context);
    return (PyObject*)result;
}

static PyObject *
GMPy_IntegerWithType_FMS(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL;

    if (!(tempx = (PyObject*)GMPy_MPZ_From_IntegerWithType(x, xtype, context)) ||
        !(tempy = (PyObject*)GMPy_MPZ_From_IntegerWithType(y, ytype, context)) ||
        !(tempz = (PyObject*)GMPy_MPZ_From_IntegerWithType(z, ztype, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPZ_FMS(tempx, tempy, tempz, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    return result;
}

static PyObject *
_GMPy_MPQ_FMS(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    MPQ_Object *result;

    if (!(result = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
    mpq_mul(result->q, MPQ(x), MPQ(y));
    mpq_sub(result->q, result->q, MPQ(z));
    GMPY_MAYBE_END_ALLOW_THREADS(context);
    return (PyObject*)result;
}

static PyObject *
GMPy_RationalWithType_FMS(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL;

    if (!(tempx = (PyObject*)GMPy_MPQ_From_RationalWithType(x, xtype, context)) ||
        !(tempy = (PyObject*)GMPy_MPQ_From_RationalWithType(y, ytype, context)) ||
        !(tempz = (PyObject*)GMPy_MPQ_From_RationalWithType(z, ztype, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPQ_FMS(tempx, tempy, tempz, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    return result;
}

static PyObject *
_GMPy_MPFR_FMS(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    MPFR_Object *result;

    if (!(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_clear_flags();

    result->rc = mpfr_fms(result->f, MPFR(x), MPFR(y), MPFR(z), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_RealWithType_FMS(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL;

    if (!(tempx = (PyObject*)GMPy_MPFR_From_RealWithType(x, xtype, 1, context)) ||
        !(tempy = (PyObject*)GMPy_MPFR_From_RealWithType(y, ytype, 1, context)) ||
        !(tempz = (PyObject*)GMPy_MPFR_From_RealWithType(z, ztype, 1, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPFR_FMS(tempx, tempy, tempz, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    return result;
}

static PyObject *
_GMPy_MPC_FMS(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    MPC_Object *result = NULL, *temp = NULL;

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    temp = GMPy_MPC_From_ComplexWithTypeAndCopy(z, OBJ_TYPE_MPC, 1, 1, context);

    mpc_neg(MPC(temp), MPC(temp), GET_MPC_ROUND(context));
    result->rc = mpc_fma(result->c, MPC(x), MPC(y), MPC(temp), GET_MPC_ROUND(context));
    Py_DECREF(temp);
    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_ComplexWithType_FMS(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, int ztype, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL;

    if (!(tempx = (PyObject*)GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context)) ||
        !(tempy = (PyObject*)GMPy_MPC_From_ComplexWithType(y, ytype, 1, 1, context)) ||
        !(tempz = (PyObject*)GMPy_MPC_From_ComplexWithType(z, ztype, 1, 1, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPC_FMS(tempx, tempy, tempz, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    return result;
}
PyDoc_STRVAR(GMPy_doc_context_fms,
"fms($self, x, y, z, /)\n--\n\n"
"Return correctly rounded result of (x * y) - z.");

PyDoc_STRVAR(GMPy_doc_function_fms,
"fms($module, x, y, z, /)\n--\n\n"
"Return correctly rounded result of (x * y) - z.");

GMPY_MPFR_MPC_TRIOP_TEMPLATEWT(FMS, fms)

/* Add support for new fmma and fmms functions from MPFR 4. */\

#if MPFR_VERSION_MAJOR > 3

static PyObject *
_GMPy_MPZ_FMMA(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    MPZ_Object *result = NULL, *temp = NULL;

    if (!(result = GMPy_MPZ_New(context)) ||
        !(temp = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(result);
        Py_XDECREF(temp);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
    mpz_mul(result->z, MPZ(x), MPZ(y));
    mpz_mul(temp->z, MPZ(z), MPZ(t));
    mpz_add(result->z, result->z, temp->z);
    GMPY_MAYBE_END_ALLOW_THREADS(context);
    Py_DECREF(temp);
    return (PyObject*)result;
}

static PyObject *
GMPy_IntegerWithType_FMMA(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, 
                          int ztype, PyObject *t, int ttype, CTXT_Object *context)
{
    PyObject *result = NULL, *tempx = NULL, *tempy = NULL, *tempz = NULL, *tempt = NULL;

    if (!(tempx = (PyObject*)GMPy_MPZ_From_IntegerWithType(x, xtype, context)) ||
        !(tempy = (PyObject*)GMPy_MPZ_From_IntegerWithType(y, ytype, context)) ||
        !(tempz = (PyObject*)GMPy_MPZ_From_IntegerWithType(z, ztype, context)) ||
        !(tempt = (PyObject*)GMPy_MPZ_From_IntegerWithType(t, ttype, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        Py_XDECREF(tempt);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPZ_FMMA(tempx, tempy, tempz, tempt, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    Py_DECREF(tempt);
    return result;
}

static PyObject *
_GMPy_MPQ_FMMA(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    MPQ_Object *result = NULL, *temp = NULL;

    if (!(result = GMPy_MPQ_New(context)) ||
        !(temp = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(result);
        Py_XDECREF(temp);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
    mpq_mul(result->q, MPQ(x), MPQ(y));
    mpq_mul(temp->q, MPQ(z), MPQ(t));
    mpq_add(result->q, result->q, temp->q);
    GMPY_MAYBE_END_ALLOW_THREADS(context);
    Py_DECREF(temp);
    return (PyObject*)result;
}

static PyObject *
GMPy_RationalWithType_FMMA(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, 
                           int ztype, PyObject *t, int ttype, CTXT_Object *context)
{
    PyObject *result = NULL, *tempx = NULL, *tempy = NULL, *tempz = NULL, *tempt = NULL;

    if (!(tempx = (PyObject*)GMPy_MPQ_From_RationalWithType(x, xtype, context)) ||
        !(tempy = (PyObject*)GMPy_MPQ_From_RationalWithType(y, ytype, context)) ||
        !(tempz = (PyObject*)GMPy_MPQ_From_RationalWithType(z, ztype, context)) ||
        !(tempt = (PyObject*)GMPy_MPQ_From_RationalWithType(t, ttype, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        Py_XDECREF(tempt);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPQ_FMMA(tempx, tempy, tempz, tempt, context);
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)tempz);
    Py_DECREF((PyObject*)tempt);
    return (result);
}

static PyObject *
_GMPy_MPFR_FMMA(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    MPFR_Object *result;

    if (!(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_clear_flags();
    
    result->rc = mpfr_fmma(result->f, MPFR(x), MPFR(y), MPFR(z), MPFR(t), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_RealWithType_FMMA(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, 
                       int ztype, PyObject *t, int ttype, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL, *tempt = NULL;

    if (!(tempx = (PyObject*)GMPy_MPFR_From_RealWithType(x, xtype, 1, context)) ||
        !(tempy = (PyObject*)GMPy_MPFR_From_RealWithType(y, ytype, 1, context)) ||
        !(tempz = (PyObject*)GMPy_MPFR_From_RealWithType(z, ztype, 1, context)) ||
        !(tempt = (PyObject*)GMPy_MPFR_From_RealWithType(t, ttype, 1, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        Py_XDECREF(tempt);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPFR_FMMA(tempx, tempy, tempz, tempt, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    Py_DECREF(tempt);
    return result;
}

PyDoc_STRVAR(GMPy_doc_context_fmma,
"fmma($self, x, y, z, t, /)\n--\n\n"
"Return correctly rounded result of (x * y) + (z * t).");

PyDoc_STRVAR(GMPy_doc_function_fmma,
"fmma($module, x, y, z, t, /)\n--\n\n"
"Return correctly rounded result of (x * y) + (z + t).");

GMPY_MPFR_QUADOP_TEMPLATEWT(FMMA, fmma)

static PyObject *
_GMPy_MPZ_FMMS(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    MPZ_Object *result = NULL, *temp = NULL;

    if (!(result = GMPy_MPZ_New(context)) ||
        !(temp = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(result);
        Py_XDECREF(temp);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
    mpz_mul(result->z, MPZ(x), MPZ(y));
    mpz_mul(temp->z, MPZ(z), MPZ(t));
    mpz_sub(result->z, result->z, temp->z);
    GMPY_MAYBE_END_ALLOW_THREADS(context);
    Py_DECREF(temp);
    return (PyObject*)result;
}

static PyObject *
GMPy_IntegerWithType_FMMS(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, 
                          int ztype, PyObject *t, int ttype, CTXT_Object *context)
{
    PyObject *result = NULL, *tempx = NULL, *tempy = NULL, *tempz = NULL, *tempt = NULL;

    if (!(tempx = (PyObject*)GMPy_MPZ_From_IntegerWithType(x, xtype, context)) ||
        !(tempy = (PyObject*)GMPy_MPZ_From_IntegerWithType(y, ytype, context)) ||
        !(tempz = (PyObject*)GMPy_MPZ_From_IntegerWithType(z, ztype, context)) ||
        !(tempt = (PyObject*)GMPy_MPZ_From_IntegerWithType(t, ttype, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        Py_XDECREF(tempt);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPZ_FMMS(tempx, tempy, tempz, tempt, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    Py_DECREF(tempt);
    return result;
}

static PyObject *
_GMPy_MPQ_FMMS(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    MPQ_Object *result = NULL, *temp = NULL;

    if (!(result = GMPy_MPQ_New(context)) ||
        !(temp = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(result);
        Py_XDECREF(temp);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
    mpq_mul(result->q, MPQ(x), MPQ(y));
    mpq_mul(temp->q, MPQ(z), MPQ(t));
    mpq_sub(result->q, result->q, temp->q);
    GMPY_MAYBE_END_ALLOW_THREADS(context);
    Py_DECREF(temp);
    return (PyObject*)result;
}

static PyObject *
GMPy_RationalWithType_FMMS(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, 
                           int ztype, PyObject *t, int ttype, CTXT_Object *context)
{
    PyObject *result = NULL, *tempx = NULL, *tempy = NULL, *tempz = NULL, *tempt = NULL;

    if (!(tempx = (PyObject*)GMPy_MPQ_From_RationalWithType(x, xtype, context)) ||
        !(tempy = (PyObject*)GMPy_MPQ_From_RationalWithType(y, ytype, context)) ||
        !(tempz = (PyObject*)GMPy_MPQ_From_RationalWithType(z, ztype, context)) ||
        !(tempt = (PyObject*)GMPy_MPQ_From_RationalWithType(t, ttype, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        Py_XDECREF(tempt);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPQ_FMMS(tempx, tempy, tempz, tempt, context);
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)tempz);
    Py_DECREF((PyObject*)tempt);
    return (result);
}

static PyObject *
_GMPy_MPFR_FMMS(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    MPFR_Object *result;

    if (!(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_clear_flags();

    result->rc = mpfr_fmms(result->f, MPFR(x), MPFR(y), MPFR(z), MPFR(t), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_RealWithType_FMMS(PyObject *x, int xtype, PyObject *y, int ytype, PyObject *z, 
                       int ztype, PyObject *t, int ttype, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL, *tempt = NULL;

    if (!(tempx = (PyObject*)GMPy_MPFR_From_RealWithType(x, xtype, 1, context)) ||
        !(tempy = (PyObject*)GMPy_MPFR_From_RealWithType(y, ytype, 1, context)) ||
        !(tempz = (PyObject*)GMPy_MPFR_From_RealWithType(z, ztype, 1, context)) ||
        !(tempt = (PyObject*)GMPy_MPFR_From_RealWithType(t, ttype, 1, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        Py_XDECREF(tempt);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPFR_FMMS(tempx, tempy, tempz, tempt, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    Py_DECREF(tempt);
    return result;
}

PyDoc_STRVAR(GMPy_doc_context_fmms,
"fmms($self, x, y, z, t, /)\n--\n\n"
"Return correctly rounded result of (x * y) - (z * t).");

PyDoc_STRVAR(GMPy_doc_function_fmms,
"fmms($module, x, y, z, t, /)\n--\n\n"
"Return correctly rounded result of (x * y) - (z + t).");

GMPY_MPFR_QUADOP_TEMPLATEWT(FMMS, fmms)

#endif
