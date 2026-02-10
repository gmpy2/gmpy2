/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_sign.c                                                            *
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
GMPy_Integer_Sign(PyObject *x, CTXT_Object *context)
{
    long res;
    MPZ_Object *tempx;

    if (!(tempx = GMPy_MPZ_From_Integer(x, context))) {
        return NULL;
    }
    else {
        res = mpz_sgn(tempx->z);
        Py_DECREF((PyObject*)tempx);
        return PyLong_FromLong(res);
    }
}

static PyObject *
GMPy_Rational_Sign(PyObject *x, CTXT_Object *context)
{
    long res;
    MPQ_Object *tempx;

    if (!(tempx = GMPy_MPQ_From_Rational(x, context))) {
        return NULL;
    }
    else {
        res = mpq_sgn(tempx->q);
        Py_DECREF((PyObject*)tempx);
        return PyLong_FromLong(res);
    }
}

static PyObject *
GMPy_Real_Sign(PyObject *x, CTXT_Object *context)
{
    long sign;
    MPFR_Object *tempx;
    PyObject *result;

    CHECK_CONTEXT(context);

    if (!(tempx = GMPy_MPFR_From_Real(x, 1, context))) {
        return NULL;
    }
    else {
        mpfr_clear_flags();
        sign = mpfr_sgn(tempx->f);
        Py_DECREF((PyObject*)tempx);
        result = PyLong_FromLong(sign);
        GMPY_CHECK_ERANGE(result, context, "sign() of invalid value (NaN)");
        return result;
    }
}

static PyObject *
GMPy_Number_Sign(PyObject *x, CTXT_Object *context)
{
    if (IS_INTEGER(x))
        return GMPy_Integer_Sign(x, context);
    else if (IS_RATIONAL_ONLY(x))
        return GMPy_Rational_Sign(x, context);
    else if (IS_REAL_ONLY(x))
        return GMPy_Real_Sign(x, context);

    TYPE_ERROR("sign() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_function_sign,
"sign($module, x, /)\n--\n\n"
"Return -1 if x < 0, 0 if x == 0, or +1 if x >0.");

static PyObject *
GMPy_Context_Sign(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Sign(other, context);
}
