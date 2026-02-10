/* * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * * * * * *
 * gmpy2_xmpz_inplace.c                                                     *
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


/* Provides inplace mutating operations for xmpz. */

#include <math.h>

/* Inplace xmpz addition. */

static PyObject *
GMPy_XMPZ_IAdd_Slot(PyObject *self, PyObject *other)
{
    /* Try to make mpz + small_int faster */

    CTXT_Object *context = NULL;
    CHECK_CONTEXT(context);

    int ytype = GMPy_ObjectType(other);

    if (IS_TYPE_PyInteger(ytype)) {
        int error;
        long temp = PyLong_AsLongAndOverflow(other, &error);

        if (!error) {
            if (temp >= 0) {
                mpz_add_ui(MPZ(self), MPZ(self), temp);
            }
            else {
                mpz_sub_ui(MPZ(self), MPZ(self), -temp);
            }
        }
        else {
            mpz_t tempz;
            mpz_init(tempz);
            if (mpz_set_PyLong(tempz, other)) {
                /* LCOV_EXCL_START */
                mpz_clear(tempz);
                return NULL;
                /* LCOV_EXCL_STOP */
            }
            GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
            mpz_add(MPZ(self), MPZ(self), tempz);
            GMPY_MAYBE_END_ALLOW_THREADS(context);
            mpz_clear(tempz);
        }
        Py_INCREF(self);
        return self;
    }

    if (IS_TYPE_MPZANY(ytype)) {
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_add(MPZ(self), MPZ(self), MPZ(other));
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace mpz subtraction.
 */

static PyObject *
GMPy_XMPZ_ISub_Slot(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;
    CHECK_CONTEXT(context);

    int ytype = GMPy_ObjectType(other);

    if (IS_TYPE_PyInteger(ytype)) {
        int error;
        long temp = PyLong_AsLongAndOverflow(other, &error);

        if (!error) {
            if (temp >= 0) {
                mpz_sub_ui(MPZ(self), MPZ(self), temp);
            }
            else {
                mpz_add_ui(MPZ(self), MPZ(self), -temp);
            }
        }
        else {
            mpz_t tempz;
            mpz_init(tempz);
            if (mpz_set_PyLong(tempz, other)) {
                /* LCOV_EXCL_START */
                mpz_clear(tempz);
                return NULL;
                /* LCOV_EXCL_STOP */
            }
            GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
            mpz_sub(MPZ(self), MPZ(self), tempz);
            GMPY_MAYBE_END_ALLOW_THREADS(context);
            mpz_clear(tempz);
        }
        Py_INCREF(self);
        return self;
    }

    if (IS_TYPE_MPZANY(ytype)) {
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_sub(MPZ(self), MPZ(self), MPZ(other));
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz multiplication.
 */

static PyObject *
GMPy_XMPZ_IMul_Slot(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;
    CHECK_CONTEXT(context);

    int ytype = GMPy_ObjectType(other);

    if (IS_TYPE_PyInteger(ytype)) {
        int error;
        long temp = PyLong_AsLongAndOverflow(other, &error);

        if (!error) {
            mpz_mul_si(MPZ(self), MPZ(self), temp);
        }
        else {
            mpz_t tempz;
            mpz_init(tempz);
            if (mpz_set_PyLong(tempz, other)) {
                /* LCOV_EXCL_START */
                mpz_clear(tempz);
                return NULL;
                /* LCOV_EXCL_STOP */
            }
            GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
            mpz_mul(MPZ(self), MPZ(self), tempz);
            GMPY_MAYBE_END_ALLOW_THREADS(context);
            mpz_clear(tempz);
        }
        Py_INCREF(self);
        return self;
    }

    if (IS_TYPE_MPZANY(ytype)) {
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_mul(MPZ(self), MPZ(self), MPZ(other));
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Pympany_floordiv follows the // semantics from Python 3.x. The result is
 * an mpz when the arguments are mpz or mpq, but the result is an mpf when
 * the arguments are mpf.
 */

static PyObject *
GMPy_XMPZ_IFloorDiv_Slot(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;
    CHECK_CONTEXT(context);

    int ytype = GMPy_ObjectType(other);

    if (IS_TYPE_PyInteger(ytype)) {
        int error;
        long temp = PyLong_AsLongAndOverflow(other, &error);

        if (!error) {
            if (temp == 0) {
                ZERO_ERROR("xmpz division by zero");
                return NULL;
            }
            else if (temp > 0) {
                mpz_fdiv_q_ui(MPZ(self), MPZ(self), temp);
            }
            else {
                mpz_cdiv_q_ui(MPZ(self), MPZ(self), -temp);
                mpz_neg(MPZ(self), MPZ(self));
            }
        }
        else {
            mpz_t tempz;
            mpz_init(tempz);
            if (mpz_set_PyLong(tempz, other)) {
                /* LCOV_EXCL_START */
                mpz_clear(tempz);
                return NULL;
                /* LCOV_EXCL_STOP */
            }
            GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
            mpz_fdiv_q(MPZ(self), MPZ(self), tempz);
            GMPY_MAYBE_END_ALLOW_THREADS(context);
            mpz_clear(tempz);
        }
        Py_INCREF(self);
        return self;
    }

    if (IS_TYPE_MPZANY(ytype)) {
        if (mpz_sgn(MPZ(other)) == 0) {
            ZERO_ERROR("xmpz division by zero");
            return NULL;
        }
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_fdiv_q(MPZ(self), MPZ(self), MPZ(other));
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz remainder.
 */

static PyObject *
GMPy_XMPZ_IRem_Slot(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;
    CHECK_CONTEXT(context);

    int ytype = GMPy_ObjectType(other);

    if (IS_TYPE_PyInteger(ytype)) {
        int error;
        long temp = PyLong_AsLongAndOverflow(other, &error);

        if (!error) {
            if (temp > 0) {
                mpz_fdiv_r_ui(MPZ(self), MPZ(self), temp);
            }
            else if (temp == 0) {
                ZERO_ERROR("xmpz modulo by zero");
                return NULL;
            }
            else {
                mpz_cdiv_r_ui(MPZ(self), MPZ(self), -temp);
            }
        }
        else {
            mpz_t tempz;
            mpz_init(tempz);
            if (mpz_set_PyLong(tempz, other)) {
                /* LCOV_EXCL_START */
                mpz_clear(tempz);
                return NULL;
                /* LCOV_EXCL_STOP */
            }
            GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
            mpz_fdiv_r(MPZ(self), MPZ(self), tempz);
            GMPY_MAYBE_END_ALLOW_THREADS(context);
            mpz_clear(tempz);
        }
        Py_INCREF(self);
        return self;
    }

    if (IS_TYPE_MPZANY(ytype)) {
        if(mpz_sgn(MPZ(other)) == 0) {
            ZERO_ERROR("xmpz modulo by zero");
            return NULL;
        }
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_fdiv_r(MPZ(self), MPZ(self), MPZ(other));
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz rshift.
 */

static PyObject *
GMPy_XMPZ_IRshift_Slot(PyObject *self, PyObject *other)
{
    mp_bitcnt_t shift = GMPy_Integer_AsMpBitCnt(other);
    if (shift == (mp_bitcnt_t)(-1) && PyErr_Occurred())
        return NULL;

    mpz_fdiv_q_2exp(MPZ(self), MPZ(self), shift);
    Py_INCREF(self);
    return self;

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz lshift.
 */

static PyObject *
GMPy_XMPZ_ILshift_Slot(PyObject *self, PyObject *other)
{
    mp_bitcnt_t shift = GMPy_Integer_AsMpBitCnt(other);
    if (shift == (mp_bitcnt_t)(-1) && PyErr_Occurred())
        return NULL;

    mpz_mul_2exp(MPZ(self), MPZ(self), shift);
    Py_INCREF(self);
    return self;

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz_pow.
 */

static PyObject *
GMPy_XMPZ_IPow_Slot(PyObject *self, PyObject *other, PyObject *mod)
{
    unsigned long exp = GMPy_Integer_AsUnsignedLong(other);
    if (exp == (unsigned long)(-1) && PyErr_Occurred())
        return NULL;

    mpz_pow_ui(MPZ(self), MPZ(self), exp);
    Py_INCREF((PyObject*)self);
    return self;

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz and.
 */

static PyObject *
GMPy_XMPZ_IAnd_Slot(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;
    CHECK_CONTEXT(context);

    if (CHECK_MPZANY(other)) {
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_and(MPZ(self), MPZ(self), MPZ(other));
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_INCREF(self);
        return self;
    }

    if (PyLong_Check(other)) {
        mpz_t tempz;
        mpz_init(tempz);
        if (mpz_set_PyLong(tempz, other)) {
            /* LCOV_EXCL_START */
            mpz_clear(tempz);
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_and(MPZ(self), MPZ(self), tempz);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        mpz_clear(tempz);
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz xor.
 */

static PyObject *
GMPy_XMPZ_IXor_Slot(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;
    CHECK_CONTEXT(context);

    if(CHECK_MPZANY(other)) {
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_xor(MPZ(self), MPZ(self), MPZ(other));
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_INCREF(self);
        return self;
    }

    if(PyLong_Check(other)) {
        mpz_t tempz;
        mpz_init(tempz);
        if (mpz_set_PyLong(tempz, other)) {
            /* LCOV_EXCL_START */
            mpz_clear(tempz);
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_xor(MPZ(self), MPZ(self), tempz);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        mpz_clear(tempz);
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz or.
 */

static PyObject *
GMPy_XMPZ_IIor_Slot(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;
    CHECK_CONTEXT(context);

    if(CHECK_MPZANY(other)) {
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_ior(MPZ(self), MPZ(self), MPZ(other));
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_INCREF(self);
        return self;
    }

    if(PyLong_Check(other)) {
        mpz_t tempz;
        mpz_init(tempz);
        if (mpz_set_PyLong(tempz, other)) {
            /* LCOV_EXCL_START */
            mpz_clear(tempz);
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_ior(MPZ(self), MPZ(self), tempz);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        mpz_clear(tempz);
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}
