/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_pow.c                                                             *
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


#include "pythoncapi_compat.h"

/* This file implements the ** operator, Python's pow() function,
 * gmpy2.powmod(), and context.pow().
 */


/* Pympz_Pow_Integer is called by GMPy_Number_Pow() after verifying that the
 * first two arguments are integers, but not necessarily mpz. The third
 * argument must either be an integer or Py_None. The context argument is not
 * currently used but may be used in the future.
 */

static PyObject *
GMPy_Integer_PowWithType(PyObject *b, int btype, PyObject *e, int etype,
                         PyObject *m, CTXT_Object *context)
{
    MPZ_Object *result = NULL, *tempb = NULL, *tempe = NULL, *tempm = NULL;
    int has_mod, mtype;

    CHECK_CONTEXT(context);

    /* Try to parse the modulus value first. */

    if (Py_IsNone(m)) {
        has_mod = 0;
    }
    else {
        has_mod = 1;
        mtype = GMPy_ObjectType(m);
        if (!IS_TYPE_INTEGER(mtype)) {
            TYPE_ERROR("pow() modulus must be an integer");
            return NULL;
        }
        else {
            if (!(tempm = GMPy_MPZ_From_IntegerWithType(m, mtype, context))) {
                goto err;
            }
        }
    }

    if (!(result = GMPy_MPZ_New(context)) ||
        !(tempb = GMPy_MPZ_From_IntegerWithType(b, btype, context)) ||
        !(tempe = GMPy_MPZ_From_IntegerWithType(e, etype, context))) {
        goto err;
    }

    if (!has_mod) {
        /* When no modulo is present, the exponent must fit in unsigned long.
         */
        unsigned long el;

        if (mpz_sgn(tempe->z) < 0) {
            Py_DECREF((PyObject*)result);
            Py_DECREF((PyObject*)tempb);
            Py_DECREF((PyObject*)tempe);
            /* This should return an mpfr result. */
            return GMPy_Real_PowWithType(b, btype, e, etype, m, context);
        }

        /* Catch -1, 0, 1 getting raised to large exponents. */

        if (mpz_cmp_si(tempb->z, 0) == 0) {
            if (mpz_cmp_si(tempe->z, 0) == 0) {
                mpz_set_ui(result->z, 1);
            }
            else {
                mpz_set_ui(result->z, 0);
            }
            goto done;
        }

        if (mpz_cmp_si(tempb->z, 1) == 0) {
            mpz_set_ui(result->z, 1);
            goto done;
        }

        if (mpz_cmp_si(tempb->z, -1) == 0) {
            if (mpz_even_p(tempe->z)) {
                mpz_set_ui(result->z, 1);
            }
            else {
                mpz_set_si(result->z, -1);
            }
            goto done;
        }

        if (!mpz_fits_ulong_p(tempe->z)) {
            VALUE_ERROR("pow() outrageous exponent");
            goto err;
        }

        el = (unsigned long) mpz_get_ui(tempe->z);
        mpz_pow_ui(result->z, tempb->z, el);
        goto done;
    }
    else {
        /* Modulo is present. */
        int sign, has_inverse;
        mpz_t mm, base, exp;

        sign = mpz_sgn(tempm->z);
        if (sign == 0) {
            VALUE_ERROR("pow() 3rd argument cannot be 0");
            goto err;
        }

        mpz_init(mm);
        mpz_abs(mm, tempm->z);

        /* A negative exponent is allowed if inverse exists. */
        if (mpz_sgn(tempe->z) < 0) {
            mpz_init(base);
            mpz_init(exp);


            has_inverse = mpz_invert(base, tempb->z, mm);
            if (has_inverse) {
                mpz_abs(exp, tempe->z);
                mpz_powm(result->z, base, exp, mm);
            }
            mpz_clear(base);
            mpz_clear(exp);
            mpz_clear(mm);

            /* Python uses a rather peculiar convention for negative modulos
            * If the modulo is negative, result should be in the interval
            * m < r <= 0 .
            */
            if ((sign < 0) && (mpz_sgn(result->z) > 0)) {
                mpz_add(result->z, result->z, tempm->z);
            }


            if (!has_inverse) {
                VALUE_ERROR("pow() base not invertible");
                goto err;
            }
        }
        else {
            GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
            mpz_powm(result->z, tempb->z, tempe->z, mm);
            mpz_clear(mm);

            /* Python uses a rather peculiar convention for negative modulos
            * If the modulo is negative, result should be in the interval
            * m < r <= 0 .
            */
            if ((sign < 0) && (mpz_sgn(result->z) > 0))
                mpz_add(result->z, result->z, tempm->z);

            GMPY_MAYBE_END_ALLOW_THREADS(context);
        }
    }

  done:
    Py_XDECREF((PyObject*)tempb);
    Py_XDECREF((PyObject*)tempe);
    Py_XDECREF((PyObject*)tempm);
    return (PyObject*)result;

  err:
    Py_XDECREF((PyObject*)tempb);
    Py_XDECREF((PyObject*)tempe);
    Py_XDECREF((PyObject*)tempm);
    Py_DECREF((PyObject*)result);
    return NULL;
}

static PyObject *
GMPy_Integer_PowModBaseListWithType(PyObject *base_lst,
                                    PyObject *e, int etype,
                                    PyObject *m, int mtype)
{
    MPZ_Object *tempe = NULL, *tempm = NULL, *tempres = NULL;
    PyObject *result = NULL, *temp = NULL;
    Py_ssize_t i, seq_length;

    if (!(tempm = GMPy_MPZ_From_IntegerWithType(m, mtype, NULL)) ||
        !(tempe = GMPy_MPZ_From_IntegerWithType(e, etype, NULL))) {
        return NULL;
    }

    if (mpz_sgn(tempm->z) < 1) {
        VALUE_ERROR("powmod_base_list() 'mod' must be > 0");
        Py_DECREF((PyObject*)tempm);
        Py_DECREF((PyObject*)tempe);
        return NULL;
    }

    /* Convert base_lst to a true list. */

    if (!(base_lst = PySequence_Fast(base_lst, "argument must be an iterable"))) {
        return NULL;
    }

    /* Try to convert all items of base_lst to an mpz.
     * 
     * Note: MUST USE GMPy_MPZ_From_IntegerAndCopy()
     *       since the value is changed in-place.
     */

    seq_length = PySequence_Fast_GET_SIZE(base_lst);
    if (!(result = PyList_New(seq_length))) {
        Py_DECREF((PyObject*)tempe);
        Py_DECREF((PyObject*)tempm);
        Py_DECREF((PyObject*)base_lst);
        return NULL;
    }

    for (i=0; i < seq_length; i++) {
        if (!(tempres = GMPy_MPZ_From_IntegerAndCopy(PySequence_Fast_GET_ITEM(base_lst, i), NULL))) {
            Py_DECREF((PyObject*)tempe);
            Py_DECREF((PyObject*)tempm);
            Py_DECREF(base_lst);
            Py_DECREF(result);
            TYPE_ERROR("all items in iterable must be integers");
            return NULL;
        }

        if (PyList_SetItem(result, i, (PyObject*)tempres) < 0) {
            Py_DECREF((PyObject*)tempe);
            Py_DECREF((PyObject*)tempm);
            Py_DECREF(base_lst);
            Py_DECREF(result);
            return NULL;
        }
    }

    Py_BEGIN_ALLOW_THREADS;
    for (i=0; i < seq_length; i++) {
        temp = PySequence_Fast_GET_ITEM(result, i);
        mpz_powm(MPZ(temp), MPZ(temp), tempe->z, tempm->z);
    }
    Py_END_ALLOW_THREADS;

    Py_DECREF((PyObject*)tempe);
    Py_DECREF((PyObject*)tempm);
    Py_DECREF((PyObject*)base_lst);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_integer_powmod_base_list,
"powmod_base_list(base_lst, exp, mod, /) -> list[mpz, ...]\n\n"
"Returns list(powmod(i, exp, mod) for i in base_lst). Will always release\n"
"the GIL. (Experimental in gmpy2 2.1.x).");

static PyObject *
GMPy_Integer_PowMod_Base_List(PyObject *self, PyObject *args)
{
    if (PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("powmod_base_list requires 3 arguments");
        return NULL;
    }

    if (!PySequence_Check(PyTuple_GET_ITEM(args, 0))) {
        TYPE_ERROR("the first argument to powmod_base_list must be a sequence");
        return NULL;
    }

    int etype = GMPy_ObjectType(PyTuple_GET_ITEM(args, 1));
    int mtype = GMPy_ObjectType(PyTuple_GET_ITEM(args, 2));

    if (IS_TYPE_INTEGER(etype) && IS_TYPE_INTEGER(mtype))
        return GMPy_Integer_PowModBaseListWithType(PyTuple_GET_ITEM(args, 0),
                                                   PyTuple_GET_ITEM(args, 1), etype,
                                                   PyTuple_GET_ITEM(args, 2), mtype);

    TYPE_ERROR("powmod_base_list() requires integer arguments");
    return NULL;
}

static PyObject *
GMPy_Integer_PowModExpListWithType(PyObject *b, int btype,
                                   PyObject *exp_lst,
                                   PyObject *m, int mtype)
{
    MPZ_Object *tempb = NULL, *tempm = NULL, *tempres = NULL;
    PyObject *result = NULL, *temp = NULL;
    Py_ssize_t i, seq_length;

    if (!(tempm = GMPy_MPZ_From_IntegerWithType(m, mtype, NULL)) ||
        !(tempb = GMPy_MPZ_From_IntegerWithType(b, btype, NULL))) {
        return NULL;
    }

    if (mpz_sgn(tempm->z) < 1) {
        VALUE_ERROR("powmod_exp_list() 'mod' must be > 0");
        Py_DECREF((PyObject*)tempm);
        Py_DECREF((PyObject*)tempb);
        return NULL;
    }

    /* Convert exp_lst to a true list. */

    if (!(exp_lst = PySequence_Fast(exp_lst, "argument must be an iterable"))) {
        return NULL;
    }

    /* Try to convert all items of exp_lst to an mpz.
     *
     * Note: MUST USE GMPy_MPZ_From_IntegerAndCopy()
     *       since the value it changes in-place.
     */

    seq_length = PySequence_Fast_GET_SIZE(exp_lst);
    if (!(result = PyList_New(seq_length))) {
        Py_DECREF((PyObject*)tempb);
        Py_DECREF((PyObject*)tempm);
        Py_DECREF((PyObject*)exp_lst);
        return NULL;
    }

    for (i=0; i < seq_length; i++) {
        if (!(tempres = GMPy_MPZ_From_IntegerAndCopy(PySequence_Fast_GET_ITEM(exp_lst, i), NULL))) {
            Py_DECREF((PyObject*)tempb);
            Py_DECREF((PyObject*)tempm);
            Py_DECREF(exp_lst);
            Py_DECREF(result);
            TYPE_ERROR("all items in iterable must be integers");
            return NULL;
        }

        if (PyList_SetItem(result, i, (PyObject*)tempres) < 0) {
            Py_DECREF((PyObject*)tempb);
            Py_DECREF((PyObject*)tempm);
            Py_DECREF(exp_lst);
            Py_DECREF(result);
            return NULL;
        }
    }

    Py_BEGIN_ALLOW_THREADS;
    for (i=0; i < seq_length; i++) {
        temp = PySequence_Fast_GET_ITEM(result, i);
        mpz_powm(MPZ(temp), tempb->z, MPZ(temp), tempm->z);
    }
    Py_END_ALLOW_THREADS;

    Py_DECREF((PyObject*)tempb);
    Py_DECREF((PyObject*)tempm);
    Py_DECREF((PyObject*)exp_lst);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_integer_powmod_exp_list,
"powmod_exp_list(base, exp_lst, mod, /) -> list[mpz, ...]\n\n"
"Returns list(powmod(base, i, mod) for i in exp_lst). Will always release\n"
"the GIL. (Experimental in gmpy2 2.1.x).");

static PyObject *
GMPy_Integer_PowMod_Exp_List(PyObject *self, PyObject *args)
{
    if (PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("powmod_exp_list requires 3 arguments");
        return NULL;
    }

    if (!PySequence_Check(PyTuple_GET_ITEM(args, 1))) {
        TYPE_ERROR("the second argument to powmod_exp_list must be a sequence");
        return NULL;
    }

    int btype = GMPy_ObjectType(PyTuple_GET_ITEM(args, 0));
    int mtype = GMPy_ObjectType(PyTuple_GET_ITEM(args, 2));

    if (IS_TYPE_INTEGER(btype) && IS_TYPE_INTEGER(mtype))
        return GMPy_Integer_PowModExpListWithType(PyTuple_GET_ITEM(args, 0), btype,
                                                  PyTuple_GET_ITEM(args, 1),
                                                  PyTuple_GET_ITEM(args, 2), mtype);

    TYPE_ERROR("powmod_exp_list() requires integer arguments");
    return NULL;
}

static PyObject *
GMPy_Rational_PowWithType(PyObject *base, int btype, PyObject *exp, int etype,
                         PyObject *mod, CTXT_Object *context)
{
    MPQ_Object *tempbq = NULL, *resultq = NULL;
    MPZ_Object *tempez = NULL;
    int bsign;
    long tempexp;

    if (!Py_IsNone(mod)) {
        TYPE_ERROR("pow() 3rd argument not allowed unless all arguments are integers");
        return NULL;
    }

    /* Only support mpq**int. Everything else gets converted to mpf. */
    if (IS_TYPE_RATIONAL(btype) && IS_TYPE_INTEGER(etype)) {

        if (!(resultq = GMPy_MPQ_New(context)) ||
            !(tempbq = GMPy_MPQ_From_RationalWithType(base, btype, context)) ||
            !(tempez = GMPy_MPZ_From_IntegerWithType(exp, etype, context))) {
            Py_XDECREF((PyObject*)resultq);
            Py_XDECREF((PyObject*)tempbq);
            Py_XDECREF((PyObject*)tempez);
            return NULL;
        }

        if (!mpz_fits_slong_p(tempez->z)) {
            VALUE_ERROR("mpq.pow() outrageous exponent");
            Py_DECREF((PyObject*)resultq);
            Py_DECREF((PyObject*)tempbq);
            Py_DECREF((PyObject*)tempez);
            return NULL;
        }
        tempexp = (long) mpz_get_si(tempez->z);

        if (tempexp == 0) {
            mpq_set_si(resultq->q, 1, 1);
            Py_DECREF((PyObject*)tempbq);
            Py_DECREF((PyObject*)tempez);
            return (PyObject*)resultq;
        }

        bsign = mpq_sgn(tempbq->q);
        if (tempexp < 0) {
            if (bsign == 0) {
                ZERO_ERROR("pow() 0 base to negative exponent");
                Py_DECREF((PyObject*)resultq);
                Py_DECREF((PyObject*)tempbq);
                Py_DECREF((PyObject*)tempez);
                return NULL;
            }
            if (bsign < 0) {
                mpz_neg(mpq_numref(resultq->q), mpq_denref(tempbq->q));
            }
            else {
                mpz_set(mpq_numref(resultq->q), mpq_denref(tempbq->q));
            }
            mpz_abs(mpq_denref(resultq->q), mpq_numref(tempbq->q));
            tempexp = -tempexp;
        }
        else {
            mpq_set(resultq->q, tempbq->q);
        }

        if (tempexp > 1) {
            mpz_pow_ui(mpq_numref(resultq->q), mpq_numref(resultq->q), tempexp);
            mpz_pow_ui(mpq_denref(resultq->q), mpq_denref(resultq->q), tempexp);
        }
        Py_DECREF((PyObject*)tempbq);
        Py_DECREF((PyObject*)tempez);
        return (PyObject*)resultq;
    }
    else {
        return GMPy_Real_PowWithType(base, btype, exp, etype, Py_None, context);
    }
}

static PyObject *
GMPy_Real_PowWithType(PyObject *base, int btype, PyObject *exp, int etype, 
                      PyObject *mod, CTXT_Object *context)
{
    MPFR_Object *tempb = NULL, *tempe = NULL, *result = NULL;
    MPZ_Object *tempz = NULL;
    MPC_Object *mpc_result = NULL;

    if (!Py_IsNone(mod)) {
        TYPE_ERROR("pow() 3rd argument not allowed unless all arguments are integers");
        return NULL;
    }

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPFR_New(0, context)) ||
        !(tempb = GMPy_MPFR_From_RealWithType(base, btype, 1, context))) {
        goto err;
    }

    mpfr_clear_flags();

    if (IS_TYPE_PyInteger(etype)) {
        int error;
        unsigned long intb;
        long temp;

        if (mpfr_fits_ulong_p(tempb->f, MPFR_RNDF)) {
            /* Need to check the inexact flag to verify that tempb is an integer. */
            intb = mpfr_get_ui(tempb->f, MPFR_RNDF);
            if (mpfr_inexflag_p()) {
                mpfr_clear_inexflag();
            }
            else {
                temp = PyLong_AsLongAndOverflow(exp, &error);
                if (!error) {
                    if (temp >= 0) {
                        result->rc = mpfr_ui_pow_ui(result->f, intb, temp, GET_MPFR_ROUND(context));
                        goto done;
                    }
                }
            }
        }
    }

    if (IS_TYPE_INTEGER(etype)) {
        if (!(tempz = GMPy_MPZ_From_IntegerWithType(exp, etype, context))) {
            goto err;
        }
        result->rc = mpfr_pow_z(result->f, tempb->f, tempz->z, GET_MPFR_ROUND(context));
        goto done;
    }
    else if (IS_TYPE_REAL(etype)) {
        if (!(tempe = GMPy_MPFR_From_RealWithType(exp, etype, 1, context))) {
            goto err;
        }
        result->rc = mpfr_pow(result->f, tempb->f, tempe->f, GET_MPFR_ROUND(context));
        goto done;
    }

    /* If the result is NaN, check if a complex result works. */
    if (mpfr_nanflag_p() && context->ctx.allow_complex) {
        mpc_result = (MPC_Object*)GMPy_Complex_PowWithType(base, btype, exp, etype, Py_None, context);
        if (!mpc_result || MPC_IS_NAN_P(mpc_result)) {
            Py_XDECREF((PyObject*)mpc_result);
            context->ctx.invalid = 1;
            GMPY_INVALID("pow() invalid operation");
            goto err;
        }
        /* return a valid complex result */
        Py_XDECREF((PyObject*)tempe);
        Py_XDECREF((PyObject*)tempz);
        Py_XDECREF((PyObject*)tempb);
        Py_XDECREF((PyObject*)result);
        return (PyObject*)mpc_result;
    }

  done:
    _GMPy_MPFR_Cleanup(&result, context);
    Py_XDECREF((PyObject*)tempz);
    Py_XDECREF((PyObject*)tempe);
    Py_XDECREF((PyObject*)tempb);
    return (PyObject*)result;

  err:
    Py_XDECREF((PyObject*)result);
    Py_XDECREF((PyObject*)tempz);
    Py_XDECREF((PyObject*)tempe);
    Py_XDECREF((PyObject*)tempb);
    return NULL;
}

static PyObject *
GMPy_Complex_PowWithType(PyObject *base, int btype, PyObject *exp, int etype,
                         PyObject *mod, CTXT_Object *context)
{
    MPC_Object *tempb = NULL, *tempe = NULL, *result= NULL;
    MPFR_Object *tempf = NULL;
    MPZ_Object *tempz = NULL;

    if (!Py_IsNone(mod)) {
        TYPE_ERROR("pow() 3rd argument not allowed unless all arguments are integers");
        return NULL;
    }

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPC_New(0, 0, context)) ||
        !(tempb = GMPy_MPC_From_ComplexWithType(base, btype, 1, 1, context))) {
        goto err;
    }

    mpfr_clear_flags();

    if (IS_TYPE_INTEGER(etype)) {
        if (!(tempz = GMPy_MPZ_From_IntegerWithType(exp, etype, context))) {
            goto err;
        }
        result->rc = mpc_pow_z(result->c, tempb->c, tempz->z, GET_MPC_ROUND(context));
        goto done;
    }

    if (IS_TYPE_REAL(etype)) {
        if (!(tempf = GMPy_MPFR_From_RealWithType(exp, etype, 1, context))) {
            goto err;
        }

        result->rc = mpc_pow_fr(result->c, tempb->c, tempf->f, GET_MPC_ROUND(context));
        goto done;
    }

    if (IS_TYPE_COMPLEX(etype)) {
        if (!(tempe = GMPy_MPC_From_ComplexWithType(exp, etype, 1, 1, context))) {
            goto err;
        }

        result->rc = mpc_pow(result->c, tempb->c, tempe->c, GET_MPC_ROUND(context));
        goto done;
    }

    TYPE_ERROR("pow() argument types not supported");
    goto err;

  done:
    _GMPy_MPC_Cleanup(&result, context);
    Py_XDECREF((PyObject*)tempz);
    Py_XDECREF((PyObject*)tempf);
    Py_XDECREF((PyObject*)tempe);
    Py_XDECREF((PyObject*)tempb);
    return (PyObject*)result;

  err:
    Py_XDECREF((PyObject*)result);
    Py_XDECREF((PyObject*)tempz);
    Py_XDECREF((PyObject*)tempf);
    Py_XDECREF((PyObject*)tempe);
    Py_XDECREF((PyObject*)tempb);
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_integer_powmod,
"powmod(x, y, m, /) -> mpz\n\n"
"Return (x**y) mod m. Same as the three argument version of Python's\n"
"built-in `pow`, but converts all three arguments to `mpz`.");

static PyObject *
GMPy_Integer_PowMod(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *m;
    int xtype, ytype, mtype;

    if (PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("powmod() requires 3 arguments.");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    m = PyTuple_GET_ITEM(args, 2);

    xtype = GMPy_ObjectType(x);
    ytype = GMPy_ObjectType(y);
    mtype = GMPy_ObjectType(m);

    if (IS_TYPE_INTEGER(xtype) &&
        IS_TYPE_INTEGER(ytype) &&
        IS_TYPE_INTEGER(mtype)) {
        return GMPy_Integer_PowWithType(x, xtype, y, ytype, m, NULL);
    }

    TYPE_ERROR("powmod() argument types not supported");
    return NULL;
}


PyDoc_STRVAR(GMPy_doc_integer_powmod_sec,
"powmod_sec(x, y, m, /) -> mpz\n\n"
"Return (x**y) mod m. Calculates x ** y (mod m) but using a constant\n"
"time algorithm to reduce the risk of side channel attacks. y must be\n"
"an integer >0. m must be an odd integer.");

static PyObject *
GMPy_Integer_PowMod_Sec(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *m;
    int xtype, ytype, mtype;
    MPZ_Object *tempx = NULL, *tempy = NULL, *tempm = NULL, *result = NULL;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("powmod_sec() requires 3 arguments.");
        goto err;
    }

    if (!(result = GMPy_MPZ_New(NULL))) {
        /* LCOV_EXCL_START */
        goto err;
        /* LCOV_EXCL_STOP */
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    m = PyTuple_GET_ITEM(args, 2);

    xtype = GMPy_ObjectType(x);
    ytype = GMPy_ObjectType(y);
    mtype = GMPy_ObjectType(m);

    /* Validate base. */

    if (!IS_TYPE_INTEGER(xtype)) {
        TYPE_ERROR("powmod_sec() base must be an integer.");
        goto err;
    }

    if (!(tempx = GMPy_MPZ_From_IntegerWithType(x, xtype, NULL))) {
        /* LCOV_EXCL_START */
        goto err;
        /* LCOV_EXCL_STOP */
    }

    /* Validate exponent. It must be > 0. */

    if (!IS_TYPE_INTEGER(ytype)) {
        TYPE_ERROR("powmod_sec() exponent must be an integer.");
        goto err;
    }

    if (!(tempy = GMPy_MPZ_From_IntegerWithType(y, ytype, NULL))) {
        /* LCOV_EXCL_START */
        goto err;
        /* LCOV_EXCL_STOP */
    }

    if (!(mpz_sgn(tempy->z) == 1)) {
        VALUE_ERROR("powmod_sec() exponent must be > 0.");
        goto err;
    }
    /* Validate modulus. It must be odd.*/

    if (!IS_TYPE_INTEGER(mtype)) {
        TYPE_ERROR("powmod_sec() modulus must be an integer.");
        goto err;
    }

    if (!(tempm = GMPy_MPZ_From_IntegerWithType(m, mtype, NULL))) {
        /* LCOV_EXCL_START */
        goto err;
        /* LCOV_EXCL_STOP */
    }

    if (mpz_even_p(tempm->z)) {
        VALUE_ERROR("powmod_sec() modulus must be odd.");
        goto err;
    }

    GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
    mpz_powm_sec(result->z, tempx->z, tempy->z, tempm->z);
    GMPY_MAYBE_END_ALLOW_THREADS(context);

    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempm);
    return (PyObject*)result;

  err:
    Py_XDECREF(tempx);
    Py_XDECREF(tempy);
    Py_XDECREF(tempm);
    return NULL;
}

static PyObject *
GMPy_Number_Pow(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    int xtype = GMPy_ObjectType(x);
    int ytype = GMPy_ObjectType(y);

    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype))
        return GMPy_Integer_PowWithType(x, xtype, y, ytype, z, context);

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype))
        return GMPy_Rational_PowWithType(x, xtype, y, ytype, z, context);

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype))
        return GMPy_Real_PowWithType(x, xtype, y, ytype, z, context);

    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype))
        return GMPy_Complex_PowWithType(x, xtype, y, ytype, z, context);

    TYPE_ERROR("pow() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_context_pow,
"context.pow(x, y, /) -> mpz | mpq | mpfr | mpc\n\n"
"Return x ** y.");

static PyObject *
GMPy_Context_Pow(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("pow() requires 2 arguments.");
        return NULL;
    }
    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Pow(PyTuple_GET_ITEM(args, 0),
                           PyTuple_GET_ITEM(args, 1),
                           Py_None, context);
}

static PyObject *
GMPy_Number_Pow_Slot(PyObject *base, PyObject *exp, PyObject *mod)
{
    int btype = GMPy_ObjectType(base);
    int etype = GMPy_ObjectType(exp);

    if (IS_TYPE_INTEGER(btype) && IS_TYPE_INTEGER(etype))
        return GMPy_Integer_PowWithType(base, btype, exp, etype, mod, NULL);

    if (IS_TYPE_RATIONAL(btype) && IS_TYPE_RATIONAL(etype))
        return GMPy_Rational_PowWithType(base, btype, exp, etype, mod, NULL);

    if (IS_TYPE_REAL(btype) && IS_TYPE_REAL(etype))
        return GMPy_Real_PowWithType(base, btype, exp, etype, mod, NULL);

    if (IS_TYPE_COMPLEX(btype) && IS_TYPE_COMPLEX(etype))
        return GMPy_Complex_PowWithType(base, btype, exp, etype, mod, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}
