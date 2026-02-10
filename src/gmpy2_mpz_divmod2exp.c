/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpz_divmod2exp.c                                                  *
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


/* This file contains functions related to division and remainder by a power
 * of two.
 */

/*
 **************************************************************************
 * Ceiling division and remainder by power of two.
 **************************************************************************
 */

PyDoc_STRVAR(doc_c_divmod_2exp,
"c_divmod_2exp($module, x, n, /)\n--\n\n"
"Return the quotient and remainder of x divided by 2**n. The quotient\n"
"is rounded towards +Inf (ceiling rounding) and the remainder will\n"
"be negative. x must be an integer. n must be >0.");

static PyObject *
GMPy_MPZ_c_divmod_2exp(PyObject *self, PyObject *args)
{
    mp_bitcnt_t nbits;
    PyObject *result = NULL;
    MPZ_Object *q = NULL, *r = NULL, *tempx = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("c_divmod_2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = GMPy_Integer_AsMpBitCnt(PyTuple_GET_ITEM(args, 1));
    if (nbits == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL)) ||
        !(q = GMPy_MPZ_New(NULL)) ||
        !(r = GMPy_MPZ_New(NULL)) ||
        !(result = PyTuple_New(2))) {

        Py_XDECREF(result);
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)q);
        Py_XDECREF((PyObject*)r);
        return NULL;
    }

    mpz_cdiv_q_2exp(q->z, tempx->z, nbits);
    mpz_cdiv_r_2exp(r->z, tempx->z, nbits);

    Py_DECREF((PyObject*)tempx);
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

PyDoc_STRVAR(doc_c_div_2exp,
"c_div_2exp($module, x, n, /)\n--\n\n"
"Returns the quotient of x divided by 2**n. The quotient is rounded\n"
"towards +Inf (ceiling rounding). x must be an integer. n must be >0.");

static PyObject *
GMPy_MPZ_c_div_2exp(PyObject *self, PyObject *args)
{
    mp_bitcnt_t nbits;
    MPZ_Object *result, *tempx;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("c_div_2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = GMPy_Integer_AsMpBitCnt(PyTuple_GET_ITEM(args, 1));
    if (nbits == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    result = GMPy_MPZ_New(NULL);
    if (!tempx || !result) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        return NULL;
    }

    mpz_cdiv_q_2exp(result->z, tempx->z, nbits);
    Py_DECREF((PyObject*)tempx);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_c_mod_2exp,
"c_mod_2exp($module, x, n, /)\n--\n\n"
"Return the remainder of x divided by 2**n. The remainder will be\n"
"negative. x must be an integer. n must be >0.");

static PyObject *
GMPy_MPZ_c_mod_2exp(PyObject *self, PyObject *args)
{
    mp_bitcnt_t nbits;
    MPZ_Object *result, *tempx;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("c_mod_2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = GMPy_Integer_AsMpBitCnt(PyTuple_GET_ITEM(args, 1));
    if (nbits == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    result = GMPy_MPZ_New(NULL);
    if (!tempx || !result) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        return NULL;
    }

    mpz_cdiv_r_2exp(result->z, tempx->z, nbits);
    Py_DECREF((PyObject*)tempx);
    return (PyObject*)result;
}

/*
 **************************************************************************
 * Floor division and remainder by power of two.
 **************************************************************************
 */

PyDoc_STRVAR(doc_f_divmod_2exp,
"f_divmod_2exp($module, x, n, /)\n--\n\n"
"Return quotient and remainder after dividing x by 2**n. The quotient\n"
"is rounded towards -Inf (floor rounding) and the remainder will be\n"
"positive. x must be an integer. n must be >0.");

static PyObject *
GMPy_MPZ_f_divmod_2exp(PyObject *self, PyObject *args)
{
    mp_bitcnt_t nbits;
    PyObject *result;
    MPZ_Object *q, *r, *tempx;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("f_divmod_2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = GMPy_Integer_AsMpBitCnt(PyTuple_GET_ITEM(args, 1));
    if (nbits == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    q = GMPy_MPZ_New(NULL);
    r = GMPy_MPZ_New(NULL);
    result = PyTuple_New(2);
    if (!tempx || !q || !r || !result) {
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)q);
        Py_XDECREF((PyObject*)r);
        return NULL;
    }

    mpz_fdiv_q_2exp(q->z, tempx->z, nbits);
    mpz_fdiv_r_2exp(r->z, tempx->z, nbits);

    Py_DECREF((PyObject*)tempx);
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

PyDoc_STRVAR(doc_f_div_2exp,
"f_div_2exp($module, x, n, /)\n--\n\n"
"Return the quotient of x divided by 2**n. The quotient is rounded\n"
"towards -Inf (floor rounding). x must be an integer. n must be >0.");

static PyObject *
GMPy_MPZ_f_div_2exp(PyObject *self, PyObject *args)
{
    mp_bitcnt_t nbits;
    MPZ_Object *result, *tempx;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("f_div_2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = GMPy_Integer_AsMpBitCnt(PyTuple_GET_ITEM(args, 1));
    if (nbits == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    result = GMPy_MPZ_New(NULL);
    if (!tempx || !result) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        return NULL;
    }

    mpz_fdiv_q_2exp(result->z, tempx->z, nbits);
    Py_DECREF((PyObject*)tempx);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_f_mod_2exp,
"f_mod_2exp($module, x, n, /)\n--\n\n"
"Return remainder of x divided by 2**n. The remainder will be\n"
"positive. x must be an integer. n must be >0.");

static PyObject *
GMPy_MPZ_f_mod_2exp(PyObject *self, PyObject *args)
{
    mp_bitcnt_t nbits;
    MPZ_Object *result, *tempx;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("f_mod_2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = GMPy_Integer_AsMpBitCnt(PyTuple_GET_ITEM(args, 1));
    if (nbits == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    result = GMPy_MPZ_New(NULL);
    if (!tempx || !result) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        return NULL;
    }

    mpz_fdiv_r_2exp(result->z, tempx->z, nbits);
    Py_DECREF((PyObject*)tempx);
    return (PyObject*)result;
}

/*
 **************************************************************************
 * Truncating division and remainder by power of two.
 **************************************************************************
 */

PyDoc_STRVAR(doc_t_divmod_2exp,
"t_divmod_2exp($module, x, n, /)\n--\n\n"
"Return the quotient and remainder of x divided by 2**n. The quotient\n"
"is rounded towards zero (truncation) and the remainder will have the\n"
"same sign as x. x must be an integer. n must be >0.");

static PyObject *
GMPy_MPZ_t_divmod_2exp(PyObject *self, PyObject *args)
{
    mp_bitcnt_t nbits;
    MPZ_Object *q, *r, *tempx;
    PyObject *result;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("t_divmod_2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = GMPy_Integer_AsMpBitCnt(PyTuple_GET_ITEM(args, 1));
    if (nbits == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    q = GMPy_MPZ_New(NULL);
    r = GMPy_MPZ_New(NULL);
    result = PyTuple_New(2);
    if (!tempx || !q || !r || !result) {
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)q);
        Py_XDECREF((PyObject*)r);
        return NULL;
    }

    mpz_tdiv_q_2exp(q->z, tempx->z, nbits);
    mpz_tdiv_r_2exp(r->z, tempx->z, nbits);

    Py_DECREF((PyObject*)tempx);
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

PyDoc_STRVAR(doc_t_div_2exp,
"t_div_2exp($module, x, n, /)\n--\n\n"
"Return the quotient of x divided by 2**n. The quotient is rounded\n"
"towards zero (truncation). n must be >0.");

static PyObject *
GMPy_MPZ_t_div_2exp(PyObject *self, PyObject *args)
{
    mp_bitcnt_t nbits;
    MPZ_Object *result, *tempx;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("t_div_2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = GMPy_Integer_AsMpBitCnt(PyTuple_GET_ITEM(args, 1));
    if (nbits == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    result = GMPy_MPZ_New(NULL);
    if (!tempx || !result) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        return NULL;
    }

    mpz_tdiv_q_2exp(result->z, tempx->z, nbits);
    Py_DECREF((PyObject*)tempx);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_t_mod_2exp,
"t_mod_2exp($module, x, n, /)\n--\n\n"
"Return the remainder of x divided by 2**n. The remainder will have\n"
"the same sign as x. x must be an integer. n must be >0.");

static PyObject *
GMPy_MPZ_t_mod_2exp(PyObject *self, PyObject *args)
{
    mp_bitcnt_t nbits;
    MPZ_Object *result, *tempx;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("t_mod_2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = GMPy_Integer_AsMpBitCnt(PyTuple_GET_ITEM(args, 1));
    if (nbits == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    result = GMPy_MPZ_New(NULL);
    if (!tempx || !result) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        return NULL;
    }

    mpz_tdiv_r_2exp(result->z, tempx->z, nbits);
    Py_DECREF((PyObject*)tempx);
    return (PyObject*)result;
}
