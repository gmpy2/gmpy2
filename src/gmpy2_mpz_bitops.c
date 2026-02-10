/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpz_bitops.c                                                      *
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

PyDoc_STRVAR(doc_bit_length_method,
"bit_length($self)\n--\n\n"
"Return the number of significant bits in the radix-2\n"
"representation of self. Note: mpz(0).bit_length() returns 0.");

static PyObject *
GMPy_MPZ_bit_length_method(PyObject *self, PyObject *other)
{
    mp_bitcnt_t n = 0;

    if (mpz_size(MPZ(self)))
        n = mpz_sizeinbase(MPZ(self), 2);

    return GMPy_PyLong_FromMpBitCnt(n);
}

PyDoc_STRVAR(doc_bit_length_function,
"bit_length($module, x, /)\n--\n\n"
"Return the number of significant bits in the radix-2\n"
"representation of x. Note: bit_length(0) returns 0.");

static PyObject *
GMPy_MPZ_bit_length_function(PyObject *self, PyObject *other)
{
    mp_bitcnt_t n = 0;
    MPZ_Object* tempx;

    if (!(tempx = GMPy_MPZ_From_Integer(other, NULL))) {
        TYPE_ERROR("bit_length() requires 'mpz' argument");
        return NULL;
    }
    if (mpz_size(MPZ(tempx)))
        n = mpz_sizeinbase(tempx->z, 2);

    Py_DECREF((PyObject*)tempx);
    return GMPy_PyLong_FromMpBitCnt(n);
}

PyDoc_STRVAR(doc_bit_mask,
"bit_mask($module, n, /)\n--\n\n"
"Return an `mpz` exactly n bits in length with all bits set.\n");

static PyObject *
GMPy_MPZ_bit_mask(PyObject *self, PyObject *other)
{
    mp_bitcnt_t n = 0;
    MPZ_Object* result;

    n = GMPy_Integer_AsMpBitCnt(other);
    if (n == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(NULL)))
        return NULL;

    mpz_set_ui(result->z, 1);
    mpz_mul_2exp(result->z, result->z, n);
    mpz_sub_ui(result->z, result->z, 1);

    return (PyObject*)result;
}

/* return scan0/scan1 for an mpz */
PyDoc_STRVAR(doc_bit_scan0_method,
"bit_scan0($self, n=0, /)\n--\n\n"
"Return the index of the first 0-bit of self with index >= n. n >= 0.\n"
"If there are no more 0-bits in self at or above index n (which can\n"
"only happen for negatives, assuming an infinitely long 2's complement\n"
"format), then `None` is returned.");

static PyObject *
GMPy_MPZ_bit_scan0_method(PyObject *self, PyObject *const *args,
                          Py_ssize_t nargs)
{
    mp_bitcnt_t index, starting_bit = 0;

    if (nargs == 1) {
        starting_bit = GMPy_Integer_AsMpBitCnt(args[0]);
        if (starting_bit == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
            return NULL;
        }
    }

    index = mpz_scan0(MPZ(self), starting_bit);

    if (index == (mp_bitcnt_t)(-1)) {
        Py_RETURN_NONE;
    }
    else {
        return GMPy_PyLong_FromMpBitCnt(index);
    }
}

PyDoc_STRVAR(doc_bit_scan0_function,
"bit_scan0($module, x, n=0, /)\n--\n\n"
"Return the index of the first 0-bit of x with index >= n. n >= 0.\n"
"If there are no more 0-bits in x at or above index n (which can\n"
"only happen for negatives, assuming an infinitely long 2's complement\n"
"format), then `None` is returned.");

static PyObject *
GMPy_MPZ_bit_scan0_function(PyObject *self, PyObject *const *args,
                            Py_ssize_t nargs)
{
    mp_bitcnt_t index, starting_bit = 0;
    MPZ_Object *tempx = NULL;

    if (nargs == 0 || nargs > 2) {
        goto err;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(args[0], NULL))) {
        goto err;
    }

    if (nargs == 2) {
        starting_bit = GMPy_Integer_AsMpBitCnt(args[1]);
        if (starting_bit == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
            goto err_index;
        }
    }

    index = mpz_scan0(tempx->z, starting_bit);

    Py_DECREF((PyObject*)tempx);
    if (index == (mp_bitcnt_t)(-1)) {
        Py_RETURN_NONE;
    }
    else {
        return GMPy_PyLong_FromMpBitCnt(index);
    }

  err:
    TYPE_ERROR("bit_scan0() requires 'mpz',['int'] arguments");
  err_index:
    Py_XDECREF((PyObject*)tempx);
    return NULL;
}

PyDoc_STRVAR(doc_bit_scan1_method,
"bit_scan1($self, n=0, /)\n--\n\n"
"Return the index of the first 1-bit of self with index >= n. n >= 0.\n"
"If there are no more 1-bits in self at or above index n (which can\n"
"only happen for non-negatives, assuming an infinitely long 2's complement\n"
"format), then `None` is returned.");

static PyObject *
GMPy_MPZ_bit_scan1_method(PyObject *self, PyObject *const *args,
                          Py_ssize_t nargs)
{
    mp_bitcnt_t index, starting_bit = 0;

    if (nargs == 1) {
        starting_bit = GMPy_Integer_AsMpBitCnt(args[0]);
        if (starting_bit == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
            return NULL;
        }
    }

    index = mpz_scan1(MPZ(self), starting_bit);

    if (index == (mp_bitcnt_t)(-1)) {
        Py_RETURN_NONE;
    }
    else {
        return GMPy_PyLong_FromMpBitCnt(index);
    }
}

PyDoc_STRVAR(doc_bit_scan1_function,
"bit_scan1($module, x, n=0, /)\n--\n\n"
"Return the index of the first 1-bit of x with index >= n. n >= 0.\n"
"If there are no more 1-bits in x at or above index n (which can\n"
"only happen for x>=0, assuming an infinitely long 2's complement\n"
"format), then `None` is returned.");

static PyObject *
GMPy_MPZ_bit_scan1_function(PyObject *self, PyObject *const *args,
                            Py_ssize_t nargs)
{
    mp_bitcnt_t index, starting_bit = 0;
    MPZ_Object *tempx = NULL;

    if (nargs == 0 || nargs > 2) {
        goto err;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(args[0], NULL))) {
        goto err;
    }

    if (nargs == 2) {
        starting_bit = GMPy_Integer_AsMpBitCnt(args[1]);
        if (starting_bit == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
            goto err_index;
        }
    }

    index = mpz_scan1(tempx->z, starting_bit);

    Py_DECREF((PyObject*)tempx);
    if (index == (mp_bitcnt_t)(-1)) {
        Py_RETURN_NONE;
    }
    else {
        return GMPy_PyLong_FromMpBitCnt(index);
    }

  err:
    TYPE_ERROR("bit_scan1() requires 'mpz',['int'] arguments");
  err_index:
    Py_XDECREF((PyObject*)tempx);
    return NULL;
}

/* get & return one bit from an mpz */
PyDoc_STRVAR(doc_bit_test_function,
"bit_test($module, x, n, /)\n--\n\n"
"Return the value of the n-th bit of x.");

static PyObject *
GMPy_MPZ_bit_test_function(PyObject *self, PyObject *const *args,
                           Py_ssize_t nargs)
{
    mp_bitcnt_t bit_index;
    int temp;
    MPZ_Object *tempx = NULL;

    if (nargs != 2) {
        goto err;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(args[0], NULL))) {
        goto err;
    }

    bit_index = GMPy_Integer_AsMpBitCnt(args[1]);
    if (bit_index == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        goto err_index;
    }

    temp = mpz_tstbit(tempx->z, bit_index);
    Py_DECREF((PyObject*)tempx);

    if (temp)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;

  err:
    TYPE_ERROR("bit_test() requires 'mpz','int' arguments");
  err_index:
    Py_XDECREF((PyObject*)tempx);
    return NULL;
}

PyDoc_STRVAR(doc_bit_test_method,
"bit_test(n, /)\n--\n\n"
"Return the value of the n-th bit of self.");

static PyObject *
GMPy_MPZ_bit_test_method(PyObject *self, PyObject *other)
{
    mp_bitcnt_t bit_index;

    bit_index = GMPy_Integer_AsMpBitCnt(other);
    if (bit_index == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if (mpz_tstbit(MPZ(self), bit_index))
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(doc_bit_clear_function,
"bit_clear($module, x, n, /)\n--\n\n"
"Return a copy of x with the n-th bit cleared.");

static PyObject *
GMPy_MPZ_bit_clear_function(PyObject *self, PyObject *args)
{
    mp_bitcnt_t bit_index;
    MPZ_Object *result = NULL, *tempx = NULL;

    if (PyTuple_GET_SIZE(args) != 2)
        goto err;

    if (!(result = GMPy_MPZ_New(NULL)))
        return NULL;

    if (!(tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL)))
        goto err;

    bit_index = GMPy_Integer_AsMpBitCnt(PyTuple_GET_ITEM(args, 1));
    if (bit_index == (mp_bitcnt_t)(-1) && PyErr_Occurred())
        goto err_index;

    mpz_set(result->z, tempx->z);
    mpz_clrbit(result->z, bit_index);

    Py_DECREF((PyObject*)tempx);
    return (PyObject*)result;

  err:
    TYPE_ERROR("bit_clear() requires 'mpz','int' arguments");
  err_index:
    Py_XDECREF((PyObject*)result);
    Py_XDECREF((PyObject*)tempx);
    return NULL;
}

PyDoc_STRVAR(doc_bit_clear_method,
"bit_clear($self, n, /)\n--\n\n"
"Return a copy of self with the n-th bit cleared.");

static PyObject *
GMPy_MPZ_bit_clear_method(PyObject *self, PyObject *other)
{
    mp_bitcnt_t bit_index;
    MPZ_Object *result = NULL;

    if (!(result = GMPy_MPZ_New(NULL)))
        return NULL;

    bit_index = GMPy_Integer_AsMpBitCnt(other);
    if (bit_index == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        Py_DECREF(result);
        return NULL;
    }

    mpz_set(result->z, MPZ(self));
    mpz_clrbit(result->z, bit_index);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_bit_set_function,
"bit_set($module, x, n, /)\n--\n\n"
"Return a copy of x with the n-th bit set.");

static PyObject *
GMPy_MPZ_bit_set_function(PyObject *self, PyObject *args)
{
    mp_bitcnt_t bit_index;
    MPZ_Object *result = NULL, *tempx = NULL;

    if (PyTuple_GET_SIZE(args) != 2)
        goto err;

    if (!(result = GMPy_MPZ_New(NULL)))
        return NULL;

    if (!(tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL)))
        goto err;

    bit_index = GMPy_Integer_AsMpBitCnt(PyTuple_GET_ITEM(args, 1));
    if (bit_index == (mp_bitcnt_t)(-1) && PyErr_Occurred())
        goto err_index;

    mpz_set(result->z, tempx->z);
    mpz_setbit(result->z, bit_index);

    Py_DECREF((PyObject*)tempx);
    return (PyObject*)result;

  err:
    TYPE_ERROR("bit_set() requires 'mpz','int' arguments");
  err_index:
    Py_XDECREF((PyObject*)result);
    Py_XDECREF((PyObject*)tempx);
    return NULL;
}

PyDoc_STRVAR(doc_bit_set_method,
"bit_set($self, n, /)\n--\n\n"
"Return a copy of self with the n-th bit set.");

static PyObject *
GMPy_MPZ_bit_set_method(PyObject *self, PyObject *other)
{
    mp_bitcnt_t bit_index;
    MPZ_Object *result = NULL;

    if (!(result = GMPy_MPZ_New(NULL)))
        return NULL;

    bit_index = GMPy_Integer_AsMpBitCnt(other);
    if (bit_index == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        Py_DECREF(result);
        return NULL;
    }

    mpz_set(result->z, MPZ(self));
    mpz_setbit(result->z, bit_index);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_bit_flip_function,
"bit_flip($module, x, n, /)\n--\n\n"
"Return a copy of x with the n-th bit inverted.");

static PyObject *
GMPy_MPZ_bit_flip_function(PyObject *self, PyObject *args)
{
    mp_bitcnt_t bit_index;
    MPZ_Object *result = NULL, *tempx = NULL;

    if (PyTuple_GET_SIZE(args) != 2)
        goto err;

    if (!(result = GMPy_MPZ_New(NULL)))
        return NULL;

    if (!(tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL)))
        goto err;

    bit_index = GMPy_Integer_AsMpBitCnt(PyTuple_GET_ITEM(args, 1));
    if (bit_index == (mp_bitcnt_t)(-1) && PyErr_Occurred())
        goto err_index;

    mpz_set(result->z, tempx->z);
    mpz_combit(result->z, bit_index);

    Py_DECREF((PyObject*)tempx);
    return (PyObject*)result;

  err:
    TYPE_ERROR("bit_flip() requires 'mpz','int' arguments");
  err_index:
    Py_XDECREF((PyObject*)result);
    Py_XDECREF((PyObject*)tempx);
    return NULL;
}

PyDoc_STRVAR(doc_bit_flip_method,
"bit_flip($self, n, /)\n--\n\n"
"Return a copy of self with the n-th bit inverted.");

static PyObject *
GMPy_MPZ_bit_flip_method(PyObject *self, PyObject *other)
{
    mp_bitcnt_t bit_index;
    MPZ_Object *result = NULL;

    if (!(result = GMPy_MPZ_New(NULL)))
        return NULL;

    bit_index = GMPy_Integer_AsMpBitCnt(other);
    if (bit_index == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        Py_DECREF(result);
        return NULL;
    }

    mpz_set(result->z, MPZ(self));
    mpz_combit(result->z, bit_index);
    return (PyObject*)result;
}

static PyObject *
GMPy_MPZ_Invert_Slot(MPZ_Object *self)
{
    MPZ_Object *result;

    if ((result = GMPy_MPZ_New(NULL)))
        mpz_com(result->z, MPZ(self));

    return (PyObject*)result;
}

static PyObject *
GMPy_MPZ_And_Slot(PyObject *self, PyObject *other)
{
    MPZ_Object *result;

    if (CHECK_MPZANY(self)) {
        if (CHECK_MPZANY(other)) {
            if (!(result = GMPy_MPZ_New(NULL)))
                return NULL;
            mpz_and(result->z, MPZ(self), MPZ(other));
        }
        else {
            if (!(result = GMPy_MPZ_From_Integer(other, NULL)))
                return NULL;
            mpz_and(result->z, MPZ(self), result->z);
        }
    }
    else if (CHECK_MPZANY(other)) {
        if (!(result = GMPy_MPZ_From_Integer(self, NULL)))
            return NULL;
        mpz_and(result->z, result->z, MPZ(other));
    }
    else {
        /* LCOV_EXCL_START */
        Py_RETURN_NOTIMPLEMENTED;
        /* LCOV_EXCL_STOP */
    }
    return (PyObject*)result;
}

static PyObject *
GMPy_MPZ_Ior_Slot(PyObject *self, PyObject *other)
{
    MPZ_Object *result;

    if (CHECK_MPZANY(self)) {
        if (CHECK_MPZANY(other)) {
            if (!(result = GMPy_MPZ_New(NULL)))
                return NULL;
            mpz_ior(result->z, MPZ(self), MPZ(other));
        }
        else {
            if (!(result = GMPy_MPZ_From_Integer(other, NULL)))
                return NULL;
            mpz_ior(result->z, MPZ(self), result->z);
        }
    }
    else if (CHECK_MPZANY(other)) {
        if (!(result = GMPy_MPZ_From_Integer(self, NULL)))
            return NULL;
        mpz_ior(result->z, result->z, MPZ(other));
    }
    else {
        /* LCOV_EXCL_START */
        Py_RETURN_NOTIMPLEMENTED;
        /* LCOV_EXCL_STOP */
    }
    return (PyObject*)result;
}

static PyObject *
GMPy_MPZ_Xor_Slot(PyObject *self, PyObject *other)
{
    MPZ_Object *result;

    if (CHECK_MPZANY(self)) {
        if (CHECK_MPZANY(other)) {
            if (!(result = GMPy_MPZ_New(NULL)))
                return NULL;
            mpz_xor(result->z, MPZ(self), MPZ(other));
        }
        else {
            if (!(result = GMPy_MPZ_From_Integer(other, NULL)))
                return NULL;
            mpz_xor(result->z, MPZ(self), result->z);
        }
    }
    else if (CHECK_MPZANY(other)) {
        if (!(result = GMPy_MPZ_From_Integer(self, NULL)))
            return NULL;
        mpz_xor(result->z, result->z, MPZ(other));
    }
    else {
        /* LCOV_EXCL_START */
        Py_RETURN_NOTIMPLEMENTED;
        /* LCOV_EXCL_STOP */
    }
    return (PyObject*)result;
}

static PyObject *
GMPy_MPZ_Rshift_Slot(PyObject *self, PyObject *other)
{
    mp_bitcnt_t count;
    MPZ_Object *result, *tempx;

    count = GMPy_Integer_AsMpBitCnt(other);
    if ((count == (mp_bitcnt_t)(-1)) && PyErr_Occurred()) {
        if (!PyErr_ExceptionMatches(PyExc_OverflowError)) {
            /* LCOV_EXCL_START */
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        PyErr_Clear();

        PyObject *tmp = PyNumber_Long(other);

        if (!tmp) {
            /* LCOV_EXCL_START */
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        if (PyLong_IsNegative(tmp)) {
            VALUE_ERROR("negative shift count");
            Py_DECREF(tmp);
            return NULL;
        }
        Py_DECREF(tmp);
        if (!(result = GMPy_MPZ_New(NULL))) {
            /* LCOV_EXCL_START */
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        mpz_set_si(result->z, mpz_sgn(MPZ(self)) < 0 ? -1 : 0);
        return (PyObject*)result;
    }

    if (!(result = GMPy_MPZ_New(NULL)))
        return NULL;

    if (CHECK_MPZANY(self)) {
        mpz_fdiv_q_2exp(result->z, MPZ(self), count);
        return (PyObject*)result;
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(self, NULL))) {
            Py_XDECREF((PyObject*)result);
            Py_XDECREF((PyObject*)tempx);
            return NULL;
        }

        mpz_fdiv_q_2exp(result->z, tempx->z, count);
        Py_DECREF((PyObject*)tempx);
        return (PyObject*)result;
    }
}

static PyObject *
GMPy_MPZ_Lshift_Slot(PyObject *self, PyObject *other)
{
    mp_bitcnt_t count;
    MPZ_Object *result, *tempx;

    count = GMPy_Integer_AsMpBitCnt(other);
    if ((count == (mp_bitcnt_t)(-1)) && PyErr_Occurred())
        return NULL;

    if (!(result = GMPy_MPZ_New(NULL)))
        return NULL;

    if (CHECK_MPZANY(self)) {
        mpz_mul_2exp(result->z, MPZ(self), count);
        return (PyObject*)result;
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(self, NULL))) {
            Py_XDECREF((PyObject*)result);
            Py_XDECREF((PyObject*)tempx);
            return NULL;
        }

        mpz_mul_2exp(result->z, tempx->z, count);
        Py_DECREF((PyObject*)tempx);
        return (PyObject*)result;
    }
}

PyDoc_STRVAR(doc_popcount,
"popcount($module, x, /)\n--\n\n"
"Return the number of 1-bits set in x. If x<0, the number of\n"
"1-bits is infinite so -1 is returned in that case.");

static PyObject *
GMPy_MPZ_popcount(PyObject *self, PyObject *other)
{
    mp_bitcnt_t n;
    MPZ_Object *tempx;

    if ((tempx = GMPy_MPZ_From_Integer(other, NULL))) {
        n = mpz_popcount(tempx->z);
        Py_DECREF((PyObject*)tempx);
        if (n == (mp_bitcnt_t)(-1))
            return PyLong_FromLong(-1);
        else
            return GMPy_PyLong_FromMpBitCnt(n);
    }
    else {
        TYPE_ERROR("popcount() requires 'mpz' argument");
        return NULL;
    }
}

PyDoc_STRVAR(doc_bit_count_method,
"bit_count($self)\n--\n\n"
"Return the number of 1-bits set in abs(self).");

static PyObject *
GMPy_MPZ_bit_count_method(PyObject *self, PyObject *other)
{
    mp_bitcnt_t n = 0;
    MPZ_Object *temppos = NULL;

    /* Since tempx could just be an additional reference to
     * an existing mpz, we need to make a copy to change the
     * sign.
     */

    if (mpz_sgn(MPZ(self)) == -1) {
        if (!(temppos = GMPy_MPZ_New(NULL))) {
            /* LCOV_EXCL_START */
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        mpz_abs(MPZ(temppos), MPZ(self));
        n = mpz_popcount(MPZ(temppos));
        Py_DECREF((PyObject*)temppos);
    }
    else {
        n = mpz_popcount(MPZ(self));
    }
    return GMPy_PyLong_FromMpBitCnt(n);
}

PyDoc_STRVAR(doc_bit_count,
"bit_count($module, x, /)\n--\n\n"
"Return the number of 1-bits set in abs(x).");

static PyObject *
GMPy_MPZ_bit_count(PyObject *self, PyObject *other)
{
    mp_bitcnt_t n;
    MPZ_Object *tempx = NULL, *temppos = NULL;


    if ((tempx = GMPy_MPZ_From_Integer(other, NULL))) {
        /* Since tempx could just be an additional reference to
         * an existing mpz, we need to make a copy to change the
         * sign.
         */
        if (mpz_sgn(MPZ(tempx)) == -1) {
            if (!(temppos = GMPy_MPZ_New(NULL))) {
                /* LCOV_EXCL_START */
                return NULL;
                /* LCOV_EXCL_STOP */
            }
            mpz_abs(MPZ(temppos), MPZ(tempx));
            n = mpz_popcount(MPZ(temppos));
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)temppos);
        }
        else {
            n = mpz_popcount(MPZ(tempx));
            Py_DECREF((PyObject*)tempx);
        }
        return GMPy_PyLong_FromMpBitCnt(n);
    }
    else {
        TYPE_ERROR("bit_count() requires 'mpz' argument");
        return NULL;
    }
}

PyDoc_STRVAR(doc_hamdist,
"hamdist($module, x, y, /)\n--\n\n"
"Return the Hamming distance (number of bit-positions where the\n"
"bits differ) between integers x and y.");

static PyObject *
GMPy_MPZ_hamdist(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;
    MPZ_Object *tempx = NULL, *tempy = NULL;

    if (PyTuple_GET_SIZE(args) != 2)
        goto err;

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    tempy = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    if (!tempx || !tempy)
        goto err;

    result = GMPy_PyLong_FromMpBitCnt(mpz_hamdist(tempx->z, tempy->z));
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    return result;

  err:
    TYPE_ERROR("hamdist() requires 'mpz','mpz' arguments");
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    return NULL;
}
