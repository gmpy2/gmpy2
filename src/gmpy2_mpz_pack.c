/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpz_pack.c                                                        *
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

/*
 **************************************************************************
 * pack and unpack methods
 *
 * Both pack and unpack use a devious trick when stuffing values into the
 * internal structure of an mpz_t. By setting a bit beyond the range of
 * interest, we are guaranteed that memory will remain available. When
 * the bit is cleared, it also triggers normalization of the value by
 * accounting for leading bits that are zero.
 **************************************************************************
 */

PyDoc_STRVAR(doc_pack,
"pack($module, lst, n, /)\n--\n\n"
"Pack a list of integers lst into a single `mpz` by concatenating\n"
"each integer element of lst after padding to length n bits. Raises\n"
"an error if any integer is negative or greater than n bits in\n"
"length.");

static PyObject *
GMPy_MPZ_pack(PyObject *self, PyObject *args)
{
    mp_bitcnt_t nbits, total_bits, tempx_bits;
    mp_size_t i;
    Py_ssize_t index, lst_count, temp_bits, limb_count;
    PyObject *lst;
    mpz_t temp, temp1;
    MPZ_Object *result, *tempx = 0;
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("pack() requires 'list','int' arguments");
        return NULL;
    }

    nbits = GMPy_Integer_AsMpBitCnt(PyTuple_GET_ITEM(args, 1));
    if (nbits == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if (!PyList_Check(PyTuple_GET_ITEM(args, 0))) {
        TYPE_ERROR("pack() requires 'list','int' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(context)))
        return NULL;

    lst = PyTuple_GET_ITEM(args, 0);
    lst_count = PyList_GET_SIZE(lst);
    total_bits = nbits * lst_count;

    if ((total_bits / lst_count) != nbits) {
        VALUE_ERROR("result too large to store in an 'mpz'");
        return NULL;
    }

    mpz_set_ui(result->z, 0);
    mpz_setbit(result->z, total_bits + (mp_bits_per_limb * 2));

    mpz_init(temp);
    mpz_init(temp1);
    mpz_set_ui(temp, 0);
    limb_count = 0;
    tempx_bits = 0;

    for (index = 0; index < lst_count; index++) {
        if (!(tempx = GMPy_MPZ_From_Integer(PyList_GetItem(lst, index), context))
            || (mpz_sgn(tempx->z) < 0)
            || (mpz_sizeinbase(tempx->z,2) > (size_t)nbits)) {
            TYPE_ERROR("pack() requires list elements be positive integers < 2^n bits");
            mpz_clear(temp);
            Py_XDECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        mpz_mul_2exp(temp1, tempx->z, tempx_bits);
        mpz_add(temp, temp, temp1);
        tempx_bits += nbits;
        i = 0;
        temp_bits = mpz_sizeinbase(temp, 2) * mpz_sgn(temp);
        while (tempx_bits >= (mp_bitcnt_t)mp_bits_per_limb) {
            if (temp_bits > 0) {
                result->z->_mp_d[limb_count] = mpz_getlimbn(temp, i);
            }
            i += 1;
            tempx_bits -= mp_bits_per_limb;
            limb_count += 1;
            temp_bits -= mp_bits_per_limb;
        }
        if (temp_bits > 0) {
            mpz_tdiv_q_2exp(temp, temp, mp_bits_per_limb * (mp_bitcnt_t)i);
        }
        else {
            mpz_set_ui(temp, 0);
        }
        Py_DECREF((PyObject*)tempx);
    }
    result->z->_mp_d[limb_count] = mpz_getlimbn(temp, 0);
    mpz_clrbit(result->z, total_bits + (mp_bits_per_limb * 2));
    mpz_clear(temp);
    mpz_clear(temp1);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_unpack,
"unpack($module, x, n, /)\n--\n\n"
"Unpack an integer x into a list of n-bit values. Equivalent to\n"
"repeated division by 2**n. Raises error if x is negative.");

static PyObject *
GMPy_MPZ_unpack(PyObject *self, PyObject *args)
{
    mp_bitcnt_t nbits, total_bits, guard_bit, extra_bits, temp_bits;
    mp_size_t index = 0;
    Py_ssize_t lst_count, i, lst_ptr = 0;
    PyObject *result;
    mpz_t temp;
    mp_limb_t extra = 0;
    MPZ_Object *item, *tempx = NULL;
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("unpack() requires 'int','int' arguments");
        return NULL;
    }

    nbits = GMPy_Integer_AsMpBitCnt(PyTuple_GET_ITEM(args, 1));
    if (nbits == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), context))) {
        TYPE_ERROR("unpack() requires 'int','int' arguments");
        return NULL;
    }

    if (mpz_sgn(tempx->z) < 0) {
        VALUE_ERROR("unpack() requires x >= 0");
        return NULL;
    }

    if (mpz_sgn(tempx->z) == 0) {
        total_bits = 0;
    }
    else {
        total_bits = mpz_sizeinbase(tempx->z, 2);
    }

    lst_count = total_bits / nbits;
    if ((total_bits % nbits) || !lst_count) {
        lst_count += 1;
    }

    if (!(result = PyList_New(lst_count))) {
        /* LCOV_EXCL_START */
        Py_DECREF((PyObject*)tempx);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (mpz_sgn(tempx->z) == 0) {
        if (!(item = GMPy_MPZ_New(context))) {
            /* LCOV_EXCL_START */
            Py_DECREF((PyObject*)tempx);
            Py_DECREF(result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        mpz_set_ui(item->z, 0);
        PyList_SET_ITEM(result, 0, (PyObject*)item);
        Py_DECREF((PyObject*)tempx);
        return result;
    }

    mpz_init(temp);
    guard_bit = nbits + (2 * mp_bits_per_limb);
    extra_bits = 0;
    index = 0;

    while (lst_ptr < lst_count) {
        i = 0;
        temp_bits = 0;
        mpz_set_ui(temp, 0);
        mpz_setbit(temp, guard_bit);
        while (temp_bits + extra_bits < nbits) {
            temp->_mp_d[i++] = mpz_getlimbn(tempx->z, index++);
            temp_bits += mp_bits_per_limb;
        }
        mpz_clrbit(temp, guard_bit);
        mpz_mul_2exp(temp, temp, extra_bits);
        if (mpz_sgn(temp) == 0 && extra != 0) {
            mpz_set_ui(temp, 1);
            temp->_mp_d[0] = extra;
        }
        else {
           mpn_add_1(temp->_mp_d, temp->_mp_d, (mp_size_t)mpz_size(temp), extra);
        }
        temp_bits += extra_bits;

        while ((lst_ptr < lst_count) && (temp_bits >= nbits)) {
            if(!(item = GMPy_MPZ_New(context))) {
                /* LCOV_EXCL_START */
                mpz_clear(temp);
                Py_DECREF((PyObject*)tempx);
                Py_DECREF(result);
                return NULL;
                /* LCOV_EXCL_STOP */
            }
            mpz_tdiv_r_2exp(item->z, temp, nbits);
            PyList_SET_ITEM(result, lst_ptr++, (PyObject*)item);
            mpz_tdiv_q_2exp(temp, temp, nbits);
            temp_bits -= nbits;
        }
        extra = mpz_getlimbn(temp, 0);
        extra_bits = temp_bits;
    }
    Py_DECREF((PyObject*)tempx);
    mpz_clear(temp);
    return result;
}
