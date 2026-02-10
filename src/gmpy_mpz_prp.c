/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz_prp.c                                                          *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP, MPFR, and MPC multiple precision           *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2011 David Cleaver                                            *
 *                                                                         *
 * Copyright 2012 - 2025 Case Van Horsen                                   *
 *                                                                         *
 * The original file is available at:                                      *
 *   <http://sourceforge.net/projects/mpzprp/files/>                       *
 *                                                                         *
 * Modified by Case Van Horsen for inclusion into GMPY2.                   *
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

/* ******************************************************************
 * mpz_prp: (also called a Fermat probable prime)
 * A "probable prime" to the base a is a number n such that,
 * (a,n)=1 and a^(n-1) = 1 mod n
 * ******************************************************************/

PyDoc_STRVAR(doc_mpz_is_fermat_prp,
"is_fermat_prp($module, n, a, /)\n--\n\n"
"Return `True` if n is a Fermat probable prime to the base a.\n"
"Assuming:\n\n"
"    gcd(n,a) == 1\n\n"
"Then a Fermat probable prime requires:\n\n"
"    a**(n-1) == 1 (mod n)");

static PyObject *
GMPY_mpz_is_fermat_prp(PyObject *self, PyObject *args)
{
    MPZ_Object *a = NULL, *n = NULL;
    PyObject *result = NULL;
    mpz_t res, nm1;

    if (PyTuple_Size(args) != 2) {
        TYPE_ERROR("is_fermat_prp() requires 2 integer arguments");
        return NULL;
    }

    mpz_init(res);
    mpz_init(nm1);

    n = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    a = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    if (!a || !n) {
        TYPE_ERROR("is_fermat_prp() requires 2 integer arguments");
        goto cleanup;
    }

    /* Require a >= 2. */
    if (mpz_cmp_ui(a->z, 2) < 0) {
        VALUE_ERROR("is_fermat_prp() requires 'a' greater than or equal to 2");
        goto cleanup;
    }

    /* Require n > 0. */
    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("is_fermat_prp() requires 'n' be greater than 0");
        goto cleanup;
    }

    /* Check for n == 1 */
    if (mpz_cmp_ui(n->z, 1) == 0) {
        result = Py_False;
        goto cleanup;
    }

    /* Handle n even. */
    /* Should n even raise an exception? */
    if (mpz_divisible_ui_p(n->z, 2)) {
        if (mpz_cmp_ui(n->z, 2) == 0)
            result = Py_True;
        else
            result = Py_False;
        goto cleanup;
    }

    /* Check gcd(a,n) */
    mpz_gcd(res, n->z, a->z);
    if (mpz_cmp_ui(res, 1) > 0) {
        VALUE_ERROR("is_fermat_prp() requires gcd(n,a) == 1");
        goto cleanup;
    }

    mpz_set(nm1, n->z);
    mpz_sub_ui(nm1, nm1, 1);
    mpz_powm(res, a->z, nm1, n->z);

    if (mpz_cmp_ui(res, 1) == 0)
        result = Py_True;
    else
        result = Py_False;

  cleanup:
    Py_XINCREF(result);
    mpz_clear(res);
    mpz_clear(nm1);
    Py_XDECREF((PyObject*)a);
    Py_XDECREF((PyObject*)n);
    return result;
}

/* *************************************************************************
 * mpz_euler_prp: (also called a Solovay-Strassen probable prime)
 * An "Euler probable prime" to the base a is an odd composite number n with,
 * (a,n)=1 such that a^((n-1)/2)=(a/n) mod n [(a/n) is the Jacobi symbol]
 * *************************************************************************/

PyDoc_STRVAR(doc_mpz_is_euler_prp,
"is_euler_prp($module, n, a, /)\n--\n\n"
"Return `True` if n is an Euler (also known as Solovay-Strassen)\n"
"probable prime to the base a.\n"
"Assuming:\n\n"
"    gcd(n,a) == 1\n"
"    n is odd\n\n"
"Then an Euler probable prime requires:\n\n"
"    a**((n-1)/2) == (a/n) (mod n)\n\n"
"where (a/n) is the Jacobi symbol.");

static PyObject *
GMPY_mpz_is_euler_prp(PyObject *self, PyObject *args)
{
    MPZ_Object *a = NULL, *n = NULL;
    PyObject *result = NULL;
    mpz_t res, exp;
    int ret;

    if (PyTuple_Size(args) != 2) {
        TYPE_ERROR("is_euler_prp() requires 2 integer arguments");
        return NULL;
    }

    mpz_init(res);
    mpz_init(exp);

    n = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    a = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    if (!a || !n) {
        TYPE_ERROR("is_euler_prp() requires 2 integer arguments");
        goto cleanup;
    }

    /* Require a >= 2. */
    if (mpz_cmp_ui(a->z, 2) < 0) {
        VALUE_ERROR("is_euler_prp() requires 'a' greater than or equal to 2");
        goto cleanup;
    }

    /* Require n > 0. */
    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("is_euler_prp() requires 'n' be greater than 0");
        goto cleanup;
    }

    /* Check for n == 1 */
    if (mpz_cmp_ui(n->z, 1) == 0) {
        result = Py_False;
        goto cleanup;
    }

    /* Handle n even. */
    if (mpz_divisible_ui_p(n->z, 2)) {
        if (mpz_cmp_ui(n->z, 2) == 0)
            result = Py_True;
        else
            result = Py_False;
        goto cleanup;
    }

    /* Check gcd(a,b) */
    mpz_gcd(res, n->z, a->z);
    if (mpz_cmp_ui(res, 1) > 0) {
        VALUE_ERROR("is_euler_prp() requires gcd(n,a) == 1");
        goto cleanup;
    }

    mpz_set(exp, n->z);
    mpz_sub_ui(exp, exp, 1);
    mpz_divexact_ui(exp, exp, 2);
    mpz_powm(res, a->z, exp, n->z);

    /* reuse exp to calculate jacobi(a,n) mod n */
    ret = mpz_jacobi(a->z,n->z);
    mpz_set(exp, n->z);
    if (ret == -1)
        mpz_sub_ui(exp, exp, 1);
    else if (ret == 1)
        mpz_add_ui(exp, exp, 1);
    mpz_mod(exp, exp, n->z);

    if (mpz_cmp(res, exp) == 0)
        result = Py_True;
    else
        result = Py_False;

  cleanup:
    Py_XINCREF(result);
    mpz_clear(res);
    mpz_clear(exp);
    Py_XDECREF((PyObject*)a);
    Py_XDECREF((PyObject*)n);
    return result;
}

/* *********************************************************************************************
 * mpz_sprp: (also called a Miller-Rabin probable prime)
 * A "strong probable prime" to the base a is an odd composite n = (2^r)*s+1 with s odd such that
 * either a^s == 1 mod n, or a^((2^t)*s) == -1 mod n, for some integer t, with 0 <= t < r.
 * *********************************************************************************************/

PyDoc_STRVAR(doc_mpz_is_strong_prp,
"is_strong_prp($module, n, a, /)\n--\n\n"
"Return `True` if n is a strong (also known as Miller-Rabin)\n"
"probable prime to the base a.\n"
"Assuming:\n\n"
"    gcd(n,a) == 1\n"
"    n is odd\n"
"    n = s*(2**r) + 1, with s odd\n\n"
"Then a strong probable prime requires one of the following is true:\n\n"
"    a**s == 1 (mod n)\n"
"    or\n"
"    a**(s*(2**t)) == -1 (mod n) for some t, 0 <= t < r.");

static PyObject *
GMPY_mpz_is_strong_prp(PyObject *self, PyObject *args)
{
    MPZ_Object *a = NULL, *n = NULL;
    PyObject *result = NULL;
    mpz_t s, nm1, mpz_test;
    mp_bitcnt_t r = 0;

    if (PyTuple_Size(args) != 2) {
        TYPE_ERROR("is_strong_prp() requires 2 integer arguments");
        return NULL;
    }

    mpz_init(s);
    mpz_init(nm1);
    mpz_init(mpz_test);

    n = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    a = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    if (!a || !n) {
        TYPE_ERROR("is_strong_prp() requires 2 integer arguments");
        goto cleanup;
    }

    /* Require a >= 2. */
    if (mpz_cmp_ui(a->z, 2) < 0) {
        VALUE_ERROR("is_strong_prp() requires 'a' greater than or equal to 2");
        goto cleanup;
    }

    /* Require n > 0. */
    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("is_strong_prp() requires 'n' be greater than 0");
        goto cleanup;
    }

    /* Check for n == 1 */
    if (mpz_cmp_ui(n->z, 1) == 0) {
        result = Py_False;
        goto cleanup;
    }

    /* Handle n even. */
    if (mpz_divisible_ui_p(n->z, 2)) {
        if (mpz_cmp_ui(n->z, 2) == 0)
            result = Py_True;
        else
            result = Py_False;
        goto cleanup;
    }

    /* Check gcd(a,b) */
    mpz_gcd(s, n->z, a->z);
    if (mpz_cmp_ui(s, 1) > 0) {
        VALUE_ERROR("is_strong_prp() requires gcd(n,a) == 1");
        goto cleanup;
    }

    mpz_set(nm1, n->z);
    mpz_sub_ui(nm1, nm1, 1);

    /* Find s and r satisfying: n-1=(2^r)*s, s odd */
    r = mpz_scan1(nm1, 0);
    mpz_fdiv_q_2exp(s, nm1, r);


    /* Check a^((2^t)*s) mod n for 0 <= t < r */
    mpz_powm(mpz_test, a->z, s, n->z);
    if ((mpz_cmp_ui(mpz_test, 1) == 0) || (mpz_cmp(mpz_test, nm1) == 0)) {
        result = Py_True;
        goto cleanup;
    }

    while (--r) {
        /* mpz_test = mpz_test^2%n */
        mpz_mul(mpz_test, mpz_test, mpz_test);
        mpz_mod(mpz_test, mpz_test, n->z);

        if (mpz_cmp(mpz_test, nm1) == 0) {
            result = Py_True;
            goto cleanup;
        }
    }

    result = Py_False;
  cleanup:
    Py_XINCREF(result);
    mpz_clear(s);
    mpz_clear(nm1);
    mpz_clear(mpz_test);
    Py_XDECREF((PyObject*)a);
    Py_XDECREF((PyObject*)n);
    return result;
}

/* *************************************************************************
 * mpz_fibonacci_prp:
 * A "Fibonacci probable prime" with parameters (P,Q), P > 0, Q=+/-1, is a
 * composite n for which V_n == P mod n
 * [V is the Lucas V sequence with parameters P,Q]
 * *************************************************************************/

PyDoc_STRVAR(doc_mpz_is_fibonacci_prp,
"is_fibonacci_prp($module, n, p, q, /)\n--\n\n"
"Return `True` if n is a Fibonacci probable prime with parameters (p,q).\n"
"Assuming:\n\n"
"    n is odd\n"
"    p > 0, q = +/-1\n"
"    p*p - 4*q != 0\n\n"
"Then a Fibonacci probable prime requires:\n\n"
"    lucasv(p,q,n) == p (mod n).");

static PyObject *
GMPY_mpz_is_fibonacci_prp(PyObject *self, PyObject *args)
{
    MPZ_Object *n = NULL, *p = NULL, *q = NULL;
    PyObject *result = NULL;
    mpz_t pmodn, zP;
    /* used for calculating the Lucas V sequence */
    mpz_t vl, vh, ql, qh, tmp;
    mp_bitcnt_t s = 0, j = 0;

    if (PyTuple_Size(args) != 3) {
        TYPE_ERROR("is_fibonacci_prp() requires 3 integer arguments");
        return NULL;
    }

    mpz_init(pmodn);
    mpz_init(zP);
    mpz_init(vl);
    mpz_init(vh);
    mpz_init(ql);
    mpz_init(qh);
    mpz_init(tmp);

    n = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    p = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    q = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 2), NULL);
    if (!n || !p || !q) {
        TYPE_ERROR("is_fibonacci_prp() requires 3 integer arguments");
        goto cleanup;
    }

    /* Check if p*p - 4*q == 0. */

    mpz_mul(tmp, p->z, p->z);
    mpz_mul_ui(qh, q->z, 4);
    mpz_sub(tmp, tmp, qh);
    if (mpz_sgn(tmp) == 0) {
        VALUE_ERROR("invalid values for p,q in is_fibonacci_prp()");
        goto cleanup;
    }

    /* Verify q = +/-1 */

    if ((mpz_cmp_si(q->z, 1) && mpz_cmp_si(q->z, -1)) || (mpz_sgn(p->z) <= 0)) {
        VALUE_ERROR("invalid values for p,q in is_fibonacci_prp()");
        goto cleanup;
    }

    /* Require n > 0. */
    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("is_fibonacci_prp() requires 'n' be greater than 0");
        goto cleanup;
    }

    /* Check for n == 1 */
    if (mpz_cmp_ui(n->z, 1) == 0) {
        result = Py_False;
        goto cleanup;
    }

    /* Handle n even. */
    if (mpz_divisible_ui_p(n->z, 2)) {
        if (mpz_cmp_ui(n->z, 2) == 0)
            result = Py_True;
        else
            result = Py_False;
        goto cleanup;
    }

    mpz_set(zP, p->z);
    mpz_mod(pmodn, zP, n->z);

    /* mpz_lucasvmod(res, p, q, n, n); */
    mpz_set_si(vl, 2);
    mpz_set(vh, p->z);
    mpz_set_si(ql, 1);
    mpz_set_si(qh, 1);
    mpz_set_si(tmp,0);

    s = mpz_scan1(n->z, 0);
    for (j = mpz_sizeinbase(n->z,2)-1; j >= s+1; j--) {
        /* ql = ql*qh (mod n) */
        mpz_mul(ql, ql, qh);
        mpz_mod(ql, ql, n->z);
        if (mpz_tstbit(n->z,j) == 1) {
            /* qh = ql*q */
            mpz_mul(qh, ql, q->z);

            /* vl = vh*vl - p*ql (mod n) */
            mpz_mul(vl, vh, vl);
            mpz_mul(tmp, ql, p->z);
            mpz_sub(vl, vl, tmp);
            mpz_mod(vl, vl, n->z);

            /* vh = vh*vh - 2*qh (mod n) */
            mpz_mul(vh, vh, vh);
            mpz_mul_si(tmp, qh, 2);
            mpz_sub(vh, vh, tmp);
            mpz_mod(vh, vh, n->z);
        }
        else {
            /* qh = ql */
            mpz_set(qh, ql);

            /* vh = vh*vl - p*ql (mod n) */
            mpz_mul(vh, vh, vl);
            mpz_mul(tmp, ql, p->z);
            mpz_sub(vh, vh, tmp);
            mpz_mod(vh, vh, n->z);

            /* vl = vl*vl - 2*ql (mod n) */
            mpz_mul(vl, vl, vl);
            mpz_mul_si(tmp, ql, 2);
            mpz_sub(vl, vl, tmp);
            mpz_mod(vl, vl, n->z);
        }
    }
    /* ql = ql*qh */
    mpz_mul(ql, ql, qh);

    /* qh = ql*q */
    mpz_mul(qh, ql, q->z);

    /* vl = vh*vl - p*ql */
    mpz_mul(vl, vh, vl);
    mpz_mul(tmp, ql, p->z);
    mpz_sub(vl, vl, tmp);

    /* ql = ql*qh */
    mpz_mul(ql, ql, qh);

    for (j = 1; j <= s; j++) {
        /* vl = vl*vl - 2*ql (mod n) */
        mpz_mul(vl, vl, vl);
        mpz_mul_si(tmp, ql, 2);
        mpz_sub(vl, vl, tmp);
        mpz_mod(vl, vl, n->z);

        /* ql = ql*ql (mod n) */
        mpz_mul(ql, ql, ql);
        mpz_mod(ql, ql, n->z);
    }

    /* vl contains our return value */
    mpz_mod(vl, vl, n->z);

    if (mpz_cmp(vl, pmodn) == 0)
        result = Py_True;
    else
        result = Py_False;

  cleanup:
    Py_XINCREF(result);
    mpz_clear(pmodn);
    mpz_clear(zP);
    mpz_clear(vl);
    mpz_clear(vh);
    mpz_clear(ql);
    mpz_clear(qh);
    mpz_clear(tmp);
    Py_XDECREF((PyObject*)p);
    Py_XDECREF((PyObject*)q);
    Py_XDECREF((PyObject*)n);
    return result;
}


/* *******************************************************************************
 * mpz_lucas_prp:
 * A "Lucas probable prime" with parameters (P,Q) is a composite n with D=P^2-4Q,
 * (n,2QD)=1 such that U_(n-(D/n)) == 0 mod n [(D/n) is the Jacobi symbol]
 * *******************************************************************************/

PyDoc_STRVAR(doc_mpz_is_lucas_prp,
"is_lucas_prp($module, n, p, q, /)\n--\n\n"
"Return `True` if n is a Lucas probable prime with parameters (p,q).\n"
"Assuming:\n\n"
"    n is odd\n"
"    D = p*p - 4*q, D != 0\n"
"    gcd(n, 2*q*D) == 1\n\n"
"Then a Lucas probable prime requires:\n\n"
"    lucasu(p,q,n - Jacobi(D,n)) == 0 (mod n)");

static PyObject *
GMPY_mpz_is_lucas_prp(PyObject *self, PyObject *args)
{
    MPZ_Object *n = NULL, *p = NULL, *q = NULL;
    PyObject *result = NULL;
    mpz_t zD, res, index;
    /* used for calculating the Lucas U sequence */
    mpz_t uh, vl, vh, ql, qh, tmp;
    mp_bitcnt_t s = 0, j = 0;
    int ret;

    if (PyTuple_Size(args) != 3) {
        TYPE_ERROR("is_lucas_prp() requires 3 integer arguments");
        return NULL;
    }

    mpz_init(zD);
    mpz_init(res);
    mpz_init(index);
    mpz_init(uh);
    mpz_init(vl);
    mpz_init(vh);
    mpz_init(ql);
    mpz_init(qh);
    mpz_init(tmp);

    n = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    p = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    q = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 2), NULL);
    if (!n || !p || !q) {
        TYPE_ERROR("is_lucas_prp() requires 3 integer arguments");
        goto cleanup;
    }

    /* Check if p*p - 4*q == 0. */
    mpz_mul(zD, p->z, p->z);
    mpz_mul_ui(tmp, q->z, 4);
    mpz_sub(zD, zD, tmp);
    if (mpz_sgn(zD) == 0) {
        VALUE_ERROR("invalid values for p,q in is_lucas_prp()");
        goto cleanup;
    }

    /* Require n > 0. */
    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("is_lucas_prp() requires 'n' be greater than 0");
        goto cleanup;
    }

    /* Check for n == 1 */
    if (mpz_cmp_ui(n->z, 1) == 0) {
        result = Py_False;
        goto cleanup;
    }

    /* Handle n even. */
    if (mpz_divisible_ui_p(n->z, 2)) {
        if (mpz_cmp_ui(n->z, 2) == 0)
            result = Py_True;
        else
            result = Py_False;
        goto cleanup;
    }

    /* Check GCD */
    mpz_mul(res, zD, q->z);
    mpz_mul_ui(res, res, 2);
    mpz_gcd(res, res, n->z);
    if ((mpz_cmp(res, n->z) != 0) && (mpz_cmp_ui(res, 1) > 0)) {
        VALUE_ERROR("is_lucas_prp() requires gcd(n,2*q*D) == 1");
        goto cleanup;
    }

    /* index = n-(D/n), where (D/n) is the Jacobi symbol */
    mpz_set(index, n->z);
    ret = mpz_jacobi(zD, n->z);
    if (ret == -1)
        mpz_add_ui(index, index, 1);
    else if (ret == 1)
        mpz_sub_ui(index, index, 1);

    /* mpz_lucasumod(res, p, q, index, n); */
    mpz_set_si(uh, 1);
    mpz_set_si(vl, 2);
    mpz_set(vh, p->z);
    mpz_set_si(ql, 1);
    mpz_set_si(qh, 1);
    mpz_set_si(tmp,0);

    s = mpz_scan1(index, 0);
    for (j = mpz_sizeinbase(index,2)-1; j >= s+1; j--) {
        /* ql = ql*qh (mod n) */
        mpz_mul(ql, ql, qh);
        mpz_mod(ql, ql, n->z);
        if (mpz_tstbit(index,j) == 1) {
            /* qh = ql*q */
            mpz_mul(qh, ql, q->z);

            /* uh = uh*vh (mod n) */
            mpz_mul(uh, uh, vh);
            mpz_mod(uh, uh, n->z);

            /* vl = vh*vl - p*ql (mod n) */
            mpz_mul(vl, vh, vl);
            mpz_mul(tmp, ql, p->z);
            mpz_sub(vl, vl, tmp);
            mpz_mod(vl, vl, n->z);

            /* vh = vh*vh - 2*qh (mod n) */
            mpz_mul(vh, vh, vh);
            mpz_mul_si(tmp, qh, 2);
            mpz_sub(vh, vh, tmp);
            mpz_mod(vh, vh, n->z);
        }
        else {
            /* qh = ql */
            mpz_set(qh, ql);

            /* uh = uh*vl - ql (mod n) */
            mpz_mul(uh, uh, vl);
            mpz_sub(uh, uh, ql);
            mpz_mod(uh, uh, n->z);

            /* vh = vh*vl - p*ql (mod n) */
            mpz_mul(vh, vh, vl);
            mpz_mul(tmp, ql, p->z);
            mpz_sub(vh, vh, tmp);
            mpz_mod(vh, vh, n->z);

            /* vl = vl*vl - 2*ql (mod n) */
            mpz_mul(vl, vl, vl);
            mpz_mul_si(tmp, ql, 2);
            mpz_sub(vl, vl, tmp);
            mpz_mod(vl, vl, n->z);
        }
    }
    /* ql = ql*qh */
    mpz_mul(ql, ql, qh);

    /* qh = ql*q */
    mpz_mul(qh, ql, q->z);

    /* uh = uh*vl - ql */
    mpz_mul(uh, uh, vl);
    mpz_sub(uh, uh, ql);

    /* vl = vh*vl - p*ql */
    mpz_mul(vl, vh, vl);
    mpz_mul(tmp, ql, p->z);
    mpz_sub(vl, vl, tmp);

    /* ql = ql*qh */
    mpz_mul(ql, ql, qh);

    for (j = 1; j <= s; j++) {
        /* uh = uh*vl (mod n) */
        mpz_mul(uh, uh, vl);
        mpz_mod(uh, uh, n->z);

        /* vl = vl*vl - 2*ql (mod n) */
        mpz_mul(vl, vl, vl);
        mpz_mul_si(tmp, ql, 2);
        mpz_sub(vl, vl, tmp);
        mpz_mod(vl, vl, n->z);

        /* ql = ql*ql (mod n) */
        mpz_mul(ql, ql, ql);
        mpz_mod(ql, ql, n->z);
    }

    /* uh contains our return value */
    mpz_mod(res, uh, n->z);
    if (mpz_cmp_ui(res, 0) == 0)
        result = Py_True;
    else
        result = Py_False;

  cleanup:
    Py_XINCREF(result);
    mpz_clear(zD);
    mpz_clear(res);
    mpz_clear(index);
    mpz_clear(uh);
    mpz_clear(vl);
    mpz_clear(vh);
    mpz_clear(ql);
    mpz_clear(qh);
    mpz_clear(tmp);
    Py_XDECREF((PyObject*)p);
    Py_XDECREF((PyObject*)q);
    Py_XDECREF((PyObject*)n);
    return result;
}

/* *********************************************************************************************
 * mpz_stronglucas_prp:
 * A "strong Lucas probable prime" with parameters (P,Q) is a composite n = (2^r)*s+(D/n), where
 * s is odd, D=P^2-4Q, and (n,2QD)=1 such that either U_s == 0 mod n or V_((2^t)*s) == 0 mod n
 * for some t, 0 <= t < r. [(D/n) is the Jacobi symbol]
 * *********************************************************************************************/

PyDoc_STRVAR(doc_mpz_is_stronglucas_prp,
"is_strong_lucas_prp($module, n, p, q, /)\n--\n\n"
"Return `True` if n is a strong Lucas probable prime with parameters (p,q).\n"
"Assuming:\n\n"
"    n is odd\n"
"    D = p*p - 4*q, D != 0\n"
"    gcd(n, 2*q*D) == 1\n"
"    n = s*(2**r) + Jacobi(D,n), s odd\n\n"
"Then a strong Lucas probable prime requires:\n\n"
"    lucasu(p,q,s) == 0 (mod n)\n"
"    or\n"
"    lucasv(p,q,s*(2**t)) == 0 (mod n) for some t, 0 <= t < r");

static PyObject *
GMPY_mpz_is_stronglucas_prp(PyObject *self, PyObject *args)
{
    MPZ_Object *n = NULL, *p= NULL, *q = NULL;
    PyObject *result = NULL;
    mpz_t zD, s, nmj, res;
    /* these are needed for the LucasU and LucasV part of this function */
    mpz_t uh, vl, vh, ql, qh, tmp;
    mp_bitcnt_t r = 0, j = 0;
    int ret = 0;

    if (PyTuple_Size(args) != 3) {
        TYPE_ERROR("is_strong_lucas_prp() requires 3 integer arguments");
        return NULL;
    }

    mpz_init(zD);
    mpz_init(s);
    mpz_init(nmj);
    mpz_init(res);
    mpz_init(uh);
    mpz_init(vl);
    mpz_init(vh);
    mpz_init(ql);
    mpz_init(qh);
    mpz_init(tmp);

    n = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    p = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    q = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 2), NULL);
    if (!n || !p || !q) {
        TYPE_ERROR("is_strong_lucas_prp() requires 3 integer arguments");
        goto cleanup;
    }

    /* Check if p*p - 4*q == 0. */
    mpz_mul(zD, p->z, p->z);
    mpz_mul_ui(tmp, q->z, 4);
    mpz_sub(zD, zD, tmp);
    if (mpz_sgn(zD) == 0) {
        VALUE_ERROR("invalid values for p,q in is_strong_lucas_prp()");
        goto cleanup;
    }

    /* Require n > 0. */
    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("is_strong_lucas_prp() requires 'n' be greater than 0");
        goto cleanup;
    }

    /* Check for n == 1 */
    if (mpz_cmp_ui(n->z, 1) == 0) {
        result = Py_False;
        goto cleanup;
    }

    /* Handle n even. */
    if (mpz_divisible_ui_p(n->z, 2)) {
        if (mpz_cmp_ui(n->z, 2) == 0)
            result = Py_True;
        else
            result = Py_False;
        goto cleanup;
    }

    /* Check GCD */
    mpz_mul(res, zD, q->z);
    mpz_mul_ui(res, res, 2);
    mpz_gcd(res, res, n->z);
    if ((mpz_cmp(res, n->z) != 0) && (mpz_cmp_ui(res, 1) > 0)) {
        VALUE_ERROR("is_strong_lucas_prp() requires gcd(n,2*q*D) == 1");
        goto cleanup;
    }

    /* nmj = n - (D/n), where (D/n) is the Jacobi symbol */
    mpz_set(nmj, n->z);
    ret = mpz_jacobi(zD, n->z);
    if (ret == -1)
        mpz_add_ui(nmj, nmj, 1);
    else if (ret == 1)
        mpz_sub_ui(nmj, nmj, 1);

    r = mpz_scan1(nmj, 0);
    mpz_fdiv_q_2exp(s, nmj, r);

    /* make sure U_s == 0 mod n or V_((2^t)*s) == 0 mod n, for some t, 0 <= t < r */
    mpz_set_si(uh, 1);
    mpz_set_si(vl, 2);
    mpz_set(vh, p->z);
    mpz_set_si(ql, 1);
    mpz_set_si(qh, 1);
    mpz_set_si(tmp,0);

    for (j = mpz_sizeinbase(s,2)-1; j >= 1; j--) {
        /* ql = ql*qh (mod n) */
        mpz_mul(ql, ql, qh);
        mpz_mod(ql, ql, n->z);
        if (mpz_tstbit(s,j) == 1) {
            /* qh = ql*q */
            mpz_mul(qh, ql, q->z);

            /* uh = uh*vh (mod n) */
            mpz_mul(uh, uh, vh);
            mpz_mod(uh, uh, n->z);

            /* vl = vh*vl - p*ql (mod n) */
            mpz_mul(vl, vh, vl);
            mpz_mul(tmp, ql, p->z);
            mpz_sub(vl, vl, tmp);
            mpz_mod(vl, vl, n->z);

            /* vh = vh*vh - 2*qh (mod n) */
            mpz_mul(vh, vh, vh);
            mpz_mul_si(tmp, qh, 2);
            mpz_sub(vh, vh, tmp);
            mpz_mod(vh, vh, n->z);
        }
        else {
            /* qh = ql */
            mpz_set(qh, ql);

            /* uh = uh*vl - ql (mod n) */
            mpz_mul(uh, uh, vl);
            mpz_sub(uh, uh, ql);
            mpz_mod(uh, uh, n->z);

            /* vh = vh*vl - p*ql (mod n) */
            mpz_mul(vh, vh, vl);
            mpz_mul(tmp, ql, p->z);
            mpz_sub(vh, vh, tmp);
            mpz_mod(vh, vh, n->z);

            /* vl = vl*vl - 2*ql (mod n) */
            mpz_mul(vl, vl, vl);
            mpz_mul_si(tmp, ql, 2);
            mpz_sub(vl, vl, tmp);
            mpz_mod(vl, vl, n->z);
        }
    }
    /* ql = ql*qh */
    mpz_mul(ql, ql, qh);

    /* qh = ql*q */
    mpz_mul(qh, ql, q->z);

    /* uh = uh*vl - ql */
    mpz_mul(uh, uh, vl);
    mpz_sub(uh, uh, ql);

    /* vl = vh*vl - p*ql */
    mpz_mul(vl, vh, vl);
    mpz_mul(tmp, ql, p->z);
    mpz_sub(vl, vl, tmp);

    /* ql = ql*qh */
    mpz_mul(ql, ql, qh);

    mpz_mod(uh, uh, n->z);
    mpz_mod(vl, vl, n->z);

    /* uh contains LucasU_s and vl contains LucasV_s */
    if ((mpz_cmp_ui(uh, 0) == 0) || (mpz_cmp_ui(vl, 0) == 0)) {
        result = Py_True;
        goto cleanup;
    }

    for (j = 1; j < r; j++) {
        /* vl = vl*vl - 2*ql (mod n) */
        mpz_mul(vl, vl, vl);
        mpz_mul_si(tmp, ql, 2);
        mpz_sub(vl, vl, tmp);
        mpz_mod(vl, vl, n->z);

        /* ql = ql*ql (mod n) */
        mpz_mul(ql, ql, ql);
        mpz_mod(ql, ql, n->z);

        if (mpz_cmp_ui(vl, 0) == 0) {
            result = Py_True;
            goto cleanup;
        }
    }

    result = Py_False;
  cleanup:
    Py_XINCREF(result);
    mpz_clear(zD);
    mpz_clear(s);
    mpz_clear(nmj);
    mpz_clear(res);
    mpz_clear(uh);
    mpz_clear(vl);
    mpz_clear(vh);
    mpz_clear(ql);
    mpz_clear(qh);
    mpz_clear(tmp);
    Py_XDECREF((PyObject*)p);
    Py_XDECREF((PyObject*)q);
    Py_XDECREF((PyObject*)n);
    return result;
}

/* *******************************************************************************************
 * mpz_extrastronglucas_prp:
 * Let U_n = LucasU(p,1), V_n = LucasV(p,1), and D=p^2-4.
 * An "extra strong Lucas probable prime" to the base p is a composite n = (2^r)*s+(D/n), where
 * s is odd and (n,2D)=1, such that either U_s == 0 mod n and V_s == +/-2 mod n, or
 * V_((2^t)*s) == 0 mod n for some t with 0 <= t < r-1 [(D/n) is the Jacobi symbol]
 * *******************************************************************************************/

PyDoc_STRVAR(doc_mpz_is_extrastronglucas_prp,
"is_extra_strong_lucas_prp($module, n, p, /)\n--\n\n"
"Return `True` if n is an extra strong Lucas probable prime with parameters\n"
"(p,1). Assuming:\n\n"
"    n is odd\n"
"    D = p*p - 4, D != 0\n"
"    gcd(n, 2*D) == 1\n"
"    n = s*(2**r) + Jacobi(D,n), s odd\n\n"
"Then an extra strong Lucas probable prime requires:\n\n"
"    lucasu(p,1,s) == 0 (mod n) and lucasv(p,1,s) == +/-2 (mod n)\n"
"    or\n"
"    lucasv(p,1,s*(2**t)) == 0 (mod n) for some t, 0 <= t < r");

static PyObject *
GMPY_mpz_is_extrastronglucas_prp(PyObject *self, PyObject *args)
{
    MPZ_Object *n = NULL, *p = NULL;
    PyObject *result = NULL;
    mpz_t zD, s, nmj, nm2, res;
    /* these are needed for the LucasU and LucasV part of this function */
    mpz_t uh, vl, vh, ql, qh, tmp;
    mp_bitcnt_t r = 0, j = 0;
    long int q = 1;
    int ret = 0;

    if (PyTuple_Size(args) != 2) {
        TYPE_ERROR("is_extra_strong_lucas_prp() requires 2 integer arguments");
        return NULL;
    }

    mpz_init(zD);
    mpz_init(s);
    mpz_init(nmj);
    mpz_init(nm2);
    mpz_init(res);
    mpz_init(uh);
    mpz_init(vl);
    mpz_init(vh);
    mpz_init(ql);
    mpz_init(qh);
    mpz_init(tmp);

    n = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    p = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    if (!n || !p) {
        TYPE_ERROR("is_extra_strong_lucas_prp() requires 2 integer arguments");
        goto cleanup;
    }

    /* Check if p*p - 4 == 0. */
    mpz_mul(zD, p->z, p->z);
    mpz_sub_ui(zD, zD, 4);
    if (mpz_sgn(zD) == 0) {
        VALUE_ERROR("invalid value for p in is_extra_strong_lucas_prp()");
        goto cleanup;
    }

    /* Require n > 0. */
    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("is_extra_strong_lucas_prp() requires 'n' be greater than 0");
        goto cleanup;
    }

    /* Check for n == 1 */
    if (mpz_cmp_ui(n->z, 1) == 0) {
        result = Py_False;
        goto cleanup;
    }

    /* Handle n even. */
    if (mpz_divisible_ui_p(n->z, 2)) {
        if (mpz_cmp_ui(n->z, 2) == 0)
            result = Py_True;
        else
            result = Py_False;
        goto cleanup;
    }

    /* Check GCD */
    mpz_mul_ui(res, zD, 2);
    mpz_gcd(res, res, n->z);
    if ((mpz_cmp(res, n->z) != 0) && (mpz_cmp_ui(res, 1) > 0)) {
        VALUE_ERROR("is_extra_strong_lucas_prp() requires gcd(n,2*D) == 1");
        goto cleanup;
    }

    /* nmj = n - (D/n), where (D/n) is the Jacobi symbol */
    mpz_set(nmj, n->z);
    ret = mpz_jacobi(zD, n->z);
    if (ret == -1)
        mpz_add_ui(nmj, nmj, 1);
    else if (ret == 1)
        mpz_sub_ui(nmj, nmj, 1);

    r = mpz_scan1(nmj, 0);
    mpz_fdiv_q_2exp(s, nmj, r);

    mpz_set(nm2, n->z);
    mpz_sub_ui(nm2, nm2, 2);

    /* make sure that either U_s == 0 mod n and V_s == +/-2 mod n, or */
    /* V_((2^t)*s) == 0 mod n for some t with 0 <= t < r-1           */
    mpz_set_si(uh, 1);
    mpz_set_si(vl, 2);
    mpz_set(vh, p->z);
    mpz_set_si(ql, 1);
    mpz_set_si(qh, 1);
    mpz_set_si(tmp,0);

    for (j = mpz_sizeinbase(s,2)-1; j >= 1; j--) {
        /* ql = ql*qh (mod n) */
        mpz_mul(ql, ql, qh);
        mpz_mod(ql, ql, n->z);
        if (mpz_tstbit(s,j) == 1) {
            /* qh = ql*q */
            mpz_mul_si(qh, ql, q);

            /* uh = uh*vh (mod n) */
            mpz_mul(uh, uh, vh);
            mpz_mod(uh, uh, n->z);

            /* vl = vh*vl - p*ql (mod n) */
            mpz_mul(vl, vh, vl);
            mpz_mul(tmp, ql, p->z);
            mpz_sub(vl, vl, tmp);
            mpz_mod(vl, vl, n->z);

            /* vh = vh*vh - 2*qh (mod n) */
            mpz_mul(vh, vh, vh);
            mpz_mul_si(tmp, qh, 2);
            mpz_sub(vh, vh, tmp);
            mpz_mod(vh, vh, n->z);
        }
        else {
            /* qh = ql */
            mpz_set(qh, ql);

            /* uh = uh*vl - ql (mod n) */
            mpz_mul(uh, uh, vl);
            mpz_sub(uh, uh, ql);
            mpz_mod(uh, uh, n->z);

            /* vh = vh*vl - p*ql (mod n) */
            mpz_mul(vh, vh, vl);
            mpz_mul(tmp, ql, p->z);
            mpz_sub(vh, vh, tmp);
            mpz_mod(vh, vh, n->z);

            /* vl = vl*vl - 2*ql (mod n) */
            mpz_mul(vl, vl, vl);
            mpz_mul_si(tmp, ql, 2);
            mpz_sub(vl, vl, tmp);
            mpz_mod(vl, vl, n->z);
        }
    }
    /* ql = ql*qh */
    mpz_mul(ql, ql, qh);

    /* qh = ql*q */
    mpz_mul_si(qh, ql, q);

    /* uh = uh*vl - ql */
    mpz_mul(uh, uh, vl);
    mpz_sub(uh, uh, ql);

    /* vl = vh*vl - p*ql */
    mpz_mul(vl, vh, vl);
    mpz_mul(tmp, ql, p->z);
    mpz_sub(vl, vl, tmp);

    /* ql = ql*qh */
    mpz_mul(ql, ql, qh);

    mpz_mod(uh, uh, n->z);
    mpz_mod(vl, vl, n->z);

    /* tmp = n-2, for the following comparison */
    mpz_sub_ui(tmp, n->z, 2);

    /* uh contains LucasU_s and vl contains LucasV_s */
    if (  ((mpz_cmp_ui(uh, 0) == 0) 
            && 
           ((mpz_cmp(vl, tmp) == 0) || (mpz_cmp_si(vl, 2) == 0))
          )
          ||
          (mpz_cmp_ui(vl, 0) == 0)
       ) {
        result = Py_True;
        goto cleanup;
    }

    for (j = 1; j < r-1; j++) {
        /* vl = vl*vl - 2*ql (mod n) */
        mpz_mul(vl, vl, vl);
        mpz_mul_si(tmp, ql, 2);
        mpz_sub(vl, vl, tmp);
        mpz_mod(vl, vl, n->z);

        /* ql = ql*ql (mod n) */
        mpz_mul(ql, ql, ql);
        mpz_mod(ql, ql, n->z);

        if (mpz_cmp_ui(vl, 0) == 0) {
            result = Py_True;
            goto cleanup;
        }
    }

    result = Py_False;
  cleanup:
    Py_XINCREF(result);
    mpz_clear(zD);
    mpz_clear(s);
    mpz_clear(nmj);
    mpz_clear(nm2);
    mpz_clear(res);
    mpz_clear(uh);
    mpz_clear(vl);
    mpz_clear(vh);
    mpz_clear(ql);
    mpz_clear(qh);
    mpz_clear(tmp);
    Py_XDECREF((PyObject*)p);
    Py_XDECREF((PyObject*)n);
    return result;
}

/* ***********************************************************************************************
 * mpz_selfridge_prp:
 * A "Lucas-Selfridge probable prime" n is a "Lucas probable prime" using Selfridge parameters of:
 * Find the first element D in the sequence {5, -7, 9, -11, 13, ...} such that Jacobi(D,n) = -1
 * Then use P=1 and Q=(1-D)/4 in the Lucas probable prime test.
 * Make sure n is not a perfect square, otherwise the search for D will only stop when D=n.
 * ***********************************************************************************************/

PyDoc_STRVAR(doc_mpz_is_selfridge_prp,
"is_selfridge_prp($module, n, /)\n--\n\n"
"Return `True` if n is a Lucas probable prime with Selfidge parameters\n"
"(p,q). The Selfridge parameters are chosen by finding the first\n"
"element D in the sequence {5, -7, 9, -11, 13, ...} such that\n"
"Jacobi(D,n) == -1. Then let p=1 and q = (1-D)/4. Then perform\n"
"a Lucas probable prime test.");

static PyObject *
GMPY_mpz_is_selfridge_prp(PyObject *self, PyObject *args)
{
    MPZ_Object *n = NULL;
    PyObject *result = NULL, *temp = NULL;
    long d = 5, p = 1, q = 0, max_d = 1000000;
    int jacobi = 0;
    mpz_t zD;

    if (PyTuple_Size(args) != 1) {
        TYPE_ERROR("is_selfridge_prp() requires 1 integer argument");
        return NULL;
    }

    mpz_init(zD);

    n = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    if (!n) {
        TYPE_ERROR("is_selfridge_prp() requires 1 integer argument");
        goto cleanup;
    }

    /* Require n > 0. */
    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("is_selfridge_prp() requires 'n' be greater than 0");
        goto cleanup;
    }

    /* Check for n == 1 */
    if (mpz_cmp_ui(n->z, 1) == 0) {
        result = Py_False;
        goto cleanup;
    }

    /* Handle n even. */
    if (mpz_divisible_ui_p(n->z, 2)) {
        if (mpz_cmp_ui(n->z, 2) == 0)
            result = Py_True;
        else
            result = Py_False;
        goto cleanup;
    }

    mpz_set_ui(zD, d);

    while (1) {
        jacobi = mpz_jacobi(zD, n->z);

        /* if jacobi == 0, d is a factor of n, therefore n is composite... */
        /* if d == n, then either n is either prime or 9... */
        if (jacobi == 0) {
            if ((mpz_cmpabs(zD, n->z) == 0) && (mpz_cmp_ui(zD, 9) != 0)) {
                result = Py_True;
                goto cleanup;
            }
            else {
                result = Py_False;
                goto cleanup;
            }
        }
        if (jacobi == -1)
            break;

        /* if we get to the 5th d, make sure we aren't dealing with a square... */
        if (d == 13) {
            if (mpz_perfect_square_p(n->z)) {
                result = Py_False;
                goto cleanup;
            }
        }

        if (d < 0) {
            d *= -1;
            d += 2;
        }
        else {
            d += 2;
            d *= -1;
        }

        /* make sure we don't search forever */
        if (d >= max_d) {
            VALUE_ERROR("appropriate value for D cannot be found in is_selfridge_prp()");
            goto cleanup;
        }

        mpz_set_si(zD, d);
    }

    q = (1-d)/4;

    /* Since "O" is used, the refcount for n is incremented so deleting
     * temp will not delete n.
     */
    temp = Py_BuildValue("Oll", n, p, q);
    if (!temp)
        goto cleanup;
    result = GMPY_mpz_is_lucas_prp(NULL, temp);
    Py_DECREF(temp);
    goto return_result;

  cleanup:
    Py_XINCREF(result);
  return_result:
    mpz_clear(zD);
    Py_DECREF((PyObject*)n);
    return result;
}

/* *********************************************************************************************************
 * mpz_strongselfridge_prp:
 * A "strong Lucas-Selfridge probable prime" n is a "strong Lucas probable prime" using Selfridge parameters of:
 * Find the first element D in the sequence {5, -7, 9, -11, 13, ...} such that Jacobi(D,n) = -1
 * Then use P=1 and Q=(1-D)/4 in the strong Lucas probable prime test.
 * Make sure n is not a perfect square, otherwise the search for D will only stop when D=n.
 * **********************************************************************************************************/

PyDoc_STRVAR(doc_mpz_is_strongselfridge_prp,
"is_strong_selfridge_prp($module, n, /)\n--\n\n"
"Return `True` if n is a strong Lucas probable prime with Selfidge\n"
"parameters (p,q). The Selfridge parameters are chosen by finding\n"
"the first element D in the sequence {5, -7, 9, -11, 13, ...} such\n"
"that Jacobi(D,n) == -1. Then let p=1 and q = (1-D)/4. Then perform\n"
"a strong Lucas probable prime test.");

static PyObject *
GMPY_mpz_is_strongselfridge_prp(PyObject *self, PyObject *args)
{
    MPZ_Object *n = NULL;
    PyObject *result = NULL, *temp = NULL;
    long d = 5, p = 1, q = 0, max_d = 1000000;
    int jacobi = 0;
    mpz_t zD;

    if (PyTuple_Size(args) != 1) {
        TYPE_ERROR("is_strong_selfridge_prp() requires 1 integer argument");
        return NULL;
    }

    mpz_init(zD);

    n = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    if (!n) {
        TYPE_ERROR("is_strong_selfridge_prp() requires 1 integer argument");
        goto cleanup;
    }

    /* Require n > 0. */
    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("is_strong_selfridge_prp() requires 'n' be greater than 0");
        goto cleanup;
    }

    /* Check for n == 1 */
    if (mpz_cmp_ui(n->z, 1) == 0) {
        result = Py_False;
        goto cleanup;
    }

    /* Handle n even. */
    if (mpz_divisible_ui_p(n->z, 2)) {
        if (mpz_cmp_ui(n->z, 2) == 0)
            result = Py_True;
        else
            result = Py_False;
        goto cleanup;
    }


    mpz_set_ui(zD, d);

    while (1) {
        jacobi = mpz_jacobi(zD, n->z);

        /* if jacobi == 0, d is a factor of n, therefore n is composite... */
        /* if d == n, then either n is either prime or 9... */
        if (jacobi == 0) {
            if ((mpz_cmpabs(zD, n->z) == 0) && (mpz_cmp_ui(zD, 9) != 0)) {
                result = Py_True;
                goto cleanup;
            }
            else {
                result = Py_False;
                goto cleanup;
            }
        }
        if (jacobi == -1)
            break;

        /* if we get to the 5th d, make sure we aren't dealing with a square... */
        if (d == 13) {
            if (mpz_perfect_square_p(n->z)) {
                result = Py_False;
                goto cleanup;
            }
        }

        if (d < 0) {
            d *= -1;
            d += 2;
        }
        else {
            d += 2;
            d *= -1;
        }

        /* make sure we don't search forever */
        if (d >= max_d) {
            VALUE_ERROR("appropriate value for D cannot be found in is_strong_selfridge_prp()");
            goto cleanup;
        }

        mpz_set_si(zD, d);
    }

    q = (1-d)/4;

    /* Since "O" is used, the refcount for n is incremented so deleting
     * temp will not delete n.
     */
    temp = Py_BuildValue("Oll", n, p, q);
    if (!temp)
        goto cleanup;
    result = GMPY_mpz_is_stronglucas_prp(NULL, temp);
    Py_DECREF(temp);
    goto return_result;

  cleanup:
    Py_XINCREF(result);
  return_result:
    mpz_clear(zD);
    Py_DECREF((PyObject*)n);
    return result;
}

/* **********************************************************************************
 * mpz_bpsw_prp:
 * A "Baillie-Pomerance-Selfridge-Wagstaff probable prime" is a composite n such that
 * n is a strong probable prime to the base 2 and
 * n is a Lucas probable prime using the Selfridge parameters.
 * **********************************************************************************/

PyDoc_STRVAR(doc_mpz_is_bpsw_prp,
"is_bpsw_prp($module, n, /)\n--\n\n"
"Return `True` if n is a Baillie-Pomerance-Selfridge-Wagstaff probable\n"
"prime. A BPSW probable prime passes the `is_strong_prp()` test with base\n"
"2 and the `is_selfridge_prp()` test.\n");

static PyObject *
GMPY_mpz_is_bpsw_prp(PyObject *self, PyObject *args)
{
    MPZ_Object *n = NULL;
    PyObject *result = NULL, *temp = NULL;

    if (PyTuple_Size(args) != 1) {
        TYPE_ERROR("is_bpsw_prp() requires 1 integer argument");
        return NULL;
    }

    n = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    if (!n) {
        TYPE_ERROR("is_bpsw_prp() requires 1 integer argument");
        goto cleanup;
    }

    /* Require n > 0. */
    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("is_bpsw_prp() requires 'n' be greater than 0");
        goto cleanup;
    }

    /* Check for n == 1 */
    if (mpz_cmp_ui(n->z, 1) == 0) {
        result = Py_False;
        goto cleanup;
    }

    /* Handle n even. */
    if (mpz_divisible_ui_p(n->z, 2)) {
        if (mpz_cmp_ui(n->z, 2) == 0)
            result = Py_True;
        else
            result = Py_False;
        goto cleanup;
    }

    /* "O" is used to increment the reference to n so deleting temp won't
     * delete n.
     */
    temp = Py_BuildValue("Oi", n, 2);
    if (!temp)
        goto cleanup;
    result = GMPY_mpz_is_strong_prp(NULL, temp);
    Py_DECREF(temp);
    if (Py_IsFalse(result))
        goto return_result;
    /* Remember to ignore the preceding result */
    Py_DECREF(result);

    temp = PyTuple_Pack(1, n);
    if (!temp)
        goto cleanup;
    result = GMPY_mpz_is_selfridge_prp(NULL, temp);
    Py_DECREF(temp);
    goto return_result;

  cleanup:
    Py_XINCREF(result);
  return_result:
    Py_DECREF((PyObject*)n);
    return result;
}


/* ****************************************************************************************
 * mpz_strongbpsw_prp:
 * A "strong Baillie-Pomerance-Selfridge-Wagstaff probable prime" is a composite n such that
 * n is a strong probable prime to the base 2 and
 * n is a strong Lucas probable prime using the Selfridge parameters.
 * ****************************************************************************************/

PyDoc_STRVAR(doc_mpz_is_strongbpsw_prp,
"is_strong_bpsw_prp($module, n, /)\n--\n\n"
"Return `True` if n is a strong Baillie-Pomerance-Selfridge-Wagstaff\n"
"probable prime. A strong BPSW probable prime passes the `is_strong_prp()`\n"
"test with base and the `is_strong_selfridge_prp()` test.\n");

static PyObject *
GMPY_mpz_is_strongbpsw_prp(PyObject *self, PyObject *args)
{
    MPZ_Object *n = NULL;
    PyObject *result = NULL, *temp = NULL;

    if (PyTuple_Size(args) != 1) {
        TYPE_ERROR("is_strong_bpsw_prp() requires 1 integer argument");
        return NULL;
    }

    n = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    if (!n) {
        TYPE_ERROR("is_strong_bpsw_prp() requires 1 integer argument");
        goto cleanup;
    }

    /* Require n > 0. */
    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("is_strong_bpsw_prp() requires 'n' be greater than 0");
        goto cleanup;
    }

    /* Check for n == 1 */
    if (mpz_cmp_ui(n->z, 1) == 0) {
        result = Py_False;
        goto cleanup;
    }

    /* Handle n even. */
    if (mpz_divisible_ui_p(n->z, 2)) {
        if (mpz_cmp_ui(n->z, 2) == 0)
            result = Py_True;
        else
            result = Py_False;
        goto cleanup;
    }

    /* "O" is used to increment the reference to n so deleting temp won't
     * delete n.
     */
    temp = Py_BuildValue("Oi", n, 2);
    if (!temp)
        goto cleanup;
    result = GMPY_mpz_is_strong_prp(NULL, temp);
    Py_DECREF(temp);
    if (Py_IsFalse(result))
        goto return_result;
    /* Remember to ignore the preceding result */
    Py_DECREF(result);

    temp = PyTuple_Pack(1, n);
    if (!temp)
        goto cleanup;
    result = GMPY_mpz_is_strongselfridge_prp(NULL, temp);
    Py_DECREF(temp);
    goto return_result;

  cleanup:
    Py_XINCREF(result);
  return_result:
    Py_DECREF((PyObject*)n);
    return result;
}
