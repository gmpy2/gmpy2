/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz_lucas.c                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP, MPFR, and MPC multiple precision           *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2011 David Cleaver                                            *
 *                                                                         *
 * Copyright 2012 - 2025 Case Van Horsen                                   *
 *                                                                         *
 * The original file is available at:                                      *
 *   <http://sourceforge.net/projects/mpzlucas/files/>                     *
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

PyDoc_STRVAR(doc_mpz_lucasu,
"lucasu($module, p, q, k, /)\n--\n\n"
"Return the k-th element of the Lucas U sequence defined by p,q.\n"
"p*p - 4*q must not equal 0; k must be greater than or equal to 0.");

static PyObject *
GMPY_mpz_lucasu(PyObject *self, PyObject *args)
{
    /* Adaptation of algorithm found in
     * https://web.archive.org/web/20150202074230/http://joye.site88.net/papers/JQ96lucas.pdf
     * calculate u[k] of Lucas U sequence for p,q.
     * Note: p^2-4q=0 is not tested, not a proper Lucas sequence!!
     */

    MPZ_Object *result = NULL, *p = NULL, *q = NULL, *k = NULL;
    size_t s = 0, j = 0;
    mpz_t uh, vl, vh, ql, qh, tmp;

    if (PyTuple_Size(args) != 3) {
        TYPE_ERROR("lucasu() requires 3 integer arguments");
        return NULL;
    }

    mpz_init(uh);
    mpz_init(vl);
    mpz_init(vh);
    mpz_init(ql);
    mpz_init(qh);
    mpz_init(tmp);

    p = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    q = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    k = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 2), NULL);
    if (!p || !q || !k) {
        TYPE_ERROR("lucasu() requires 3 integer arguments");
        goto cleanup;
    }

    /* Check if p*p - 4*q == 0. */

    mpz_mul(tmp, p->z, p->z);
    mpz_mul_ui(qh, q->z, 4);
    mpz_sub(tmp, tmp, qh);
    if (mpz_sgn(tmp) == 0) {
        VALUE_ERROR("invalid values for p,q in lucasu()");
        goto cleanup;
    }

    /* Check if k < 0. */

    if (mpz_sgn(k->z) < 0) {
        VALUE_ERROR("invalid value for k in lucasu()");
        goto cleanup;
    }

    mpz_set_si(uh, 1);
    mpz_set_si(vl, 2);
    mpz_set(vh, p->z);
    mpz_set_si(ql, 1);
    mpz_set_si(qh, 1);
    mpz_set_si(tmp, 0);

    if (mpz_sgn(k->z) == 0) {
        mpz_set_si(uh, 0);
        goto end;
    }

    s = mpz_scan1(k->z, 0);
    for (j = mpz_sizeinbase(k->z,2)-1; j >= s+1; j--) {
        /* ql = ql*qh */
        mpz_mul(ql, ql, qh);
        if (mpz_tstbit(k->z,j) == 1) {
            /* qh = ql*q */
            mpz_mul(qh, ql, q->z);

            /* uh = uh*vh */
            mpz_mul(uh, uh, vh);

            /* vl = vh*vl - p*ql */
            mpz_mul(vl, vh, vl);
            mpz_mul(tmp, ql, p->z);
            mpz_sub(vl, vl, tmp);

            /* vh = vh*vh - 2*qh */
            mpz_mul(vh, vh, vh);
            mpz_mul_si(tmp, qh, 2);
            mpz_sub(vh, vh, tmp);
        }
        else {
            /* qh = ql */
            mpz_set(qh, ql);

            /* uh = uh*vl - ql */
            mpz_mul(uh, uh, vl);
            mpz_sub(uh, uh, ql);

            /* vh = vh*vl - p*ql */
            mpz_mul(vh, vh, vl);
            mpz_mul(tmp, ql, p->z);
            mpz_sub(vh, vh, tmp);

            /* vl = vl*vl - 2*ql */
            mpz_mul(vl, vl, vl);
            mpz_mul_si(tmp, ql, 2);
            mpz_sub(vl, vl, tmp);
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
        /* uh = uh*vl */
        mpz_mul(uh, uh, vl);

        /* vl = vl*vl - 2*ql */
        mpz_mul(vl, vl, vl);
        mpz_mul_si(tmp, ql, 2);
        mpz_sub(vl, vl, tmp);

        /* ql = ql*ql */
        mpz_mul(ql, ql, ql);
    }

  end:
    if (!(result = GMPy_MPZ_New(NULL)))
        goto cleanup;

    /* uh contains our return value */
    mpz_set(result->z, uh);

  cleanup:
    mpz_clear(uh);
    mpz_clear(vl);
    mpz_clear(vh);
    mpz_clear(ql);
    mpz_clear(qh);
    mpz_clear(tmp);
    Py_XDECREF((PyObject*)p);
    Py_XDECREF((PyObject*)q);
    Py_XDECREF((PyObject*)k);

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpz_lucasu_mod,
"lucasu_mod($module, p, q, k, n, /)\n--\n\n"
"Return the k-th element of the Lucas U sequence defined by p,q (mod n).\n"
"p*p - 4*q must not equal 0; k must be greater than or equal to 0;\n"
"n must be greater than 0.");

static PyObject *
GMPY_mpz_lucasu_mod(PyObject *self, PyObject *args)
{
    /* Adaptation of algorithm found in http://joye.site88.net/papers/JQ96lucas.pdf
     * calculate u[k] (modulo n) of Lucas U sequence for p,q.
     * Note: p^2-4q=0 is not tested, not a proper Lucas sequence!!
     */

    MPZ_Object *result = NULL, *p = NULL, *q = NULL, *k = NULL, *n = NULL;

    size_t s = 0, j = 0;
    mpz_t uh, vl, vh, ql, qh, tmp;

    if (PyTuple_Size(args) != 4) {
        TYPE_ERROR("lucasu_mod() requires 4 integer arguments");
        return NULL;
    }

    mpz_init(uh);
    mpz_init(vl);
    mpz_init(vh);
    mpz_init(ql);
    mpz_init(qh);
    mpz_init(tmp);

    p = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    q = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    k = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 2), NULL);
    n = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 3), NULL);
    if (!p || !q || !k || !n) {
        TYPE_ERROR("lucasu_mod() requires 4 integer arguments");
        goto cleanup;
    }

    /* Check if p*p - 4*q == 0. */

    mpz_mul(tmp, p->z, p->z);
    mpz_mul_ui(qh, q->z, 4);
    mpz_sub(tmp, tmp, qh);
    if (mpz_sgn(tmp) == 0) {
        VALUE_ERROR("invalid values for p,q in lucasu_mod()");
        goto cleanup;
    }

    /* Check if k < 0. */

    if (mpz_sgn(k->z) < 0) {
        VALUE_ERROR("invalid value for k in lucasu_mod()");
        goto cleanup;
    }

    /* Check if n > 0. */

    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("invalid value for n in lucasu_mod()");
        goto cleanup;
    }

    mpz_set_si(uh, 1);
    mpz_set_si(vl, 2);
    mpz_set(vh, p->z);
    mpz_set_si(ql, 1);
    mpz_set_si(qh, 1);
    mpz_set_si(tmp, 0);

    if (mpz_sgn(k->z) == 0) {
        mpz_set_si(uh, 0);
        goto end;
    }

    s = mpz_scan1(k->z, 0);
    for (j = mpz_sizeinbase(k->z,2)-1; j >= s+1; j--) {
        /* ql = ql*qh (mod n) */
        mpz_mul(ql, ql, qh);
        mpz_mod(ql, ql, n->z);
        if (mpz_tstbit(k->z,j) == 1) {
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

  end:
    if (!(result = GMPy_MPZ_New(NULL)))
        goto cleanup;

    /* uh contains our return value */
    mpz_mod(result->z, uh, n->z);

  cleanup:
    mpz_clear(uh);
    mpz_clear(vl);
    mpz_clear(vh);
    mpz_clear(ql);
    mpz_clear(qh);
    mpz_clear(tmp);
    Py_XDECREF((PyObject*)p);
    Py_XDECREF((PyObject*)q);
    Py_XDECREF((PyObject*)k);
    Py_XDECREF((PyObject*)n);

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpz_lucasv,
"lucasv($module, p, q, k, /)\n--\n\n"
"Return the k-th element of the Lucas V sequence defined by p,q.\n"
"p*p - 4*q must not equal 0; k must be greater than or equal to 0.");

static PyObject *
GMPY_mpz_lucasv(PyObject *self, PyObject *args)
{
    /* Adaptation of algorithm found in http://joye.site88.net/papers/JQ96lucas.pdf
     * calculate v[k] of Lucas V sequence for p,q.
     * Note: p^2-4q=0 is not tested, not a proper Lucas sequence!!
     */

    MPZ_Object *result = NULL, *p = NULL, *q = NULL, *k = NULL;
    size_t s = 0, j = 0;
    mpz_t vl, vh, ql, qh, tmp;

    if (PyTuple_Size(args) != 3) {
        TYPE_ERROR("lucasv() requires 3 integer arguments");
        return NULL;
    }

    mpz_init(vl);
    mpz_init(vh);
    mpz_init(ql);
    mpz_init(qh);
    mpz_init(tmp);

    p = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    q = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    k = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 2), NULL);
    if (!p || !q || !k) {
        TYPE_ERROR("lucasv() requires 3 integer arguments");
        goto cleanup;
    }

    /* Check if p*p - 4*q == 0. */

    mpz_mul(tmp, p->z, p->z);
    mpz_mul_ui(qh, q->z, 4);
    mpz_sub(tmp, tmp, qh);
    if (mpz_sgn(tmp) == 0) {
        VALUE_ERROR("invalid values for p,q in lucasv()");
        goto cleanup;
    }

    /* Check if k < 0. */

    if (mpz_sgn(k->z) < 0) {
        VALUE_ERROR("invalid value for k in lucasv()");
        goto cleanup;
    }

    mpz_set_si(vl, 2);
    mpz_set(vh, p->z);
    mpz_set_si(ql, 1);
    mpz_set_si(qh, 1);
    mpz_set_si(tmp,0);

    if (mpz_sgn(k->z) == 0) {
        goto end;
    }

    s = mpz_scan1(k->z, 0);
    for (j = mpz_sizeinbase(k->z,2)-1; j >= s+1; j--) {
        /* ql = ql*qh */
        mpz_mul(ql, ql, qh);
        if (mpz_tstbit(k->z,j) == 1) {
            /* qh = ql*q */
            mpz_mul(qh, ql, q->z);

            /* vl = vh*vl - p*ql */
            mpz_mul(vl, vh, vl);
            mpz_mul(tmp, ql, p->z);
            mpz_sub(vl, vl, tmp);

            /* vh = vh*vh - 2*qh */
            mpz_mul(vh, vh, vh);
            mpz_mul_si(tmp, qh, 2);
            mpz_sub(vh, vh, tmp);
        }
        else {
            /* qh = ql */
            mpz_set(qh, ql);

            /* vh = vh*vl - p*ql */
            mpz_mul(vh, vh, vl);
            mpz_mul(tmp, ql, p->z);
            mpz_sub(vh, vh, tmp);

            /* vl = vl*vl - 2*ql */
            mpz_mul(vl, vl, vl);
            mpz_mul_si(tmp, ql, 2);
            mpz_sub(vl, vl, tmp);
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
        /* vl = vl*vl - 2*ql */
        mpz_mul(vl, vl, vl);
        mpz_mul_si(tmp, ql, 2);
        mpz_sub(vl, vl, tmp);

        /* ql = ql*ql */
        mpz_mul(ql, ql, ql);
    }

  end:
    if (!(result = GMPy_MPZ_New(NULL)))
        goto cleanup;

    /* vl contains our return value */
    mpz_set(result->z, vl);

  cleanup:
    mpz_clear(vl);
    mpz_clear(vh);
    mpz_clear(ql);
    mpz_clear(qh);
    mpz_clear(tmp);
    Py_XDECREF((PyObject*)p);
    Py_XDECREF((PyObject*)q);
    Py_XDECREF((PyObject*)k);

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpz_lucasv_mod,
"lucasv_mod($module, p, q, k, n, /)\n--\n\n"
"Return the k-th element of the Lucas V sequence defined by p,q (mod n).\n"
"p*p - 4*q must not equal 0; k must be greater than or equal to 0;\n"
"n must be greater than 0.");

static PyObject *
GMPY_mpz_lucasv_mod(PyObject *self, PyObject *args)
{
    /* Adaptation of algorithm found in http://joye.site88.net/papers/JQ96lucas.pdf
     * calculate v[k] (modulo n) of Lucas V sequence for p,q.
     * Note: p^2-4q=0 is not tested, not a proper Lucas sequence!!
     */

    MPZ_Object *result = NULL, *p = NULL, *q = NULL, *k = NULL, *n = NULL;

    size_t s = 0, j = 0;
    mpz_t vl, vh, ql, qh, tmp;

    if (PyTuple_Size(args) != 4) {
        TYPE_ERROR("lucasv_mod() requires 4 integer arguments");
        return NULL;
    }

    mpz_init(vl);
    mpz_init(vh);
    mpz_init(ql);
    mpz_init(qh);
    mpz_init(tmp);

    p = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    q = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    k = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 2), NULL);
    n = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 3), NULL);
    if (!p || !q || !k || !n) {
        TYPE_ERROR("lucasv_mod() requires 4 integer arguments");
        goto cleanup;
    }

    /* Check if p*p - 4*q == 0. */

    mpz_mul(tmp, p->z, p->z);
    mpz_mul_ui(qh, q->z, 4);
    mpz_sub(tmp, tmp, qh);
    if (mpz_sgn(tmp) == 0) {
        VALUE_ERROR("invalid values for p,q in lucasv_mod()");
        goto cleanup;
    }

    /* Check if k < 0. */

    if (mpz_sgn(k->z) < 0) {
        VALUE_ERROR("invalid value for k in lucasv_mod()");
        goto cleanup;
    }

    /* Check if n > 0. */

    if (mpz_sgn(n->z) <= 0) {
        VALUE_ERROR("invalid value for n in lucasv_mod()");
        goto cleanup;
    }

    mpz_set_si(vl, 2);
    mpz_set(vh, p->z);
    mpz_set_si(ql, 1);
    mpz_set_si(qh, 1);
    mpz_set_si(tmp,0);

    if (mpz_sgn(k->z) == 0) {
        goto end;
    }

    s = mpz_scan1(k->z, 0);
    for (j = mpz_sizeinbase(k->z,2)-1; j >= s+1; j--) {
        /* ql = ql*qh (mod n) */
        mpz_mul(ql, ql, qh);
        mpz_mod(ql, ql, n->z);
        if (mpz_tstbit(k->z,j) == 1) {
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

  end:
    if (!(result = GMPy_MPZ_New(NULL)))
        goto cleanup;

    /* vl contains our return value */
    mpz_mod(result->z, vl, n->z);

  cleanup:
    mpz_clear(vl);
    mpz_clear(vh);
    mpz_clear(ql);
    mpz_clear(qh);
    mpz_clear(tmp);
    Py_XDECREF((PyObject*)p);
    Py_XDECREF((PyObject*)q);
    Py_XDECREF((PyObject*)k);
    Py_XDECREF((PyObject*)n);

    return (PyObject*)result;
}
