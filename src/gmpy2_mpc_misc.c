/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpc_misc.c                                                        *
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

PyDoc_STRVAR(GMPy_doc_context_phase,
"phase($self, x, /)\n--\n\n"
"Return the phase angle, also known as argument, of a complex x.");

PyDoc_STRVAR(GMPy_doc_function_phase,
"phase($module, x, /)\n--\n\n"
"Return the phase angle, also known as argument, of a complex x.");

static PyObject *
GMPy_Complex_Phase(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result = NULL;
    MPC_Object *tempx = NULL;

    CHECK_CONTEXT(context)

    if (!(IS_COMPLEX_ONLY(x))) {
        TYPE_ERROR("phase() argument type not supported");
        return NULL;
    }

    result = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPC_From_Complex(x, 1, 1, context);
    if (!result || !tempx) {
        Py_XDECREF(result);
        Py_XDECREF(tempx);
        return NULL;
    }

    mpfr_clear_flags();

    result->rc = mpc_arg(result->f, tempx->c, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Context_Phase(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Complex_Phase(other, context);
}

PyDoc_STRVAR(GMPy_doc_context_root_of_unity,
"root_of_unity($self, n, k, /)\n--\n\n"
"Return the n-th root of mpc(1) raised to the k-th power..");

PyDoc_STRVAR(GMPy_doc_function_root_of_unity,
"root_of_unity($module, n, k, /)\n--\n\n"
"Return the n-th root of mpc(1) raised to the k-th power..");

static PyObject *
GMPy_Complex_Root_Of_Unity(PyObject *n, PyObject *k, CTXT_Object *context)
{
    MPC_Object *result;
    unsigned long n_val, k_val;

    CHECK_CONTEXT(context)

    result = GMPy_MPC_New(0, 0, context);
    if (!result) {
        return NULL;
    }

    n_val = GMPy_Integer_AsUnsignedLong(n);
    k_val = GMPy_Integer_AsUnsignedLong(k);
    if ((n_val == (unsigned long)(-1) && PyErr_Occurred()) ||
        (k_val == (unsigned long)(-1) && PyErr_Occurred())) {
        Py_DECREF((PyObject*)result);
        VALUE_ERROR("root_of_unity() requires positive integer arguments.");
        return NULL;
    }

    result->rc = mpc_rootofunity(result->c, n_val, k_val, GET_MPC_ROUND(context));

    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Context_Root_Of_Unity(PyObject *self, PyObject *args)
{
    PyObject *n, *k;
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("root_of_unity() requires 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    n = PyTuple_GET_ITEM(args, 0);
    k = PyTuple_GET_ITEM(args, 1);

    if (IS_INTEGER(n) && IS_INTEGER(k)) {
        return GMPy_Complex_Root_Of_Unity(n, k, context);
    }
    else {
        TYPE_ERROR("root_of_unity() requires integer arguments");
        return NULL;
    }
}

PyDoc_STRVAR(GMPy_doc_context_norm,
"norm($self, x, /)\n--\n\n"
"Return the norm of a complex x. The norm(x) is defined as\n"
"x.real**2 + x.imag**2. abs(x) is the square root of norm(x).\n");

PyDoc_STRVAR(GMPy_doc_function_norm,
"norm($module, x, /)\n--\n\n"
"Return the norm of a complex x. The norm(x) is defined as\n"
"x.real**2 + x.imag**2. abs(x) is the square root of norm(x).\n");

static PyObject *
GMPy_Complex_Norm(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result = NULL;
    MPC_Object *tempx = NULL;

    CHECK_CONTEXT(context);

    if (!(IS_COMPLEX_ONLY(x))) {
        TYPE_ERROR("norm() argument type not supported");
        return NULL;
    }

    result = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPC_From_Complex(x, 1, 1, context);
    if (!result || !tempx) {
        Py_XDECREF(result);
        Py_XDECREF(tempx);
        return NULL;
    }

    mpfr_clear_flags();

    result->rc = mpc_norm(result->f, tempx->c, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Context_Norm(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Complex_Norm(other, context);
}

PyDoc_STRVAR(GMPy_doc_context_polar,
"polar($self, x, /)\n--\n\n"
"Return the polar coordinate form of a complex x that is in\n"
"rectangular form.");

PyDoc_STRVAR(GMPy_doc_function_polar,
"polar($module, x, /)\n--\n\n"
"Return the polar coordinate form of a complex x that is in\n"
"rectangular form.");

static PyObject *
GMPy_Complex_Polar(PyObject *x, CTXT_Object *context)
{
    PyObject *tempx, *abs, *phase, *result;

    CHECK_CONTEXT(context);

    if (!(IS_COMPLEX_ONLY(x))) {
        TYPE_ERROR("polar() argument type not supported");
        return NULL;
    }

    if (!(tempx = (PyObject*)GMPy_MPC_From_Complex(x, 1, 1, context))) {
        return NULL;
    }

    abs = GMPy_Complex_AbsWithType(tempx, OBJ_TYPE_MPC, context);
    phase = GMPy_Complex_Phase(tempx, context);
    Py_DECREF(tempx);
    result = PyTuple_New(2);
    if (!abs || !phase || !result) {
        Py_XDECREF(abs);
        Py_XDECREF(phase);
        Py_XDECREF(result);
        return NULL;
    }

    PyTuple_SET_ITEM(result, 0, abs);
    PyTuple_SET_ITEM(result, 1, phase);
    return result;
}

static PyObject *
GMPy_Context_Polar(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Complex_Polar(other, context);
}

PyDoc_STRVAR(GMPy_doc_context_rect,
"rect($self, r, phi, /)\n--\n\n"
"Return the rectangular coordinate form of a complex number that is\n"
"given in polar form.");

PyDoc_STRVAR(GMPy_doc_function_rect,
"rect($module, r, phi, /)\n--\n\n"
"Return the rectangular coordinate form of a complex number that is\n"
"given in polar form.");

static PyObject *
GMPy_Complex_Rect(PyObject *r, PyObject *phi, CTXT_Object *context)
{
    MPFR_Object *tempx, *tempy;
    MPC_Object *result;

    CHECK_CONTEXT(context);

    tempx = GMPy_MPFR_From_Real(r, 1, context);
    tempy = GMPy_MPFR_From_Real(phi, 1, context);
    result = GMPy_MPC_New(0, 0, context);
    if (!tempx || !tempy || !result) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)tempy);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    mpfr_cos(mpc_realref(result->c), tempy->f, GET_REAL_ROUND(context));
    mpfr_mul(mpc_realref(result->c), mpc_realref(result->c), tempx->f, GET_REAL_ROUND(context));
    mpfr_sin(mpc_imagref(result->c), tempy->f, GET_IMAG_ROUND(context));
    mpfr_mul(mpc_imagref(result->c), mpc_imagref(result->c), tempx->f, GET_IMAG_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);

    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Context_Rect(PyObject *self, PyObject *args)
{
    PyObject *r, *phi;
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("rect() requires 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    r = PyTuple_GET_ITEM(args, 0);
    phi = PyTuple_GET_ITEM(args, 1);

    if (IS_REAL(r) && IS_REAL(phi)) {
        return GMPy_Complex_Rect(r, phi, context);
    }
    else {
        TYPE_ERROR("rect() argument type not supported");
        return NULL;
    }
}

PyDoc_STRVAR(GMPy_doc_context_proj,
"proj($self, x, /)\n--\n\n"
"Returns the projection of a complex x on to the Riemann sphere.");

PyDoc_STRVAR(GMPy_doc_function_proj,
"proj($module, x, /)\n--\n\n"
"Returns the projection of a complex x on to the Riemann sphere.");

static PyObject *
GMPy_Complex_Proj(PyObject *x, CTXT_Object *context)
{
    MPC_Object *result = NULL, *tempx = NULL;

    CHECK_CONTEXT(context);

    if (!(IS_COMPLEX_ONLY(x))) {
        TYPE_ERROR("proj() argument type not supported");
        return NULL;
    }

    result = GMPy_MPC_New(0, 0, context);
    tempx = GMPy_MPC_From_Complex(x, 1, 1, context);
    if (!result || !tempx) {
        Py_XDECREF(result);
        Py_XDECREF(tempx);
        return NULL;
    }

    mpfr_clear_flags();

    result->rc = mpc_proj(result->c, tempx->c, GET_MPC_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Context_Proj(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Complex_Proj(other, context);
}
/* Implement the conjugate() method. */

PyDoc_STRVAR(GMPy_doc_mpc_conjugate_method,
"conjugate($self)\n--\n\n"
"Returns the conjugate of self.");

static PyObject *
GMPy_MPC_Conjugate_Method(PyObject *self, PyObject *args)
{
    MPC_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        return NULL;
    }

    result->rc = mpc_conj(result->c, MPC(self), GET_MPC_ROUND(context));

    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

/* Implement the .precision attribute of an mpc. */

static PyObject *
GMPy_MPC_GetPrec_Attrib(MPC_Object *self, void *closure)
{
    mpfr_prec_t rprec = 0, iprec = 0;

    mpc_get_prec2(&rprec, &iprec, self->c);
    return Py_BuildValue("(nn)", rprec, iprec);
}

/* Implement the .rc attribute of an mpc. */

static PyObject *
GMPy_MPC_GetRc_Attrib(MPC_Object *self, void *closure)
{
    return Py_BuildValue("(ii)", MPC_INEX_RE(self->rc), MPC_INEX_IM(self->rc));
}

/* Implement the .imag attribute of an mpc. */

static PyObject *
GMPy_MPC_GetImag_Attrib(MPC_Object *self, void *closure)
{
    MPFR_Object *result = NULL;
    CTXT_Object *context = NULL;
    mpfr_prec_t rprec = 0, iprec = 0;
    mpc_get_prec2(&rprec, &iprec, self->c);

    CHECK_CONTEXT(context);

    if ((result = GMPy_MPFR_New(iprec, context))) {
        result->rc = mpc_imag(result->f, self->c, GET_MPFR_ROUND(context));
        _GMPy_MPFR_Cleanup(&result, context);
    }
    return (PyObject*)result;
}

/* Implement the .real attribute of an mpc. */

static PyObject *
GMPy_MPC_GetReal_Attrib(MPC_Object *self, void *closure)
{
    MPFR_Object *result = NULL;
    CTXT_Object *context = NULL;
    mpfr_prec_t rprec = 0, iprec = 0;
    mpc_get_prec2(&rprec, &iprec, self->c);

    CHECK_CONTEXT(context);

    if ((result = GMPy_MPFR_New(rprec, context))) {
        result->rc = mpc_real(result->f, self->c, context->ctx.mpfr_round);
        _GMPy_MPFR_Cleanup(&result, context);
    }
    return (PyObject*)result;
}

/* Implement the nb_bool slot. */

static int
GMPy_MPC_NonZero_Slot(MPC_Object *self)
{
    return !MPC_IS_ZERO_P(self->c);
}

PyDoc_STRVAR(GMPy_doc_mpc_sizeof_method,
"__sizeof__($self)\n--\n\n"
"Returns the amount of memory consumed by self.");

static PyObject *
GMPy_MPC_SizeOf_Method(PyObject *self, PyObject *other)
{
    return PyLong_FromSize_t(sizeof(MPC_Object) + \
        (((mpc_realref(MPC(self))->_mpfr_prec + mp_bits_per_limb - 1) / \
        mp_bits_per_limb) * sizeof(mp_limb_t)) + \
        (((mpc_imagref(MPC(self))->_mpfr_prec + mp_bits_per_limb - 1) / \
        mp_bits_per_limb) * sizeof(mp_limb_t)));
}
