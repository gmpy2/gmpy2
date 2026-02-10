/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_const.c                                                           *
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

PyDoc_STRVAR(GMPy_doc_function_const_pi,
"const_pi($module, precision=0)\n--\n\n"
"Return the constant pi using the specified precision. If no\n"
"precision is specified, the default precision is used.");

static PyObject *
GMPy_Function_Const_Pi(PyObject *self, PyObject *args, PyObject *keywds)
{
    MPFR_Object *result = NULL;
    mpfr_prec_t bits = 0;
    static char *kwlist[] = {"precision", NULL};
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "|l", kwlist, &bits)) {
        return NULL;
    }

    if ((result = GMPy_MPFR_New(bits, context))) {
        mpfr_clear_flags();
        result->rc = mpfr_const_pi(result->f, GET_MPFR_ROUND(context));
        _GMPy_MPFR_Cleanup(&result, context);
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_context_const_pi,
"const_pi($self)\n--\n\n"
"Return the constant pi using the context's precision.");

static PyObject *
GMPy_Context_Const_Pi(PyObject *self, PyObject *args)
{
    MPFR_Object *result = NULL;
    CTXT_Object *context = (CTXT_Object*)self;

    if ((result = GMPy_MPFR_New(0, context))) {
        mpfr_clear_flags();
        result->rc = mpfr_const_pi(result->f, GET_MPFR_ROUND(context));
        _GMPy_MPFR_Cleanup(&result, context);
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_function_const_euler,
"const_euler($module, precision=0)\n--\n\n"
"Return the euler constant using the specified precision. If no\n"
"precision is specified, the default precision is used.");

static PyObject *
GMPy_Function_Const_Euler(PyObject *self, PyObject *args, PyObject *keywds)
{
    MPFR_Object *result = NULL;
    mpfr_prec_t bits = 0;
    static char *kwlist[] = {"precision", NULL};
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "|l", kwlist, &bits)) {
        return NULL;
    }

    if ((result = GMPy_MPFR_New(bits, context))) {
        mpfr_clear_flags();
        result->rc = mpfr_const_euler(result->f, GET_MPFR_ROUND(context));
        _GMPy_MPFR_Cleanup(&result, context);
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_context_const_euler,
"const_euler($self)\n--\n\n"
"Return the euler constant using the context's precision.");

static PyObject *
GMPy_Context_Const_Euler(PyObject *self, PyObject *args)
{
    MPFR_Object *result = NULL;
    CTXT_Object *context = (CTXT_Object*)self;

    if ((result = GMPy_MPFR_New(0, context))) {
        mpfr_clear_flags();
        result->rc = mpfr_const_euler(result->f, GET_MPFR_ROUND(context));
        _GMPy_MPFR_Cleanup(&result, context);
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_function_const_log2,
"const_log2($module, precision=0)\n--\n\n"
"Return the log2 constant using the specified precision. If no\n"
"precision is specified, the default precision is used.");

static PyObject *
GMPy_Function_Const_Log2(PyObject *self, PyObject *args, PyObject *keywds)
{
    MPFR_Object *result = NULL;
    mpfr_prec_t bits = 0;
    static char *kwlist[] = {"precision", NULL};
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "|l", kwlist, &bits)) {
        return NULL;
    }

    if ((result = GMPy_MPFR_New(bits, context))) {
        mpfr_clear_flags();
        result->rc = mpfr_const_log2(result->f, GET_MPFR_ROUND(context));
        _GMPy_MPFR_Cleanup(&result, context);
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_context_const_log2,
"const_log2($self)\n--\n\n"
"Return the log2 constant using the context's precision.");

static PyObject *
GMPy_Context_Const_Log2(PyObject *self, PyObject *args)
{
    MPFR_Object *result = NULL;
    CTXT_Object *context = (CTXT_Object*)self;

    if ((result = GMPy_MPFR_New(0, context))) {
        mpfr_clear_flags();
        result->rc = mpfr_const_log2(result->f, GET_MPFR_ROUND(context));
        _GMPy_MPFR_Cleanup(&result, context);
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_function_const_catalan,
"const_catalan($module, precision=0)\n--\n\n"
"Return the catalan constant using the specified precision. If no\n"
"precision is specified, the default precision is used.");

static PyObject *
GMPy_Function_Const_Catalan(PyObject *self, PyObject *args, PyObject *keywds)
{
    MPFR_Object *result = NULL;
    mpfr_prec_t bits = 0;
    static char *kwlist[] = {"precision", NULL};
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "|l", kwlist, &bits)) {
        return NULL;
    }

    if ((result = GMPy_MPFR_New(bits, context))) {
        mpfr_clear_flags();
        result->rc = mpfr_const_catalan(result->f, GET_MPFR_ROUND(context));
        _GMPy_MPFR_Cleanup(&result, context);
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_context_const_catalan,
"const_catalan($self)\n--\n\n"
"Return the catalan constant using the context's precision.");

static PyObject *
GMPy_Context_Const_Catalan(PyObject *self, PyObject *args)
{
    MPFR_Object *result = NULL;
    CTXT_Object *context = (CTXT_Object*)self;

    if ((result = GMPy_MPFR_New(0, context))) {
        mpfr_clear_flags();
        result->rc = mpfr_const_catalan(result->f, GET_MPFR_ROUND(context));
        _GMPy_MPFR_Cleanup(&result, context);
    }

    return (PyObject*)result;
}
