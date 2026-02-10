/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_xmpz_misc.c                                                       *
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

PyDoc_STRVAR(GMPy_doc_xmpz_function_xbit_mask,
"xbit_mask($module, n, /)\n--\n\n"
"Return an 'xmpz' exactly n bits in length with all bits set.\n");

static PyObject *
GMPy_XMPZ_Function_XbitMask(PyObject *self, PyObject *other)
{
    Py_ssize_t i = 0;
    XMPZ_Object* result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    i = GMPy_Integer_AsSsize_t(other);
    if (i == -1 && PyErr_Occurred()) {
        TYPE_ERROR("xbit_mask() requires 'int' argument");
        return NULL;
    }

    if (i < 0) {
        VALUE_ERROR("mask length must be >= 0");
        return NULL;
    }

    if (!(result = GMPy_XMPZ_New(context))) {
        return NULL;
    }

    mpz_set_ui(result->z, 1);
    mpz_mul_2exp(result->z, result->z, i);
    mpz_sub_ui(result->z, result->z, 1);

    return (PyObject*)result;
}

static PyObject *
GMPy_XMPZ_Abs_Slot(XMPZ_Object *x)
{
    mpz_abs(x->z, x->z);
    Py_RETURN_NONE;
}

static PyObject *
GMPy_XMPZ_Neg_Slot(XMPZ_Object *x)
{
    mpz_neg(x->z, x->z);
    Py_RETURN_NONE;
}

static PyObject *
GMPy_XMPZ_Pos_Slot(XMPZ_Object *x)
{
    Py_RETURN_NONE;
}

static int
GMPy_XMPZ_NonZero_Slot(XMPZ_Object *x)
{
    return mpz_sgn(x->z) != 0;
}

/* BIT OPERATIONS */

static PyObject *
GMPy_XMPZ_Com_Slot(XMPZ_Object *x)
{
    mpz_com(x->z, x->z);
    Py_RETURN_NONE;
}

PyDoc_STRVAR(GMPy_doc_xmpz_method_make_mpz,
"make_mpz($self)\n--\n\n"
"Return an `mpz` by converting self as quickly as possible.\n\n"
"NOTE: Optimized for speed so the original `xmpz` value is set to 0!");

static PyObject *
GMPy_XMPZ_Method_MakeMPZ(PyObject *self, PyObject *other)
{
    MPZ_Object* result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPZ_New(context))) {
        return NULL;
    }
    mpz_swap(result->z, MPZ(self));
    mpz_set_ui(MPZ(self), 0);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_xmpz_method_copy,
"copy($self)\n--\n\n"
"Return a copy of a self.");

static PyObject *
GMPy_XMPZ_Method_Copy(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    return (PyObject*)GMPy_XMPZ_From_XMPZ((XMPZ_Object*)self, context);
}

/*
 * Add mapping support to xmpz objects.
 */

static Py_ssize_t
GMPy_XMPZ_Method_Length(XMPZ_Object *obj)
{
    return mpz_sizeinbase(obj->z, 2);
}

static PyObject *
GMPy_XMPZ_Method_SubScript(XMPZ_Object* self, PyObject* item)
{
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyIndex_Check(item)) {
        Py_ssize_t i;

        i = PyLong_AsSsize_t(item);
        if (i == -1 && PyErr_Occurred()) {
            INDEX_ERROR("argument too large to be converted to an index");
            return NULL;
        }
        if (i < 0) {
            i += mpz_sizeinbase(self->z, 2);
        }
        return PyLong_FromLong(mpz_tstbit(self->z, i));
    }
    else if (PySlice_Check(item)) {
        Py_ssize_t start, stop, step, slicelength, cur, i;
        MPZ_Object *result;

        if (PySlice_GetIndicesEx(item,
                         mpz_sizeinbase(self->z, 2),
                         &start, &stop, &step, &slicelength) < 0) {
            return NULL;
        }

        if ((step < 0 && start < stop) || (step > 0 && start > stop)) {
            stop = start;
        }

        if (!(result = GMPy_MPZ_New(context))) {
            return NULL;
        }

        mpz_set_ui(result->z, 0);
        if (slicelength > 0) {
            for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                if (mpz_tstbit(self->z, cur)) {
                    mpz_setbit(result->z, i);
                }
            }
        }
        return (PyObject*)result;
    }
    else {
        TYPE_ERROR("bit positions must be integers");
        return NULL;
    }
}

static int
GMPy_XMPZ_Method_AssignSubScript(XMPZ_Object* self, PyObject* item, PyObject* value)
{
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_M1(context);

    if (PyIndex_Check(item)) {
        Py_ssize_t bit_value, i;

        i = PyLong_AsSsize_t(item);
        if (i == -1 && PyErr_Occurred()) {
            INDEX_ERROR("argument too large to be converted to an index");
            return -1;
        }
        if (i < 0) {
            i += mpz_sizeinbase(self->z, 2);
        }

        bit_value = PyLong_AsSsize_t(value);
        if (bit_value == -1 && PyErr_Occurred()) {
            VALUE_ERROR("bit value must be 0 or 1");
            return -1;
        }
        if (bit_value == 1) {
            mpz_setbit(self->z, i);
            return 0;
        }
        else if (bit_value == 0) {
            mpz_clrbit(self->z, i);
            return 0;
        }
        else {
            VALUE_ERROR("bit value must be 0 or 1");
            return -1;
        }
    }
    else if (PySlice_Check(item)) {
        Py_ssize_t cur, i, seq_len, start, stop, step, slicelength, temp;

        seq_len = mpz_sizeinbase(self->z, 2);
        if (!Py_IsNone(((PySliceObject*)item)->stop)) {
            /* If a fixed endpoint is specified, and the endpoint is greater
             * than the length of the xmpz object, allow the underlying xmpz
             * object to be made larger.
             */
            temp = PyLong_AsSsize_t(((PySliceObject*)item)->stop);
            if (temp == -1  && PyErr_Occurred()) {
                return 0;
            }
            if (temp > seq_len) {
                seq_len = temp;
            }
        }

        if (PySlice_GetIndicesEx(item,
                        seq_len,
                        &start, &stop, &step, &slicelength) < 0) {
            return -1;
        }

        if ((step < 0 && start < stop) || (step > 0 && start > stop)) {
            stop = start;
        }

        if (value == NULL) {
            TYPE_ERROR("deleting bits not supported");
            return -1;
        }

        else {
            int bit;
            MPZ_Object *tempx;

            if (!(tempx = GMPy_MPZ_From_Integer(value, context))) {
                VALUE_ERROR("must specify bit sequence as an integer");
                return -1;
            }
            if (mpz_sgn(tempx->z) == 0) {
                for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                    mpz_clrbit(self->z, cur);
                }
            }
            else if (!(mpz_cmp_si(tempx->z, -1))) {
                for (cur = start + (slicelength-1) * step, i = 0;
                     i < slicelength;
                     cur -= step, i++) {
                    mpz_setbit(self->z, cur);
                }
            }
            else {
                for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                    bit = mpz_tstbit(tempx->z, i);
                    if (bit)
                        mpz_setbit(self->z, cur);
                    else
                        mpz_clrbit(self->z, cur);
                }
            }
            Py_DECREF((PyObject*)tempx);
        }
        return 0;
    }
    else {
        TYPE_ERROR("bit positions must be integers");
        return -1;
    }
    return -1;
}

/* Implement a multi-purpose iterator object that iterates over the bits in
 * an xmpz. Three different iterators can be created:
 *   1) xmpz.iter_bits(start=0, stop=-1) will return True/False for each bit
 *      position in the xmpz, beginning with bit 'start'. If stop is specified,
 *      the xmpz will be padded with 0-bits (False) until stop is reached.
 *   2) xmpz.iter_set(start=0, stop=-1, scale=1, offset=0) will return
 *      (scale*bit_position + offset) when bit_position is set, beginning at
 *      'start', ending at 'stop'.
 *   3) xmpz.iter_clear(start=0, stop=-1, scale=1, offset=0) will return
 *      (scale*bit_position + offset) when bit_position is clear, beginning at
 *      'start', ending at 'stop'.
 *
 */

static GMPy_Iter_Object *
GMPy_Iter_New(void)
{
    GMPy_Iter_Object *result;

    if ((result = PyObject_New(GMPy_Iter_Object,
                               &GMPy_Iter_Type))) {
        result->bitmap = NULL;
        result->start = 0;
        result->stop = (mp_bitcnt_t)-1;
        result->iter_type = 1;
    }
    return result;
}

static void
GMPy_Iter_Dealloc(GMPy_Iter_Object *self)
{
    Py_XDECREF((PyObject*)self->bitmap);
    PyObject_Free(self);
}

static PyObject *
GMPy_Iter_Next(GMPy_Iter_Object *self) {
    PyObject *result = 0;
    mp_bitcnt_t temp, current_stop;

    if (self->stop == (mp_bitcnt_t)(-1))
        current_stop = mpz_sizeinbase(self->bitmap->z, 2);
    else
        current_stop = self->stop;

    switch (self->iter_type) {
        case 1:
            if (self->start >= current_stop)
                PyErr_SetNone(PyExc_StopIteration);
            else {
                temp = mpz_tstbit(self->bitmap->z, self->start);
                self->start += 1;
                result = temp ? Py_True : Py_False;
                Py_INCREF(result);
            }
            break;
        case 2:
            if (self->start >= current_stop)
                PyErr_SetNone(PyExc_StopIteration);
            else {
                temp = mpz_scan1(self->bitmap->z, self->start);
                if (temp == (mp_bitcnt_t)(-1))
                    PyErr_SetNone(PyExc_StopIteration);
                else {
                    self->start = temp + 1;
                    result = PyLong_FromSsize_t(temp);
                }
            }
            break;
        case 3:
            if (self->start >= current_stop)
                PyErr_SetNone(PyExc_StopIteration);
            else {
                temp = mpz_scan0(self->bitmap->z, self->start);
                if (temp >= current_stop)
                    PyErr_SetNone(PyExc_StopIteration);
                else {
                    self->start = temp + 1;
                    result = PyLong_FromSsize_t(temp);
                }
            }
            break;
        default:
            SYSTEM_ERROR("Illegal iter_type in gmpy2.Iterator.");
    }
    return result;
}

static PyObject *
GMPy_Iter_Repr(GMPy_Iter_Object *self)
{
    return PyUnicode_FromString("<gmpy2.Iterator>");
}

PyDoc_STRVAR(GMPy_doc_xmpz_method_iter_bits,
"iter_bits($self, start=0, stop=-1)\n--\n\n"
"Return `True` or `False` for each bit position in self beginning at\n"
"'start'. If a positive value is specified for 'stop', iteration is\n"
"continued until 'stop' is reached. If a negative value is specified,\n"
"iteration is continued until the last 1-bit. Note: the value of the\n"
"underlying `xmpz` object can change during iteration.");

static PyObject *
GMPy_XMPZ_Method_IterBits(PyObject *self, PyObject *args, PyObject *kwargs)
{
    GMPy_Iter_Object *result;
    Py_ssize_t start = 0, stop = -1;

    static char *kwlist[] = {"start", "stop", NULL };

    if (!(result = GMPy_Iter_New())) {
        return NULL;
    }

    if (!(PyArg_ParseTupleAndKeywords(args, kwargs, "|nn", kwlist, &start, &stop))) {
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    result->iter_type = 1;
    result->bitmap = (XMPZ_Object*)self;
    Py_INCREF(self);
    result->start = start;
    result->stop = stop;
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_xmpz_method_iter_set,
"iter_set($self, start=0, stop=-1)\n--\n\n"
"Return an iterator yielding the bit position for every bit that\n"
"is set in self, beginning at 'start'. If a positive value is\n"
"specified for 'stop', iteration is continued until 'stop' is\n"
"reached. To match the behavior of slicing, 'stop' is not included.\n"
"If a negative value is specified, iteration is continued until\n"
"the last 1-bit. Note: the value of the underlying `xmpz` object can\n"
"change during iteration.");

static PyObject *
GMPy_XMPZ_Method_IterSet(PyObject *self, PyObject *args, PyObject *kwargs)
{
    GMPy_Iter_Object *result;
    Py_ssize_t start = 0, stop = -1;

    static char *kwlist[] = {"start", "stop", NULL };

    if (!(result = GMPy_Iter_New())) {
        return NULL;
    }

    if (!(PyArg_ParseTupleAndKeywords(args, kwargs, "|nn", kwlist, &start, &stop))) {
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    result->iter_type = 2;
    result->bitmap = (XMPZ_Object*)self;
    Py_INCREF(self);
    result->start = start;
    result->stop = stop;
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_xmpz_method_iter_clear,
"iter_clear($self, start=0, stop=-1)\n--\n\n"
"Return every bit position that is clear in self, beginning at\n"
"'start'. If a positive value is specified for 'stop', iteration\n"
"is continued until 'stop' is reached. If a negative value is specified,\n"
"iteration is continued until the last 1-bit. Note: the value of the\n"
"underlying `xmpz` object can change during iteration.");

static PyObject *
GMPy_XMPZ_Method_IterClear(PyObject *self, PyObject *args, PyObject *kwargs)
{
    GMPy_Iter_Object *result;
    Py_ssize_t start = 0, stop = -1;

    static char *kwlist[] = {"start", "stop", NULL };

    if (!(result = GMPy_Iter_New())) {
        return NULL;
    }

    if (!(PyArg_ParseTupleAndKeywords(args, kwargs, "|nn", kwlist, &start, &stop))) {
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    result->iter_type = 3;
    result->bitmap = (XMPZ_Object*)self;
    Py_INCREF(self);
    result->start = start;
    result->stop = stop;
    return (PyObject*)result;
}

static PyObject *
GMPy_XMPZ_Attrib_GetNumer(XMPZ_Object *self, void *closure)
{
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
}

static PyObject *
GMPy_XMPZ_Attrib_GetReal(XMPZ_Object *self, void *closure)
{
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
}

static PyObject *
GMPy_XMPZ_Attrib_GetDenom(XMPZ_Object *self, void *closure)
{
    XMPZ_Object *result;

    if ((result = GMPy_XMPZ_New(NULL))) {
        mpz_set_ui(result->z, 1);
    }
    return (PyObject*)result;
}

static PyObject *
GMPy_XMPZ_Attrib_GetImag(XMPZ_Object *self, void *closure)
{
    XMPZ_Object *result;

    if ((result = GMPy_XMPZ_New(NULL))) {
        mpz_set_ui(result->z, 0);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_xmpz_method_sizeof,
"__sizeof__($self)\n--\n\n"
"Returns the amount of memory consumed by self. Note: deleted xmpz objects\n"
"are reused and may or may not be resized when a new value is assigned.");

static PyObject *
GMPy_XMPZ_Method_SizeOf(PyObject *self, PyObject *other)
{
    return PyLong_FromSize_t(sizeof(XMPZ_Object) + \
        (MPZ(self)->_mp_alloc * sizeof(mp_limb_t)));
}

static PyTypeObject GMPy_Iter_Type =
{
    PyVarObject_HEAD_INIT(0, 0)
    "gmpy2 iterator",                       /* tp_name          */
    sizeof(GMPy_Iter_Object),               /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    (destructor) GMPy_Iter_Dealloc,         /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) GMPy_Iter_Repr,              /* tp_repr          */
        0,                                  /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
        0,                                  /* tp_hash          */
        0,                                  /* tp_call          */
        0,                                  /* tp_str           */
        0,                                  /* tp_getattro      */
        0,                                  /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
    "GMPY2 Iterator Object",                /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
        0,                                  /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
    PyObject_SelfIter,                      /* tp_iter          */
    (iternextfunc)GMPy_Iter_Next,           /* tp_iternext      */
};
