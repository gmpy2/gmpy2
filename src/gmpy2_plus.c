/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_plus.c                                                            *
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

/* This file implements __pos__ and context.plus().
 *
 * Public API
 * ==========
 * The following function is available as part of GMPY2's C API. If the value
 * of context is NULL, then the function should use the currently active
 * context.
 *
 *   GMPy_Number_Plus(Number, context)
 *
 * Private API
 * ===========
 *   GMPy_MPZ_Plus_Slot
 *   GMPy_MPQ_Plus_Slot
 *   GMPy_MPFR_Plus_Slot
 *   GMPy_MPC_Plus_Slot
 *
 *   GMPy_Integer_Plus(Integer, context|NULL)
 *   GMPy_Rational_Plus(Rational, context|NULL)
 *   GMPy_Real_Plus(Real, context|NULL)
 *   GMPy_Complex_Plus(Complex, context|NULL)
 *
 *   GMPy_Context_Plus(context, args)
 */

static PyObject *
GMPy_Integer_PlusWithType(PyObject *x, int xtype, CTXT_Object *context)
{
    return (PyObject*)GMPy_MPZ_From_IntegerWithType(x, xtype, context);
}

static PyObject *
GMPy_MPZ_Plus_Slot(MPZ_Object *x)
{
    Py_INCREF((PyObject*)x);
    return (PyObject*)x;
}

static PyObject *
GMPy_Rational_PlusWithType(PyObject *x, int xtype, CTXT_Object *context)
{
    return (PyObject*)GMPy_MPQ_From_RationalWithType(x, xtype, context);
}

static PyObject *
GMPy_MPQ_Plus_Slot(MPQ_Object *x)
{
    Py_INCREF((PyObject*)x);
    return (PyObject*)x;
}

static PyObject *
GMPy_Real_PlusWithType(PyObject *x, int xtype, CTXT_Object *context)
{
    return (PyObject*)GMPy_MPFR_From_RealWithType(x, xtype, 0, context);
}

static PyObject *
GMPy_MPFR_Plus_Slot(MPFR_Object *x)
{
    return (PyObject*)GMPy_MPFR_From_MPFR(x, 0, NULL);
}

static PyObject *
GMPy_Complex_PlusWithType(PyObject *x, int xtype, CTXT_Object *context)
{
    return (PyObject*)GMPy_MPC_From_ComplexWithType(x, xtype, 0, 0, context);
}

static PyObject *
GMPy_MPC_Plus_Slot(MPC_Object *x)
{
    return (PyObject*)GMPy_MPC_From_MPC(x, 0, 0, NULL);
}

static PyObject *
GMPy_Number_Plus(PyObject *x, CTXT_Object *context)
{
    int xtype = GMPy_ObjectType(x);

    if (IS_TYPE_INTEGER(xtype))
        return GMPy_Integer_PlusWithType(x, xtype, context);

    if (IS_TYPE_RATIONAL(xtype))
        return GMPy_Rational_PlusWithType(x, xtype, context);

    if (IS_TYPE_REAL(xtype))
        return GMPy_Real_PlusWithType(x, xtype, context);

    if (IS_TYPE_COMPLEX(xtype))
        return GMPy_Complex_PlusWithType(x, xtype, context);

    TYPE_ERROR("plus() argument type not supported");
    return NULL;
}

/* Implement context.plus(). The following code assumes it used a as method of
 * a context. */

PyDoc_STRVAR(GMPy_doc_context_plus,
"plus($self, x, /)\n--\n\n"
"Return +x, the context is applied to the result.");

static PyObject *
GMPy_Context_Plus(PyObject *self, PyObject *args)
{
    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("plus() requires 1 argument.");
        return NULL;
    }

    return GMPy_Number_Plus(PyTuple_GET_ITEM(args, 0), (CTXT_Object*)self);
}
