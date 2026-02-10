/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpfr.c                                                            *
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


/* General purpose functions are defined here. They are used as an alternative
 * to code bloat via macro overuse.
 */

static inline void
_GMPy_MPFR_Cleanup(MPFR_Object **v, CTXT_Object *ctext)
{
    /* GMPY_MPFR_CHECK_RANGE(V, CTX) */
    if (mpfr_regular_p((*v)->f) &&
        (!(((*v)->f->_mpfr_exp >= ctext->ctx.emin) &&
           ((*v)->f->_mpfr_exp <= ctext->ctx.emax)))) {
        mpfr_exp_t _oldemin, _oldemax;
        _oldemin = mpfr_get_emin();
        _oldemax = mpfr_get_emax();
        mpfr_set_emin(ctext->ctx.emin);
        mpfr_set_emax(ctext->ctx.emax);
        (*v)->rc = mpfr_check_range((*v)->f, (*v)->rc, GET_MPFR_ROUND(ctext));
        mpfr_set_emin(_oldemin);
        mpfr_set_emax(_oldemax);
    }

    /* GMPY_MPFR_SUBNORMALIZE(V, CTX) */
    if (ctext->ctx.subnormalize &&
        (*v)->f->_mpfr_exp >= ctext->ctx.emin &&
        (*v)->f->_mpfr_exp <= ctext->ctx.emin + mpfr_get_prec((*v)->f) - 2) {
        mpfr_exp_t _oldemin, _oldemax;
        _oldemin = mpfr_get_emin();
        _oldemax = mpfr_get_emax();
        mpfr_set_emin(ctext->ctx.emin);
        mpfr_set_emax(ctext->ctx.emax);
        (*v)->rc = mpfr_subnormalize((*v)->f, (*v)->rc, GET_MPFR_ROUND(ctext));
        mpfr_set_emin(_oldemin);
        mpfr_set_emax(_oldemax);
    }

    /* GMPY_MPFR_EXCEPTIONS(V, CTX) */
    ctext->ctx.underflow |= mpfr_underflow_p();
    ctext->ctx.overflow |= mpfr_overflow_p();
    ctext->ctx.invalid |= mpfr_nanflag_p();
    ctext->ctx.inexact |= mpfr_inexflag_p();
    ctext->ctx.divzero |= mpfr_divby0_p();
    if (ctext->ctx.traps) {
        if ((ctext->ctx.traps & TRAP_UNDERFLOW) && mpfr_underflow_p()) {
            PyErr_SetString(GMPyExc_Underflow, "underflow");
            Py_XDECREF((PyObject*)(*v));
            (*v) = NULL;
        }
        if ((ctext->ctx.traps & TRAP_OVERFLOW) && mpfr_overflow_p()) {
            PyErr_SetString(GMPyExc_Overflow, "overflow");
            Py_XDECREF((PyObject*)(*v));
            (*v) = NULL;
        }
        if ((ctext->ctx.traps & TRAP_INEXACT) && mpfr_inexflag_p()) {
            PyErr_SetString(GMPyExc_Inexact, "inexact result");
            Py_XDECREF((PyObject*)(*v));
            (*v) = NULL;
        }
        if ((ctext->ctx.traps & TRAP_INVALID) && mpfr_nanflag_p()) {
            PyErr_SetString(GMPyExc_Invalid, "invalid operation");
            Py_XDECREF((PyObject*)(*v));
            (*v) = NULL;
        }
        if ((ctext->ctx.traps & TRAP_DIVZERO) && mpfr_divby0_p()) {
            PyErr_SetString(GMPyExc_DivZero, "division by zero");
            Py_XDECREF((PyObject*)(*v));
            (*v) = NULL;
        }
    }
}

static inline void
_GMPy_MPFR_CleanupV2(MPFR_Object **v, CTXT_Object *ctext)
{
    /* GMPY_MPFR_EXCEPTIONS(V, CTX) */
    ctext->ctx.underflow |= mpfr_underflow_p();
    ctext->ctx.overflow |= mpfr_overflow_p();
    ctext->ctx.invalid |= mpfr_nanflag_p();
    ctext->ctx.inexact |= mpfr_inexflag_p();
    ctext->ctx.divzero |= mpfr_divby0_p();
    if (ctext->ctx.traps) {
        if ((ctext->ctx.traps & TRAP_UNDERFLOW) && mpfr_underflow_p()) {
            PyErr_SetString(GMPyExc_Underflow, "underflow");
            Py_XDECREF((PyObject*)(*v));
            (*v) = NULL;
        }
        if ((ctext->ctx.traps & TRAP_OVERFLOW) && mpfr_overflow_p()) {
            PyErr_SetString(GMPyExc_Overflow, "overflow");
            Py_XDECREF((PyObject*)(*v));
            (*v) = NULL;
        }
        if ((ctext->ctx.traps & TRAP_INEXACT) && mpfr_inexflag_p()) {
            PyErr_SetString(GMPyExc_Inexact, "inexact result");
            Py_XDECREF((PyObject*)(*v));
            (*v) = NULL;
        }
        if ((ctext->ctx.traps & TRAP_INVALID) && mpfr_nanflag_p()) {
            PyErr_SetString(GMPyExc_Invalid, "invalid operation");
            Py_XDECREF((PyObject*)(*v));
            (*v) = NULL;
        }
        if ((ctext->ctx.traps & TRAP_DIVZERO) && mpfr_divby0_p()) {
            PyErr_SetString(GMPyExc_DivZero, "division by zero");
            Py_XDECREF((PyObject*)(*v));
            (*v) = NULL;
        }
    }
}

PyDoc_STRVAR(GMPy_doc_mpfr,
"mpfr(n=0, /, precision=0)\n"
"mpfr(n, /, precision, context)\n"
"mpfr(s, /, precision=0, base=0)\n"
"mpfr(s, /, precision, base, context)\n\n"
"Return a floating-point number after converting a numeric value n or\n"
"a string s made of digits in the given base.\n\n"
"A string can be with fraction-part (with a period as a separator)\n"
"and/or exponent-part with an exponent marker 'e' or 'E' for bases up to\n"
"10, else '@' in any base.  In bases 2 and 16, the exponent prefix can also\n"
"be 'p' or 'P', in which case the exponent indicates\n"
"a multiplication by a power of 2 instead of the base.  The value of\n"
"an exponent is always written in base 10.  The fractional-part digits\n"
"are parsed the same as the `mpz` type constructor\n"
"does and both the whole number and exponent-part optionally can be\n"
"preceded by ‘+’ or ‘-’.  Every input, accepted by the `float` type\n"
"constructor or the `float.fromhex` method is also accepted.\n\n"
"If a precision greater than or equal to 2 is specified, then it\n"
"is used.  A precision of 0 (the default) implies the precision of either\n"
"the specified context or the current context is used.\n"
"A precision of 1 minimizes the loss of precision by following\n"
"these rules:\n\n"
"    1) If n is a radix-2 floating-point number, then the full\n"
"       precision of n is retained.\n"
"    2) If n is an integer, then the precision is the bit length\n"
"       of the integer.\n");

static PyNumberMethods mpfr_number_methods =
{
    .nb_add = (binaryfunc) GMPy_Number_Add_Slot,
    .nb_subtract = (binaryfunc) GMPy_Number_Sub_Slot,
    .nb_multiply = (binaryfunc) GMPy_Number_Mul_Slot,
    .nb_remainder = (binaryfunc) GMPy_Number_Mod_Slot,
    .nb_divmod = (binaryfunc) GMPy_Number_DivMod_Slot,
    .nb_power = (ternaryfunc) GMPy_Number_Pow_Slot,
    .nb_negative = (unaryfunc) GMPy_MPFR_Minus_Slot,
    .nb_positive = (unaryfunc) GMPy_MPFR_Plus_Slot,
    .nb_absolute = (unaryfunc) GMPy_MPFR_Abs_Slot,
    .nb_bool = (inquiry) GMPy_MPFR_NonZero_Slot,
    .nb_int = (unaryfunc) GMPy_MPFR_Int_Slot,
    .nb_float = (unaryfunc) GMPy_MPFR_Float_Slot,
    .nb_floor_divide = (binaryfunc) GMPy_Number_FloorDiv_Slot,
    .nb_true_divide = (binaryfunc) GMPy_Number_TrueDiv_Slot,
};

static PyGetSetDef Pympfr_getseters[] =
{
    {"precision", (getter)GMPy_MPFR_GetPrec_Attrib, NULL, "precision in bits", NULL},
    {"rc", (getter)GMPy_MPFR_GetRc_Attrib, NULL, "return code", NULL},
    {"imag", (getter)GMPy_MPFR_GetImag_Attrib, NULL, "imaginary component", NULL},
    {"real", (getter)GMPy_MPFR_GetReal_Attrib, NULL, "real component", NULL},
    {"_mpf_", (getter)GMPy_MPFR_Get_Mpmath_MPF_Tuple, NULL, "return raw mpf tuple", NULL},
    {NULL}
};

static PyMethodDef Pympfr_methods [] =
{
    { "__ceil__", GMPy_MPFR_Method_Ceil, METH_NOARGS, GMPy_doc_mpfr_ceil_method },
    { "__floor__", GMPy_MPFR_Method_Floor, METH_NOARGS, GMPy_doc_mpfr_floor_method },
    { "__format__", GMPy_MPFR_Format, METH_VARARGS, GMPy_doc_mpfr_format },
    { "__round__", GMPy_MPFR_Method_Round10, METH_VARARGS, GMPy_doc_method_round10 },
    { "__sizeof__", GMPy_MPFR_SizeOf_Method, METH_NOARGS, GMPy_doc_mpfr_sizeof_method },
    { "__trunc__", GMPy_MPFR_Method_Trunc, METH_NOARGS, GMPy_doc_mpfr_trunc_method },
    { "as_integer_ratio", GMPy_MPFR_Integer_Ratio_Method, METH_NOARGS, GMPy_doc_method_integer_ratio },
    { "as_mantissa_exp", GMPy_MPFR_Mantissa_Exp_Method, METH_NOARGS, GMPy_doc_method_mantissa_exp },
    { "as_simple_fraction", (PyCFunction)GMPy_MPFR_Simple_Fraction_Method, METH_VARARGS | METH_KEYWORDS, GMPy_doc_method_simple_fraction },
    { "conjugate", GMPy_MP_Method_Conjugate, METH_NOARGS, GMPy_doc_mp_method_conjugate },
    { "digits", GMPy_MPFR_Digits_Method, METH_VARARGS, GMPy_doc_mpfr_digits_method },
    { "is_finite", GMPy_Number_Method_Is_Finite, METH_NOARGS, GMPy_doc_method_is_finite },
    { "is_infinite", GMPy_Number_Method_Is_Infinite, METH_NOARGS, GMPy_doc_method_is_infinite },
    { "is_integer", GMPy_MPFR_Is_Integer_Method, METH_NOARGS, GMPy_doc_method_is_integer },
    { "is_nan", GMPy_Number_Method_Is_NAN, METH_NOARGS, GMPy_doc_method_is_nan },
    { "is_regular", GMPy_MPFR_Is_Regular_Method, METH_NOARGS, GMPy_doc_method_is_regular },
    { "is_signed", GMPy_MPFR_Is_Signed_Method, METH_NOARGS, GMPy_doc_method_is_signed },
    { "is_zero", GMPy_Number_Method_Is_Zero, METH_NOARGS, GMPy_doc_method_is_zero },
    { NULL }
};

static PyTypeObject MPFR_Type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "gmpy2.mpfr",
    .tp_basicsize = sizeof(MPFR_Object),
    .tp_dealloc = (destructor) GMPy_MPFR_Dealloc,
    .tp_repr = (reprfunc) GMPy_MPFR_Repr_Slot,
    .tp_as_number = &mpfr_number_methods,
    .tp_hash = (hashfunc) GMPy_MPFR_Hash_Slot,
    .tp_str = (reprfunc) GMPy_MPFR_Str_Slot,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = GMPy_doc_mpfr,
    .tp_richcompare = (richcmpfunc)&GMPy_RichCompare_Slot,
    .tp_methods = Pympfr_methods,
    .tp_getset = Pympfr_getseters,
    .tp_new = GMPy_MPFR_NewInit,
};
