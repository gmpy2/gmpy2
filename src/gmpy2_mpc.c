/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpc.c                                                             *
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

static void
_GMPy_MPC_Cleanup(MPC_Object **v, CTXT_Object *ctext)
{
    /* GMPY_MPC_CHECK_RANGE(V, CTX) */
    {
        int rcr, rci;
        rcr = MPC_INEX_RE((*v)->rc);
        rci = MPC_INEX_IM((*v)->rc);
        if (mpfr_regular_p(mpc_realref((*v)->c)) &&
            (!((mpc_realref((*v)->c)->_mpfr_exp >= ctext->ctx.emin) &&
               (mpc_realref((*v)->c)->_mpfr_exp <= ctext->ctx.emax)))) {
            mpfr_exp_t _oldemin, _oldemax;
            _oldemin = mpfr_get_emin();
            _oldemax = mpfr_get_emax();
            mpfr_set_emin(ctext->ctx.emin);
            mpfr_set_emax(ctext->ctx.emax);
            rcr = mpfr_check_range(mpc_realref((*v)->c), rcr, GET_REAL_ROUND(ctext));
            mpfr_set_emin(_oldemin);
            mpfr_set_emax(_oldemax);
        }
        if (mpfr_regular_p(mpc_imagref((*v)->c)) &&
            (!((mpc_imagref((*v)->c)->_mpfr_exp >= ctext->ctx.emin) &&
               (mpc_imagref((*v)->c)->_mpfr_exp <= ctext->ctx.emax)))) {
            mpfr_exp_t _oldemin, _oldemax;
            _oldemin = mpfr_get_emin();
            _oldemax = mpfr_get_emax();
            mpfr_set_emin(ctext->ctx.emin);
            mpfr_set_emax(ctext->ctx.emax);
            rci = mpfr_check_range(mpc_imagref((*v)->c), rci, GET_IMAG_ROUND(ctext));
            mpfr_set_emin(_oldemin);
            mpfr_set_emax(_oldemax);
        }
        (*v)->rc = MPC_INEX(rcr, rci);
    }

    /* GMPY_MPC_SUBNORMALIZE(V, CTX)  */
    {
        int rcr, rci;
        rcr = MPC_INEX_RE((*v)->rc);
        rci = MPC_INEX_IM((*v)->rc);
        if (ctext->ctx.subnormalize &&
            (!((mpc_realref((*v)->c)->_mpfr_exp >= ctext->ctx.emin) &&
               (mpc_realref((*v)->c)->_mpfr_exp <= ctext->ctx.emin + mpfr_get_prec(mpc_realref((*v)->c)) - 2)))) {
            mpfr_exp_t _oldemin, _oldemax;
            _oldemin = mpfr_get_emin();
            _oldemax = mpfr_get_emax();
            mpfr_set_emin(ctext->ctx.emin);
            mpfr_set_emax(ctext->ctx.emax);
            rcr = mpfr_subnormalize(mpc_realref((*v)->c), rcr, GET_REAL_ROUND(ctext));
            mpfr_set_emin(_oldemin);
            mpfr_set_emax(_oldemax);
        }
        if (ctext->ctx.subnormalize &&
            (!((mpc_imagref((*v)->c)->_mpfr_exp >= ctext->ctx.emin) &&
               (mpc_imagref((*v)->c)->_mpfr_exp <= ctext->ctx.emin + mpfr_get_prec(mpc_imagref((*v)->c)) - 2)))) {
            mpfr_exp_t _oldemin, _oldemax;
            _oldemin = mpfr_get_emin();
            _oldemax = mpfr_get_emax();
            mpfr_set_emin(ctext->ctx.emin);
            mpfr_set_emax(ctext->ctx.emax);
            rci = mpfr_check_range(mpc_imagref((*v)->c), rci, GET_IMAG_ROUND(ctext));
            mpfr_set_emin(_oldemin);
            mpfr_set_emax(_oldemax);
        }
        (*v)->rc = MPC_INEX(rcr, rci);
    }

    /* GMPY_MPC_EXCEPTIONS(V, CTX) */
    {
        int _invalid = 0, _underflow = 0, _overflow = 0, _inexact = 0;
        int rcr, rci;
        rcr = MPC_INEX_RE((*v)->rc);
        rci = MPC_INEX_IM((*v)->rc);
        if (MPC_IS_NAN_P(*v)) {
            ctext->ctx.invalid = 1;
            _invalid = 1;
        }
        if ((*v)->rc) {
            ctext->ctx.inexact = 1;
            _inexact = 1;
        }
        if ((rcr && mpfr_zero_p(mpc_realref((*v)->c))) || (rci && mpfr_zero_p(mpc_imagref((*v)->c)))) {
            ctext->ctx.underflow = 1;
            _underflow = 1;
        }
        if ((rcr && mpfr_inf_p(mpc_realref((*v)->c))) || (rci && mpfr_inf_p(mpc_imagref((*v)->c)))) {
            ctext->ctx.overflow = 1;
            _overflow = 1;
        }
        if (ctext->ctx.traps) {
            if ((ctext->ctx.traps & TRAP_UNDERFLOW) && _underflow) { \
                GMPY_UNDERFLOW("underflow");
                Py_XDECREF((PyObject*)(*v));
                (*v) = NULL;
            }
            if ((ctext->ctx.traps & TRAP_OVERFLOW) && _overflow) {
                GMPY_OVERFLOW("overflow");
                Py_XDECREF((PyObject*)(*v));
                (*v) = NULL;
            }
            if ((ctext->ctx.traps & TRAP_INEXACT) && _inexact) {
                GMPY_INEXACT("inexact result");
                Py_XDECREF((PyObject*)(*v));
                (*v) = NULL;
            }
            if ((ctext->ctx.traps & TRAP_INVALID) && _invalid) {
                GMPY_INVALID("invalid operation");
                Py_XDECREF((PyObject*)(*v));
                (*v) = NULL;
            }
        }
    }
}

PyDoc_STRVAR(GMPy_doc_mpc,
"mpc(c=0, /, precision=0)\n"
"mpc(c=0, /, precision, context)\n"
"mpc(real, /, imag=0, precision=0)\n"
"mpc(real, /, imag, precision, context)\n"
"mpc(s, /, precision=0, base=10)\n"
"mpc(s, /, precision, base, context)\n\n"
"Return a complex floating-point number constructed from a numeric value\n"
"c or from a pair of two non-complex numbers real and imag or from a\n"
"string s made of digits in the given base.\n\n"
"A string can be possibly with real-part and/or imaginary-part (that\n"
"have 'j' as a suffix), separated by '+' and parsed the same as the\n"
"`mpfr` constructor does (but the base must be up to 36).\n\n"
"The precision can be specified by either a single number that is used\n"
"for both the real and imaginary components, or as a pair of different\n"
"precisions for the real and imaginary components.  For every component,\n"
"the meaning of its precision value is the same as in the `mpfr`\n"
"type constructor.");

static PyMethodDef Pympc_methods[] =
{
    { "__complex__", GMPy_PyComplex_From_MPC, METH_NOARGS, GMPy_doc_mpc_complex },
    { "__format__", GMPy_MPC_Format, METH_VARARGS, GMPy_doc_mpc_format },
    { "__sizeof__", GMPy_MPC_SizeOf_Method, METH_NOARGS, GMPy_doc_mpc_sizeof_method },
    { "conjugate", GMPy_MPC_Conjugate_Method, METH_NOARGS, GMPy_doc_mpc_conjugate_method },
    { "digits", GMPy_MPC_Digits_Method, METH_VARARGS, GMPy_doc_mpc_digits_method },
    { "is_finite", GMPy_Number_Method_Is_Finite, METH_NOARGS, GMPy_doc_method_is_finite },
    { "is_infinite", GMPy_Number_Method_Is_Infinite, METH_NOARGS, GMPy_doc_method_is_infinite },
    { "is_nan", GMPy_Number_Method_Is_NAN, METH_NOARGS, GMPy_doc_method_is_nan },
    { "is_zero", GMPy_Number_Method_Is_Zero, METH_NOARGS, GMPy_doc_method_is_zero },
    { NULL }
};


static PyNumberMethods mpc_number_methods =
{
    .nb_add = (binaryfunc) GMPy_Number_Add_Slot,
    .nb_subtract = (binaryfunc) GMPy_Number_Sub_Slot,
    .nb_multiply = (binaryfunc) GMPy_Number_Mul_Slot,
    .nb_remainder = (binaryfunc) GMPy_Number_Mod_Slot,
    .nb_divmod = (binaryfunc) GMPy_Number_DivMod_Slot,
    .nb_power = (ternaryfunc) GMPy_Number_Pow_Slot,
    .nb_negative = (unaryfunc) GMPy_MPC_Minus_Slot,
    .nb_positive = (unaryfunc) GMPy_MPC_Plus_Slot,
    .nb_absolute = (unaryfunc) GMPy_MPC_Abs_Slot,
    .nb_bool = (inquiry) GMPy_MPC_NonZero_Slot,
    .nb_int = (unaryfunc) GMPy_MPC_Int_Slot,
    .nb_float = (unaryfunc) GMPy_MPC_Float_Slot,
    .nb_floor_divide = (binaryfunc) GMPy_Number_FloorDiv_Slot,
    .nb_true_divide = (binaryfunc) GMPy_Number_TrueDiv_Slot,
};

static PyGetSetDef Pympc_getseters[] =
{
    {"precision", (getter)GMPy_MPC_GetPrec_Attrib, NULL, "precision in bits", NULL},
    {"rc", (getter)GMPy_MPC_GetRc_Attrib, NULL, "return code", NULL},
    {"imag", (getter)GMPy_MPC_GetImag_Attrib, NULL, "imaginary component", NULL},
    {"real", (getter)GMPy_MPC_GetReal_Attrib, NULL, "real component", NULL},
    {NULL}
};

static PyTypeObject MPC_Type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "mpc",
    .tp_basicsize = sizeof(MPC_Object),
    .tp_dealloc = (destructor) GMPy_MPC_Dealloc,
    .tp_repr = (reprfunc) GMPy_MPC_Repr_Slot,
    .tp_as_number = &mpc_number_methods,
    .tp_hash = (hashfunc) GMPy_MPC_Hash_Slot,
    .tp_str = (reprfunc) GMPy_MPC_Str_Slot,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = GMPy_doc_mpc,
    .tp_richcompare = (richcmpfunc)&GMPy_RichCompare_Slot,
    .tp_methods = Pympc_methods,
    .tp_getset = Pympc_getseters,
    .tp_new = GMPy_MPC_NewInit,
};
