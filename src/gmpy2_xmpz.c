/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_xmpz.c                                                            *
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

PyDoc_STRVAR(GMPy_doc_xmpz,
"xmpz(n=0, /)\n"
"xmpz(s, /, base=0)\n\n"
"Return a mutable integer constructed from a numeric value n\n"
"or a string s made of digits in the given base.  Every input,\n"
"that is accepted by the `mpz` type constructor is also accepted.\n\n"
"Note: This type can be faster when used for augmented assignment\n"
"(+=, -=, etc), but `xmpz` objects cannot be used as dictionary keys.");


static PyNumberMethods GMPy_XMPZ_number_methods =
{
    .nb_add = (binaryfunc) GMPy_Number_Add_Slot,
    .nb_subtract = (binaryfunc) GMPy_Number_Sub_Slot,
    .nb_multiply = (binaryfunc) GMPy_Number_Mul_Slot,
    .nb_remainder = (binaryfunc) GMPy_Number_Mod_Slot,
    .nb_divmod = (binaryfunc) GMPy_Number_DivMod_Slot,
    .nb_power = (ternaryfunc) GMPy_Number_Pow_Slot,
    .nb_negative = (unaryfunc) GMPy_XMPZ_Neg_Slot,
    .nb_positive = (unaryfunc) GMPy_XMPZ_Pos_Slot,
    .nb_absolute = (unaryfunc) GMPy_XMPZ_Abs_Slot,
    .nb_bool = (inquiry) GMPy_XMPZ_NonZero_Slot,
    .nb_invert = (unaryfunc) GMPy_XMPZ_Com_Slot,
    .nb_lshift = (binaryfunc) GMPy_MPZ_Lshift_Slot,
    .nb_rshift = (binaryfunc) GMPy_MPZ_Rshift_Slot,
    .nb_and = (binaryfunc) GMPy_MPZ_And_Slot,
    .nb_xor = (binaryfunc) GMPy_MPZ_Xor_Slot,
    .nb_or = (binaryfunc) GMPy_MPZ_Ior_Slot,
    .nb_int = (unaryfunc) GMPy_MPZ_Int_Slot,
    .nb_float = (unaryfunc) GMPy_MPZ_Float_Slot,
    .nb_inplace_add = (binaryfunc) GMPy_XMPZ_IAdd_Slot,
    .nb_inplace_subtract = (binaryfunc) GMPy_XMPZ_ISub_Slot,
    .nb_inplace_multiply = (binaryfunc) GMPy_XMPZ_IMul_Slot,
    .nb_inplace_remainder = (binaryfunc) GMPy_XMPZ_IRem_Slot,
    .nb_inplace_power = (ternaryfunc) GMPy_XMPZ_IPow_Slot,
    .nb_inplace_lshift = (binaryfunc) GMPy_XMPZ_ILshift_Slot,
    .nb_inplace_rshift = (binaryfunc) GMPy_XMPZ_IRshift_Slot,
    .nb_inplace_and = (binaryfunc) GMPy_XMPZ_IAnd_Slot,
    .nb_inplace_xor = (binaryfunc) GMPy_XMPZ_IXor_Slot,
    .nb_inplace_or = (binaryfunc) GMPy_XMPZ_IIor_Slot,
    .nb_floor_divide = (binaryfunc) GMPy_Number_FloorDiv_Slot,
    .nb_true_divide = (binaryfunc) GMPy_Number_TrueDiv_Slot,
    .nb_inplace_floor_divide = (binaryfunc) GMPy_XMPZ_IFloorDiv_Slot,
    .nb_index = (unaryfunc) GMPy_MPZ_Int_Slot,
};

static PyMappingMethods GMPy_XMPZ_mapping_methods = {
    (lenfunc)GMPy_XMPZ_Method_Length,
    (binaryfunc)GMPy_XMPZ_Method_SubScript,
    (objobjargproc)GMPy_XMPZ_Method_AssignSubScript
};

static PyGetSetDef GMPy_XMPZ_getseters[] =
{
    { "numerator", (getter)GMPy_XMPZ_Attrib_GetNumer, NULL,
        "the numerator of a rational number in lowest terms", NULL },
    { "denominator", (getter)GMPy_XMPZ_Attrib_GetDenom, NULL,
        "the denominator of a rational number in lowest terms", NULL },
    { "real", (getter)GMPy_XMPZ_Attrib_GetReal, NULL,
        "the real part of a complex number", NULL },
    { "imag", (getter)GMPy_XMPZ_Attrib_GetImag, NULL,
        "the imaginary part of a complex number", NULL },
    {NULL}
};

static PyMethodDef GMPy_XMPZ_methods [] =
{
    { "__format__", GMPy_MPZ_Format, METH_VARARGS, GMPy_doc_mpz_format },
    { "__sizeof__", GMPy_XMPZ_Method_SizeOf, METH_NOARGS, GMPy_doc_xmpz_method_sizeof },
    { "bit_clear", GMPy_MPZ_bit_clear_method, METH_O, doc_bit_clear_method },
    { "bit_flip", GMPy_MPZ_bit_flip_method, METH_O, doc_bit_flip_method },
    { "bit_length", GMPy_MPZ_bit_length_method, METH_NOARGS, doc_bit_length_method },
    { "bit_count", GMPy_MPZ_bit_count_method, METH_NOARGS, doc_bit_count_method },
    { "bit_scan0", (PyCFunction)GMPy_MPZ_bit_scan0_method, METH_FASTCALL, doc_bit_scan0_method },
    { "bit_scan1", (PyCFunction)GMPy_MPZ_bit_scan1_method, METH_FASTCALL, doc_bit_scan1_method },
    { "bit_set", GMPy_MPZ_bit_set_method, METH_O, doc_bit_set_method },
    { "bit_test", GMPy_MPZ_bit_test_method, METH_O, doc_bit_test_method },
    { "conjugate", GMPy_MP_Method_Conjugate, METH_NOARGS, GMPy_doc_mp_method_conjugate },
    { "copy", GMPy_XMPZ_Method_Copy, METH_NOARGS, GMPy_doc_xmpz_method_copy },
    { "digits", GMPy_XMPZ_Digits_Method, METH_VARARGS, GMPy_doc_mpz_digits_method },
    { "iter_bits", (PyCFunction)GMPy_XMPZ_Method_IterBits, METH_VARARGS | METH_KEYWORDS, GMPy_doc_xmpz_method_iter_bits },
    { "iter_clear", (PyCFunction)GMPy_XMPZ_Method_IterClear, METH_VARARGS | METH_KEYWORDS, GMPy_doc_xmpz_method_iter_clear },
    { "iter_set", (PyCFunction)GMPy_XMPZ_Method_IterSet, METH_VARARGS | METH_KEYWORDS, GMPy_doc_xmpz_method_iter_set },
    { "make_mpz", GMPy_XMPZ_Method_MakeMPZ, METH_NOARGS, GMPy_doc_xmpz_method_make_mpz },
    { "num_digits", (PyCFunction)GMPy_MPZ_Method_NumDigits, METH_FASTCALL, GMPy_doc_mpz_method_num_digits },
    { "num_limbs", GMPy_XMPZ_Method_NumLimbs, METH_NOARGS, GMPy_doc_xmpz_method_num_limbs },
    { "limbs_read", GMPy_XMPZ_Method_LimbsRead, METH_NOARGS, GMPy_doc_xmpz_method_limbs_read },
    { "limbs_write", GMPy_XMPZ_Method_LimbsWrite, METH_O, GMPy_doc_xmpz_method_limbs_write },
    { "limbs_modify", GMPy_XMPZ_Method_LimbsModify, METH_O, GMPy_doc_xmpz_method_limbs_modify },
    { "limbs_finish", GMPy_XMPZ_Method_LimbsFinish, METH_O, GMPy_doc_xmpz_method_limbs_finish },
    { NULL }
};

static PyTypeObject XMPZ_Type =
{
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "gmpy2.xmpz",
    .tp_basicsize = sizeof(XMPZ_Object),
    .tp_dealloc = (destructor) GMPy_XMPZ_Dealloc,
    .tp_repr = (reprfunc) GMPy_XMPZ_Repr_Slot,
    .tp_as_number = &GMPy_XMPZ_number_methods,
    .tp_as_mapping = &GMPy_XMPZ_mapping_methods,
    .tp_str = (reprfunc) GMPy_XMPZ_Str_Slot,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = GMPy_doc_xmpz,
    .tp_richcompare = (richcmpfunc)&GMPy_RichCompare_Slot,
    .tp_methods = GMPy_XMPZ_methods,
    .tp_getset = GMPy_XMPZ_getseters,
    .tp_new = GMPy_XMPZ_NewInit,
};
