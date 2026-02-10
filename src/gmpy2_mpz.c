/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpz.c                                                             *
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

PyDoc_STRVAR(GMPy_doc_mpz,
"mpz(n=0, /)\n"
"mpz(s, /, base=0)\n\n"
"Return an immutable integer constructed from a numeric value n (truncating\n"
"n to its integer part) or a string s made of digits in the given base.\n"
"Every input, that is accepted by the `int` type constructor is also accepted.\n\n"
"The base may vary from 2 to 62, or if base is 0, then binary, octal, or\n"
"hexadecimal strings are recognized by leading '0b', '0o', or '0x'\n"
"characters (case is ignored), otherwise the string is assumed to be\n"
"decimal.  For bases up to 36, digits case is ignored.  For bases 37\n"
"to 62, upper-case letter represent the usual 10..35 range, while\n"
"lower-case letter represent 36..61.  Optionally the string can be\n"
"preceded by '+' or '-'.  White space and underscore is simply ignored.");

/* Since `gmpy2.mpz` is now a type and no longer a factory function, see
 * gmpy2_cache.c/GMPy_MPZ_NewInit for details on creation.
 */

static PyNumberMethods GMPy_MPZ_number_methods = {
    .nb_add = (binaryfunc) GMPy_Number_Add_Slot,
    .nb_subtract = (binaryfunc) GMPy_Number_Sub_Slot,
    .nb_multiply = (binaryfunc) GMPy_Number_Mul_Slot,
    .nb_remainder = (binaryfunc) GMPy_Number_Mod_Slot,
    .nb_divmod = (binaryfunc) GMPy_Number_DivMod_Slot,
    .nb_power = (ternaryfunc) GMPy_Number_Pow_Slot,
    .nb_negative = (unaryfunc) GMPy_MPZ_Minus_Slot,
    .nb_positive = (unaryfunc) GMPy_MPZ_Plus_Slot,
    .nb_absolute = (unaryfunc) GMPy_MPZ_Abs_Slot,
    .nb_bool = (inquiry) GMPy_MPZ_NonZero_Slot,
    .nb_invert = (unaryfunc) GMPy_MPZ_Invert_Slot,
    .nb_lshift = (binaryfunc) GMPy_MPZ_Lshift_Slot,
    .nb_rshift = (binaryfunc) GMPy_MPZ_Rshift_Slot,
    .nb_and = (binaryfunc) GMPy_MPZ_And_Slot,
    .nb_xor = (binaryfunc) GMPy_MPZ_Xor_Slot,
    .nb_or = (binaryfunc) GMPy_MPZ_Ior_Slot,
    .nb_int = (unaryfunc) GMPy_MPZ_Int_Slot,
    .nb_float = (unaryfunc) GMPy_MPZ_Float_Slot,
    .nb_floor_divide = (binaryfunc) GMPy_Number_FloorDiv_Slot,
    .nb_true_divide = (binaryfunc) GMPy_Number_TrueDiv_Slot,
    .nb_index = (unaryfunc) GMPy_MPZ_Int_Slot,
};

static PyMappingMethods GMPy_MPZ_mapping_methods = {
    (lenfunc)GMPy_MPZ_Method_Length,
    (binaryfunc)GMPy_MPZ_Method_SubScript,
    NULL
};

static PyGetSetDef GMPy_MPZ_getseters[] = {
    { "numerator", (getter)GMPy_MPZ_Attrib_GetNumer, NULL,
        "the numerator of a rational number in lowest terms", NULL },
    { "denominator", (getter)GMPy_MPZ_Attrib_GetDenom, NULL,
        "the denominator of a rational number in lowest terms", NULL },
    { "real", (getter)GMPy_MPZ_Attrib_GetReal, NULL,
        "the real part of a complex number", NULL },
    { "imag", (getter)GMPy_MPZ_Attrib_GetImag, NULL,
        "the imaginary part of a complex number", NULL },
    {NULL}
};

static PyMethodDef GMPy_MPZ_methods[] = {
    { "__format__", GMPy_MPZ_Format, METH_VARARGS, GMPy_doc_mpz_format },
    { "__ceil__", GMPy_MPZ_Method_Ceil, METH_NOARGS, GMPy_doc_mpz_method_ceil },
    { "__floor__", GMPy_MPZ_Method_Floor, METH_NOARGS, GMPy_doc_mpz_method_floor },
    { "__round__", (PyCFunction)GMPy_MPZ_Method_Round, METH_FASTCALL, GMPy_doc_mpz_method_round },
    { "__sizeof__", GMPy_MPZ_Method_SizeOf, METH_NOARGS, GMPy_doc_mpz_method_sizeof },
    { "__trunc__", GMPy_MPZ_Method_Trunc, METH_NOARGS, GMPy_doc_mpz_method_trunc },
    { "__array__", (PyCFunction)GMPy_MPZ_Method_Array, METH_FASTCALL | METH_KEYWORDS, GMPy_doc_mpz_method_array },
    { "bit_clear", GMPy_MPZ_bit_clear_method, METH_O, doc_bit_clear_method },
    { "bit_count", GMPy_MPZ_bit_count_method, METH_NOARGS, doc_bit_count_method },
    { "bit_flip", GMPy_MPZ_bit_flip_method, METH_O, doc_bit_flip_method },
    { "bit_length", GMPy_MPZ_bit_length_method, METH_NOARGS, doc_bit_length_method },
    { "bit_scan0", (PyCFunction)GMPy_MPZ_bit_scan0_method, METH_FASTCALL, doc_bit_scan0_method },
    { "bit_scan1", (PyCFunction)GMPy_MPZ_bit_scan1_method, METH_FASTCALL, doc_bit_scan1_method },
    { "bit_set", GMPy_MPZ_bit_set_method, METH_O, doc_bit_set_method },
    { "bit_test", GMPy_MPZ_bit_test_method, METH_O, doc_bit_test_method },
    { "conjugate", GMPy_MP_Method_Conjugate, METH_NOARGS, GMPy_doc_mp_method_conjugate },
    { "digits", GMPy_MPZ_Digits_Method, METH_VARARGS, GMPy_doc_mpz_digits_method },
    { "is_congruent", (PyCFunction)GMPy_MPZ_Method_IsCongruent, METH_FASTCALL, GMPy_doc_mpz_method_is_congruent },
    { "is_divisible", GMPy_MPZ_Method_IsDivisible, METH_O, GMPy_doc_mpz_method_is_divisible },
    { "is_even", GMPy_MPZ_Method_IsEven, METH_NOARGS, GMPy_doc_mpz_method_is_even },
    { "is_odd", GMPy_MPZ_Method_IsOdd, METH_NOARGS, GMPy_doc_mpz_method_is_odd },
    { "is_power", GMPy_MPZ_Method_IsPower, METH_NOARGS, GMPy_doc_mpz_method_is_power },
    { "is_prime", (PyCFunction)GMPy_MPZ_Method_IsPrime, METH_FASTCALL, GMPy_doc_mpz_method_is_prime },
    { "is_probab_prime", (PyCFunction)GMPy_MPZ_Method_IsProbabPrime, METH_FASTCALL, GMPy_doc_mpz_method_is_probab_prime },
    { "is_square", GMPy_MPZ_Method_IsSquare, METH_NOARGS, GMPy_doc_mpz_method_is_square },
    { "is_integer", GMPy_MPZ_Method_IsInteger, METH_NOARGS, GMPy_doc_mpz_method_is_integer },
    { "num_digits", (PyCFunction)GMPy_MPZ_Method_NumDigits, METH_FASTCALL, GMPy_doc_mpz_method_num_digits },
    { "as_integer_ratio", GMPy_MPZ_Method_As_Integer_Ratio, METH_NOARGS, GMPy_doc_mpz_method_as_integer_ratio },
    { "to_bytes", (PyCFunction)GMPy_MPZ_Method_To_Bytes, METH_FASTCALL | METH_KEYWORDS, GMPy_doc_mpz_method_to_bytes },
    { "from_bytes", (PyCFunction)GMPy_MPZ_Method_From_Bytes, METH_FASTCALL | METH_KEYWORDS | METH_CLASS, GMPy_doc_mpz_method_from_bytes },
    { NULL }
};

static PyTypeObject MPZ_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "gmpy2.mpz",
    .tp_basicsize = sizeof(MPZ_Object),
    .tp_dealloc = (destructor) GMPy_MPZ_Dealloc,
    .tp_repr = (reprfunc) GMPy_MPZ_Repr_Slot,
    .tp_as_number = &GMPy_MPZ_number_methods,
    .tp_as_mapping = &GMPy_MPZ_mapping_methods,
    .tp_hash = (hashfunc) GMPy_MPZ_Hash_Slot,
    .tp_str = (reprfunc) GMPy_MPZ_Str_Slot,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = GMPy_doc_mpz,
    .tp_richcompare = (richcmpfunc)&GMPy_RichCompare_Slot,
    .tp_methods = GMPy_MPZ_methods,
    .tp_getset = GMPy_MPZ_getseters,
    .tp_new = GMPy_MPZ_NewInit,
};
