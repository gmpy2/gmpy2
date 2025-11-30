/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2.c                                                                 *
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

 /* Todo list
  * ---------
  * Add all MPFR and MPC functions as context methods.
  * All MPFR and MPC functions need to set exponent range on entry. The
  *    current approach where only set_context() and context.__enter__ set
  *    the exponent range fails for context methods.
  * Should a read-only (or template) context prevent the setting of
  *    exception flags?
  * Add context option to control the result of integer division:
  *    integer (mpz), exact (mpq), or true (mpfr).
  * Add modular arithmetic functions.
  * Implement Chinese Remainder Theorem.
  * Update PRP code.
  */

/*
 * originally written for GMP-2.0 (by AMK...?)
 * Rewritten by Niels Möller, May 1996
 *
 * Version for GMP-4, Python 2.X, with support for MSVC++6,
 * addition of mpf's, &c: Alex Martelli (now aleaxit@gmail.com, Nov 2000).
 * cleanups & reorgs leading to 1.0: Alex Martelli (until Aug 2003)
 * further cleanups and bugfixes leading to 1.01, Alex Martelli (Nov 2005)
 * minor bugfixes+new decimal (&c) support to 1.02, Alex Martelli (Feb 2006)
 * various bugfixes for 64-bit platforms, 1.03, aleaxit and casevh (Jun 2008)
 * rich comparisons, 1.04, aleaxit and casevh (Jan 2009)
 * support for Python 3.x, 1.10, casevh (Oct 2009)
 *
 * Some hacks by Gustavo Niemeyer <niemeyer@conectiva.com>.
 *
 ************************************************************************
 *
 * Discussion on sizes of C integer types.
 *
 * GMP, MPFR, and MPC use typedef to create integer objects with
 * different sizes. It can become confusing to map the different types
 * onto the standard C types used by Python's C API. Below are external
 * types and how they map to C types. The assumptions are verified when
 * the module is initialized.
 *
 * mp_limb_t: This is usually an 'unsigned long' but is an 'unsigned
 *     long long' on certain 64-bit Windows builds.
 *
 * mp_bitcnt_t: This is usually an 'unsigned long' but is an 'unsigned
 *     long long' on MPIR/64-bit Windows. 'size_t' is the best match.
 *
 * mp_size_t: This is a 'long'.
 *
 * mpfr_rnd_t: This is an 'int'.
 *
 * mpfr_prec_t: This is a 'long'.
 *
 * mpfr_sign_t: This is an 'int'.
 *
 * mpfr_exp_t: This is currently the same as mp_exp_t but will change
 *     to a signed 64-bit integer in the future.
 *
 * mpc_rnd_t: This is an 'int'.
 *
 * mpc_prec_t: See mpfr_exp_t.
 *
 ************************************************************************
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <signal.h>

#define GMPY2_MODULE
#include "gmpy2.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Global data declarations begin here.                                    *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* The following global strings are used by gmpy_misc.c. */

#define GMPY_VERSION "2.3.0a2"

char gmpy_version[] = GMPY_VERSION;

char gmpy_license[] = "\
The GMPY2 source code is licensed under LGPL 3 or later. The supported \
versions of the GMP, MPFR, and MPC libraries are also licensed under \
LGPL 3 or later.";

/* The following global structures are used by gmpy_cache.c.
 */

#if !defined(PYPY_VERSION)
#define CACHE_SIZE (100)
#else
#define CACHE_SIZE (0)
#endif
#define MAX_CACHE_MPZ_LIMBS (64)
#define MAX_CACHE_MPFR_BITS (1024)

typedef struct {
    MPZ_Object *gmpympzcache[CACHE_SIZE+1];
    int in_gmpympzcache;

    XMPZ_Object *gmpyxmpzcache[CACHE_SIZE+1];
    int in_gmpyxmpzcache;

    MPQ_Object *gmpympqcache[CACHE_SIZE+1];
    int in_gmpympqcache;

    MPFR_Object *gmpympfrcache[CACHE_SIZE+1];
    int in_gmpympfrcache;

    MPC_Object *gmpympccache[CACHE_SIZE+1];
    int in_gmpympccache;
} gmpy_global;

#if !defined(_MSC_VER)
#  define _Py_thread_local _Thread_local
#else
#  define _Py_thread_local __declspec(thread)
#endif

_Py_thread_local gmpy_global global = {
    .in_gmpympzcache = 0,
    .in_gmpyxmpzcache = 0,
    .in_gmpympqcache = 0,
    .in_gmpympfrcache = 0,
    .in_gmpympccache = 0,
};

/* Support for context manager using context vars.
 */

static PyObject *current_context_var = NULL;

/* Define gmpy2 specific errors for mpfr and mpc data types. No change will
 * be made the exceptions raised by mpz, xmpz, and mpq.
 */

static PyObject *GMPyExc_GmpyError = NULL;
static PyObject *GMPyExc_DivZero = NULL;
static PyObject *GMPyExc_Inexact = NULL;
static PyObject *GMPyExc_Invalid = NULL;
static PyObject *GMPyExc_Overflow = NULL;
static PyObject *GMPyExc_Underflow = NULL;
static PyObject *GMPyExc_Erange = NULL;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * End of global data declarations.                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* The code for object creation, deletion, and caching is in gmpy_cache.c. */

#include "gmpy2_cache.c"

/* Miscellaneous helper functions and simple methods are in gmpy_misc.c. */

#include "gmpy2_misc.c"

/* Support for conversion to/from binary representation. */

#include "gmpy2_binary.c"

/* Support for conversions to/from numeric types. */

#include "gmpy2_convert.c"
#include "gmpy2_convert_utils.c"
#include "gmpy2_convert_gmp.c"
#include "gmpy2_convert_mpfr.c"
#include "gmpy2_convert_mpc.c"

/* Support for random numbers. */

#include "gmpy2_random.c"

/* Support for Lucas sequences. */

#include "gmpy_mpz_lucas.c"

/* Support for probable-prime tests. */

#include "gmpy_mpz_prp.c"

/* Include helper functions for mpmath. */

#include "gmpy2_mpmath.c"

#include "gmpy2_mpz_divmod.c"
#include "gmpy2_mpz_divmod2exp.c"
#include "gmpy2_mpz_pack.c"
#include "gmpy2_mpz_bitops.c"
#include "gmpy2_xmpz_inplace.c"

/* Begin includes of refactored code. */

#include "gmpy2_abs.c"
#include "gmpy2_add.c"
#include "gmpy2_divmod.c"
#include "gmpy2_floordiv.c"
#include "gmpy2_minus.c"
#include "gmpy2_mod.c"
#include "gmpy2_mul.c"
#include "gmpy2_plus.c"
#include "gmpy2_pow.c"
#include "gmpy2_sub.c"
#include "gmpy2_truediv.c"
#include "gmpy2_math.c"
#include "gmpy2_const.c"
#include "gmpy2_square.c"
#include "gmpy2_format.c"
#include "gmpy2_hash.c"
#include "gmpy2_fused.c"
#include "gmpy2_muldiv_2exp.c"
#include "gmpy2_predicate.c"
#include "gmpy2_sign.c"
#include "gmpy2_richcompare.c"
#include "gmpy2_cmp.c"
#include "gmpy2_mpc_misc.c"
#include "gmpy2_mpfr_misc.c"
#include "gmpy2_mpq_misc.c"
#include "gmpy2_mpz_misc.c"
#include "gmpy2_xmpz_misc.c"
#include "gmpy2_xmpz_limbs.c"

#ifdef VECTOR
#include "gmpy2_vector.c"
#endif

/* Include gmpy_context last to avoid adding doc names to .h files. */

#include "gmpy2_mpz.c"
#include "gmpy2_xmpz.c"
#include "gmpy2_mpq.c"
#include "gmpy2_mpfr.c"
#include "gmpy2_mpc.c"

#include "gmpy2_context.c"

static PyMethodDef Pygmpy_methods [] =
{
    { "add", GMPy_Context_Add, METH_VARARGS, GMPy_doc_function_add },
    { "bit_clear", GMPy_MPZ_bit_clear_function, METH_VARARGS, doc_bit_clear_function },
    { "bit_count", GMPy_MPZ_bit_count, METH_O, doc_bit_count },
    { "bit_flip", GMPy_MPZ_bit_flip_function, METH_VARARGS, doc_bit_flip_function },
    { "bit_length", GMPy_MPZ_bit_length_function, METH_O, doc_bit_length_function },
    { "bit_mask", GMPy_MPZ_bit_mask, METH_O, doc_bit_mask },
    { "bit_scan0", (PyCFunction)GMPy_MPZ_bit_scan0_function, METH_FASTCALL, doc_bit_scan0_function },
    { "bit_scan1", (PyCFunction)GMPy_MPZ_bit_scan1_function, METH_FASTCALL, doc_bit_scan1_function },
    { "bit_set", GMPy_MPZ_bit_set_function, METH_VARARGS, doc_bit_set_function },
    { "bit_test", (PyCFunction)GMPy_MPZ_bit_test_function, METH_FASTCALL, doc_bit_test_function },
    { "bincoef", (PyCFunction)GMPy_MPZ_Function_Bincoef, METH_FASTCALL, GMPy_doc_mpz_function_bincoef },
    { "cmp", GMPy_MPANY_cmp, METH_VARARGS, GMPy_doc_mpany_cmp },
    { "cmp_abs", GMPy_MPANY_cmp_abs, METH_VARARGS, GMPy_doc_mpany_cmp_abs },
    { "comb", (PyCFunction)GMPy_MPZ_Function_Bincoef, METH_FASTCALL, GMPy_doc_mpz_function_comb },
    { "c_div", GMPy_MPZ_c_div, METH_VARARGS, doc_c_div },
    { "c_div_2exp", GMPy_MPZ_c_div_2exp, METH_VARARGS, doc_c_div_2exp },
    { "c_divmod", GMPy_MPZ_c_divmod, METH_VARARGS, doc_c_divmod },
    { "c_divmod_2exp", GMPy_MPZ_c_divmod_2exp, METH_VARARGS, doc_c_divmod_2exp },
    { "c_mod", GMPy_MPZ_c_mod, METH_VARARGS, doc_c_mod },
    { "c_mod_2exp", GMPy_MPZ_c_mod_2exp, METH_VARARGS, doc_c_mod_2exp },
    { "denom", GMPy_MPQ_Function_Denom, METH_O, GMPy_doc_mpq_function_denom },
    { "digits", GMPy_Context_Digits, METH_VARARGS, GMPy_doc_context_digits },
    { "div", GMPy_Context_TrueDiv, METH_VARARGS, GMPy_doc_truediv },
    { "divexact", (PyCFunction)GMPy_MPZ_Function_Divexact, METH_FASTCALL, GMPy_doc_mpz_function_divexact },
    { "divm", (PyCFunction)GMPy_MPZ_Function_Divm, METH_FASTCALL, GMPy_doc_mpz_function_divm },
    { "double_fac", GMPy_MPZ_Function_DoubleFac, METH_O, GMPy_doc_mpz_function_double_fac },
    { "fac", GMPy_MPZ_Function_Fac, METH_O, GMPy_doc_mpz_function_fac },
    { "fib", GMPy_MPZ_Function_Fib, METH_O, GMPy_doc_mpz_function_fib },
    { "fib2", GMPy_MPZ_Function_Fib2, METH_O, GMPy_doc_mpz_function_fib2 },
    { "floor_div", GMPy_Context_FloorDiv, METH_VARARGS, GMPy_doc_floordiv },
    { "from_binary", GMPy_MPANY_From_Binary, METH_O, doc_from_binary },
    { "f_div", GMPy_MPZ_f_div, METH_VARARGS, doc_f_div },
    { "f_div_2exp", GMPy_MPZ_f_div_2exp, METH_VARARGS, doc_f_div_2exp },
    { "f_divmod", GMPy_MPZ_f_divmod, METH_VARARGS, doc_f_divmod },
    { "f_divmod_2exp", GMPy_MPZ_f_divmod_2exp, METH_VARARGS, doc_f_divmod_2exp },
    { "f_mod", GMPy_MPZ_f_mod, METH_VARARGS, doc_f_mod },
    { "f_mod_2exp", GMPy_MPZ_f_mod_2exp, METH_VARARGS, doc_f_mod_2exp },
    { "gcd", (PyCFunction)GMPy_MPZ_Function_GCD, METH_FASTCALL, GMPy_doc_mpz_function_gcd },
    { "gcdext", (PyCFunction)GMPy_MPZ_Function_GCDext, METH_FASTCALL, GMPy_doc_mpz_function_gcdext },
    { "hamdist", GMPy_MPZ_hamdist, METH_VARARGS, doc_hamdist },
    { "invert", (PyCFunction)GMPy_MPZ_Function_Invert, METH_FASTCALL, GMPy_doc_mpz_function_invert },
    { "iroot", (PyCFunction)GMPy_MPZ_Function_Iroot, METH_FASTCALL, GMPy_doc_mpz_function_iroot },
    { "iroot_rem", (PyCFunction)GMPy_MPZ_Function_IrootRem, METH_FASTCALL, GMPy_doc_mpz_function_iroot_rem },
    { "isqrt", GMPy_MPZ_Function_Isqrt, METH_O, GMPy_doc_mpz_function_isqrt },
    { "isqrt_rem", GMPy_MPZ_Function_IsqrtRem, METH_O, GMPy_doc_mpz_function_isqrt_rem },
    { "is_bpsw_prp", GMPY_mpz_is_bpsw_prp, METH_VARARGS, doc_mpz_is_bpsw_prp },
    { "is_congruent", (PyCFunction)GMPy_MPZ_Function_IsCongruent, METH_FASTCALL, GMPy_doc_mpz_function_is_congruent },
    { "is_divisible", (PyCFunction)GMPy_MPZ_Function_IsDivisible, METH_FASTCALL, GMPy_doc_mpz_function_is_divisible },
    { "is_even", GMPy_MPZ_Function_IsEven, METH_O, GMPy_doc_mpz_function_is_even },
    { "is_euler_prp", GMPY_mpz_is_euler_prp, METH_VARARGS, doc_mpz_is_euler_prp },
    { "is_extra_strong_lucas_prp", GMPY_mpz_is_extrastronglucas_prp, METH_VARARGS, doc_mpz_is_extrastronglucas_prp },
    { "is_fermat_prp", GMPY_mpz_is_fermat_prp, METH_VARARGS, doc_mpz_is_fermat_prp },
    { "is_fibonacci_prp", GMPY_mpz_is_fibonacci_prp, METH_VARARGS, doc_mpz_is_fibonacci_prp },
    { "is_lucas_prp", GMPY_mpz_is_lucas_prp, METH_VARARGS, doc_mpz_is_lucas_prp },
    { "is_odd", GMPy_MPZ_Function_IsOdd, METH_O, GMPy_doc_mpz_function_is_odd },
    { "is_power", GMPy_MPZ_Function_IsPower, METH_O, GMPy_doc_mpz_function_is_power },
    { "is_prime", (PyCFunction)GMPy_MPZ_Function_IsPrime, METH_FASTCALL, GMPy_doc_mpz_function_is_prime },
    { "is_probab_prime", (PyCFunction)GMPy_MPZ_Function_IsProbabPrime, METH_FASTCALL, GMPy_doc_mpz_function_is_probab_prime },
    { "is_selfridge_prp", GMPY_mpz_is_selfridge_prp, METH_VARARGS, doc_mpz_is_selfridge_prp },
    { "is_square", GMPy_MPZ_Function_IsSquare, METH_O, GMPy_doc_mpz_function_is_square },
    { "is_strong_prp", GMPY_mpz_is_strong_prp, METH_VARARGS, doc_mpz_is_strong_prp },
    { "is_strong_bpsw_prp", GMPY_mpz_is_strongbpsw_prp, METH_VARARGS, doc_mpz_is_strongbpsw_prp },
    { "is_strong_lucas_prp", GMPY_mpz_is_stronglucas_prp, METH_VARARGS, doc_mpz_is_stronglucas_prp },
    { "is_strong_selfridge_prp", GMPY_mpz_is_strongselfridge_prp, METH_VARARGS, doc_mpz_is_strongselfridge_prp },
    { "jacobi", (PyCFunction)GMPy_MPZ_Function_Jacobi, METH_FASTCALL, GMPy_doc_mpz_function_jacobi },
    { "kronecker", (PyCFunction)GMPy_MPZ_Function_Kronecker, METH_FASTCALL, GMPy_doc_mpz_function_kronecker },
    { "lcm", (PyCFunction)GMPy_MPZ_Function_LCM, METH_FASTCALL, GMPy_doc_mpz_function_lcm },
    { "legendre", (PyCFunction)GMPy_MPZ_Function_Legendre, METH_FASTCALL, GMPy_doc_mpz_function_legendre },
    { "license", GMPy_get_license, METH_NOARGS, GMPy_doc_license },
    { "lucas", GMPy_MPZ_Function_Lucas, METH_O, GMPy_doc_mpz_function_lucas },
    { "lucasu", GMPY_mpz_lucasu, METH_VARARGS, doc_mpz_lucasu },
    { "lucasu_mod", GMPY_mpz_lucasu_mod, METH_VARARGS, doc_mpz_lucasu_mod },
    { "lucasv", GMPY_mpz_lucasv, METH_VARARGS, doc_mpz_lucasv },
    { "lucasv_mod", GMPY_mpz_lucasv_mod, METH_VARARGS, doc_mpz_lucasv_mod },
    { "lucas2", GMPy_MPZ_Function_Lucas2, METH_O, GMPy_doc_mpz_function_lucas2 },
    { "mod", GMPy_Context_Mod, METH_VARARGS, GMPy_doc_mod },
    { "mp_version", GMPy_get_mp_version, METH_NOARGS, GMPy_doc_mp_version },
    { "mp_limbsize", GMPy_get_mp_limbsize, METH_NOARGS, GMPy_doc_mp_limbsize },
    { "mpc_version", GMPy_get_mpc_version, METH_NOARGS, GMPy_doc_mpc_version },
    { "mpfr_version", GMPy_get_mpfr_version, METH_NOARGS, GMPy_doc_mpfr_version },
    { "mpq_from_old_binary", GMPy_MPQ_From_Old_Binary, METH_O, doc_mpq_from_old_binary },
    { "mpz_from_old_binary", GMPy_MPZ_From_Old_Binary, METH_O, doc_mpz_from_old_binary },
    { "mpz_random", GMPy_MPZ_random_Function, METH_VARARGS, GMPy_doc_mpz_random_function },
    { "mpz_rrandomb", GMPy_MPZ_rrandomb_Function, METH_VARARGS, GMPy_doc_mpz_rrandomb_function },
    { "mpz_urandomb", GMPy_MPZ_urandomb_Function, METH_VARARGS, GMPy_doc_mpz_urandomb_function },
    { "mul", GMPy_Context_Mul, METH_VARARGS, GMPy_doc_function_mul },
    { "multi_fac", (PyCFunction)GMPy_MPZ_Function_MultiFac, METH_FASTCALL, GMPy_doc_mpz_function_multi_fac },
    { "next_prime", GMPy_MPZ_Function_NextPrime, METH_O, GMPy_doc_mpz_function_next_prime },
#if (__GNU_MP_VERSION > 6) || (__GNU_MP_VERSION == 6 &&  __GNU_MP_VERSION_MINOR >= 3)
    { "prev_prime", GMPy_MPZ_Function_PrevPrime, METH_O, GMPy_doc_mpz_function_prev_prime },
#endif
    { "numer", GMPy_MPQ_Function_Numer, METH_O, GMPy_doc_mpq_function_numer },
    { "num_digits", (PyCFunction)GMPy_MPZ_Function_NumDigits, METH_FASTCALL, GMPy_doc_mpz_function_num_digits },
    { "pack", GMPy_MPZ_pack, METH_VARARGS, doc_pack },
    { "popcount", GMPy_MPZ_popcount, METH_O, doc_popcount },
    { "powmod", GMPy_Integer_PowMod, METH_VARARGS, GMPy_doc_integer_powmod },
    { "powmod_base_list", GMPy_Integer_PowMod_Base_List, METH_VARARGS, GMPy_doc_integer_powmod_base_list },
    { "powmod_exp_list", GMPy_Integer_PowMod_Exp_List, METH_VARARGS, GMPy_doc_integer_powmod_exp_list },
    { "powmod_sec", GMPy_Integer_PowMod_Sec, METH_VARARGS, GMPy_doc_integer_powmod_sec },
    { "primorial", GMPy_MPZ_Function_Primorial, METH_O, GMPy_doc_mpz_function_primorial },
    { "qdiv", GMPy_MPQ_Function_Qdiv, METH_VARARGS, GMPy_doc_function_qdiv },
    { "remove", (PyCFunction)GMPy_MPZ_Function_Remove, METH_FASTCALL, GMPy_doc_mpz_function_remove },
    { "random_state", GMPy_RandomState_Factory, METH_VARARGS, GMPy_doc_random_state_factory },
    { "sign", GMPy_Context_Sign, METH_O, GMPy_doc_function_sign },
    { "square", GMPy_Context_Square, METH_O, GMPy_doc_function_square },
    { "sub", GMPy_Context_Sub, METH_VARARGS, GMPy_doc_sub },
    { "to_binary", GMPy_MPANY_To_Binary, METH_O, doc_to_binary },
    { "t_div", GMPy_MPZ_t_div, METH_VARARGS, doc_t_div },
    { "t_div_2exp", GMPy_MPZ_t_div_2exp, METH_VARARGS, doc_t_div_2exp },
    { "t_divmod", GMPy_MPZ_t_divmod, METH_VARARGS, doc_t_divmod },
    { "t_divmod_2exp", GMPy_MPZ_t_divmod_2exp, METH_VARARGS, doc_t_divmod_2exp },
    { "t_mod", GMPy_MPZ_t_mod, METH_VARARGS, doc_t_mod },
    { "t_mod_2exp", GMPy_MPZ_t_mod_2exp, METH_VARARGS, doc_t_mod_2exp },
    { "unpack", GMPy_MPZ_unpack, METH_VARARGS, doc_unpack },
    { "version", GMPy_get_version, METH_NOARGS, GMPy_doc_version },
    { "xbit_mask", GMPy_XMPZ_Function_XbitMask, METH_O, GMPy_doc_xmpz_function_xbit_mask },
    { "_mpmath_normalize", (PyCFunction)Pympz_mpmath_normalize_fast, METH_FASTCALL, doc_mpmath_normalizeg },
    { "_mpmath_create", (PyCFunction)Pympz_mpmath_create_fast, METH_FASTCALL, doc_mpmath_create },

    { "acos", GMPy_Context_Acos, METH_O, GMPy_doc_function_acos },
    { "acosh", GMPy_Context_Acosh, METH_O, GMPy_doc_function_acosh },
    { "ai", GMPy_Context_Ai, METH_O, GMPy_doc_function_ai },
    { "agm", GMPy_Context_AGM, METH_VARARGS, GMPy_doc_function_agm },
    { "asin", GMPy_Context_Asin, METH_O, GMPy_doc_function_asin },
    { "asinh", GMPy_Context_Asinh, METH_O, GMPy_doc_function_asinh },
    { "atan", GMPy_Context_Atan, METH_O, GMPy_doc_function_atan },
    { "atanh", GMPy_Context_Atanh, METH_O, GMPy_doc_function_atanh },
    { "atan2", GMPy_Context_Atan2, METH_VARARGS, GMPy_doc_function_atan2 },
    { "can_round", GMPy_MPFR_Can_Round, METH_VARARGS, GMPy_doc_mpfr_can_round },
    { "cbrt", GMPy_Context_Cbrt, METH_O, GMPy_doc_function_cbrt },
    { "ceil", GMPy_Context_Ceil, METH_O, GMPy_doc_function_ceil },
    { "check_range", GMPy_Context_CheckRange, METH_O, GMPy_doc_function_check_range },
    { "const_catalan", (PyCFunction)GMPy_Function_Const_Catalan, METH_VARARGS | METH_KEYWORDS, GMPy_doc_function_const_catalan },
    { "const_euler", (PyCFunction)GMPy_Function_Const_Euler, METH_VARARGS | METH_KEYWORDS, GMPy_doc_function_const_euler },
    { "const_log2", (PyCFunction)GMPy_Function_Const_Log2, METH_VARARGS | METH_KEYWORDS, GMPy_doc_function_const_log2 },
    { "const_pi", (PyCFunction)GMPy_Function_Const_Pi, METH_VARARGS | METH_KEYWORDS, GMPy_doc_function_const_pi },
    { "copy_sign", GMPy_MPFR_copy_sign, METH_VARARGS, GMPy_doc_mpfr_copy_sign },
    { "cos", GMPy_Context_Cos, METH_O, GMPy_doc_function_cos },
    { "cosh", GMPy_Context_Cosh, METH_O, GMPy_doc_function_cosh },
    { "cot", GMPy_Context_Cot, METH_O, GMPy_doc_function_cot },
    { "coth", GMPy_Context_Coth, METH_O, GMPy_doc_function_coth },
    { "csc", GMPy_Context_Csc, METH_O, GMPy_doc_function_csc },
    { "csch", GMPy_Context_Csch, METH_O, GMPy_doc_function_csch },
    { "degrees", GMPy_Context_Degrees, METH_O, GMPy_doc_function_degrees },
    { "digamma", GMPy_Context_Digamma, METH_O, GMPy_doc_function_digamma },
    { "div_2exp", GMPy_Context_Div_2exp, METH_VARARGS, GMPy_doc_function_div_2exp },
    { "eint", GMPy_Context_Eint, METH_O, GMPy_doc_function_eint },
    { "erf", GMPy_Context_Erf, METH_O, GMPy_doc_function_erf },
    { "erfc", GMPy_Context_Erfc, METH_O, GMPy_doc_function_erfc },
    { "exp", GMPy_Context_Exp, METH_O, GMPy_doc_function_exp },
    { "expm1", GMPy_Context_Expm1, METH_O, GMPy_doc_function_expm1 },
    { "exp10", GMPy_Context_Exp10, METH_O, GMPy_doc_function_exp10 },
    { "exp2", GMPy_Context_Exp2, METH_O, GMPy_doc_function_exp2 },
    { "f2q", GMPy_Context_F2Q, METH_VARARGS, GMPy_doc_function_f2q },
    { "factorial", GMPy_Context_Factorial, METH_O, GMPy_doc_function_factorial },
    { "floor", GMPy_Context_Floor, METH_O, GMPy_doc_function_floor },
    { "fma", GMPy_Context_FMA, METH_VARARGS, GMPy_doc_function_fma },
    { "fms", GMPy_Context_FMS, METH_VARARGS, GMPy_doc_function_fms },
    { "fmma", GMPy_Context_FMMA, METH_VARARGS, GMPy_doc_function_fmma },
    { "fmms", GMPy_Context_FMMS, METH_VARARGS, GMPy_doc_function_fmms },
    { "fmod", GMPy_Context_Fmod, METH_VARARGS, GMPy_doc_function_fmod },
    { "frac", GMPy_Context_Frac, METH_O, GMPy_doc_function_frac },
    { "free_cache", GMPy_MPFR_Free_Cache, METH_NOARGS, GMPy_doc_mpfr_free_cache },
    { "frexp", GMPy_Context_Frexp, METH_O, GMPy_doc_function_frexp },
    { "fsum", GMPy_Context_Fsum, METH_O, GMPy_doc_function_fsum },
    { "gamma", GMPy_Context_Gamma, METH_O, GMPy_doc_function_gamma },
    { "gamma_inc", GMPy_Context_Gamma_Inc, METH_VARARGS, GMPy_doc_function_gamma_inc },
    { "get_context", GMPy_CTXT_Get, METH_NOARGS, GMPy_doc_get_context },
    { "get_emax_max", GMPy_MPFR_get_emax_max, METH_NOARGS, GMPy_doc_mpfr_get_emax_max },
    { "get_emin_min", GMPy_MPFR_get_emin_min, METH_NOARGS, GMPy_doc_mpfr_get_emin_min },
    { "get_exp", GMPy_MPFR_get_exp, METH_O, GMPy_doc_mpfr_get_exp },
    { "get_max_precision", GMPy_MPFR_get_max_precision, METH_NOARGS, GMPy_doc_mpfr_get_max_precision },
    { "hypot", GMPy_Context_Hypot, METH_VARARGS, GMPy_doc_function_hypot },
    { "ieee", (PyCFunction)GMPy_CTXT_ieee, METH_VARARGS | METH_KEYWORDS, GMPy_doc_context_ieee },
    { "inf", GMPy_MPFR_set_inf, METH_VARARGS, GMPy_doc_mpfr_set_inf },
    { "is_finite", GMPy_Context_Is_Finite, METH_O, GMPy_doc_function_is_finite },
    { "is_infinite", GMPy_Context_Is_Infinite, METH_O, GMPy_doc_function_is_infinite },
    { "is_integer", GMPy_Context_Is_Integer, METH_O, GMPy_doc_function_is_integer },
    { "is_lessgreater", GMPy_Context_Is_LessGreater, METH_VARARGS, GMPy_doc_function_is_lessgreater },
    { "is_nan", GMPy_Context_Is_NAN, METH_O, GMPy_doc_function_is_nan },
    { "is_regular", GMPy_Context_Is_Regular, METH_O, GMPy_doc_function_is_regular },
    { "is_signed", GMPy_Context_Is_Signed, METH_O, GMPy_doc_function_is_signed },
    { "is_unordered", GMPy_Context_Is_Unordered, METH_VARARGS, GMPy_doc_function_is_unordered },
    { "is_zero", GMPy_Context_Is_Zero, METH_O, GMPy_doc_function_is_zero },
    { "jn", GMPy_Context_Jn, METH_VARARGS, GMPy_doc_function_jn },
    { "j0", GMPy_Context_J0, METH_O, GMPy_doc_function_j0 },
    { "j1", GMPy_Context_J1, METH_O, GMPy_doc_function_j1 },
    { "lgamma", GMPy_Context_Lgamma, METH_O, GMPy_doc_function_lgamma },
    { "li2", GMPy_Context_Li2, METH_O, GMPy_doc_function_li2 },
    { "lngamma", GMPy_Context_Lngamma, METH_O, GMPy_doc_function_lngamma },
    { "local_context", (PyCFunction)GMPy_CTXT_Local, METH_VARARGS | METH_KEYWORDS, GMPy_doc_local_context },
    { "log", GMPy_Context_Log, METH_O, GMPy_doc_function_log },
    { "log1p", GMPy_Context_Log1p, METH_O, GMPy_doc_function_log1p },
    { "log10", GMPy_Context_Log10, METH_O, GMPy_doc_function_log10 },
    { "log2", GMPy_Context_Log2, METH_O, GMPy_doc_function_log2 },
    { "maxnum", GMPy_Context_Maxnum, METH_VARARGS, GMPy_doc_function_maxnum },
    { "minnum", GMPy_Context_Minnum, METH_VARARGS, GMPy_doc_function_minnum },
    { "modf", GMPy_Context_Modf, METH_O, GMPy_doc_function_modf },
    { "mpfr_from_old_binary", GMPy_MPFR_From_Old_Binary, METH_O, doc_mpfr_from_old_binary },
    { "mpfr_random", GMPy_MPFR_random_Function, METH_VARARGS, GMPy_doc_mpfr_random_function },
    { "mpfr_grandom", GMPy_MPFR_grandom_Function, METH_VARARGS, GMPy_doc_mpfr_grandom_function },
    { "mpfr_nrandom", GMPy_MPFR_nrandom_Function, METH_VARARGS, GMPy_doc_mpfr_nrandom_function },
    { "mul_2exp", GMPy_Context_Mul_2exp, METH_VARARGS, GMPy_doc_function_mul_2exp },
    { "nan", GMPy_MPFR_set_nan, METH_NOARGS, GMPy_doc_mpfr_set_nan },
    { "next_above", GMPy_Context_NextAbove, METH_O, GMPy_doc_function_next_above },
    { "next_below", GMPy_Context_NextBelow, METH_O, GMPy_doc_function_next_below },
    { "next_toward", GMPy_Context_NextToward, METH_VARARGS, GMPy_doc_function_next_toward },
    { "radians", GMPy_Context_Radians, METH_O, GMPy_doc_function_radians },
    { "rec_sqrt", GMPy_Context_RecSqrt, METH_O, GMPy_doc_function_rec_sqrt },
    { "reldiff", GMPy_Context_RelDiff, METH_VARARGS, GMPy_doc_function_reldiff },
    { "remainder", GMPy_Context_Remainder, METH_VARARGS, GMPy_doc_function_remainder },
    { "remquo", GMPy_Context_RemQuo, METH_VARARGS, GMPy_doc_function_remquo },
    { "rint", GMPy_Context_Rint, METH_O, GMPy_doc_function_rint },
    { "rint_ceil", GMPy_Context_RintCeil, METH_O, GMPy_doc_function_rint_ceil },
    { "rint_floor", GMPy_Context_RintFloor, METH_O, GMPy_doc_function_rint_floor },
    { "rint_round", GMPy_Context_RintRound, METH_O, GMPy_doc_function_rint_round },
    { "rint_trunc", GMPy_Context_RintTrunc, METH_O, GMPy_doc_function_rint_trunc },
    { "root", GMPy_Context_Root, METH_VARARGS, GMPy_doc_function_root },
    { "rootn", GMPy_Context_Rootn, METH_VARARGS, GMPy_doc_function_rootn },
    { "round_away", GMPy_Context_RoundAway, METH_O, GMPy_doc_function_round_away },
    { "round2", GMPy_Context_Round2, METH_VARARGS, GMPy_doc_function_round2 },
    { "sec", GMPy_Context_Sec, METH_O, GMPy_doc_function_sec },
    { "sech", GMPy_Context_Sech, METH_O, GMPy_doc_function_sech },
    { "set_context", GMPy_CTXT_Set, METH_O, GMPy_doc_set_context },
    { "set_exp", GMPy_MPFR_set_exp, METH_VARARGS, GMPy_doc_mpfr_set_exp },
    { "set_sign", GMPy_MPFR_set_sign, METH_VARARGS, GMPy_doc_mpfr_set_sign },
    { "sin", GMPy_Context_Sin, METH_O, GMPy_doc_function_sin },
    { "sin_cos", GMPy_Context_Sin_Cos, METH_O, GMPy_doc_function_sin_cos },
    { "sinh", GMPy_Context_Sinh, METH_O, GMPy_doc_function_sinh },
    { "sinh_cosh", GMPy_Context_Sinh_Cosh, METH_O, GMPy_doc_function_sinh_cosh },
    { "sqrt", GMPy_Context_Sqrt, METH_O, GMPy_doc_function_sqrt },
    { "tan", GMPy_Context_Tan, METH_O, GMPy_doc_function_tan },
    { "tanh", GMPy_Context_Tanh, METH_O, GMPy_doc_function_tanh },
    { "trunc", GMPy_Context_Trunc, METH_O, GMPy_doc_function_trunc},
#ifdef VECTOR
    { "vector", GMPy_Context_Vector, METH_O, GMPy_doc_function_vector},
    { "vector2", GMPy_Context_Vector2, METH_VARARGS, GMPy_doc_function_vector2},
#endif
    { "yn", GMPy_Context_Yn, METH_VARARGS, GMPy_doc_function_yn },
    { "y0", GMPy_Context_Y0, METH_O, GMPy_doc_function_y0 },
    { "y1", GMPy_Context_Y1, METH_O, GMPy_doc_function_y1 },
    { "zero", GMPy_MPFR_set_zero, METH_VARARGS, GMPy_doc_mpfr_set_zero },
    { "zeta", GMPy_Context_Zeta, METH_O, GMPy_doc_function_zeta },

    { "mpc_random", GMPy_MPC_random_Function, METH_VARARGS, GMPy_doc_mpc_random_function },
    { "norm", GMPy_Context_Norm, METH_O, GMPy_doc_function_norm },
    { "polar", GMPy_Context_Polar, METH_O, GMPy_doc_function_polar },
    { "phase", GMPy_Context_Phase, METH_O, GMPy_doc_function_phase },
    { "proj", GMPy_Context_Proj, METH_O, GMPy_doc_function_proj },
    { "root_of_unity", GMPy_Context_Root_Of_Unity, METH_VARARGS, GMPy_doc_function_root_of_unity },
    { "rect", GMPy_Context_Rect, METH_VARARGS, GMPy_doc_function_rect },
    { NULL, NULL, 1}
};

static char _gmpy_docs[] = "gmpy2 " GMPY_VERSION
" - General Multiple-precision arithmetic for Python\n"
"\n"
"gmpy2 supports several multiple-precision libraries. Integer and\n"
"rational arithmetic is provided by the GMP library. Real floating-\n"
"point arithmetic is provided by the MPFR library. Complex floating-\n"
"point arithmetic is provided by the MPC library.\n"
"\n"
"The integer type 'mpz' has comparable functionality to Python's\n"
"builtin integers, but is faster for operations on large numbers.\n"
"A wide variety of additional functions are provided:\n"
"      - bit manipulations\n"
"      - GCD, Extended GCD, LCM\n"
"      - Fibonacci and Lucas sequences\n"
"      - primality testing\n"
"      - powers and integer Nth roots\n"
"\n"
"The rational type 'mpq' is equivalent to Python's fractions\n"
"module, but is faster.\n"
"\n"
"The real type 'mpfr' and complex type 'mpc' provide multiple-\n"
"precision real and complex numbers with user-definable precision,\n"
"rounding, and exponent range. All the advanced functions from the\n"
"MPFR and MPC libraries are available.\n\
";

void
gmp_abort_handler(int i)
{
    if (__GNU_MP_VERSION == 6 &&  __GNU_MP_VERSION_MINOR >= 3) {
        if (gmp_errno & GMP_ERROR_DIVISION_BY_ZERO) {
            printf("gmp: divide by zero\n");
        }
        else if (gmp_errno & GMP_ERROR_SQRT_OF_NEGATIVE) {
            printf("gmp: root of negative\n");
        }
        else {
            printf("gmp: overflow in mpz type\n");
        }
    }
    abort();
}

static int
gmpy_exec(PyObject *gmpy_module)
{
    PyObject *result = NULL;
    PyObject *namespace = NULL;
    PyObject *copy_reg_module = NULL;
    PyObject *temp = NULL;
    PyObject *numbers_module = NULL;
    PyObject* xmpz = NULL;
    PyObject* limb_size = NULL;

#ifndef STATIC
    static void *GMPy_C_API[GMPy_API_pointers];
    PyObject *c_api_object;
#endif

    /* Validate the sizes of the various typedef'ed integer types. */

    if (sizeof(mpfr_prec_t) != sizeof(long)) {
        /* LCOV_EXCL_START */
        SYSTEM_ERROR("Size of mpfr_prec_t and long not compatible");
        return -1;;
        /* LCOV_EXCL_STOP */
    }

    if (sizeof(mpfr_exp_t) != sizeof(long)) {
        /* LCOV_EXCL_START */
        SYSTEM_ERROR("Size of mpfr_exp_t and long not compatible");
        return -1;;
        /* LCOV_EXCL_STOP */
    }

    /* Initialize the types. */
    if (PyType_Ready(&MPZ_Type) < 0) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    if (PyType_Ready(&MPQ_Type) < 0) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    if (PyType_Ready(&XMPZ_Type) < 0) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    if (PyType_Ready(&GMPy_Iter_Type) < 0) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    if (PyType_Ready(&MPFR_Type) < 0) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    if (PyType_Ready(&CTXT_Type) < 0) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    if (PyType_Ready(&MPC_Type) < 0) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    if (PyType_Ready(&RandomState_Type) < 0) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }

    /* Initialize exceptions. */
    GMPyExc_GmpyError = PyErr_NewException("gmpy2.gmpy2Error", PyExc_ArithmeticError, NULL);
    if (!GMPyExc_GmpyError) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }

    GMPyExc_Erange = PyErr_NewException("gmpy2.RangeError", GMPyExc_GmpyError, NULL);
    if (!GMPyExc_Erange) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }

    GMPyExc_Inexact = PyErr_NewException("gmpy2.InexactResultError", GMPyExc_GmpyError, NULL);
    if (!GMPyExc_Inexact) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }

    GMPyExc_Overflow = PyErr_NewException("gmpy2.OverflowResultError", GMPyExc_Inexact, NULL);
    if (!GMPyExc_Overflow) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }

    GMPyExc_Underflow = PyErr_NewException("gmpy2.UnderflowResultError", GMPyExc_Inexact, NULL);
    if (!GMPyExc_Underflow) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }

    temp = PyTuple_Pack(2, GMPyExc_GmpyError, PyExc_ValueError);
    if (!temp) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    GMPyExc_Invalid = PyErr_NewException("gmpy2.InvalidOperationError", temp, NULL);
    Py_DECREF(temp);
    if (!GMPyExc_Invalid) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }

    temp = PyTuple_Pack(2, GMPyExc_GmpyError, PyExc_ZeroDivisionError);
    if (!temp) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    GMPyExc_DivZero = PyErr_NewException("gmpy2.DivisionByZeroError", temp, NULL);
    Py_DECREF(temp);
    if (!GMPyExc_DivZero) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }

    /* Add the context type to the module namespace. */

    Py_INCREF(&CTXT_Type);
    PyModule_AddObject(gmpy_module, "context", (PyObject*)&CTXT_Type);

    /* Add the mpz type to the module namespace. */

    Py_INCREF(&MPZ_Type);
    PyModule_AddObject(gmpy_module, "mpz", (PyObject*)&MPZ_Type);

    /* Add the xmpz type to the module namespace. */

    Py_INCREF(&XMPZ_Type);
    PyModule_AddObject(gmpy_module, "xmpz", (PyObject*)&XMPZ_Type);

    xmpz = XMPZ_Type.tp_dict;
    limb_size = PyLong_FromSize_t(sizeof(mp_limb_t));
    PyDict_SetItemString(xmpz, "limb_size", limb_size);
    Py_DECREF(limb_size);

    /* Add the MPQ type to the module namespace. */

    Py_INCREF(&MPQ_Type);
    PyModule_AddObject(gmpy_module, "mpq", (PyObject*)&MPQ_Type);

    /* Add the MPFR type to the module namespace. */

    Py_INCREF(&MPFR_Type);
    PyModule_AddObject(gmpy_module, "mpfr", (PyObject*)&MPFR_Type);

    /* Add the MPC type to the module namespace. */

    Py_INCREF(&MPC_Type);
    PyModule_AddObject(gmpy_module, "mpc", (PyObject*)&MPC_Type);

    /* Initialize context var. */
    if (!(current_context_var = PyContextVar_New("gmpy2_context", NULL))) {
        return -1;
    }

    /* Add the constants for defining rounding modes. */
    if (PyModule_AddIntConstant(gmpy_module, "RoundToNearest", MPFR_RNDN) < 0) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    if (PyModule_AddIntConstant(gmpy_module, "RoundToZero", MPFR_RNDZ) < 0) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    if (PyModule_AddIntConstant(gmpy_module, "RoundUp", MPFR_RNDU) < 0) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    if (PyModule_AddIntConstant(gmpy_module, "RoundDown", MPFR_RNDD) < 0) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    if (PyModule_AddIntConstant(gmpy_module, "RoundAwayZero", MPFR_RNDA) < 0) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    if (PyModule_AddIntConstant(gmpy_module, "Default", GMPY_DEFAULT) < 0) {
        /* LCOV_EXCL_START */
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    if (PyModule_AddStringConstant(gmpy_module, "__version__", gmpy_version) < 0) {
        /* LCOV_EXCL_START */
        return -1;
        /* LCOV_EXCL_STOP */
    }

    /* Add the exceptions. */
    Py_INCREF(GMPyExc_DivZero);
    if (PyModule_AddObject(gmpy_module, "DivisionByZeroError", GMPyExc_DivZero) < 0) {
        /* LCOV_EXCL_START */
        Py_DECREF(GMPyExc_DivZero);
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    Py_INCREF(GMPyExc_Inexact);
    if (PyModule_AddObject(gmpy_module, "InexactResultError", GMPyExc_Inexact) < 0) {
        /* LCOV_EXCL_START */
        Py_DECREF(GMPyExc_Inexact);
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    Py_INCREF(GMPyExc_Invalid);
    if (PyModule_AddObject(gmpy_module, "InvalidOperationError", GMPyExc_Invalid) < 0 ) {
        /* LCOV_EXCL_START */
        Py_DECREF(GMPyExc_Invalid);
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    Py_INCREF(GMPyExc_Overflow);
    if (PyModule_AddObject(gmpy_module, "OverflowResultError", GMPyExc_Overflow) < 0) {
        /* LCOV_EXCL_START */
        Py_DECREF(GMPyExc_Overflow);
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    Py_INCREF(GMPyExc_Underflow);
    if (PyModule_AddObject(gmpy_module, "UnderflowResultError", GMPyExc_Underflow) < 0) {
        /* LCOV_EXCL_START */
        Py_DECREF(GMPyExc_Underflow);
        return -1;;
        /* LCOV_EXCL_STOP */
    }
    Py_INCREF(GMPyExc_Erange);
    if (PyModule_AddObject(gmpy_module, "RangeError", GMPyExc_Erange) < 0) {
        /* LCOV_EXCL_START */
        Py_DECREF(GMPyExc_Erange);
        return -1;;
        /* LCOV_EXCL_STOP */
    }

#ifdef SHARED
    /* Create the Capsule for the C-API. */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
    GMPy_C_API[MPZ_Type_NUM] = (void*)&MPZ_Type;
    GMPy_C_API[XMPZ_Type_NUM] = (void*)&XMPZ_Type;
    GMPy_C_API[MPQ_Type_NUM] = (void*)&MPQ_Type;
    GMPy_C_API[XMPQ_Type_NUM] = (void*)&MPQ_Type;
    GMPy_C_API[MPFR_Type_NUM] = (void*)&MPFR_Type;
    GMPy_C_API[XMPFR_Type_NUM] = (void*)&MPFR_Type;
    GMPy_C_API[MPC_Type_NUM] = (void*)&MPC_Type;
    GMPy_C_API[XMPC_Type_NUM] = (void*)&MPC_Type;
    GMPy_C_API[CTXT_Type_NUM] = (void*)&CTXT_Type;
    GMPy_C_API[RandomState_Type_NUM] = (void*)&RandomState_Type;

    GMPy_C_API[GMPy_MPZ_New_NUM] = (void*)GMPy_MPZ_New;
    GMPy_C_API[GMPy_MPZ_NewInit_NUM] = (void*)GMPy_MPZ_NewInit;
    GMPy_C_API[GMPy_MPZ_Dealloc_NUM] = (void*)GMPy_MPZ_Dealloc;
    GMPy_C_API[GMPy_MPZ_ConvertArg_NUM] = (void*)GMPy_MPZ_ConvertArg;

    GMPy_C_API[GMPy_XMPZ_New_NUM] = (void*)GMPy_XMPZ_New;
    GMPy_C_API[GMPy_XMPZ_NewInit_NUM] = (void*)GMPy_XMPZ_NewInit;
    GMPy_C_API[GMPy_XMPZ_Dealloc_NUM] = (void*)GMPy_XMPZ_Dealloc;

    GMPy_C_API[GMPy_MPQ_New_NUM] = (void*)GMPy_MPQ_New;
    GMPy_C_API[GMPy_MPQ_NewInit_NUM] = (void*)GMPy_MPQ_NewInit;
    GMPy_C_API[GMPy_MPQ_Dealloc_NUM] = (void*)GMPy_MPQ_Dealloc;
    GMPy_C_API[GMPy_MPQ_ConvertArg_NUM] = (void*)GMPy_MPQ_ConvertArg;

    GMPy_C_API[GMPy_MPFR_New_NUM] = (void*)GMPy_MPFR_New;
    GMPy_C_API[GMPy_MPFR_NewInit_NUM] = (void*)GMPy_MPFR_NewInit;
    GMPy_C_API[GMPy_MPFR_Dealloc_NUM] = (void*)GMPy_MPFR_Dealloc;
    GMPy_C_API[GMPy_MPFR_ConvertArg_NUM] = (void*)GMPy_MPFR_ConvertArg;

    GMPy_C_API[GMPy_MPC_New_NUM] = (void*)GMPy_MPC_New;
    GMPy_C_API[GMPy_MPC_NewInit_NUM] = (void*)GMPy_MPC_NewInit;
    GMPy_C_API[GMPy_MPC_Dealloc_NUM] = (void*)GMPy_MPC_Dealloc;
    GMPy_C_API[GMPy_MPC_ConvertArg_NUM] = (void*)GMPy_MPC_ConvertArg;
#pragma GCC diagnostic pop
    c_api_object = PyCapsule_New((void *)GMPy_C_API, "gmpy2._C_API", NULL);

    if (c_api_object != NULL) {
        PyModule_AddObject(gmpy_module, "_C_API", c_api_object);
    }
#endif

    /* Add support for pickling. */
    copy_reg_module = PyImport_ImportModule("copyreg");
    if (copy_reg_module) {
        char* enable_pickle =
            "def gmpy2_reducer(x): return (gmpy2.from_binary, (gmpy2.to_binary(x),))\n"
            "copyreg.pickle(gmpy2.mpz, gmpy2_reducer)\n"
            "copyreg.pickle(gmpy2.xmpz, gmpy2_reducer)\n"
            "copyreg.pickle(gmpy2.mpq, gmpy2_reducer)\n"
            "copyreg.pickle(gmpy2.mpfr, gmpy2_reducer)\n"
            "copyreg.pickle(gmpy2.mpc, gmpy2_reducer)\n";

        namespace = PyDict_New();
        PyDict_SetItemString(namespace, "copyreg", copy_reg_module);
        PyDict_SetItemString(namespace, "gmpy2", gmpy_module);
        result = PyRun_String(enable_pickle, Py_file_input, namespace, namespace);
        if (!result) {
            /* LCOV_EXCL_START */
            PyErr_Clear();
            /* LCOV_EXCL_STOP */
        }
        Py_DECREF(namespace);
        Py_DECREF(copy_reg_module);
        Py_XDECREF(result);
    }
    else {
        /* LCOV_EXCL_START */
        PyErr_Clear();
        /* LCOV_EXCL_STOP */
    }

    /* Register the gmpy2 types with the numeric tower. */

    numbers_module = PyImport_ImportModule("numbers");
    if (numbers_module) {
        char* register_numbers =
            "numbers.Integral.register(gmpy2.mpz)\n"
            "numbers.Rational.register(gmpy2.mpq)\n"
            "numbers.Real.register(gmpy2.mpfr)\n"
            "numbers.Complex.register(gmpy2.mpc)\n"
        ;
        namespace = PyDict_New();
        PyDict_SetItemString(namespace, "numbers", numbers_module);
        PyDict_SetItemString(namespace, "gmpy2", gmpy_module);
        result = PyRun_String(register_numbers, Py_file_input,
                              namespace, namespace);
        if (!result) {
            /* LCOV_EXCL_START */
            PyErr_Clear();
            /* LCOV_EXCL_STOP */
        }
        Py_DECREF(namespace);
        Py_DECREF(numbers_module);
        Py_XDECREF(result);
    }
    else {
        /* LCOV_EXCL_START */
        PyErr_Clear();
        /* LCOV_EXCL_STOP */
    }

    PyOS_setsig(SIGFPE, gmp_abort_handler);

    return 0;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
static PyModuleDef_Slot Pygmpy_slots[] = {
    {Py_mod_exec, gmpy_exec},
#  if PY_VERSION_HEX >= 0x030C0000
    {Py_mod_multiple_interpreters,
     Py_MOD_MULTIPLE_INTERPRETERS_NOT_SUPPORTED},
#  endif
#  if PY_VERSION_HEX >= 0x030D0000
#    ifdef MPFR_USE_THREAD_SAFE
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#    else
    {Py_mod_gil, Py_MOD_GIL_USED},
#    endif
#  endif
    {0, NULL}
};
#pragma GCC diagnostic pop

/* Notes on Python 3.x support: Full support for PEP-3121 has not been
 * implemented. No per-module state has been defined.
 */

static struct PyModuleDef gmpy_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "gmpy2",
    .m_doc = _gmpy_docs,
    .m_size = 0,
    .m_methods = Pygmpy_methods,
    .m_slots = Pygmpy_slots,
};

PyMODINIT_FUNC PyInit_gmpy2(void)
{
    return PyModuleDef_Init(&gmpy_module);
}
