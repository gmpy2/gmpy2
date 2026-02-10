/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpfr.h                                                            *
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


/* gmpy_mpfr C API extension header file.
 *
 * Provide interface to the MPFR (Multiple Precision Floating-point with
 * Rounding) library.
 *
 * Version 2.00, April 2011 (created) casevh
 */

#ifndef GMPY_MPFR_H
#define GMPY_MPFR_H

#ifdef __cplusplus
extern "C" {
#endif


#if !defined(FLT_RADIX) || (FLT_RADIX!=2)
#   error "FLT_RADIX undefined or != 2, GMPY2 is confused. :("
#endif

#if defined(MS_WIN32) && defined(_MSC_VER)
#  pragma comment(lib,"mpfr.lib")
#endif

#ifdef FAST

/* Very bad code ahead. I've copied portions of mpfr-impl.h and
 * hacked them so they work. This code will only be enabled if you
 * specify the --fast option.
 */

#define MPFR_THREAD_ATTR __thread
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR unsigned int __gmpfr_flags;
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_exp_t   __gmpfr_emin;
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_exp_t   __gmpfr_emax;

/* Flags of __gmpfr_flags */
#define MPFR_FLAGS_UNDERFLOW 1
#define MPFR_FLAGS_OVERFLOW 2
#define MPFR_FLAGS_NAN 4
#define MPFR_FLAGS_INEXACT 8
#define MPFR_FLAGS_ERANGE 16
#define MPFR_FLAGS_DIVBY0 32
#define MPFR_FLAGS_ALL (MPFR_FLAGS_UNDERFLOW | \
                        MPFR_FLAGS_OVERFLOW  | \
                        MPFR_FLAGS_NAN       | \
                        MPFR_FLAGS_INEXACT   | \
                        MPFR_FLAGS_ERANGE    | \
                        MPFR_FLAGS_DIVBY0)

/* Replace some common functions for direct access to the global vars */
#define mpfr_get_emin() (__gmpfr_emin + 0)
#define mpfr_get_emax() (__gmpfr_emax + 0)

#define mpfr_clear_flags()      ((void) (__gmpfr_flags = 0))
#define mpfr_clear_underflow()  ((void) (__gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_UNDERFLOW))
#define mpfr_clear_overflow()   ((void) (__gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_OVERFLOW))
#define mpfr_clear_nanflag()    ((void) (__gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_NAN))
#define mpfr_clear_inexflag()   ((void) (__gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_INEXACT))
#define mpfr_clear_erangeflag() ((void) (__gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_ERANGE))
#define mpfr_clear_divby0()     ((void) (__gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_DIVBY0))
#define mpfr_underflow_p()      ((int) (__gmpfr_flags & MPFR_FLAGS_UNDERFLOW))
#define mpfr_overflow_p()       ((int) (__gmpfr_flags & MPFR_FLAGS_OVERFLOW))
#define mpfr_nanflag_p()        ((int) (__gmpfr_flags & MPFR_FLAGS_NAN))
#define mpfr_inexflag_p()       ((int) (__gmpfr_flags & MPFR_FLAGS_INEXACT))
#define mpfr_erangeflag_p()     ((int) (__gmpfr_flags & MPFR_FLAGS_ERANGE))
#define mpfr_divby0_p()         ((int) (__gmpfr_flags & MPFR_FLAGS_DIVBY0))

#define mpfr_check_range(x,t,r) \
 ((((x)->_mpfr_exp) >= __gmpfr_emin && ((x)->_mpfr_exp) <= __gmpfr_emax) \
  ? ((t) ? (__gmpfr_flags |= MPFR_FLAGS_INEXACT, (t)) : 0)               \
  : mpfr_check_range(x,t,r))

/* End of the really bad code. Hopefully.
 */

#endif

static PyTypeObject MPFR_Type;
#define MPFR_Check(v) (((PyObject*)v)->ob_type == &MPFR_Type)

#define GMPY_DIVZERO(msg) PyErr_SetString(GMPyExc_DivZero, msg)
#define GMPY_INEXACT(msg) PyErr_SetString(GMPyExc_Inexact, msg)
#define GMPY_INVALID(msg) PyErr_SetString(GMPyExc_Invalid, msg)
#define GMPY_OVERFLOW(msg) PyErr_SetString(GMPyExc_Overflow, msg)
#define GMPY_UNDERFLOW(msg) PyErr_SetString(GMPyExc_Underflow, msg)
#define GMPY_ERANGE(msg) PyErr_SetString(GMPyExc_Erange, msg)

#define GMPY_MPFR_CHECK_RANGE(V, CTX) \
    if (mpfr_regular_p(V->f) && \
        (!((V->f->_mpfr_exp >= CTX->ctx.emin) && \
            (V->f->_mpfr_exp <= CTX->ctx.emax)))) { \
        mpfr_exp_t _oldemin, _oldemax; \
        _oldemin = mpfr_get_emin(); \
        _oldemax = mpfr_get_emax(); \
        mpfr_set_emin(CTX->ctx.emin); \
        mpfr_set_emax(CTX->ctx.emax); \
        V->rc = mpfr_check_range(V->f, V->rc, GET_MPFR_ROUND(CTX)); \
        mpfr_set_emin(_oldemin); \
        mpfr_set_emax(_oldemax); \
    }

#define GMPY_MPFR_SUBNORMALIZE(V, CTX) \
    if (CTX->ctx.subnormalize && \
        V->f->_mpfr_exp >= CTX->ctx.emin && \
        V->f->_mpfr_exp <= CTX->ctx.emin + mpfr_get_prec(V->f) - 2) { \
        mpfr_exp_t _oldemin, _oldemax; \
        _oldemin = mpfr_get_emin(); \
        _oldemax = mpfr_get_emax(); \
        mpfr_set_emin(CTX->ctx.emin); \
        mpfr_set_emax(CTX->ctx.emax); \
        V->rc = mpfr_subnormalize(V->f, V->rc, GET_MPFR_ROUND(CTX)); \
        mpfr_set_emin(_oldemin); \
        mpfr_set_emax(_oldemax); \
    }

/* Exceptions should be checked in order of least important to most important.
 */

#define GMPY_MPFR_EXCEPTIONS(V, CTX) \
    CTX->ctx.underflow |= mpfr_underflow_p(); \
    CTX->ctx.overflow |= mpfr_overflow_p(); \
    CTX->ctx.invalid |= mpfr_nanflag_p(); \
    CTX->ctx.inexact |= mpfr_inexflag_p(); \
    CTX->ctx.divzero |= mpfr_divby0_p(); \
    if (CTX->ctx.traps) { \
        if ((CTX->ctx.traps & TRAP_UNDERFLOW) && mpfr_underflow_p()) { \
            GMPY_UNDERFLOW("underflow"); \
            Py_XDECREF((PyObject*)V); \
            V = NULL; \
        } \
        if ((CTX->ctx.traps & TRAP_OVERFLOW) && mpfr_overflow_p()) { \
            GMPY_OVERFLOW("overflow"); \
            Py_XDECREF((PyObject*)V); \
            V = NULL; \
        } \
        if ((CTX->ctx.traps & TRAP_INEXACT) && mpfr_inexflag_p()) { \
            GMPY_INEXACT("inexact result"); \
            Py_XDECREF((PyObject*)V); \
            V = NULL; \
        } \
        if ((CTX->ctx.traps & TRAP_INVALID) && mpfr_nanflag_p()) { \
            GMPY_INVALID("invalid operation"); \
            Py_XDECREF((PyObject*)V); \
            V = NULL; \
        } \
        if ((CTX->ctx.traps & TRAP_DIVZERO) && mpfr_divby0_p()) { \
            GMPY_DIVZERO("division by zero"); \
            Py_XDECREF((PyObject*)V); \
            V = NULL; \
        } \
    }

#define GMPY_MPFR_CLEANUP(V, CTX) \
    GMPY_MPFR_CHECK_RANGE(V, CTX); \
    GMPY_MPFR_SUBNORMALIZE(V, CTX); \
    GMPY_MPFR_EXCEPTIONS(V, CTX);

#define GMPY_CHECK_ERANGE(V, CTX, MSG) \
    CTX->ctx.erange |= mpfr_erangeflag_p(); \
    if (CTX->ctx.traps) { \
        if ((CTX->ctx.traps & TRAP_ERANGE) && mpfr_erangeflag_p()) { \
            GMPY_ERANGE(MSG); \
            Py_XDECREF((PyObject*)V); \
            V = NULL; \
        } \
    } \

static void _GMPy_MPFR_Cleanup(MPFR_Object **v, CTXT_Object *ctext);
static void _GMPy_MPFR_CleanupV2(MPFR_Object **v, CTXT_Object *ctext);

#ifdef __cplusplus
}
#endif
#endif
