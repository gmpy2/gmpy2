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

static PyTypeObject MPFR_Type;
#define MPFR_Check(v) (Py_TYPE(v) == &MPFR_Type)

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

#ifdef __cplusplus
}
#endif
#endif
