/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_cache.h                                                           *
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


/* gmpy2 caches objects so they can be reused quickly without involving a new
 * memory allocation or object construction.
 */

#ifndef GMPY_CACHE_H
#define GMPY_CACHE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Private functions */

/* C-API functions */

/* static MPZ_Object *  GMPy_MPZ_New(CTXT_Object *context); */
/* static PyObject *    GMPy_MPZ_NewInit(PyTypeObject *type, PyObject *args, PyObject *keywds); */
/* static void          GMPy_MPZ_Dealloc(MPZ_Object *self); */
static GMPy_MPZ_New_RETURN     GMPy_MPZ_New     GMPy_MPZ_New_PROTO;
static GMPy_MPZ_NewInit_RETURN GMPy_MPZ_NewInit GMPy_MPZ_NewInit_PROTO;
static GMPy_MPZ_Dealloc_RETURN GMPy_MPZ_Dealloc GMPy_MPZ_Dealloc_PROTO;

/* static XMPZ_Object *  GMPy_XMPZ_New(CTXT_Object *context); */
/* static PyObject *     GMPy_XMPZ_NewInit(PyTypeObject *type, PyObject *args, PyObject *keywds); */
/* static void           GMPy_XMPZ_Dealloc(XMPZ_Object *self); */
static GMPy_XMPZ_New_RETURN     GMPy_XMPZ_New     GMPy_XMPZ_New_PROTO;
static GMPy_XMPZ_NewInit_RETURN GMPy_XMPZ_NewInit GMPy_XMPZ_NewInit_PROTO;
static GMPy_XMPZ_Dealloc_RETURN GMPy_XMPZ_Dealloc GMPy_XMPZ_Dealloc_PROTO;

/* static MPQ_Object *  GMPy_MPQ_New(CTXT_Object *context); */
/* static PyObject *    GMPy_MPQ_NewInit(PyTypeObject *type, PyObject *args, PyObject *keywds); */
/* static void          GMPy_MPQ_Dealloc(MPQ_Object *self); */
static GMPy_MPQ_New_RETURN     GMPy_MPQ_New     GMPy_MPQ_New_PROTO;
static GMPy_MPQ_NewInit_RETURN GMPy_MPQ_NewInit GMPy_MPQ_NewInit_PROTO;
static GMPy_MPQ_Dealloc_RETURN GMPy_MPQ_Dealloc GMPy_MPQ_Dealloc_PROTO;

/* static MPFR_Object * GMPy_MPFR_New(CTXT_Object *context); */
/* static PyObject *    GMPy_MPFR_NewInit(PyTypeObject *type, PyObject *args, PyObject *keywds); */
/* static void          GMPy_MPFR_Dealloc(MPFR_Object *self); */
static GMPy_MPFR_New_RETURN     GMPy_MPFR_New     GMPy_MPFR_New_PROTO;
static GMPy_MPFR_NewInit_RETURN GMPy_MPFR_NewInit GMPy_MPFR_NewInit_PROTO;
static GMPy_MPFR_Dealloc_RETURN GMPy_MPFR_Dealloc GMPy_MPFR_Dealloc_PROTO;

static MPC_Object *  GMPy_MPC_New(mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static void          GMPy_MPC_Dealloc(MPC_Object *self);

#ifdef __cplusplus
}
#endif
#endif
