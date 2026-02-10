/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_random.h                                                           *
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

#ifndef GMPY_RANDOM_H
#define GMPY_RANDOM_H

#ifdef __cplusplus
extern "C" {
#endif

/* gmpy_random C API extension header file.
 *
 * Provide support random number state.
 *
 * Version 2.00, December 2011 (created) casevh
 *
 * This file is expected to be included from gmpy.h
 */

static PyTypeObject RandomState_Type;
#define RANDOM_STATE(obj) (((RandomState_Object *)(obj))->state)
#define RandomState_Check(v) (((PyObject*)v)->ob_type == &RandomState_Type)

static RandomState_Object * GMPy_RandomState_New(void);
static void                 GMPy_RandomState_Dealloc(RandomState_Object *self);

static PyObject * GMPy_RandomState_Repr(RandomState_Object *self);
static PyObject * GMPy_RandomState_Factory(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_urandomb_Function(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_rrandomb_Function(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_random_Function(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_random_Function(PyObject *self, PyObject *args);
#if MPFR_VERSION_MAJOR > 3
static PyObject * GMPy_MPFR_nrandom_Function(PyObject *self, PyObject *args);
#endif
static PyObject * GMPy_MPFR_grandom_Function(PyObject *self, PyObject *args);
static PyObject * GMPy_MPC_random_Function(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif
