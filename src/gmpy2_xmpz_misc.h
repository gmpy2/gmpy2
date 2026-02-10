/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_xmpz_misc.h                                                       *
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

#ifndef GMPY_XMPZ_MISC_H
#define GMPY_XMPZ_MISC_H

#ifdef __cplusplus
extern "C" {
#endif

static int        GMPy_XMPZ_NonZero_Slot(XMPZ_Object *x);
static PyObject * GMPy_XMPZ_Function_XbitMask(PyObject *self, PyObject *other);
static PyObject * GMPy_XMPZ_Abs_Slot(XMPZ_Object *x);
static PyObject * GMPy_XMPZ_Neg_Slot(XMPZ_Object *x);
static PyObject * GMPy_XMPZ_Pos_Slot(XMPZ_Object *x);
static PyObject * GMPy_XMPZ_Com_Slot(XMPZ_Object *x);
static PyObject * GMPy_XMPZ_Method_MakeMPZ(PyObject *self, PyObject *other);
static PyObject * GMPy_XMPZ_Method_Copy(PyObject *self, PyObject *other);
static Py_ssize_t GMPy_XMPZ_Method_Length(XMPZ_Object *obj);
static PyObject * GMPy_XMPZ_Method_SubScript(XMPZ_Object* self, PyObject* item);
static int        GMPy_XMPZ_Method_AssignSubScript(XMPZ_Object* self, PyObject* item, PyObject* value);

static PyObject * GMPy_XMPZ_Attrib_GetNumer(XMPZ_Object *self, void *closure);
static PyObject * GMPy_XMPZ_Attrib_GetDenom(XMPZ_Object *self, void *closure);
static PyObject * GMPy_XMPZ_Attrib_GetReal(XMPZ_Object *self, void *closure);
static PyObject * GMPy_XMPZ_Attrib_GetImag(XMPZ_Object *self, void *closure);

static GMPy_Iter_Object * GMPy_Iter_New(void);
static void               GMPy_Iter_Dealloc(GMPy_Iter_Object *self);
static PyObject *         GMPy_Iter_Next(GMPy_Iter_Object *self);
static PyObject *         GMPy_Iter_Repr(GMPy_Iter_Object *self);
static PyObject *         GMPy_XMPZ_Method_IterBits(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *         GMPy_XMPZ_Method_IterSet(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *         GMPy_XMPZ_Method_IterClear(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *         GMPy_XMPZ_Method_SizeOf(PyObject *self, PyObject *other);


#ifdef __cplusplus
}
#endif
#endif
