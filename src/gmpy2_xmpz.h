/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_xmpz.h                                                            *
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

#ifndef GMPY_XMPZ_H
#define GMPY_XMPZ_H

#ifdef __cplusplus
extern "C" {
#endif

static PyTypeObject XMPZ_Type;
#define XMPZ(obj) (((XMPZ_Object*)(obj))->z)
#define XMPZ_Check(v) (((PyObject*)v)->ob_type == &XMPZ_Type)
#define CHECK_MPZANY(v) (MPZ_Check(v) || XMPZ_Check(v))

typedef struct {
    PyObject_HEAD
    XMPZ_Object *bitmap;
    mp_bitcnt_t start, stop;
    int iter_type;
} GMPy_Iter_Object;

static PyTypeObject GMPy_Iter_Type;
#define GMPy_Iter_Check(v) (((PyObject*)v)->ob_type == &GMPy_Iter_Type)

#ifdef __cplusplus
}
#endif
#endif
