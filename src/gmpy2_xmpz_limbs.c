/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_xmpz_limbs.c                                                      *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP, MPFR, and MPC multiple precision           *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2020 Tyler Lanphear                                           *
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

PyDoc_STRVAR(GMPy_doc_xmpz_method_num_limbs,
"num_limbs($self)\n--\n\n"
"Return the number of limbs of self.");
static PyObject* GMPy_XMPZ_Method_NumLimbs(PyObject* obj, PyObject* other)
{
    return PyLong_FromSize_t(mpz_size(XMPZ(obj)));
}

PyDoc_STRVAR(GMPy_doc_xmpz_method_limbs_read,
"limbs_read($self)\n--\n\n"
"Returns the address of the immutable buffer representing the\n"
"limbs of self.");
static PyObject* GMPy_XMPZ_Method_LimbsRead(PyObject* obj, PyObject* args)
{
    const mp_limb_t* limbs = mpz_limbs_read(XMPZ(obj));
    return PyLong_FromVoidPtr((void *) limbs);
}

PyDoc_STRVAR(GMPy_doc_xmpz_method_limbs_write,
"limbs_write($self, n, /)\n--\n\n"
"Returns the address of a mutable buffer representing the limbs\n"
"of self, resized so that it may hold at least n limbs.\n"
"Must be followed by a call to self.limbs_finish(n) after writing to\n"
"the returned address in order for the changes to take effect.\n"
"WARNING: this operation is destructive and may destroy the old\n"
"value of self.");
static PyObject* GMPy_XMPZ_Method_LimbsWrite(PyObject* obj, PyObject* other)
{

    if (!PyLong_Check(other)) {
        TYPE_ERROR("number of limbs must be an int or a long");
        return NULL;
    }
    else {
        size_t num_limbs = (size_t) PyLong_AsSsize_t(other);
        mp_limb_t * limbs = mpz_limbs_write(XMPZ(obj), (mp_size_t) num_limbs);
        return PyLong_FromVoidPtr((void *) limbs);
    }
}

PyDoc_STRVAR(GMPy_doc_xmpz_method_limbs_modify,
"limbs_modify($self, n, /)\n--\n\n"
"Returns the address of a mutable buffer representing the limbs\n"
"of self, resized so that it may hold at least n limbs.\n"
"Must be followed by a call to self.limbs_finish(n) after writing to\n"
"the returned address in order for the changes to take effect.");
static PyObject* GMPy_XMPZ_Method_LimbsModify(PyObject* obj, PyObject* other)
{
    if (!PyLong_Check(other)) {
        TYPE_ERROR("number of limbs must be an int or a long");
        return NULL;
    }
    else {
        size_t num_limbs = (size_t) PyLong_AsSsize_t(other);
        mp_limb_t * limbs = mpz_limbs_modify(XMPZ(obj), (mp_size_t) num_limbs);
        return PyLong_FromVoidPtr((void *) limbs);
    }
}

PyDoc_STRVAR(GMPy_doc_xmpz_method_limbs_finish,
"limbs_finish($self, n, /)\n--\n\n"
"Must be called after writing to the address returned by\n"
"self.limbs_write(n) or self.limbs_modify(n) to update\n"
"the limbs of self.");
static PyObject* GMPy_XMPZ_Method_LimbsFinish(PyObject* obj, PyObject* other)
{
    if (!PyLong_Check(other)) {
        TYPE_ERROR("number of limbs must be an int or long");
        return NULL;
    }
    else {
        size_t num_limbs = (size_t) PyLong_AsSsize_t(other);
        mpz_limbs_finish(XMPZ(obj), num_limbs);
        Py_RETURN_NONE;
    }
}
