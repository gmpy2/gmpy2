/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_misc.c                                                            *
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

/* Miscellaneous module-level functions and helper functions. */

PyDoc_STRVAR(GMPy_doc_license,
"license() -> str\n\n"
"Return string giving license information.");

static PyObject *
GMPy_get_license(PyObject *self, PyObject *args)
{
    return PyUnicode_FromString(gmpy_license);
}

PyDoc_STRVAR(GMPy_doc_version,
"version() -> str\n\n"
"Return string giving current GMPY2 version.");

static PyObject *
GMPy_get_version(PyObject *self, PyObject *args)
{
    PyObject *gmpy2 = PyImport_ImportModule("gmpy2");

    if (!gmpy2) {
        return NULL; /* LCOV_EXCL_LINE */
    }

    PyObject *version = PyObject_GetAttrString(gmpy2, "__version__");

    Py_DECREF(gmpy2);
    return version;
}

PyDoc_STRVAR(GMPy_doc_mp_version,
"mp_version() -> str\n\n"
"Return string giving current GMP version.");

static PyObject *
GMPy_get_mp_version(PyObject *self, PyObject *args)
{
    return PyUnicode_FromFormat("GMP %s", gmp_version);
}

PyDoc_STRVAR(GMPy_doc_mpfr_version,
"mpfr_version() -> str\n\n"
"Return string giving current MPFR version.");

static PyObject *
GMPy_get_mpfr_version(PyObject *self, PyObject *args)
{
    return PyUnicode_FromFormat("MPFR %s", MPFR_VERSION_STRING);
}

PyDoc_STRVAR(GMPy_doc_mpc_version,
"mpc_version() -> str\n\n"
"Return string giving current MPC version.");

static PyObject *
GMPy_get_mpc_version(PyObject *self, PyObject *args)
{
    return PyUnicode_FromFormat("MPC %s", MPC_VERSION_STRING);
}

PyDoc_STRVAR(GMPy_doc_mp_limbsize,
"mp_limbsize() -> int\n\n\
Return the number of bits per limb.");

static PyObject *
GMPy_get_mp_limbsize(PyObject *self, PyObject *args)
{
    return PyLong_FromLong(mp_bits_per_limb);
}
