/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_random.c                                                           *
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

static RandomState_Object *
GMPy_RandomState_New(void)
{
    RandomState_Object *result;

    if ((result = PyObject_New(RandomState_Object, &RandomState_Type))) {
        gmp_randinit_default(result->state);
    }
    return result;
}

static void
GMPy_RandomState_Dealloc(RandomState_Object *self)
{
    gmp_randclear(self->state);
    PyObject_Free(self);
}

static PyObject *
GMPy_RandomState_Repr(RandomState_Object *self)
{
    return PyUnicode_FromString("<gmpy2.RandomState>");
}

PyDoc_STRVAR(GMPy_doc_random_state_factory,
"random_state($module, seed=0, /)\n--\n\n"
"Return new object containing state information for the random number\n"
"generator. An optional integer can be specified as the seed value.");

static PyObject *
GMPy_RandomState_Factory(PyObject *self, PyObject *args)
{
    RandomState_Object *result;
    MPZ_Object *temp;

    if (!(result = GMPy_RandomState_New())) {
        return NULL;
    }

    if (PyTuple_GET_SIZE(args) == 0) {
        gmp_randseed_ui(result->state, 0);
    }
    else if (PyTuple_GET_SIZE(args) == 1) {
        if (!(temp = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args,0), NULL))) {
            Py_DECREF((PyObject*)result);
            TYPE_ERROR("seed must be an integer");
            return NULL;
        }
        gmp_randseed(result->state, temp->z);
        Py_DECREF((PyObject*)temp);
    }
    else {
        Py_DECREF((PyObject*)result);
        TYPE_ERROR("random_state() requires 0 or 1 integer arguments");
        return NULL;
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_urandomb_function,
"mpz_urandomb($module, random_state, bit_count, /)\n--\n\n"
"Return uniformly distributed random integer between 0 and\n"
"2**bit_count-1.");

static PyObject *
GMPy_MPZ_urandomb_Function(PyObject *self, PyObject *args)
{
    MPZ_Object *result;
    PyObject *temp0, *temp1;
    unsigned long len;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("mpz_urandomb() requires 2 arguments");
        return NULL;
    }

    temp0 = PyTuple_GET_ITEM(args, 0);
    temp1 = PyTuple_GET_ITEM(args, 1);

    if (!RandomState_Check(temp0)) {
        TYPE_ERROR("mpz_urandomb() requires 'random_state' and 'bit_count' arguments");
        return NULL;
    }

    len = GMPy_Integer_AsUnsignedLongWithType(temp1, GMPy_ObjectType(temp1));
    if (len == (unsigned long)(-1) && PyErr_Occurred()) {
        TYPE_ERROR("mpz_urandomb() requires 'random_state' and 'bit_count' arguments");
        return NULL;
    }

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_urandomb(result->z, RANDOM_STATE(temp0), len);
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_rrandomb_function,
"mpz_rrandomb($module, random_state, bit_count, /)\n--\n\n"
"Return a random integer between 0 and 2**bit_count-1 with long\n"
"sequences of zeros and one in its binary representation.");

static PyObject *
GMPy_MPZ_rrandomb_Function(PyObject *self, PyObject *args)
{
    MPZ_Object *result;
    PyObject *temp0, *temp1;
    unsigned long len;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("mpz_rrandomb() requires 2 arguments");
        return NULL;
    }

    temp0 = PyTuple_GET_ITEM(args, 0);
    temp1 = PyTuple_GET_ITEM(args, 1);

    if (!RandomState_Check(temp0)) {
        TYPE_ERROR("mpz_rrandomb() requires 'random_state' and 'bit_count' arguments");
        return NULL;
    }

    len = GMPy_Integer_AsUnsignedLongWithType(temp1, GMPy_ObjectType(temp1));
    if (len == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        TYPE_ERROR("mpz_rrandomb() requires 'random_state' and 'bit_count' arguments");
        return NULL;
    }

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_rrandomb(result->z, RANDOM_STATE(temp0), len);
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_random_function,
"mpz_random($module, random_state, int, /)\n--\n\n"
"Return uniformly distributed random integer between 0 and n-1.");

static PyObject *
GMPy_MPZ_random_Function(PyObject *self, PyObject *args)
{
    MPZ_Object *result, *temp;
    PyObject *temp0, *temp1;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("mpz_random() requires 2 arguments");
        return NULL;
    }

    temp0 = PyTuple_GET_ITEM(args, 0);
    temp1 = PyTuple_GET_ITEM(args, 1);

    if (!RandomState_Check(temp0)) {
        TYPE_ERROR("mpz_random() requires 'random_state' and 'int' arguments");
        return NULL;
    }

    if (!(temp = GMPy_MPZ_From_IntegerWithType(temp1, GMPy_ObjectType(temp1), NULL))) {
        TYPE_ERROR("mpz_random() requires 'random_state' and 'int' arguments");
        return NULL;
    }

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_urandomm(result->z, RANDOM_STATE(PyTuple_GET_ITEM(args, 0)), temp->z);
    }

    Py_DECREF((PyObject*)temp);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpfr_random_function,
"mpfr_random($module, random_state, /)\n--\n\n"
"Return uniformly distributed number between [0,1].");

static PyObject *
GMPy_MPFR_random_Function(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("mpfr_random() requires 1 argument");
        return NULL;
    }

    if (!RandomState_Check(PyTuple_GET_ITEM(args, 0))) {
        TYPE_ERROR("mpfr_random() requires 'random_state' argument");
        return NULL;
    }

    if ((result = GMPy_MPFR_New(0, context))) {
        mpfr_urandom(result->f, RANDOM_STATE(PyTuple_GET_ITEM(args, 0)), GET_MPFR_ROUND(context));
    }

    return (PyObject*)result;
}

#if MPFR_VERSION_MAJOR > 3

PyDoc_STRVAR(GMPy_doc_mpfr_nrandom_function,
"mpfr_nrandom($module, random_state, /)\n--\n\n"
"Return a random number with gaussian distribution.");

static PyObject *
GMPy_MPFR_nrandom_Function(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("mpfr_nrandom() requires 1 argument");
        return NULL;
    }

    if (!RandomState_Check(PyTuple_GET_ITEM(args, 0))) {
        TYPE_ERROR("mpfr_nrandom() requires 'random_state' argument");
        return NULL;
    }

    if ((result = GMPy_MPFR_New(0, context))) {
        mpfr_nrandom(result->f,
                    RANDOM_STATE(PyTuple_GET_ITEM(args, 0)),
                    GET_MPFR_ROUND(context));
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpfr_grandom_function,
"mpfr_grandom($module, random_state, /)\n--\n\n"
"Return two random numbers with gaussian distribution.");

static PyObject *
GMPy_MPFR_grandom_Function(PyObject *self, PyObject *args)
{
    MPFR_Object *result1, *result2;
    PyObject *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("mpfr_grandom() requires 1 argument");
        return NULL;
    }

    if (!RandomState_Check(PyTuple_GET_ITEM(args, 0))) {
        TYPE_ERROR("mpfr_grandom() requires 'random_state' argument");
        return NULL;
    }

    result1 = GMPy_MPFR_New(0, context);
    result2 = GMPy_MPFR_New(0, context);
    if (!result1 || !result2) {
        Py_XDECREF((PyObject*)result1);
        Py_XDECREF((PyObject*)result2);
        return NULL;
    }

    mpfr_nrandom(result1->f,
                 RANDOM_STATE(PyTuple_GET_ITEM(args, 0)),
                 GET_MPFR_ROUND(context));

    mpfr_nrandom(result2->f,
                 RANDOM_STATE(PyTuple_GET_ITEM(args, 0)),
                 GET_MPFR_ROUND(context));

    result = Py_BuildValue("(NN)", (PyObject*)result1, (PyObject*)result2);
    if (!result) {
        Py_DECREF((PyObject*)result1);
        Py_DECREF((PyObject*)result2);
    }
    return result;
}
#else
PyDoc_STRVAR(GMPy_doc_mpfr_grandom_function,
"mpfr_grandom($module, random_state, /)\n--\n\n"
"Return two random numbers with gaussian distribution.");

static PyObject *
GMPy_MPFR_grandom_Function(PyObject *self, PyObject *args)
{
    MPFR_Object *result1, *result2;
    PyObject *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("mpfr_grandom() requires 1 argument");
        return NULL;
    }

    if (!RandomState_Check(PyTuple_GET_ITEM(args, 0))) {
        TYPE_ERROR("mpfr_grandom() requires 'random_state' argument");
        return NULL;
    }

    result1 = GMPy_MPFR_New(0, context);
    result2 = GMPy_MPFR_New(0, context);
    if (!result1 || !result2) {
        Py_XDECREF((PyObject*)result1);
        Py_XDECREF((PyObject*)result2);
        return NULL;
    }

    mpfr_grandom(result1->f, result2->f,
                 RANDOM_STATE(PyTuple_GET_ITEM(args, 0)),
                 GET_MPFR_ROUND(context));

    result = Py_BuildValue("(NN)", (PyObject*)result1, (PyObject*)result2);
    if (!result) {
        Py_DECREF((PyObject*)result1);
        Py_DECREF((PyObject*)result2);
    }
    return result;
}
#endif


PyDoc_STRVAR(GMPy_doc_mpc_random_function,
"mpc_random($module, random_state, /)\n--\n\n"
"Return uniformly distributed number in the unit square [0,1]x[0,1].");

static PyObject *
GMPy_MPC_random_Function(PyObject *self, PyObject *args)
{
    MPC_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("mpfc_random() requires 1 argument");
        return NULL;
    }

    if (!RandomState_Check(PyTuple_GET_ITEM(args, 0))) {
        TYPE_ERROR("mpc_random() requires 'random_state' argument");
        return NULL;
    }

    if ((result = GMPy_MPC_New(0, 0, context))) {
        mpc_urandom(result->c, RANDOM_STATE(PyTuple_GET_ITEM(args, 0)));
    }

    return (PyObject*)result;
}

static PyTypeObject RandomState_Type =
{
    PyVarObject_HEAD_INIT(0, 0)
    .tp_name = "gmpy2.random_state",       
    .tp_basicsize = sizeof(RandomState_Object),     
    .tp_dealloc = (destructor) GMPy_RandomState_Dealloc, 
    .tp_repr = (reprfunc) GMPy_RandomState_Repr,
    .tp_flags = Py_TPFLAGS_DEFAULT, 
    .tp_doc = "GMPY2 Random number generator state", 
};
