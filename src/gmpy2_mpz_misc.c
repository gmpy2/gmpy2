/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpz_misc.c                                                        *
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

/* return number-of-digits for an mpz in requested base, default 10 */
PyDoc_STRVAR(GMPy_doc_mpz_method_num_digits,
"num_digits($self, base=10, /)\n--\n\n"
"Return length of string representing the absolute value of self in\n"
"the given base. Values  for base can range between 2 and 62. The\n"
"value returned may be 1 too large.");

PyDoc_STRVAR(GMPy_doc_mpz_function_num_digits,
"num_digits($module, x, base=10, /)\n--\n\n"
"Return length of string representing the absolute value of x in\n"
"the given base. Values  for base can range between 2 and 62. The\n"
"value returned may be 1 too large.");

static PyObject *
GMPy_MPZ_Method_NumDigits(PyObject *self, PyObject *const *args,
                          Py_ssize_t nargs)
{
    long base = 10;
    PyObject *result;

    if (nargs == 1) {
        base = PyLong_AsLong(args[0]);
        if (base == -1 && PyErr_Occurred()) {
            return NULL;
        }
    }

    if ((base < 2) || (base > 62)) {
        VALUE_ERROR("base must be in the interval [2, 62]");
        return NULL;
    }

    result = PyLong_FromSize_t(mpz_sizeinbase(MPZ(self), (int)base));
    return result;
}

static PyObject *
GMPy_MPZ_Function_NumDigits(PyObject *self, PyObject *const *args,
                            Py_ssize_t nargs)
{
    long base = 10;
    MPZ_Object *temp;
    PyObject *result;

    if (nargs == 0 || nargs > 2) {
        TYPE_ERROR("num_digits() requires 'mpz',['int'] arguments");
        return NULL;
    }

    if (nargs == 2) {
        base = PyLong_AsLong(args[1]);
        if (base == -1 && PyErr_Occurred()) {
            return NULL;
        }
    }

    if ((base < 2) || (base > 62)) {
        VALUE_ERROR("base must be in the interval [2, 62]");
        return NULL;
    }

    if (!(temp = GMPy_MPZ_From_Integer(args[0], NULL))) {
        return NULL;
    }

    result = PyLong_FromSize_t(mpz_sizeinbase(temp->z, (int)base));
    Py_DECREF((PyObject*)temp);
    return result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_iroot,
"iroot($module, x, n, /)\n--\n\n"
"Return the integer n-th root of x and boolean value that is `True`\n"
"iff the root is exact. x >= 0. n > 0.");

static PyObject *
GMPy_MPZ_Function_Iroot(PyObject *self, PyObject *const *args,
                        Py_ssize_t nargs)
{
    unsigned long n;
    int exact, is_signed = 0;
    MPZ_Object *root = NULL, *tempx = NULL;
    PyObject *result = NULL;

    if (nargs != 2 || !IS_INTEGER(args[0]) || !IS_INTEGER(args[1])) {
        TYPE_ERROR("iroot() requires 'int','int' arguments");
        return NULL;
    }

    n = (unsigned long)GMPy_Integer_AsUnsignedLongOrLong(args[1], &is_signed);
    if ((n == (unsigned long)(-1)) && PyErr_Occurred()) {
        return NULL;
    }
    if (is_signed || !n) {
        VALUE_ERROR("n must be > 0");
        return NULL;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(args[0], NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (mpz_sgn(tempx->z) < 0) {
        VALUE_ERROR("iroot() of negative number");
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    if (!(result = PyTuple_New(2)) ||
        !(root = GMPy_MPZ_New(NULL))) {

        /* LCOV_EXCL_START */
        Py_DECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)root);
        Py_XDECREF(result);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    exact = mpz_root(root->z, tempx->z, n);
    Py_DECREF((PyObject*)tempx);

    PyTuple_SET_ITEM(result, 0, (PyObject*)root);
    PyTuple_SET_ITEM(result, 1, (PyObject*)PyBool_FromLong(exact));
    return result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_iroot_rem,
"iroot_rem($module, x, n, /)\n--\n\n"
"Return a 2-element tuple (y,r), such that y is the integer n-th\n"
"root of x and x=y**n + r. x >= 0. n > 0.");

static PyObject *
GMPy_MPZ_Function_IrootRem(PyObject *self, PyObject *const *args,
                           Py_ssize_t nargs)
{
    unsigned long n;
    MPZ_Object *root = NULL, *rem = NULL, *tempx = NULL;
    PyObject *result = NULL;

    if (nargs != 2 || !IS_INTEGER(args[0]) || !IS_INTEGER(args[1])) {
        TYPE_ERROR("iroot_rem() requires 'int','int' arguments");
        return NULL;
    }

    n = GMPy_Integer_AsUnsignedLong(args[1]);
    if ((n == 0) || ((n == (unsigned long)(-1)) && PyErr_Occurred())) {
        VALUE_ERROR("n must be > 0");
        return NULL;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(args[0], NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (mpz_sgn(tempx->z) < 0) {
        VALUE_ERROR("iroot_rem() of negative number");
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    if (!(result = PyTuple_New(2)) ||
        !(root = GMPy_MPZ_New(NULL)) ||
        !(rem = GMPy_MPZ_New(NULL))) {

        /* LCOV_EXCL_START */
        Py_DECREF((PyObject*)tempx);
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)root);
        Py_XDECREF((PyObject*)rem);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_rootrem(root->z, rem->z, tempx->z, n);
    Py_DECREF((PyObject*)tempx);

    PyTuple_SET_ITEM(result, 0, (PyObject*)root);
    PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
    return result;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_ceil, "Ceiling of an mpz returns itself.");

static PyObject *
GMPy_MPZ_Method_Ceil(PyObject *self, PyObject *other)
{
    Py_INCREF(self);
    return self;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_floor, "Floor of an mpz returns itself.");

static PyObject *
GMPy_MPZ_Method_Floor(PyObject *self, PyObject *other)
{
    Py_INCREF(self);
    return self;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_trunc, "Truncating an mpz returns itself.");

static PyObject *
GMPy_MPZ_Method_Trunc(PyObject *self, PyObject *other)
{
    Py_INCREF(self);
    return self;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_round,
"__round__($self, ndigits=0, /)\n--\n\nRound an mpz to power of 10.");

static PyObject *
GMPy_MPZ_Method_Round(PyObject *self, PyObject *const *args,
                      Py_ssize_t nargs)
{
    Py_ssize_t round_digits;
    MPZ_Object *result;

    if (nargs == 0) {
        Py_INCREF(self);
        return self;
    }

    round_digits = GMPy_Integer_AsSsize_t(args[0]);
    if (round_digits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("__round__() requires 'int' argument");
        return NULL;
    }

    if (round_digits >= 0) {
        Py_INCREF(self);
        return self;
    }

    /* We can now assume round_digits > 0. */
    round_digits = -round_digits;

    if ((result = GMPy_MPZ_New(NULL))) {
        if ((unsigned)round_digits > mpz_sizeinbase(MPZ(self), 10)) {
            mpz_set_ui(result->z, 0);
        }
        else {
            mpz_t temp;
            CTXT_Object *context = NULL;
            mpfr_rnd_t ctx_round;

            CHECK_CONTEXT(context);
            ctx_round = GET_MPFR_ROUND(context);
            mpz_init(temp);
            mpz_ui_pow_ui(temp, 10, (unsigned long)round_digits);
            if (ctx_round == MPFR_RNDA) {
                if (mpz_sgn(MPZ(self)) < 0) {
                    ctx_round = MPFR_RNDD;
                }
                else {
                    ctx_round = MPFR_RNDU;
                }
            }
            switch (ctx_round) {
                case MPFR_RNDZ:
                    mpz_tdiv_q(result->z, MPZ(self), temp);
                    break;
                case MPFR_RNDD:
                    mpz_fdiv_q(result->z, MPZ(self), temp);
                    break;
                case MPFR_RNDU:
                    mpz_cdiv_q(result->z, MPZ(self), temp);
                    break;
                case MPFR_RNDN:
                default:
                    {
                        mpz_t rem;

                        mpz_init(rem);
                        mpz_fdiv_qr(result->z, rem, MPZ(self), temp);
                        mpz_mul_2exp(rem, rem, 1);
                        if (mpz_cmp(rem, temp) > 0) {
                            mpz_add_ui(result->z, result->z, 1);
                        }
                        else if (mpz_cmp(rem, temp) == 0) {
                            if (mpz_odd_p(result->z)) {
                                mpz_add_ui(result->z, result->z, 1);
                            }
                        }
                        mpz_clear(rem);
                    }
            }
            mpz_mul(result->z, result->z, temp);
            mpz_clear(temp);
        }
    }

    return (PyObject*)result;
}

static int
GMPy_MPZ_NonZero_Slot(MPZ_Object *self)
{
    return mpz_sgn(self->z) != 0;
}

/* Miscellaneous gmpy functions */

PyDoc_STRVAR(GMPy_doc_mpz_function_gcd,
"gcd($module, /, *integers)\n--\n\n"
"Return the greatest common divisor of integers.");

static PyObject *
GMPy_MPZ_Function_GCD(PyObject *self, PyObject * const *args, Py_ssize_t nargs)
{
    MPZ_Object *arg, *result = NULL;
    CTXT_Object *context = NULL;
    Py_ssize_t i;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    for (i = 0; i < nargs; i++) {
        if (!(arg = GMPy_MPZ_From_Integer(args[i], context))) {
            TYPE_ERROR("gcd() requires 'mpz' arguments");
            Py_XDECREF((PyObject*)arg);
            Py_DECREF((PyObject*)result);
            return NULL;
        }

        if (mpz_cmp_si(MPZ(result), 1) != 0) {
            GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
            mpz_gcd(MPZ(result), MPZ(arg), MPZ(result));
            GMPY_MAYBE_END_ALLOW_THREADS(context);
        }

        Py_DECREF((PyObject*)arg);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_lcm,
"lcm($module, /, *integers)\n--\n\n"
"Return the lowest common multiple of integers.");

static PyObject *
GMPy_MPZ_Function_LCM(PyObject *self, PyObject * const *args, Py_ssize_t nargs)
{
    MPZ_Object *arg, *result = NULL;
    CTXT_Object *context = NULL;
    Py_ssize_t i;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_set_ui(result->z, 1);

    for (i = 0; i < nargs; i++) {
        if (!(arg = GMPy_MPZ_From_Integer(args[i], context))) {
            TYPE_ERROR("lcm() requires 'mpz' arguments");
            Py_XDECREF((PyObject*)arg);
            Py_DECREF((PyObject*)result);
            return NULL;
        }

        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_lcm(MPZ(result), MPZ(arg), MPZ(result));
        GMPY_MAYBE_END_ALLOW_THREADS(context);

        Py_DECREF((PyObject*)arg);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_gcdext,
"gcdext($module, a, b, /)\n--\n\n"
"Return a 3-element tuple (g,s,t) such that g == gcd(a,b)\n"
"and g == a*s + b*t.");

static PyObject *
GMPy_MPZ_Function_GCDext(PyObject *self, PyObject * const *args,
                         Py_ssize_t nargs)
{
    PyObject *arg0, *arg1, *result = NULL;
    MPZ_Object *g = NULL, *s = NULL, *t = NULL, *tempa = NULL, *tempb = NULL;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if(nargs != 2) {
        TYPE_ERROR("gcdext() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = PyTuple_New(3)) ||
        !(g = GMPy_MPZ_New(NULL)) ||
        !(s = GMPy_MPZ_New(NULL)) ||
        !(t = GMPy_MPZ_New(NULL))) {

        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)g);
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)t);
        Py_XDECREF(result);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    arg0 = args[0];
    arg1 = args[1];

    if (MPZ_Check(arg0) && MPZ_Check(arg1)) {
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_gcdext(g->z, s->z, t->z, MPZ(arg0), MPZ(arg1));
        GMPY_MAYBE_END_ALLOW_THREADS(context);
    }
    else {
        if(!(tempa = GMPy_MPZ_From_Integer(arg0, NULL)) ||
           !(tempb = GMPy_MPZ_From_Integer(arg1, NULL))) {

            TYPE_ERROR("gcdext() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempa);
            Py_XDECREF((PyObject*)tempb);
            Py_DECREF((PyObject*)g);
            Py_DECREF((PyObject*)s);
            Py_DECREF((PyObject*)t);
            Py_DECREF(result);
            return NULL;
        }
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_gcdext(g->z, s->z, t->z, tempa->z, tempb->z);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_DECREF((PyObject*)tempa);
        Py_DECREF((PyObject*)tempb);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)g);
    PyTuple_SET_ITEM(result, 1, (PyObject*)s);
    PyTuple_SET_ITEM(result, 2, (PyObject*)t);
    return result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_divm,
"divm($module, a, b, m, /)\n--\n\n"
"Return x such that b*x == a mod m. Raises a `ZeroDivisionError`\n"
"exception if no such value x exists.");

static PyObject *
GMPy_MPZ_Function_Divm(PyObject *self, PyObject * const *args,
                       Py_ssize_t nargs)
{
    MPZ_Object *result = NULL, *num = NULL, *den = NULL, *mod = NULL;
    mpz_t numz, denz, modz, gcdz;
    int ok = 0;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (nargs != 3) {
        TYPE_ERROR("divm() requires 'mpz','mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (!(num = GMPy_MPZ_From_Integer(args[0], NULL)) ||
        !(den = GMPy_MPZ_From_Integer(args[1], NULL)) ||
        !(mod = GMPy_MPZ_From_Integer(args[2], NULL))) {

        TYPE_ERROR("divm() requires 'mpz','mpz','mpz' arguments");
        Py_XDECREF((PyObject*)num);
        Py_XDECREF((PyObject*)den);
        Py_XDECREF((PyObject*)mod);
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    /* Make copies so we don't destroy the input. */
    mpz_init(numz);
    mpz_init(denz);
    mpz_init(modz);
    mpz_set(numz, num->z);
    mpz_set(denz, den->z);
    mpz_set(modz, mod->z);
    Py_DECREF((PyObject*)num);
    Py_DECREF((PyObject*)den);
    Py_DECREF((PyObject*)mod);


    GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
    ok = mpz_invert(result->z, denz, modz);
    GMPY_MAYBE_END_ALLOW_THREADS(context);

    if (!ok) {
        /* last-ditch attempt: do num, den AND mod have a gcd>1 ? */
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_init(gcdz);
        mpz_gcd(gcdz, numz, denz);
        mpz_gcd(gcdz, gcdz, modz);
        mpz_divexact(numz, numz, gcdz);
        mpz_divexact(denz, denz, gcdz);
        mpz_divexact(modz, modz, gcdz);
        mpz_clear(gcdz);
        ok = mpz_invert(result->z, denz, modz);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
    }

    if (ok) {
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_mul(result->z, result->z, numz);
        mpz_mod(result->z, result->z, modz);
        mpz_clear(numz);
        mpz_clear(denz);
        mpz_clear(modz);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        return (PyObject*)result;
    }
    else {
        ZERO_ERROR("not invertible");
        mpz_clear(numz);
        mpz_clear(denz);
        mpz_clear(modz);
        Py_DECREF((PyObject*)result);
        return NULL;
    }
}

PyDoc_STRVAR(GMPy_doc_mpz_function_fac,
"fac($module, n, /)\n--\n\n"
"Return the exact factorial of n.\n\n"
"See factorial(n) to get the floating-point approximation.");

static PyObject *
GMPy_MPZ_Function_Fac(PyObject *self, PyObject *other)
{
    MPZ_Object *result = NULL;
    unsigned long n;

    n = GMPy_Integer_AsUnsignedLong(other);
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_fac_ui(result->z, n);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_double_fac,
"double_fac($module, n, /)\n--\n\n"
"Return the exact double factorial (n!!) of n. The double\n"
"factorial is defined as n*(n-2)*(n-4)...");

static PyObject *
GMPy_MPZ_Function_DoubleFac(PyObject *self, PyObject *other)
{
    MPZ_Object *result = NULL;
    unsigned long n;

    n = GMPy_Integer_AsUnsignedLong(other);
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_2fac_ui(result->z, n);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_primorial,
"primorial($module, n, /)\n--\n\n"
"Return the product of all positive prime numbers less than or\n"
"equal to n.");

static PyObject *
GMPy_MPZ_Function_Primorial(PyObject *self, PyObject *other)
{
    MPZ_Object *result = NULL;
    unsigned long n;

    n = GMPy_Integer_AsUnsignedLong(other);
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_primorial_ui(result->z, n);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_multi_fac,
"multi_fac($module, n, m, /)\n--\n\n"
"Return the exact m-multi factorial of n. The m-multi"
"factorial is defined as n*(n-m)*(n-2m)...");

static PyObject *
GMPy_MPZ_Function_MultiFac(PyObject *self, PyObject *const *args,
                           Py_ssize_t nargs)
{
    MPZ_Object *result = NULL;
    unsigned long n, m;

    if (nargs != 2) {
        TYPE_ERROR("multi_fac() requires 2 integer arguments");
        return NULL;
    }

    n = GMPy_Integer_AsUnsignedLong(args[0]);
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    m = GMPy_Integer_AsUnsignedLong(args[1]);
    if (m == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_mfac_uiui(result->z, n, m);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_fib,
"fib($module, n, /)\n--\n\n"
"Return the n-th Fibonacci number.");

static PyObject *
GMPy_MPZ_Function_Fib(PyObject *self, PyObject *other)
{
    MPZ_Object *result = NULL;
    unsigned long  n;

    n = GMPy_Integer_AsUnsignedLong(other);
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }
    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_fib_ui(result->z, n);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_fib2,
"fib2($module, n, /)\n--\n\n"
"Return a 2-tuple with the (n-1)-th and n-th Fibonacci numbers.");

static PyObject *
GMPy_MPZ_Function_Fib2(PyObject *self, PyObject *other)
{
    PyObject *result = NULL;
    MPZ_Object *fib1 = NULL, *fib2 = NULL;
    unsigned long n;

    n = GMPy_Integer_AsUnsignedLong(other);
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if (!(result = PyTuple_New(2)) ||
        !(fib1 = GMPy_MPZ_New(NULL)) ||
        !(fib2 = GMPy_MPZ_New(NULL))) {

        /* LCOV_EXCL_START */
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)fib1);
        Py_XDECREF((PyObject*)fib2);
        result = NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_fib2_ui(fib1->z, fib2->z, n);
    PyTuple_SET_ITEM(result, 0, (PyObject*)fib1);
    PyTuple_SET_ITEM(result, 1, (PyObject*)fib2);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_lucas,
"lucas($module, n, /)\n--\n\n"
"Return the n-th Lucas number.");

static PyObject *
GMPy_MPZ_Function_Lucas(PyObject *self, PyObject *other)
{
    MPZ_Object *result = NULL;
    unsigned long n;

    n = GMPy_Integer_AsUnsignedLong(other);
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_lucnum_ui(result->z, n);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_lucas2,
"lucas2($module, n, /)\n--\n\n"
"Return a 2-tuple with the (n-1)-th and n-th Lucas numbers.");

static PyObject *
GMPy_MPZ_Function_Lucas2(PyObject *self, PyObject *other)
{
    PyObject *result = NULL;
    MPZ_Object *luc1 = NULL, *luc2 = NULL;
    unsigned long n;

    n = GMPy_Integer_AsUnsignedLong(other);
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if (!(result = PyTuple_New(2)) ||
        !(luc1 = GMPy_MPZ_New(NULL)) ||
        !(luc2 = GMPy_MPZ_New(NULL))) {

        /* LCOV_EXCL_START */
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)luc1);
        Py_XDECREF((PyObject*)luc2);
        result = NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_lucnum2_ui(luc1->z, luc2->z, n);
    PyTuple_SET_ITEM(result, 0, (PyObject*)luc1);
    PyTuple_SET_ITEM(result, 1, (PyObject*)luc2);
    return result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_bincoef,
"bincoef($module, n, k, /)\n--\n\n"
"Return the binomial coefficient ('n choose k'). k >= 0.");

PyDoc_STRVAR(GMPy_doc_mpz_function_comb,
"comb($module, n, k, /)\n--\n\n"
"Return the number of combinations of 'n things, taking k at a\n"
"time'. k >= 0. Same as bincoef(n, k)");

static PyObject *
GMPy_MPZ_Function_Bincoef(PyObject *self, PyObject * const *args,
                          Py_ssize_t nargs)
{
    MPZ_Object *result = NULL, *tempx;
    unsigned long n, k;

    if (nargs != 2) {
        TYPE_ERROR("bincoef() requires two integer arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    k = GMPy_Integer_AsUnsignedLong(args[1]);
    if (k == (unsigned long)(-1) && PyErr_Occurred()) {
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    n = GMPy_Integer_AsUnsignedLong(args[0]);
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        /* Since we plan to skip the else clause and continue,
         * we need to clear the error since we aren't acting on it.
         */
        PyErr_Clear();
    }
    else {
        /* Use mpz_bin_uiui which should be faster. */
        mpz_bin_uiui(result->z, n, k);
        return (PyObject*)result;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(args[0], NULL))) {
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    mpz_bin_ui(result->z, tempx->z, k);
    Py_DECREF((PyObject*)tempx);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_isqrt,
"isqrt($module, x, /)\n--\n\n"
"Return the integer square root of a non-negative integer x.");

static PyObject *
GMPy_MPZ_Function_Isqrt(PyObject *self, PyObject *other)
{
    MPZ_Object *result;

    if (CHECK_MPZANY(other)) {
        if (mpz_sgn(MPZ(other)) < 0) {
            VALUE_ERROR("isqrt() of negative number");
            return NULL;
        }
        if ((result = GMPy_MPZ_New(NULL))) {
            mpz_sqrt(result->z, MPZ(other));
        }
    }
    else {
        if (!(result = GMPy_MPZ_From_Integer(other, NULL))) {
            TYPE_ERROR("isqrt() requires 'mpz' argument");
            return NULL;
        }
        if (mpz_sgn(result->z) < 0) {
            VALUE_ERROR("isqrt() of negative number");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        mpz_sqrt(result->z, result->z);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_isqrt_rem,
"isqrt_rem($module, x, /)\n--\n\n"
"Return a 2-element tuple (s,t) such that s=isqrt(x) and t=x-s*s.\n"
"x >=0.");

static PyObject *
GMPy_MPZ_Function_IsqrtRem(PyObject *self, PyObject *other)
{
    MPZ_Object *root = NULL, *rem = NULL, *temp = NULL;
    PyObject *result;

    if (!(temp = GMPy_MPZ_From_Integer(other, NULL))) {
        TYPE_ERROR("isqrt_rem() requires 'mpz' argument");
        return NULL;
    }

    if (mpz_sgn(temp->z) < 0) {
        VALUE_ERROR("isqrt_rem() of negative number");
        Py_DECREF((PyObject*)temp);
        return NULL;
    }

    if (!(result = PyTuple_New(2)) ||
        !(root = GMPy_MPZ_New(NULL)) ||
        !(rem = GMPy_MPZ_New(NULL))) {

        /* LCOV_EXCL_START */
        Py_DECREF((PyObject*)temp);
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)root);
        Py_XDECREF((PyObject*)rem);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_sqrtrem(root->z, rem->z, temp->z);
    Py_DECREF((PyObject*)temp);
    PyTuple_SET_ITEM(result, 0, (PyObject*)root);
    PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
    return result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_remove,
"remove($module, x, f, /)\n--\n\n"
"Return a 2-element tuple (y,m) such that x=y*(f**m) and f does\n"
"not divide y. Remove the factor f from x as many times as\n"
"possible. m is the multiplicity f in x. f > 1.");

static PyObject *
GMPy_MPZ_Function_Remove(PyObject *self, PyObject * const *args,
                         Py_ssize_t nargs)
{
    MPZ_Object *result = NULL, *tempx = NULL, *tempf = NULL;
    PyObject *x, *f;
    size_t multiplicity;

    if (nargs != 2) {
        TYPE_ERROR("remove() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    x = args[0];
    f = args[1];

    if (MPZ_Check(x) && MPZ_Check(f)) {
        if (mpz_cmp_si(MPZ(f), 2) < 0) {
            VALUE_ERROR("factor must be > 1");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        multiplicity = mpz_remove(result->z, MPZ(x), MPZ(f));
        return Py_BuildValue("(Nk)", result, multiplicity);
    }
    else {


        if (!(tempx = GMPy_MPZ_From_Integer(x, NULL)) ||
            !(tempf = GMPy_MPZ_From_Integer(f, NULL))) {

            TYPE_ERROR("remove() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempf);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        if (mpz_cmp_si(MPZ(tempf), 2) < 0) {
            VALUE_ERROR("factor must be > 1");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempf);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        multiplicity = mpz_remove(result->z, tempx->z, tempf->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempf);
        return Py_BuildValue("(Nk)", result, multiplicity);
    }
}

PyDoc_STRVAR(GMPy_doc_mpz_function_invert,
"invert($module, x, m, /)\n--\n\n"
"Return y such that x*y == 1 modulo m. Raises `ZeroDivisionError` if no\n"
"inverse exists.");

static PyObject *
GMPy_MPZ_Function_Invert(PyObject *self, PyObject * const *args,
                         Py_ssize_t nargs)
{
    PyObject *x, *y;
    MPZ_Object *result = NULL, *tempx = NULL, *tempy = NULL;
    int success;

    if (nargs != 2) {
        TYPE_ERROR("invert() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    x = args[0];
    y = args[1];

    if (MPZ_Check(x) && MPZ_Check(y)) {
        if (mpz_sgn(MPZ(y)) == 0) {
            ZERO_ERROR("invert() division by 0");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        success = mpz_invert(result->z, MPZ(x), MPZ(y));
        if (!success) {
            ZERO_ERROR("invert() no inverse exists");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
    }
    else {


        if (!(tempx = GMPy_MPZ_From_Integer(x, NULL)) ||
            !(tempy = GMPy_MPZ_From_Integer(y, NULL))) {

            TYPE_ERROR("invert() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        if (mpz_sgn(tempy->z) == 0) {
            ZERO_ERROR("invert() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF(result);
            return NULL;
        }
        success = mpz_invert(result->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        if (!success) {
            ZERO_ERROR("invert() no inverse exists");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_divexact,
"divexact($module, x, y, /)\n--\n\n"
"Return the quotient of x divided by y. Faster than standard\n"
"division but requires the remainder is zero!");

static PyObject *
GMPy_MPZ_Function_Divexact(PyObject *self, PyObject * const *args,
                           Py_ssize_t nargs)
{
    PyObject *x, *y;
    MPZ_Object *result, *tempx= NULL, *tempy = NULL;

    if (nargs != 2) {
        TYPE_ERROR("divexact() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    x = args[0];
    y = args[1];

    if (MPZ_Check(x) && MPZ_Check(y)) {
        if (mpz_sgn(MPZ(y)) == 0) {
            ZERO_ERROR("divexact() division by 0");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        mpz_divexact(result->z, MPZ(x), MPZ(y));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(x, NULL)) ||
            !(tempy = GMPy_MPZ_From_Integer(y, NULL))) {

            TYPE_ERROR("divexact() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        if (mpz_sgn(MPZ(tempy)) == 0) {
            ZERO_ERROR("divexact() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        mpz_divexact(result->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_is_square,
"is_square($module, x, /)\n--\n\n"
"Returns `True` if x is a perfect square, else return `False`.");

static PyObject *
GMPy_MPZ_Function_IsSquare(PyObject *self, PyObject *other)
{
    int res;
    MPZ_Object *tempx;

    if (MPZ_Check(other)) {
        res = mpz_perfect_square_p(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(other, NULL))) {
            TYPE_ERROR("is_square() requires 'mpz' argument");
            return NULL;
        }
        else {
            res = mpz_perfect_square_p(tempx->z);
            Py_DECREF((PyObject*)tempx);
        }
    }

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_is_integer,
"is_integer($self)\n--\n\n"
"Returns `True`.");

static PyObject *
GMPy_MPZ_Method_IsInteger(PyObject *self, PyObject *other)
{
    Py_RETURN_TRUE;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_is_square,
"is_square($self)\n--\n\n"
"Returns `True` if self is a perfect square, else return `False`.");

static PyObject *
GMPy_MPZ_Method_IsSquare(PyObject *self, PyObject *other)
{
    int res;

    res = mpz_perfect_square_p(MPZ(self));

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_is_divisible,
"is_divisible($module, x, d, /)\n--\n\n"
"Returns `True` if x is divisible by d, else return `False`.");

static PyObject *
GMPy_MPZ_Function_IsDivisible(PyObject *self, PyObject * const * args,
                              Py_ssize_t nargs)
{
    unsigned long temp;
    int res = 0;
    MPZ_Object *tempx, *tempd;

    if (nargs != 2) {
        TYPE_ERROR("is_divisible() requires 2 integer arguments");
        return NULL;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(args[0], NULL))) {
        return NULL;
    }

    temp = GMPy_Integer_AsUnsignedLong(args[1]);
    if (temp == (unsigned long)-1 && PyErr_Occurred()) {
        PyErr_Clear();
        /* Implement mpz_divisible_p here. */

        if (!(tempd = GMPy_MPZ_From_Integer(args[1], NULL))) {
            TYPE_ERROR("is_divisible() requires 2 integer arguments");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }

        res = mpz_divisible_p(tempx->z, tempd->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempd);
    }
    else {
        /* Implement mpz_divisible_ui_p here. */

        res = mpz_divisible_ui_p(tempx->z, temp);
        Py_DECREF((PyObject*)tempx);
    }
    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_is_divisible,
"is_divisible($self, d, /)\n--\n\n"
"Returns `True` if self is divisible by d, else return `False`.");

static PyObject *
GMPy_MPZ_Method_IsDivisible(PyObject *self, PyObject *other)
{
    unsigned long temp;
    int res;
    MPZ_Object *tempd;

    temp = GMPy_Integer_AsUnsignedLong(other);
    if (temp == (unsigned long)-1 && PyErr_Occurred()) {
        PyErr_Clear();
        /* Implement mpz_divisible_p here. */

        if (!(tempd = GMPy_MPZ_From_Integer(other, NULL))) {
            TYPE_ERROR("is_divisible() requires 2 integer arguments");
            return NULL;
        }

        res = mpz_divisible_p(MPZ(self), tempd->z);
        Py_DECREF((PyObject*)tempd);
    }
    else {
        /* Implement mpz_divisible_ui_p here. */

        res = mpz_divisible_ui_p(MPZ(self), temp);
    }
    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_is_congruent,
"is_congruent($module, x, y, m, /)\n--\n\n"
"Returns `True` if x is congruent to y modulo m, else return `False`.");

static PyObject *
GMPy_MPZ_Function_IsCongruent(PyObject *self, PyObject * const *args,
                              Py_ssize_t nargs)
{
    int res;
    MPZ_Object *tempx = NULL, *tempy = NULL, *tempm = NULL;

    if (nargs != 3) {
        TYPE_ERROR("is_congruent() requires 3 integer arguments");
        return NULL;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(args[0], NULL)) ||
        !(tempy = GMPy_MPZ_From_Integer(args[1], NULL)) ||
        !(tempm = GMPy_MPZ_From_Integer(args[2], NULL))) {

        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)tempy);
        Py_XDECREF((PyObject*)tempm);
        TYPE_ERROR("is_congruent() requires 3 integer arguments");
        return NULL;
    }

    res = mpz_congruent_p(tempx->z, tempy->z, tempm->z);
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)tempm);
    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_is_congruent,
"is_congruent($self, y, m, /)\n--\n\n"
"Returns `True` if self is congruent to y modulo m, else return `False`.");

static PyObject *
GMPy_MPZ_Method_IsCongruent(PyObject *self, PyObject *const *args,
                            Py_ssize_t nargs)
{
    int res;
    MPZ_Object *tempy = NULL, *tempm = NULL;

    if (nargs != 2) {
        TYPE_ERROR("is_congruent() requires 2 integer arguments");
        return NULL;
    }

    if (!(tempy = GMPy_MPZ_From_Integer(args[0], NULL)) ||
        !(tempm = GMPy_MPZ_From_Integer(args[1], NULL))) {

        Py_XDECREF((PyObject*)tempy);
        Py_XDECREF((PyObject*)tempm);
        TYPE_ERROR("is_congruent() requires 2 integer arguments");
        return NULL;
    }

    res = mpz_congruent_p(MPZ(self), tempy->z, tempm->z);
    Py_DECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)tempm);
    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_is_power,
"is_power($module, x, /)\n--\n\n"
"Return `True` if x is a perfect power (there exists a y and an\n"
"n > 1, such that x=y**n), else return `False`.");

static PyObject *
GMPy_MPZ_Function_IsPower(PyObject *self, PyObject *other)
{
    int res;
    MPZ_Object* tempx;

    if (MPZ_Check(other)) {
        res = mpz_perfect_power_p(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(other, NULL))) {
            TYPE_ERROR("is_power() requires 'mpz' argument");
            return NULL;
        }
        else {
            res = mpz_perfect_power_p(tempx->z);
            Py_DECREF((PyObject*)tempx);
        }
    }

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_is_power,
"is_power($self)\n--\n\n"
"Return `True` if self is a perfect power (there exists a y and an\n"
"n > 1, such that self=y**n), else return `False`.");

static PyObject *
GMPy_MPZ_Method_IsPower(PyObject *self, PyObject *other)
{
    int res;

    res = mpz_perfect_power_p(MPZ(self));

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_is_prime,
"is_prime($module, x, n=25, /)\n--\n\n"
"Return `True` if x is *probably* prime, else `False` if x is\n"
"definitely composite. x is checked for small divisors and up\n"
"to n Miller-Rabin tests are performed.");

static PyObject *
GMPy_MPZ_Function_IsPrime(PyObject *self, PyObject * const *args,
                          Py_ssize_t nargs)
{
    int i;
    unsigned long reps = 25;
    MPZ_Object* tempx;

    if (nargs == 0 || nargs > 2) {
        TYPE_ERROR("is_prime() requires 'mpz'[,'int'] arguments");
        return NULL;
    }

    if (nargs == 2) {
        reps = GMPy_Integer_AsUnsignedLong(args[1]);
        if (reps == (unsigned long)(-1) && PyErr_Occurred()) {
            return NULL;
        }
        /* Silently limit n to a reasonable value. */
        if (reps > 1000) {
            reps = 1000;
        }
    }

    if (!(tempx = GMPy_MPZ_From_Integer(args[0], NULL))) {
        return NULL;
    }

    if (mpz_sgn(tempx->z) == -1) {
        Py_DECREF((PyObject*)tempx);
        Py_RETURN_FALSE;
    }

    i = mpz_probab_prime_p(tempx->z, (int)reps);
    Py_DECREF((PyObject*)tempx);

    if (i)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_is_prime,
"is_prime($self, n=25, /)\n--\n\n"
"Return `True` if self is *probably* prime, else `False` if self is\n"
"definitely composite. self is checked for small divisors and up\n"
"to n Miller-Rabin tests are performed.");

static PyObject *
GMPy_MPZ_Method_IsPrime(PyObject *self, PyObject * const *args,
                        Py_ssize_t nargs)
{
    int i;
    unsigned long reps = 25;

    if (nargs > 1) {
        TYPE_ERROR("is_prime() takes at most 1 argument");
        return NULL;
    }

    if (nargs == 1) {
        reps = GMPy_Integer_AsUnsignedLong(args[0]);
        if (reps == (unsigned long)(-1) && PyErr_Occurred()) {
            return NULL;
        }
        /* Silently limit n to a reasonable value. */
        if (reps > 1000) {
            reps = 1000;
        }
    }

    if (mpz_sgn(MPZ(self)) == -1) {
        Py_RETURN_FALSE;
    }

    i = mpz_probab_prime_p(MPZ(self), (int)reps);

    if (i)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_is_probab_prime,
"is_probab_prime($module, x, n=25, /)\n--\n\n"
"Return 2 if x is definitely prime, 1 if x is probably prime,\n"
"or return 0 if x is definitely non-prime.  x is checked for small\n"
"divisors and up to n Miller-Rabin tests are performed.  Reasonable\n"
"values of n are between 15 and 50.");

static PyObject *
GMPy_MPZ_Function_IsProbabPrime(PyObject *module, PyObject *const *args,
                                Py_ssize_t nargs)
{
    int ret;
    unsigned long reps = 25;
    MPZ_Object* tempx;

    if (nargs == 0 || nargs > 2) {
        TYPE_ERROR("is_probab_prime() requires 'mpz'[,'int'] arguments");
        return NULL;
    }

    if (nargs == 2) {
        reps = PyLong_AsUnsignedLong(args[1]);
        if (reps == (unsigned long)(-1) && PyErr_Occurred()) {
            return NULL;
        }
    }

    if (!(tempx = GMPy_MPZ_From_Integer(args[0], NULL))) {
        return NULL;
    }

    if (mpz_sgn(MPZ(tempx)) == -1) {
        Py_DECREF((PyObject*)tempx);
        return PyLong_FromLong(0);
    }

    ret = mpz_probab_prime_p(MPZ(tempx), reps);
    Py_DECREF((PyObject*)tempx);

    return PyLong_FromLong(ret);
}

PyDoc_STRVAR(GMPy_doc_mpz_method_is_probab_prime,
"is_probab_prime($self, n=25, /)\n--\n\n"
"Return 2 if self is definitely prime, 1 if self is probably prime,\n"
"or return 0 if self is definitely non-prime.  self is checked for small\n"
"divisors and up to n Miller-Rabin tests are performed.  Reasonable\n"
"values of n are between 15 and 50.");

static PyObject *
GMPy_MPZ_Method_IsProbabPrime(PyObject *self, PyObject *const *args,
                              Py_ssize_t nargs)
{
    int ret;
    unsigned long reps = 25;

    if (nargs > 1) {
        TYPE_ERROR("is_probab_prime() takes at most 1 argument");
        return NULL;
    }

    if (nargs == 1) {
        reps = PyLong_AsUnsignedLong(args[0]);
        if (reps == (unsigned long)(-1) && PyErr_Occurred()) {
            return NULL;
        }
    }

    if (mpz_sgn(MPZ(self)) == -1) {
        return PyLong_FromLong(0);
    }

    ret = mpz_probab_prime_p(MPZ(self), (int)reps);

    return PyLong_FromLong(ret);
}

PyDoc_STRVAR(GMPy_doc_mpz_function_next_prime,
"next_prime($module, x, /)\n--\n\n"
"Return the next *probable* prime number > x.");

static PyObject *
GMPy_MPZ_Function_NextPrime(PyObject *self, PyObject *other)
{
    MPZ_Object *result;

    if(MPZ_Check(other)) {
        if(!(result = GMPy_MPZ_New(NULL))) {
            /* LCOV_EXCL_START */
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        mpz_nextprime(result->z, MPZ(other));
    }
    else {
        if (!(result = GMPy_MPZ_From_Integer(other, NULL))) {
            TYPE_ERROR("next_prime() requires 'mpz' argument");
            return NULL;
        }
        else {
            mpz_nextprime(result->z, result->z);
        }
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_prev_prime,
"prev_prime($module, x, /)\n--\n\n"
"Return the previous *probable* prime number < x.\n"
"Only present when compiled with GMP 6.3.0 or later.");

#if (__GNU_MP_VERSION > 6) || (__GNU_MP_VERSION == 6 &&  __GNU_MP_VERSION_MINOR >= 3)
static PyObject *
GMPy_MPZ_Function_PrevPrime(PyObject *self, PyObject *other)
{
        MPZ_Object *result;

        if(MPZ_Check(other)) {
            if(!(result = GMPy_MPZ_New(NULL))) {
                /* LCOV_EXCL_START */
                return NULL;
                /* LCOV_EXCL_STOP */
            }
            if (!mpz_prevprime(result->z, MPZ(other))) {
                /* no previous prime, raise value error. */
                VALUE_ERROR("x must be >= 3");
                return NULL;
            }
        }
        else {
            if (!(result = GMPy_MPZ_From_Integer(other, NULL))) {
                TYPE_ERROR("prev_prime() requires 'mpz' argument");
                return NULL;
            }
            else {
                if (!mpz_prevprime(result->z, result->z)) {
                    /* no previous prime, raise value error. */
                    VALUE_ERROR("x must be >= 3");
                    return NULL;
                }
            }
        }
        return (PyObject*)result;
}
#endif


PyDoc_STRVAR(GMPy_doc_mpz_function_jacobi,
"jacobi($module, x, y, /)\n--\n\n"
"Return the Jacobi symbol (x|y). y must be odd and >0.");

static PyObject *
GMPy_MPZ_Function_Jacobi(PyObject *self, PyObject *const *args,
                         Py_ssize_t nargs)
{
    MPZ_Object *tempx = NULL, *tempy = NULL;
    long res;

    if (nargs != 2) {
        TYPE_ERROR("jacobi() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(args[0], NULL)) ||
        !(tempy = GMPy_MPZ_From_Integer(args[1], NULL))) {

        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)tempy);
        return NULL;
    }

    if (mpz_sgn(tempy->z) <= 0 || mpz_even_p(tempy->z)) {
        VALUE_ERROR("y must be odd and >0");
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return NULL;
    }

    res = (long)(mpz_jacobi(tempx->z, tempy->z));
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    return PyLong_FromLong(res);
}

PyDoc_STRVAR(GMPy_doc_mpz_function_legendre,
"legendre($module, x, y, /)\n--\n\n"
"Return the Legendre symbol (x|y). y is assumed to be an odd prime.");

static PyObject *
GMPy_MPZ_Function_Legendre(PyObject *self, PyObject * const *args,
                           Py_ssize_t nargs)
{
    MPZ_Object *tempx = NULL, *tempy = NULL;
    long res;

    if (nargs != 2) {
        TYPE_ERROR("legendre() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(args[0], NULL)) ||
        !(tempy = GMPy_MPZ_From_Integer(args[1], NULL))) {

        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)tempy);
        return NULL;
    }

    if (mpz_sgn(tempy->z) <= 0 || mpz_even_p(tempy->z)) {
        VALUE_ERROR("y must be odd, prime, and >0");
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return NULL;
    }

    res = (long)(mpz_legendre(tempx->z, tempy->z));
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    return PyLong_FromLong(res);
}

PyDoc_STRVAR(GMPy_doc_mpz_function_kronecker,
"kronecker($module, x, y, /)\n--\n\n"
"Return the Kronecker-Jacobi symbol (x|y).");

static PyObject *
GMPy_MPZ_Function_Kronecker(PyObject *self, PyObject * const *args,
                            Py_ssize_t nargs)
{
    MPZ_Object *tempx = NULL, *tempy = NULL;
    long res;

    if (nargs != 2) {
        TYPE_ERROR("kronecker() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(args[0], NULL)) ||
        !(tempy = GMPy_MPZ_From_Integer(args[1], NULL))) {

        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)tempy);
        return NULL;
    }

    res = (long)(mpz_kronecker(tempx->z, tempy->z));
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    return PyLong_FromLong(res);
}

PyDoc_STRVAR(GMPy_doc_mpz_function_is_even,
"is_even($module, x, /)\n--\n\n"
"Return `True` if x is even, `False` otherwise.");

static PyObject *
GMPy_MPZ_Function_IsEven(PyObject *self, PyObject *other)
{
    int res;
    MPZ_Object *tempx;

    if (MPZ_Check(other)) {
        res = mpz_even_p(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(other, NULL))) {
            TYPE_ERROR("is_even() requires 'mpz' argument");
            return NULL;
        }
        else {
            res = mpz_even_p(tempx->z);
            Py_DECREF((PyObject*)tempx);
        }
    }

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_is_even,
"is_even($self)\n--\n\n"
"Return `True` if self is even, `False` otherwise.");

static PyObject *
GMPy_MPZ_Method_IsEven(PyObject *self, PyObject *other)
{
    int res;

    res = mpz_even_p(MPZ(self));

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_is_odd,
"is_odd($module, x, /)\n--\n\n"
"Return `True` if x is odd, `False` otherwise.");

static PyObject *
GMPy_MPZ_Function_IsOdd(PyObject *self, PyObject *other)
{
    int res;
    MPZ_Object *tempx;

    if (CHECK_MPZANY(other)) {
        res = mpz_odd_p(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(other, NULL))) {
            TYPE_ERROR("is_odd() requires 'mpz' argument");
            return NULL;
        }
        else {
            res = mpz_odd_p(tempx->z);
            Py_DECREF((PyObject*)tempx);
        }
    }

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_is_odd,
"is_odd($self)\n--\n\n"
"Return `True` if self is odd, `False` otherwise.");

static PyObject *
GMPy_MPZ_Method_IsOdd(PyObject *self, PyObject *other)
{
    int res;

    res = mpz_odd_p(MPZ(self));

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

/*
 * Add mapping support to mpz objects.
 */

static Py_ssize_t
GMPy_MPZ_Method_Length(MPZ_Object *self)
{
    return mpz_sizeinbase(self->z, 2);
}

static PyObject *
GMPy_MPZ_Method_SubScript(MPZ_Object *self, PyObject *item)
{
    if (PyIndex_Check(item)) {
        Py_ssize_t i;

        i = PyLong_AsSsize_t(item);
        if (i == -1 && PyErr_Occurred()) {
            INDEX_ERROR("argument too large to convert to an index");
            return NULL;
        }
        if (i < 0) {
            i += mpz_sizeinbase(self->z, 2);
        }
        return PyLong_FromLong(mpz_tstbit(self->z, i));
    }
    else if (PySlice_Check(item)) {
        Py_ssize_t start, stop, step, slicelength, cur, i;
        MPZ_Object *result;

        if (PySlice_GetIndicesEx(item,
                        mpz_sizeinbase(self->z, 2),
                        &start, &stop, &step, &slicelength) < 0) {
            return NULL;
        }

        if ((step < 0 && start < stop) || (step > 0 && start > stop)) {
            stop = start;
        }

        if (!(result = GMPy_MPZ_New(NULL))) {
            return NULL;
        }

        mpz_set_ui(result->z, 0);
        if (slicelength > 0) {
            for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                if(mpz_tstbit(self->z, cur)) {
                    mpz_setbit(result->z, i);
                }
            }
        }
        return (PyObject*)result;
    }
    else {
        TYPE_ERROR("bit positions must be integers");
        return NULL;
    }
}

static PyObject *
GMPy_MPZ_Attrib_GetNumer(MPZ_Object *self, void *closure)
{
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
}

static PyObject *
GMPy_MPZ_Attrib_GetReal(MPZ_Object *self, void *closure)
{
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
}

static PyObject *
GMPy_MPZ_Attrib_GetDenom(MPZ_Object *self, void *closure)
{
    MPZ_Object *result;

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_set_ui(result->z, 1);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_as_integer_ratio,
"as_integer_ratio($self)\n--\n\n"
"Return a pair of integers, whose ratio is exactly equal to the\n"
"original number.  The ratio is in lowest terms and has a\n"
"positive denominator.");
static PyObject *
GMPy_MPZ_Method_As_Integer_Ratio(PyObject *self, PyObject *args)
{
    return PyTuple_Pack(2,
                        GMPy_MPZ_Attrib_GetNumer((MPZ_Object*)self, NULL),
                        GMPy_MPZ_Attrib_GetDenom((MPZ_Object*)self, NULL));
}

PyDoc_STRVAR(GMPy_doc_mpz_method_to_bytes,
"to_bytes($self, length=1, byteorder=\'big\', *, signed=False)\n--\n\n\
Return an array of bytes representing an integer.\n\n\
  length\n\
    Length of bytes object to use.  An `OverflowError` is raised if the\n\
    integer is not representable with the given number of bytes.\n\
  byteorder\n\
    The byte order used to represent the integer.  If byteorder is\n\
    \'big\', the most significant byte is at the beginning of the byte\n\
    array.  If byteorder is \'little\', the most significant byte is at\n\
    the end of the byte array.  To request the native byte order of the\n\
    host system, use `sys.byteorder` as the byte order value.\n\
  signed\n\
    Determines whether two\'s complement is used to represent the\n\
    integer.  If signed is `False` and a negative integer is given,\n\
    an `OverflowError` is raised.");
static PyObject *
GMPy_MPZ_Method_To_Bytes(PyObject *self, PyObject *const *args,
                         Py_ssize_t nargs, PyObject *kwnames)
{
    Py_ssize_t i, nkws = 0, size, gap, length = 1;
    PyObject *bytes, *arg;
    mpz_t tmp, *px = &MPZ(self);
    char *buffer;
    int sign, is_signed = 0, is_negative, is_big, argidx[2] = {-1, -1};
    const char *byteorder = NULL, *kwname;

    if (nargs > 2) {
        TYPE_ERROR("to_bytes() takes at most 2 positional arguments");
        return NULL;
    }
    if (nargs >= 1) {
        argidx[0] = 0;
    }
    if (nargs == 2) {
        argidx[1] = 1;
    }

    if (kwnames) {
        nkws = PyTuple_GET_SIZE(kwnames);
    }
    if (nkws > 3) {
        TYPE_ERROR("to_bytes() takes at most 3 keyword arguments");
        return NULL;
    }
    for (i = 0; i < nkws; i++) {
        kwname = PyUnicode_AsUTF8(PyTuple_GET_ITEM(kwnames, i));
        if (strcmp(kwname, "signed") == 0) {
            is_signed = PyObject_IsTrue(args[nargs + i]);
        }
        else if (strcmp(kwname, "length") == 0) {
            if (nargs == 0) {
                argidx[0] = (int)(nargs + i);
            }
            else {
                TYPE_ERROR("argument for to_bytes() given by name ('length') and position (1)");
                return NULL;
            }
        }
        else if (strcmp(kwname, "byteorder") == 0) {
            if (nargs <= 1) {
                argidx[1] = (int)(nargs + i);
            }
            else {
                TYPE_ERROR("argument for to_bytes() given by name ('byteorder') and position (2)");
                return NULL;
            }
        }
        else {
            TYPE_ERROR("got an invalid keyword argument for to_bytes()");
            return NULL;
        }
    }

    if (argidx[0] >= 0) {
        arg = args[argidx[0]];
        if (PyLong_Check(arg)) {
            length = PyLong_AsSsize_t(arg);
        }
        else {
            TYPE_ERROR("to_bytes() takes an integer argument 'length'");
            return NULL;
        }
    }
    if (argidx[1] >= 0) {
        arg = args[argidx[1]];
        if (PyUnicode_Check(arg)) {
            byteorder = PyUnicode_AsUTF8(arg);
        }
        else {
            TYPE_ERROR("to_bytes() argument 'byteorder' must be str");
            return NULL;
        }
    }

    if (length < 0) {
        VALUE_ERROR("length argument must be non-negative");
        return NULL;
    }

    if (byteorder == NULL || strcmp(byteorder, "big") == 0) {
        is_big = 1;
    }
    else if (strcmp(byteorder, "little") == 0) {
        is_big = 0;
    }
    else {
        VALUE_ERROR("byteorder must be either 'little' or 'big'");
        return NULL;
    }

    is_negative = mpz_sgn(*px) < 0;

    if (is_negative) {
        if (!is_signed) {
            OVERFLOW_ERROR("can't convert negative mpz to unsigned");
            return NULL;
        }
        mpz_init(tmp);
        mpz_ui_pow_ui(tmp, 256, length);
        mpz_add(tmp, tmp, *px);
        px = &tmp;
    }

    sign = mpz_sgn(*px);
    size = mpz_sizeinbase(*px, 256) - !sign;
    gap = length - size;

    if (gap < 0 || sign < 0 ||
        (is_signed &&
#if (PY_VERSION_HEX < 0x030D08F0 || (PY_VERSION_HEX >= 0x030E0000 \
                                     && PY_VERSION_HEX < 0x030E00C3))
         length &&
#endif
         mpz_tstbit(*px, 8*length - 1) == !is_negative))
    {
        OVERFLOW_ERROR("mpz too big to convert");
        return NULL;
    }

    bytes = PyBytes_FromStringAndSize(NULL, length);
    if (bytes == NULL) {
        return NULL;
    }
    buffer = PyBytes_AS_STRING(bytes);
    memset(buffer, is_negative ? 0xFF : 0, gap);

    if ((*px)->_mp_size) {
        mpn_get_str((unsigned char *)(buffer + gap), 256,
                    (*px)->_mp_d, (*px)->_mp_size);
    }
    if (!is_big && length) {
        revstr(buffer, 0, length - 1);
    }
    if (is_negative) {
        mpz_clear(tmp);
    }

    return bytes;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_from_bytes,
"from_bytes($type, /, bytes, byteorder=\'big\', *, signed=False)\n--\n\n\
Return the integer represented by the given array of bytes.\n\n\
  bytes\n\
    Holds the array of bytes to convert.  The argument must either\n\
    support the buffer protocol or be an iterable object producing bytes.\n\
    `bytes` and `bytearray` are examples of built-in objects that support\n\
    the buffer protocol.\n\
  byteorder\n\
    The byte order used to represent the integer.  If byteorder is \'big\',\n\
    the most significant byte is at the beginning of the byte array.  If\n\
    byteorder is \'little\', the most significant byte is at the end of the\n\
    byte array.  To request the native byte order of the host system, use\n\
    `sys.byteorder` as the byte order value.\n\
  signed\n\
    Indicates whether two\'s complement is used to represent the integer.");
static PyObject *
GMPy_MPZ_Method_From_Bytes(PyTypeObject *type, PyObject *const *args, Py_ssize_t nargs, PyObject *kwnames)
{
    Py_ssize_t i, nkws = 0, length;
    PyObject *arg, *bytes;
    char *buffer;
    int is_signed = 0, endian, argidx[2] = {-1, -1};
    const char *byteorder = NULL, *kwname;
    mpz_t tmp;
    MPZ_Object *result;

    if (nargs > 2) {
        TYPE_ERROR("from_bytes() takes at most 2 positional arguments");
        return NULL;
    }
    if (nargs >= 1) {
        argidx[0] = 0;
    }
    if (nargs == 2) {
        argidx[1] = 1;
    }

    if (kwnames) {
        nkws = PyTuple_GET_SIZE(kwnames);
    }
    if (nkws > 3) {
        TYPE_ERROR("from_bytes() takes at most 3 keyword arguments");
        return NULL;
    }
    if (nkws + nargs < 1) {
        TYPE_ERROR("from_bytes() missing required argument 'bytes' (pos 1)");
        return NULL;
    }
    for (i = 0; i < nkws; i++) {
        kwname = PyUnicode_AsUTF8(PyTuple_GET_ITEM(kwnames, i));
        if (strcmp(kwname, "signed") == 0) {
            is_signed = PyObject_IsTrue(args[nargs + i]);
        }
        else if (strcmp(kwname, "bytes") == 0) {
            if (nargs == 0) {
                argidx[0] = (int)(nargs + i);
            }
            else {
                TYPE_ERROR("argument for from_bytes() given by name ('bytes') and position (1)");
                return NULL;
            }
        }
        else if (strcmp(kwname, "byteorder") == 0) {
            if (nargs <= 1) {
                argidx[1] = (int)(nargs + i);
            }
            else {
                TYPE_ERROR("argument for from_bytes() given by name ('byteorder') and position (2)");
                return NULL;
            }
        }
        else {
            TYPE_ERROR("got an invalid keyword argument for from_bytes()");
            return NULL;
        }
    }

    if (argidx[1] >= 0) {
        arg = args[argidx[1]];
        if (PyUnicode_Check(arg)) {
            byteorder = PyUnicode_AsUTF8(arg);
        }
        else {
            TYPE_ERROR("from_bytes() argument 'byteorder' must be str");
            return NULL;
        }
    }

    if (byteorder == NULL || strcmp(byteorder, "big") == 0) {
        endian = 1;
    }
    else if (strcmp(byteorder, "little") == 0) {
        endian = -1;
    }
    else {
        VALUE_ERROR("byteorder must be either 'little' or 'big'");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(NULL))) {
        return NULL;
    }

    bytes = PyObject_Bytes(args[argidx[0]]);
    if (bytes == NULL) {
        return NULL;
    }
    if (PyBytes_AsStringAndSize(bytes, &buffer, &length) == -1) {
        return NULL;
    }

    mpz_import(MPZ(result), length, endian, sizeof(char), 0, 0, buffer);
    Py_DECREF(bytes);

    if (is_signed && mpz_tstbit(MPZ(result), 8*length - 1)) {
        mpz_init(tmp);
        mpz_ui_pow_ui(tmp, 256, (mp_size_t)length);
        mpz_sub(MPZ(result), tmp, MPZ(result));
        mpz_clear(tmp);
        mpz_neg(MPZ(result), MPZ(result));
    }

    return (PyObject*)result;
}

static PyObject *
GMPy_MPZ_Attrib_GetImag(MPZ_Object *self, void *closure)
{
    MPZ_Object *result;

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_set_ui(result->z, 0);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_sizeof,
"__sizeof__($self)\n--\n\n"
"Returns the amount of memory consumed by self. Note: deleted mpz objects\n"
"are reused and may or may not be resized when a new value is assigned.");

static PyObject *
GMPy_MPZ_Method_SizeOf(PyObject *self, PyObject *other)
{
    return PyLong_FromSize_t(sizeof(MPZ_Object) + \
        (MPZ(self)->_mp_alloc * sizeof(mp_limb_t)));
}

/* Note: this particular function is also used for xmpz, mpq, and mpfr. Only
 * mpc.conjugate() does more that just return another reference to the original
 * object.
 */

PyDoc_STRVAR(GMPy_doc_mp_method_conjugate,
"conjugate($self)\n--\n\n"
"Return the conjugate of self (which is just a new reference to self since self is\n"
"not a complex number).");

static PyObject *
GMPy_MP_Method_Conjugate(PyObject *self, PyObject *args)
{
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_array,
"__array__($self, dtype=None, copy=None)\n--\n\n");

static PyObject *
GMPy_MPZ_Method_Array(PyObject *self, PyObject *const *args,
                      Py_ssize_t nargs, PyObject *kwnames)
{
    Py_ssize_t i, nkws = 0;
    int argidx[2] = {-1, -1};
    const char* kwname;
    PyObject *dtype = Py_None;

    if (nargs > 2) {
        TYPE_ERROR("__array__() takes at most 2 positional arguments");
        return NULL;
    }
    if (nargs >= 1) {
        argidx[0] = 0;
    }
    if (nargs == 2) {
        argidx[1] = 1;
    }

    if (kwnames) {
        nkws = PyTuple_GET_SIZE(kwnames);
    }
    if (nkws > 2) {
        TYPE_ERROR("__array__() takes at most 2 keyword arguments");
        return NULL;
    }
    for (i = 0; i < nkws; i++) {
        kwname = PyUnicode_AsUTF8(PyTuple_GET_ITEM(kwnames, i));
        if (strcmp(kwname, "dtype") == 0) {
            if (nargs == 0) {
                argidx[0] = (int)(nargs + i);
            }
            else {
                TYPE_ERROR("argument for __array__() given by name ('dtype') and position (1)");
                return NULL;
            }
        }
        else if (strcmp(kwname, "copy") == 0) {
            if (nargs <= 1) {
                argidx[1] = (int)(nargs + i);
            }
            else {
                TYPE_ERROR("argument for __array__() given by name ('copy') and position (2)");
                return NULL;
            }
        }
        else {
            TYPE_ERROR("got an invalid keyword argument for __array__()");
            return NULL;
        }
    }

    if (argidx[0] >= 0) {
        dtype = args[argidx[0]];
    }

    PyObject *mod = PyImport_ImportModule("numpy");

    if (!mod) {
        return NULL;
    }

    PyObject *tmp_long = GMPy_PyLong_From_MPZ((MPZ_Object *)self, NULL);

    if (!tmp_long) {
        Py_DECREF(mod);
        return NULL;
    }

    PyObject *result = PyObject_CallMethod(mod, "array",
                                           "OO", tmp_long, dtype);

    Py_DECREF(mod);
    Py_DECREF(tmp_long);

    return result;
}
