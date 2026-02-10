/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_convert.c                                                          *
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

/* This file contains all the conversion functions for gmpy2.
 *
 * Overview
 * --------
 * gmpy2 tries to optimize the performance and accuracy of conversions from
 * other numeric types. gmpy2 uses a LBYL (Look Before You Leap) approach and
 * identifies the numeric type before conversion before conversion to a gmpy2
 * type. The basic operations (+, -, *, /) are optimized to directly work with
 * some basic types such as C longs or doubles.
 */

/* ======================================================================== *
 * Miscellaneous helper functions.                                          *
 * ======================================================================== */

/* Since macros are used in gmpy2's codebase, these functions are skipped
 * until they are needed for the C API in the future.
 */

#if 0
static int GMPy_isInteger(PyObject *obj)
{
    return IS_INTEGER(obj) ? 1 : 0;
}

static int GMPy_isFraction(PyObject *obj)
{
    return (!strcmp(Py_TYPE(obj)->tp_name, "Fraction")) ? 1 : 0;
}

static int GMPy_isRational(PyObject *obj)
{
    return IS_RATIONAL(obj) ? 1 : 0;
}

static int GMPy_isReal(PyObject *obj)
{
    return IS_REAL(obj) ? 1 : 0;
}

static int GMPy_isComplex(PyObject *obj)
{
    return IS_COMPLEX(obj) ? 1 : 0;
}
#endif

/* GMPy_ObjectType(PyObject *obj) returns an integer that identifies the
 * object's type. See gmpy2_convert.h for details.
 * 
 * Exceptions are never raised.
 */

static inline int GMPy_ObjectType(PyObject *obj)
{
    /* Tests are sorted by order by (best guess of) most common argument type.
     * Tests that require attribute lookups are done last.
     */

    if (MPZ_Check(obj)) return OBJ_TYPE_MPZ;

    if (MPFR_Check(obj)) return OBJ_TYPE_MPFR;

    if (MPC_Check(obj))  return OBJ_TYPE_MPC;
    
    if (MPQ_Check(obj))  return OBJ_TYPE_MPQ;

    if (XMPZ_Check(obj)) return OBJ_TYPE_XMPZ;

    if (PyLong_Check(obj)) return OBJ_TYPE_PyInteger;

    if (PyFloat_Check(obj)) return OBJ_TYPE_PyFloat;

    if (PyComplex_Check(obj)) return OBJ_TYPE_PyComplex;

    if (IS_FRACTION(obj)) return OBJ_TYPE_PyFraction;

    /* Now we look for the presence of __mpz__, __mpq__, __mpfr__, and __mpc__.
     * Since a type may define more than one of the special methods, we perform
     * the checks in reverse order.
     */

    if (HAS_MPC_CONVERSION(obj)) return OBJ_TYPE_HAS_MPC;

    if (HAS_MPFR_CONVERSION(obj)) return OBJ_TYPE_HAS_MPFR;

    if (HAS_MPQ_CONVERSION(obj)) return OBJ_TYPE_HAS_MPQ;

    if (HAS_MPZ_CONVERSION(obj)) return OBJ_TYPE_HAS_MPZ;

    return OBJ_TYPE_UNKNOWN;
}

static PyObject *
GMPy_RemoveIgnoredASCII(PyObject *s)
{
    PyObject *ascii_str = NULL, *temp = NULL, *filtered = NULL, *symbol = NULL, *blank = NULL;

    if (PyBytes_CheckExact(s)) {
        temp = PyUnicode_DecodeASCII(PyBytes_AS_STRING(s), PyBytes_GET_SIZE(s), "strict");
        if (!temp) {
            VALUE_ERROR("string contains non-ASCII characters");
            return NULL;
        }
    }
    else if (PyUnicode_Check(s)) {
        Py_INCREF(s);
        temp = s;
    }
    else {
        /* LCOV_EXCL_START */
        TYPE_ERROR("object is not string or Unicode");
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    blank = PyUnicode_FromString("");

    symbol = PyUnicode_FromString(" ");
    filtered = PyUnicode_Replace(temp, symbol, blank, -1);
    Py_XDECREF(symbol);
    Py_XDECREF(temp);

    temp = filtered;

    symbol = PyUnicode_FromString("_");
    filtered = PyUnicode_Replace(temp, symbol, blank, -1);
    Py_XDECREF(symbol);
    Py_XDECREF(temp);

    Py_XDECREF(blank);

    if (!filtered) {
        return NULL;
    }

    ascii_str = PyUnicode_AsASCIIString(filtered);
    Py_DECREF(filtered);
    if (!ascii_str) {
        VALUE_ERROR("string contains non-ASCII characters");
        return NULL;
    }

    return ascii_str;
}

/* mpz_set_PyStr converts a Python "string" into a mpz_t structure. It accepts
 * a sequence of bytes (i.e. str in Python 2, bytes in Python 3) or a Unicode
 * string (i.e. unicode in Python 3, str in Python 3). Returns -1 on error,
 * 1 if successful.
 */

static int
mpz_set_PyStr(mpz_t z, PyObject *s, int base)
{
    char *cp;
    int negative = 0;
    PyObject *ascii_str;

    ascii_str = GMPy_RemoveIgnoredASCII(s);

    if (!ascii_str) return -1;

    cp = PyBytes_AsString(ascii_str);

    if (cp[0] == '+') cp++;
    if (cp[0] == '-') {
        cp++;
        negative = 1;
    }

    /* Check for leading base indicators. */
    if (cp[0] == '0' && cp[1] != '\0') {
        if (base == 0) {
            /* GMP uses prefix '0' for octal, so set base here. */
            if (tolower(cp[1]) == 'o') {
                base = 8;
                cp += 2;
            }
            else if (tolower(cp[1]) != 'b' && tolower(cp[1]) != 'x') {
                base = 10;
            }
        }
        else {
            /* If the specified base matches the leading base indicators,
             * then we need to skip the base indicators.
             */
            if ((tolower(cp[1]) == 'b' && base ==  2) ||
                (tolower(cp[1]) == 'o' && base ==  8) ||
                (tolower(cp[1]) == 'x' && base == 16))
            {
                cp += 2;
            }
        }
    }

    while (cp[0] == '0' && cp[1] != '\0' && base != 0) cp++;

    /* delegate rest to GMP's function */
    if (-1 == mpz_set_str(z, cp, base)) {
        VALUE_ERROR("invalid digits");
        Py_DECREF(ascii_str);
        return -1;
    }
    if (negative) {
        mpz_neg(z, z);
    }
    Py_DECREF(ascii_str);
    return 1;
}

/* Format an mpz into any base (2 to 62). Bits in the option parameter
 * control various behaviors:
 *   bit 0: if set, output is wrapped with mpz(...) or xmpz(...)
 *   bit 1: if set, a '+' is included for positive numbers
 *   bit 2: if set, a ' ' is included for positive nubmers
 *   bit 3: if set, a '0b', '0o', or '0x' is included for binary, octal, hex
 *   bit 4: if set, no prefix is included for binary, octal, hex
 *
 * Note: if neither bit 3 or 4 is set, prefixes that match the platform
 * default are included.
 *
 * If base < 0, capital letters are used.
 *
 * If which = 0, then mpz formatting is used (if bit 0 set). Otherwise xmpz
 * formatting is used (if bit 0 is set).
 */

static char* _ztag = "mpz(";
static char* _xztag = "xmpz(";

static PyObject *
mpz_ascii(mpz_t z, int base, int option, int which)
{
    PyObject *result;
    char *buffer, *p;
    int negative = 0;
    size_t size;

    if (
        !(
          ((base >= -36) && (base <= -2)) ||
          ((base >= 2) && (base <= 62))
         )
       ) {
        VALUE_ERROR("base must be in the interval 2 ... 62");
        return NULL;
    }

    /* Allocate extra space for:
     *
     * minus sign and trailing NULL byte (2)
     * 'xmpz()' tag                      (6)
     * '0x' prefix                       (2)
     *                                  -----
     *                                   10
     *
     * And add one more to be sure...
     */

    size = mpz_sizeinbase(z, (base < 0 ? -base : base)) + 11;
    TEMP_ALLOC(buffer, size);

    if (mpz_sgn(z) < 0) {
        negative = 1;
        mpz_neg(z, z);
    }

    p = buffer;
    if (option & 1) {
        if (which)
            strcpy(p, _xztag);
        else
            strcpy(p, _ztag);
        p += strlen(p);
    }

    if (negative) {
        *(p++) = '-';
    }
    else {
        if (option & 2)
            *(p++) = '+';
        else if (option & 4)
            *(p++) = ' ';
    }

    if (option & 8) {
        if (base == 2)        { *(p++) = '0'; *(p++) = 'b'; }
        else if (base == 8)   { *(p++) = '0'; *(p++) = 'o'; }
        else if (base == 16)  { *(p++) = '0'; *(p++) = 'x'; }
        else if (base == -16) { *(p++) = '0'; *(p++) = 'X'; }
    }
    else if (!(option & 24)) {
        if (base == 2)        { *(p++) = '0'; *(p++) = 'b'; }
        else if (base == 8)   { *(p++) = '0'; *(p++) = 'o'; }
        else if (base == 16)  { *(p++) = '0'; *(p++) = 'x'; }
        else if (base == -16) { *(p++) = '0'; *(p++) = 'X'; }
    }

    /* Call GMP. */
    mpz_get_str(p, base, z);
    p = buffer + strlen(buffer);

    if (option & 1)
        *(p++) = ')';
    *(p++) = '\00';

    result = PyUnicode_FromString(buffer);
    if (negative == 1) {
        mpz_neg(z, z);
    }
    TEMP_FREE(buffer, size);
    return result;
}
