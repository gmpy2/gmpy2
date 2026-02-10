/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_binary.c                                                          *
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

/* Conversion routines between GMPY2 objects and a compact, portable
 * binary representation. The binary format of GMPY2 is not compatible
 * with GMPY 1.x. Methods to read the old format are provided.
 */

/* Provide functions to access the old binary formats. */

PyDoc_STRVAR(doc_mpz_from_old_binary,
"mpz_from_old_binary($module, string, /)\n--\n\n"
"Return an 'mpz' from a GMPY 1.x binary format.");

static PyObject *
GMPy_MPZ_From_Old_Binary(PyObject *self, PyObject *other)
{
    unsigned char *cp;
    Py_ssize_t len;
    int negative = 0;
    MPZ_Object *result;

    if (!(PyBytes_Check(other))) {
        TYPE_ERROR("mpz_from_old_binary() requires bytes argument");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    len = PyBytes_Size(other);
    cp = (unsigned char*)PyBytes_AsString(other);

    if (cp[len-1] == 0xFF) {
        negative = 1;
        --len;
    }
    mpz_import(result->z, len, -1, sizeof(char), 0, 0, cp);
    if (negative)
        mpz_neg(result->z, result->z);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpq_from_old_binary,
"mpq_from_old_binary($module, string, /)\n--\n\n"
"Return an 'mpq' from a GMPY 1.x binary format.");

static PyObject *
GMPy_MPQ_From_Old_Binary(PyObject *self, PyObject *other)
{
    unsigned char *cp;
    Py_ssize_t len;
    int topper, negative, numlen;
    mpz_t numerator, denominator;
    MPQ_Object *result;

    if (!(PyBytes_Check(other))) {
        TYPE_ERROR("mpq_from_old_binary() requires bytes argument");
        return NULL;
    }

    if (!(result = GMPy_MPQ_New(NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    len = PyBytes_Size(other);
    cp = (unsigned char*)PyBytes_AsString(other);

    if (len < 6) {
        VALUE_ERROR("invalid mpq binary (too short)");
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    topper = cp[3] & 0x7f;
    negative = cp[3] & 0x80;
    numlen = cp[0] + 256 * (cp[1] + 256 * (cp[2] + 256 * topper));
    if (len < (4 + numlen + 1)) {
        VALUE_ERROR("invalid mpq binary (num len)");
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    mpz_init(numerator);
    mpz_init(denominator);
    mpz_import(numerator, numlen, -1, sizeof(char), 0, 0, cp+4);
    mpz_import(denominator, len-4-numlen, -1, sizeof(char), 0, 0, cp+4+numlen);
    if (negative)
        mpz_neg(numerator, numerator);

    mpq_set_num(result->q, numerator);
    mpq_set_den(result->q, denominator);
    mpq_canonicalize(result->q);
    mpz_clear(numerator);
    mpz_clear(denominator);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpfr_from_old_binary,
"mpfr_from_old_binary($module, string, /)\n--\n\n"
"Return an `mpfr` from a GMPY 1.x binary mpf format.");

static PyObject *
GMPy_MPFR_From_Old_Binary(PyObject *self, PyObject *other)
{
    unsigned char *cp;
    Py_ssize_t len;
    MPFR_Object *result;
    mpfr_t digit;
    mpfr_prec_t prec;
    int i, codebyte, resusign, exposign, resuzero, precilen;
    unsigned int expomag = 0;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (!(PyBytes_Check(other))) {
        TYPE_ERROR("mpfr_from_old_binary() requires bytes argument");
        return NULL;
    }

    len = PyBytes_Size(other);
    cp = (unsigned char*)PyBytes_AsString(other);

    if (len == 1) {
        prec = 0;
    }
    else {
        prec = (mpfr_prec_t)(8 * (len - 5));
        if ((len>=5) && (cp[0]&8)) {
            prec = 0;
            for (i=4; i>0; --i) {
                prec = (prec << 8) | cp[i];
            }
        }
    }

    /*
     * binary format for MP floats: first, a code-byte, then, a LSB
     * 4-byte unsigned int (exponent magnitude), then the "mantissa"
     * (actually, "significand", but "mantissa" is the usual term...)
     * in MSB form.
     *
     * The codebyte encodes both the signs, exponent and result, or
     * also the zeroness of the result (in which case, nothing more).
     */
    codebyte = cp[0];
    resusign = codebyte & 1;
    exposign = codebyte & 2;
    resuzero = codebyte & 4;
    precilen = (codebyte & 8)?4:0;

    /* mpfr zero has a very compact (1-byte) binary encoding!-) */
    if (resuzero) {
        if (!(result = GMPy_MPFR_New(prec, context))) {
            /* LCOV_EXCL_START */
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        result->rc = mpfr_set_ui(result->f, 0, MPFR_RNDN);
        return (PyObject*)result;
    }

    /* all other numbers are 6+ bytes: codebyte, 4-byte exp, 1+
     * bytes for the mantissa; check this string is 6+ bytes
     */
    if (len < 6 + precilen) {
        VALUE_ERROR("invalid mpf binary encoding (too short)");
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(prec, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    /* reconstruct exponent */
    for (i = 4 + precilen; i > precilen; --i) {
        expomag = (expomag<<8) | cp[i];
    }

    /* reconstruct 'mantissa' (significand) */
    mpfr_set_si(result->f, 0, MPFR_RNDN);
    mpfr_init2(digit, prec);
    for (i = 5 + precilen; i<len; i++) {
        mpfr_set_ui(digit, cp[i], MPFR_RNDN);
        mpfr_div_2ui(digit, digit, (unsigned long)((i-4-precilen) * 8),
                     context->ctx.mpfr_round);
        mpfr_add(result->f, result->f, digit, MPFR_RNDN);
    }
    mpfr_clear(digit);
    /* apply exponent, with its appropriate sign */
    if (exposign)
        mpfr_div_2ui(result->f, result->f, 8*expomag, MPFR_RNDN);
    else
        mpfr_mul_2ui(result->f, result->f, 8*expomag, MPFR_RNDN);
    /* apply significand-sign (sign of the overall number) */
    if (resusign)
        mpfr_neg(result->f, result->f, MPFR_RNDN);

    return (PyObject*)result;
}

/* Format of the binary representation of an mpz/xmpz.
 *
 * byte[0]:     1 => mpz
 *              2 => xmpz
 *              3 => mpq  (see Pympq_To_Binary)
 *              4 => mpfr (see Pympfr_To_Binary)
 *              5 => mpc  (see Pympc_To_Binary)
 * byte[1:0-1]: 0 => value is 0
 *              1 => value is > 0
 *              2 => value is < 0
 *              3 => unassigned
 * byte[2]+: value
 */

static PyObject *
GMPy_MPZ_To_Binary(MPZ_Object *self)
{
    size_t size = 2;
    int sgn;
    char *buffer;
    PyObject *result;

    sgn = mpz_sgn(self->z);
    if (sgn == 0) {
        TEMP_ALLOC(buffer, size);
        buffer[0] = 0x01;
        buffer[1] = 0x00;
        goto done;
    }

    size = mpz_sizeinbase(self->z, 256) + 2;

    TEMP_ALLOC(buffer, size);
    buffer[0] = 0x01;
    if (sgn > 0)
        buffer[1] = 0x01;
    else
        buffer[1] = 0x02;
    mpn_get_str((unsigned char *)(buffer + 2), 256,
                self->z->_mp_d, Py_ABS(self->z->_mp_size));
    revstr(buffer, 2, size - 1);

  done:
    result = PyBytes_FromStringAndSize(buffer, size);
    TEMP_FREE(buffer, size);
    return result;
}

static PyObject *
GMPy_XMPZ_To_Binary(XMPZ_Object *self)
{
    size_t size = 2;
    int sgn;
    char *buffer;
    PyObject *result;

    sgn = mpz_sgn(self->z);
    if (sgn == 0) {
        TEMP_ALLOC(buffer, size);
        buffer[0] = 0x02;
        buffer[1] = 0x00;
        goto done;
    }

    size = ((mpz_sizeinbase(self->z, 2) + 7) / 8) + 2;

    TEMP_ALLOC(buffer, size);
    buffer[0] = 0x02;
    if (sgn > 0)
        buffer[1] = 0x01;
    else
        buffer[1] = 0x02;
    mpz_export(buffer+2, NULL, -1, sizeof(char), 0, 0, self->z);

  done:
    result = PyBytes_FromStringAndSize(buffer, size);
    TEMP_FREE(buffer, size);
    return result;
}

/* Format of the binary representation of an mpq.
 *
 * byte[0]:     1 => mpz  (see Pympz_To_Binary)
 *              2 => xmpz (see Pyxmpz_To_Binary)
 *              3 => mpq
 *              4 => mpfr (see Pympfr_To_Binary)
 *              5 => mpc  (see Pympc_To_Binary)
 * byte[1:0-1]: 0 => value is 0
 *              1 => value is > 0
 *              2 => value is < 0
 *              3 => unassigned
 * byte[1:2-2]: 0 => 32-bit length (n=4)
 *              1 => 64-bit length (n=8)
 * byte[2+]:    numerator length, using either 4 or 8 bytes
 * byte[2+n]+:  numerator, followed by denominator
 */

static PyObject *
GMPy_MPQ_To_Binary(MPQ_Object *self)
{
    size_t sizenum, sizeden, sizesize = 4, size = 2, sizetemp, i;
    size_t count = 0;
    int sgn;
    char *buffer, large = 0x00;
    PyObject *result = 0;

    sgn = mpq_sgn(self->q);
    if (sgn == 0) {
        TEMP_ALLOC(buffer, size);
        buffer[0] = 0x03;
        buffer[1] = 0x00;
        goto done;
    }

    sizenum = (mpz_sizeinbase(mpq_numref(self->q), 2) + 7) / 8;
    sizeden = (mpz_sizeinbase(mpq_denref(self->q), 2) + 7) / 8;
    size = sizenum + sizeden + 2;

    /* Check if sizenum larger than 32 bits. */
    if ((sizenum >> 16) >> 16) {
        /* Current versions of GMP do not allow values to be this large. The
         * test is left to (possibly) support future versions that support
         * larger values.
         */
        /* LCOV_EXCL_START */
        large = 0x04;
        sizesize = 8;
        /* LCOV_EXCL_STOP */
    }
    size += sizesize;

    TEMP_ALLOC(buffer, size);
    buffer[0] = 0x03;
    if (sgn > 0)
        buffer[1] = 0x01 | large;
    else
        buffer[1] = 0x02 | large;

    /* Copy sizenum to the buffer. */
    sizetemp = sizenum;
    for (i=0; i<sizesize; i++) {
        buffer[i+2] = (char)(sizetemp & 0xff);
        sizetemp >>= 8;
    }

    mpz_export(buffer+sizesize+2, &count, -1,
               sizeof(char), 0, 0, mpq_numref(self->q));
    if (count != sizenum) {
        /* LCOV_EXCL_START */
        SYSTEM_ERROR("internal error in Pympq_To_Binary");
        TEMP_FREE(buffer, size);
        return NULL;
        /* LCOV_EXCL_STOP */
    }
    count = 0;
    mpz_export(buffer+sizenum+sizesize+2, &count, -1,
               sizeof(char), 0, 0, mpq_denref(self->q));
    if (count != sizeden) {
        /* LCOV_EXCL_START */
        SYSTEM_ERROR("internal error in Pympq_To_Binary");
        TEMP_FREE(buffer, size);
        return NULL;
        /* LCOV_EXCL_STOP */
    }


  done:
    result = PyBytes_FromStringAndSize(buffer, size);
    TEMP_FREE(buffer, size);
    return result;
}

/* Format of the binary representation of an mpfr.
 *
 * byte[0]:      1 => mpz  (see Pympz_To_Binary)
 *               2 => xmpz (see Pyxmpz_To_Binary)
 *               3 => mpq  (see Pympq_To_Binary)
 *               4 => mpfr
 *               5 => mpc  (see Pympc_To_Binary)
 * byte[1:0]:    0 => value is "special"
 *               1 => value is an actual number
 * byte[1:1]:    0 => signbit is clear
 *               1 => signbit is set
 * byte[1:2-2]:  0 => 32-bit lengths (n=4)
 *               1 => 64-bit lengths (n=8)
 * byte[1:3-4]:  0 => 0 (see signbit)
 *               1 => value is NaN
 *               2 => value is Inf (see signbit)
 *               3 => unassigned
 * byte[1:5]:    0 => exponent is positive
 *               1 => exponent is negative
 * byte[1:6]:    0 => 4 byte limbs
 *               1 => 8 byte limbs
 * byte[2]:      0 => rc = 0
 *               1 => rc > 0
 *               2 => rc < 0
 * byte[3]:      mpfr.round_mode
 * byte[4]+:     precision, saved in 4 or 8 bytes
 * byte[4+n]+:   exponent, saved in 4 or 8 bytes
 * byte[4+2n]+:  mantissa
 */

static PyObject *
GMPy_MPFR_To_Binary(MPFR_Object *self)
{
    size_t sizemant = 0, sizesize = 4, size = 4, sizetemp, i;
    mp_limb_t templimb;
    mpfr_prec_t precision;
    mpfr_exp_t exponent = 0;
    int sgn;
    char *buffer, *cp, large = 0x00, expsgn = 0x00;
    PyObject *result = 0;

    /* Check if the precision, exponent and mantissa length can fit in
     * 32 bits.
     */

    sgn = mpfr_signbit(self->f);
    precision = mpfr_get_prec(self->f);

    /* Exponent and mantiss are only valid for regular numbers
     * (not 0, Nan, Inf, -Inf).
     */
    if (mpfr_regular_p(self->f)) {
        exponent = self->f->_mpfr_exp;
        if (exponent < 0) {
            exponent = -exponent;
            expsgn = 0x20;
        }
        /* Calculate the size of mantissa in limbs */
        sizemant = (self->f->_mpfr_prec + mp_bits_per_limb - 1)/mp_bits_per_limb;
    }
    if (((exponent >> 16) >> 16) ||
        ((precision >> 16) >> 16) ||
        ((sizemant >> 16) >> 16)) {
        /* This can only be tested on 64-bit platforms. lcov will report the
         * code as not tested until 64-bit specific tests are created.
         */
        sizesize = 8;
        large = 0x04;
    }

    if (!mpfr_regular_p(self->f)) {
        /* Only need to save the precision. */
        size += sizesize;
        TEMP_ALLOC(buffer, size);
        buffer[0] = 0x04;

        /* Set to all 0 since we are special. */
        buffer[1] = 0x00;

        /* Set the sign bit. */
        if (sgn) buffer[1] |= 0x02;

        /* 4 or 8 byte values. */
        buffer[1] |= large;

        /* Check if NaN. */
        if (mpfr_nan_p(self->f)) buffer[1] |= 0x08;

        /* Check if Infinity. */
        if (mpfr_inf_p(self->f)) buffer[1] |= 0x10;

        /* Save the result code */
        if (self->rc == 0)     buffer[2] = 0x00;
        else if (self->rc > 0) buffer[2] = 0x01;
        else                   buffer[2] = 0x02;

        /* Save the precision */
        sizetemp = precision;
        for (i=0; i<sizesize; i++) {
            buffer[i+4] = (char)(sizetemp & 0xff);
            sizetemp >>= 8;
        }
        goto done;
    }

    /* Now process all actual numbers. */
    size += (2 * sizesize) + (sizemant * (mp_bits_per_limb >> 3));
    TEMP_ALLOC(buffer, size);
    buffer[0] = 0x04;

    /* Set bit 0 to 1 since we are an actual number. */
    buffer[1] = 0x01;

    /* Save the sign bit. */
    if (sgn) buffer[1] |= 0x02;

    /* Save the size of the values. */
    buffer[1] |= large;

    /* Save the exponent sign. */
    buffer[1] |= expsgn;

    /* Save the limb size. */
    if ((mp_bits_per_limb >> 3) == 8)
        buffer[1] |= 0x40;
    /* This branch can only be reached on 32-bit platforms. */
    else if ((mp_bits_per_limb >> 3) != 4) {
        /* LCOV_EXCL_START */
        SYSTEM_ERROR("cannot support current limb size");
        TEMP_FREE(buffer, size);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    /* Save the result code. */
    if (self->rc == 0)     buffer[2] = 0x00;
    else if (self->rc > 0) buffer[2] = 0x01;
    else                   buffer[2] = 0x02;

    /* Rounding mode is no longer used, so just store a null byte. */
    buffer[3] = 0x00;

    /* Save the precision */
    cp = buffer + 4;
    sizetemp = precision;
    for (i=0; i<sizesize; i++) {
        cp[i] = (char)(sizetemp & 0xff);
        sizetemp >>= 8;
    }

    /* Save the exponenet */
    cp += sizesize;
    sizetemp = exponent;
    for (i=0; i<sizesize; i++) {
        cp[i] = (char)(sizetemp & 0xff);
        sizetemp >>= 8;
    }

    /* Save the actual mantissa */
    cp += sizesize;
    for (i=0; i<sizemant; i++) {
        templimb = self->f->_mpfr_d[i];
#if GMP_LIMB_BITS == 64
        cp[0] = (char)(templimb & 0xff);
        templimb >>= 8;
        cp[1] = (char)(templimb & 0xff);
        templimb >>= 8;
        cp[2] = (char)(templimb & 0xff);
        templimb >>= 8;
        cp[3] = (char)(templimb & 0xff);
        templimb >>= 8;
        cp[4] = (char)(templimb & 0xff);
        templimb >>= 8;
        cp[5] = (char)(templimb & 0xff);
        templimb >>= 8;
        cp[6] = (char)(templimb & 0xff);
        templimb >>= 8;
        cp[7] = (char)(templimb & 0xff);
        cp += 8;
#endif
#if GMP_LIMB_BITS == 32
        cp[0] = (char)(templimb & 0xff);
        templimb >>= 8;
        cp[1] = (char)(templimb & 0xff);
        templimb >>= 8;
        cp[2] = (char)(templimb & 0xff);
        templimb >>= 8;
        cp[3] = (char)(templimb & 0xff);
        cp += 4;
#endif
}

  done:
    result = PyBytes_FromStringAndSize(buffer, size);
    TEMP_FREE(buffer, size);
    return result;
}

/* Format of the binary representation of an mpc.
 *
 * The format consists of the concatenation of mpfrs (real and imaginary)
 * converted to binary format. The 0x04 leading byte of each binary string
 * is replaced by 0x05.
 */

static PyObject *
GMPy_MPC_To_Binary(MPC_Object *obj)
{
    MPFR_Object *real = NULL, *imag = NULL;
    PyObject *result = NULL, *temp = NULL;
    mpfr_prec_t rprec = 0, iprec = 0;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    mpc_get_prec2(&rprec, &iprec, obj->c);

    if (!(real = GMPy_MPFR_New(rprec, context)) ||
        !(imag = GMPy_MPFR_New(iprec, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)real);
        Py_XDECREF((PyObject*)imag);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_set(real->f, mpc_realref(obj->c), MPFR_RNDN);
    mpfr_set(imag->f, mpc_imagref(obj->c), MPFR_RNDN);
    real->rc = obj->rc;

    if (!(result = GMPy_MPFR_To_Binary(real)) ||
        !(temp = GMPy_MPFR_To_Binary(imag))) {
        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)temp);
        Py_DECREF((PyObject*)real);
        Py_DECREF((PyObject*)imag);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    Py_DECREF((PyObject*)real);
    Py_DECREF((PyObject*)imag);

    Py_ssize_t result_size, temp_size;
    char *result_str, *temp_str;

    if ((PyBytes_AsStringAndSize(result, &result_str, &result_size) < 0)
        || (PyBytes_AsStringAndSize(temp, &temp_str, &temp_size) < 0))
    {
        /* LCOV_EXCL_START */
        Py_DECREF(result);
        Py_DECREF(temp);
        return NULL;
        /* LCOV_EXCL_STOP */
    }
    result_str[0] = 0x05;
    temp_str[0] = 0x05;

    char *buf = PyMem_Malloc(result_size + temp_size);

    if (!buf) {
        /* LCOV_EXCL_START */
        Py_DECREF(result);
        Py_DECREF(temp);
        return NULL;
        /* LCOV_EXCL_STOP */
    }
    memcpy(buf, result_str, result_size);
    memcpy(buf + result_size, temp_str, temp_size);
    Py_DECREF(result);
    Py_DECREF(temp);

    PyObject *ret = PyBytes_FromStringAndSize(buf, result_size + temp_size);

    PyMem_Free(buf);
    return ret;
}

PyDoc_STRVAR(doc_from_binary,
"from_binary($module, bytes, /)\n--\n\n"
"Return a Python object from a byte sequence created by `to_binary()`.");

static PyObject *
GMPy_MPANY_From_Binary(PyObject *self, PyObject *other)
{
    unsigned char *buffer, *cp;
    Py_ssize_t len;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (!(PyBytes_Check(other))) {
        TYPE_ERROR("from_binary() requires bytes argument");
        return NULL;
    }

    len = PyBytes_Size(other);
    if (len < 2) {
        VALUE_ERROR("byte sequence too short for from_binary()");
        return NULL;
    }
    buffer = (unsigned char*)PyBytes_AsString(other);
    cp = buffer;

    switch (cp[0]) {
        case 0x01: {
            MPZ_Object *result;

            if (!(result = GMPy_MPZ_New(NULL))) {
                /* LCOV_EXCL_START */
                return NULL;
                /* LCOV_EXCL_STOP */
            }
            if (cp[1] == 0x00) {
                mpz_set_ui(result->z, 0);
                return (PyObject*)result;
            }
            mpz_import(result->z, len-2, -1, sizeof(char), 0, 0, cp+2);
            if (cp[1] == 0x02)
                mpz_neg(result->z, result->z);
            return (PyObject*)result;
            break;
        }
        case 0x02: {
            XMPZ_Object *result;

            if (!(result = GMPy_XMPZ_New(NULL))) {
                /* LCOV_EXCL_START */
                return NULL;
                /* LCOV_EXCL_STOP */
            }
            if (cp[1] == 0x00) {
                mpz_set_ui(result->z, 0);
                return (PyObject*)result;
            }
            mpz_import(result->z, len-2, -1, sizeof(char), 0, 0, cp+2);
            if (cp[1] == 0x02)
                mpz_neg(result->z, result->z);
            return (PyObject*)result;
            break;
        }
        case 0x03: {
            MPQ_Object *result;
            Py_ssize_t numlen = 0, sizesize = 4, i;
            mpz_t num, den;

            if (!(result = GMPy_MPQ_New(NULL))) {
                /* LCOV_EXCL_START */
                return NULL;
                /* LCOV_EXCL_STOP */
            }
            if (cp[1] == 0x00) {
                mpq_set_ui(result->q, 0, 1);
                return (PyObject*)result;
            }
            if (cp[1] & 0x04)
                sizesize = 8;

            if (len < 2 + sizesize) {
                VALUE_ERROR("byte sequence too short for from_binary()");
                return NULL;
            }

            for (i=sizesize; i>0; --i) {
                numlen = (numlen << 8) + cp[i+1];
            }

            if (len < 2 + sizesize + numlen + 1) {
                VALUE_ERROR("byte sequence too short for from_binary()");
                return NULL;
            }

            mpz_init(num);
            mpz_init(den);
            mpz_import(num, numlen, -1,
                       sizeof(char), 0, 0, cp+sizesize+2);
            mpz_import(den, len-numlen-sizesize-2, -1,
                       sizeof(char), 0, 0, cp+sizesize+numlen+2);
            mpq_set_num(result->q, num);
            mpq_set_den(result->q, den);
            mpq_canonicalize(result->q);
            mpz_clear(num);
            mpz_clear(den);

            if (cp[1] == 0x02)
                mpq_neg(result->q, result->q);
            return (PyObject*)result;
            break;
        }
        case 0x04: {
            MPFR_Object *result;
            Py_ssize_t sizemant = 0, sizesize = 4, i, newmant;
            mpfr_prec_t precision = 0;
            mpfr_exp_t exponent = 0;
            mp_limb_t templimb;
            int sgn = 1, expsgn = 1, limbsize = 4;
            int newlimbsize = (mp_bits_per_limb >> 3);

            if (len < 4) {
                VALUE_ERROR("byte sequence too short for from_binary()");
                return NULL;
            }

            /* Get size of values. */
            if (cp[1] & 0x04) sizesize = 8;

            /* Get the original precision. */
            for (i=sizesize; i>0; --i) {
                precision = (precision << 8) + cp[i+3];
            }

            /* Get the original sign bit. */
            if (cp[1] & 0x02) sgn = -1;

            /* Get the original exponent sign. */
            if (cp[1] & 0x20) expsgn = -1;

            /* Get the limb size of the originating system. */
            if (cp[1] & 0x40) limbsize = 8;


            if (!(result = GMPy_MPFR_New(precision, context))) {
                /* LCOV_EXCL_START */
                return NULL;
                /* LCOV_EXCL_STOP */
            }

            /* Restore the original result code and rounding mode. */

            /* Get the original result code. */
            if (cp[2] == 0)      result->rc = 0;
            else if (cp[2] == 1) result->rc = 1;
            else                 result->rc = -1;

            if (!(cp[1] & 0x01)) {
                /* Process special numbers. */
                if ((cp[1] & 0x18) == 0x00)
                    mpfr_set_zero(result->f, sgn);
                else if ((cp[1] & 0x18) == 0x08)
                    mpfr_set_nan(result->f);
                else
                    mpfr_set_inf(result->f, sgn);
                return (PyObject*)result;
            }
            /* Process actual numbers. */

            /* Calculate the number of limbs on the original system. */
            if (limbsize == 8) sizemant = ((precision + 63) / 64);
            else               sizemant = ((precision + 31) / 32);

            /* Calculate the number of limbs on the current system. */
            newmant = (precision + mp_bits_per_limb - 1) / mp_bits_per_limb;

            /* Get the original exponent. */
            cp = buffer + 4 + sizesize - 1;
            for (i=sizesize; i>0; --i) {
                exponent = (exponent << 8) + cp[i];
            }

            if (len < 2 + sizesize) {
                VALUE_ERROR("byte sequence too short for from_binary()");
                return NULL;
            }

            /* Check if the mantissa occupies the same number of bytes
             * on both the source and target system. */
            if (limbsize * sizemant == newmant * newlimbsize) {
                mpfr_set_ui(result->f, 1, MPFR_RNDN);
                cp = buffer + 4 + (2 * sizesize);
                for (i=0; i<newmant; i++) {
#if GMP_LIMB_BITS == 64
                    templimb = cp[7];
                    templimb = (templimb << 8) + cp[6];
                    templimb = (templimb << 8) + cp[5];
                    templimb = (templimb << 8) + cp[4];
                    templimb = (templimb << 8) + cp[3];
                    templimb = (templimb << 8) + cp[2];
                    templimb = (templimb << 8) + cp[1];
                    templimb = (templimb << 8) + cp[0];
#endif
#if GMP_LIMB_BITS == 32
                    templimb = cp[3];
                    templimb = (templimb << 8) + cp[2];
                    templimb = (templimb << 8) + cp[1];
                    templimb = (templimb << 8) + cp[0];
#endif
                    result->f->_mpfr_d[i] = templimb;
                    cp += newlimbsize;
                }
                result->f->_mpfr_exp = expsgn * exponent;
                if (sgn == -1)
                    mpfr_neg(result->f, result->f, MPFR_RNDN);
                return (PyObject*)result;
            }
            else if (limbsize * sizemant > newmant * newlimbsize) {
                /* Since the amount of saved data is greater than the amount of
                 * data needed on the new system, we skip the first 32 bits
                 * since they must be 0.
                 */

                /* Verify we are on a 32-bit system and the source was 64-bit. */
                if ((limbsize == 8) && (newlimbsize == 4)) {
                    VALUE_ERROR("byte sequence invalid for from_binary()");
                    return NULL;
                }

                mpfr_set_ui(result->f, 1, MPFR_RNDN);
                cp = buffer + 4 + (2 * sizesize) + 4;
                for (i=0; i<newmant; i++) {
                    templimb = cp[3];
                    templimb = (templimb << 8) + cp[2];
                    templimb = (templimb << 8) + cp[1];
                    templimb = (templimb << 8) + cp[0];
                    result->f->_mpfr_d[i] = templimb;
                    cp += newlimbsize;
                }
                result->f->_mpfr_exp = expsgn * exponent;
                if (sgn == -1)
                    mpfr_neg(result->f, result->f, MPFR_RNDN);
                return (PyObject*)result;
            }
            else {
                /* Since the amount of saved data is less than the amount of
                 * data needed on the new system, we must "add" 32 0-bits at
                 * the low end.
                 */

                /* Verify we are on a 64-bit system and the source was 32-bit. */
                if ((limbsize == 4) && (newlimbsize == 8)) {
                    VALUE_ERROR("byte sequence invalid for from_binary()");
                    return NULL;
                }

                mpfr_set_ui(result->f, 1, MPFR_RNDN);
                cp = buffer + 4 + (2 * sizesize);
                templimb = cp[3];
                templimb = (templimb << 8) + cp[2];
                templimb = (templimb << 8) + cp[1];
                templimb = (templimb << 8) + cp[0];
                result->f->_mpfr_d[i] = ((templimb << 16) << 16);
                cp += 4;
                for (i=0; i<newmant-1; i++) {
                    templimb = cp[7];
                    templimb = (templimb << 8) + cp[6];
                    templimb = (templimb << 8) + cp[5];
                    templimb = (templimb << 8) + cp[4];
                    templimb = (templimb << 8) + cp[3];
                    templimb = (templimb << 8) + cp[2];
                    templimb = (templimb << 8) + cp[1];
                    templimb = (templimb << 8) + cp[0];
                    result->f->_mpfr_d[i] = templimb;
                    cp += newlimbsize;
                }
                result->f->_mpfr_exp = expsgn * exponent;
                if (sgn == -1)
                    mpfr_neg(result->f, result->f, MPFR_RNDN);
                return (PyObject*)result;
            }
        }
        case 0x05: {
            MPC_Object *result;
            MPFR_Object *real = 0, *imag = 0;
            Py_ssize_t sizemant = 0, sizesize = 4, i, newmant;
            mpfr_prec_t precision = 0;
            mpfr_exp_t exponent = 0;
            mp_limb_t templimb;
            int sgn = 1, expsgn = 1, limbsize = 4;
            int newlimbsize = (mp_bits_per_limb >> 3);
            unsigned char *tempbuf;

            if (len < 4) {
                VALUE_ERROR("byte sequence too short for from_binary()");
                return NULL;
            }
            /* read the real part first */
            if (cp[1] & 0x04) sizesize = 8;
            for (i=sizesize; i>0; --i) {
                precision = (precision << 8) + cp[i+3];
            }
            if (cp[1] & 0x02) sgn = -1;
            if (cp[1] & 0x20) expsgn = -1;
            if (cp[1] & 0x40) limbsize = 8;
            if (!(real = GMPy_MPFR_New(precision, context))) {
                /* LCOV_EXCL_START */
                return NULL;
                /* LCOV_EXCL_STOP */
            }
            if (cp[2] == 0)      real->rc = 0;
            else if (cp[2] == 1) real->rc = 1;
            else                 real->rc = -1;
            if (!(cp[1] & 0x01)) {
                if ((cp[1] & 0x18) == 0x00)
                    mpfr_set_zero(real->f, sgn);
                else if ((cp[1] & 0x18) == 0x08)
                    mpfr_set_nan(real->f);
                else
                    mpfr_set_inf(real->f, sgn);
                cp += 4 + sizesize;
                goto readimag;
            }
            if (limbsize == 8) sizemant = ((precision + 63) / 64);
            else               sizemant = ((precision + 31) / 32);
            newmant = (precision + mp_bits_per_limb - 1) / mp_bits_per_limb;
            cp = buffer + 4 + sizesize - 1;
            for (i=sizesize; i>0; --i) {
                exponent = (exponent << 8) + cp[i];
            }
            if (limbsize * sizemant == newmant * newlimbsize) {
                mpfr_set_ui(real->f, 1, MPFR_RNDN);
                cp = buffer + 4 + (2 * sizesize);
                for (i=0; i<newmant; i++) {
#if GMP_LIMB_BITS == 64
                    templimb = cp[7];
                    templimb = (templimb << 8) + cp[6];
                    templimb = (templimb << 8) + cp[5];
                    templimb = (templimb << 8) + cp[4];
                    templimb = (templimb << 8) + cp[3];
                    templimb = (templimb << 8) + cp[2];
                    templimb = (templimb << 8) + cp[1];
                    templimb = (templimb << 8) + cp[0];
#endif
#if GMP_LIMB_BITS == 32
                    templimb = cp[3];
                    templimb = (templimb << 8) + cp[2];
                    templimb = (templimb << 8) + cp[1];
                    templimb = (templimb << 8) + cp[0];
#endif
                    real->f->_mpfr_d[i] = templimb;
                    cp += newlimbsize;
                }
                real->f->_mpfr_exp = expsgn * exponent;
                if (sgn == -1)
                    mpfr_neg(real->f, real->f, MPFR_RNDN);
            }
            else if (limbsize * sizemant > newmant * newlimbsize) {
                if ((limbsize == 8) && (newlimbsize == 4)) {
                    VALUE_ERROR("byte sequence invalid for from_binary()");
                    Py_DECREF((PyObject*)real);
                    return NULL;
                }
                mpfr_set_ui(real->f, 1, MPFR_RNDN);
                cp = buffer + 4 + (2 * sizesize) + 4;
                for (i=0; i<newmant; i++) {
                    templimb = cp[3];
                    templimb = (templimb << 8) + cp[2];
                    templimb = (templimb << 8) + cp[1];
                    templimb = (templimb << 8) + cp[0];
                    real->f->_mpfr_d[i] = templimb;
                    cp += newlimbsize;
                }
                real->f->_mpfr_exp = expsgn * exponent;
                if (sgn == -1)
                    mpfr_neg(real->f, real->f, MPFR_RNDN);
            }
            else {
                if ((limbsize == 4) && (newlimbsize == 8)) {
                    VALUE_ERROR("byte sequence invalid for from_binary()");
                    Py_DECREF((PyObject*)real);
                    return NULL;
                }
                mpfr_set_ui(real->f, 1, MPFR_RNDN);
                cp = buffer + 4 + (2 * sizesize);
                templimb = cp[3];
                templimb = (templimb << 8) + cp[2];
                templimb = (templimb << 8) + cp[1];
                templimb = (templimb << 8) + cp[0];
                real->f->_mpfr_d[i] = ((templimb << 16) << 16);
                cp += 4;
                for (i=0; i<newmant-1; i++) {
                    templimb = cp[7];
                    templimb = (templimb << 8) + cp[6];
                    templimb = (templimb << 8) + cp[5];
                    templimb = (templimb << 8) + cp[4];
                    templimb = (templimb << 8) + cp[3];
                    templimb = (templimb << 8) + cp[2];
                    templimb = (templimb << 8) + cp[1];
                    templimb = (templimb << 8) + cp[0];
                    real->f->_mpfr_d[i] = templimb;
                    cp += newlimbsize;
                }
                real->f->_mpfr_exp = expsgn * exponent;
                if (sgn == -1)
                    mpfr_neg(real->f, real->f, MPFR_RNDN);
            }
  readimag:
            /* Set all the variables back to default. */
            tempbuf = cp;

            sizemant = 0;
            sizesize = 4;
            precision = 0;
            exponent = 0;
            sgn = 1;
            expsgn = 1;
            limbsize = 4;

            /* Done reading the real part. The next byte should be 0x05. */
            if (!(cp[0] == 0x05)) {
                VALUE_ERROR("byte sequence invalid for from_binary()");
                Py_DECREF((PyObject*)real);
                return NULL;
            }
            if (cp[1] & 0x04) sizesize = 8;
            for (i=sizesize; i>0; --i) {
                precision = (precision << 8) + cp[i+3];
            }
            if (cp[1] & 0x02) sgn = -1;
            if (cp[1] & 0x20) expsgn = -1;
            if (cp[1] & 0x40) limbsize = 8;
            if (!(imag = GMPy_MPFR_New(precision, context))) {
                /* LCOV_EXCL_START */
                Py_DECREF((PyObject*)real);
                return NULL;
                /* LCOV_EXCL_STOP */
            }
            if (cp[2] == 0)      imag->rc = 0;
            else if (cp[2] == 1) imag->rc = 1;
            else                 imag->rc = -1;
            if (!(cp[1] & 0x01)) {
                if ((cp[1] & 0x18) == 0x00)
                    mpfr_set_zero(imag->f, sgn);
                else if ((cp[1] & 0x18) == 0x08)
                    mpfr_set_nan(imag->f);
                else
                    mpfr_set_inf(imag->f, sgn);
                goto alldone;
            }
            if (limbsize == 8) sizemant = ((precision + 63) / 64);
            else               sizemant = ((precision + 31) / 32);
            newmant = (precision + mp_bits_per_limb - 1) / mp_bits_per_limb;
            cp = tempbuf + 4 + sizesize - 1;
            for (i=sizesize; i>0; --i) {
                exponent = (exponent << 8) + cp[i];
            }
            if (limbsize * sizemant == newmant * newlimbsize) {
                mpfr_set_ui(imag->f, 1, MPFR_RNDN);
                cp = tempbuf + 4 + (2 * sizesize);
                for (i=0; i<newmant; i++) {
#if GMP_LIMB_BITS == 64
                    templimb = cp[7];
                    templimb = (templimb << 8) + cp[6];
                    templimb = (templimb << 8) + cp[5];
                    templimb = (templimb << 8) + cp[4];
                    templimb = (templimb << 8) + cp[3];
                    templimb = (templimb << 8) + cp[2];
                    templimb = (templimb << 8) + cp[1];
                    templimb = (templimb << 8) + cp[0];
#endif
#if GMP_LIMB_BITS == 32
                    templimb = cp[3];
                    templimb = (templimb << 8) + cp[2];
                    templimb = (templimb << 8) + cp[1];
                    templimb = (templimb << 8) + cp[0];
#endif
                    imag->f->_mpfr_d[i] = templimb;
                    cp += newlimbsize;
                }
                imag->f->_mpfr_exp = expsgn * exponent;
                if (sgn == -1)
                    mpfr_neg(imag->f, imag->f, MPFR_RNDN);
            }
            else if (limbsize * sizemant > newmant * newlimbsize) {
                if ((limbsize == 8) && (newlimbsize == 4)) {
                    VALUE_ERROR("byte sequence invalid for from_binary()");
                    Py_DECREF((PyObject*)real);
                    Py_DECREF((PyObject*)imag);
                    return NULL;
                }
                mpfr_set_ui(imag->f, 1, MPFR_RNDN);
                cp = tempbuf + 4 + (2 * sizesize) + 4;
                for (i=0; i<newmant; i++) {
                    templimb = cp[3];
                    templimb = (templimb << 8) + cp[2];
                    templimb = (templimb << 8) + cp[1];
                    templimb = (templimb << 8) + cp[0];
                    imag->f->_mpfr_d[i] = templimb;
                    cp += newlimbsize;
                }
                imag->f->_mpfr_exp = expsgn * exponent;
                if (sgn == -1)
                    mpfr_neg(imag->f, imag->f, MPFR_RNDN);
            }
            else {
                if ((limbsize == 4) && (newlimbsize == 8)) {
                    VALUE_ERROR("byte sequence invalid for from_binary()");
                    Py_DECREF((PyObject*)real);
                    Py_DECREF((PyObject*)imag);
                    return NULL;
                }
                mpfr_set_ui(imag->f, 1, MPFR_RNDN);
                cp = tempbuf + 4 + (2 * sizesize);
                templimb = cp[3];
                templimb = (templimb << 8) + cp[2];
                templimb = (templimb << 8) + cp[1];
                templimb = (templimb << 8) + cp[0];
                imag->f->_mpfr_d[i] = ((templimb << 16) << 16);
                cp += 4;
                for (i=0; i<newmant-1; i++) {
                    templimb = cp[7];
                    templimb = (templimb << 8) + cp[6];
                    templimb = (templimb << 8) + cp[5];
                    templimb = (templimb << 8) + cp[4];
                    templimb = (templimb << 8) + cp[3];
                    templimb = (templimb << 8) + cp[2];
                    templimb = (templimb << 8) + cp[1];
                    templimb = (templimb << 8) + cp[0];
                    imag->f->_mpfr_d[i] = templimb;
                    cp += newlimbsize;
                }
                imag->f->_mpfr_exp = expsgn * exponent;
                if (sgn == -1)
                    mpfr_neg(imag->f, imag->f, MPFR_RNDN);
            }
  alldone:
            if (!(result = (MPC_Object*)GMPy_MPC_New(0, 0, context))) {
                /* LCOV_EXCL_START */
                Py_DECREF((PyObject*)real);
                Py_DECREF((PyObject*)imag);
                return NULL;
                /* LCOV_EXCL_STOP */
            }
            mpfr_swap(mpc_realref(result->c), real->f);
            mpfr_swap(mpc_imagref(result->c), imag->f);
            Py_DECREF((PyObject*)real);
            Py_DECREF((PyObject*)imag);
            return (PyObject*)result;
        }
        default: {
            TYPE_ERROR("from_binary() argument type not supported");
            return NULL;
        }
    }
}

PyDoc_STRVAR(doc_to_binary,
"to_binary($module, x, /)\n--\n\n"
"Return a Python byte sequence that is a portable binary\n"
"representation of a gmpy2 object x. The byte sequence can\n"
"be passed to `from_binary()` to obtain an exact copy of\n"
"x's value.  Raises a `TypeError` if x is not a gmpy2 object.");

static PyObject *
GMPy_MPANY_To_Binary(PyObject *self, PyObject *other)
{
    if(MPZ_Check(other))
        return GMPy_MPZ_To_Binary((MPZ_Object*)other);
    else if(XMPZ_Check(other))
        return GMPy_XMPZ_To_Binary((XMPZ_Object*)other);
    else if(MPQ_Check(other))
        return GMPy_MPQ_To_Binary((MPQ_Object*)other);
    else if(MPFR_Check(other))
        return GMPy_MPFR_To_Binary((MPFR_Object*)other);
    else if(MPC_Check(other))
        return GMPy_MPC_To_Binary((MPC_Object*)other);
    TYPE_ERROR("to_binary() argument type not supported");
    return NULL;
}
