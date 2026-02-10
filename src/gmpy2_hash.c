/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_hash.c                                                            *
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


#include "pythoncapi_compat.h"

static Py_hash_t
GMPy_MPZ_Hash_Slot(MPZ_Object *self)
{
    Py_hash_t hash;

    if (self->hash_cache != -1) {
        return self->hash_cache;
    }

    hash = (Py_hash_t)mpn_mod_1(self->z->_mp_d, (mp_size_t)mpz_size(self->z), PyHASH_MODULUS);
    if (mpz_sgn(self->z) < 0) {
        hash = -hash;
    }
    if (hash == -1) {
        hash = -2;
    }
    return (self->hash_cache = hash);
}

static Py_hash_t
GMPy_MPQ_Hash_Slot(MPQ_Object *self)
{
    Py_hash_t hash = 0;
    mpz_t temp, temp1, mask;

    if (self->hash_cache != -1) {
        return self->hash_cache;
    }

    mpz_init(temp);
    mpz_init(temp1);
    mpz_init(mask);
    mpz_set_si(mask, 1);
    mpz_mul_2exp(mask, mask, PyHASH_BITS);
    mpz_sub_ui(mask, mask, 1);

    if (!mpz_invert(temp, mpq_denref(self->q), mask)) {
        mpz_clear(temp);
        mpz_clear(temp1);
        mpz_clear(mask);
        hash = PyHASH_INF;
        if (mpz_sgn(mpq_numref(self->q)) < 0) {
            hash = -hash;
        }
        self->hash_cache = hash;
        return hash;
    }
    mpz_set(temp1, mask);
    mpz_sub_ui(temp1, temp1, 2);
    mpz_powm(temp, mpq_denref(self->q), temp1, mask);

    mpz_tdiv_r(temp1, mpq_numref(self->q), mask);
    mpz_mul(temp, temp, temp1);
    hash = (Py_hash_t)mpn_mod_1(temp->_mp_d, (mp_size_t)mpz_size(temp), PyHASH_MODULUS);

    if (mpz_sgn(mpq_numref(self->q)) < 0) {
        hash = -hash;
    }
    if (hash == -1) {
        hash = -2;
    }
    mpz_clear(temp);
    mpz_clear(temp1);
    mpz_clear(mask);
    self->hash_cache = hash;
    return hash;
}

static Py_hash_t
_mpfr_hash(mpfr_t f)
{
    Py_uhash_t hash = 0;
    Py_ssize_t exp;
    size_t msize;
    int sign;

    /* Handle special cases first */
    if (!mpfr_number_p(f)) {
        if (mpfr_inf_p(f)) {
            if (mpfr_sgn(f) > 0) {
                return PyHASH_INF;
            }
            else {
                return -PyHASH_INF;
            }
        }
        else {
#if PY_VERSION_HEX >= 0x030A00A0
            return Py_HashPointer(f);
#else
            return _PyHASH_NAN;
#endif
        }
    }

    /* Calculate the number of limbs in the mantissa. */
    msize = (f->_mpfr_prec + mp_bits_per_limb - 1) / mp_bits_per_limb;

    /* Calculate the hash of the mantissa. */
    if (mpfr_sgn(f) > 0) {
        hash = mpn_mod_1(f->_mpfr_d, (mp_size_t)msize, PyHASH_MODULUS);
        sign = 1;
    }
    else if (mpfr_sgn(f) < 0) {
        hash = mpn_mod_1(f->_mpfr_d, (mp_size_t)msize, PyHASH_MODULUS);
        sign = -1;
    }
    else {
        return 0;
    }

    /* Calculate the final hash. */
    exp = f->_mpfr_exp - (msize * mp_bits_per_limb);
    exp = exp >= 0 ? exp % PyHASH_BITS : PyHASH_BITS-1-((-1-exp) % PyHASH_BITS);
    hash = ((hash << exp) & PyHASH_MODULUS) | hash >> (PyHASH_BITS - exp);

    hash *= sign;
    if (hash == (Py_uhash_t)(-1)) {
        hash = (Py_uhash_t)(-2);
    }
    return (Py_hash_t)hash;
}

static Py_hash_t
GMPy_MPFR_Hash_Slot(MPFR_Object *self)
{
    if (self->hash_cache == -1) {
        self->hash_cache = _mpfr_hash(self->f);
    }
    return self->hash_cache;
}

static Py_hash_t
GMPy_MPC_Hash_Slot(MPC_Object *self)
{
    Py_uhash_t hashreal, hashimag, combined;

    if (self->hash_cache != -1) {
        return self->hash_cache;
    }

    hashreal = (Py_uhash_t)_mpfr_hash(mpc_realref(self->c));
    hashimag = (Py_uhash_t)_mpfr_hash(mpc_imagref(self->c));
    combined = hashreal + PyHASH_IMAG * hashimag;
    if (combined == (Py_uhash_t)(-1)) {
        combined = (Py_uhash_t)(-2);
    }
    self->hash_cache = combined;
    return (Py_hash_t)combined;
}
