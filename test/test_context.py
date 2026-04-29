import warnings

import pytest

import gmpy2
from gmpy2 import (context, get_context, ieee, local_context, mpc, mpfr, mpz,
                   set_context)


def test_context_abs():
    ctx = gmpy2.context()

    assert ctx.abs(-1) == mpz(1)
    assert ctx.abs(0) == mpz(0)
    assert ctx.abs(1) == mpz(1)
    assert ctx.abs(mpz(8)) == mpz(8)
    assert ctx.abs(mpz(-8)) == mpz(8)
    assert ctx.abs(-1.0) == mpfr('1.0')
    assert ctx.abs(mpfr(-2)) == mpfr('2.0')
    assert ctx.abs(2+3j) == mpfr('3.6055512754639891')
    assert ctx.abs(mpc(2+3j)) == mpfr('3.6055512754639891')


def test_context_ieee():
    ctx = ieee(32)

    assert (ctx.precision == 24 and ctx.emax == 128 and
            ctx.emin == -148 and ctx.subnormalize)

    ctx = ieee(64)

    assert (ctx.precision == 53 and ctx.emax == 1024 and
            ctx.emin == -1073 and ctx.subnormalize)

    ctx = ieee(128)

    assert (ctx.precision == 113 and ctx.emax == 16384 and
            ctx.emin == -16493 and ctx.subnormalize)

    ctx = ieee(256)

    assert (ctx.precision == 237 and ctx.emax == 262144 and
            ctx.emin == -262377 and ctx.subnormalize)

    pytest.raises(ValueError, lambda: ieee(-1))
    pytest.raises(TypeError, lambda: ieee("a"))

    set_context(ieee(32))

    assert gmpy2.const_pi().digits(2) == ('110010010000111111011011', 2, 24)

    set_context(ieee(64))

    assert gmpy2.const_pi().digits(2) == ('11001001000011111101101010100010001000010110100011000', 2, 53)

    set_context(ieee(128))

    assert gmpy2.const_pi().digits(2) == ('11001001000011111101101010100010001000010110100011000010001101001100010011000110011000101000101110000000110111000', 2, 113)

    pytest.raises(ValueError, lambda: ieee(16*31))
    assert ieee(32*31).precision == 965


def test_context():
    ctx = context()

    assert (ctx.precision == 53 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    ctx = context(precision=100)

    assert (ctx.precision == 100 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    ctx = context(real_prec=100)

    assert (ctx.precision == 53 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize and
            ctx.real_prec == 100)

    ctx = context(real_prec=100,imag_prec=200)

    assert (ctx.precision == 53 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize and
            ctx.real_prec == 100 and ctx.imag_prec == 200)

    pytest.raises(ValueError, lambda: context(1, 2))
    pytest.raises(ValueError, lambda: context(spam=123))


def test_get_context():
    set_context(context())

    ctx = get_context()

    assert (ctx.precision == 53 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    ctx = get_context()
    ctx.precision = 100

    assert (ctx.precision == 100 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    ctx = get_context()

    assert (ctx.precision == 100 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    b = ctx.copy()
    b.precision = 200

    assert (b.precision == 200 and b.emax == 1073741823 and
            b.emin == -1073741823 and not b.subnormalize)
    assert (ctx.precision == 100 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    ctx = get_context()

    assert (ctx.precision == 100 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)


def test_context_2():
    set_context(context())

    with context() as ctx:
        assert ctx.precision == 53
        ctx.precision += 20
        assert ctx.precision == 73

    ctx = get_context()

    assert (ctx.precision == 53 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    with context(ieee(64)) as ctx:
        assert (ctx.precision == 53 and ctx.emax == 1024 and
                ctx.emin == -1073 and ctx.subnormalize)

    ctx = get_context()

    assert (ctx.precision == 53 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    with get_context() as ctx:
        assert ctx.precision == 53
        ctx.precision += 100
        assert ctx.precision == 153

    ctx = get_context()

    assert (ctx.precision == 53 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)

    with context(precision=200) as ctx:
        assert ctx.precision == 200
        ctx.precision += 100
        assert ctx.precision == 300

    ctx = get_context()

    assert (ctx.precision == 53 and ctx.emax == 1073741823 and
            ctx.emin == -1073741823 and not ctx.subnormalize)


def test_nested_context():
    set_context(context())

    r = [get_context().precision]

    with ieee(128):
        r.append(get_context().precision)
        with ieee(256):
            r.append(get_context().precision)
            with ieee(512):
                r.append(get_context().precision)
            r.append(get_context().precision)
        r.append(get_context().precision)
    r.append(get_context().precision)
    assert r == [53, 113, 237, 489, 237, 113, 53]


@pytest.mark.filterwarnings("ignore:local_context().*:DeprecationWarning")
def test_nested_local_context():
    set_context(context())

    r = [get_context().precision]

    with local_context(ieee(128)):
        r.append(get_context().precision)
        with local_context(ieee(256)):
            r.append(get_context().precision)
            with local_context(ieee(512)):
                r.append(get_context().precision)
            r.append(get_context().precision)
        r.append(get_context().precision)
    r.append(get_context().precision)
    assert r == [53, 113, 237, 489, 237, 113, 53]


def test_context_repr():
    ctx = get_context()
    assert repr(ctx) == \
"""context(precision=53, real_prec=Default, imag_prec=Default,\n\
        round=RoundToNearest, real_round=Default, imag_round=Default,\n\
        emax=1073741823, emin=-1073741823,\n        subnormalize=False,\n\
        trap_underflow=False, underflow=False,\n        trap_overflow=False,\
 overflow=False,\n        trap_inexact=False, inexact=False,\n\
        trap_invalid=False, invalid=False,\n        trap_erange=False,\
 erange=False,\n        trap_divzero=False, divzero=False,\n\
        allow_complex=False,\n        rational_division=False,\n\
        allow_release_gil=False)"""

    ctx.real_prec = 100
    ctx.imag_prec = 200
    assert repr(ctx) == \
"""context(precision=53, real_prec=100, imag_prec=200,\n\
        round=RoundToNearest, real_round=Default, imag_round=Default,\n\
        emax=1073741823, emin=-1073741823,\n        subnormalize=False,\n\
        trap_underflow=False, underflow=False,\n        trap_overflow=False,\
 overflow=False,\n        trap_inexact=False, inexact=False,\n\
        trap_invalid=False, invalid=False,\n        trap_erange=False,\
 erange=False,\n        trap_divzero=False, divzero=False,\n\
        allow_complex=False,\n        rational_division=False,\n\
        allow_release_gil=False)"""
    ctx.trap_invalid = True
    assert repr(ctx) == \
"""context(precision=53, real_prec=100, imag_prec=200,\n\
        round=RoundToNearest, real_round=Default, imag_round=Default,\n\
        emax=1073741823, emin=-1073741823,\n        subnormalize=False,\n\
        trap_underflow=False, underflow=False,\n        trap_overflow=False,\
 overflow=False,\n        trap_inexact=False, inexact=False,\n\
        trap_invalid=True, invalid=False,\n        trap_erange=False,\
 erange=False,\n        trap_divzero=False, divzero=False,\n\
        allow_complex=False,\n        rational_division=False,\n\
        allow_release_gil=False)"""
    pytest.raises(gmpy2.InvalidOperationError, lambda: mpfr('nan') % 123)
    assert repr(ctx) == \
"""context(precision=53, real_prec=100, imag_prec=200,\n\
        round=RoundToNearest, real_round=Default, imag_round=Default,\n\
        emax=1073741823, emin=-1073741823,\n        subnormalize=False,\n\
        trap_underflow=False, underflow=False,\n        trap_overflow=False,\
 overflow=False,\n        trap_inexact=False, inexact=False,\n\
        trap_invalid=True, invalid=True,\n        trap_erange=False,\
 erange=False,\n        trap_divzero=False, divzero=False,\n\
        allow_complex=False,\n        rational_division=False,\n\
        allow_release_gil=False)"""
    set_context(ieee(32))
    ctx = get_context()
    ctx.trap_underflow = True
    c = mpc(0.1 + 0.1j)
    pytest.raises(gmpy2.UnderflowResultError, lambda: c**201)
    assert repr(ctx) == \
"""context(precision=24, real_prec=Default, imag_prec=Default,\n\
        round=RoundToNearest, real_round=Default, imag_round=Default,\n\
        emax=128, emin=-148,\n        subnormalize=True,\n\
        trap_underflow=True, underflow=True,\n        trap_overflow=False,\
 overflow=False,\n        trap_inexact=False, inexact=True,\n\
        trap_invalid=False, invalid=False,\n        trap_erange=False,\
 erange=False,\n        trap_divzero=False, divzero=False,\n\
        allow_complex=False,\n        rational_division=False,\n\
        allow_release_gil=False)"""


def test_local_context_deprecated():
    with pytest.deprecated_call():
        local_context()

    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        pytest.raises(DeprecationWarning, lambda: local_context())


@pytest.mark.filterwarnings("ignore:local_context().*:DeprecationWarning")
def test_local_context():
    ctx_orig = get_context()
    ctx_orig.precision = 123
    with context() as ctx:
        assert ctx.precision == 53
    with local_context() as ctx:
        assert ctx.precision == 123
        ctx.precision = 321
    with local_context(ctx_orig) as ctx:
        assert ctx.precision == 123

    pytest.raises(ValueError, lambda: local_context(1, 2))
    pytest.raises(ValueError, lambda: local_context(spam=123))
