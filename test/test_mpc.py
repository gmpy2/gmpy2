from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import cmath
import inspect
import sys

import pytest
from hypothesis import example, given, settings
from hypothesis.strategies import complex_numbers, floats
from supportclasses import a, b, c, cx, d, q, r, z

import gmpy2
from gmpy2 import (cmp, cmp_abs, from_binary, get_context, is_nan, mpc,
                   mpc_random, mpfr, mpq, mpz, nan, random_state, to_binary)


def test_mpc_cmp():
    pytest.raises(TypeError, lambda: cmp(mpc(1,2), mpc(3,4)))
    assert cmp_abs(mpc(1,2), mpc(3,4)) == -1
    assert cmp_abs(mpc(1,2), mpc(1,2)) == 0
    assert cmp_abs(mpc(3,4), mpc(1,2)) == 1
    gmpy2.get_context().clear_flags()
    assert gmpy2.get_context().erange is False
    assert cmp_abs(mpc(nan(),1), mpc(4.5)) == 0
    assert gmpy2.get_context().erange is True


def test_mpc_comparisons():
    a = mpz(123)
    c = mpc(4.5,6)

    assert (c == a, c != a) == (False, True)

    pytest.raises(TypeError, lambda: c < a)

    cnan = mpc(mpfr('nan'), 6)

    assert (c == cnan, c != cnan) == (False, True)

    pytest.raises(TypeError, lambda: c > cnan)

    cx = complex(4.5,6)

    assert (c == cx, c != cx) == (True, False)

    pytest.raises(TypeError, lambda: c > cx)

    from supportclasses import c

    assert (c == cx) is False


def test_mpc_conversion():
    x = mpc(a)
    assert isinstance(x, mpc)
    assert x == 42+67j

    pytest.raises(TypeError, lambda: mpc(b))
    pytest.raises(TypeError, lambda: mpc(c))
    pytest.raises(TypeError, lambda: mpc(d))

    assert mpc(mpfr(5.6), precision=(0,0)) == mpc('5.5999999999999996+0.0j')
    assert mpc(Fraction(4,5)) == mpc('0.80000000000000004+0.0j')
    assert mpc(b'5+6j') == mpc('5.0+6.0j')
    assert mpc('5+6j') == mpc('5.0+6.0j')

    pytest.raises(ValueError, lambda: mpc('\xc3'))
    pytest.raises(ValueError, lambda: mpc('notanumber'))
    pytest.raises(ValueError, lambda: mpc('\0'))

    assert mpc('(5+6j)') == mpc('5.0+6.0j')
    assert mpc('   5+6j   ') == mpc('5.0+6.0j')

    pytest.raises(ValueError, lambda: mpc('5+6ji'))

    assert str(mpc(5,6)) == '5.0+6.0j'
    assert complex(mpc(4,5)) == (4+5j)
    assert complex(mpc(4,5)) == (4+5j)


def test_mpc_creation():
    ctx = gmpy2.get_context()
    ctx.clear_flags()
    a = mpc("1.2")
    assert a.rc == (-1, 0)
    assert ctx.inexact
    ctx.clear_flags()
    a = mpc("(1 2)")
    assert a == 1 + 2j
    assert a.rc == (0, 0)
    assert ctx.inexact is False
    ctx.clear_flags()
    a = mpc("1   + 2.1  j")
    assert a == 1 + 2.1j
    assert a.rc == (0, 1)
    assert ctx.inexact

    c = mpc('67+87j', precision=70)

    assert c.precision == (70, 70)
    assert c.real.precision == 70
    assert c.imag.precision == 70

    c = mpc('42e56+42.909j', precision=(45,300))

    assert c.precision == (45, 300)
    assert c.real.precision == 45
    assert c.imag.precision == 300

    x = mpc("1.3142123+4.3e-1001j", precision=(70,37))

    assert mpc(x.real, x.imag, precision=(70,37)) == x

    assert mpc('1_2+4_5j') == mpc('12.0+45.0j')

    pytest.raises(TypeError, lambda: mpc(1, base=2))
    pytest.raises(TypeError, lambda: mpc(1, s=2))
    pytest.raises(TypeError, lambda: mpc("1", s=2))
    pytest.raises(TypeError, lambda: mpc("1", imag=2))
    pytest.raises(TypeError, lambda: mpc(1j, imag=2))
    pytest.raises(TypeError, lambda: mpc(1j, base=2))
    pytest.raises(TypeError, lambda: mpc(1j, s=2))


def test_mpc_random():
    assert (mpc_random(random_state(42))
            == mpc('0.86555158787663011+0.4422082613292212j'))


def test_mpc_to_from_binary():
    x = mpc("0+0j")
    assert x == from_binary(to_binary(x))
    x = mpc("1+0j")
    assert x == from_binary(to_binary(x))
    x = mpc("-1+0j")
    assert x == from_binary(to_binary(x))
    x = mpc("0+1j")
    assert x == from_binary(to_binary(x))
    x = mpc("0-1j")
    assert x == from_binary(to_binary(x))
    x = mpc("inf+0j")
    assert x == from_binary(to_binary(x))
    x = mpc("0+infj")
    assert x == from_binary(to_binary(x))
    x = mpc("inf-infj")
    assert x == from_binary(to_binary(x))
    x = mpc("inf+nanj")
    y = from_binary(to_binary(x))
    assert x.real == y.real
    assert is_nan(y.imag)
    x = mpc("-inf+0j")
    assert x == from_binary(to_binary(x))
    x = mpc("0-infj")
    assert x == from_binary(to_binary(x))
    x = mpc("-inf-infj")
    assert x == from_binary(to_binary(x))
    x = mpc("-inf+nanj")
    y = from_binary(to_binary(x))
    assert x.real == y.real
    assert is_nan(y.imag)
    x = mpc("nan+0j")
    y = from_binary(to_binary(x))
    assert x.imag == y.imag
    assert is_nan(y.real)
    x = mpc("0+nanj")
    y = from_binary(to_binary(x))
    assert x.real == y.real
    assert is_nan(y.imag)
    x = mpc("nan-infj")
    y = from_binary(to_binary(x))
    assert x.imag == y.imag
    assert is_nan(y.real)
    x = mpc("nan+nanj")
    y = from_binary(to_binary(x))
    assert is_nan(y.real)
    assert is_nan(y.imag)
    get_context().real_prec=100
    get_context().imag_prec=110
    assert (from_binary(to_binary(mpc("1.3-4.7j"))) ==
            mpc('1.2999999999999999999999999999994-4.7000000000000000000000000000000025j',
                (100,110)))


def test_mpc_format():
    gmpy2.set_context(gmpy2.context())

    c, c1 = mpc(mpq(1/3), 5), mpc(-1, -2)

    assert '{:>30}'.format(c) == '      0.33333333333333331+5.0j'
    assert '{:<30}'.format(c) == '0.33333333333333331+5.0j      '
    assert '{:^30}'.format(c) == '   0.33333333333333331+5.0j   '

    pytest.raises(ValueError, lambda: '{:<<20}'.format(c))

    assert '{:>+30}'.format(c) == '     +0.33333333333333331+5.0j'

    pytest.raises(ValueError, lambda: '{:+^20}'.format(c))

    assert '{:.10f}'.format(c) == '0.3333333333+5.0000000000j'

    pytest.raises(ValueError, lambda: '{:Z.10f}'.format(c))
    pytest.raises(ValueError, lambda: '{:Z 10}'.format(c))

    assert '{:Z}'.format(c) == '0.33333333333333331+5.0j'
    assert '{:U}'.format(c) == '0.33333333333333332+5.0j'

    pytest.raises(ValueError, lambda: '{:PU}'.format(c))

    assert '{:UP}'.format(c) == '0.333334+5.0j'

    pytest.raises(ValueError, lambda: '{:PP}'.format(c))

    assert '{:G}'.format(c) == '0.333333+5.0j'
    assert '{:M}'.format(c) == '(0.333333 5.0)'
    assert '{:b}'.format(c) == '1.0101010101010101010101010101010101010101010101010101p-2+1.01p+2j'
    assert '{:a}'.format(c) == '0x5.5555555555554p-4+0x5p+0j'
    assert '{:e}'.format(c) == '3.333333e-01+5.000000e+00j'
    assert '{:M}'.format(c1) == '(-1.0 -2.0)'

    assert "{:#g}".format(mpc(1)) == '1.00000+0.00000j'

    # issue 666
    r = mpc('1.5707963267948966')
    assert f'{r:e}' == '1.570796e+00+0.000000e+00j'

    # issue 503
    c = mpc(2.675)
    assert f'{c:.2f}' == '2.67+0.00j'
    gmpy2.set_context(gmpy2.context(round=gmpy2.RoundUp))
    assert f'{c:.2f}' == '2.68+0.00j'


def test_mpc_repr():
    c = mpc('1.2999999999999999999999999999994-4.7000000000000000000000000000000025j',(100,110))
    assert repr(c) == "mpc('1.2999999999999999999999999999994-4.7000000000000000000000000000000025j',(100,110))"
    assert repr(mpc("1+nanj")) == "mpc('1.0+nanj')"
    assert repr(mpc("infj")) == "mpc('0.0+infj')"


def test_mpc_digits():
    c = mpc(mpq(1/3), 5)

    assert c.digits(8) == (('2525252525252525250', 0, 53), ('5000000000000000000', 1, 53))
    assert c.digits(8, 2) == (('25', 0, 53), ('50', 1, 53))

    pytest.raises(ValueError, lambda: c.digits(8, -2))
    pytest.raises(ValueError, lambda: c.digits(0))


def test_mpc_abs():
    assert abs(mpc(-1,0)) == mpfr('1.0')
    assert abs(-1+0j) == 1.0
    assert abs(mpc(1,1)) == mpfr('1.4142135623730951')

    ctx = gmpy2.get_context()
    ctx.clear_flags()

    c = mpc('nan+0.0j')

    assert is_nan(c.real) and c.imag == 0.0
    assert ctx.invalid

    ctx.clear_flags()

    c = mpc('nan+0j')

    assert is_nan(c.real) and c.imag == 0.0
    assert ctx.invalid

    ctx.clear_flags()

    assert is_nan(abs(mpc('nanj'))) and ctx.invalid

    ctx.clear_flags()

    assert abs(mpc('inf+10j')) == mpfr('inf')
    assert abs(mpc('-infj')) == mpfr('inf')

    a = mpc('nan+infj')

    assert abs(a) and not ctx.invalid

    a = mpc('-inf+nanj')

    assert abs(a) and not ctx.invalid


def test_mpc_add():
    a = mpz(123)
    aj = mpc(1+2j)
    bj = mpc(4+5j)

    assert aj + bj == mpc('5.0+7.0j')
    assert bj + aj == mpc('5.0+7.0j')
    assert aj + a == mpc('124.0+2.0j')
    assert a + aj == mpc('124.0+2.0j')
    assert aj + 1 == mpc('2.0+2.0j')
    assert 1 + aj == mpc('2.0+2.0j')
    assert aj + 0 == mpc('1.0+2.0j')
    assert 0 + aj == mpc('1.0+2.0j')
    assert -1 + aj == mpc('0.0+2.0j')
    assert aj - 1 == mpc('0.0+2.0j')
    assert aj + 1.2 == 2.2 + 2j

    assert aj + float('inf') == mpc('inf+2.0j')
    assert aj + float('-inf') == mpc('-inf+2.0j')

    x = aj + float('nan')

    assert is_nan(x.real) and x.imag == 2.0

    a = mpc('1-0j')
    b = mpfr('1')
    c = mpc('2-0j')

    assert repr(a+b) == repr(b + a) == repr(c)


def test_mpc_sub():
    pytest.raises(TypeError, lambda: mpc(1,2) - 'a')

    assert mpfr(1) - mpc(1,2) == mpc('0.0-2.0j')
    assert mpc(1,2) - mpfr(1) == mpc('0.0+2.0j')
    assert mpc(1,2) - 1+0j == mpc('0.0+2.0j')
    assert 1+0j - mpc(1,2) == mpc('0.0-2.0j')
    assert mpc(1,2) - z == mpc('-1.0+2.0j')
    assert mpc(1,2) - q == mpc('-0.5+2.0j')
    assert mpc(1,2) - r == mpc('-0.5+2.0j')
    assert mpc(1,2) - cx == mpc('-41.0-65.0j')

    a = mpz(123)
    aj = mpc(1+2j)
    bj = mpc(4+5j)

    assert aj - bj == mpc('-3.0-3.0j')
    assert bj - aj == mpc('3.0+3.0j')
    assert aj - a == mpc('-122.0+2.0j')
    assert a - aj == mpc('122.0-2.0j')
    assert aj - 1 == mpc('0.0+2.0j')
    assert 1 - aj == mpc('0.0-2.0j')
    assert 0 - aj == mpc('-1.0-2.0j')
    assert aj - 0 == mpc('1.0+2.0j')
    assert aj - -1 == mpc('2.0+2.0j')
    assert -1 - aj == mpc('-2.0-2.0j')
    assert aj - 1.2 == (1+2j) - 1.2

    assert float('inf') - aj == mpc('inf-2.0j')
    assert aj - float('inf') == mpc('-inf+2.0j')
    assert aj - float('-inf') == mpc('inf+2.0j')

    x = aj - float('nan')

    assert is_nan(x.real) and x.imag == 2.0

    a = mpc('1+0j')
    b = mpfr('1')
    c = mpc('0-0j')

    assert repr(b - a) == repr(c)


def test_mpc_mul():
    pytest.raises(TypeError, lambda: mpc(1,2) * 'a')

    assert mpfr(1) * mpc(1,2) == mpc('1.0+2.0j')
    assert mpc(1,2) * mpfr(1) == mpc('1.0+2.0j')
    assert mpc(1,2) * mpfr(-1) == mpc('-1.0-2.0j')
    assert mpc(1,2) * (1+0j) == mpc('1.0+2.0j')
    assert (1+0j) * mpc(1,2) == mpc('1.0+2.0j')
    assert mpc(1,2) * z == mpc('2.0+4.0j')
    assert mpc(1,2) * q == mpc('1.5+3.0j')
    assert mpc(1,2) * r == mpc('1.5+3.0j')
    assert mpc(1,2) * cx == mpc('-92.0+151.0j')

    a = mpz(123)
    aj = mpc(1+2j)
    bj = mpc(4+5j)

    assert aj * bj == mpc('-6.0+13.0j')
    assert bj * aj == mpc('-6.0+13.0j')
    assert aj * a == mpc('123.0+246.0j')
    assert a * aj == mpc('123.0+246.0j')
    assert aj * -1 == mpc('-1.0-2.0j')
    assert aj * (0.0+1j) == mpc('-2.0+1.0j')

    assert aj * float('inf') == mpc('inf+infj')
    assert aj * float('-inf') == mpc('-inf-infj')

    for x in [aj * float('nan'), mpc(0,0) * float('inf'),
              mpc(0,0) * float('-inf'), mpc(0,0) * float('nan')]:
        assert all(is_nan(_) for _ in [x.real, x.imag])

    a = mpc('inf-1j')
    b = mpfr('inf')
    c = mpc('inf-infj')

    assert a*b == b*a == c


def test_mpc_divmod():
    pytest.raises(TypeError, lambda: divmod(mpc(1),'a'))

    ctx = gmpy2.context()

    pytest.raises(TypeError, lambda: ctx.divmod(mpc(1,2),mpc(3,4)))
    pytest.raises(TypeError, lambda: divmod(mpc(1,2), mpc(1,2)))
    pytest.raises(TypeError, lambda: ctx.divmod(mpc(1,2),mpc(3,4)))

    aj = mpc(1+2j)
    bj = mpc(4+5j)

    with pytest.raises(TypeError):
        divmod(aj, bj)


def test_mpc_div():
    a = mpz(123)
    aj = mpc(1+2j)
    bj = mpc(4+5j)

    assert aj / bj == mpc('0.34146341463414637+0.073170731707317069j')
    assert gmpy2.div(aj, bj) == mpc('0.34146341463414637+0.073170731707317069j')

    with pytest.raises(TypeError):
        aj // bj

    assert aj / a == mpc('0.008130081300813009+0.016260162601626018j')
    assert a / aj == mpc('24.600000000000001-49.200000000000003j')
    assert aj / 0 == mpc('inf+infj')
    assert mpc('2.0+2.0j') / z == mpc('1.0+1.0j')
    assert mpc('2.0+2.0j') / q == mpc('1.3333333333333333+1.3333333333333333j')
    assert mpc('2.0+2.0j') / r == mpc('1.3333333333333333+1.3333333333333333j')
    assert mpc(15,15) / cx == mpc('0.26147449224372299-0.059971213817367662j')

    with gmpy2.context(trap_divzero=True):
        with pytest.raises(gmpy2.DivisionByZeroError):
            mpc(15, 15)/mpc(0)
        with pytest.raises(gmpy2.DivisionByZeroError):
            mpfr(1)/mpc(0)
        with pytest.raises(gmpy2.DivisionByZeroError):
            mpc(1)/mpfr(0)

    assert aj / float('inf') == mpc('0.0+0.0j')
    assert aj / float('-inf') == mpc('-0.0-0.0j')
    assert float('inf') / aj == mpc('inf-infj')
    assert float('-inf') / aj == mpc('-inf+infj')

    a = mpc('1-infj')
    b = 1

    assert a/b == a
    assert repr(b/a) == repr(mpc(0))


def test_mpc_mod():
    a = mpz(123)
    aj = mpc(1+2j)
    bj = mpc(4+5j)

    with pytest.raises(TypeError):
        aj % bj
    with pytest.raises(TypeError):
        aj % a
    with pytest.raises(TypeError):
        a % aj


def test_mpc_pow():
    c1, c2 = mpc(2,5), mpc(5,2)
    ctx = gmpy2.get_context()

    assert ctx.pow(complex(2,5), complex(5,2)) == mpc('-416.55882051164394+44.334999625388825j')
    assert pow(c1, c2) == mpc('-416.55882051164394+44.334999625388825j')
    assert ctx.pow(c1, c2) == mpc('-416.55882051164394+44.334999625388825j')
    assert ctx.pow(c1, c2) == c1 ** c2

    pytest.raises(TypeError, lambda: pow(c1, c2, 5))

    assert pow(c1, 5) == mpc('4282.0-1475.0j')
    assert c1 ** mpz(5) == mpc('4282.0-1475.0j')
    assert c1 ** mpfr(2.5) == mpc('-66.373652915897722+11.111336616269842j')


@settings(max_examples=1000)
@given(complex_numbers(allow_nan=False))
@example(complex())
@example(complex(-1))
@example(complex(-2))
def test_mpc_hash(c):
    assert hash(mpc(c)) == hash(c)


def test_mpc_exc():
    gmpy2.set_context(gmpy2.ieee(32))

    ctx = gmpy2.get_context()
    ctx.trap_overflow = True
    ctx.trap_underflow = True

    c = mpc(0.1 + 0.1j)

    pytest.raises(gmpy2.UnderflowResultError, lambda: c**201)
    pytest.raises(gmpy2.OverflowResultError, lambda: c**-201)

    ctx.trap_inexact = True

    pytest.raises(gmpy2.InexactResultError, lambda: mpc(0.25)**0.25)

    ctx.trap_invalid = True

    pytest.raises(gmpy2.InvalidOperationError, lambda: mpc(mpfr('nan')))


def test_issue_520():
    z = gmpy2.mpc(-0.0, 2)
    res = gmpy2.asinh(z)
    assert cmath.isclose(gmpy2.log(z + gmpy2.sqrt(1 + z*z)), res)


def test_mpc_thread_safe():
    def worker():
        test_mpc_creation()
    tpe = ThreadPoolExecutor(max_workers=20)
    for _ in range(1000):
        tpe.submit(worker)


def test_issue_650():
    x = -mpfr('nan')
    z = mpc(x, 1)

    assert gmpy2.copy_sign(mpfr(1), z.real) == -1
    assert str(z) == 'nan+1.0j'
    assert gmpy2.copy_sign(mpfr(1), z.real) == -1

    z = mpc(1, x)

    assert gmpy2.copy_sign(mpfr(1), z.imag) == -1
    assert str(z) == '1.0+nanj'
    assert gmpy2.copy_sign(mpfr(1), z.imag) == -1


@pytest.mark.skipif(sys.version_info < (3, 14), reason="requires v3.14+")
def test_mpc_mixed_arithmetics():
    from itertools import combinations_with_replacement
    from math import inf, nan
    from operator import add, mul, sub, truediv

    gmpy2.set_context(gmpy2.ieee(64))
    gmpy2.get_context().trap_divzero = True
    cases = [0.0, -0.0, inf, -inf, nan, -1.0, 1.0, 2.0, -3.0]

    for x in cases:
        gx = gmpy2.mpfr(x)
        for z in combinations_with_replacement(cases, 2):
            z = complex(*z)
            gz = gmpy2.mpc(z)
            for op in [add, mul, sub, truediv]:
                try:
                    r1 = op(x, z)
                except ZeroDivisionError:
                    try:
                        op(gx, gz)
                        assert False
                    except ZeroDivisionError:
                        g1 = r1 = 0j
                else:
                    g1 = op(gx, gz)
                try:
                    r2 = op(z, x)
                except ZeroDivisionError:
                    try:
                        op(gz, gx)
                        assert False
                    except ZeroDivisionError:
                        g2 = r2 = 0j
                else:
                    g2 = op(gz, gx)
                if op not in [sub, truediv]:
                    assert str(r1) == str(r2)
                    assert str(g1) == str(g2)
                if str(r1) != str(complex(g1)):
                    assert cmath.isfinite(r1) and x and all(_ for _ in [z.real,
                                                                        z.imag])
                    assert cmath.isclose(r1, complex(g1))
                assert str(r2) == str(complex(g2))


@pytest.mark.skipif(sys.version_info < (3, 13), reason="requires v3.13+")
def test_mpc_signatures():
    cls = gmpy2.mpc
    for f in dir(cls):
        a = getattr(cls, f)
        if callable(a) and f != '__class__':
            _ = inspect.signature(a)  # not raises
