import math
import pickle
import sys
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
from fractions import Fraction

import pytest
from hypothesis import example, given, settings
from hypothesis.strategies import floats, integers
from supportclasses import a, b, c, d, q, r, z

import gmpy2
from gmpy2 import (cmp, cmp_abs, from_binary, gamma_inc, get_context, inf,
                   is_nan, mpc, mpfr, mpfr_grandom, mpfr_nrandom, mpq, mpz,
                   nan, random_state, to_binary, xmpz, zero)


def test_mpfr_gamma_inc():
    assert gamma_inc(1, 1) == mpfr('0.36787944117144233')
    assert gamma_inc(1, 0) == mpfr('1.0')
    assert gamma_inc(0, 1) == mpfr('0.21938393439552029')


def test_mpfr_cmp():
    assert cmp(mpfr(0), mpfr(0)) == 0
    assert cmp(mpfr(0), mpz(0)) == 0
    assert cmp(mpfr(0), mpq(0,1)) == 0
    assert cmp(zero(-1), zero(-1)) == 0
    assert cmp(zero(1), zero(-1)) == 0
    assert cmp(zero(-1), zero(1)) == 0
    assert cmp(mpfr(1), mpfr(0)) == 1
    assert cmp(mpfr(1), mpz(0)) == 1
    assert cmp(mpfr(1), mpq(0,1)) == 1
    assert cmp(mpfr(-1), mpfr(0)) == -1
    assert cmp(mpfr(-1), mpz(0)) == -1
    assert cmp(mpfr(-1), mpq(0,1)) == -1
    assert cmp(nan(), mpfr(0)) == 0
    assert cmp(nan(), mpz(0)) == 0
    assert cmp(nan(), mpq(0,1)) == 0

    gmpy2.get_context().clear_flags()

    assert cmp(nan(), 1) == 0
    assert gmpy2.get_context().erange is True

    assert cmp_abs(mpfr(-1), mpfr(0)) == 1
    assert cmp_abs(mpfr(-1), mpz(0)) == 1
    assert cmp_abs(mpfr(-1), mpq(0,1)) == 1
    assert cmp_abs(mpfr(0), mpfr(-1)) == -1
    assert cmp_abs(mpz(0), mpfr(-1)) == -1
    assert cmp_abs(mpq(0,1), mpfr(-1)) == -1

    assert cmp(mpfr(1.5), q) == 0
    assert cmp(r, mpfr(1.5)) == 0


def test_mpfr_comparisons():
    from supportclasses import a, q, r

    assert mpfr(1.5) == q
    assert r == mpfr(1.5)
    assert (r == a) is False

    a = mpz(123)
    r = mpfr('inf')
    q = mpq('45/7')
    f = float(0.7)
    r2 = mpfr(454.6)

    assert (r == a, r != a, r > a, r >= a, r < a, r <= a) == (False, True, True, True, False, False)

    r = mpfr('-inf')

    assert (r == a, r != a, r > a, r >= a, r < a, r <= a) == (False, True, False, False, True, True)

    r = mpfr('nan')

    assert (r == a, r != a, r > a, r >= a, r < a, r <= a) == (False, True, False, False, False, False)
    assert (r == q, r != q, r > q, r >= q, r < q, r <= q) == (False, True, False, False, False, False)
    assert (r == f, r != f, r > f, r >= f, r < f, r <= f) == (False, True, False, False, False, False)
    assert (r == r2, r != r2, r > r2, r >= r2, r < r2, r <= r2) == (False, True, False, False, False, False)

    r = mpfr(126.5)

    assert (r == a, r != a, r > a, r >= a, r < a, r <= a) == (False, True, True, True, False, False)
    assert (r == q, r != q, r > q, r >= q, r < q, r <= q) == (False, True, True, True, False, False)

    f = float(126.5)

    assert (r == f, r != f, r > f, r >= f, r < f, r <= f) == (True, False, False, True, False, True)
    assert (r == r2, r != r2, r > r2, r >= r2, r < r2, r <= r2) == (False, True, False, False, True, True)


def test_mpfr_conversion():
    x = mpfr(a)
    assert isinstance(x, mpfr)
    assert x == 1.5
    pytest.raises(TypeError, lambda: mpfr(b))
    pytest.raises(TypeError, lambda: mpfr(c))
    pytest.raises(TypeError, lambda: mpfr(d))

    pytest.raises(OverflowError, lambda: mpz(mpfr('inf')))

    assert mpz(mpfr(5.51)) == mpz(6)

    pytest.raises(OverflowError, lambda: xmpz(mpfr('inf')))

    assert xmpz(mpfr(5.51)) == xmpz(6)

    pytest.raises(OverflowError, lambda: mpq(mpfr('inf')))

    assert float(mpfr('inf')) == float('inf')
    assert float(mpfr('-inf')) == float('-inf')
    assert float(mpfr('1e+400')) == float('inf')

    assert mpq(mpfr(4.5)) == mpq(9,2)
    assert mpq(mpfr(0)) == mpq(0,1)
    assert int(mpfr(5.3)) == 5
    assert float(mpfr(5.3)) == 5.3
    assert str(float('5.656')) == '5.656'


def test_mpfr_create():
    assert mpfr() == mpfr('0.0')
    assert mpfr(0) == mpfr('0.0')
    assert mpfr(1) == mpfr('1.0')
    assert mpfr(-1) == mpfr('-1.0')
    assert mpfr("123e17") == mpfr('1.23e+19')
    assert mpfr('1._3_5e4_5') == mpfr('1.3499999999999999e+45')
    assert mpfr('-0b1.1000000000001p+4') == mpfr('-24.001953125')
    assert mpfr('-0x1.9000000000000p+4') == mpfr('-25.0')
    assert mpfr('-0x1.9p+4') == mpfr('-25.0')

    pytest.raises(ValueError, lambda: mpfr("foo"))

    assert mpfr(1,100) == mpfr('1.0',100)
    assert mpfr("1",100) == mpfr('1.0',100)

    pytest.raises(TypeError, lambda: mpfr("1","hi"))

    assert mpfr("-inf") == mpfr('-inf')
    assert mpfr("inf") == mpfr('inf')
    assert is_nan(mpfr("nan"))
    assert mpfr(float("-inf")) == mpfr('-inf')
    assert mpfr(float("inf")) == mpfr('inf')
    assert is_nan(mpfr(float("nan")))
    assert mpfr(float("0")) == mpfr('0.0')
    assert mpfr(float("-0")) == mpfr('-0.0')
    assert mpfr("-0") == mpfr('-0.0')

    a = float("1.2345678901234567890")
    b = mpfr("1.2345678901234567890")
    assert a == 1.2345678901234567
    assert b == mpfr('1.2345678901234567')
    assert a == b

    c = mpfr(b)

    assert b is c
    assert mpfr(Fraction(0,1)) == mpfr('0.0')
    assert mpfr(Fraction(0,-1)) == mpfr('0.0')
    assert mpfr(Fraction(1,2)) == mpfr('0.5')
    assert mpfr(-1) == mpfr('-1.0')
    assert mpfr(12345678901234567890) == mpfr('1.2345678901234567e+19')
    assert mpfr(mpz(12345678901234567890)) == mpfr('1.2345678901234567e+19')
    assert mpfr(mpz(-1)) == mpfr('-1.0')
    assert mpfr(2**15 - 1) == mpfr('32767.0')
    assert mpfr(2**15) == mpfr('32768.0')
    assert mpfr(2**15 + 1) == mpfr('32769.0')
    assert mpfr(2**30 - 1) == mpfr('1073741823.0')
    assert mpfr(2**30 ) == mpfr('1073741824.0')
    assert mpfr(2**30 + 1) == mpfr('1073741825.0')

    ctx = gmpy2.get_context()
    ctx.clear_flags()
    a = mpfr("1.2")

    assert a.rc == -1

    ctx.clear_flags()
    a = mpfr("1.25")

    assert a.rc == 0

    ctx.clear_flags()
    a = mpfr('nan')

    assert is_nan(a)
    assert ctx.invalid

    ctx.clear_flags()

    assert is_nan(mpfr(a))
    assert not ctx.invalid

    ctx.clear_flags()

    assert is_nan(mpfr(float('nan')))
    assert ctx.invalid

    gmpy2.set_context(gmpy2.ieee(128))

    assert mpfr(1)/7 == mpfr('0.14285714285714285714285714285714285',113)
    assert mpfr(1.0/7) == mpfr('0.142857142857142849212692681248881854',113)
    assert mpfr(1.0/7).digits(2) == ('10010010010010010010010010010010010010010010010010010000000000000000000000000000000000000000000000000000000000000', -2, 113)

    gmpy2.set_context(gmpy2.ieee(32))

    assert mpfr(1)/7 == mpfr('0.142857149',24)
    assert (mpfr(1)/7).digits(2) == ('100100100100100100100101', -2, 24)
    assert mpfr(1.0/7) == mpfr('0.142857149',24)
    assert mpfr(1.0/7).digits(2) == ('100100100100100100100101', -2, 24)
    assert mpfr(1.0/7, precision=0) == mpfr('0.142857149',24)
    assert repr(mpfr(1.0/7, precision=1)) == "mpfr('0.14285714285714285')"
    assert repr(mpfr(1.0/7, precision=5)) == "mpfr('0.141',5)"

    pytest.raises(TypeError, lambda: mpfr(1, base=2))
    pytest.raises(TypeError, lambda: mpfr("1", s=2))


@settings(max_examples=1000)
@given(floats())
@example(0.0)
@example(1.0)
@example(2.0)
@example(-1.0)
@example(123.456)
@example(123.5)
@example(float('inf'))
@example(-float('inf'))
@example(float('nan'))
def test_mpfr_hash(x):
    if math.isnan(x):
        if sys.version_info < (3, 10):
            assert hash(mpfr(x)) == hash(x) == sys.hash_info.nan
        else:
            assert hash(mpfr(x)) != hash(x)
    else:
        assert hash(mpfr(x)) == hash(x)


@given(floats())
@example(0.0)
@example(1.0)
@example(-1.0)
@example(+float('inf'))
@example(-float('inf'))
@example(float('nan'))
@example(1.345)
def test_mpfr_to_from_binary_bulk(r):
    x = mpfr(r)
    y = from_binary(to_binary(x))
    assert x == y or (is_nan(x) and is_nan(y))


def test_mpfr_to_from_binary():
    x = mpfr("1.345e1000")
    assert x == from_binary(to_binary(x))
    assert -x == from_binary(to_binary(-x))
    x = mpfr("1e1234567890123456789")
    assert x == from_binary(to_binary(x))
    assert x.rc == 1
    x = mpfr("-1e1234567890123456789")
    assert x == from_binary(to_binary(x))
    assert x.rc == -1
    x = gmpy2.const_pi()
    assert x.rc == -1
    y = from_binary(to_binary(x))
    assert x == y and y.rc == -1
    -1
    x = gmpy2.const_pi(precision=104)
    assert x.rc == 1
    y = from_binary(to_binary(x))
    assert x == y and y.rc == 1
    with gmpy2.context() as ctx:
        ctx.precision = 20
        x = gmpy2.const_pi()
        assert x == from_binary(to_binary(x))
        ctx.precision = 100
        x = gmpy2.const_pi()
        assert x == from_binary(to_binary(x))
        ctx.precision = 200
        x = mpfr(gmpy2.const_pi())
        assert x == from_binary(to_binary(x))
        x = gmpy2.const_pi()
        ctx.precision = 300
        x = from_binary(to_binary(x))
        assert x.precision == 200


def test_mpfr_random():
    assert gmpy2.mpfr_random(random_state(42)) == mpfr('0.93002690534702315')


def test_mpfr_grandom():
    assert mpfr_grandom(random_state(42)) == (mpfr('-0.32898912492644183'),
                                              mpfr('0.03656576719642516'))


def test_mpfr_nrandom():
    assert mpfr_nrandom(random_state(42)) == mpfr('-0.32898912492644183')


def test_mpfr_mpmath():
    mpmath = pytest.importorskip("mpmath")
    a, b, c, d = '1.1', '-1.1', '-3.14', '0'
    assert mpfr(a)._mpf_ == (0, mpz(4953959590107546), -52, 53)
    assert mpmath.mpf(mpfr(a)) == mpmath.mpf(a)
    assert mpfr(b)._mpf_ == (1, mpz(4953959590107546), -52, 53)
    assert mpmath.mpf(mpfr(b)) == mpmath.mpf(b)
    assert mpfr(c, precision=10)._mpf_ == (1, mpz(804), -8, 10)
    assert mpmath.mpf(mpfr(c, precision=10), prec=10) == mpmath.mpf(c, prec=10)
    assert mpfr(d)._mpf_ == (0, mpz(0), 1, 1)
    pytest.raises(TypeError, lambda: gmpy2._mpmath_create(1))
    pytest.raises(ValueError, lambda: gmpy2._mpmath_create(1, 1, -1))
    pytest.raises(TypeError, lambda: mpmath.mpf(("!", 1,)))


def test_mpfr_format():
    r, r1, r2 = mpfr(5.6), mpfr(-3), mpfr(5)

    assert '{:<30}'.format(r1) == '-3.000000                     '
    assert '{:>+20}'.format(r2) == '           +5.000000'
    assert '{:>-15}'.format(r2) == '       5.000000'
    assert '{:>-15}'.format(r1) == '      -3.000000'
    assert '{:U}'.format(r) == '5.600000'

    pytest.raises(ValueError, lambda: '{:U-}'.format(r))

    assert '{:Z}'.format(r) == '5.599999'

    pytest.raises(ValueError, lambda: '{:Z+}'.format(r))

    assert '{:+Z}'.format(r) == '+5.599999'

    pytest.raises(ValueError, lambda: '{:Z }'.format(r))

    assert '{:.10f}'.format(r) == '5.6000000000'
    assert '{:.10f.}'.format(r) == '5.6000000000'

    assert "{:.0f}".format(mpfr('123')) == '123.0'

    pytest.raises(ValueError, lambda: '{:Z.}'.format(r))
    pytest.raises(ValueError, lambda: '{:->}'.format(r))
    pytest.raises(ValueError, lambda: '{:YZ}'.format(r))

    # issue 503
    r = mpfr(2.675)
    assert f'{r:.2f}' == '2.67'
    gmpy2.set_context(gmpy2.context(round=gmpy2.RoundUp))
    assert f'{r:.2f}' == '2.68'


def test_mpfr_digits():
    r, r2 = mpfr(5.6), mpfr(5)

    assert r.digits() == ('55999999999999996', 1, 53)
    assert r.digits(2) == ('10110011001100110011001100110011001100110011001100110', 3, 53)
    assert r.digits(2,54) == ('101100110011001100110011001100110011001100110011001100', 3, 53)
    assert r.digits(10,54) == ('559999999999999964472863211994990706443786621093750000', 1, 53)
    assert r2.digits(2) == ('10100000000000000000000000000000000000000000000000000', 3, 53)

    pytest.raises(TypeError, lambda: r2.digits(2, 5, 6))
    pytest.raises(ValueError, lambda: r.digits(0))


def test_mpfr_abs():
    a = mpfr(1.0)
    b = abs(a)

    assert a is not b
    assert abs(mpfr(1, precision=100)) == mpfr('1.0')

    ctx = gmpy2.get_context()
    ctx.clear_flags()

    assert is_nan(abs(mpfr('nan')))
    assert ctx.invalid

    ctx.clear_flags()

    assert abs(mpfr('inf')) == mpfr('inf')
    assert abs(mpfr('-inf')) == mpfr('inf')


def test_mpfr_not():
    assert not mpfr(0)
    assert not mpfr(1) is False


def test_mpfr_add():
    a = mpfr("12.34")
    b = mpfr("45.67")

    assert a+1 == mpfr('13.34')
    assert a+1.0 == mpfr('13.34')
    assert a+mpz(1) == mpfr('13.34')
    assert a+mpq(1,1) == mpfr('13.34')
    assert 1+a == mpfr('13.34')
    assert 1.0+a == mpfr('13.34')
    assert mpz(1)+a == mpfr('13.34')
    assert mpq(1,1)+a == mpfr('13.34')
    assert a+b == mpfr('58.010000000000005')
    assert a+0 == mpfr('12.34')
    assert 0+a == mpfr('12.34')
    assert b+a == mpfr('58.010000000000005')

    pytest.raises(TypeError, lambda: a+'b')
    pytest.raises(TypeError, lambda: 'b'+a)

    assert a == 12.34
    assert b == 45.67
    assert a+b == 12.34+45.67

    assert a + float('Inf') == mpfr('inf')
    assert float('Inf') + a == mpfr('inf')
    assert a + float('-Inf') == mpfr('-inf')
    assert float('-Inf') + a == mpfr('-inf')
    assert is_nan(a + float('nan'))
    assert is_nan(float('nan') + a)

    assert a + mpfr('Inf') == mpfr('inf')
    assert mpfr('Inf') + a == mpfr('inf')
    assert a + mpfr('-Inf') == mpfr('-inf')
    assert mpfr('-Inf') + a == mpfr('-inf')
    assert is_nan(a + mpfr('nan'))
    assert is_nan(mpfr('nan') + a)


def test_mpfr_sub():
    assert mpfr(10) - 1 == mpfr('9.0')
    assert 10 - mpfr(1) == mpfr('9.0')
    assert mpfr(10) - mpz(1) == mpfr('9.0')
    assert mpz(10) - mpfr(1) == mpfr('9.0')
    assert mpfr(10) - mpfr(1) == mpfr('9.0')
    assert mpfr(10) - mpq(1,1) == mpfr('9.0')
    assert mpq(10,1) - mpfr(1) == mpfr('9.0')
    assert mpfr(10) - Fraction(1,1) == mpfr('9.0')
    assert Fraction(10,1) - mpfr(1) == mpfr('9.0')
    assert mpfr(10) - 1.0 == mpfr('9.0')
    assert 10.0 - mpfr(1) == mpfr('9.0')
    assert mpfr(0) - (1 << 100) == mpfr('-1p100', base=2)
    assert (1 << 100) - mpfr(0) == mpfr('1p100', base=2)
    assert mpfr(10) - z == mpfr('8.0')
    assert mpfr(10) - q == mpfr('8.5')
    assert mpfr(10) - r == mpfr('8.5')

    a = mpfr("12.34")
    b = mpfr("45.67")

    assert a-1 == mpfr('11.34')
    assert a-1.0 == mpfr('11.34')
    assert a-(-1) == mpfr('13.34')
    assert a-(-1.0) == mpfr('13.34')
    assert a-b == mpfr('-33.329999999999998')
    assert b-a == mpfr('33.329999999999998')

    pytest.raises(TypeError, lambda: a-'b')
    pytest.raises(TypeError, lambda: 'b'-a)

    assert a-b == 12.34-45.67

    assert a - float('Inf') == mpfr('-inf')
    assert float('Inf') - a == mpfr('inf')
    assert a - float('-Inf') == mpfr('inf')
    assert float('-Inf') - a == mpfr('-inf')
    assert is_nan(a - float('nan'))
    assert is_nan(float('nan') - a)

    assert a - mpfr('Inf') == mpfr('-inf')
    assert mpfr('Inf') - a == mpfr('inf')
    assert a - mpfr('-Inf') == mpfr('inf')
    assert mpfr('-Inf') - a == mpfr('-inf')
    assert is_nan(a - mpfr('nan'))
    assert is_nan(mpfr('nan') - a)


def test_mpfr_mul():
    c = 12345678901234567890

    assert mpfr(10) * 1 == mpfr('10.0')
    assert 10 * mpfr(1) == mpfr('10.0')
    assert mpfr(10) * mpz(1) == mpfr('10.0')
    assert mpz(10) * mpfr(1) == mpfr('10.0')
    assert mpfr(10) * mpfr(1) == mpfr('10.0')
    assert mpfr(10) * mpq(1,1) == mpfr('10.0')
    assert mpq(10,1) * mpfr(1) == mpfr('10.0')
    assert mpfr(10) * Fraction(1,1) == mpfr('10.0')
    assert Fraction(10,1) * mpfr(1) == mpfr('10.0')
    assert mpfr(10) * 1.0 == mpfr('10.0')
    assert 10.0 * mpfr(1) == mpfr('10.0')
    assert mpfr(1) * c == mpfr(c)
    assert c * mpfr(1) == mpfr(c)
    assert mpfr(10) * z == mpfr('20.0')
    assert mpfr(10) * q == mpfr('15.0')
    assert mpfr(10) * r == mpfr('15.0')

    pytest.raises(TypeError, lambda: mpfr(10) * 'a')
    pytest.raises(TypeError, lambda: 'a' * mpfr(10))

    a = mpfr("12.34")
    b = mpfr("45.67")

    assert a*b == mpfr('563.56780000000003')
    assert a*0 == mpfr('0.0')
    assert a*0.0 == mpfr('0.0')
    assert a*1 == mpfr('12.34')
    assert a*1.0 == mpfr('12.34')
    assert a*(-1.0) == mpfr('-12.34')
    assert 0*a == mpfr('0.0')
    assert 0.0*a == mpfr('0.0')
    assert 1*a == mpfr('12.34')
    assert 1.0*a == mpfr('12.34')
    assert a*b == mpfr('563.56780000000003')
    assert a*b == 12.34*45.67

    assert a * float('Inf') == mpfr('inf')
    assert float('Inf') * a == mpfr('inf')
    assert a * float('-Inf') == mpfr('-inf')
    assert float('-Inf') * a == mpfr('-inf')
    assert -a * float('Inf') == mpfr('-inf')
    assert float('Inf') * -a == mpfr('-inf')
    assert -a * float('-Inf') == mpfr('inf')
    assert float('-Inf') * -a == mpfr('inf')
    assert is_nan(a * float('nan'))
    assert is_nan(float('nan') * a)
    assert is_nan(mpfr(0) * float('Inf'))
    assert is_nan(mpfr(0) * float('-Inf'))
    assert is_nan(float('Inf') * mpfr(0))
    assert is_nan(float('-Inf') * mpfr(0))

    assert a * mpfr('Inf') == mpfr('inf')
    assert mpfr('Inf') * a == mpfr('inf')
    assert a * mpfr('-Inf') == mpfr('-inf')
    assert mpfr('-Inf') * a == mpfr('-inf')
    assert -a * mpfr('Inf') == mpfr('-inf')
    assert mpfr('Inf') * -a == mpfr('-inf')
    assert -a * mpfr('-Inf') == mpfr('inf')
    assert mpfr('-Inf') * -a == mpfr('inf')
    assert is_nan(a * mpfr('nan'))
    assert is_nan(mpfr('nan') * a)
    assert is_nan(mpz(0) * mpfr('Inf'))
    assert is_nan(mpz(0) * mpfr('-Inf'))
    assert is_nan(mpfr('Inf') * mpfr(0))
    assert is_nan(mpfr('-Inf') * mpfr(0))


def test_mpfr_divmod():
    pytest.raises(TypeError, lambda: divmod(mpfr(1),'a'))

    ctx = gmpy2.context()

    assert ctx.divmod(mpfr(2),mpfr(1)) == (mpfr('2.0'), mpfr('0.0'))

    pytest.raises(TypeError, lambda: divmod(mpfr(1), mpc(1,2)))

    assert divmod(mpfr(1.2), mpfr(0.7)) == (mpfr('1.0'), mpfr('0.5'))

    ctx = gmpy2.ieee(64)
    gmpy2.set_context(ctx)

    assert divmod(mpfr(1.2), mpfr(0.7)) == (mpfr('1.0'), mpfr('0.5'))
    ctx.clear_flags()

    assert ctx.divzero is False
    assert all(map(ctx.is_nan, divmod(mpfr(1.2), mpfr(0))))

    with gmpy2.context(trap_divzero=True):
        pytest.raises(gmpy2.DivisionByZeroError, lambda: divmod(mpfr(1), mpfr(0)))
    with gmpy2.context(trap_invalid=True):
        pytest.raises(gmpy2.InvalidOperationError, lambda: divmod(gmpy2.nan(), mpfr(1)))
    with gmpy2.context(trap_invalid=True):
        pytest.raises(gmpy2.InvalidOperationError, lambda: divmod(mpfr(1), gmpy2.inf()))

    assert divmod(mpfr(111), mpfr(-222)) == (mpfr('-1.0'), mpfr('-111.0'))

    a = mpfr("12.34")
    b = mpfr("45.67")

    assert divmod(12.34, 45.67) == (0.0, 12.34)
    assert divmod(a,b) == (mpfr('0.0'), mpfr('12.34'))
    assert divmod(b,a) == (mpfr('3.0'), mpfr('8.6500000000000021'))
    assert divmod(45.67,12.34) == divmod(b,a)

    assert divmod(a, float('Inf')) == (mpfr('0.0'), mpfr('12.34'))
    assert divmod(a, float('-Inf')) == (mpfr('-1.0'), mpfr('-inf'))
    assert divmod(-a, float('Inf')) == (mpfr('-1.0'), mpfr('inf'))
    assert divmod(-a, float('-Inf')) == (mpfr('0.0'), mpfr('-12.34'))
    assert all(map(is_nan, divmod(a, float('nan'))))
    assert all(map(is_nan, divmod(-a, float('nan'))))
    assert divmod(mpfr(0), float('Inf')) == (mpfr('0.0'), mpfr('0.0'))
    assert divmod(mpfr(0), float('-Inf')) == (mpfr('-0.0'), mpfr('-0.0'))
    assert divmod(mpfr(0), float('nan'))
    assert all(map(is_nan, divmod(float('Inf'), a)))
    assert all(map(is_nan, divmod(float('-Inf'), a)))
    assert all(map(is_nan, divmod(float('Inf'), -a)))
    assert all(map(is_nan, divmod(float('-Inf'), -a)))
    assert all(map(is_nan, divmod(float('nan'), a)))
    assert all(map(is_nan, divmod(float('nan'), -a)))
    assert all(map(is_nan, divmod(float('Inf'), mpfr(0))))
    assert all(map(is_nan, divmod(float('-Inf'), mpfr(0))))
    assert all(map(is_nan, divmod(float('nan'), mpfr(0))))

    assert divmod(a, mpfr('Inf')) == (mpfr('0.0'), mpfr('12.34'))
    assert divmod(a, mpfr('-Inf')) == (mpfr('-1.0'), mpfr('-inf'))
    assert divmod(-a, mpfr('Inf')) == (mpfr('-1.0'), mpfr('inf'))
    assert divmod(-a, mpfr('-Inf')) == (mpfr('0.0'), mpfr('-12.34'))
    assert all(map(is_nan, divmod(a, mpfr('nan'))))
    assert all(map(is_nan, divmod(-a, mpfr('nan'))))
    assert divmod(mpfr(0), mpfr('Inf')) == (mpfr('0.0'), mpfr('0.0'))
    assert divmod(mpfr(0), mpfr('-Inf')) == (mpfr('-0.0'), mpfr('-0.0'))
    assert divmod(mpfr(0), mpfr('nan'))
    assert all(map(is_nan, divmod(mpfr('Inf'), a)))
    assert all(map(is_nan, divmod(mpfr('-Inf'), a)))
    assert all(map(is_nan, divmod(mpfr('Inf'), -a)))
    assert all(map(is_nan, divmod(mpfr('-Inf'), -a)))
    assert all(map(is_nan, divmod(mpfr('nan'), a)))
    assert all(map(is_nan, divmod(mpfr('nan'), -a)))
    assert all(map(is_nan, divmod(mpfr('Inf'), mpfr(0))))
    assert all(map(is_nan, divmod(mpfr('-Inf'), mpfr(0))))
    assert all(map(is_nan, divmod(mpfr('nan'), mpfr(0))))


def test_mpfr_floordiv():
    ctx = gmpy2.get_context()
    a, b = mpz(45), mpz(6)
    r, r2 = mpfr(45), mpfr(3.1)
    q, q2 = mpq(118,18), mpq(3,2)
    pyq, pyq2 = Fraction(118,18), Fraction(3,2)
    c, c2 = mpc(51, 65), mpc(4, 6)

    assert ctx.floor_div(r, r2) == mpfr('14.0')
    assert ctx.floor_div(r, r2) == mpfr('14.0')
    assert ctx.floor_div(r, 3.1) == mpfr('14.0')
    assert ctx.floor_div(r, 4) == mpfr('11.0')
    assert ctx.floor_div(r, b) == mpfr('7.0')
    assert ctx.floor_div(r, q2) == mpfr('30.0')
    assert ctx.floor_div(r, pyq2) == mpfr('30.0')
    assert ctx.floor_div(r, 0) == mpfr('inf')
    assert ctx.floor_div(r, mpz(0)) == mpfr('inf')
    assert ctx.floor_div(45.0, r2) == mpfr('14.0')
    assert ctx.floor_div(45.0, r2) == mpfr('14.0')
    assert ctx.floor_div(45.0, 3.1) == mpfr('14.0')
    assert ctx.floor_div(45.0, 4) == mpfr('11.0')
    assert ctx.floor_div(45.0, b) == mpfr('7.0')
    assert ctx.floor_div(45.0, q2) == mpfr('30.0')
    assert ctx.floor_div(45.0, pyq2) == mpfr('30.0')
    assert ctx.floor_div(45.0, 0) == mpfr('inf')
    assert ctx.floor_div(45.0, mpz(0)) == mpfr('inf')
    assert ctx.floor_div(45, r2) == mpfr('14.0')

    pytest.raises(TypeError, lambda: ctx.floor_div(r, 'not'))

    assert r // b == mpfr('7.0')
    assert r // q2 == mpfr('30.0')
    assert r // r2 == mpfr('14.0')

    pytest.raises(TypeError, lambda: r // c2)
    pytest.raises(TypeError, lambda: r // 'not')

    a = mpfr("12.34")
    b = mpfr("45.67")

    assert a//b == mpfr('0.0')
    assert b//a == mpfr('3.0')


def test_mpfr_truediv():
    a = mpfr("12.34")
    b = mpfr("45.67")

    assert a/b == mpfr('0.27019925552879348')
    assert gmpy2.div(a, b) == mpfr('0.27019925552879348')
    assert get_context().div(a, b) == mpfr('0.27019925552879348')
    assert b/a == mpfr('3.7009724473257699')
    assert a/b == 12.34/45.67
    assert a/0 == mpfr('inf')
    assert a/(-0.0) == mpfr('-inf')
    assert a/b == 12.34/45.67
    assert mpfr(10) / z == mpfr('5.0')
    assert mpfr(10) / q == mpfr('6.666666666666667')
    assert mpfr(10) / r == mpfr('6.666666666666667')
    assert b / mpc(4, 4) == mpc('5.7087500000000002-5.7087500000000002j')
    assert get_context().div(b, mpc(4, 4)) == mpc('5.7087500000000002-5.7087500000000002j')

    ans = b / mpc(0, 0)

    assert ans.real == mpfr('inf') and is_nan(ans.imag)

    ans = get_context().div(b, mpc(0, 0))

    assert ans.real == mpfr('inf') and is_nan(ans.imag)

    pytest.raises(TypeError, lambda: b / 'str')
    pytest.raises(TypeError, lambda: get_context().div(b, 'str'))

    assert a / float('Inf') == mpfr('0.0')
    assert -a / float('Inf') == mpfr('-0.0')
    assert float('Inf') / a == mpfr('inf')
    assert float('Inf') / -a == mpfr('-inf')
    assert a / float('-Inf') == mpfr('-0.0')
    assert -a / float('-Inf') == mpfr('0.0')
    assert float('-Inf') / a == mpfr('-inf')
    assert float('-Inf') / -a == mpfr('inf')
    assert is_nan(a / float('nan'))
    assert is_nan(float('nan') / a)
    assert is_nan(float('nan') / mpfr(0))
    assert is_nan(float('nan') / mpfr(0))

    assert a / mpfr('Inf') == mpfr('0.0')
    assert -a / mpfr('Inf') == mpfr('-0.0')
    assert mpfr('Inf') / a == mpfr('inf')
    assert mpfr('Inf') / -a == mpfr('-inf')
    assert a / mpfr('-Inf') == mpfr('-0.0')
    assert -a / mpfr('-Inf') == mpfr('0.0')
    assert mpfr('-Inf') / a == mpfr('-inf')
    assert mpfr('-Inf') / -a == mpfr('inf')
    assert is_nan(a / mpfr('nan'))
    assert is_nan(mpfr('nan') / a)
    assert is_nan(mpfr('nan') / mpfr(0))
    assert is_nan(mpfr('nan') / mpfr(0))


def test_mpfr_mod():
    r = mpfr('0.0') % mpfr('-1.0')
    assert r.is_zero() and r.is_signed()

    a = mpfr("12.34")
    b = mpfr("45.67")

    assert a % 1 == mpfr('0.33999999999999986')
    assert 12.34 % 1 == 0.33999999999999986
    assert a % z == mpfr('0.33999999999999986')
    assert is_nan(b % 0)
    assert is_nan(b % mpz(0))
    assert is_nan(get_context().mod(b, 0))
    assert is_nan(get_context().mod(b, mpz(0)))
    assert is_nan(b % mpfr(0))
    assert is_nan(get_context().mod(b, mpfr(0)))

    get_context().trap_divzero = True

    pytest.raises(gmpy2.DivisionByZeroError, lambda: b % 0)
    pytest.raises(gmpy2.DivisionByZeroError, lambda: b % mpz(0))
    pytest.raises(gmpy2.DivisionByZeroError, lambda: get_context().mod(b, 0))
    pytest.raises(gmpy2.DivisionByZeroError, lambda: get_context().mod(b, mpz(0)))
    pytest.raises(gmpy2.DivisionByZeroError, lambda: b % mpfr(0))
    pytest.raises(gmpy2.DivisionByZeroError, lambda: get_context().mod(b, mpfr(0)))

    get_context().trap_divzero = False

    pytest.raises(TypeError, lambda: b % 'str')
    pytest.raises(TypeError, lambda: get_context().mod(b,'str'))
    pytest.raises(TypeError, lambda: b % mpc(3, 5))

    assert is_nan(mpfr('nan') % a)
    assert is_nan(a % mpfr('nan'))
    assert is_nan(mpfr('inf') % a)
    assert a % mpfr('inf') == mpfr('12.34')

    get_context().trap_invalid = True

    pytest.raises(gmpy2.InvalidOperationError, lambda: mpfr('nan') % a)
    pytest.raises(gmpy2.InvalidOperationError, lambda: a % mpfr('nan'))
    pytest.raises(gmpy2.InvalidOperationError, lambda: mpfr('inf') % a)
    pytest.raises(gmpy2.InvalidOperationError, lambda: a % mpfr('inf'))

    get_context().trap_invalid = False

    pytest.raises(TypeError, lambda: get_context().mod(a, b, 5))

    assert get_context().mod(a, -mpfr('inf')) == mpfr('-inf')


def test_mpfr_pow():
    r1, r2 = mpfr(5.0), mpfr(2.5)
    ctx = gmpy2.get_context()

    assert r1 ** mpz(2) == mpfr('25.0')
    assert r2 ** mpz(2) == mpfr('6.25')
    assert r2 ** 2 == mpfr('6.25')
    assert pow(r1, r2) == mpfr('55.901699437494742')
    assert ctx.pow(r1, r2) == mpfr('55.901699437494742')
    assert ctx.pow(r1, r2) == r1 ** r2

    pytest.raises(TypeError, lambda: pow(r1, r2, 5))
    pytest.raises(TypeError, lambda: ctx.pow(r1, r2, 5))

    assert pow(r1, 4) == mpfr('625.0')
    assert ctx.pow(r1, 4) == mpfr('625.0')


def test_mpfr_subnormalize():
    gmpy2.set_context(gmpy2.ieee(64))

    zeroes = '0.' + '0' * 323

    assert mpfr(zeroes + '2470328229206232720') == mpfr('0.0')
    assert mpfr(zeroes + '2470328229206232721') == mpfr('4.9406564584124654e-324')
    assert mpfr(zeroes + '247032822920623272088') == mpfr('0.0')
    assert mpfr(zeroes + '247032822920623272089') == mpfr('4.9406564584124654e-324')

    gmpy2.set_context(gmpy2.ieee(32))

    def fmt(v):
        return '{:.23b}'.format(v)


    a = mpfr('0x1p-126')

    assert fmt(a) == '1.00000000000000000000000p-126'
    assert fmt(gmpy2.next_below(a)) == '1.11111111111111111111110p-127'
    assert fmt(gmpy2.next_above(a)) == '1.00000000000000000000001p-126'
    assert fmt(gmpy2.next_below(0)) == '-1.00000000000000000000000p-149'
    assert fmt(gmpy2.next_above(0)) == '1.00000000000000000000000p-149'
    assert fmt(gmpy2.next_toward(a, -10)) == '1.11111111111111111111110p-127'
    assert fmt(gmpy2.next_toward(a, 10)) == '1.00000000000000000000001p-126'

    gmpy2.set_context(gmpy2.context())


def test_mpfr_as_integer_ratio():
    assert mpfr('1.1e+2').as_integer_ratio() == (mpz(110), mpz(1))

    pytest.raises(ValueError, lambda: mpfr('nan').as_integer_ratio())
    pytest.raises(OverflowError, lambda: mpfr('inf').as_integer_ratio())


def test_mpfr_round():
    pytest.raises(TypeError, lambda: round(mpfr('1.0'), "spam"))

    r = round(mpfr('-0.0'), 123)
    assert r.is_zero() and r.is_signed()

    assert round(mpfr('12.34'), -1) == mpfr('10.0')

    r = mpfr(4.55)

    assert round(r, 1) == mpfr('4.5')
    assert round(mpfr('5.1')) == mpz(5)

    pytest.raises(ValueError, lambda: round(nan()))
    pytest.raises(OverflowError, lambda: round(inf()))


def test_mpfr_as_mantissa_exp():
    r = mpfr(4.55)

    assert r.as_mantissa_exp() == (mpz(5122844576133939), mpz(-50))
    assert mpfr(0).as_mantissa_exp() == (mpz(0), mpz(1))

    pytest.raises(OverflowError, lambda: inf().as_mantissa_exp())
    pytest.raises(ValueError, lambda: nan().as_mantissa_exp())


def test_mpfr_as_simple_fraction():
    r = mpfr(4.55)

    assert r.as_simple_fraction() == mpq(91,20)

    pytest.raises(ValueError, lambda: r.as_simple_fraction(precision=100))
    pytest.raises(TypeError, lambda: r.as_simple_fraction(nosense=10))


def test_mpfr_real_imag():
    r = mpfr(4.55)
    a = mpfr('12.34')

    assert r.imag == mpfr('0.0')
    assert r.real == mpfr('4.5499999999999998')
    assert a.real == mpfr('12.34')
    assert a.imag == mpfr('0.0')


def test_mpfr_conjugate():
    r = mpfr(4.55)
    a = mpfr('12.34')

    assert r.conjugate() == r
    assert a.conjugate() == a


def test_mpfr_pickle():
    a = mpfr('12.34')

    assert pickle.loads(pickle.dumps(a)) == a
    assert pickle.loads(pickle.dumps(mpfr("inf"))) == mpfr('inf')
    assert pickle.loads(pickle.dumps(mpfr("-inf"))) == mpfr('-inf')
    assert is_nan(pickle.loads(pickle.dumps(mpfr("nan"))))
    assert pickle.loads(pickle.dumps(mpfr(0))) == mpfr('0.0')


def test_mpfr_floor():
    a = mpfr('12.34')

    ctx = get_context()
    r = ctx.floor(a)

    assert r == mpz(12) and isinstance(r, mpz)


def test_mpfr_ceil():
    a = mpfr('12.34')

    ctx = get_context()
    r = ctx.ceil(a)

    assert r == mpz(13) and isinstance(r, mpz)


def test_mpfr_trunc():
    a = mpfr('12.34')

    ctx = get_context()
    r = ctx.trunc(a)

    assert r == mpz(12) and isinstance(r, mpz)


@given(floats(allow_nan=False, allow_infinity=False),
       integers(-10, 40))
def test_mpfr_round_roundtrip_bulk(x, n):
    q = mpq(*x.as_integer_ratio())
    assert float(round(q, n)) == round(x, n)
    assert mpfr(round(q, n)) == round(mpfr(x), n)


def test_issue_540():
    a = mpfr('1')
    b = mpfr('10')
    ctxD = gmpy2.context(round=gmpy2.RoundDown)

    assert ctxD.div(a, b) == mpfr('0.099999999999999992')


def test_mpfr_thread_safe():
    def worker():
        ctx = gmpy2.get_context()
        ctx.clear_flags()
        a = mpfr("1.2")
        assert a.rc == -1
        assert ctx.inexact
        ctx.clear_flags()
        a = mpfr("1")
        assert a == 1
        assert a.rc == 0
        assert ctx.inexact is False
        ctx.clear_flags()
        a = mpfr("2.1")
        assert a == 2.1
        assert a.rc == 1
        assert ctx.inexact
    tpe = ThreadPoolExecutor(max_workers=20)
    for _ in range(1000):
        tpe.submit(worker)


def test_issue_650():
    x = -mpfr('nan')

    assert gmpy2.copy_sign(mpfr(1), x) == -1
    assert str(x) == 'nan'
    assert gmpy2.copy_sign(mpfr(1), x) == -1
