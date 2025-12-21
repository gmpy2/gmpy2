import math
import numbers
import pickle
import sys
from decimal import Decimal
from fractions import Fraction

import pytest
from hypothesis import example, given, settings
from hypothesis.strategies import fractions, integers
from supportclasses import a, b, c, d, q, z

import gmpy2
from gmpy2 import (cmp, cmp_abs, from_binary, is_nan, mpc, mpfr, mpq, mpz,
                   to_binary, xmpz)


def test_mpq_constructor():
    assert mpq('1/1') == mpq(1,1)
    assert mpq('1_/_1') == mpq(1,1)

    pytest.raises(TypeError, lambda: mpq('1', s=1))


def test_mpq_as_integer_ratio():
    assert mpq(2, 3).as_integer_ratio() == (mpz(2), mpz(3))


def test_mpq_numbers_abc():
    assert isinstance(mpq(1, 3), numbers.Rational)


def test_mpq_pickling():
    for proto in range(pickle.HIGHEST_PROTOCOL + 1):
        for x in [mpq(1234, 6789), mpq(-1234, 6789), mpq(0, 1)]:
            assert pickle.loads(pickle.dumps(x, protocol=proto)) == x


def test_mpq_from_float():
    assert mpq.from_float(3.2) == mpq(3602879701896397, 1125899906842624)


def test_mpq_float():
    assert float(mpq(9, 5)) == 1.8


DBL_MAX_EXP = sys.float_info.max_exp
DBL_MIN_EXP = sys.float_info.min_exp
DBL_MANT_DIG = sys.float_info.mant_dig
DBL_MIN_OVERFLOW = 2**DBL_MAX_EXP - 2**(DBL_MAX_EXP - DBL_MANT_DIG - 1)


# Pure Python version of correctly-rounded true division,
# copied from CPython test suite (Lib/test/test_long.py).
# The CPython's true division of integers suffers from double
# rounding on x86 machines that operate with e x87 FPU set
# to 64-bit precision.  See python/cpython#142449.
def truediv(a, b):
    """Correctly-rounded true division for integers."""
    negative = a^b < 0
    a, b = abs(a), abs(b)

    # exceptions:  division by zero, overflow
    if not b:
        raise ZeroDivisionError("division by zero")
    if a >= DBL_MIN_OVERFLOW * b:
        raise OverflowError("int/int too large to represent as a float")

   # find integer d satisfying 2**(d - 1) <= a/b < 2**d
    d = a.bit_length() - b.bit_length()
    if d >= 0 and a >= 2**d * b or d < 0 and a * 2**-d >= b:
        d += 1

    # compute 2**-exp * a / b for suitable exp
    exp = max(d, DBL_MIN_EXP) - DBL_MANT_DIG
    a, b = a << max(-exp, 0), b << max(exp, 0)
    q, r = divmod(a, b)

    # round-half-to-even: fractional part is r/b, which is > 0.5 iff
    # 2*r > b, and == 0.5 iff 2*r == b.
    if 2*r > b or 2*r == b and q % 2 == 1:
        q += 1

    result = math.ldexp(q, exp)
    return -result if negative else result


@settings(max_examples=10000)
@given(fractions())
@example(Fraction(12143, 31517))
def test_mpq_float_bulk(x):
    q = mpq(x)
    try:
        fx = truediv(x.numerator, x.denominator)
    except OverflowError:
        with pytest.raises(OverflowError):
            float(q)
    else:
        assert fx == float(q)


def test_mpq_from_Decimal():
    assert mpq(Decimal("5e-3")) == mpq(5, 1000)
    assert mpq(Decimal(1)) == mpq(1)  # issue 327
    assert mpq(Decimal('0.6')) == mpq(3, 5)
    assert mpq.from_decimal(Decimal("5e-3")) == mpq(5, 1000)


def test_mpq_cmp():
    assert cmp(mpq(1,2), mpq(1,2)) == 0
    assert cmp(mpq(-1,2), mpq(1,2)) == -1
    assert cmp(mpq(1,2), mpq(-1,2)) == 1
    assert cmp(0, mpq(1,2)) == -1
    assert cmp(mpq(1,2), 0) == 1

    assert cmp_abs(mpq(1,2), mpq(1,3)) == 1
    assert cmp_abs(mpq(-1,2), mpq(1,3)) == 1
    assert cmp_abs(mpq(1,2), mpq(-1,3)) == 1
    assert cmp_abs(mpq(1,4), mpq(1,3)) == -1
    assert cmp_abs(mpq(-1,4), mpq(1,3)) == -1
    assert cmp_abs(mpq(1,4), mpq(-1,3)) == -1

    assert cmp(mpq(3,2), q) == 0
    assert cmp(q, mpq(3,5)) == 1


def test_mpq_comparisons():
    from supportclasses import q

    assert mpq(3,2) == q
    assert (q == mpq(3,5)) is False

    a = mpz(123)
    q = mpq(4, 5)

    assert (q == a, q != a, q > a, q >= a, q < a, q <= a) == (False, True, False, False, True, True)
    assert (mpq(246,2) != a) is False

    f = float(0.7)

    assert (q == f, q != f, q > f, q >= f, q < f, q <= f) == (False, True, True, True, False, False)

    f = float('nan')

    assert (q == f, q != f, q > f, q >= f, q < f, q <= f) == (False, True, False, False, False, False)

    f = float('inf')

    assert (q == f, q != f, q > f, q >= f, q < f, q <= f) == (False, True, False, False, True, True)

    f = -f

    assert (q == f, q != f, q > f, q >= f, q < f, q <= f) == (False, True, True, True, False, False)


def test_mpq_conversion():
    x = mpq(a)
    assert isinstance(x, mpq)
    assert 2*x == 3
    assert mpq(z, z) == mpq(1, 1)

    class Three:
        def __mpz__(self): return mpz(3)

    assert mpq(z, Three()) == mpq(2,3)

    pytest.raises(TypeError, lambda: mpq(b))
    pytest.raises(TypeError, lambda: mpq(c))
    pytest.raises(TypeError, lambda: mpq(d))

    assert mpq('2/3') == mpq(2,3)
    assert mpq(b'2/3') == mpq(2,3)

    pytest.raises(ValueError, lambda: mpq('2,3'))
    pytest.raises(ValueError, lambda: mpq('2/a'))

    assert mpq('2.3') == mpq(23,10)

    assert pytest.raises(ValueError, lambda: mpq('2.3/10'))

    assert mpq(4.5) == mpq(9,2)

    pytest.raises(OverflowError, lambda: mpq(float('inf')))

    assert mpq(xmpz(15)) == mpq(15,1)
    assert mpq(mpfr(4.5)) == mpq(9,2)

    pytest.raises(TypeError, lambda: mpq(dict()))

    assert float(mpq(1,2)) == 0.5
    assert int(mpq(15,2)) == 7


def test_mpq_round():
    pytest.raises(TypeError, lambda: round(mpq(7,3), 4.5))

    q = mpq('4/5')

    assert round(mpq('7/2')) == mpz(4)
    assert round(q, 4) == mpq(4,5)

    pytest.raises(TypeError, lambda: round(q, 4, 2))


@settings(max_examples=1000)
@given(integers(), integers(min_value=1))
@example(0, 1)
@example(1, 1)
@example(-1, 2)
@example(123456789123456789, 9876)
def test_mpq_to_from_binary(p, q):
    x = mpq(p,q)
    assert x == from_binary(to_binary(x))


@settings(max_examples=1000)
@given(fractions())
@example(Fraction(0))
@example(Fraction(-1))
@example(Fraction(-2))
@example(Fraction(15432, 125))
@example(Fraction(1, sys.hash_info.modulus))
@example(Fraction(-1, sys.hash_info.modulus))
def test_mpq_hash(q):
    assert hash(mpq(q)) == hash(q)


def test_mpq_digits():
    q = mpq(2/3)

    assert q.digits() == '6004799503160661/9007199254740992'
    assert q.digits(16) == '0x15555555555555/0x20000000000000'

    pytest.raises(TypeError, lambda: q.digits(16, 5))
    pytest.raises(ValueError, lambda: q.digits(0))


def test_mpq_abs():
    a = mpq(12,7)
    b = abs(a)

    assert a is b

    a = mpq(-12,7)
    b = abs(a)

    assert b == mpq(12,7)
    assert a == mpq(-12,7)


def test_mpq_add():
    a = mpq(3,11)

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


def test_mpq_sub():
    assert mpq(1,2) - Fraction(3,2) == mpq(-1,1)
    assert Fraction(1,2) - mpq(3,2) == mpq(-1,1)
    assert mpq(1,2) - mpq(3,2) == mpq(-1,1)
    assert mpq(1,2) - 0 == mpq(1,2)
    assert mpq(1,2) - mpz(1) == mpq(-1,2)
    assert mpq(1,2) + (-1) == mpq(-1,2)
    assert 1 - mpq(1,2) == mpq(1,2)
    assert mpz(1) - mpq(1,2) == mpq(1,2)
    assert mpq(1,2) - mpz(1) == mpq(-1,2)
    assert mpq(1,1) - mpc(0) == mpc('1.0+0.0j')
    assert mpc(0) - mpq(1,1) == mpc('-1.0+0.0j')
    assert mpq(1,2) - z == mpq(-3,2)
    assert mpq(1,2) - q ==mpq(-1,1)

    ctx = gmpy2.context()

    assert ctx.sub(mpq(1,2), mpq(3,2)) == mpq(-1,1)
    assert ctx.sub(mpq(1,2), Fraction(3,2)) == mpq(-1,1)
    assert ctx.sub(Fraction(1,2), mpq(3,2)) == mpq(-1,1)
    assert ctx.sub(Fraction(1,2), Fraction(3,2)) == mpq(-1,1)

    pytest.raises(TypeError, lambda: ctx.sub(1,'a'))
    pytest.raises(TypeError, lambda: mpq(1,2) - 'a')
    pytest.raises(TypeError, lambda: 'a' - mpq(1,2))

    a = mpq(3,11)
    b = mpq(1,2)
    c = Fraction(5,7)

    assert a-b == mpq(-5,22)
    assert b-a == mpq(5,22)
    assert a-1 == mpq(-8,11)
    assert 1-a == mpq(8,11)
    assert a-c == mpq(-34,77)
    assert c-a == mpq(34,77)
    assert a-a == mpq(0,1)
    assert a-mpz(123456) == mpq(-1358013,11)
    assert mpz(-123456)-a == mpq(-1358019,11)

    pytest.raises(TypeError, lambda: a-'b')
    pytest.raises(TypeError, lambda: 'b'-a)

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


def test_mpq_mul():
    assert mpq(1,2) * Fraction(3,2) == mpq(3,4)
    assert Fraction(1,2) * mpq(3,2) == mpq(3,4)
    assert mpq(1,2) * mpq(3,2) == mpq(3,4)
    assert mpq(1,2) * 0 == mpq(0,1)
    assert mpq(1,2) * mpz(1) == mpq(1,2)
    assert mpq(1,2) * (-1) == mpq(-1,2)
    assert mpq(1,1) * mpc(1,0) == mpc('1.0+0.0j')
    assert mpc(1,0) * mpq(1,1) == mpc('1.0+0.0j')
    assert mpq(1,2) * z == mpq(1,1)
    assert mpq(1,2) * q == mpq(3,4)

    ctx = gmpy2.context()

    assert ctx.mul(mpq(1,2), mpq(3,2)) == mpq(3,4)
    assert ctx.mul(mpq(1,2), Fraction(3,2)) == mpq(3,4)
    assert ctx.mul(Fraction(1,2), mpq(3,2)) == mpq(3,4)
    assert ctx.mul(Fraction(1,2), Fraction(3,2)) == mpq(3,4)

    pytest.raises(TypeError, lambda: ctx.mul(1,'a'))
    pytest.raises(TypeError, lambda: mpq(1,2) * 'a')
    pytest.raises(TypeError, lambda: 'a' * mpq(1,2))

    a = mpq(3,11)
    b = mpq(1,2)

    assert a*b == mpq(3,22)
    assert b*a == mpq(3,22)
    assert a*0 == mpq(0,1)
    assert 0*a == mpq(0,1)
    assert a*-1 == mpq(-3,11)
    assert -1*a == mpq(-3,11)
    assert a*mpz(17) == mpq(51,11)
    assert mpz(17)*a == mpq(51,11)
    assert a*a == mpq(9,121)

    pytest.raises(TypeError, lambda: a*'b')
    pytest.raises(TypeError, lambda: 'b'*a)

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
    assert is_nan(mpz(0) * float('Inf'))
    assert is_nan(mpz(0) * float('-Inf'))
    assert is_nan(float('Inf') * mpz(0))
    assert is_nan(float('-Inf') * mpz(0))

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
    assert is_nan(mpfr('Inf') * mpz(0))
    assert is_nan(mpfr('-Inf') * mpz(0))


def test_mpq_mod():
    a = mpq(3,11)
    b = mpq(1,2)
    ctx = gmpy2.get_context()

    assert a%b == mpq(3,11)
    assert b%a == mpq(5,22)
    assert a%z == mpq(3,11)
    assert mpq(3,1) % q == mpq(0,1)

    pytest.raises(ZeroDivisionError, lambda: a % mpq(0,5))

    assert ctx.mod(a, mpfr(1.5)) == mpfr('0.27272727272727271')

    pytest.raises(TypeError, lambda: ctx.mod(a, mpc(0.5, 2)))
    pytest.raises(ZeroDivisionError, lambda: ctx.mod(a, mpq(0, 2)))

    assert a % mpfr(1.5) == mpfr('0.27272727272727271')

    pytest.raises(TypeError, lambda: a % mpc(0.5, 2))
    pytest.raises(TypeError, lambda: a % 'str')
    pytest.raises(TypeError, lambda: ctx.mod(a, 'str'))


def test_mpq_divmod():
    a = mpq(3,11)
    b = mpq(1,2)
    c = Fraction(5,7)
    ctx = gmpy2.get_context()

    assert divmod(a,b) == (mpz(0), mpq(3,11))
    assert divmod(b,a) == (mpz(1), mpq(5,22))

    pytest.raises(ZeroDivisionError, lambda: divmod(a,0))

    assert divmod(17,a) == (mpz(62), mpq(1,11))
    assert divmod(mpz(17),a) == (mpz(62), mpq(1,11))
    assert divmod(a,c) == (mpz(0), mpq(3,11))

    pytest.raises(TypeError, lambda: divmod(mpq(1,2),'a'))

    ctx = gmpy2.ieee(64)
    gmpy2.set_context(ctx)

    assert ctx.divmod(mpq(3,2),mpq(3,7)) == (mpz(3), mpq(3,14))

    pytest.raises(TypeError, lambda: divmod(mpq(1,2), mpc(1,2)))

    assert divmod(a, float('Inf')) == (mpfr('0.0'), mpfr('0.27272727272727271'))
    assert divmod(a, float('-Inf')) == (mpfr('-1.0'), mpfr('-inf'))
    assert divmod(-a, float('Inf')) == (mpfr('-1.0'), mpfr('inf'))
    assert divmod(-a, float('-Inf')) == (mpfr('0.0'), mpfr('-0.27272727272727271'))
    assert all(map(is_nan, divmod(a, float('nan'))))
    assert all(map(is_nan, divmod(-a, float('nan'))))
    assert divmod(mpz(0), float('Inf')) == (mpfr('0.0'), mpfr('0.0'))
    assert divmod(mpz(0), float('-Inf')) == (mpfr('-0.0'), mpfr('-0.0'))
    assert all(map(is_nan, divmod(mpz(0), float('nan'))))
    assert all(map(is_nan, divmod(float('Inf'), a)))
    assert all(map(is_nan, divmod(float('-Inf'), a)))
    assert all(map(is_nan, divmod(float('Inf'), -a)))
    assert all(map(is_nan, divmod(float('-Inf'), -a)))
    assert all(map(is_nan, divmod(float('nan'), a)))
    assert all(map(is_nan, divmod(float('nan'), -a)))
    assert all(map(is_nan, divmod(float('Inf'), mpz(0))))
    assert all(map(is_nan, divmod(float('-Inf'), mpz(0))))
    assert all(map(is_nan, divmod(float('nan'), mpz(0))))

    assert divmod(a, mpfr('Inf')) == (mpfr('0.0'), mpfr('0.27272727272727271'))
    assert divmod(a, mpfr('-Inf')) == (mpfr('-1.0'), mpfr('-inf'))
    assert divmod(-a, mpfr('Inf')) == (mpfr('-1.0'), mpfr('inf'))
    assert divmod(-a, mpfr('-Inf')) == (mpfr('0.0'), mpfr('-0.27272727272727271'))
    assert all(map(is_nan, divmod(a, mpfr('nan'))))
    assert all(map(is_nan, divmod(-a, mpfr('nan'))))
    assert divmod(mpz(0), mpfr('Inf')) == (mpfr('0.0'), mpfr('0.0'))
    assert divmod(mpz(0), mpfr('-Inf')) == (mpfr('-0.0'), mpfr('-0.0'))
    assert divmod(mpz(0), mpfr('nan'))
    assert divmod(mpfr('Inf'), a)
    assert divmod(mpfr('-Inf'), a)
    assert divmod(mpfr('Inf'), -a)
    assert divmod(mpfr('-Inf'), -a)
    assert divmod(mpfr('nan'), a)
    assert divmod(mpfr('nan'), -a)
    assert divmod(mpfr('Inf'), mpz(0))
    assert divmod(mpfr('-Inf'), mpz(0))
    assert divmod(mpfr('nan'), mpz(0))


def test_mpq_floordiv():
    ctx = gmpy2.get_context()
    a, b = mpz(45), mpz(6)
    r, r2 = mpfr(45), mpfr(3.1)
    q, q2 = mpq(118,18), mpq(3,2)
    pyq, pyq2 = Fraction(118,18), Fraction(3,2)
    c, c2 = mpc(51, 65), mpc(4, 6)

    assert ctx.floor_div(q, q2) == mpz(4)
    assert ctx.floor_div(q, pyq2) == mpz(4)
    assert ctx.floor_div(q, pyq) == mpz(1)
    assert ctx.floor_div(pyq, q2) == mpz(4)
    assert ctx.floor_div(pyq, pyq2) == mpz(4)
    assert ctx.floor_div(pyq, q) == mpz(1)

    pytest.raises(ZeroDivisionError, lambda: ctx.floor_div(q, mpq(0,1)))
    pytest.raises(ZeroDivisionError, lambda: ctx.floor_div(q, Fraction(0)))
    pytest.raises(ZeroDivisionError, lambda: ctx.floor_div(pyq, mpq(0,1)))
    pytest.raises(ZeroDivisionError, lambda: ctx.floor_div(pyq, Fraction(0)))

    assert q // b == mpz(1)
    assert q // q2 == mpz(4)
    assert q // r2 == mpfr('2.0')

    pytest.raises(TypeError, lambda: q // c2)
    pytest.raises(TypeError, lambda: q // 'not')

    a = mpq(3,11)
    b = mpq(1,2)

    assert a//b == mpz(0)
    assert b//a == mpz(1)

    pytest.raises(ZeroDivisionError, lambda: a//0)

    assert 0//a == mpz(0)
    assert mpq(355, 113) // 2 == mpz(1)
    assert mpq(355, 113) // mpz(2) == mpz(1)
    assert mpq(355, 113) // z == mpz(1)


def test_mpq_truediv():
    a = mpq(3,11)
    b = mpq(1,2)
    ctx = gmpy2.get_context()

    assert a/b == mpq(6,11)
    assert gmpy2.div(a, b) == mpq(6,11)
    assert ctx.div(a, b) == mpq(6,11)
    assert b/a == mpq(11,6)

    pytest.raises(ZeroDivisionError, lambda: a/0)

    assert 0/a == mpq(0,1)
    assert a / z == mpq(3,22)
    assert mpq(3,11) / q == mpq(2,11)
    assert a / mpc(5, 4) == mpc('0.03325942350332594-0.02660753880266075j')

    pytest.raises(ZeroDivisionError, lambda: a / mpq(0,1))
    pytest.raises(TypeError, lambda: a / 'str')

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
    assert is_nan(float('nan') / mpz(0))

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
    assert is_nan(mpfr('nan') / mpz(0))
    assert is_nan(mpfr('nan') / mpz(0))


def test_mpq_pow():
    q = mpq(2,3)
    ctx = gmpy2.get_context()

    assert q ** 2 == mpq(4,9)
    assert q ** 0 == mpq(1,1)
    assert q ** -5 == mpq(243,32)
    assert ctx.pow(Fraction(2,3),2) == q ** 2
    assert mpq(-5,8) ** 5 == mpq(-3125,32768)
    assert q ** mpq(4,5) == mpfr('0.72298118079846574')

    pytest.raises(TypeError, lambda: pow(q, 5, 2))


def test_mpq_attributes():
    q = mpq('4/5')
    a = mpq(3,11)
    pyq = Fraction(4, 5)

    assert q.numerator == mpz(4)
    assert q.denominator == mpz(5)
    assert gmpy2.numer(q) == mpz(4)
    assert gmpy2.denom(q) == mpz(5)

    assert a.numerator == mpz(3)
    assert a.denominator == mpz(11)
    assert a.real == mpq(3,11)
    assert a.imag == mpz(0)

    pytest.raises(TypeError, lambda: gmpy2.numer(6.2))
    pytest.raises(TypeError, lambda: gmpy2.denom(5.6))
    pytest.raises(TypeError, lambda: gmpy2.denom(mpfr(5)))

    assert gmpy2.numer(pyq) == mpz(4)
    assert gmpy2.denom(pyq) == mpz(5)


def test_mpq_is_integer():
    assert mpq(0).is_integer()
    assert mpq(123).is_integer()
    assert not mpq(1, 2).is_integer()


def test_mpq_conjugate():
    a = mpq(3, 11)

    assert a.conjugate() == a


def test_mpq_floor():
    assert math.floor(mpq('7/2')) == mpz(3)


def test_mpq_ceil():
    assert math.ceil(mpq('7/2')) == mpz(4)


def test_mpq_trunc():
    assert math.trunc(mpq('7/2')) == mpz(3)


def test_mpq_not():
    q = mpq('4/5')

    assert q
    assert not mpq('0/5')


def test_mpq_qdiv():
    q = mpq('4/5')
    pyq = Fraction(4, 5)

    assert gmpy2.qdiv(q) == mpq(4,5)
    assert gmpy2.qdiv(pyq) == mpq(4,5)
    assert gmpy2.qdiv(5) == mpz(5)

    pytest.raises(TypeError, lambda: gmpy2.qdiv(mpc(4, 5)))
    pytest.raises(TypeError, lambda: gmpy2.qdiv(4, 5, 4))
    pytest.raises(TypeError, lambda: gmpy2.qdiv(4, 5.6))

    assert gmpy2.qdiv(q, 2) == mpq(2,5)
    assert gmpy2.qdiv(10, q) == mpq(25,2)
    assert gmpy2.qdiv(1) == mpz(1)

    assert gmpy2.qdiv(8,1) == mpz(8)
    assert gmpy2.qdiv(8,2) == mpz(4)
    assert gmpy2.qdiv(8,3) == mpq(8,3)
    assert gmpy2.qdiv(8,4) == mpz(2)
    assert gmpy2.qdiv(mpq(3,4), mpq(1,3)) == mpq(9,4)
    assert gmpy2.qdiv(mpq(3,4), mpq(1,4)) == mpz(3)

    args = mpq(2), 1/gmpy2.mpq(2)

    assert gmpy2.qdiv(*args) == mpz(4)
    assert args == (mpq(2), 1/mpq(2))


def test_mpq_repr():
    assert repr(mpq(11,13)) == 'mpq(11,13)'


def test_mpq_str():
    assert str(mpq(11,13)) == '11/13'


def test_issue_334():
    x = mpq(3,2)
    y = mpq(x,2)

    assert x == mpq(3,2)
    assert y == mpq(3,4)
    assert id(x) is not id(y)


def test_mpq_limit_denominator():
    x = mpq(1, 2)

    pytest.raises(TypeError, lambda: x.limit_denominator(1, 2))
    pytest.raises(TypeError, lambda: x.limit_denominator(x=1, y=2))
    pytest.raises(TypeError, lambda: x.limit_denominator(x=1))
    pytest.raises(TypeError, lambda: x.limit_denominator(1, max_denominator=2))
    pytest.raises(TypeError, lambda: x.limit_denominator(1.1))
    pytest.raises(ValueError, lambda: x.limit_denominator(-10))

    x = mpq('3.141592653589793')

    assert x.limit_denominator(10) == mpq(22, 7)
    assert x.limit_denominator(100) == mpq(311, 99)
    assert mpq(4321, 8765).limit_denominator(10000) == mpq(4321, 8765)

    x = mpq('3.1415926535897932')
    assert x.limit_denominator(10000) == mpq(355, 113)
    assert -x.limit_denominator(10000) == mpq(-355, 113)
    assert x.limit_denominator(113) == mpq(355, 113)
    assert x.limit_denominator(112) == mpq(333, 106)
    assert mpq(201, 200).limit_denominator(100) == mpq(1)
    assert mpq(201, 200).limit_denominator(101) == mpq(102, 101)
    assert mpq(0).limit_denominator(10000) == mpq(0)
