import inspect
import sys
from ctypes import memmove

import pytest
from hypothesis import example, given, settings
from hypothesis.strategies import integers

import gmpy2
from gmpy2 import from_binary, mpfr, mpq, mpz, to_binary, xmpz


def test_xmpz_digits():
    x = xmpz(6)
    assert x.digits() == '6'
    assert x.digits(2) == '0b110'
    assert x.digits(8) == '0o6'
    assert x.digits(16) == '0x6'
    assert xmpz(30).digits(16) == '0x1e'

    pytest.raises(ValueError, lambda: x.digits(0))
    pytest.raises(ValueError, lambda: x.digits(63))


def test_xmpz_limbs():
    x, y = xmpz(123456789), xmpz(0);
    x_limbs, num_limbs = x.limbs_read(), x.num_limbs()

    assert x_limbs != 0
    assert num_limbs > 0

    y_limbs = y.limbs_write(num_limbs)

    assert y_limbs != 0
    assert memmove(y_limbs, x_limbs, num_limbs * xmpz.limb_size) > 0

    y.limbs_finish(num_limbs)

    assert int(y) == 123456789

    x, y = xmpz(987654321), xmpz(0)
    x_limbs, num_limbs = x.limbs_read(), x.num_limbs()

    assert x_limbs != 0
    assert num_limbs > 0

    y_limbs = y.limbs_modify(num_limbs);

    assert y_limbs != 0
    assert memmove(y_limbs, x_limbs, num_limbs * xmpz.limb_size) > 0

    y.limbs_finish(num_limbs)

    assert int(y) == 987654321


def test_xmpz_attributes():
    x = xmpz(10)

    assert x.numerator == xmpz(10)
    assert x.denominator == xmpz(1)
    assert x.real == xmpz(10)


def test_xmpz_misc():
    z = xmpz(5)

    assert gmpy2.xbit_mask(7) == xmpz(127)
    assert gmpy2.xbit_mask(z) == xmpz(31)

    pytest.raises(TypeError, lambda: gmpy2.xbit_mask(4.5))

    assert not xmpz(0)
    assert z
    -z
    assert z == xmpz(-5)
    -z
    assert z == xmpz(5)
    +z
    assert z == xmpz(5)

    z2 = z.copy()

    assert z2 == z
    assert z2 is not z
    assert z.copy().make_mpz() == mpz(5)
    assert xmpz(-100).make_mpz() == mpz(-100)
    assert len(z) == 3
    assert len(xmpz(3000)) == 12

    assert str(xmpz(42)) == '42'
    assert repr(xmpz(42)) == 'xmpz(42)'


def test_xmpz_create():
    pytest.raises(TypeError, lambda: xmpz(s=1))
    pytest.raises(TypeError, lambda: xmpz(1, s=2))


def test_xmpz_subscripts():
    x = xmpz(10)

    assert x[0] == x[2] == 0
    assert x[1] == x[3] == 1

    x[0] = 1

    assert x == xmpz(11)
    assert x[0:2] == mpz(3)
    assert x[0:] == mpz(11)
    assert x[:4] == mpz(11)

    x[:4] = 14

    assert x == xmpz(14)
    assert x[0] == 0

    pytest.raises(TypeError, lambda: x[mpfr('inf')])

    assert x[3:0] == mpz(0)
    assert x[-1] == 1

    x[-1] = 0

    assert x == xmpz(6)

    x[0:4] = 15

    assert x == xmpz(15)

    with pytest.raises(ValueError):
        x[0] = 16

    x[0:4] = 16

    assert x == xmpz(0)

    x[0:4] = 15

    assert x == xmpz(15)

    x[0:5] = 16

    assert x == xmpz(16)


def test_xmpz_iterators():
    x = xmpz(16)

    assert [b for b in x.iter_bits()] == [False, False, False, False, True]

    x = xmpz(30)

    assert [b for b in x.iter_bits()] == [False, True, True, True, True]
    assert [b for b in x.iter_bits(-1, 2)] == []
    assert [b for b in x.iter_bits(1, 3)] == [True, True]
    assert [b for b in x.iter_set()] == [1, 2, 3, 4]
    assert [b for b in x.iter_clear()] == [0]

    x = xmpz(10)

    assert [b for b in x.iter_clear()] == [0, 2]


def test_xmpz_conversion():
    assert xmpz('5') == xmpz(5)
    assert xmpz('5') == xmpz(5)

    pytest.raises(ValueError, lambda: xmpz('not'))

    assert xmpz(-3.5) == xmpz(-3)

    pytest.raises(OverflowError, lambda: xmpz(float('inf')))

    assert xmpz(mpz(100)) == xmpz(100)
    assert xmpz(xmpz(100)) == xmpz(100)
    assert xmpz(mpq(30,2)) == xmpz(15)
    assert str(xmpz(100)) == '100'


def test_xmpz_abs():
    a = xmpz(123)
    b = abs(a)

    assert a is not b

    a = xmpz(-123)
    b = abs(a)

    assert a == xmpz(123)
    assert b is None

    b = abs(a)

    assert b is None


def test_xmpz_iadd():
    x = xmpz(5)
    x += mpz(6)

    assert x == xmpz(11)

    x += 7

    assert x == xmpz(18)

    x += -6

    assert x == xmpz(12)

    x += mpfr(2.5)

    assert x == mpfr('14.5')


def test_xmpz_isub():
    x = xmpz(7)
    x -= xmpz(1)

    assert x == xmpz(6)

    x -= 1

    assert x == xmpz(5)

    x -= mpz(7)

    assert x == xmpz(-2)

    x -= -5

    assert x == xmpz(3)

    x -= -mpfr(5)

    assert x == mpfr('8.0')


def test_xmpz_imul():
    x = xmpz(2)
    x *= xmpz(2)

    assert x == xmpz(4)

    x *= 2

    assert x == xmpz(8)

    x *= mpz(3)

    assert x == xmpz(24)

    x *= -1

    assert x == xmpz(-24)

    x *= mpfr(-0.5)

    assert x == mpfr('12.0')

def test_xmpz_ifloordiv():
    x = xmpz(49)
    x //= xmpz(3)

    assert x == xmpz(16)

    x //= mpz(3)

    assert x == xmpz(5)

    x //= 2

    assert x == xmpz(2)

    with pytest.raises(ZeroDivisionError):
        x //= 0

    x == xmpz(2)
    x //= mpfr(-0.5)

    assert x == mpfr('-4.0')


def test_xmpz_imod():
    x = xmpz(45)
    x %= xmpz(18)

    assert x == xmpz(9)

    x %= mpz(2)

    assert x == xmpz(1)

    x = xmpz(40)
    x %= 21

    assert x == xmpz(19)

    with pytest.raises(ZeroDivisionError):
        x %= 0

    x == xmpz(19)
    x %= mpfr(10)

    assert x == mpfr('9.0')


def test_xmpz_ishifts():
    x = xmpz(63)
    x >>= xmpz(63)

    assert x == xmpz(0)

    x = xmpz(63)
    x >>= xmpz(1)

    assert x == xmpz(31)

    x >>= mpz(2)

    assert x == xmpz(7)

    x >>= 1

    assert x == xmpz(3)

    x <<= xmpz(2)

    assert x == xmpz(12)

    x <<= mpz(1)

    assert x == xmpz(24)

    x <<= 0

    assert x == xmpz(24)

    with pytest.raises(TypeError):
        x >>= mpfr(2)
    with pytest.raises(TypeError):
        x <<= mpfr(2)
    with pytest.raises(OverflowError):
        x >>= -1
    with pytest.raises(OverflowError):
        x <<= -5


def test_xmpz_ipow():
    x = xmpz(5)
    x **= xmpz(2)

    assert x == xmpz(25)

    x **= mpz(2)

    assert x == xmpz(625)

    with pytest.raises(OverflowError):
        x **= -2

    x **= 2

    assert x == xmpz(390625)

    with pytest.raises(TypeError):
        x **= mpfr(2)


def test_xmpz_iand():
    x = xmpz(7)
    x &= xmpz(5)

    assert x == xmpz(5)

    x &= mpz(4)

    assert x == xmpz(4)

    x &= 9

    assert x == xmpz(0)

    x = mpz(4)
    x &= 12

    assert x == mpz(4)

    with pytest.raises(TypeError):
        x &= mpfr(4)


def test_xmpz_ior():
    x = xmpz(0)
    x |= xmpz(1)

    assert x == xmpz(1)

    x |= xmpz(0)

    assert x == xmpz(1)

    x = xmpz(0)
    x |= xmpz(0)

    assert x == xmpz(0)

    x |= 5

    assert x == xmpz(5)

    with pytest.raises(TypeError):
        x |= mpfr(3)

def test_xmpz_ixor():
    x = xmpz(1)
    x ^= xmpz(0)

    assert x == xmpz(1)

    x ^= mpz(1)

    assert x == xmpz(0)

    x ^= 1

    assert x == xmpz(1)

    with pytest.raises(TypeError):
        x ^= mpfr(0)


@settings(max_examples=1000)
@example(0)
@example(1)
@example(-1)
@example(1234567890123456789)
@given(integers())
def test_xmpz_to_from_binary(x):
    x = xmpz(x)

    assert from_binary(to_binary(x)) == x


@pytest.mark.skipif(sys.version_info < (3, 13), reason="requires v3.13+")
def test_xmpz_signatures():
    cls = gmpy2.xmpz
    for f in dir(cls):
        a = getattr(cls, f)
        if callable(a) and f != '__class__':
            _ = inspect.signature(a)  # not raises
