import inspect
import math
import numbers
import pickle
import sys
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction

import pytest
from hypothesis import assume, example, given, settings
from hypothesis.strategies import booleans, integers, sampled_from
from pytest import mark, raises
from supportclasses import a, b, c, d, q, z

import gmpy2
from gmpy2 import (cmp, cmp_abs, from_binary, is_nan, is_prime, mp_version,
                   mpc, mpfr, mpq, mpz, mpz_random, mpz_rrandomb, mpz_urandomb,
                   next_prime, pack, random_state, to_binary, unpack, xmpz)


def test_mpz_to_bytes_interface():
    x = mpz(1)
    with raises(TypeError):
        x.to_bytes(1, 2, 3)
    with raises(TypeError):
        x.to_bytes(1, 2)
    with raises(TypeError):
        x.to_bytes('spam')
    with raises(TypeError):
        x.to_bytes(a=1, b=2, c=3, d=4)
    with raises(TypeError):
        x.to_bytes(2, length=2)
    with raises(TypeError):
        x.to_bytes(2, 'big', byteorder='big')
    with raises(TypeError):
        x.to_bytes(spam=1)

    with raises(ValueError):
        x.to_bytes(2, 'spam')
    with raises(ValueError):
        x.to_bytes(-1)

    assert x.to_bytes(2) == x.to_bytes(length=2)
    assert x.to_bytes(2, byteorder='little') == x.to_bytes(2, 'little')

    assert x.to_bytes() == x.to_bytes(1)
    assert x.to_bytes() == x.to_bytes(1, 'big')
    assert x.to_bytes() == x.to_bytes(signed=False)


@settings(max_examples=1000)
@given(integers(), integers(min_value=0, max_value=10000),
       sampled_from(['big', 'little']), booleans())
@example(0, 0, 'big', False)
@example(0, 0, 'little', False)
@example(0, 1, 'big', False)
@example(128, 1, 'big', True)
@example(128, 1, 'little', True)
@example(-129, 1, 'big', True)
@example(-129, 1, 'little', True)
@example(-1, 0, 'big', True)
@example(-2, 0, 'big', True)
@example(-2, 0, 'little', True)
@example(42, 1, 'big', False)
@example(42, 1, 'little', False)
@example(42, 3, 'big', False)
@example(42, 3, 'little', False)
@example(1000, 2, 'big', False)
@example(1000, 4, 'big', False)
@example(-2049, 1, 'big', True)
@example(-65281, 3, 'big', True)
@example(-65281, 3, 'little', True)
def test_mpz_to_bytes(x, length, byteorder, signed):
    try:
        rx = x.to_bytes(length, byteorder, signed=signed)
    except OverflowError:
        with raises(OverflowError):
            mpz(x).to_bytes(length, byteorder, signed=signed)
    else:
        assert rx == mpz(x).to_bytes(length, byteorder, signed=signed)


def test_mpz_from_bytes_interface():
    with raises(TypeError):
        mpz.from_bytes()
    with raises(TypeError):
        mpz.from_bytes(1, 2, 3)
    with raises(TypeError):
        mpz.from_bytes(b'', 2)
    with raises(TypeError):
        mpz.from_bytes(1)
    with raises(TypeError):
        mpz.from_bytes(b'', bytes=b'')
    with raises(TypeError):
        mpz.from_bytes(b'', 'big', byteorder='big')
    with raises(TypeError):
        mpz.from_bytes(b'', spam=1)
    with raises(TypeError):
        mpz.from_bytes(a=1, b=2, c=3, d=4)

    with raises(ValueError):
        mpz.from_bytes(b'', 'spam')

    assert mpz.from_bytes(b'\x01', byteorder='little') == mpz.from_bytes(b'\x01', 'little')

    assert mpz.from_bytes(b'\x01') == mpz.from_bytes(bytes=b'\x01')
    assert mpz.from_bytes(b'\x01') == mpz.from_bytes(b'\x01', 'big')
    assert mpz.from_bytes(b'\x01') == mpz.from_bytes(b'\x01', signed=False)


@settings(max_examples=1000)
@given(integers(), integers(min_value=0, max_value=10000),
       sampled_from(['big', 'little']), booleans())
@example(0, 0, 'big', False)
@example(0, 0, 'little', False)
@example(0, 1, 'big', False)
@example(128, 1, 'big', True)
@example(128, 1, 'little', True)
@example(-129, 1, 'big', True)
@example(-129, 1, 'little', True)
@example(-1, 0, 'big', True)
@example(-1, 1, 'big', True)
@example(-1, 1, 'little', True)
@example(-2, 0, 'big', True)
@example(-2, 0, 'little', True)
@example(-1, 3, 'big', True)
@example(-2, 3, 'big', True)
@example(-2, 5, 'little', True)
def test_mpz_from_bytes(x, length, byteorder, signed):
    try:
        bytes = x.to_bytes(length, byteorder, signed=signed)
    except OverflowError:
        assume(False)
    else:
        rx = int.from_bytes(bytes, byteorder, signed=signed)
        assert rx == mpz.from_bytes(bytes, byteorder, signed=signed)
        assert rx == mpz.from_bytes(bytearray(bytes), byteorder, signed=signed)
        assert rx == mpz.from_bytes(list(bytes), byteorder, signed=signed)


def test_mpz_as_integer_ratio():
    assert mpz(3).as_integer_ratio() == (mpz(3), mpz(1))


def test_mpz_numbers_abc():
    assert isinstance(mpz(2), numbers.Integral)


def test_mpz_pickling():
    for proto in range(pickle.HIGHEST_PROTOCOL + 1):
        for x in [mpz(12346789), mpz(-12346789), mpz(0)]:
            assert pickle.loads(pickle.dumps(x, protocol=proto)) == x


@settings(max_examples=1000)
@given(integers(), integers())
@example(0, 0)
@example(0, 1)
@example(-11, 75)
@example(14, 105)
@example(64, 123456789012345678901234567890)
@example(321, -123)
def test_mpz_arithmetics(i, z):
    i, z = map(mpz, [i, z])

    assert int(i) + int(z) == i + z
    assert int(z) + int(i) == z + i

    assert int(i) - int(z) == i - z
    assert int(z) - int(i) == z - i

    assert int(i) * int(z) == i * z
    assert int(z) * int(i) == z * i

    # Test all permutations of floor division
    if z:
        assert int(i) // int(z) == i // z
        assert int(i) % int(z) == i % z
        assert i % int(z) == i % z
        assert int(i) % z == i % z
        assert divmod(int(i), int(z)) == divmod(i, z)
        assert divmod(i, int(z)) == divmod(i, z)
        assert divmod(int(i), z) == divmod(i, z)

    if i:
        assert int(z) // int(i) == z // i
        assert int(z) % int(i) == z % i
        assert z % int(i) == z % i
        assert int(z) % i == z % i
        assert divmod(int(z), int(i)) == divmod(z, i)
        assert divmod(z, int(i)) == divmod(z, i)
        assert divmod(int(z), i) == divmod(z, i)


@settings(max_examples=1000)
@given(integers(min_value=0),
       integers(min_value=1, max_value=100000))
def test_mpz_pack_unpack_bulk(x, n):
    lst = unpack(x, n)
    assert pack(lst, n) == x


def test_mpz_pack_unpack():
    x = mpz(0)
    assert all((x == pack(unpack(x,i),i) for i in range(1,100)))
    x = mpz(1)
    assert all((x == pack(unpack(x,i),i) for i in range(1,100)))
    x = mpz(2)
    assert all((x == pack(unpack(x,i),i) for i in range(1,100)))
    x = mpz(3141592635)
    assert all((x == pack(unpack(x,i),i) for i in range(1,100)))
    x = mpz(1234567891234567890000000000000000000000000000000000000123)
    assert all((x == pack(unpack(x,i),i) for i in range(1,100)))
    x = mpz(1) << 500
    assert all((x == pack(unpack(x,i),i) for i in range(1,200)))
    x -= 1
    assert all((x == pack(unpack(x,i),i) for i in range(1,200)))

    raises(TypeError, lambda: pack(x))
    raises(TypeError, lambda: pack(1, 1))
    raises(TypeError, lambda: pack([mpz(1), mpz(-1)], 2))

    raises(TypeError, lambda: unpack(x))
    raises(TypeError, lambda: unpack([], 1))
    raises(ValueError, lambda: unpack(-1, 1))


def test_mpz_cmp():
    assert cmp(0, mpz(0)) == 0
    assert cmp(1, mpz(0)) == 1
    assert cmp(0, mpz(1)) == -1
    assert cmp(-1, mpz(0)) == -1
    assert cmp(0, mpz(-1)) == 1

    assert cmp_abs(mpz(0), 0) == 0
    assert cmp_abs(mpz(1), 0) == 1
    assert cmp_abs(mpz(0), 1) == -1
    assert cmp_abs(mpz(-1), 0) == 1
    assert cmp_abs(mpz(0), -1) == -1

    a = mpz(-10)
    assert cmp_abs(a, 0) == 1
    assert a == mpz(-10)
    assert cmp_abs(100, a) == 1
    assert a == mpz(-10)

    assert cmp(mpz(2), z) == 0
    assert cmp(z, mpz(3)) == -1
    assert cmp(mpz(1), q) == -1
    assert cmp(mpz(1), mpz(q)) == 0


def test_mpz_comparisons():
    from supportclasses import q

    assert mpz(2) == z
    assert (z == mpz(3)) is False
    assert (mpz(1) == q) is False
    assert mpz(1) == mpz(q)

    a = mpz(123)
    b = mpz(456)
    c = mpz(a)
    q = mpq(4, 5)

    assert a == mpz(123)
    assert b == mpz(456)
    assert c is a
    assert c==a
    assert (c>a) is False
    assert (c<a) is False
    assert (a>b) is False
    assert a<b
    assert not mpz(0)
    assert not a is False
    assert (mpz(1) == None) is False
    assert (mpz(1) == '1') is False
    assert (mpz(1) == 'abc') is False
    assert [mpz(23), None].count(None) == 1
    assert (a == q, a != q, a > q, a >= q, a < q, a <= q) == (False, True, True, True, False, False)

    q = mpq(123, 1)

    assert (a == q, a != q, a > q, a >= q, a < q, a <= q) == (True, False, False, True, False, True)

    gmpy2.context().trap_divzero == False
    f = float('inf')

    assert (a == f, a != f, a > f, a >= f, a < f, a <= f) == (False, True, False, False, True, True)
    assert (f == a, f != a, f > a, f >= a, f < a, f <= a) == (False, True, True, True, False, False)

    f = float('-inf')

    assert (a == f, a != f, a > f, a >= f, a < f, a <= f) == (False, True, True, True, False, False)
    assert (f == a, f != a, f > a, f >= a, f < a, f <= a) == (False, True, False, False, True, True)

    f = float('nan')

    assert (a == f, a != f, a > f, a >= f, a < f, a <= f) == (False, True, False, False, False, False)
    assert (f == a, f != a, f > a, f >= a, f < a, f <= a) == (False, True, False, False, False, False)

    r = mpfr('inf')

    assert (a == r, a != r, a > r, a >= r, a < r, a <= r) == (False, True, False, False, True, True)

    r = mpfr('-inf')

    assert (a == r, a != r, a > r, a >= r, a < r, a <= r) == (False, True, True, True, False, False)

    r = mpfr('nan')

    assert (a == r, a != r, a > r, a >= r, a < r, a <= r) == (False, True, False, False, False, False)


def test_mpz_conversion():
    x = mpz(a)
    assert isinstance(x, mpz)
    assert x == 1
    raises(TypeError, lambda: mpz(b))
    raises(TypeError, lambda: mpz(c))
    raises(TypeError, lambda: mpz(d))

    assert float(mpz(1)) == 1.0
    raises(OverflowError, lambda: float(mpz(99**199)))
    assert mpz(xmpz(1)) == mpz(1)

    assert int(mpz(-3)) == -3

    assert int(mpz(11)) is int(mpz(11))


def test_mpz_create():
    assert mpz() == mpz(0)
    assert mpz(0) == mpz(0)
    assert mpz(1) == mpz(1)
    assert mpz(-1) == mpz(-1)
    assert mpz(2**15-2) == mpz(32766)
    assert mpz(2**15-1) == mpz(32767)
    assert mpz(2**15) == mpz(32768)
    assert mpz(2**15+1) == mpz(32769)
    assert mpz(2**15+2) == mpz(32770)
    assert mpz(2**30-2) == mpz(1073741822)
    assert mpz(2**30-1) == mpz(1073741823)
    assert mpz(2**30) == mpz(1073741824)
    assert mpz(2**30+1) == mpz(1073741825)
    assert mpz(2**30+2) == mpz(1073741826)
    assert mpz(2**16-2) == mpz(65534)
    assert mpz(2**16-1) == mpz(65535)
    assert mpz(2**16) == mpz(65536)
    assert mpz(2**16+1) == mpz(65537)
    assert mpz(2**16+2) == mpz(65538)
    assert mpz(1000000000000) == mpz(1000000000000)
    assert mpz(-1000000000000) == mpz(-1000000000000)

    raises(ValueError, lambda: mpz(''))
    raises(ValueError, lambda: mpz('a'))

    assert mpz('a',16) == mpz(10)

    raises(ValueError, lambda: mpz('z',16))

    assert mpz('0b1101') == mpz(13)
    assert mpz('0b1101',2) == mpz(13)
    assert mpz('1101',2) == mpz(13)
    assert mpz('0b0010') == mpz(2)
    assert mpz('0b0010',2) == mpz(2)

    raises(ValueError, lambda: mpz('0b0b10',2))
    raises(ValueError, lambda: mpz('0b0b10'))
    raises(ValueError, lambda: mpz('0b0012'))

    assert mpz('0o0012') == mpz(10)
    assert mpz('0o0012',8) == mpz(10)
    assert mpz('12',8) == mpz(10)
    assert mpz('0x12') == mpz(18)
    assert mpz('0x12',16) == mpz(18)
    assert mpz('12',16) == mpz(18)
    assert mpz('-1') == mpz(-1)
    assert mpz('+1') == mpz(1)
    assert mpz('  0xA', base=0) == mpz(10)

    raises(ValueError, lambda: mpz(float('nan')))
    raises(OverflowError, lambda: mpz(float('inf')))
    raises(OverflowError, lambda: mpz(float('-inf')))
    raises(TypeError, lambda: mpz(12, base=16))

    assert mpz('12', base=16) == mpz(18)

    raises(ValueError, lambda: mpz('\xff'))
    raises(ValueError, lambda: mpz('\x0cf'))
    raises(ValueError, lambda: mpz('\0xff'))

    assert mpz(b'12') == mpz(12)

    raises(TypeError, lambda: mpz(None))
    raises(TypeError, lambda: mpz(None,base=10))
    raises(ValueError, lambda: mpz('99',base=100))
    raises(TypeError, lambda: mpz('99',base='a'))

    assert mpz('99',base=10) == mpz(99)
    assert mpz(xmpz(5)) == mpz(5)

    raises(ValueError, lambda: mpz('ы'))
    raises(ValueError, lambda: mpz(bytes('ы', encoding='utf-8')))

    assert mpz(3.14) == mpz(3)
    assert mpz(mpq(17,3)) == mpz(5)
    assert mpz(23) == mpz(23)
    assert mpz(-23) == mpz(-23)

    x = 1000*1000*1000*1000*1000*1000*1000

    assert mpz(x) == 1000000000000000000000
    assert mpz(0.0) == mpz(0)
    assert mpz(-0.0) == mpz(0)

    raises(ValueError, lambda: mpz(float("nan")))
    raises(OverflowError, lambda: mpz(float("inf")))
    raises(OverflowError, lambda: mpz(float("-inf")))

    assert mpz("0") == mpz(0)
    assert mpz("-0") == mpz(0)

    raises(ValueError, lambda: mpz("hi"))

    assert mpz("123456", 7) == mpz(22875)

    raises(ValueError, lambda: mpz("123456", base=3))

    assert mpz() == mpz(0)
    assert mpz(Fraction(1,2)) == mpz(0)
    assert mpz(Fraction(-3,2)) == mpz(-1)
    assert mpz(Fraction(3,2)) == mpz(1)
    assert mpz('043') == mpz(43)
    assert mpz('43',0) == mpz(43)
    assert mpz('0o43') == mpz(35)
    assert mpz('0x43') == mpz(67)

    raises(ValueError, lambda: mpz('0x43',10))

    assert mpz('43') == mpz(43)
    assert mpz('1_2') == mpz(12)
    assert mpz('_1_2') == mpz(12)
    assert mpz('1 2') == mpz(12)
    assert mpz(' 1 2') == mpz(12)

    raises(TypeError, lambda: mpz(s=1))
    raises(TypeError, lambda: mpz(1, s=2))


@settings(max_examples=10000)
@given(integers())
@example(0)
@example(-3)
def test_mpz_conversion_bulk(n):
    m = mpz(n)
    assert int(m) == n
    assert str(m) == str(n)


@settings(max_examples=1000)
@given(integers())
@example(0)
@example(1)
@example(-1)
@example(123456789123456789)
def test_mpz_to_from_binary(n):
    x = mpz(n)
    assert x == from_binary(to_binary(x))


@settings(max_examples=1000)
@given(integers())
@example(0)
@example(1)
@example(-1)
@example(-2)
@example(123)
def test_mpz_hash(n):
    assert hash(mpz(n)) == hash(n)


def test_mpz_ceil():
    a = mpz(123)
    assert math.ceil(a) == a
    assert math.ceil(a) is a


def test_mpz_floor():
    a = mpz(123)
    assert math.floor(a) == a
    assert math.floor(a) is a


def test_mpz_trunc():
    a = mpz(123)
    assert math.trunc(a) == a
    assert math.trunc(a) is a


def test_mpz_round():
    assert round(mpz(123456), 2) == mpz(123456)
    assert round(mpz(123456), -22) == mpz(0)
    assert round(mpz(123456), -2) == mpz(123500)
    assert round(mpz(123456), -1) == mpz(123460)
    assert round(mpz(123455), -1) == mpz(123460)
    assert round(mpz(123454), -1) == mpz(123450)
    assert round(mpz(123445), -1) == mpz(123440)
    assert round(mpz(123445)) == mpz(123445)

    raises(TypeError, lambda: round(mpz(123456),'a'))
    raises(TypeError, lambda: round(mpz(123456),'a',4))

    # issue 552
    assert round(mpz(501), -3) == mpz(1000)

    # issue 567
    gmpy2.set_context(gmpy2.ieee(64))
    ctx = gmpy2.get_context()
    x = mpz(123456)

    assert round(x, -1) == mpz(123460)
    ctx.round = gmpy2.RoundUp
    assert round(x, -1) == mpz(123460)
    ctx.round = gmpy2.RoundAwayZero
    assert round(x, -1) == mpz(123460)
    ctx.round = gmpy2.RoundToZero
    assert round(x, -1) == mpz(123450)
    ctx.round = gmpy2.RoundDown
    assert round(x, -1) == mpz(123450)


@settings(max_examples=10000)
@given(integers())
@example(38732858750156021)
@example(225188150488381457)
def test_mpz_float_bulk(n):
    m = mpz(n)
    try:
        fn = float(n)
    except OverflowError:
        with raises(OverflowError):
            float(m)
    else:
        assert fn == float(m)


def test_mpz_bool():
    assert bool(mpz(100))
    assert not bool(mpz(0))
    assert bool(mpz(-100))


def test_mpz_random():
    r1 = random_state(42)
    r2 = random_state(42)

    assert mpz_random(r1, 2**88) == mpz(171378365038768291737094841)
    assert mpz_random(r2, 2**88) == mpz(171378365038768291737094841)
    assert mpz_random(r1, 2**88) == mpz(62749575961297098445301393)
    assert mpz_random(r2, 2**88) == mpz(62749575961297098445301393)


def test_mpz_urandomb():
    assert (mpz_urandomb(random_state(42), 64).digits(2) ==
            '1100100011011011101100101001100100111110010111011100101010111001')


def test_mpz_rrandomb():
    assert (mpz_rrandomb(random_state(42), 64).digits(2) ==
            '1111111111111111111111111100000000111111111111111111000000000000')

@mark.skipif(mp_version() < "GMP 6.3.0", reason="requires GMP 6.3.0 or higher")
def test_prev_prime():
    # Imported here as symbol won't exist if mp_version() < 6.3.0
    from gmpy2 import prev_prime
    assert prev_prime(3) == mpz(2)
    assert prev_prime(4) == mpz(3)
    assert prev_prime(5) == mpz(3)
    assert prev_prime(6) == mpz(5)
    assert prev_prime(1000004) == mpz(1000003)
    assert prev_prime(1000033) == mpz(1000003)

    with raises(ValueError):
        prev_prime(-100)

    with raises(ValueError):
        prev_prime(2)

    with raises(TypeError):
        prev_prime('a')

    with raises(TypeError):
        prev_prime(4.5)


def test_mpz_format():
    z1, z2 = mpz(-3), mpz(5)

    assert '{:<5}'.format(z1) == '-3   '
    assert '{:>+5}'.format(z2) == '   +5'

    raises(ValueError, lambda: '{:5+}'.format(z1))

    assert '{:>-4}'.format(z2) == '   5'
    assert '{:<-4}'.format(z1) == '-3  '

    raises(ValueError, lambda: '{:>4-}'.format(z1))
    raises(ValueError, lambda: '{:<4 }'.format(z1))

    assert '{:#x}'.format(z1) == '-0x3'
    assert '{:#o}'.format(z1) == '-0o3'

    raises(ValueError, lambda: '{:>5#}'.format(z1))
    raises(ValueError, lambda: '{:~}'.format(z1))

    a = mpz(123)

    assert str(a) == '123'
    assert repr(a) == 'mpz(123)'
    assert hex(a) == '0x7b'
    assert oct(a) == '0o173'
    assert mpz('1001001011',2) == mpz(587)
    assert bin(mpz('1001001011',2)) == '0b1001001011'
    assert '1001001011' == mpz('1001001011',2).digits(2)
    assert [a.digits(i) for i in range(2,63)] == ['1111011', '11120', '1323',
                                                  '443', '323', '234', '173',
                                                  '146', '123', '102', 'a3',
                                                  '96', '8b', '83', '7b', '74',
                                                  '6f', '69', '63', '5i', '5d',
                                                  '58', '53', '4n', '4j', '4f',
                                                  '4b', '47', '43', '3u', '3r',
                                                  '3o', '3l', '3i', '3f', '3C',
                                                  '39', '36', '33', '30', '2d',
                                                  '2b', '2Z', '2X', '2V', '2T',
                                                  '2R', '2P', '2N', '2L', '2J',
                                                  '2H', '2F', '2D', '2B', '29',
                                                  '27', '25', '23', '21', '1z']

    raises(ValueError, lambda: a.digits(63))

    assert '{}'.format(a) == '123'
    assert '{:d}'.format(a) == '123'
    assert '{:b}'.format(a) == '1111011'
    assert '{:o}'.format(a) == '173'
    assert '{:x}'.format(a) == '7b'
    assert '{:#x}'.format(a) == '0x7b'
    assert '{:#X}'.format(a) == '0X7B'
    assert '{:#o}'.format(a) == '0o173'
    assert '{:#15o}'.format(a) == '          0o173'
    assert '{:<#15o}'.format(a) == '0o173          '
    assert '{:^#15o}'.format(a) == '     0o173     '
    assert '{:>#15o}'.format(a) == '          0o173'
    assert '{:^ #15o}'.format(a) == '     0o173     '
    assert '{:^#15o}'.format(a) == '     0o173     '
    assert '{:^ #16o}'.format(a) == '      0o173     '

    raises(ValueError, lambda: '{:#^16o}'.format(a))

    assert '{:^#16o}'.format(a) == '     0o173      '

    # floating-point formats
    assert '{:f}'.format(a) == '123.000000'
    assert '{:g}'.format(a) == '123.0'

    gmpy2.set_context(gmpy2.ieee(64))
    assert '{:f}'.format(a<<1024) == 'inf'


def test_mpz_digits():
    z1, z2 = mpz(-3), mpz(15)

    assert z1.digits() == '-3'
    assert z1.digits(2) == '-11'
    assert z1.digits(8) == '-3'
    assert z2.digits(16) == 'f'

    raises(ValueError, lambda: z1.digits(0))
    raises(ValueError, lambda: z1.digits(1))


def test_mpz_abs():
    a = mpz(123)
    b = abs(a)

    assert a is b
    assert abs(-a) == a

    a = mpz(-123)
    b = abs(a)

    assert b == mpz(123)
    assert a is not b
    assert a == mpz(-123)


def test_mpz_add():
    a = mpz(123)
    b = mpz(456)
    c = 12345678901234567890

    assert a+1 == mpz(124)
    assert a+(-1) == mpz(122)
    assert 1+a == mpz(124)
    assert (-1)+a == mpz(122)
    assert a+b == mpz(579)
    assert b+a == mpz(579)

    raises(TypeError, lambda: a+'b')
    raises(TypeError, lambda: 'b'+a)
    assert a+c == 12345678901234568013
    assert c+a == 12345678901234568013

    assert a+True == mpz(124)
    assert a+False == mpz(123)

    assert a + float('Inf') == mpfr('inf')
    assert float('Inf') + a == mpfr('inf')
    assert a + float('-Inf') == mpfr('-inf')
    assert float('-Inf') + a == mpfr('-inf')
    assert is_nan(a + float('nan'))
    assert is_nan(float('nan') + a)


def test_mpz_iadd():
    x = mpz(5)
    x += mpz(6)

    assert x == mpz(11)

    x += 7

    assert x == mpz(18)

    x += -6

    assert x == mpz(12)

    x += 0

    assert x == mpz(12)

    x += mpfr(2.5)

    assert x == mpfr('14.5')


def test_mpz_sub():
    a, b = mpz(123), mpz(456)
    c = 12345678901234567890

    assert a-1 == mpz(122)
    assert a-(-1) == mpz(124)
    assert 1-a == mpz(-122)
    assert (-1)-a == mpz(-124)
    assert a-c == -12345678901234567767
    assert c-a == 12345678901234567767
    assert a-(-c) == 12345678901234568013
    assert (-c)-a == -12345678901234568013
    assert a-b == mpz(-333)
    assert b-a == mpz(333)
    assert a-(-b) == mpz(579)
    assert (-b)-a == mpz(-579)
    assert a-z == mpz(121)
    assert gmpy2.sub(2,1) == mpz(1)

    ctx = gmpy2.context()

    assert ctx.sub(a,b) == a-b
    assert ctx.sub(c,c) == c-c
    assert ctx.sub(1, 1) == mpz(0)
    assert ctx.sub(a, 1) == mpz(122)
    assert ctx.sub(1, a) == mpz(-122)
    assert ctx.sub(a, mpq(0)) == mpq(123,1)
    assert ctx.sub(a, mpfr(0)) == mpfr('123.0')
    assert ctx.sub(a, mpc(0)) == mpc('123.0+0.0j')

    raises(TypeError, lambda: ctx.sub(1))
    raises(TypeError, lambda: ctx.sub(1,2,3))
    raises(TypeError, lambda: a-'b')
    raises(TypeError, lambda: 'b'-a)

    assert a-1 == mpz(122)
    assert a-(-1) == mpz(124)
    assert 1-a == mpz(-122)
    assert (-1)-a == mpz(-124)
    assert a-b == mpz(-333)
    assert b-a == mpz(333)

    raises(TypeError, lambda: a-'b')
    raises(TypeError, lambda: 'b'-a)

    assert a-c == -12345678901234567767
    assert c-a == 12345678901234567767

    assert a - float('Inf') == mpfr('-inf')
    assert float('Inf') - a == mpfr('inf')
    assert a - float('-Inf') == mpfr('inf')
    assert float('-Inf') - a == mpfr('-inf')
    assert is_nan(a - float('nan'))
    assert is_nan(float('nan') - a)


def test_mpz_isub():
    x = mpz(7)
    x -= mpz(1)

    assert x == mpz(6)

    x -= 1

    assert x == mpz(5)

    x -= xmpz(7)

    assert x == mpz(-2)

    x -= -5

    assert x == mpz(3)

    x -= -mpfr(5)

    assert x == mpfr('8.0')


def test_mpz_mul():
    a = mpz(123)
    b = mpz(456)
    c = 12345678901234567890

    assert mpz(2)*z == mpz(4)
    assert gmpy2.mul(2,1) == mpz(2)

    ctx = gmpy2.context()

    assert ctx.mul(a,b) == a*b
    assert ctx.mul(c,c) == c*c
    assert ctx.mul(a, mpq(1)) == mpq(123,1)
    assert ctx.mul(a, mpfr(1)) == mpfr('123.0')
    assert ctx.mul(a, mpc(1)) == mpc('123.0+0.0j')

    raises(TypeError, lambda: ctx.mul(1))
    raises(TypeError, lambda: ctx.mul(1,2,3))

    assert a*b == mpz(56088)
    assert b*a == mpz(56088)
    assert a*0 == mpz(0)
    assert 0*a == mpz(0)
    assert a*123 == mpz(15129)
    assert 123*a == mpz(15129)
    assert a*c == 1518518504851851850470
    assert c*a == 1518518504851851850470

    a = mpz(3)

    assert a*'b' == 'bbb'
    assert 'b'*a == 'bbb'

    assert a*False == mpz(0)

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


def test_mul_imul():
    x = mpz(2)
    x *= mpz(2)

    assert x == mpz(4)

    x *= 2

    assert x == mpz(8)

    x *= xmpz(3)

    assert x == mpz(24)

    x *= -1

    assert x == mpz(-24)

    x *= mpfr(-0.5)

    assert x == mpfr('12.0')


def test_mpz_divmod():
    a = mpz(123)
    b = mpz(456)
    c = 12345678901234567890
    ctx = gmpy2.get_context()

    raises(TypeError, lambda: divmod(mpz(123),'a'))

    ctx = gmpy2.context()

    raises(TypeError, lambda: ctx.divmod('a',456))
    raises(TypeError, lambda: ctx.divmod(1,2,3))
    raises(ZeroDivisionError, lambda: ctx.divmod(456, 0))
    raises(TypeError, lambda: ctx.divmod(a,mpc(456)))
    raises(TypeError, lambda: divmod(mpz(1), mpc(1,2)))

    assert ctx.divmod(a,b) == (mpz(0), mpz(123))
    assert ctx.divmod(123,456) == (mpz(0), mpz(123))
    assert divmod(mpz(3), z) == (mpz(1), mpz(1))
    assert divmod(z, mpz(3)) == (mpz(0), mpz(2))

    assert divmod(a,b) == (mpz(0), mpz(123))
    assert divmod(b,a) == (mpz(3), mpz(87))

    raises(ZeroDivisionError, lambda: divmod(a,0))
    raises(ZeroDivisionError, lambda: divmod(a, mpz(0)))
    raises(ZeroDivisionError, lambda: divmod(123, mpz(0)))

    assert divmod(b,123) == (mpz(3), mpz(87))
    assert divmod(a,c) == (mpz(0), mpz(123))
    assert divmod(a,int(c)) == (mpz(0), mpz(123))
    assert divmod(a*(c-1),c) == (122, 12345678901234567767)
    assert divmod(a*(c-1),int(c)) == (122, 12345678901234567767)
    assert divmod(a*(c-1),-c) == (mpz(-123), mpz(-123))
    assert divmod(a*(c-1),-int(c)) == (mpz(-123), mpz(-123))
    assert divmod(int(a*(c-1)),-int(c)) == (-123, -123)

    assert divmod(a, mpfr('Inf')) == (mpfr('0.0'), mpfr('123.0'))
    assert divmod(a, mpfr('-Inf')) == (mpfr('-1.0'), mpfr('-inf'))
    assert divmod(-a, mpfr('Inf')) == (mpfr('-1.0'), mpfr('inf'))
    assert divmod(-a, mpfr('-Inf')) == (mpfr('0.0'), mpfr('-123.0'))
    assert all(is_nan(_) for _ in divmod(a, mpfr('nan')))
    assert all(is_nan(_) for _ in divmod(-a, mpfr('nan')))
    assert divmod(mpz(0), mpfr('Inf')) == (mpfr('0.0'), mpfr('0.0'))
    assert divmod(mpz(0), mpfr('-Inf')) == (mpfr('-0.0'), mpfr('-0.0'))
    assert divmod(mpz(0), mpfr('nan'))
    assert all(is_nan(_) for _ in divmod(mpfr('Inf'), a))
    assert all(is_nan(_) for _ in divmod(mpfr('-Inf'), a))
    assert all(is_nan(_) for _ in divmod(mpfr('Inf'), -a))
    assert all(is_nan(_) for _ in divmod(mpfr('-Inf'), -a))
    assert all(is_nan(_) for _ in divmod(mpfr('nan'), a))
    assert all(is_nan(_) for _ in divmod(mpfr('nan'), -a))
    assert all(is_nan(_) for _ in divmod(mpfr('Inf'), mpz(0)))
    assert all(is_nan(_) for _ in divmod(mpfr('-Inf'), mpz(0)))
    assert all(is_nan(_) for _ in divmod(mpfr('nan'), mpz(0)))

    assert divmod(a, mpfr('Inf')) == (mpfr('0.0'), mpfr('123.0'))
    assert divmod(a, mpfr('-Inf')) == (mpfr('-1.0'), mpfr('-inf'))
    assert divmod(-a, mpfr('Inf')) == (mpfr('-1.0'), mpfr('inf'))
    assert divmod(-a, mpfr('-Inf')) == (mpfr('0.0'), mpfr('-123.0'))
    assert all(is_nan(_) for _ in divmod(a, mpfr('nan')))
    assert all(is_nan(_) for _ in divmod(-a, mpfr('nan')))
    assert divmod(mpz(0), mpfr('Inf')) == (mpfr('0.0'), mpfr('0.0'))
    assert divmod(mpz(0), mpfr('-Inf')) == (mpfr('-0.0'), mpfr('-0.0'))
    assert divmod(mpz(0), mpfr('nan'))
    assert all(is_nan(_) for _ in divmod(mpfr('Inf'), a))
    assert all(is_nan(_) for _ in divmod(mpfr('-Inf'), a))
    assert all(is_nan(_) for _ in divmod(mpfr('Inf'), -a))
    assert all(is_nan(_) for _ in divmod(mpfr('-Inf'), -a))
    assert all(is_nan(_) for _ in divmod(mpfr('nan'), a))
    assert all(is_nan(_) for _ in divmod(mpfr('nan'), -a))
    assert all(is_nan(_) for _ in divmod(mpfr('Inf'), mpz(0)))
    assert all(is_nan(_) for _ in divmod(mpfr('-Inf'), mpz(0)))
    assert all(is_nan(_) for _ in divmod(mpfr('nan'), mpz(0)))


def test_mpz_floordiv():
    ctx = gmpy2.get_context()
    a, b = mpz(45), mpz(6)
    r, r2 = mpfr(45), mpfr(3.1)
    q, q2 = mpq(118,18), mpq(3,2)
    c, c2 = mpc(51, 65), mpc(4, 6)

    assert ctx.floor_div(a, 6) == mpz(7)
    assert ctx.floor_div(a, b) == mpz(7)

    raises(ZeroDivisionError, lambda: ctx.floor_div(a, 0))
    raises(ZeroDivisionError, lambda: ctx.floor_div(a, mpz(0)))

    assert ctx.floor_div(45, b) == mpz(7)

    raises(ZeroDivisionError, lambda: ctx.floor_div(45, 0))
    raises(ZeroDivisionError, lambda: ctx.floor_div(45, mpz(0)))
    raises(TypeError, lambda: ctx.floor_div())
    raises(TypeError, lambda: gmpy2.floor_div(4,5,6))

    assert a // b == mpz(7)

    raises(ZeroDivisionError, lambda: a // 0)
    raises(ZeroDivisionError, lambda: a // mpz(0))

    assert ctx.floor_div(a, q2) == mpz(30)
    assert ctx.floor_div(a, r2) == mpfr('14.0')

    raises(TypeError, lambda: ctx.floor_div(a, c))

    assert a // b == mpz(7)
    assert a // q == mpz(6)
    assert a // r2 == mpfr('14.0')

    raises(TypeError, lambda: a // c2)
    raises(TypeError, lambda: a // 'not')

    a = mpz(123)
    b = mpz(456)
    c = 12345678901234567890

    assert a//b == mpz(0)
    assert b//a == mpz(3)
    assert (a*b)//b == mpz(123)
    assert (a*b)//a == mpz(456)

    raises(ZeroDivisionError, lambda: a//0)

    assert c//a == 100371373180768844
    assert a**10//c == mpz(64)
    assert a // z == mpz(61)

    assert a//True == mpz(123)


def test_mpz_ifloordiv():
    x = mpz(49)
    x //= mpz(3)

    assert x == mpz(16)

    x //= xmpz(3)

    assert x == mpz(5)

    x //= 2

    assert x == mpz(2)

    with raises(ZeroDivisionError):
        x //= mpz(0)
    with raises(ZeroDivisionError):
        x //= 0

    assert x == mpz(2)

    x //= mpfr(-0.5)

    assert x == mpfr('-4.0')

    x = mpz(11)
    x //= -5

    assert x == mpz(-3)


def test_mpz_mod():
    a = mpz(123)
    b = mpz(456)
    c = 12345678901234567890
    ctx = gmpy2.get_context()

    assert a % b == mpz(123)
    assert b % a == mpz(87)
    assert gmpy2.mod(b, a) == mpz(87)
    assert ctx.mod(b, a) == mpz(87)
    assert a % z == mpz(1)

    raises(ZeroDivisionError, lambda: a % mpz(0))
    raises(ZeroDivisionError, lambda: a % 0)
    raises(ZeroDivisionError, lambda: 14 % mpz(0))
    raises(ZeroDivisionError, lambda: gmpy2.mod(14, mpz(0)))
    raises(ZeroDivisionError, lambda: gmpy2.mod(124, 0))
    raises(ZeroDivisionError, lambda: gmpy2.mod(b, mpz(0)))
    raises(ZeroDivisionError, lambda: ctx.mod(b, mpz(0)))
    raises(TypeError, lambda: gmpy2.mod(124, 'str'))
    raises(TypeError, lambda: a % 'str')

    assert gmpy2.mod(124, mpz(5)) == mpz(4)
    assert z % mpq(1,2) == mpq(0,1)
    assert a % mpq(2,3) == mpq(1,3)


def test_mpz_imod():
    x = mpz(45)
    x %= mpz(18)

    assert x == mpz(9)

    x %= xmpz(2)

    assert x == mpz(1)

    x = mpz(40)
    x %= 21

    assert x == mpz(19)

    with raises(ZeroDivisionError):
        x %= 0

    assert x == mpz(19)

    with raises(ZeroDivisionError):
        x %= mpz(0)

    x %= -9

    assert x == mpz(-8)

    x %= mpfr(10)

    assert x == mpfr('2.0')


def test_mpz_truediv():
    a = mpz(123)
    b = mpz(456)
    c = 12345678901234567890
    ctx = gmpy2.get_context()

    assert a/b == mpfr('0.26973684210526316')
    assert gmpy2.div(a, b) == mpfr('0.26973684210526316')
    assert ctx.div(a, b) == mpfr('0.26973684210526316')
    assert b/a == mpfr('3.7073170731707319')

    raises(ZeroDivisionError, lambda: a/0)

    assert a/0.0 == mpfr('inf')
    assert a / z == mpfr('61.5')

    raises(TypeError, lambda: ctx.div(a, b, 5))
    raises(TypeError, lambda: ctx.div(a, 'str'))
    raises(TypeError, lambda: a / 'str')

    with gmpy2.context(rational_division=True):
        assert mpz(1)/mpz(2) == mpq(1, 2)

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


def test_mpz_pow():
    z1, z2 = mpz(5), mpz(2)
    ctx = gmpy2.get_context()

    assert z1 ** z2 == mpz(25)
    assert ctx.pow(z1, z2) == mpz(25)
    assert z1 ** -z2 == mpfr('0.040000000000000001')
    assert z1 ** 0 == mpz(1)
    assert mpz(0) ** 32 == mpz(0)
    assert mpz(-1) ** 32 == mpz(1)
    assert mpz(1) ** mpz(324) == mpz(1)
    assert mpz(0) ** 0 == mpz(1)
    assert mpz(-1) ** 3 == mpz(-1)
    assert z1 ** 2 == pow(z1, 2)
    assert pow(z1, 2, 19) == mpz(6)
    assert pow(z1, -2, 19) == mpz(16)

    raises(ValueError, lambda: pow(mpz(0), -2, 19))

    assert pow(z1, 2, -19) == mpz(-13)

    raises(ValueError, lambda: pow(5, 2, 0))
    raises(TypeError, lambda: ctx.pow(z1, 'invalid'))

    a = mpz(123)
    b = mpz(456)

    assert pow(a,10) == 792594609605189126649
    assert pow(a,7,b) == mpz(99)


def test_mpz_ipow():
    x = mpz(5)
    x **= mpz(2)

    assert x == mpz(25)

    x **= xmpz(2)

    assert x == mpz(625)

    x **= -2

    assert x == mpfr('2.5600000000000001e-06')

    x = mpz(625)
    x **= 2

    assert x == mpz(390625)

    x **= mpfr(2)

    assert x == mpfr('152587890625.0')

    x = mpz(390625)
    x **= mpfr(-2)

    assert x == mpfr('6.5535999999999999e-12')


def test_lucasu():
    assert gmpy2.lucasu(2,4,0) == mpz(0)
    assert gmpy2.lucasu(2,4,1) == mpz(1)

    raises(ValueError, lambda: gmpy2.lucasu(2,1,1))

    assert gmpy2.lucasu(2,4,8) == mpz(128)

    raises(TypeError, lambda: gmpy2.lucasu('a',4,8))
    raises(TypeError, lambda: gmpy2.lucasu(2,'b',8))
    raises(TypeError, lambda: gmpy2.lucasu(2,4,None))
    raises(ValueError, lambda: gmpy2.lucasu(mpz(2), mpz(1), mpz(7)))

    assert gmpy2.lucasu_mod(3,2,0,7) == mpz(0)
    assert gmpy2.lucasu_mod(3,2,5,7) == mpz(3)
    assert gmpy2.lucasu_mod(3,2,555,777777777) == mpz(387104641)

    raises(ValueError, lambda: gmpy2.lucasu_mod(2,1,555,777777777))

    raises(ValueError, lambda: gmpy2.lucasv(2,1,4))
    raises(TypeError, lambda: gmpy2.lucasv('a',1,2))
    raises(TypeError, lambda: gmpy2.lucasv(4,'b',2))
    raises(TypeError, lambda: gmpy2.lucasv(4,3,'c'))

    assert gmpy2.lucasv(4,3,0) == mpz(2)
    assert gmpy2.lucasv(4,3,7) == mpz(2188)
    assert gmpy2.lucasv(4,3,8) == mpz(6562)

    assert gmpy2.lucasv_mod(4,3,0,2) == mpz(0)
    assert gmpy2.lucasv_mod(4,3,0,3) == mpz(2)
    assert gmpy2.lucasv_mod(4,3,55,123456) == mpz(35788)
    assert gmpy2.lucasv_mod(4,3,56,123456) == mpz(107362)
    assert gmpy2.lucasv_mod(4,3,57,123456) == mpz(75172)


def test_mpz_attributes():
    a = mpz(123)

    assert a.numerator == mpz(123)
    assert a.denominator == mpz(1)

    assert a.real == mpz(123)
    assert a.imag == mpz(0)


def test_mpz_conjugate():
    a = mpz(123)

    assert a.conjugate() == a


def test_mpz_invert():
    a = mpz(123)

    assert ~a == mpz(-124)


def test_mpz_and():
    a = mpz(123)
    b = mpz(456)

    assert a&b == mpz(72)
    assert a&int(b) == mpz(72)
    assert int(a)&b == mpz(72)

    raises(TypeError, lambda: a&mpq(1))


def test_mpz_and_operator():
    assert (mpz(0) and mpz(7)) == mpz(0)
    assert (mpz(7) and mpz(0)) == mpz(0)
    assert (mpz(7) and mpz(5)) == mpz(5)
    assert (mpz(0) and 5) == mpz(0)
    assert (mpz(7) and 5.2) == 5.2
    assert (mpz(7) and None) is None


def test_mpz_iand():
    x = mpz(7)
    x &= mpz(5)

    assert x == mpz(5)

    x &= xmpz(4)

    assert x == mpz(4)

    x &= 9

    assert x == mpz(0)

    x = mpz(4)
    x &= 12

    assert x == mpz(4)

    with raises(TypeError):
        x &= mpfr(4)


def test_mpz_or():
    a = mpz(123)
    b = mpz(456)

    assert a|b == mpz(507)
    assert a|int(b) == mpz(507)
    assert int(a)|b == mpz(507)

    raises(TypeError, lambda: a|mpq(1))


def test_mpz_ior():
    x = mpz(0)
    x |= mpz(1)

    assert x == mpz(1)

    x |= mpz(0)

    assert x == mpz(1)

    x = mpz(0)

    x |= mpz(0)

    assert x == mpz(0)

    x |= 5

    assert x == mpz(5)

    with raises(TypeError):
        x |= mpfr(3)


def test_mpz_xor():
    a = mpz(123)
    b = mpz(456)

    assert a^b == mpz(435)
    assert a^int(b) == mpz(435)
    assert int(a)^b == mpz(435)

    raises(TypeError, lambda: a^mpq(1))


def test_mpz_ixor():
    x = mpz(1)
    x ^= mpz(0)

    assert x == mpz(1)

    x ^= xmpz(1)

    assert x == mpz(0)

    x ^= 1

    assert x == mpz(1)

    with raises(TypeError):
        x ^= mpfr(0)


def test_mpz_lshift():
    a = mpz(123)

    assert a<<1 == mpz(246)
    assert int(a)<<mpz(1) == mpz(246)

    raises(OverflowError, lambda: a<<-1)

    assert a<<0 == mpz(123)

    raises(TypeError, lambda: "a" << a)


def test_mpz_rshift():
    a = mpz(123)

    assert a>>1 == mpz(61)
    assert int(a)>>mpz(1) == mpz(61)
    assert a>>111111111111111111111 == mpz(0)
    assert (-a)>>111111111111111111111 == mpz(-1)

    raises(ValueError, lambda: a>>-2)

    assert a>>0 == mpz(123)

    raises(TypeError, lambda: "a" >> a)


def test_mpz_ilshift_irshift():
    x = mpz(63)
    x >>= mpz(63)

    assert x == mpz(0)

    x = mpz(63)
    x >>= mpz(1)

    assert x == mpz(31)

    x >>= xmpz(2)

    assert x == mpz(7)

    x >>= 1

    assert x == mpz(3)

    x <<= mpz(2)

    assert x == mpz(12)

    x <<= mpz(1)

    assert x == mpz(24)

    x <<= 0

    assert x == mpz(24)

    with raises(TypeError):
        x >>= mpfr(2)
    with raises(TypeError):
        x <<= mpfr(2)
    with raises(ValueError):
        x >>= -1
    with raises(OverflowError):
        x <<= -5


def test_mpz_index():
    a = mpz(123)
    b = mpz(456)

    assert range(333)[a] == 123

    raises(IndexError, lambda: range(333)[b])


def test_mpz_seq():
    a = mpz(10)

    assert a[1] == 1
    assert a[-1] == 1
    assert a[-2] == 0
    assert a[0:2] == mpz(2)
    assert a[0:3] == mpz(2)
    assert a[0:3:-1] == mpz(0)

    raises(IndexError, lambda: a[111111111111111111111])
    raises(TypeError, lambda: a["spam"])


def test_mpz_len():
    assert len(mpz(0)) == 1
    assert len(mpz(1)) == 1
    assert len(mpz(17)) == 5


def test_mpz_num_digits():
    assert mpz(123456).num_digits() == 6
    assert mpz(123456).num_digits(2) == 17

    raises(ValueError, lambda: mpz(123456).num_digits(-1))
    raises(OverflowError, lambda: mpz(123456).num_digits(9999999999999999999999999999999999))
    raises(ValueError, lambda: mpz(123456).num_digits(100))
    raises(TypeError, lambda: gmpy2.num_digits(123456,-1,7))
    raises(ValueError, lambda: gmpy2.num_digits(123456,-1))
    raises(TypeError, lambda: gmpy2.num_digits('123456'))
    raises(ValueError, lambda: gmpy2.num_digits(123456,100))
    raises(TypeError, lambda: gmpy2.num_digits(123456,'a'))

    assert gmpy2.num_digits(123456) == 6
    assert gmpy2.num_digits(123456,2) == 17


def test_mpz_is_square():
    raises(TypeError, lambda: gmpy2.is_square('a'))

    assert gmpy2.is_square(mpz(9))
    assert not gmpy2.is_square(10)
    assert mpz(16).is_square()
    assert not mpz(17).is_square()


def test_mpz_is_integer():
    assert mpz(0).is_integer()
    assert mpz(123).is_integer()


def test_mpz_is_divisible():
    raises(TypeError, lambda: gmpy2.is_divisible())
    raises(TypeError, lambda: gmpy2.is_divisible('a',2))
    raises(TypeError, lambda: gmpy2.is_divisible(2,'a'))

    assert gmpy2.is_divisible(12,2)
    assert not gmpy2.is_divisible(12,7)
    assert mpz(12).is_divisible(2)
    assert not mpz(12).is_divisible(7)
    assert gmpy2.is_divisible(mpz(123456789123456789123456789),
                              123456789123456789123456789)
    assert not gmpy2.is_divisible(mpz(1234567891234567891234567897),
                                  123456789123456789123456789)

    raises(TypeError, lambda: mpz(12).is_divisible('a'))

    assert mpz(123456789123456789123456789).is_divisible(123456789123456789123456789)
    assert not mpz(1234567891234567891234567897).is_divisible(123456789123456789123456789)


def test_mpz_is_congruent():
    raises(TypeError, lambda: gmpy2.is_congruent(1))
    raises(TypeError, lambda: gmpy2.is_congruent(1,'a',3))

    assert gmpy2.is_congruent(7*3+2, 7*11+2, 7)
    assert not gmpy2.is_congruent(7*3+2, 7*11+5, 7)

    raises(TypeError, lambda: mpz(7*3+2).is_congruent(1))
    raises(TypeError, lambda: mpz(7*3+2).is_congruent('a',7))

    assert mpz(7*3+2).is_congruent(7*11+2,7)
    assert not mpz(7*3+2).is_congruent(7*11+5,7)


def test_mpz_is_power():
    a = mpz(123)

    raises(TypeError, lambda: gmpy2.is_power())
    raises(TypeError, lambda: gmpy2.is_power('a'))

    assert not a.is_power()
    assert mpz(123**11).is_power()
    assert not gmpy2.is_power(a)
    assert gmpy2.is_power(99*99*99)
    assert not gmpy2.is_power(99*98)


def test_mpz_is_prime():
    raises(OverflowError, lambda: gmpy2.is_prime(3,-3))
    raises(TypeError, lambda: gmpy2.is_prime())
    raises(TypeError, lambda: gmpy2.is_prime(1,2,3))
    raises(TypeError, lambda: gmpy2.is_prime('a'))

    assert not gmpy2.is_prime(12345)
    assert gmpy2.is_prime(80**81 + 81**80)
    assert gmpy2.is_prime(80**81 + 81**80, 10000)

    raises(TypeError, lambda: mpz(129).is_prime(1,2))
    raises(OverflowError, lambda: mpz(129).is_prime(-7))

    assert not mpz(129).is_prime(10000)
    assert mpz(80**81 + 81**80).is_prime()
    assert not mpz(1234567890).is_prime()
    assert not mpz(-3).is_prime()
    assert not gmpy2.is_prime(-3)


def test_mpz_is_probab_prime():
    raises(OverflowError, lambda: gmpy2.is_probab_prime(3,-3))
    raises(TypeError, lambda: gmpy2.is_probab_prime())
    raises(TypeError, lambda: gmpy2.is_probab_prime(1,2,3))
    raises(TypeError, lambda: gmpy2.is_probab_prime('a'))

    assert gmpy2.is_probab_prime(71) == 2
    assert gmpy2.is_probab_prime(12345) == 0
    assert gmpy2.is_probab_prime(80**81 + 81**80) == 1
    assert gmpy2.is_probab_prime(80**81 + 81**80, 10000) == 1

    raises(TypeError, lambda: mpz(129).is_probab_prime(1,2))
    raises(OverflowError, lambda: mpz(129).is_probab_prime(-7))

    assert mpz(71).is_probab_prime(71) == 2
    assert mpz(129).is_probab_prime(10000) == 0
    assert mpz(80**81 + 81**80).is_probab_prime() == 1
    assert mpz(1234567890).is_probab_prime() == 0
    assert mpz(-3).is_probab_prime() == 0
    assert gmpy2.is_probab_prime(-3) == 0


def test_mpz_is_even():
    a = mpz(123)
    b = mpz(456)

    raises(TypeError, lambda: gmpy2.is_even('a'))

    assert not gmpy2.is_even(a)
    assert gmpy2.is_even(b)
    assert not a.is_even()
    assert b.is_even()
    assert not gmpy2.is_even(11)
    assert gmpy2.is_even(14)


def test_mpz_is_odd():
    a = mpz(123)
    b = mpz(456)

    raises(TypeError, lambda: gmpy2.is_odd('a'))

    assert gmpy2.is_odd(a)
    assert not gmpy2.is_odd(b)
    assert a.is_odd()
    assert not b.is_odd()
    assert gmpy2.is_odd(11)
    assert not gmpy2.is_odd(14)


def test_mpz_bit_length():
    assert mpz(0).bit_length() == 0
    assert mpz(1).bit_length() == 1
    assert mpz(5).bit_length() == 3
    assert mpz(8).bit_length() == 4
    assert gmpy2.bit_length(mpz(10**30)) == 100
    assert gmpy2.bit_length(56) == 6

    raises(TypeError, lambda: gmpy2.bit_length(mpfr(4.0)))


def test_mpz_bit_mask():
    assert gmpy2.bit_mask(mpz(0)) == mpz(0)
    assert gmpy2.bit_mask(mpz(4)) == mpz(15)
    assert gmpy2.bit_mask(mpz(3)) == mpz(7)
    assert gmpy2.bit_mask(mpz(16)) == mpz(65535)
    assert gmpy2.bit_mask(8) == mpz(255)

    raises(OverflowError, lambda: gmpy2.bit_mask(-1))


def test_mpz_bit_scan0():
    assert mpz(6).bit_scan0() == 0
    assert mpz(7).bit_scan0() == 3
    assert mpz(8).bit_scan0(2) == 2
    assert mpz(7).bit_scan0(2) == 3

    raises(OverflowError, lambda: mpz(7).bit_scan0(-2))

    assert gmpy2.bit_scan0(mpz(7), 2) == 3
    assert gmpy2.bit_scan0(mpz(8), 2) == 2
    assert gmpy2.bit_scan0(8) == 0

    raises(TypeError, lambda: gmpy2.bit_scan0())
    raises(TypeError, lambda: gmpy2.bit_scan0(mpz(7), 2.5))
    raises(TypeError, lambda: gmpy2.bit_scan0(mpz(7), 2, 5))
    raises(TypeError, lambda: gmpy2.bit_scan0(7.5, 0))
    raises(OverflowError, lambda: gmpy2.bit_scan0(8, -2))

    assert gmpy2.bit_scan0(mpz(-1), 1) is None
    assert mpz(-1).bit_scan0(1) is None


def test_mpz_bit_scan1():
    assert mpz(7).bit_scan1() == 0
    assert mpz(8).bit_scan1() == 3
    assert mpz(7).bit_scan1(2) == 2

    raises(OverflowError, lambda: mpz(7).bit_scan1(-2))

    assert gmpy2.bit_scan1(7) == 0
    assert gmpy2.bit_scan1(8) == 3
    assert gmpy2.bit_scan1(7, 2) == 2

    raises(TypeError, lambda: gmpy2.bit_scan1(mpz(7), 2, 5))
    raises(TypeError, lambda: gmpy2.bit_scan1())
    raises(TypeError, lambda: gmpy2.bit_scan1(mpz(6), 2.5))
    raises(TypeError, lambda: gmpy2.bit_scan1(7.5, 0))
    raises(OverflowError, lambda: gmpy2.bit_scan1(8, -1))

    assert gmpy2.bit_scan1(mpz(1), 1) is None
    assert mpz(1).bit_scan1(1) is None


def test_mpz_bit_test():
    assert mpz(7).bit_test(2)
    assert not mpz(8).bit_test(2)
    assert not mpz(-8).bit_test(2)

    raises(OverflowError, lambda: mpz(8).bit_test(-2))

    assert gmpy2.bit_test(mpz(7), 2)
    assert not gmpy2.bit_test(mpz(8), 2)

    raises(TypeError, lambda: gmpy2.bit_test())
    raises(TypeError, lambda: gmpy2.bit_test(mpz(7), 2.5))
    raises(TypeError, lambda: gmpy2.bit_test(7.5, 2))
    raises(OverflowError, lambda: gmpy2.bit_test(8, -2))


def test_mpz_bit_clear():
    assert mpz(7).bit_clear(0) == mpz(6)
    assert mpz(7).bit_clear(2) == mpz(3)
    assert mpz(8).bit_clear(2) == mpz(8)

    raises(OverflowError, lambda: mpz(8).bit_clear(-1))

    assert gmpy2.bit_clear(4, 2) == mpz(0)

    raises(TypeError, lambda: gmpy2.bit_clear())
    raises(TypeError, lambda: gmpy2.bit_clear(7.2, 2))
    raises(TypeError, lambda: gmpy2.bit_clear(mpz(4), 2.5))
    raises(OverflowError, lambda: gmpy2.bit_clear(4, -2))


def test_mpz_bit_set():
    assert mpz(4).bit_set(0) == mpz(5)
    assert mpz(7).bit_set(3) == mpz(15)
    assert mpz(0).bit_set(2) == mpz(4)

    raises(OverflowError, lambda: mpz(0).bit_set(-2))

    assert gmpy2.bit_set(8, 1) == mpz(10)

    raises(TypeError, lambda: gmpy2.bit_set(0))
    raises(TypeError, lambda: gmpy2.bit_set())
    raises(TypeError, lambda: gmpy2.bit_set(8.5, 1))
    raises(TypeError, lambda: gmpy2.bit_set(8, 1.5))
    raises(OverflowError, lambda: gmpy2.bit_set(8, -1))


def test_mpz_bit_flip():
    assert mpz(4).bit_flip(2) == mpz(0)
    assert mpz(4).bit_flip(1) == mpz(6)
    assert mpz(0).bit_flip(3) == mpz(8)

    raises(OverflowError, lambda: mpz(5).bit_flip(-3))

    assert gmpy2.bit_flip(mpz(7), mpz(1)) == mpz(5)
    assert gmpy2.bit_flip(mpz(7), 2) == mpz(3)

    raises(TypeError, lambda: gmpy2.bit_flip())
    raises(TypeError, lambda: gmpy2.bit_flip(4.5, 2))
    raises(TypeError, lambda: gmpy2.bit_flip(4, 2.5))
    raises(OverflowError, lambda: gmpy2.bit_flip(mpz(7), -2))


def test_mpz_bit_count():
    a = mpz(10009)

    assert a.bit_count() == 7
    assert gmpy2.popcount(a) == 7
    assert gmpy2.bit_count(a) == 7

    a = -a

    assert a == mpz(-10009)
    assert a.bit_count() == 7
    assert gmpy2.bit_count(a) == 7
    assert gmpy2.popcount(a) == -1

    raises(TypeError, lambda: gmpy2.bit_count('spam'))


def test_mpz_popcount():
    assert gmpy2.popcount(-65) == -1
    assert gmpy2.popcount(7) == 3
    assert gmpy2.popcount(8) == 1
    assert gmpy2.popcount(15) == 4
    assert gmpy2.popcount(mpz(0)) == 0

    raises(TypeError, lambda: gmpy2.popcount(4.5))


def test_mpz_hamdist():
    assert gmpy2.hamdist(mpz(5), mpz(7)) == 1
    assert gmpy2.hamdist(mpz(0), mpz(7)) == 3
    assert gmpy2.hamdist(mpz(0), 7) == 3

    raises(TypeError, lambda: gmpy2.hamdist(mpq(14,2), 5))
    raises(TypeError, lambda: gmpy2.hamdist(5,6,5))


def test_issue_339():
    samples = map(mpz, [13157547707030902665, 1070317427780135395,
                        18019609787501108695, 3978762157568107671,
                        14444587867185512177])
    assert all((2*q).is_divisible(q) for q in samples)


def test_issue_312():
    assert not is_prime(-7)
    assert not is_prime(mpz(-7))
    assert not is_prime(1 - 2**4423)
    assert all(not is_prime(-a) for a in range(8))
    assert next_prime(-8) == 2


def test_mpz_array():
    numpy = pytest.importorskip('numpy')
    i = 5579686107214117131790972086716881
    m = gmpy2.mpz(i)
    assert numpy.longdouble(m) == numpy.longdouble(i)
    assert m.__array__(dtype=numpy.longdouble) == numpy.longdouble(i)

    raises(TypeError, lambda: m.__array__(1, 2, 3))
    raises(TypeError, lambda: m.__array__(dtype=None, copy=None, spam=123))
    raises(TypeError, lambda: m.__array__(int, dtype=None))
    raises(TypeError, lambda: m.__array__(int, None, copy=None))
    raises(TypeError, lambda: m.__array__(spam=123))


def test_mpz_thread_safe():
    def worker():
        ctx = gmpy2.get_context()
        a = mpz(-1)
        b = abs(a) + 2
        a = b - 2
        assert str(a/b) == '0.33333333333333331'
        del a
        ctx.rational_division = True
        a = b - 2
        assert str(a/b) == '1/3'
        del b
    tpe = ThreadPoolExecutor(max_workers=20)
    for _ in range(1000):
        tpe.submit(worker)


@pytest.mark.skipif(sys.version_info < (3, 13), reason="requires v3.13+")
def test_mpz_signatures():
    cls = gmpy2.mpz
    for f in dir(cls):
        a = getattr(cls, f)
        if callable(a) and f != '__class__':
            _ = inspect.signature(a)  # not raises
